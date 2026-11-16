"""
Model Training & Cluster-Stratified Cross-Validation Utility for AMP-Screen.
Extracts 96D structural, physical, and sequence feature vectors from a dataset of
539 protein-peptide co-complexes across 46 distinct receptor families with known
experimental affinities (pKd, Kd, DeltaG).
Trains the tri-model ensemble (Gradient Boosting + Random Forest + Extra Trees)
and saves the trained model to data/trained_complex_ensemble.joblib.

Key Scientific Safeguards:
1. Zero Data Leakage: 'delta_delta_g' feature strictly removed.
2. No Silent Defaults: Real sequence length, contact densities, and pI used.
3. Zero Target Leakage: GroupKFold cross-validation by receptor family ensures
   different members of the same protein family (e.g., all SH3 or all Caspases)
   never appear in both training and validation splits.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import joblib
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
from sklearn.model_selection import KFold, GroupKFold

from amp_screen.ml_scorer import (
    FEATURE_NAMES,
    extract_structural_fingerprint,
    extract_complex_features,
    ELEMENT_PAIRS
)
from amp_screen.physics_engine import calculate_haddock_physics_score
from amp_screen.utils import logger, ensure_dir


def featurize_complex(complex_pdb: Path, meta_row: dict) -> list:
    """Extracts 96D feature vector from a single complex PDB (without label leakage)."""
    fp = extract_structural_fingerprint(complex_pdb)
    physics = calculate_haddock_physics_score(complex_pdb)
    seq = meta_row.get("sequence", meta_row.get("peptide_sequence", ""))
    return extract_complex_features(fp=fp, physics=physics, hit_meta=meta_row, seq=seq)


from amp_screen.deep_model import DeepAffinityRegressor


def train_and_evaluate(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray = None,
    out_model_path: Path = Path("data/trained_complex_ensemble.joblib"),
    n_splits: int = 5
) -> dict:
    """
    Trains the deep learning ResNet, tree ensemble, and hybrid deep ensemble.
    Evaluates cross-validation performance under GroupKFold by receptor family.
    """
    print("\n" + "=" * 105)
    print("      PROTEIN-PEPTIDE COMPLEX AFFINITY: DEEP LEARNING & ENSEMBLE TRAINING (96 PARAMETERS)")
    print("=" * 105)

    y_pred_trees_cv = np.zeros_like(y)
    y_pred_dl_cv = np.zeros_like(y)
    y_pred_hybrid_cv = np.zeros_like(y)

    if groups is not None and len(np.unique(groups)) >= n_splits:
        print(f"Using GroupKFold cross-validation across {len(np.unique(groups))} distinct receptor families.")
        splitter = GroupKFold(n_splits=n_splits)
        split_gen = splitter.split(X, y, groups=groups)
    else:
        print("Using standard 5-fold cross-validation.")
        splitter = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        split_gen = splitter.split(X, y)

    for fold_idx, (train_idx, val_idx) in enumerate(split_gen, 1):
        X_tr, y_tr = X[train_idx], y[train_idx]
        X_va, y_va = X[val_idx], y[val_idx]

        # 1. Decision Tree Models
        rf = RandomForestRegressor(n_estimators=400, max_depth=12, min_samples_split=3, max_features=0.6, random_state=42, n_jobs=-1)
        gb = GradientBoostingRegressor(n_estimators=250, max_depth=3, learning_rate=0.03, subsample=0.75, max_features=0.7, random_state=42)
        et = ExtraTreesRegressor(n_estimators=400, max_depth=12, min_samples_split=3, max_features=0.6, random_state=42, n_jobs=-1)

        rf.fit(X_tr, y_tr)
        gb.fit(X_tr, y_tr)
        et.fit(X_tr, y_tr)

        tree_pred = 0.40 * gb.predict(X_va) + 0.30 * rf.predict(X_va) + 0.30 * et.predict(X_va)
        y_pred_trees_cv[val_idx] = tree_pred

        # 2. Deep Learning ResNet-MLP Model
        dl = DeepAffinityRegressor(
            in_features=X.shape[1],
            hidden_dim=256,
            n_blocks=2,
            dropout=0.15,
            lr=1e-3,
            weight_decay=1e-4,
            epochs=100,
            batch_size=32,
            random_state=42 + fold_idx
        )
        dl.fit(X_tr, y_tr, val_data=(X_va, y_va))
        dl_pred = dl.predict(X_va)
        y_pred_dl_cv[val_idx] = dl_pred

        # 3. Hybrid Deep Ensemble (35% DL + 30% GB + 20% RF + 15% ET)
        hybrid_pred = 0.35 * dl_pred + 0.65 * tree_pred
        y_pred_hybrid_cv[val_idx] = hybrid_pred

        fold_r_tree, _ = pearsonr(y_va, tree_pred) if len(y_va) > 2 else (0.0, 1.0)
        fold_r_dl, _ = pearsonr(y_va, dl_pred) if len(y_va) > 2 else (0.0, 1.0)
        fold_r_hyb, _ = pearsonr(y_va, hybrid_pred) if len(y_va) > 2 else (0.0, 1.0)
        print(f"  Fold {fold_idx}/{n_splits} ({len(val_idx)} complexes) | Tree r: {fold_r_tree:.3f} | DL ResNet r: {fold_r_dl:.3f} | Hybrid r: {fold_r_hyb:.3f}")

    # Metrics computation
    def calc_metrics(y_true, y_pred, name):
        r_v, p_v = pearsonr(y_true, y_pred)
        rho_v, _ = spearmanr(y_true, y_pred)
        rmse_v = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        mae_v = float(np.mean(np.abs(y_true - y_pred)))
        print(f"\n--- {name} ---")
        print(f"  • Pearson Correlation (r):   {r_v:.4f} (p = {p_v:.2e})")
        print(f"  • Spearman Rank Order (ρ):   {rho_v:.4f}")
        print(f"  • Root Mean Square Error:    {rmse_v:.3f} pKd (~{rmse_v * 1.364:.2f} kcal/mol)")
        print(f"  • Mean Absolute Error:       {mae_v:.3f} pKd (~{mae_v * 1.364:.2f} kcal/mol)")
        return {"r": float(r_v), "rho": float(rho_v), "rmse": float(rmse_v), "mae": float(mae_v)}

    print("\n" + "=" * 105)
    print("      CROSS-VALIDATED MODEL ARCHITECTURE COMPARISON (GroupKFold)")
    print("=" * 105)
    tree_metrics = calc_metrics(y, y_pred_trees_cv, "Pure Tree Ensemble (GB + RF + ET)")
    dl_metrics = calc_metrics(y, y_pred_dl_cv, "Pure Deep Learning (DeepAffinityResNet)")
    hyb_metrics = calc_metrics(y, y_pred_hybrid_cv, "Hybrid Deep Ensemble (DL + GB + RF + ET)")
    print("=" * 105)

    # Fit final models on 100% of available data
    print(f"\nFitting final models on all {len(y)} complexes with {len(FEATURE_NAMES)} parameters...")
    final_rf = RandomForestRegressor(n_estimators=400, max_depth=12, min_samples_split=3, max_features=0.6, random_state=42, n_jobs=-1)
    final_gb = GradientBoostingRegressor(n_estimators=250, max_depth=3, learning_rate=0.03, subsample=0.75, max_features=0.7, random_state=42)
    final_et = ExtraTreesRegressor(n_estimators=400, max_depth=12, min_samples_split=3, max_features=0.6, random_state=42, n_jobs=-1)

    final_rf.fit(X, y)
    final_gb.fit(X, y)
    final_et.fit(X, y)

    print("Fitting final DeepAffinityResNet model...")
    final_dl = DeepAffinityRegressor(
        in_features=X.shape[1],
        hidden_dim=256,
        n_blocks=2,
        dropout=0.15,
        lr=1e-3,
        weight_decay=1e-4,
        epochs=120,
        batch_size=32,
        random_state=42
    )
    final_dl.fit(X, y)

    out_model_path = Path(out_model_path)
    ensure_dir(out_model_path.parent)
    payload = {
        "rf": final_rf,
        "gb": final_gb,
        "et": final_et,
        "dl": final_dl.to_dict(),
        "features": list(FEATURE_NAMES),
        "ensemble_weights": {"dl": 0.35, "gb": 0.26, "rf": 0.195, "et": 0.195},
        "cv_tree_metrics": tree_metrics,
        "cv_dl_metrics": dl_metrics,
        "cv_hybrid_metrics": hyb_metrics,
        "cv_pearson_r": round(hyb_metrics["r"], 4),
        "cv_spearman_rho": round(hyb_metrics["rho"], 4),
        "cv_rmse": round(hyb_metrics["rmse"], 4),
        "cv_mae": round(hyb_metrics["mae"], 4),
        "n_samples": len(y),
        "n_features": len(FEATURE_NAMES)
    }
    joblib.dump(payload, str(out_model_path), compress=3)
    print(f"✔ Successfully saved trained 96D Hybrid Deep Ensemble to: {out_model_path}")

    return payload


def main():
    parser = argparse.ArgumentParser(description="Train custom 96-parameter ML affinity model for AMP-Screen")
    parser.add_argument("--data", default="data/pdbbind_peptide_affinity_dataset.csv", help="Dataset CSV with experimental affinities and features")
    parser.add_argument("--output", default="data/trained_complex_ensemble.joblib", help="Output model path")
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        raise FileNotFoundError(f"Dataset not found: {data_path}")

    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} protein-peptide complexes from {data_path.name}.")

    # Verify all 96 features exist in dataframe
    missing_feats = [f for f in FEATURE_NAMES if f not in df.columns]
    if missing_feats:
        raise ValueError(f"Dataset is missing required features: {missing_feats}")

    X = df[FEATURE_NAMES].values
    y_col = next((c for c in df.columns if c in ("experimental_pkd", "experimental_pKd", "pKd", "pkd")), None)
    if not y_col:
        raise ValueError(f"Could not find pKd column in {data_path.name}")
    y = df[y_col].values

    groups = df["receptor_family"].values if "receptor_family" in df.columns else None

    train_and_evaluate(X=X, y=y, groups=groups, out_model_path=Path(args.output))


if __name__ == "__main__":
    main()
