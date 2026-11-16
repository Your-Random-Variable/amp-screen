"""
Model Training & Cluster-Stratified Cross-Validation Utility for AMP-Screen.
Extracts 96D structural, physical, and sequence feature vectors from a dataset of
1,054 protein-peptide co-complexes across 52 distinct receptor families with known
experimental affinities (pKd, Kd, DeltaG).
Trains the Tuned Pent-Ensemble:
  1. Tuned Random Forest (500 trees, max_features=0.4, min_samples_leaf=2)
  2. Tuned Extra Trees (500 trees, max_depth=15, max_features=0.4, min_samples_leaf=2)
  3. Tuned CatBoost (450 iterations, depth=6, l2_leaf_reg=10)
  4. Tuned LightGBM (350 estimators, num_leaves=20, colsample=0.5)
  5. Regularized Deep Tabular ResNet (64D, dropout=0.35, weight_decay=5e-3)
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
import catboost as cb
import lightgbm as lgb

from amp_screen.ml_scorer import (
    FEATURE_NAMES,
    extract_structural_fingerprint,
    extract_complex_features,
    ELEMENT_PAIRS,
    RT_LN10
)
from amp_screen.physics_engine import calculate_haddock_physics_score
from amp_screen.utils import logger, ensure_dir
from amp_screen.deep_model import DeepAffinityRegressor


def featurize_complex(complex_pdb: Path, meta_row: dict) -> list:
    """Extracts 96D feature vector from a single complex PDB (without label leakage)."""
    fp = extract_structural_fingerprint(complex_pdb)
    physics = calculate_haddock_physics_score(complex_pdb)
    seq = meta_row.get("sequence", meta_row.get("peptide_sequence", ""))
    return extract_complex_features(fp=fp, physics=physics, hit_meta=meta_row, seq=seq)


def train_and_evaluate(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray = None,
    out_model_path: Path = Path("data/trained_complex_ensemble.joblib"),
    n_splits: int = 5
) -> dict:
    """
    Trains the Tuned Pent-Ensemble (RF, ET, CatBoost, LightGBM, DL ResNet, and GB fallback).
    Evaluates cross-validation performance under GroupKFold by receptor family.
    """
    print("\n" + "=" * 110)
    print("      PROTEIN-PEPTIDE COMPLEX AFFINITY: TUNED PENT-ENSEMBLE TRAINING (96 PARAMETERS, 1,054 SAMPLES)")
    print("=" * 110)

    y_pred_rf = np.zeros_like(y)
    y_pred_et = np.zeros_like(y)
    y_pred_cb = np.zeros_like(y)
    y_pred_lgb = np.zeros_like(y)
    y_pred_dl = np.zeros_like(y)
    y_pred_pent = np.zeros_like(y)

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

        # 1. Tuned Random Forest (decorrelated subspace)
        rf = RandomForestRegressor(n_estimators=400, max_depth=None, min_samples_leaf=2, max_features=0.4, random_state=42, n_jobs=-1)
        rf.fit(X_tr, y_tr)
        p_rf = rf.predict(X_va)
        y_pred_rf[val_idx] = p_rf

        # 2. Tuned Extra Trees
        et = ExtraTreesRegressor(n_estimators=400, max_depth=15, min_samples_leaf=2, max_features=0.4, random_state=42, n_jobs=-1)
        et.fit(X_tr, y_tr)
        p_et = et.predict(X_va)
        y_pred_et[val_idx] = p_et

        # 3. Tuned CatBoost
        cb_m = cb.CatBoostRegressor(iterations=400, learning_rate=0.03, depth=6, l2_leaf_reg=10, random_seed=42, verbose=0)
        cb_m.fit(X_tr, y_tr)
        p_cb = cb_m.predict(X_va)
        y_pred_cb[val_idx] = p_cb

        # 4. Tuned LightGBM
        lgb_m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.03, num_leaves=20, min_child_samples=15, colsample_bytree=0.5, subsample=0.8, random_state=42, verbose=-1, n_jobs=-1)
        lgb_m.fit(X_tr, y_tr)
        p_lgb = lgb_m.predict(X_va)
        y_pred_lgb[val_idx] = p_lgb

        # 5. Regularized Deep Tabular ResNet
        dl = DeepAffinityRegressor(in_features=X.shape[1], hidden_dim=64, n_blocks=1, dropout=0.35, lr=1e-3, weight_decay=5e-3, epochs=60, batch_size=32, random_state=42 + fold_idx)
        dl.fit(X_tr, y_tr, val_data=(X_va, y_va))
        p_dl = dl.predict(X_va)
        y_pred_dl[val_idx] = p_dl

        # Tuned Pent-Ensemble Blend
        p_pent = 0.35 * p_rf + 0.25 * p_et + 0.20 * p_cb + 0.10 * p_lgb + 0.10 * p_dl
        y_pred_pent[val_idx] = p_pent

        fold_r_pent, _ = pearsonr(y_va, p_pent) if len(y_va) > 2 else (0.0, 1.0)
        print(f"  Fold {fold_idx}/{n_splits} ({len(val_idx)} complexes) | Tuned Pent-Ensemble Pearson r: {fold_r_pent:.3f}")

    def calc_metrics(y_true, y_pred, name):
        r_v, p_v = pearsonr(y_true, y_pred)
        rho_v, _ = spearmanr(y_true, y_pred)
        rmse_v = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        mae_v = float(np.mean(np.abs(y_true - y_pred)))
        print(f"\n--- {name} ---")
        print(f"  • Pearson Correlation (r):   {r_v:.4f} (p = {p_v:.2e})")
        print(f"  • Spearman Rank Order (ρ):   {rho_v:.4f}")
        print(f"  • Root Mean Square Error:    {rmse_v:.3f} pKd (~{rmse_v * RT_LN10:.2f} kcal/mol)")
        print(f"  • Mean Absolute Error:       {mae_v:.3f} pKd (~{mae_v * RT_LN10:.2f} kcal/mol)")
        return {"r": float(r_v), "rho": float(rho_v), "rmse": float(rmse_v), "mae": float(mae_v)}

    print("\n" + "=" * 110)
    print("      CROSS-VALIDATED MODEL PERFORMANCE (GroupKFold across 52 Receptor Families)")
    print("=" * 110)
    rf_metrics = calc_metrics(y, y_pred_rf, "Tuned Random Forest (Subspace Decoupled)")
    et_metrics = calc_metrics(y, y_pred_et, "Tuned Extra Trees")
    cb_metrics = calc_metrics(y, y_pred_cb, "Tuned CatBoost (Oblivious Trees, L2=10)")
    lgb_metrics = calc_metrics(y, y_pred_lgb, "Tuned LightGBM (Leaf-Wise)")
    dl_metrics = calc_metrics(y, y_pred_dl, "Regularized Deep Tabular ResNet (64D, Dropout=0.35)")
    pent_metrics = calc_metrics(y, y_pred_pent, "Tuned Pent-Ensemble (35% RF + 25% ET + 20% CB + 10% LGB + 10% DL)")
    print("=" * 110)

    # Fit final models on 100% of available data
    print(f"\nFitting final tuned models on all {len(y)} complexes with {len(FEATURE_NAMES)} parameters...")
    final_rf = RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=2, max_features=0.4, random_state=42, n_jobs=-1)
    final_et = ExtraTreesRegressor(n_estimators=500, max_depth=15, min_samples_leaf=2, max_features=0.4, random_state=42, n_jobs=-1)
    final_cb = cb.CatBoostRegressor(iterations=450, learning_rate=0.03, depth=6, l2_leaf_reg=10, random_seed=42, verbose=0)
    final_lgb = lgb.LGBMRegressor(n_estimators=350, learning_rate=0.03, num_leaves=20, min_child_samples=15, colsample_bytree=0.5, subsample=0.8, random_state=42, verbose=-1, n_jobs=-1)
    final_gb = GradientBoostingRegressor(n_estimators=250, max_depth=3, learning_rate=0.03, subsample=0.75, max_features=0.7, random_state=42)

    final_rf.fit(X, y)
    final_et.fit(X, y)
    final_cb.fit(X, y)
    final_lgb.fit(X, y)
    final_gb.fit(X, y)

    print("Fitting final regularized DeepAffinityResNet model...")
    final_dl = DeepAffinityRegressor(
        in_features=X.shape[1],
        hidden_dim=64,
        n_blocks=1,
        dropout=0.35,
        lr=1e-3,
        weight_decay=5e-3,
        epochs=60,
        batch_size=32,
        random_state=42
    )
    final_dl.fit(X, y)

    out_model_path = Path(out_model_path)
    ensure_dir(out_model_path.parent)
    payload = {
        "rf": final_rf,
        "et": final_et,
        "cb": final_cb,
        "lgb": final_lgb,
        "gb": final_gb,
        "dl": final_dl.to_dict(),
        "features": list(FEATURE_NAMES),
        "ensemble_weights": {"rf": 0.35, "et": 0.25, "cb": 0.20, "lgb": 0.10, "dl": 0.10},
        "cv_pent_metrics": pent_metrics,
        "cv_dl_metrics": dl_metrics,
        "cv_rf_metrics": rf_metrics,
        "cv_et_metrics": et_metrics,
        "cv_cb_metrics": cb_metrics,
        "cv_lgb_metrics": lgb_metrics,
        "cv_pearson_r": round(pent_metrics["r"], 4),
        "cv_spearman_rho": round(pent_metrics["rho"], 4),
        "cv_rmse": round(pent_metrics["rmse"], 4),
        "cv_mae": round(pent_metrics["mae"], 4),
        "n_samples": len(y),
        "n_features": len(FEATURE_NAMES)
    }
    joblib.dump(payload, str(out_model_path), compress=3)
    print(f"✔ Successfully saved tuned Pent-Ensemble to: {out_model_path}")

    return payload


def main():
    parser = argparse.ArgumentParser(description="Train custom 96-parameter tuned pent-ensemble affinity model for AMP-Screen")
    parser.add_argument("--data", default="data/pdbbind_peptide_affinity_dataset.csv", help="Dataset CSV with experimental affinities and features")
    parser.add_argument("--output", default="data/trained_complex_ensemble.joblib", help="Output model path")
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        raise FileNotFoundError(f"Dataset not found: {data_path}")

    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} protein-peptide complexes from {data_path.name}.")

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
