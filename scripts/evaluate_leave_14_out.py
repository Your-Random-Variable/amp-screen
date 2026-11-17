#!/usr/bin/env python3
"""
AMP-Screen: Leave-14-Out Reference Cross-Validation Procedure.
Evaluates out-of-sample generalization when all 14 standardized experimental
co-crystal complexes are strictly held out from ensemble training.
Saves row-level predictions to benchmarks/leave_14_out_results.csv.
"""

import sys
from pathlib import Path

# Add project root to path for standalone execution
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
import catboost as cb
import lightgbm as lgb

from amp_screen.ml_scorer import FEATURE_NAMES, RT_LN10
from amp_screen.utils import get_default_dataset_path


BENCHMARK_14_PDBS = [
    "1AWR", "1B07", "1BBZ", "1BE9", "1CKA", "1DKX", "1ER8",
    "1GBQ", "1LCJ", "1PRL", "1SEM", "1SSH", "1VJ6", "1YCR"
]


def evaluate_leave_14_out(
    dataset_path: Path = None,
    out_csv: Path = Path("benchmarks/leave_14_out_results.csv")
) -> pd.DataFrame:
    if dataset_path is None:
        dataset_path = get_default_dataset_path()
        if dataset_path is None or not Path(dataset_path).exists():
            raise FileNotFoundError("Could not find default dataset path.")

    df = pd.read_csv(dataset_path)
    print("=" * 100)
    print("      AMP-SCREEN: LEAVE-14-OUT REFERENCE CROSS-VALIDATION EVALUATION")
    print("=" * 100)
    print(f"Loaded training corpus from: {dataset_path} ({len(df)} samples)")

    is_holdout = df["pdb_id"].str.upper().isin(set(BENCHMARK_14_PDBS))
    train_df = df[~is_holdout].copy()
    test_df = df[is_holdout].copy()

    features = [f for f in FEATURE_NAMES if f in df.columns]
    assert len(features) == 96, f"Expected 96 features, found {len(features)}"

    X_train = train_df[features].values
    y_train = train_df["experimental_pkd"].values

    X_test = test_df[features].values
    y_test = test_df["experimental_pkd"].values

    print(f"Training set: {len(train_df)} non-reference complexes")
    print(f"Held-out test set: {len(test_df)} reference co-crystals\n")

    # Fit sub-models with standard ensemble hyperparameters
    rf = RandomForestRegressor(n_estimators=300, max_features=0.4, min_samples_leaf=2, random_state=42, n_jobs=-1).fit(X_train, y_train)
    et = ExtraTreesRegressor(n_estimators=300, max_depth=15, max_features=0.4, min_samples_leaf=2, random_state=42, n_jobs=-1).fit(X_train, y_train)
    cat = cb.CatBoostRegressor(iterations=350, depth=6, l2_leaf_reg=5.0, learning_rate=0.03, verbose=0, random_seed=42).fit(X_train, y_train)
    lgbm = lgb.LGBMRegressor(n_estimators=250, num_leaves=20, colsample_bytree=0.5, learning_rate=0.03, random_state=42, verbose=-1).fit(X_train, y_train)

    # Weighted ensemble prediction
    pred_pkd = 0.35 * rf.predict(X_test) + 0.25 * et.predict(X_test) + 0.20 * cat.predict(X_test) + 0.20 * lgbm.predict(X_test)

    # Thermodynamics conversion at 298.15 K
    pred_dg = -RT_LN10 * pred_pkd
    exp_dg = -RT_LN10 * y_test

    test_df["predicted_pkd"] = np.round(pred_pkd, 2)
    test_df["predicted_delta_g"] = np.round(pred_dg, 2)
    test_df["exp_delta_g"] = np.round(exp_dg, 2)
    test_df["abs_error_kcal"] = np.round(np.abs(pred_dg - exp_dg), 2)

    # Filter only quantitative comparisons (12 systems with known exp Kd)
    quant_df = test_df.dropna(subset=["experimental_pkd", "predicted_pkd"]).copy()
    mae = float(np.mean(quant_df["abs_error_kcal"]))
    rmse = float(np.sqrt(np.mean(quant_df["abs_error_kcal"] ** 2)))
    r_val, p_val = pearsonr(quant_df["predicted_delta_g"], quant_df["exp_delta_g"])
    rho_val, p_rho = spearmanr(quant_df["predicted_delta_g"], quant_df["exp_delta_g"])

    fmt = "%-6s | %-20s | %-12s | %-8s | %-8s | %-8s | %-8s"
    print(fmt % ("PDB", "Receptor Family", "Peptide Seq", "Exp pKd", "Pred pKd", "Exp dG", "Pred dG"))
    print("-" * 90)
    for _, r in test_df.iterrows():
        print(fmt % (
            r["pdb_id"], str(r["receptor_family"])[:20], str(r["peptide_sequence"])[:12],
            f"{r['experimental_pkd']:.2f}", f"{r['predicted_pkd']:.2f}",
            f"{r['exp_delta_g']:.2f}", f"{r['predicted_delta_g']:.2f}"
        ))
    print("=" * 90)

    print(f"\nLeave-14-Out Independent Generalization Metrics (N={len(quant_df)}):")
    print(f"  • MAE:  {mae:.2f} kcal/mol")
    print(f"  • RMSE: {rmse:.2f} kcal/mol")
    print(f"  • Pearson r:  {r_val:.3f} (p = {p_val:.2e})")
    print(f"  • Spearman ρ: {rho_val:.3f} (p = {p_rho:.2e})")
    print("=" * 90 + "\n")

    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cols_to_save = [
        "pdb_id", "receptor_family", "peptide_sequence", "assay_type",
        "experimental_pkd", "predicted_pkd", "exp_delta_g", "predicted_delta_g",
        "abs_error_kcal"
    ]
    test_df[cols_to_save].to_csv(out_csv, index=False)
    print(f"[OK] Leave-14-out predictions saved to: {out_csv}")
    return test_df


if __name__ == "__main__":
    evaluate_leave_14_out()
