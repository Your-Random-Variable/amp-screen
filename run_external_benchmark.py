#!/usr/bin/env python3
"""
AMP-Screen External Benchmark Evaluation Suite.
Strictly evaluates the physics engine and 56D ML ensemble on an external,
prospective test set of 8 high-resolution protein-peptide co-crystal structures
that have ZERO overlap with the 381 PDBbind training dataset.
"""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import gemmi

from amp_screen.physics_engine import calculate_haddock_physics_score
from amp_screen.amp_analyzer import (
    analyze_cocomplex_contacts,
    estimate_solution_affinity,
    calculate_net_charge,
    calculate_hydrophobic_moment,
    calculate_boman_index
)
from amp_screen.utils import AA_3TO1, logger
from amp_screen.ml_scorer import MLPeptideScorer

EXTERNAL_BENCHMARKS = {
    "1X2R": {
        "target": "Keap1 Kelch Domain",
        "pep_name": "Nrf2 ETGE (9-mer)",
        "exp_kd_molar": 5.0e-9,
        "exp_kd_str": "5.0 nM",
        "exp_delta_g": -11.32,
        "citation": "Tong et al., Mol. Cell. Biol. 26:2887 (2006)"
    },
    "2FLU": {
        "target": "Keap1 Kelch Domain",
        "pep_name": "Nrf2 Neh2 (16-mer)",
        "exp_kd_molar": 9.0e-9,
        "exp_kd_str": "9.0 nM",
        "exp_delta_g": -10.97,
        "citation": "Lo et al., EMBO J. 25:3605 (2006)"
    },
    "2NL9": {
        "target": "Mcl-1 Oncoprotein",
        "pep_name": "Bim BH3 (23-mer)",
        "exp_kd_molar": 2.0e-9,
        "exp_kd_str": "2.0 nM",
        "exp_delta_g": -11.87,
        "citation": "Lee et al., Cell Death Differ. 15:888 (2008)"
    },
    "6QG8": {
        "target": "Bcl-2 Oncoprotein",
        "pep_name": "PUMA BH3 (20-mer)",
        "exp_kd_molar": 2.4e-9,
        "exp_kd_str": "2.4 nM",
        "exp_delta_g": -11.75,
        "citation": "Stewart et al., Nat. Chem. Biol. 15:1020 (2019)"
    },
    "1QZ7": {
        "target": "Beta-Catenin",
        "pep_name": "Axin (17-mer)",
        "exp_kd_molar": 7.2e-8,
        "exp_kd_str": "72 nM",
        "exp_delta_g": -9.75,
        "citation": "Xing et al., Structure 11:1391 (2003)"
    },
    "1TW6": {
        "target": "cIAP1 BIR3 Domain",
        "pep_name": "Smac (4-mer AVPI)",
        "exp_kd_molar": 2.8e-7,
        "exp_kd_str": "280 nM",
        "exp_delta_g": -8.94,
        "citation": "Sun et al., J. Med. Chem. 47:4149 (2004)"
    },
    "1SFI": {
        "target": "Bovine Trypsin",
        "pep_name": "SFTI-1 (14-mer Cyclic)",
        "exp_kd_molar": 1.0e-10,
        "exp_kd_str": "0.1 nM",
        "exp_delta_g": -13.64,
        "citation": "Korsinczky et al., JMB 311:579 (2001)"
    },
    "3D9L": {
        "target": "Nrd1 CTD Domain",
        "pep_name": "RNAPII CTD (7-mer)",
        "exp_kd_molar": 2.0e-5,
        "exp_kd_str": "20.0 µM",
        "exp_delta_g": -6.41,
        "citation": "Vassal et al., Mol. Cell 31:479 (2008)"
    }
}


def run_external_benchmark(
    benchmarks_dir: Path = Path("benchmarks/external_test_set"),
    out_csv: Path = Path("benchmarks/external_test_set/external_benchmark_results.csv")
) -> pd.DataFrame:
    benchmarks_dir = Path(benchmarks_dir).resolve()
    ml_model_file = Path("data/trained_complex_ensemble.joblib").resolve()
    ml_scorer = MLPeptideScorer(model_path=ml_model_file) if ml_model_file.exists() else None

    results = []

    print("\n" + "=" * 115)
    print("      AMP-SCREEN: FULLY EXTERNAL BENCHMARK EVALUATION (0% OVERLAP WITH 381 TRAINING SET)")
    print("=" * 115)
    print(f"Benchmarking directory: {benchmarks_dir} ({len(EXTERNAL_BENCHMARKS)} external complexes)\n")

    for pdb_id, meta in EXTERNAL_BENCHMARKS.items():
        pdb_path = benchmarks_dir / f"{pdb_id.lower()}.pdb"
        if not pdb_path.exists():
            logger.warning(f"File not found: {pdb_path}")
            continue

        st = gemmi.read_structure(str(pdb_path))
        seq = "".join(
            AA_3TO1.get(res.name.upper(), "X")
            for c in st[0] if c.name == "B"
            for res in c if res.name.upper() in AA_3TO1
        )

        phys = calculate_haddock_physics_score(pdb_path, receptor_chain="A", peptide_chain="B")
        contacts = analyze_cocomplex_contacts(pdb_path, receptor_chain="A", peptide_chain="B")

        z = calculate_net_charge(seq) if seq else 0.0
        mu_h = calculate_hydrophobic_moment(seq) if seq else 0.0
        boman = calculate_boman_index(seq) if seq else 0.0

        grid_pot = phys["haddock_score"] * 0.28
        sol_aff = estimate_solution_affinity(grid_pot, sequence=seq, conformation="auto")

        ml_delta_g = None
        ml_kd = "N/A"
        ml_pkd = None

        if ml_scorer and ml_scorer.is_calibrated:
            ml_res = ml_scorer.score_complex(
                pdb_path,
                hit_meta={
                    "pep_id": pdb_id,
                    "sequence": seq,
                    "length": len(seq),
                    "affinity": sol_aff["solution_delta_g"],
                    "docking_energy": sol_aff["solution_delta_g"],
                    "hbonds": contacts["hbonds"],
                    "salt_bridges": contacts["salt_bridges"],
                    "hydrophobic_contacts": contacts["hydrophobic_contacts"],
                    "clashes": contacts["clashes"],
                    "net_charge_z": z,
                    "hydrophobic_moment_muH": mu_h,
                    "boman_index": boman,
                }
            )
            ml_delta_g = ml_res["ml_delta_g"]
            ml_kd = ml_res["ml_kd"]
            ml_pkd = ml_res["ml_pkd"]

        results.append({
            "pdb_id": pdb_id,
            "target": meta["target"],
            "peptide_name": meta["pep_name"],
            "sequence": seq,
            "length": len(seq),
            "haddock_score": round(phys["haddock_score"], 1),
            "bsa_angstrom2": round(phys["bsa"], 1),
            "hbonds": contacts["hbonds"],
            "salt_bridges": contacts["salt_bridges"],
            "hydrophobic_contacts": contacts["hydrophobic_contacts"],
            "phys_delta_g": round(sol_aff["solution_delta_g"], 2),
            "phys_kd": sol_aff["solution_kd"],
            "ml_delta_g": round(ml_delta_g, 2) if ml_delta_g is not None else None,
            "ml_kd": ml_kd,
            "ml_pkd": round(ml_pkd, 2) if ml_pkd is not None else None,
            "exp_delta_g": meta["exp_delta_g"],
            "exp_kd": meta["exp_kd_str"],
            "citation": meta["citation"]
        })

    df = pd.DataFrame(results)
    df["err_phys"] = np.abs(df["phys_delta_g"] - df["exp_delta_g"])
    df["err_ml"] = np.abs(df["ml_delta_g"] - df["exp_delta_g"])

    # Pretty print summary table
    fmt = "%-6s | %-20s | %-16s | %-8s | %-8s | %-8s | %-8s | %-10s | %-10s"
    print(fmt % ("PDB", "Target Receptor", "Peptide (L)", "HADDOCK", "Exp dG", "Phys dG", "ML dG", "Exp Kd", "ML Kd"))
    print("-" * 115)
    for _, r in df.iterrows():
        pep_str = f"{r['sequence'][:8]}..({r['length']})" if r['length'] > 8 else f"{r['sequence']} ({r['length']})"
        print(fmt % (
            r["pdb_id"], r["target"][:20], pep_str,
            f"{r['haddock_score']:.1f}", f"{r['exp_delta_g']:.2f}",
            f"{r['phys_delta_g']:.2f}", f"{r['ml_delta_g']:.2f}",
            r["exp_kd"], r["ml_kd"]
        ))
    print("=" * 115)

    mae_phys = float(np.mean(df["err_phys"]))
    rmse_phys = float(np.sqrt(np.mean(df["err_phys"] ** 2)))
    r_phys, p_phys = pearsonr(df["phys_delta_g"], df["exp_delta_g"])
    rho_phys, _ = spearmanr(df["phys_delta_g"], df["exp_delta_g"])

    mae_ml = float(np.mean(df["err_ml"]))
    rmse_ml = float(np.sqrt(np.mean(df["err_ml"] ** 2)))
    r_ml, p_ml = pearsonr(df["ml_delta_g"], df["exp_delta_g"])
    rho_ml, _ = spearmanr(df["ml_delta_g"], df["exp_delta_g"])

    print("\nDescriptive Performance on Completely Unseen External Test Set (N=8):")
    print(f"  • Physical Force-Field: MAE = {mae_phys:.2f} kcal/mol, RMSE = {rmse_phys:.2f} kcal/mol, r = {r_phys:.3f} (p = {p_phys:.4f}), ρ = {rho_phys:.3f}")
    print(f"  • 56D ML Ensemble:      MAE = {mae_ml:.2f} kcal/mol, RMSE = {rmse_ml:.2f} kcal/mol, r = {r_ml:.3f} (p = {p_ml:.4f}), ρ = {rho_ml:.3f}")
    print("=" * 115 + "\n")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[OK] External benchmark results saved to: {out_csv}\n")
    return df


if __name__ == "__main__":
    run_external_benchmark()
