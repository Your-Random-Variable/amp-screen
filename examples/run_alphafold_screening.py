#!/usr/bin/env python3
"""
AlphaFold-Driven Virtual Screening Example.
===========================================
End-to-end AlphaFold workflow for target-based peptide screening:
1. Target Retrieval: Fetches high-confidence predicted structure from AlphaFold Database (AFDB).
2. Input Preparation: Prepares AlphaFold-Multimer & AlphaFold 3 co-folding inputs.
3. Complex Scoring: Evaluates co-folded complexes for interface pLDDT, HADDOCK physics,
   calibrated ML binding affinity (Kd, ΔG), and biophysical safety profiling.

Usage:
    python examples/run_alphafold_screening.py
"""

from pathlib import Path
import pandas as pd

from amp_screen.alphafold import (
    fetch_alphafold_target,
    prepare_alphafold_multimer_inputs,
    score_alphafold_complex
)
from amp_screen.utils import ensure_dir

# 1. Target Definition (Specified by UniProt accession)
TARGET_UNIPROT = "P04637"  # Cellular tumor antigen p53 (or any target accession)
TARGET_NAME = "p53"

# 2. Candidate Peptides under evaluation
CANDIDATE_PEPTIDES = {
    "Lead_1": {
        "sequence": "ETFSDLWKLLPEN",
        "description": "Wild-type interaction motif"
    },
    "Lead_2": {
        "sequence": "LTFEHWWAQLTS",
        "description": "High-affinity optimized variant"
    },
    "Lead_3": {
        "sequence": "TSFAEYWNLLSP",
        "description": "Stapled peptide mimetic sequence"
    }
}


def main():
    base_out = Path("alphafold_screening_workspace")
    ensure_dir(base_out)

    print("\n" + "=" * 100)
    print("      ALPHAFOLD-DRIVEN VIRTUAL SCREENING WORKFLOW")
    print(f"      Target: {TARGET_NAME} (UniProt: {TARGET_UNIPROT}) | AlphaFold Database + AF-Multimer Suite")
    print("=" * 100 + "\n")

    # Step 1: Fetch target receptor structure from AlphaFold Database (AFDB)
    target_out_pdb = base_out / f"{TARGET_NAME}_AF.pdb"
    print(f"[Step 1] Fetching receptor structure from AlphaFold Database for UniProt: {TARGET_UNIPROT}...")
    target_info = fetch_alphafold_target(
        uniprot_id=TARGET_UNIPROT,
        out_pdb=target_out_pdb,
        clean=True
    )
    print(f"  ✓ Organism          : {target_info['organism']}")
    print(f"  ✓ Sequence Length   : {target_info['length']} residues")
    print(f"  ✓ Global pLDDT      : {target_info['global_plddt']:.1f} / 100")
    print(f"  ✓ High-Confidence   : {target_info['confident_residues_pct']:.1f}% residues (pLDDT >= 70)")
    print(f"  ✓ Receptor PDB      : {target_info['pdb_path']}\n")

    # Step 2: Prepare AlphaFold-Multimer & AlphaFold 3 Co-Folding Inputs
    print("[Step 2] Preparing AlphaFold-Multimer & AlphaFold 3 co-folding inputs for candidates...")
    af_inputs_dir = base_out / "af_multimer_inputs"
    ensure_dir(af_inputs_dir)

    for pep_id, info in CANDIDATE_PEPTIDES.items():
        pep_seq = info["sequence"]
        inputs = prepare_alphafold_multimer_inputs(
            target_seq=target_info["sequence"],
            peptide_seq=pep_seq,
            job_name=f"{TARGET_NAME}_{pep_id}",
            out_dir=af_inputs_dir / pep_id,
            target_name=TARGET_NAME,
            peptide_name=pep_id
        )
        print(f"  ✓ [{pep_id:<8}] FASTA: {inputs['colabfold_fasta'].name} | AF3 JSON: {inputs['af3_json'].name}")

    print("\n  [Tip] Inputs can be submitted to:")
    print("        * AlphaFold 3 Server: https://alphafoldserver.com (upload the _af3.json)")
    print("        * Google ColabFold: upload the _colabfold.fasta\n")

    # Step 3: Example of scoring a co-folded or benchmark complex PDB
    benchmark_complex = Path("benchmarks/1dkx.pdb")
    if benchmark_complex.exists():
        print(f"[Step 3] Scoring reference complex structure ({benchmark_complex.name})...")
        score_res = score_alphafold_complex(
            complex_pdb=benchmark_complex,
            target_chain="A",
            peptide_chain="B"
        )
        print(f"  ✓ Decision          : {score_res['decision']}")
        print(f"  ✓ Rating            : {score_res['star_rating']}")
        print(f"  ✓ Net Charge (pH 7) : {score_res['net_charge']:+.1f}")
        print(f"  ✓ Predicted ML Kd   : {score_res['ml_kd']}")
        print(f"  ✓ Predicted ML ΔG   : {score_res['ml_delta_g']:.2f} kcal/mol")
        print(f"  ✓ Buried SASA (BSA) : {score_res['bsa_angstrom2']:.1f} Å²")
        print(f"  ✓ Hemolysis Risk    : {score_res['hemolysis_risk']}\n")

    print("[Done] AlphaFold workflow completed successfully.\n")


if __name__ == "__main__":
    main()
