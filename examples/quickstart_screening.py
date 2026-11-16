#!/usr/bin/env python3
"""
Quickstart Virtual Screening Example.
=====================================
Demonstrates fast, high-throughput protein-peptide virtual screening:
1. Takes a target receptor structure (PDB) and a candidate peptide library (FASTA).
2. Performs automated surface pocket detection.
3. Rapidly screens peptides using native 3D geometric surface docking (2-5 s / peptide).
4. Predicts solution binding affinity (Kd, ΔG) with calibrated 96D/56D ML models.
5. Generates ranked leaderboard and docked co-complex PDB structures.

Usage:
    python examples/quickstart_screening.py
"""

from pathlib import Path
import subprocess
import sys


def main():
    repo_root = Path(__file__).resolve().parent.parent
    receptor_pdb = repo_root / "benchmarks" / "1dkx.pdb"
    peptide_fasta = repo_root / "examples" / "peptides" / "benchmark_amps.fasta"
    output_dir = repo_root / "screening_quickstart_results"

    print("\n" + "=" * 90)
    print("                AMP-SCREEN QUICKSTART VIRTUAL SCREENING")
    print("=" * 90)
    print(f"  • Target Receptor : {receptor_pdb}")
    print(f"  • Peptide Library : {peptide_fasta}")
    print(f"  • Output Directory: {output_dir}")
    print("=" * 90 + "\n")

    cmd = [
        sys.executable, "-m", "amp_screen.cli", "run",
        "--receptor", str(receptor_pdb),
        "--rec-chain", "A",
        "--library", str(peptide_fasta),
        "--output", str(output_dir),
        "--auto-pocket",
        "--refine-top", "3",
        "--top", "3",
        "--threads", "4"
    ]

    print(f"Executing: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, check=True)

    leaderboard_csv = output_dir / "screening_leaderboard.csv"
    if leaderboard_csv.exists():
        print("\n" + "=" * 90)
        print("                        TOP SCREENING HITS")
        print("=" * 90)
        with open(leaderboard_csv, "r") as f:
            for i, line in enumerate(f):
                if i < 10:  # Print header + top hits
                    print(line.strip())
        print("=" * 90)
        print(f"\n[Success] Screen completed! Docked co-complexes saved to: {output_dir / 'top_complexes'}\n")


if __name__ == "__main__":
    main()
