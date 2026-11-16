"""
Automated Comprehensive Integration & Feature Test for AMP-Screen.
Validates:
1. Auto-Pocket geometric cavity detector
2. Disulfide bond and macrocycle locking
3. Multi-conformation ensemble generation
4. End-to-end Screening funnel (HTVS + Refinement + HADDOCK/ML rescoring)
"""

import sys
import subprocess
from pathlib import Path
import tempfile
import shutil
import pytest
import pandas as pd

from amp_screen.receptor_prep import detect_surface_cavities, calculate_grid_box
from amp_screen.peptide_prep import detect_and_annotate_disulfides, prepare_peptide_pdbqt
from amp_screen.structure_builder import batch_fold_peptides

ROOT = Path(__file__).resolve().parent.parent
PYTHON = sys.executable

OBABEL = shutil.which("obabel") or (str(Path(sys.prefix) / "bin" / "obabel") if (Path(sys.prefix) / "bin" / "obabel").exists() else None)
HAS_OBABEL = OBABEL is not None

try:
    import meeko
    HAS_MEEKO = True
except ImportError:
    HAS_MEEKO = False


def test_auto_pocket():
    target_pdb = ROOT / "examples" / "receptors" / "target.pdb"
    if not target_pdb.exists():
        pytest.skip("target.pdb not found in examples/receptors")
    pockets = detect_surface_cavities(target_pdb, max_pockets=3)
    assert len(pockets) >= 1, "Failed to detect surface cavities in target.pdb"

    grid = calculate_grid_box(target_pdb, auto_pocket=True)
    assert "center_x" in grid and "detected_pocket_residues" in grid


def test_ensemble_folding():
    test_fa = {"MiniAMP": "RRWWRF"}
    with tempfile.TemporaryDirectory() as tmp_dir:
        folded = batch_fold_peptides(test_fa, Path(tmp_dir), conformation="ensemble")
        assert len(folded) == 3, f"Expected 3 ensemble conformers, got {len(folded)}"


@pytest.mark.skipif(not (HAS_OBABEL and HAS_MEEKO), reason="OpenBabel binary or Meeko library not available")
def test_disulfide_locking():
    pdb_text = """ATOM      1  N   CYS A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  CYS A   1       1.450   0.000   0.000  1.00 20.00           C
ATOM      3  C   CYS A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   CYS A   1       1.300   2.400   0.000  1.00 20.00           O
ATOM      5  CB  CYS A   1       2.000  -0.750   1.200  1.00 20.00           C
ATOM      6  SG  CYS A   1       1.500  -2.500   1.200  1.00 20.00           S
ATOM      7  N   CYS A   2       3.300   1.400   0.000  1.00 20.00           N
ATOM      8  CA  CYS A   2       4.000   2.700   0.000  1.00 20.00           C
ATOM      9  C   CYS A   2       5.500   2.500   0.000  1.00 20.00           C
ATOM     10  O   CYS A   2       6.200   3.500   0.000  1.00 20.00           O
ATOM     11  CB  CYS A   2       3.600   3.500   1.200  1.00 20.00           C
ATOM     12  SG  CYS A   2       1.500  -0.500   1.200  1.00 20.00           S
END
"""
    with tempfile.TemporaryDirectory() as tmp_dir:
        t_dir = Path(tmp_dir)
        test_pdb = t_dir / "test_cyclic.pdb"
        test_ann = t_dir / "test_cyclic.annotated.pdb"
        test_pdbqt = t_dir / "test_cyclic.pdbqt"

        test_pdb.write_text(pdb_text)
        disulfides, _ = detect_and_annotate_disulfides(test_pdb, test_ann)
        assert len(disulfides) == 1, "Failed to perceive CYS1-CYS2 disulfide bond!"

        out_pdbqt, torsions = prepare_peptide_pdbqt(
            test_pdb,
            test_pdbqt,
            flexibility="sidechain",
            obabel_bin=OBABEL,
            lock_disulfides=True
        )
        assert out_pdbqt.exists(), "PDBQT generation failed for cyclic peptide"


@pytest.mark.end_to_end
def test_end_to_end_native_screening():
    target_pdb = ROOT / "benchmarks" / "1dkx.pdb"
    peptides_fa = ROOT / "examples" / "peptides" / "benchmark_amps.fasta"

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir) / "screening_out"
        cmd = [
            PYTHON, "-m", "amp_screen.cli", "run",
            "--receptor", str(target_pdb),
            "--library", str(peptides_fa),
            "--output", str(out_dir),
            "--auto-pocket",
            "--flexibility", "rigid",
            "--max-workers", "2",
            "--exhaustiveness", "2",
            "--refine-top", "2",
            "--refine-exhaustiveness", "2",
            "--top", "2",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        assert res.returncode == 0, f"Screening pipeline failed: {res.stderr}"

        csv_report = out_dir / "screening_leaderboard.csv"
        assert csv_report.exists(), "Missing screening_leaderboard.csv"
        df = pd.read_csv(csv_report)
        assert "haddock_score" in df.columns
        assert "ml_pkd" in df.columns
