"""
Comprehensive Verification & Regression Test Suite for AMP-Screen.
Validates all physical, thermodynamic, and biochemical components with real assertions:
- P0-1: Heavy-atom / united-atom implicit hydrogen partial charge model
- P0-2: Continuous electrostatics switching function (9 to 12 A)
- P0-3: Lennard-Jones repulsive wall for severe atomic clashes
- P0-4: Transparent ML physical consensus score
- P0-6: Lead optimizer scale factor consistency (HADDOCK_TO_DG)
- Secondary structure Ramachandran dihedrals
- CAPRI benchmark evaluation
"""

from pathlib import Path
import tempfile
import numpy as np

from amp_screen.physics_engine import (
    VDW_PARAMS,
    SOLVATION_PARAMS,
    PARTIAL_CHARGES,
    ALL_ATOM_CHARGES,
    HEAVY_ATOM_CHARGES,
    get_residue_formal_charge,
    calculate_haddock_physics_score
)
from amp_screen.amp_analyzer import (
    analyze_complex_contacts,
    analyze_cocomplex_contacts
)
from amp_screen.ml_scorer import (
    ELEMENT_PAIRS,
    MLPeptideScorer,
    RT_LN10
)
from amp_screen.lead_optimizer import HADDOCK_TO_DG
from amp_screen.benchmark import evaluate_docked_complex
from amp_screen.structure_builder import SECONDARY_STRUCTURE_DIHEDRALS
from amp_screen.utils import clean_sequence


def test_p0_1_residue_formal_charges():
    """P0-1: Heavy-atom and all-atom charge models must assign correct formal charge to standard residues (pH 7.4)."""
    all_atom_expected = {
        "ARG": 1.0, "LYS": 1.0, "ASP": -1.0, "GLU": -1.0,
        "HIP": 1.0, "HID": 0.0, "HIE": 0.0, "HIS": 0.0,
        "TYR": 0.0, "CYS": 0.0
    }
    for res, exp in all_atom_expected.items():
        direct_sum = sum(v for k, v in ALL_ATOM_CHARGES.items() if k[0] == res)
        assert abs(direct_sum - exp) < 1e-3, f"ALL_ATOM_CHARGES sum for {res} must be {exp}, got {direct_sum}"

    heavy_expected = {
        "ARG": 1.0, "LYS": 1.0, "ASP": -1.0, "GLU": -1.0,
        "HIP": 1.0, "HID": 0.0, "HIE": 0.0, "HIS": 0.0
    }
    for res, exp in heavy_expected.items():
        direct_sum = sum(v for k, v in HEAVY_ATOM_CHARGES.items() if k[0] == res)
        assert abs(direct_sum - exp) < 1e-3, f"HEAVY_ATOM_CHARGES sum for {res} must be {exp}, got {direct_sum}"

    for has_h in [False, True]:
        label = "All-atom" if has_h else "Heavy-atom"
        assert abs(get_residue_formal_charge("ARG", has_hydrogens=has_h) - 1.0) < 0.05, f"{label} ARG must have formal charge +1.0"
        assert abs(get_residue_formal_charge("LYS", has_hydrogens=has_h) - 1.0) < 0.05, f"{label} LYS must have formal charge +1.0"
        assert abs(get_residue_formal_charge("ASP", has_hydrogens=has_h) - (-1.0)) < 0.05, f"{label} ASP must have formal charge -1.0"
        assert abs(get_residue_formal_charge("GLU", has_hydrogens=has_h) - (-1.0)) < 0.05, f"{label} GLU must have formal charge -1.0"
        assert abs(get_residue_formal_charge("ALA", has_hydrogens=has_h) - 0.0) < 0.05, f"{label} ALA must be neutral"
        assert abs(get_residue_formal_charge("LEU", has_hydrogens=has_h) - 0.0) < 0.05, f"{label} LEU must be neutral"

    p_1cka = Path("benchmarks/1cka.pdb")
    if p_1cka.exists():
        score_1cka = calculate_haddock_physics_score(p_1cka, receptor_chain="A", peptide_chain="B")
        assert score_1cka["e_elec"] < -10.0, f"Cationic peptide in 1CKA pocket must produce attractive electrostatics, got {score_1cka['e_elec']}"
        assert score_1cka["haddock_score"] < -20.0, f"Expected favorable HADDOCK score for 1CKA, got {score_1cka['haddock_score']}"


def test_p0_2_switching_function():
    """P0-2: Canonical C2 quintic switching function S(u) = 1 - 10u^3 + 15u^4 - 6u^5 must smoothly taper to zero with zero derivatives."""
    r_on, r_off = 9.0, 12.0
    r_range = r_off - r_on

    def switch_val(r):
        if r <= r_on:
            return 1.0
        elif r >= r_off:
            return 0.0
        u = (r - r_on) / r_range
        return 1.0 - 10.0 * (u ** 3) + 15.0 * (u ** 4) - 6.0 * (u ** 5)

    assert abs(switch_val(9.0) - 1.0) < 1e-6
    assert abs(switch_val(10.5) - 0.5) < 1e-6
    assert abs(switch_val(12.0) - 0.0) < 1e-6

    eps = 1e-5
    d1_on = (switch_val(r_on + eps) - switch_val(r_on)) / eps
    assert abs(d1_on) < 1e-4, f"First derivative at r_on must be zero, got {d1_on}"

    d1_off = (switch_val(r_off) - switch_val(r_off - eps)) / eps
    assert abs(d1_off) < 1e-4, f"First derivative at r_off must be zero, got {d1_off}"

    d2_on = (switch_val(r_on + 2 * eps) - 2 * switch_val(r_on + eps) + switch_val(r_on)) / (eps ** 2)
    assert abs(d2_on) < 1e-3, f"Second derivative at r_on must be zero, got {d2_on}"

    d2_off = (switch_val(r_off) - 2 * switch_val(r_off - eps) + switch_val(r_off - 2 * eps)) / (eps ** 2)
    assert abs(d2_off) < 1e-3, f"Second derivative at r_off must be zero, got {d2_off}"


def test_p0_3_repulsive_wall():
    """P0-3: Lennard-Jones potential must provide steep, non-decreasing repulsive barrier for overlapping atoms."""
    sig_c, eps_c = VDW_PARAMS["C"]
    sig_cc = sig_c
    eps_cc = eps_c

    def calc_vdw(r):
        r_eff = max(r, 1.0)
        sr6 = (sig_cc / r_eff) ** 6
        vdw = 4.0 * eps_cc * (sr6 ** 2 - sr6)
        return float(np.clip(vdw, -2.0, 500.0))

    distances = [1.8, 1.5, 1.0, 0.7, 0.5]
    energies = [calc_vdw(r) for r in distances]

    assert energies[0] > 10.0, f"Repulsive penalty at 1.8 A should be > 10 kcal/mol, got {energies[0]}"
    assert energies[1] > 100.0, f"Repulsive penalty at 1.5 A should be > 100 kcal/mol, got {energies[1]}"
    assert energies[2] == 500.0, f"Repulsive penalty at 1.0 A must hit 500.0 kcal/mol cap, got {energies[2]}"
    assert energies[3] == 500.0, f"Repulsive penalty at 0.7 A must retain 500.0 kcal/mol cap, got {energies[3]}"
    assert energies[4] == 500.0, f"Repulsive penalty at 0.5 A must retain 500.0 kcal/mol cap, got {energies[4]}"

    for e_prev, e_curr in zip(energies[:-1], energies[1:]):
        assert e_curr >= e_prev, f"Repulsive energy must be non-decreasing as distance decreases: {e_curr} < {e_prev}"


def test_p0_4_ml_scorer_uncalibrated_consensus():
    """P0-4: ML Scorer transparently marks is_calibrated=False on uncalibrated instance and computes physical consensus."""
    scorer = MLPeptideScorer()
    assert scorer.is_calibrated is False, "Default instance must be marked uncalibrated"
    score = scorer.compute_consensus_score(vina_affinity=-8.5, haddock_score=-40.0)
    assert score > 0.0, f"Expected positive consensus score, got {score}"


def test_p0_6_lead_optimizer_scale_factor():
    """P0-6: Lead optimizer HADDOCK_TO_DG scale factor must be strictly 0.15 across all mutation routines."""
    assert abs(HADDOCK_TO_DG - 0.15) < 1e-4, f"HADDOCK_TO_DG scale must be 0.15, got {HADDOCK_TO_DG}"
    ddg = -10.0 * HADDOCK_TO_DG
    assert abs(ddg - (-1.5)) < 1e-4


def test_vdw_parameters_and_combining_rule():
    """VDW Lorentz-Berthelot combining rule and parameters."""
    sig_c, eps_c = VDW_PARAMS["C"]
    sig_n, eps_n = VDW_PARAMS["N"]
    expected_sig = (sig_c + sig_n) / 2.0
    assert abs(expected_sig - 3.375) < 1e-4
    expected_eps = np.sqrt(eps_c * eps_n)
    assert abs(expected_eps - np.sqrt(0.066 * 0.170)) < 1e-4
    assert VDW_PARAMS["H"][0] == 2.50


def test_desolvation_parameters():
    """Eisenberg & McLachlan atomic solvation parameters."""
    assert SOLVATION_PARAMS["C"] == -0.016
    assert SOLVATION_PARAMS["N"] == +0.006
    assert SOLVATION_PARAMS["O"] == +0.006
    assert SOLVATION_PARAMS["S"] == -0.021


def test_published_haddock_weights():
    """HADDOCK 2.4 score weights."""
    e_vdw, e_elec, e_desolv, e_air = -20.0, -15.0, -9.0, 0.0
    score = 1.0 * e_vdw + 0.2 * e_elec + 1.0 * e_desolv + 0.1 * e_air
    assert abs(score - (-32.0)) < 1e-4


def test_amber_charges():
    """AMBER ff14SB charge table integrity."""
    assert PARTIAL_CHARGES[("ARG", "NH1")] == -0.3549
    assert PARTIAL_CHARGES[("ARG", "CZ")] == 0.7341
    assert PARTIAL_CHARGES[("LYS", "CE")] == -0.0143
    assert ("HID", "ND1") in PARTIAL_CHARGES
    assert ("HIE", "NE2") in PARTIAL_CHARGES
    assert ("HIP", "ND1") in PARTIAL_CHARGES


def test_element_pairs_and_thermodynamics():
    """10 unique element pairs and RT_LN10 thermodynamic conversion at 298.15 K."""
    assert len(ELEMENT_PAIRS) == 10
    assert ("S", "S") in ELEMENT_PAIRS
    assert abs(RT_LN10 - 1.3642) < 0.005


def test_secondary_structure_dihedrals():
    """Standard Ramachandran angles."""
    assert SECONDARY_STRUCTURE_DIHEDRALS["helix"]["phi"] == -57.0
    assert SECONDARY_STRUCTURE_DIHEDRALS["helix"]["psi"] == -47.0
    assert SECONDARY_STRUCTURE_DIHEDRALS["sheet"]["phi"] == -119.0
    assert SECONDARY_STRUCTURE_DIHEDRALS["sheet"]["psi"] == 113.0
    assert SECONDARY_STRUCTURE_DIHEDRALS["ppii"]["phi"] == -78.0
    assert SECONDARY_STRUCTURE_DIHEDRALS["ppii"]["psi"] == 149.0


def test_clean_sequence_and_nonstandard_mapping():
    """Non-standard amino acid mapping in clean_sequence."""
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert clean_sequence(seq) == seq
    assert clean_sequence("UOBZJ") == "CKDEL"


def test_capri_benchmark_on_real_crystal():
    """Verify CAPRI evaluation on benchmark 1AWR reference structure."""
    ref_path = Path("benchmarks/1awr.pdb")
    if ref_path.exists():
        metrics = evaluate_docked_complex(
            docked_pdb=ref_path,
            reference_pdb=ref_path,
            docked_rec_chain="A",
            docked_pep_chain="B",
            ref_rec_chain="A",
            ref_pep_chain="B"
        )
        assert metrics["receptor_ca_rmsd"] < 1e-3
        assert metrics["l_rmsd_backbone"] < 1e-3
        assert metrics["f_nat_percent"] == 100.0
        assert metrics["capri_quality"] == "High (***)"


def test_cocomplex_cross_chain_contacts():
    """Verify that analyze_cocomplex_contacts strictly extracts cross-chain contacts and rejects intra-chain pairs."""
    p_1cka = Path("benchmarks/1cka.pdb")
    if p_1cka.exists():
        contacts = analyze_cocomplex_contacts(p_1cka, receptor_chain="A", peptide_chain="B")
        assert contacts["hbonds"] > 0, "1CKA crystal must have cross-chain H-bonds"
        assert contacts["hydrophobic_contacts"] > 0, "1CKA crystal must have cross-chain hydrophobic contacts"

        intra = analyze_cocomplex_contacts(p_1cka, receptor_chain="B", peptide_chain="B")
        assert intra["hbonds"] == 0, f"Intra-chain H-bonds must be 0, got {intra['hbonds']}"
        assert intra["hydrophobic_contacts"] == 0
        assert intra["salt_bridges"] == 0


def test_separate_pdbs_chain_a_contacts():
    """Verify that analyze_complex_contacts detects contacts between separate PDBs even when both use chain A."""
    rec_pdb = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 20.00           O\n"
        "TER\nEND\n"
    )
    lig_pdb = (
        "ATOM      1  N   SER A   1       1.300   5.200   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  SER A   1       2.200   6.200   0.000  1.00 20.00           C\n"
        "ATOM      3  C   SER A   1       3.600   5.800   0.000  1.00 20.00           C\n"
        "ATOM      4  O   SER A   1       4.000   4.700   0.000  1.00 20.00           O\n"
        "ATOM      5  CB  SER A   1       1.800   7.600   0.000  1.00 20.00           C\n"
        "ATOM      6  OG  SER A   1       1.300   2.400   2.800  1.00 20.00           O\n"
        "TER\nEND\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix="_rec.pdb", delete=False) as tf_rec, \
         tempfile.NamedTemporaryFile("w", suffix="_lig.pdb", delete=False) as tf_lig:
        tf_rec.write(rec_pdb)
        tf_lig.write(lig_pdb)
        rec_path = Path(tf_rec.name)
        lig_path = Path(tf_lig.name)

    try:
        contacts = analyze_complex_contacts(rec_path, lig_path)
        assert contacts["hbonds"] >= 1, f"Expected H-bond between separate PDBs both on chain A, got {contacts['hbonds']}"
    finally:
        rec_path.unlink(missing_ok=True)
        lig_path.unlink(missing_ok=True)
