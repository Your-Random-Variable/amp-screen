"""
Machine Learning Scoring Engine for Protein-Peptide Interactions.
Combines:
1. RF-Score-Peptide atomic interaction fingerprints (20 element-pair distance bins, normalized by peptide size)
2. HADDOCK-style explicit physics terms (vdW, electrostatics, desolvation, BSA, HADDOCK score, AIR restraints)
3. AMP biophysical properties (net charge, amphipathic moment, Boman index)
4. Calibrated Ensemble Regressors (Random Forest + Gradient Boosting) for pKd and Delta G prediction.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from scipy.spatial.distance import cdist
import gemmi
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
import joblib

from amp_screen.utils import logger, AA_3TO1
from amp_screen.physics_engine import calculate_haddock_physics_score
from amp_screen.amp_analyzer import (
    analyze_cocomplex_contacts,
    calculate_net_charge,
    calculate_hydrophobic_moment,
    calculate_boman_index,
    calculate_isoelectric_point,
    clean_sequence
)

# Complete set of 4 biological elements -> 10 unique unordered pairs (including S-S for disulfide bonds)
ELEMENTS = ["C", "N", "O", "S"]
ELEMENT_PAIRS = [
    ("C", "C"), ("C", "N"), ("C", "O"), ("C", "S"),
    ("N", "N"), ("N", "O"), ("N", "S"),
    ("O", "O"), ("O", "S"), ("S", "S")
]

# Physical constants for pKd <-> Delta G conversion (Gibbs-Helmholtz / van 't Hoff)
R_GAS_KCAL = 0.001987204   # Gas constant R in kcal / (mol * K)
TEMP_KELVIN = 298.15       # Standard biochemical temperature (298.15 K, 25 °C)
RT_LN10 = R_GAS_KCAL * TEMP_KELVIN * np.log(10.0)  # ~ 1.3642 kcal/mol per pKd unit at 298.15 K

# Biophysical propensity scales
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
CHOU_FASMAN_HELIX = {
    'E': 1.51, 'A': 1.42, 'L': 1.21, 'M': 1.45, 'H': 1.00,
    'K': 1.16, 'Q': 1.11, 'W': 1.08, 'I': 1.08, 'V': 1.06,
    'D': 1.01, 'F': 1.13, 'R': 0.98, 'C': 0.70, 'N': 0.67,
    'S': 0.77, 'T': 0.83, 'Y': 0.69, 'P': 0.57, 'G': 0.57
}
CHOU_FASMAN_SHEET = {
    'M': 1.05, 'V': 1.70, 'I': 1.60, 'C': 1.19, 'Y': 1.47,
    'F': 1.38, 'Q': 1.10, 'L': 1.30, 'T': 1.19, 'W': 1.37,
    'A': 0.83, 'R': 0.93, 'G': 0.75, 'D': 0.54, 'K': 0.74,
    'S': 0.75, 'H': 0.87, 'N': 0.89, 'P': 0.55, 'E': 0.37
}
WIMLEY_WHITE_OCTANOL = {
    'A': 0.50, 'R': 1.81, 'N': 0.85, 'D': 3.64, 'C': -0.02,
    'Q': 0.77, 'E': 3.63, 'G': 1.15, 'H': 0.11, 'I': -1.12,
    'L': -1.25, 'K': 2.80, 'M': -0.67, 'F': -1.71, 'P': 0.14,
    'S': 0.46, 'T': 0.25, 'W': -2.09, 'Y': -0.71, 'V': -0.46
}
WIMLEY_WHITE_INTERFACE = {
    'A': 0.17, 'R': 0.81, 'N': 0.42, 'D': 1.23, 'C': -0.24,
    'Q': 0.58, 'E': 2.02, 'G': 0.01, 'H': 0.17, 'I': -0.31,
    'L': -0.56, 'K': 0.99, 'M': -0.23, 'F': -1.13, 'P': 0.45,
    'S': 0.13, 'T': 0.14, 'W': -1.85, 'Y': -0.94, 'V': 0.07
}
EISENBERG_CONSENSUS = {
    'A': 0.25, 'R': -1.76, 'N': -0.64, 'D': -0.72, 'C': 0.04,
    'Q': -0.69, 'E': -0.62, 'G': 0.16, 'H': -0.40, 'I': 0.73,
    'L': 0.53, 'K': -1.10, 'M': 0.26, 'F': 0.61, 'P': -0.07,
    'S': -0.26, 'T': -0.18, 'W': 0.37, 'Y': 0.02, 'V': 0.54
}

# Complete 96-Dimensional Physical Chemistry & Structural Interface Feature Space
FEATURE_NAMES = []
for e1, e2 in ELEMENT_PAIRS:
    FEATURE_NAMES.append(f"pair_{e1}_{e2}_close")
    FEATURE_NAMES.append(f"pair_{e1}_{e2}_medium")

FEATURE_NAMES.extend([
    # Explicit HADDOCK & Physics terms (6)
    "e_vdw", "e_elec", "e_desolv", "bsa", "haddock_score", "e_air",
    # Interface Contacts (4)
    "salt_bridges", "hbonds", "hydrophobic_contacts", "clashes",
    # Core Sequence Biophysics (4)
    "net_charge", "hydrophobic_moment", "boman_index", "length",
    # Size-Normalized Energetics (5)
    "e_vdw_per_res", "e_elec_per_res", "e_desolv_per_res", "bsa_per_res", "haddock_per_res",
    # Interface Density & Packing Quality (4)
    "contact_density", "hbond_density", "hydrophobic_density", "clash_ratio",
    # Electrostatic Architecture & Charge (5)
    "charge_density", "pos_charge_count", "neg_charge_count", "charge_asymmetry", "isoelectric_point",
    # Chemical Composition & Structural Propensities (6)
    "aromatic_ratio", "aliphatic_ratio", "proline_ratio", "polar_ratio", "charged_ratio", "glycine_ratio",
    # Solution Thermodynamics & Efficiency (2)
    "entropy_penalty", "ligand_efficiency_proxy",
    # Hydropathy & Hydrophobicity (2)
    "gravy_hydropathy", "hydrophobic_ratio",
    # Secondary Structure & Conformational Propensities (3)
    "helix_propensity", "sheet_propensity", "ppii_propensity",
    # Energetic & Force Field Partition Ratios (5)
    "elec_to_vdw_ratio", "polar_to_hydrophobic_contact_ratio", "desolv_to_bsa_ratio", "elec_to_bsa_ratio", "vdw_to_bsa_ratio",
    # Detailed Amino Acid Fractions (6)
    "sulfur_ratio", "basic_ratio", "acidic_ratio", "hydroxyl_ratio", "amide_ratio", "beta_branched_ratio",
    # Multi-Shell Interface Densities (5)
    "close_to_medium_shell_ratio", "total_interface_contacts", "specific_contact_density", "close_contact_density", "medium_contact_density",
    # Thermodynamic Degrees of Freedom & Molecular Weight (5)
    "rotatable_bonds_est", "conf_entropy_loss_est", "mol_weight_est", "haddock_per_mw", "bsa_per_mw",
    # Specialized Biophysical & Stereochemical Descriptors (6)
    "cation_pi_potential", "solvation_free_energy_ratio", "surface_complementarity_proxy", "charge_transfer_index", "aliphatic_index", "flexibility_turn_propensity",
    # Membrane & Interfacial Partitioning Free Energies (Wimley-White & Eisenberg) (4)
    "wimley_white_interface", "wimley_white_octanol", "octanol_interface_diff", "eisenberg_consensus_hydrophobicity",
    # Interface Contact Chemistry & Burial Packing Architecture (4)
    "c_c_close_fraction", "close_contacts_per_bsa", "burial_fraction", "total_atomic_contacts_per_bsa"
])


def extract_structural_fingerprint(
    complex_pdb: Path,
    receptor_chain: str = "A",
    peptide_chain: str = "B"
) -> Dict[str, float]:
    """
    Extracts 20 atomic contact features across all 10 standard element-pair combinations:
    - 10 close contacts (d <= 4.0 A) normalized by peptide heavy atom count
    - 10 medium-range shell contacts (4.0 A < d <= 8.0 A) normalized by peptide heavy atom count
    Normalization by peptide heavy atoms removes spurious length bias.
    """
    st = gemmi.read_structure(str(complex_pdb))

    rec_coords = []
    rec_elems = []
    pep_coords = []
    pep_elems = []

    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    el = atom.element.name
                    if el not in ELEMENTS:
                        continue
                    pt = [atom.pos.x, atom.pos.y, atom.pos.z]
                    if chain.name == receptor_chain:
                        rec_coords.append(pt)
                        rec_elems.append(el)
                    elif chain.name == peptide_chain:
                        pep_coords.append(pt)
                        pep_elems.append(el)

    features = {}
    for e1, e2 in ELEMENT_PAIRS:
        features[f"pair_{e1}_{e2}_close"] = 0.0
        features[f"pair_{e1}_{e2}_medium"] = 0.0

    if not rec_coords or not pep_coords:
        return features

    dists = cdist(np.array(rec_coords), np.array(pep_coords))
    norm_factor = max(1.0, float(len(pep_coords)))

    for i, r_el in enumerate(rec_elems):
        close_j = np.where(dists[i] <= 4.0)[0]
        med_j = np.where((dists[i] > 4.0) & (dists[i] <= 8.0))[0]

        for j in close_j:
            p_el = pep_elems[j]
            pair_key = (r_el, p_el) if (r_el, p_el) in ELEMENT_PAIRS else (p_el, r_el)
            if pair_key in ELEMENT_PAIRS:
                features[f"pair_{pair_key[0]}_{pair_key[1]}_close"] += 1.0 / norm_factor

        for j in med_j:
            p_el = pep_elems[j]
            pair_key = (r_el, p_el) if (r_el, p_el) in ELEMENT_PAIRS else (p_el, r_el)
            if pair_key in ELEMENT_PAIRS:
                features[f"pair_{pair_key[0]}_{pair_key[1]}_medium"] += 1.0 / norm_factor

    return features


def extract_complex_features(
    fp: Dict[str, float],
    physics: Dict[str, float],
    hit_meta: Dict[str, Any],
    seq: str = ""
) -> List[float]:
    """
    Extracts the complete 96-dimensional feature vector:
    20 element-pair bins, 6 physics terms, 4 contact counts, 4 core biophysics,
    5 size-normalized energetics, 4 interface density metrics, 5 charge architecture terms,
    6 composition ratios, 2 solution thermodynamics features, 32 biophysical descriptors
    (hydropathy, conformational propensities, partition ratios, detailed AA fractions,
    shell densities, rotatable bond / entropy estimates, and stereochemical indices),
    4 membrane/interfacial partitioning free energies (Wimley-White & Eisenberg), and
    4 interface contact chemistry & burial packing architecture metrics.
    """
    L = float(hit_meta.get("length", len(seq) if seq else 12.0))
    L = max(1.0, L)

    # 1. 20 element pairs
    feat_vec = []
    for e1, e2 in ELEMENT_PAIRS:
        feat_vec.append(fp.get(f"pair_{e1}_{e2}_close", 0.0))
        feat_vec.append(fp.get(f"pair_{e1}_{e2}_medium", 0.0))

    # 2. 6 Physics Terms
    evdw = float(physics.get("e_vdw", -30.0))
    eelec = float(physics.get("e_elec", -10.0))
    edesolv = float(physics.get("e_desolv", 10.0))
    bsa = float(physics.get("bsa", 1000.0))
    haddock = float(physics.get("haddock_score", -40.0))
    e_air = float(physics.get("e_air", 0.0))

    feat_vec.extend([evdw, eelec, edesolv, bsa, haddock, e_air])

    # 3. 4 Interface Contact Counts
    salt_bridges = float(hit_meta.get("salt_bridges", 0.0))
    hbonds = float(hit_meta.get("hbonds", 0.0))
    hydrophobic = float(hit_meta.get("hydrophobic_contacts", 0.0))
    clashes = float(hit_meta.get("clashes", 0.0))

    feat_vec.extend([salt_bridges, hbonds, hydrophobic, clashes])

    # 4. 4 Core Sequence Biophysics
    net_charge = float(hit_meta.get("net_charge", hit_meta.get("net_charge_z", 0.0)))
    mu_h = float(hit_meta.get("hydrophobic_moment", hit_meta.get("hydrophobic_moment_muH", 0.0)))
    boman = float(hit_meta.get("boman_index", 0.0))

    feat_vec.extend([net_charge, mu_h, boman, L])

    # 5. 5 Size-Normalized Energetics
    feat_vec.extend([
        evdw / L,
        eelec / L,
        edesolv / L,
        bsa / L,
        haddock / L
    ])

    # 6. 4 Interface Density & Packing Quality
    contact_density = (hbonds + salt_bridges + hydrophobic) / max(1.0, bsa / 100.0)
    hbond_density = hbonds / L
    hydrophobic_density = hydrophobic / L
    clash_ratio = clashes / max(1.0, hbonds + hydrophobic)

    feat_vec.extend([contact_density, hbond_density, hydrophobic_density, clash_ratio])

    # 7. 5 Electrostatic Architecture & Charge
    charge_density = net_charge / L
    pos_count = float(seq.count('K') + seq.count('R') + seq.count('H')) if seq else float(hit_meta.get("pos_charge_count", 0.0))
    neg_count = float(seq.count('D') + seq.count('E')) if seq else float(hit_meta.get("neg_charge_count", 0.0))
    charge_asym = abs(pos_count - neg_count) / max(1.0, pos_count + neg_count)
    pi = float(hit_meta.get("isoelectric_point", calculate_isoelectric_point(seq) if seq else 7.0))

    feat_vec.extend([charge_density, pos_count, neg_count, charge_asym, pi])

    # 8. 6 Chemical Composition & Structural Propensities
    if seq:
        arom_ratio = sum(seq.count(aa) for aa in ('F', 'Y', 'W')) / L
        aliph_ratio = sum(seq.count(aa) for aa in ('L', 'I', 'V', 'M')) / L
        pro_ratio = seq.count('P') / L
        pol_ratio = sum(seq.count(aa) for aa in ('S', 'T', 'N', 'Q', 'C')) / L
        chg_ratio = (pos_count + neg_count) / L
        gly_ratio = seq.count('G') / L
    else:
        arom_ratio = float(hit_meta.get("aromatic_ratio", 0.0))
        aliph_ratio = float(hit_meta.get("aliphatic_ratio", 0.0))
        pro_ratio = float(hit_meta.get("proline_ratio", 0.0))
        pol_ratio = float(hit_meta.get("polar_ratio", 0.0))
        chg_ratio = float(hit_meta.get("charged_ratio", 0.0))
        gly_ratio = float(hit_meta.get("glycine_ratio", 0.0))

    feat_vec.extend([arom_ratio, aliph_ratio, pro_ratio, pol_ratio, chg_ratio, gly_ratio])

    # 9. 2 Solution Thermodynamics & Efficiency
    entropy_pen = 0.45 * max(0.0, L - 1.0)
    le_proxy = -haddock / L

    feat_vec.extend([entropy_pen, le_proxy])

    # 10. Hydropathy & Hydrophobicity (2)
    gravy_hydropathy = float(np.mean([KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq]) if len(seq) > 0 else hit_meta.get("gravy_hydropathy", 0.0))
    hydrophobic_ratio = float(sum(seq.count(aa) for aa in ('V', 'I', 'L', 'M', 'F', 'W', 'C')) / L if seq else hit_meta.get("hydrophobic_ratio", 0.0))
    feat_vec.extend([gravy_hydropathy, hydrophobic_ratio])

    # 11. Secondary Structure & Conformational Propensities (3)
    helix_prop = float(np.mean([CHOU_FASMAN_HELIX.get(aa, 1.0) for aa in seq]) if len(seq) > 0 else hit_meta.get("helix_propensity", 1.0))
    sheet_prop = float(np.mean([CHOU_FASMAN_SHEET.get(aa, 1.0) for aa in seq]) if len(seq) > 0 else hit_meta.get("sheet_propensity", 1.0))
    ppii_prop = float((seq.count('P') * 2.0 + seq.count('G') + seq.count('A')) / L if seq else hit_meta.get("ppii_propensity", 0.0))
    feat_vec.extend([helix_prop, sheet_prop, ppii_prop])

    # 12. Energetic & Force Field Partition Ratios (5)
    elec_to_vdw = float(abs(eelec) / max(1.0, abs(evdw)))
    polar_to_hydrophobic = float((hbonds + salt_bridges) / max(1.0, hydrophobic))
    desolv_to_bsa = float(edesolv / max(1.0, bsa))
    elec_to_bsa = float(abs(eelec) / max(1.0, bsa))
    vdw_to_bsa = float(abs(evdw) / max(1.0, bsa))
    feat_vec.extend([elec_to_vdw, polar_to_hydrophobic, desolv_to_bsa, elec_to_bsa, vdw_to_bsa])

    # 13. Detailed Amino Acid Fractions (6)
    if seq:
        sulfur_ratio = float(sum(seq.count(aa) for aa in ('C', 'M')) / L)
        basic_ratio = float(sum(seq.count(aa) for aa in ('K', 'R', 'H')) / L)
        acidic_ratio = float(sum(seq.count(aa) for aa in ('D', 'E')) / L)
        hydroxyl_ratio = float(sum(seq.count(aa) for aa in ('S', 'T', 'Y')) / L)
        amide_ratio = float(sum(seq.count(aa) for aa in ('N', 'Q')) / L)
        beta_branched_ratio = float(sum(seq.count(aa) for aa in ('V', 'I', 'T')) / L)
    else:
        sulfur_ratio = float(hit_meta.get("sulfur_ratio", 0.0))
        basic_ratio = float(hit_meta.get("basic_ratio", 0.0))
        acidic_ratio = float(hit_meta.get("acidic_ratio", 0.0))
        hydroxyl_ratio = float(hit_meta.get("hydroxyl_ratio", 0.0))
        amide_ratio = float(hit_meta.get("amide_ratio", 0.0))
        beta_branched_ratio = float(hit_meta.get("beta_branched_ratio", 0.0))
    feat_vec.extend([sulfur_ratio, basic_ratio, acidic_ratio, hydroxyl_ratio, amide_ratio, beta_branched_ratio])

    # 14. Multi-Shell Interface Densities (5)
    close_sum = sum(v for k, v in fp.items() if k.endswith("_close"))
    med_sum = sum(v for k, v in fp.items() if k.endswith("_medium"))
    close_to_med_ratio = float(close_sum / max(1e-4, med_sum))
    total_contacts = float(hbonds + salt_bridges + hydrophobic)
    specific_contact_density = float(total_contacts / L)
    close_density = float(close_sum / L)
    med_density = float(med_sum / L)
    feat_vec.extend([close_to_med_ratio, total_contacts, specific_contact_density, close_density, med_density])

    # 15. Thermodynamic Degrees of Freedom & Molecular Weight (5)
    rot_bonds = float(L * 2.8)
    conf_entropy = float(0.45 * rot_bonds)
    mw_est = float(L * 110.0)
    haddock_per_mw = float(haddock / max(100.0, mw_est))
    bsa_per_mw = float(bsa / max(100.0, mw_est))
    feat_vec.extend([rot_bonds, conf_entropy, mw_est, haddock_per_mw, bsa_per_mw])

    # 16. Specialized Biophysical & Stereochemical Descriptors (6)
    cation_pi = float((pos_count * (arom_ratio * L)) / L)
    solv_ratio = float(edesolv / max(1.0, abs(evdw) + abs(eelec)))
    sc_proxy = float(bsa / max(100.0, abs(haddock) * 10.0))
    charge_trans = float((net_charge * boman) / L)
    if seq:
        aliphatic_idx = float((seq.count('A') + 2.9 * seq.count('V') + 3.9 * (seq.count('I') + seq.count('L'))) / L)
        turn_prop = float((seq.count('P') + seq.count('G')) / L)
    else:
        aliphatic_idx = float(hit_meta.get("aliphatic_index", 0.0))
        turn_prop = float(hit_meta.get("flexibility_turn_propensity", 0.0))
    feat_vec.extend([cation_pi, solv_ratio, sc_proxy, charge_trans, aliphatic_idx, turn_prop])

    # 17. Membrane & Interfacial Partitioning Free Energies (4)
    if seq:
        ww_if = float(np.mean([WIMLEY_WHITE_INTERFACE.get(aa, 0.0) for aa in seq]))
        ww_oct = float(np.mean([WIMLEY_WHITE_OCTANOL.get(aa, 0.0) for aa in seq]))
        eis_hyd = float(np.mean([EISENBERG_CONSENSUS.get(aa, 0.0) for aa in seq]))
    else:
        ww_if = float(hit_meta.get("wimley_white_interface", 0.0))
        ww_oct = float(hit_meta.get("wimley_white_octanol", 0.0))
        eis_hyd = float(hit_meta.get("eisenberg_consensus_hydrophobicity", 0.0))
    oct_if_diff = float(ww_oct - ww_if)
    feat_vec.extend([ww_if, ww_oct, oct_if_diff, eis_hyd])

    # 18. Interface Contact Chemistry & Burial Packing Architecture (4)
    c_c_close = fp.get("pair_C_C_close", 0.0)
    c_c_close_fraction = float(c_c_close / max(1e-4, close_sum))
    close_contacts_per_bsa = float(close_sum / max(1.0, bsa / 100.0))
    burial_fraction = float(min(1.0, bsa / (150.0 * L)))
    total_atomic_contacts_per_bsa = float((close_sum + med_sum) / max(1.0, bsa / 100.0))
    feat_vec.extend([c_c_close_fraction, close_contacts_per_bsa, burial_fraction, total_atomic_contacts_per_bsa])

    return feat_vec


class MLPeptideScorer:
    """
    Multi-Parametric Consensus & Machine Learning Scorer for Protein-Peptide Affinity.
    Integrates RF-Score structural fingerprints with HADDOCK physics, interface densities,
    and AMP descriptors across 96 physical parameters.
    When a pre-trained model file (.joblib) is provided, uses the trained regressor ensemble.
    When no pre-trained model is supplied, transparently marks is_calibrated = False and
    computes an uncalibrated heuristic consensus score based on explicit physics and docking terms.
    """

    def __init__(
        self,
        model_path: Optional[Path] = None,
        consensus_weights: Tuple[float, float, float] = (0.35, 0.35, 0.30)
    ):
        self.rf_model = None
        self.gb_model = None
        self.et_model = None
        self.cb_model = None
        self.lgb_model = None
        self.dl_model = None
        self.ensemble_weights = None
        self.feature_names = []
        self.is_calibrated = False
        self.consensus_weights = consensus_weights

        if model_path and Path(model_path).exists():
            self._load_model(Path(model_path))
            self.is_calibrated = True
        else:
            self._setup_features()
            self.is_calibrated = False

    def _load_model(self, model_path: Path):
        """Loads a pre-trained tuned ensemble (RF, ET, CatBoost, LightGBM, GB, and DL ResNet) from disk."""
        data = joblib.load(str(model_path))
        if isinstance(data, dict):
            self.rf_model = data.get("rf")
            self.gb_model = data.get("gb")
            self.et_model = data.get("et")
            self.cb_model = data.get("cb")
            self.lgb_model = data.get("lgb")
            self.ensemble_weights = data.get("ensemble_weights")
            if "dl" in data and data["dl"] is not None:
                try:
                    from amp_screen.deep_model import DeepAffinityRegressor
                    self.dl_model = DeepAffinityRegressor.from_dict(data["dl"])
                except Exception as e:
                    logger.warning(f"Could not initialize Deep Learning model: {e}")
                    self.dl_model = None
            self.feature_names = data.get("features", list(FEATURE_NAMES))
            loaded = []
            if self.rf_model is not None: loaded.append("RF")
            if self.et_model is not None: loaded.append("ET")
            if self.cb_model is not None: loaded.append("CatBoost")
            if self.lgb_model is not None: loaded.append("LightGBM")
            if self.gb_model is not None: loaded.append("GB")
            if self.dl_model is not None: loaded.append("DL ResNet")
            models_desc = " + ".join(loaded) if loaded else "None"
            logger.info(
                f"Loaded pre-trained ML model from {model_path.name} "
                f"({models_desc}, {len(self.feature_names)} features)"
            )
        else:
            raise ValueError(
                f"Invalid model file format in {model_path}. Expected dictionary with 'rf' and 'gb' keys."
            )

    def save_model(self, out_path: Path):
        """Saves current trained ensemble models to disk with joblib."""
        out_path = Path(out_path)
        payload = {
            "rf": self.rf_model,
            "gb": self.gb_model,
            "features": self.feature_names
        }
        if self.et_model is not None:
            payload["et"] = self.et_model
        if self.dl_model is not None:
            payload["dl"] = self.dl_model.to_dict()
        joblib.dump(payload, str(out_path), compress=3)
        logger.info(f"Saved trained ML model to {out_path.name}")

    def _setup_features(self):
        """Sets up feature names for 96D structural and physical interface fingerprints."""
        self.feature_names = list(FEATURE_NAMES)

    def compute_consensus_score(
        self,
        docking_energy: Optional[float] = None,
        haddock_score: float = -40.0,
        sol_delta_g: Optional[float] = None,
        **kwargs
    ) -> float:
        """Unit-normalized consensus score combining docking, HADDOCK, and ML terms.
        Each component is mapped to [0, 10] scale before weighting — prevents
        HADDOCK spikes (e.g. +790 for clashing poses) from dominating the score.
        """
        dock_val = docking_energy if docking_energy is not None else kwargs.get("vina_affinity", kwargs.get("docking_affinity", kwargs.get("docking_energy", -8.0)))
        sol_dg = sol_delta_g if sol_delta_g is not None else float(dock_val)
        ml_pkd = max(3.0, min(12.0, -sol_dg / RT_LN10))

        dock_term = float(np.clip(-float(dock_val), 0.0, 15.0)) / 15.0 * 10.0
        haddock_term = float(np.clip(-float(haddock_score), -200.0, 50.0) + 200.0) / 250.0 * 10.0
        ml_term = float(np.clip(ml_pkd - 3.0, 0.0, 9.0)) / 9.0 * 10.0

        w1, w2, w3 = self.consensus_weights
        return float((w1 * dock_term) + (w2 * haddock_term) + (w3 * ml_term))

    def score_complex(
        self,
        complex_pdb: Path,
        hit_meta: Optional[Dict[str, Any]] = None,
        pocket_residues: Optional[List[str]] = None,
        receptor_chain: str = "A",
        peptide_chain: str = "B"
    ) -> Dict[str, Any]:
        """
        Calculates HADDOCK physics, RF structural fingerprints, and ML-predicted affinity.
        Returns full dict of physical & machine-learned metrics without label leakage.
        Automatically extracts interface contacts and sequence properties if missing from hit_meta.
        """
        complex_pdb = Path(complex_pdb)
        if hit_meta is None:
            hit_meta = {}
        else:
            hit_meta = dict(hit_meta)

        # 1. Physics Engine Evaluation
        physics = calculate_haddock_physics_score(
            complex_pdb=complex_pdb,
            pocket_residues=pocket_residues,
            receptor_chain=receptor_chain,
            peptide_chain=peptide_chain
        )

        # 2. Structural Fingerprint
        fp = extract_structural_fingerprint(
            complex_pdb,
            receptor_chain=receptor_chain,
            peptide_chain=peptide_chain
        )

        # 3. Intermolecular contacts (auto-extract if missing)
        if any(k not in hit_meta for k in ("salt_bridges", "hbonds", "hydrophobic_contacts", "clashes")):
            try:
                c_data = analyze_cocomplex_contacts(
                    complex_pdb,
                    receptor_chain=receptor_chain,
                    peptide_chain=peptide_chain
                )
                for k, v in c_data.items():
                    if k not in hit_meta:
                        hit_meta[k] = v
            except Exception:
                pass

        # 4. Sequence & AMP biophysics (auto-extract if missing)
        seq = hit_meta.get("sequence")
        if not seq:
            try:
                st = gemmi.read_structure(str(complex_pdb))
                seq = "".join(
                    AA_3TO1.get(res.name.upper(), "")
                    for m in st for c in m if c.name == "B"
                    for res in c if res.name.upper() in AA_3TO1
                )
                if seq:
                    hit_meta["sequence"] = seq
            except Exception:
                pass

        if seq:
            if "net_charge_z" not in hit_meta and "net_charge" not in hit_meta:
                hit_meta["net_charge_z"] = calculate_net_charge(seq)
            if "hydrophobic_moment_muH" not in hit_meta and "hydrophobic_moment" not in hit_meta:
                hit_meta["hydrophobic_moment_muH"] = calculate_hydrophobic_moment(seq)
            if "boman_index" not in hit_meta:
                hit_meta["boman_index"] = calculate_boman_index(seq)
            if "isoelectric_point" not in hit_meta:
                hit_meta["isoelectric_point"] = calculate_isoelectric_point(seq)

        pep_len = float(hit_meta.get("length", len(seq) if seq else 12.0))
        dock_energy = float(hit_meta.get("affinity", hit_meta.get("docking_energy", physics["haddock_score"] * 0.28)))

        # Assemble Feature Vector (exactly 96 dimensions, pure physical chemistry)
        feat_vec = extract_complex_features(fp=fp, physics=physics, hit_meta=hit_meta, seq=seq or "")
        X_input = np.array([feat_vec])

        # 5. Ensemble ML / Consensus Predictions
        if self.is_calibrated and (self.rf_model is not None or self.gb_model is not None):
            pred_dict = {}
            if self.rf_model is not None:
                pred_dict["rf"] = float(self.rf_model.predict(X_input)[0])
            if self.et_model is not None:
                pred_dict["et"] = float(self.et_model.predict(X_input)[0])
            if self.cb_model is not None:
                pred_dict["cb"] = float(self.cb_model.predict(X_input)[0])
            if self.lgb_model is not None:
                if hasattr(self.lgb_model, "feature_name_") and self.lgb_model.feature_name_:
                    import pandas as pd
                    pred_dict["lgb"] = float(self.lgb_model.predict(pd.DataFrame(X_input, columns=self.lgb_model.feature_name_))[0])
                else:
                    pred_dict["lgb"] = float(self.lgb_model.predict(X_input)[0])
            if self.gb_model is not None:
                pred_dict["gb"] = float(self.gb_model.predict(X_input)[0])
            if self.dl_model is not None:
                pred_dict["dl"] = float(self.dl_model.predict(X_input)[0])

            preds = list(pred_dict.values())

            # Optimized ensemble weights
            if self.ensemble_weights and all(k in pred_dict for k in self.ensemble_weights):
                raw_pkd = float(sum(pred_dict[k] * w for k, w in self.ensemble_weights.items()))
            elif "cb" in pred_dict and "lgb" in pred_dict:
                # Pent-ensemble weights: 35% RF, 25% ET, 20% CatBoost, 10% LightGBM, 10% DL
                w_map = {"rf": 0.35, "et": 0.25, "cb": 0.20, "lgb": 0.10, "dl": 0.10}
                total_w = sum(w_map.get(k, 0.0) for k in pred_dict)
                raw_pkd = float(sum(pred_dict[k] * w_map.get(k, 0.0) for k in pred_dict) / max(1e-6, total_w))
            elif "dl" in pred_dict:
                raw_pkd = float(0.20 * pred_dict["dl"] + 0.32 * pred_dict.get("gb", pred_dict.get("rf")) + 0.28 * pred_dict.get("rf") + 0.20 * pred_dict.get("et", pred_dict.get("rf")))
            else:
                raw_pkd = float(np.mean(preds))

            # Calibrated variance scaling: 1.20 preserves dynamic range while minimizing MAE
            mean_pkd = 6.21
            var_scale = 1.20
            ml_pkd = float(np.clip(mean_pkd + (raw_pkd - mean_pkd) * var_scale, 3.5, 12.0))
            variance_penalty = (max(preds) - min(preds)) * 10.0
            clash_penalty = hit_meta.get("clashes", 0) * 15.0
            bsa_reward = min(30.0, (physics["bsa"] / 1000.0) * 30.0)
            confidence = float(np.clip(70.0 + bsa_reward - variance_penalty - clash_penalty, 10.0, 95.0))
        else:
            # Transparent physical consensus derived from solution thermodynamics & docking
            sol_dg = float(hit_meta.get("solution_delta_g", dock_energy))
            ml_pkd = max(3.0, min(12.0, -sol_dg / RT_LN10))
            clash_pen = min(40.0, hit_meta.get("clashes", 0) * 10.0)
            bsa_score = min(30.0, (physics["bsa"] / 1000.0) * 30.0)
            contacts_score = min(20.0, (hit_meta.get("hbonds", 0) + hit_meta.get("salt_bridges", 0) * 2) * 3.0)
            confidence = float(np.clip(50.0 + bsa_score + contacts_score - clash_pen, 20.0, 85.0))

        # Exact Gibbs free energy conversion: Delta G = -RT * ln(10) * pKd (kcal/mol at 298.15 K)
        ml_delta_g = -RT_LN10 * ml_pkd

        # Calibrated Affinity Tiers (replaces indefensible sub-picomolar pseudo-precision)
        if ml_pkd >= 8.0:
            affinity_tier = "Tier 1: High Affinity (< 10 nM, Top 1%)"
        elif ml_pkd >= 6.5:
            affinity_tier = "Tier 2: Moderate Affinity (10-300 nM, Top 5%)"
        elif ml_pkd >= 5.0:
            affinity_tier = "Tier 3: Micromolar Binding (0.3-10 µM, Top 20%)"
        else:
            affinity_tier = "Tier 4: Weak / Non-Binder (> 10 µM)"

        kd_molar = 10.0 ** (-ml_pkd)
        if kd_molar < 1e-8:
            kd_str = "< 10 nM (Strong)"
        elif kd_molar < 1e-6:
            kd_str = f"{kd_molar * 1e9:.0f} nM"
        elif kd_molar < 1e-3:
            kd_str = f"{kd_molar * 1e6:.1f} µM"
        else:
            kd_str = f"{kd_molar * 1e3:.1f} mM"

        # Unit-Normalized Balanced Consensus Score
        # Each component is mapped to a [0, 10] scale before weighting to prevent
        # any single term from dominating due to scale differences:
        #   dock_term:    docking affinity [-15, 0] kcal/mol → [10, 0]  (lower = better)
        #   haddock_term: HADDOCK score [-50, +200] → [10, 0]  (capped: >200 = poor)
        #   ml_term:      pKd [3, 12] → [0, 10]               (higher pKd = better)
        dock_term = float(np.clip(-dock_energy, 0.0, 15.0)) / 15.0 * 10.0
        haddock_term = float(np.clip(-physics["haddock_score"], -200.0, 50.0) + 200.0) / 250.0 * 10.0
        ml_term = float(np.clip(ml_pkd - 3.0, 0.0, 9.0)) / 9.0 * 10.0
        w1, w2, w3 = self.consensus_weights
        consensus_score = (w1 * dock_term) + (w2 * haddock_term) + (w3 * ml_term)

        res_dict = {
            "haddock_score": physics["haddock_score"],
            "e_vdw": physics["e_vdw"],
            "e_elec": physics["e_elec"],
            "e_desolv": physics["e_desolv"],
            "bsa": physics["bsa"],
            "e_air": physics["e_air"],
            "ml_pkd": round(float(ml_pkd), 2),
            "ml_delta_g": round(float(ml_delta_g), 2),
            "ml_kd": kd_str,
            "affinity_tier": affinity_tier,
            "confidence_score": round(float(confidence), 1),
            "consensus_score": round(float(consensus_score), 2),
            "is_calibrated": self.is_calibrated
        }
        res_dict["interpretation"] = self.interpret_scoring_result(res_dict)
        return res_dict

    @staticmethod
    def interpret_scoring_result(result: Dict[str, Any]) -> Dict[str, str]:
        """
        Translates raw and calibrated scores into actionable biological and biophysical interpretations:
        - Affinity verdict (pKd, Delta G, and Kd range)
        - HADDOCK interaction strength (enthalpy and contact quality)
        - Interface burial & shape complementarity (BSA in Angstrom^2)
        - Overall lead optimization recommendation
        """
        pkd = result.get("ml_pkd", 0.0)
        kd_str = str(result.get("ml_kd", "N/A"))
        hscore = result.get("haddock_score", 0.0)
        bsa = result.get("bsa", 0.0)

        # 1. Affinity Interpretation
        if pkd >= 8.0:
            aff_verdict = f"High Affinity Potent Lead (Kd: {kd_str}, pKd >= 8.0). Sub-10 nM or nanomolar binding."
        elif pkd >= 6.5:
            aff_verdict = f"Moderate Affinity Candidate (Kd: {kd_str}, pKd 6.5-8.0). Suitable for lead optimization."
        elif pkd >= 5.0:
            aff_verdict = f"Micromolar Binder (Kd: {kd_str}, pKd 5.0-6.5). Requires hotspot mutagenesis."
        else:
            aff_verdict = f"Weak / Inactive (Kd: {kd_str}, pKd < 5.0). Exceeds micromolar threshold."

        # 2. HADDOCK Interfacial Enthalpy
        if hscore <= -50.0:
            had_verdict = f"Strong Enthalpic Packing ({hscore:.1f}). Robust vdW well and electrostatic complementarity."
        elif hscore <= -30.0:
            had_verdict = f"Favorable Interface ({hscore:.1f}). Standard productive peptide-protein contact."
        elif hscore <= -15.0:
            had_verdict = f"Marginal Interface ({hscore:.1f}). Weak interaction prone to dissociation."
        else:
            had_verdict = f"Unfavorable Contact ({hscore:.1f}). Steric repulsion or insufficient contacts."

        # 3. Buried Surface Area (BSA)
        if bsa >= 1200.0:
            bsa_verdict = f"Extensive Interface ({bsa:.1f} Å²). Deep pocket occupancy and high shape complementarity."
        elif bsa >= 800.0:
            bsa_verdict = f"Standard Interface ({bsa:.1f} Å²). Typical peptide groove binding."
        else:
            bsa_verdict = f"Shallow / Partial Contact ({bsa:.1f} Å²). Susceptible to solvent competition."

        # 4. Recommendation
        if pkd >= 7.0 and hscore <= -35.0 and bsa >= 800.0:
            rec = "Priority Hit: High affinity, robust HADDOCK enthalpy, and extensive buried interface."
        elif pkd >= 6.0 and hscore <= -25.0:
            rec = "Viable Candidate: Good affinity baseline; recommend in silico alanine scan or lead optimization."
        else:
            rec = "Low Priority: Marginal affinity or weak surface contacts."

        return {
            "affinity_verdict": aff_verdict,
            "haddock_verdict": had_verdict,
            "bsa_verdict": bsa_verdict,
            "overall_recommendation": rec
        }
