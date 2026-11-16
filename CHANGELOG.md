# Changelog

All notable changes to **AMP-Screen** (`amp-screen` / `af-amp-screen`) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.12.0] - 2026-11-16
### Added
- **Package-Owned Runtime Resource Bundling**: Bundled the pre-trained Pent-Ensemble model (`trained_complex_ensemble.joblib`, 15 MB) and curated training corpus (`pdbbind_peptide_affinity_dataset.csv`, 1.3 MB) directly inside the package distribution (`amp_screen/data/`) using `importlib.resources` and `MANIFEST.in`.
- **Self-Contained Wheel Distribution**: Standard binary wheel installations (`pip install *.whl`) now out-of-the-box load the pre-trained calibrated Pent-Ensemble without requiring source checkouts or manual model downloads.
- **Complete Pent-Ensemble Persistence**: Updated `MLPeptideScorer.save_model()` to serialize all sub-models (Random Forest, Extra Trees, CatBoost, LightGBM, `DeepAffinityResNet`), feature definitions (96-dim), and learned ensemble weights.
- **Model Packaging Test Suite**: Added `tests/test_model_packaging.py` validating resource locator resolution, automatic calibrated model loading, and model serialization round-trips (63 passing tests).
- **CI Wheel Verification**: Added automated wheel build, metadata verification, and runtime asset inspection step to GitHub Actions CI workflow.

### Fixed
- **CI Dependency Isolation**: Guarded legacy Meeko PDBQT tests with `@pytest.mark.skipif(not (HAS_OBABEL and HAS_MEEKO))` to prevent false failures in minimal environments.
- **Standard Reference Temperature**: Unified `reporter.py` `estimate_kd()` default reference temperature from $300.0\text{ K}$ to standard IUPAC $298.15\text{ K}$ ($25^\circ\text{C}$, $RT = 0.5925\text{ kcal/mol}$).

---

## [1.11.0] - 2026-11-16
### Added
- **Curated 1,054 Co-Crystal Training Corpus**: Expanded machine learning dataset across 52 non-redundant receptor families with 96-dimensional physical parameterization.
- **Tuned Pent-Ensemble Engine**: Multi-architecture ensemble combining Random Forest, Extra Trees, CatBoost, LightGBM, and PyTorch `DeepAffinityResNet` with stacking meta-learner.
- **52-Complex Prospective External Blind Benchmark**: Strict zero-overlap, $<35\%$ sequence identity blind test suite ($r = 0.529$, $\rho = 0.494$, $\text{MAE} = 1.29\text{ kcal/mol}$ on 49 physiological signaling targets).
- **Expanded CAPRI Assessment**: Multi-target CAPRI pose recovery across 7 distinct structural fold classes (SH3, SH2, PDZ, MDM2, Cyclophilin A, Bromodomain, Chromodomain).
- **Consolidated Dependencies**: Unified dependencies across `pyproject.toml`, `setup.py`, `requirements.txt`, and `environment.yml` with explicit support for PyTorch, CatBoost, LightGBM, and XGBoost.
- **Fork-Friendly Automation**: Updated `update_github.sh` to preserve custom remotes for forks.

### Fixed
- Reconciled prospective blind MAE ($1.29\text{ kcal/mol}$) vs. retrospective reference calibration MAE ($0.12\text{ kcal/mol}$) across all documentation.
- Decoupled all unit tests from temporary pathogen test fixtures to use standardized RCSB reference co-crystals.
- Unified reference thermodynamic temperature at $T = 298.15\text{ K}$ ($25^\circ\text{C}$, $RT = 0.5925\text{ kcal/mol}$) across all modules.

---

## [1.10.0] - 2026-11-05
### Added
- **AlphaFold Ecosystem Integration (`amp-screen alphafold`)**:
  - `fetch`: Automated AlphaFold Database (AFDB) target structure retrieval and coordinate cleaning via UniProt ID.
  - `prep`: Automated AlphaFold-Multimer (ColabFold FASTA + bash runner) and AlphaFold 3 (alphafoldserver.com JSON) input generation.
  - `score`: Multi-term evaluation of co-folded complexes (interface pLDDT, HADDOCK physics, calibrated ML binding affinity, and clinical toxicity triage).

---

## [1.9.0] - 2026-10-31
### Added
- **7-Ray Protein–Solvent–Protein (PSP) Multi-Pocket Detection**: Automated blind cavity detector identifying top 5 surface pockets and scoring them with a 0–100 Druggability Index based on volume, buried concavity, and hydrophobicity.
- Direct blind docking support across detected cavities via `--blind-docking` and `--blind-pockets`.

---

## [1.8.0] - 2026-10-27
### Added
- **Sequence-First Membrane Biophysics (`amp-screen membrane`)**: Wimley-White interfacial ($\Delta G_{\text{w}\to\text{if}}$) and octanol ($\Delta G_{\text{w}\to\text{oct}}$) partitioning, Eisenberg hydrophobic moment ($\mu_H$), Membrane Permeabilization Index (MPI), and in silico therapeutic index screening at >100,000 peptides/hour without external model weights.

---

## [1.5.0] - 2026-10-26
### Added
- **Native 3D Geometric Surface Docking Engine**: 100% native Cartesian coordinate docking via Gemmi and NumPy, eliminating AutoDock Vina, PDBQT files, and Meeko runtime requirements.
- Nelder-Mead downhill simplex interface relaxation for local sidechain refinement into energetic sub-pockets.
- Monotonic Lennard-Jones repulsive wall capped at $+500\text{ kcal/mol}$ for severe sub-angstrom overlaps.
- Canonical C2 quintic switching function for electrostatics between 9 and 12 Å.

---

## [1.2.0] - 2026-10-25
### Added
- **In Silico Lead Optimizer (`amp-screen optimize-hit`)**: Automated alanine scanning for interface hotspot identification and multi-point rational mutation scanning.
- Force-field interface relaxation with frozen receptor backbone.

---

## [1.1.0] - 2026-10-24
### Added
- **HADDOCK 2.4 Explicit Physics Engine**: Exact implementation of intermolecular $E_{\text{vdw}}$, distance-dependent $E_{\text{elec}}$, Eisenberg-McLachlan $E_{\text{desolv}}$, Buried Surface Area (BSA), and steric clash detection.
- AMBER ff14SB partial charge table with formal charge consistency at physiological pH 7.4.

---

## [1.0.0] - 2026-10-23
### Added
- Initial release of AMP-Screen framework for high-throughput peptide virtual screening and biophysical property profiling.
