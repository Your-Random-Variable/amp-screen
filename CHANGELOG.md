# Changelog

All notable changes to **AMP-Screen** (`amp-screen` / `af-amp-screen`) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.12.2] - 2026-11-17
### Added
- **Explicit Non-Redistribution Rights & Intellectual Property Notice**: Added formal legal declaration in `DATA_PROVENANCE_AND_LICENSE.md` certifying zero redistribution of third-party proprietary database archives, affirming CC0 public domain provenance of wwPDB structural coordinates, public factual status of literature binding constants under US 17 U.S.C. § 102(b) and EU Directive 96/9/EC, and 100% de novo algorithmic calculation of all 96 feature parameters.
- **Enhanced Wheel Smoke Test Assertions**: Updated `scripts/verify_wheel.py` and `tests/test_model_packaging.py` with explicit assertions verifying initialization of all Pent-Ensemble models: Random Forest, Extra Trees, CatBoost, LightGBM, Gradient Boosting fallback, PyTorch DeepAffinityResNet, and ensemble weights.
- **Clarified Isolation Scope**: Documented the scope and rationale of the clean-venv wheel validation testing workflow in `scripts/verify_wheel.py`.
- **Supported Platforms & CI Validation Documentation**: Explicitly documented formal Linux (Ubuntu x86_64, Python 3.10 and 3.11) continuous integration coverage alongside cross-platform macOS/Windows architectural support in `README.md`.

### Fixed
- **Retrained Ensemble Model Serialization**: Retrained and regenerated `trained_complex_ensemble.joblib` with current scikit-learn / joblib versions, eliminating unpickling warnings and ensuring clean cross-environment loading.
- **Corrected Legacy Script Docstring**: Fixed stale dataset size numbers in `data/curate_training_data.py` (539 complexes / 46 families → 1,054 complexes / 52 families) and documented its legacy status relative to `data/generate_expanded_dataset.py`.
- **Fixed Stale Paths in README**: Updated retrain command instructions in `README.md` to reference canonical `amp_screen/data/` paths.
- **Pytest Third-Party Warning Filters**: Configured pytest `filterwarnings` in `pyproject.toml` to suppress noisy NumPy deprecation warnings during automated test suite execution.
- **Missing Optional Import**: Added missing `Optional` typing import in `train_ml_model.py`.

---

## [1.12.1] - 2026-11-17
### Added
- **Dataset Provenance & Legal Rights Document**: Added [`DATA_PROVENANCE_AND_LICENSE.md`](DATA_PROVENANCE_AND_LICENSE.md) establishing explicit licensing distinctions between the permissive MIT software code, wwPDB CC0 primary coordinate files, factual literature binding constants, and ODC-By v1.0 / CC-BY 4.0 derived 96-dimensional physical descriptors.
- **Academic Third-Party Attribution**: Added formal citations and terms-of-use compliance notices for PDBbind (Wang et al., 2004, 2005; Liu et al., 2017), PepBDB (Wen et al., 2019), and wwPDB (Berman et al., 2000).
- **Automated Wheel & Clean Virtualenv Smoke Test**: Added [`scripts/verify_wheel.py`](scripts/verify_wheel.py) which verifies package data contents, creates an isolated temporary virtualenv, installs the built wheel, and asserts calibrated scorer loading in a sterile execution context.
- **Multi-Python CI Matrix**: Configured GitHub Actions CI across Python 3.10 and 3.11 with automated clean-venv wheel verification.

### Fixed
- **PyPI Availability Statement Reconciled**: Updated manuscript (`Bioinformatics_Application_Note.md` and `.tex`) and documentation to clarify that installation is currently available via GitHub releases / git checkout, with PyPI distribution scheduled upon peer-reviewed acceptance.
- **Eliminated Duplicate 16 MB Binary Assets**: Replaced duplicate copies of `trained_complex_ensemble.joblib` and `pdbbind_peptide_affinity_dataset.csv` in `data/` with relative symlinks pointing to canonical `amp_screen/data/`, reducing source distribution archive size from ~32 MB to ~16 MB.
- **Updated Documentation Test Count**: Synchronized all documentation (`README.md` and `CONTRIBUTING.md`) to reflect the complete 63 passing unit and integration tests.
- **Refined Core-vs-Optional Wording**: Clarified that AutoDock Vina, PDBQT, and Meeko are completely eliminated from the core native Cartesian docking engine, while preserving optional legacy interop utilities (`amp-screen-legacy`).

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
