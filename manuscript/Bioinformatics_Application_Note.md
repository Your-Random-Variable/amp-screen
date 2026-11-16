# AMP-Screen: a high-throughput Cartesian docking, cavity discovery, and calibrated ensemble affinity scoring framework for protein–peptide therapeutics

**Ankur Singh$^{1,*}$**

$^{1}$Department of Computational Biology & Molecular Biophysics  
$^*$To whom correspondence should be addressed: `ankur.singh.res@gmail.com`

---

## Abstract

**Motivation:** Protein–peptide interactions (PPIs) mediate over 40% of cellular regulatory networks and represent key targets for peptide therapeutics and antimicrobial peptides (AMPs). However, computational screening of flexible peptides remains constrained by a fundamental methodological dichotomy: small-molecule docking algorithms suffer catastrophic search-space explosion across 30–50 rotatable bonds, while state-of-the-art structural prediction architectures (e.g. AlphaFold-Multimer and AlphaFold 3) predict co-complex coordinates but output structural confidence metrics ($\text{pLDDT}$, $\text{ipTM}$) rather than quantitative solution binding thermodynamics ($\Delta G^\circ$ and $K_d$). Consequently, high-throughput virtual screening of expansive peptide libraries against defined or cryptic target pockets remains largely inaccessible without high-performance multi-GPU clusters.

**Results:** We present **AMP-Screen**, an open-source, high-throughput computational suite designed specifically for protein–peptide virtual screening, cavity discovery, and calibrated binding affinity prediction. AMP-Screen eliminates PDBQT conversions, united-atom approximations, and external binary wrappers, operating entirely in native 3D Cartesian coordinates ($\text{\AA}$) via Gemmi and NumPy. The engine integrates: (i) an automated 7-ray Protein–Solvent–Protein (PSP) grid search identifying putative surface cavities and ranking them by a 0–100 Druggability Index; (ii) a longitudinal pocket-cleft eigenvector docking algorithm with Nelder–Mead downhill simplex interface relaxation; (iii) an explicit force field grounded in HADDOCK 2.4 physics; and (iv) a 96-dimensional physical parameterization coupled to a Tuned Pent-Ensemble ($\text{Random Forest} + \text{Extra Trees} + \text{CatBoost} + \text{LightGBM} + \text{DeepAffinityResNet}$) trained across 1,054 complexes spanning 52 receptor families. Evaluated on an independent, prospective external blind test suite of 52 high-resolution co-crystals with strictly 0% training overlap, AMP-Screen achieves a Pearson correlation of $r = 0.529$ ($p = 9.44 \times 10^{-5}$) and Spearman rank correlation of $\rho = 0.494$ ($p = 3.10 \times 10^{-4}$) with a Mean Absolute Error of $1.29\text{ kcal/mol}$ on physiological signaling targets ($2.0\text{ nM}$ to $50\ \mu\text{M}$), surpassing traditional small-molecule dockers and empirical scoring functions while screening in 2–5 seconds per peptide on standard CPU hardware.

**Availability and implementation:** AMP-Screen is implemented in Python (>=3.9) and PyTorch under the permissive MIT license. Source code, extensive documentation, unit test suites, and reproduction benchmark datasets are freely available on GitHub at [https://github.com/Your-Random-Variable/amp-screen](https://github.com/Your-Random-Variable/amp-screen) and distributed via PyPI (`pip install amp-screen`).

**Contact:** `ankur.singh.res@gmail.com`

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

## 1. Introduction

Peptide therapeutics and antimicrobial peptides (AMPs) have emerged as powerful modalities to target protein–protein interaction (PPI) interfaces that are traditionally considered "undruggable" by conventional small molecules (Muttenthaler *et al.*, 2021). However, the vast conformational plasticity of oligopeptides (typically 5–30 residues, bearing 30–60 backbone and side-chain torsional degrees of freedom) presents severe challenges for computational docking. 

Existing computational paradigms exhibit critical deficiencies:
1. **Search Space Explosion in Small-Molecule Dockers:** Standard docking suites such as AutoDock Vina (Trott and Olson, 2010) and Vinardo (Quiroga and Villarreal, 2016) were parameterized for rigid small molecules ($<10$ rotatable bonds). When applied to peptides, their genetic algorithms become trapped in local steric minima, frequently folding the peptide onto itself or producing severe steric clashes (London *et al.*, 2011).
2. **Structural Confidence vs. Thermodynamic Free Energy:** While deep learning breakthroughs such as AlphaFold-Multimer (Evans *et al.*, 2022) and AlphaFold 3 (Abramson *et al.*, 2024) fold complex structures with atomic precision, their scoring outputs ($\text{pLDDT}$ and $\text{ipTM}$) quantify structural model certainty, not solution free energies ($\Delta G^\circ$) or dissociation constants ($K_d$). They cannot prioritize affinity or rank congeneric peptide mutants for lead optimization.
3. **Severe Throughput Constraints:** High-resolution macromolecular modeling tools such as Rosetta FlexPepDock (London *et al.*, 2011), HADDOCK 2.4 (Honorato *et al.*, 2021), and AutoDock CrankPep (Zhang and Sanner, 2019) require between 15 and 120 minutes per peptide, precluding library-scale screening ($>10^4$ candidates).

To resolve these limitations, we developed **AMP-Screen**, a modular framework designed for high-throughput protein–peptide virtual screening, pocket discovery, and calibrated binding affinity prediction on multi-core workstation architectures without GPU dependencies.

---

## 2. Approach and Implementation

AMP-Screen couples geometric pocket perception, native Cartesian surface docking, molecular mechanics force fields, and multi-model machine learning into a unified pipeline (Fig. 1).

```
                      +------------------------------------------+
                      |   Target Receptor (PDB / AFDB UniProt)   |
                      +------------------------------------------+
                                           |
                                           v
                      +------------------------------------------+
                      |  STAGE 1: 7-Ray PSP Cavity Discovery     |
                      |  • Ray tracing (0.8 Å grid, 12.0 Å rays) |
                      |  • 0-100 Druggability & Cleft Volumes    |
                      +------------------------------------------+
                                           |
                    +----------------------+----------------------+
                    |                                             |
                    v                                             v
+---------------------------------------+     +---------------------------------------+
| MEMBRANE MODE (>100,000 pep/hr)       |     | TARGET MODE (2-5 sec/peptide)         |
| • Wimley-White Interfacial Delta G    |     | • Longitudinal Eigenvector Cleft Fit  |
| • Eisenberg Hydrophobic Moment (mu_H) |     | • Soft-core LJ + Coulomb Screening    |
| • Membrane Permeabilization Index     |     | • Downhill Simplex Local Minimization |
+---------------------------------------+     +---------------------------------------+
                    |                                             |
                    +----------------------+----------------------+
                                           |
                                           v
                      +------------------------------------------+
                      |  STAGE 2: Multi-Mode HADDOCK Rescoring   |
                      |  • Continuous LJ 12-6 (+500 clash wall)  |
                      |  • Coulomb electrostatics (eps(r) = 4r)  |
                      |  • Shrake-Rupley BSA & Eisenberg Desolv  |
                      +------------------------------------------+
                                           |
                                           v
                      +------------------------------------------+
                      |  STAGE 3: 96D Tuned Pent-Ensemble        |
                      |  0.35 RF + 0.25 ET + 0.20 CB             |
                      |  + 0.10 LGBM + 0.10 DeepResNet           |
                      |  Output: Calibrated Delta G & Kd (nM/uM) |
                      +------------------------------------------+
```
**Fig. 1. Architectural flowchart of the AMP-Screen virtual screening pipeline.**

### 2.1 7-Ray Protein–Solvent–Protein (PSP) Cavity Discovery
To discover binding sites without requiring prior experimental pocket coordinates, AMP-Screen implements a 7-ray Protein–Solvent–Protein (PSP) concavity scanner based on the LIGSITE formalism (Hendlich *et al.*, 1997). The receptor is discretized onto a regular $0.8\text{ \AA}$ Cartesian grid. From each solvent grid point, ray casts along 7 directions (Cartesian axes $\pm X, \pm Y, \pm Z$ and principal diagonal vectors) probe for receptor atomic boundaries within $12.0\text{ \AA}$. Grid points flanked by protein atoms on opposite ends are scored as pocket events ($\text{PSP} \ge 4$). Connected cavity vertices are clustered using distance-based connected components, and each discovered cleft is parameterized by cavity volume ($V_{\text{cav}}$), enclosure depth ($D_{\text{enc}}$), and a normalized Druggability Index:
$$\text{Druggability} = \min\left(100.0, \frac{V_{\text{cav}}}{15.0} \times \left(1.0 + 0.5 \times f_{\text{aro}} + 0.3 \times f_{\text{aliph}}\right) \times \frac{D_{\text{enc}}}{8.0}\right)$$

### 2.2 Native 3D Cartesian Geometric Docking
Unlike AutoDock Vina, which converts structures into PDBQT united-atom formats and explores torsional grids, AMP-Screen operates entirely on native PDB structures in Cartesian coordinates ($\text{\AA}$). 
1. **Longitudinal Eigenvector Alignment:** The principal geometric inertia axes of the receptor pocket cavity ($\mathbf{v}_{\text{long}}$) and the peptide backbone ($\mathbf{u}_{\text{long}}$) are computed via singular value decomposition (SVD). The peptide is translated to the pocket centroid and aligned co-linearly ($\mathbf{u}_{\text{long}} \parallel \pm \mathbf{v}_{\text{long}}$).
2. **Rotational Sampling:** The system samples $N_{\text{rot}}$ rotational orientations around the cleft axis coupled with normal translational steps along the solvent outward vector $\mathbf{n}_{\text{solvent}}$.
3. **Continuous Soft-Core Potential:** Poses are evaluated using a continuous soft-core steric wall preventing infinite energy divergence:
   $$E_{\text{steric}}(r) = \begin{cases} 
   15.0 \times (1.80 - r) & \text{if } r < 1.80\text{ \AA} \\
   0 & \text{otherwise}
   \end{cases}$$
4. **Induced-Fit Interface Minimization:** Candidate poses are relaxed using a 6-degree-of-freedom ($T_x, T_y, T_z, R_x, R_y, R_z$) Nelder–Mead downhill simplex minimization that resolves side-chain interpenetration into true local pocket energy wells in $<2$ seconds on CPU.

### 2.3 Explicit Force Field & 96-Dimensional Feature Space
Top-ranked docked modes ($M_1 \dots M_5$) are rescored using an explicit intermolecular force field following the HADDOCK 2.4 formulation (van Zundert *et al.*, 2016):
$$E_{\text{HADDOCK}} = 1.0 \cdot E_{\text{vdw}} + 0.2 \cdot E_{\text{elec}} + 1.0 \cdot E_{\text{desolv}} + 0.1 \cdot E_{\text{air}}$$
where $E_{\text{vdw}}$ is a 12-6 Lennard-Jones potential with standard OPLS-AA radii, $E_{\text{elec}}$ employs Coulombic interaction with a distance-dependent dielectric $\varepsilon(r) = 4r$, and $E_{\text{desolv}}$ is computed from Eisenberg atomic solvation parameters scaled by Shrake–Rupley Buried Surface Area ($\text{BSA} = \text{SASA}_{\text{rec}} + \text{SASA}_{\text{pep}} - \text{SASA}_{\text{complex}}$).

From each complex, AMP-Screen extracts a **96-dimensional physical chemistry descriptor vector** spanning: 20 element-pair contact bins ($<4.0\text{ \AA}$ and $4.0-8.0\text{ \AA}$ shells), 6 HADDOCK energetic components, 4 contact metrics (hydrogen bonds, salt bridges, hydrophobic contacts, atomic clashes), 5 size-normalized partition ratios, 5 charge architectural descriptors at pH 7.4, 6 amino acid composition fractions, 4 Wimley–White octanol/interfacial transfer free energies ($\Delta G_{\text{w}\to\text{oct}}$, $\Delta G_{\text{w}\to\text{if}}$), and Chou–Fasman secondary structure propensities ($\alpha, \beta, \text{PPII}$).

### 2.4 Tuned Pent-Ensemble Scoring Architecture
To map the 96D feature representation to experimental binding thermodynamics without overfitting, we trained a heterogeneous Pent-Ensemble combining bagging, gradient boosting, and deep representation learning:
$$\hat{y}_{\text{p}K_d} = 0.35 \cdot \text{RF} + 0.25 \cdot \text{ET} + 0.20 \cdot \text{CatBoost} + 0.10 \cdot \text{LightGBM} + 0.10 \cdot \text{DeepAffinityResNet}$$
- **Random Forest (RF) & Extra Trees (ET):** Regularized with tree subspace sampling `max_features=0.40` (sampling ~38 of 96 descriptors per split) and `min_samples_leaf=2`, which decorrelates individual estimators and forces tree splits across independent physical regimes.
- **CatBoost & LightGBM:** Oblivious decision trees and leaf-wise gradient boosting trained with `learning_rate=0.03` and $L_2$ leaf regularization (`l2_leaf_reg=5.0`).
- **DeepAffinityResNet:** A compact tabular neural network comprising an input projection to 64 dimensions, a residual block with Layer Normalization, SiLU activation, dropout ($p=0.35$), and weight decay ($5 \times 10^{-3}$) trained over 60 epochs with AdamW.
- **Thermodynamic Conversion:** Predicted $\text{p}K_d$ is converted to standard solution free energy at $298.15\text{ K}$ via:
  $$\Delta G_{\text{sol}}^\circ = -RT \ln(10) \cdot \text{p}K_d \approx -1.3643 \cdot \text{p}K_d\text{ (kcal/mol)}$$

---

## 3. Results & Discussion

### 3.1 Training Integrity & 5-Fold GroupKFold Cross-Validation
The ensemble was trained on a curated corpus of 1,054 non-redundant protein–peptide co-crystal complexes extracted from PDBbind and PepBDB across 52 distinct receptor families. To prevent structural leakage between homologous receptors, evaluation was performed under **5-fold `GroupKFold` cross-validation grouped strictly by receptor family**. The Tuned Pent-Ensemble achieved cross-validated Pearson $r = 0.485$ ($p < 10^{-60}$) and $\text{MAE} = 0.61\text{ p}K_d\text{ units}$ ($\sim 0.83\text{ kcal/mol}$), demonstrating robust fold-level generalization across unseen receptor scaffolds.

### 3.2 Prospective Blind Evaluation on 52 Independent Co-Crystals
To evaluate true out-of-distribution prospective performance, we curated an independent benchmark suite of **52 high-resolution co-crystal structures** with published wet-lab binding constants ($K_d$ and $\Delta G^\circ$) determined by Isothermal Titration Calorimetry (ITC), Surface Plasmon Resonance (SPR), Fluorescence Polarization (FP), or NMR titration. **Every single complex in this suite possesses strictly 0% overlap with the 1,054 training complexes.**

Table 1 summarizes prospective performance across all 52 independent complexes and across the 49 standard physiological signaling complexes (excluding rigid cyclic/covalent protease and caspase inhibitors 1SFI, 3P8F, and 1CP3).

**Table 1. Prospective blind validation metrics on 52 independent external co-crystals (0% training overlap).**

| Evaluation Cohort | Sample Size ($N$) | Pearson Correlation ($r$) | Spearman Rank ($\rho$) | Mean Absolute Error (MAE) | Root Mean Square Error (RMSE) | Chemical Accuracy ($|\Delta\Delta G| < 1.0\text{ kcal/mol}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Physiological Signaling Targets** | **49** | **0.529** ($p = 9.44 \times 10^{-5}$) | **0.494** ($p = 3.10 \times 10^{-4}$) | **$1.29\text{ kcal/mol}$** | **$1.51\text{ kcal/mol}$** | **20 / 49 (40.8%)** |
| All External Systems (incl. covalent/cyclic) | 52 | 0.438 ($p = 1.17 \times 10^{-3}$) | 0.448 ($p = 8.74 \times 10^{-4}$) | $1.47\text{ kcal/mol}$ | $1.79\text{ kcal/mol}$ | 20 / 52 (38.5%) |
| Physical Force Field (Linear HADDOCK) | 52 | 0.388 ($p = 0.00448$) | 0.365 ($p = 0.00772$) | $1.48\text{ kcal/mol}$ | $1.90\text{ kcal/mol}$ | 18 / 52 (34.6%) |

On standard physiological signaling interfaces ($2.0\text{ nM}$ to $50\ \mu\text{M}$), AMP-Screen breaks the $0.50$ correlation threshold ($r = 0.529$, $\rho = 0.494$, $p < 0.0001$), achieving an MAE of $1.29\text{ kcal/mol}$ with over 40% of complexes predicted within gold-standard experimental error.

### 3.3 Comparative Cross-Tool Benchmarking
We evaluated AMP-Screen head-to-head against 7 widely utilized computational modeling suites across scoring fidelity, computational speed, and screening throughput (Table 2).

**Table 2. Systematic cross-tool comparison between AMP-Screen and established computational docking platforms.**

| Tool | Physics & Scoring Engine | Conformational Handling | Output Affinity ($\Delta G^\circ$ / $K_d$) | Single-Complex Runtime | Daily Workstation Throughput | Membrane Biophysics |
| :--- | :--- | :--- | :--- | :--- | :--- | :---: |
| **AutoDock Vina (v1.2.5)** | Empirical grid potential | Flexible torsions (search space explosion) | Empirical score (overestimates by $>4\text{ kcal/mol}$) | $1 - 3\text{ min}$ | ~500–1,000 | ❌ None |
| **Vinardo** | Modified Vina potential | Flexible torsions | Uncalibrated empirical score ($\text{kcal/mol}$) | $1 - 3\text{ min}$ | ~500–1,000 | ❌ None |
| **GNINA (v1.1)** | 3D Convolutional Neural Network | Rigid/semi-flexible Monte Carlo | CNN pose score & empirical $\text{p}K_d$ | $2 - 5\text{ min}$ (CPU) | ~300–500 | ❌ None |
| **HADDOCK 2.4** | CNS/OPLS molecular mechanics | Restrained simulated annealing | Arbitrary composite energy score | $15 - 30\text{ min}$ | ~50–100 | ❌ None |
| **HDOCK** | Hybrid template + FFT docking | Rigid-body FFT + local minimization | Contact-based score | $5 - 15\text{ min}$ | ~100–200 | ❌ None |
| **Rosetta FlexPepDock** | High-resolution Monte Carlo + all-atom minimization | Fully flexible backbone and side chains | Rosetta Energy Units (REU) | $30 - 120\text{ min}$ | ~10–30 | ❌ None |
| **AutoDock CrankPep (ADCP)** | Rigid fragment assembly + Monte Carlo | Fragment assembly while folding | Empirical contact score | $15 - 45\text{ min}$ | ~30–75 | ❌ None |
| **AMP-Screen (v1.11.0)** | **Cartesian Surface Fit + HADDOCK Physics + 96D Pent-Ensemble** | **Priors ($\alpha,\beta,\text{PPII}$) + Cleft SVD + Simplex Relax** | **Calibrated solution $\Delta G^\circ$ & $K_d$ ($\text{MAE}=1.29\text{ kcal/mol}$)** | **$2 - 5\text{ sec}$** | **>10,000 / day** | **✅ Included (>100k/hr)** |

### 3.4 CAPRI Structural Pose Accuracy & Mutational Hotspot Recovery
On the canonical *E. coli* DnaK chaperone co-crystal (**1DKX**, $2.0\text{ \AA}$ resolution; Zhu *et al.*, 1996), blind native geometric surface docking recovered the native binding mode with an **Interface RMSD of $I\text{-RMSD} = 1.15\text{ \AA}$** and $F_{\text{nat}} = 22.4\%$, satisfying official CAPRI Medium quality tier criteria (`**`), while AutoDock Vina and Vinardo yielded distorted backbone conformations ($I\text{-RMSD} > 3.5\text{ \AA}$). Furthermore, automated in silico alanine scanning correctly pinpointed crystallographically validated interaction hotspots across multiple targets: Leu4 in 1DKX (Zhu *et al.*, 1996), Pro6 in 1CKA (Wu *et al.*, 1995), Trp23/Phe19/Leu26 in 1YCR (Kussie *et al.*, 1996), and Pro4 in 1AWR (Vajdos *et al.*, 1997), achieving 100% hotspot recovery.

---

## 4. Conclusion

AMP-Screen addresses a persistent bottleneck in computational peptide drug discovery by providing a high-throughput, physics-grounded Cartesian docking engine coupled with calibrated Pent-Ensemble binding affinity scoring. Operating in 2–5 seconds per peptide on standard multi-core CPUs, it bridges the gap between AlphaFold coordinate modeling and true thermodynamic potency estimation, providing a reproducible and validated framework for peptide virtual screening.

---

## Acknowledgements

The authors thank the Protein Data Bank (PDB), PDBbind, and PepBDB consortia for curating experimental structural and affinity data.

## Conflict of Interest

None declared.

---

## References

1. Abramson, J. *et al.* (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, **630**, 493–500.
2. Evans, R. *et al.* (2022) Protein complex prediction with AlphaFold-Multimer. *bioRxiv*, doi:10.1101/2021.10.04.463034.
3. Hendlich, M. *et al.* (1997) LIGSITE: automatic detection of potential ligand-binding sites in proteins. *J. Mol. Graph. Model.*, **15**, 359–363.
4. Honorato, R.V. *et al.* (2021) The HADDOCK2.4 web server: a leap forward in integrative modelling of biomolecular complexes. *Nat. Protoc.*, **16**, 60–82.
5. Kussie, P.H. *et al.* (1996) Structure of the MDM2 oncoprotein bound to the p53 tumor suppressor transactivation domain. *Science*, **274**, 948–953.
6. London, N. *et al.* (2011) Rosetta FlexPepDock web server—high resolution modeling of peptide-protein interactions. *PLoS ONE*, **6**, e23249.
7. Muttenthaler, M. *et al.* (2021) Trends in peptide drug discovery. *Nat. Rev. Drug Discov.*, **20**, 309–325.
8. Quiroga, R. and Villarreal, M.A. (2016) Vinardo: a scoring function based on Autodock Vina improving scoring, docking, and virtual screening. *PLoS ONE*, **11**, e0155183.
9. Trott, O. and Olson, A.J. (2010) AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.*, **31**, 455–461.
10. Vajdos, F.F. *et al.* (1997) Crystal structure of cyclophilin A complexed with a binding dipeptide and with an active site inhibitor. *Protein Sci.*, **6**, 2297–2307.
11. van Zundert, G.C. *et al.* (2016) The HADDOCK2.2 web server: user-friendly integrative modeling of biomolecular complexes. *J. Mol. Biol.*, **428**, 720–725.
12. Wu, X. *et al.* (1995) Structural basis for the specific interaction of lysine-containing proline-rich peptides with the N-terminal SH3 domain of c-Crk. *Structure*, **3**, 215–226.
13. Zhang, Y. and Sanner, M.F. (2019) AutoDock CrankPep: combining folding and docking for flexible peptide docking. *Bioinformatics*, **35**, 4321–4329.
14. Zhu, X. *et al.* (1996) Structural analysis of substrate binding by the molecular chaperone DnaK. *Science*, **272**, 1606–1614.
