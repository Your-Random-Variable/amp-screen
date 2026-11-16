# Supplementary Information: AMP-Screen

**Title:** AMP-Screen: a high-throughput Cartesian docking, cavity discovery, and calibrated ensemble affinity scoring framework for protein–peptide therapeutics  
**Author:** Ankur Singh  
**Journal:** *Bioinformatics* (Oxford Academic)  

---

## Supplementary Note S1: 96-Dimensional Physical Feature Space

The 96-dimensional physical descriptor vector is categorized into the following biophysical, physical chemistry, and geometric feature blocks:

1. **Intermolecular Element-Pair Contacts ($2 \times 10 = 20$ features):**
   - Distance bins for 10 atomic element pairs: $(\text{C-C}, \text{C-N}, \text{C-O}, \text{C-S}, \text{N-N}, \text{N-O}, \text{N-S}, \text{O-O}, \text{O-S}, \text{S-S})$.
   - Shell 1 (Close): $r \le 4.0\text{ \AA}$ (normalized by peptide length $L$).
   - Shell 2 (Medium): $4.0 < r \le 8.0\text{ \AA}$ (normalized by $L$).

2. **HADDOCK 2.4 Force Field Terms (6 features):**
   - $E_{\text{vdw}}$ (van der Waals Lennard-Jones 12-6).
   - $E_{\text{elec}}$ (Coulomb electrostatic energy with $\varepsilon(r) = 4r$).
   - $E_{\text{desolv}}$ (Eisenberg atomic solvation energy).
   - $E_{\text{HADDOCK}}$ (composite weighted linear combination).
   - $\text{BSA}$ (Buried Surface Area via Shrake-Rupley, $\text{\AA}^2$).
   - Total non-bonded intermolecular interaction energy ($E_{\text{vdw}} + E_{\text{elec}}$).

3. **Interface Contact Chemistry (4 features):**
   - Hydrogen bond count ($d_{\text{D-A}} \le 3.5\text{ \AA}$).
   - Salt bridge count ($d_{\text{N-O}} \le 4.0\text{ \AA}$ between charged Arg/Lys and Asp/Glu).
   - Hydrophobic contact count ($d \le 4.5\text{ \AA}$ between non-polar carbons).
   - Atomic clash count ($d < 2.0\text{ \AA}$).

4. **Normalized Energetics & Contact Densities (9 features):**
   - $E_{\text{vdw}} / L$, $E_{\text{elec}} / L$, $\text{BSA} / L$.
   - Contact density ($\text{Total Contacts} / L$).
   - Energetic partition ratios: $|E_{\text{elec}}| / \max(1, |E_{\text{vdw}}|)$, $\text{Polar Contacts} / \max(1, \text{Hydrophobic Contacts})$, $E_{\text{desolv}} / \max(1, \text{BSA})$.
   - Close-to-medium contact shell ratio.

5. **Charge & Electrostatic Architecture at pH 7.4 (5 features):**
   - Formal net charge $Z$ calculated via Henderson-Hasselbalch equation at physiological pH 7.4.
   - Positive charge count ($K+R$).
   - Negative charge count ($D+E$).
   - Charge density ($Z / L$).
   - Salt bridge potential ratio ($\min(N_{\text{pos}}, N_{\text{neg}}) / L$).

6. **Amino Acid Composition & Stereochemical Classes (12 features):**
   - Hydrophobic ratio ($\text{V, I, L, M, F, W, Y}$).
   - Net hydropathy index (Kyte-Doolittle GRAVY).
   - Sulfur-containing fraction ($C+M$).
   - Basic fraction ($K+R+H$).
   - Acidic fraction ($D+E$).
   - Hydroxyl fraction ($S+T+Y$).
   - Amide fraction ($N+Q$).
   - $\beta$-branched aliphatic fraction ($V+I+T$).
   - Proline fraction ($P / L$).
   - Glycine flexibility fraction ($G / L$).
   - Aromatic interaction density ($F+W+Y$).
   - Molecular weight ($M_w$) in Daltons.

7. **Membrane & Interfacial Partitioning Biophysics (4 features):**
   - Wimley-White water-to-octanol transfer free energy ($\Delta G_{\text{w}\to\text{oct}}$).
   - Wimley-White water-to-interface transfer free energy ($\Delta G_{\text{w}\to\text{if}}$).
   - Eisenberg hydrophobic moment ($\mu_H$).
   - Boman protein-binding index (kcal/mol).

8. **Secondary Structure Propensities (3 features):**
   - Chou-Fasman $\alpha$-helical propensity.
   - Chou-Fasman $\beta$-sheet propensity.
   - Polyproline-II (PPII) extended conformation propensity.

9. **Backbone Flexibility & Torsional Estimates (6 features):**
   - Estimated rotatable bonds ($N_{\text{rot}}$).
   - Conformational entropy penalty estimate ($\Delta S_{\text{conf}} \sim N_{\text{rot}} \cdot R \ln 3$).
   - Packing burial ratio ($\text{BSA} / \text{SASA}_{\text{pep}}$).
   - Interface packing density ($\text{Heavy Atoms} / \text{BSA}$).
   - Polar surface area fraction ($\text{PSA} / \text{BSA}$).
   - Aromatic burial density.

---

## Supplementary Table S1: Complete Prospective Blind External Benchmark ($N=52$)

Full tabular results for all 52 co-crystal complexes evaluated under blind conditions (0% training overlap) are deposited in `benchmarks/external_test_set/external_benchmark_results.csv` and documented in `BENCHMARKS.md`.
