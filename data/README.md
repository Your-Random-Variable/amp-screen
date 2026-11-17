# Data Directory & Canonical Asset Organization

## Canonical Asset Location
To comply with Python packaging standards (PEP 517/621) and enable self-contained wheel installations, the primary runtime data assets are maintained directly within the Python package:

- **Canonical Pre-trained Ensemble Model:** [`amp_screen/data/trained_complex_ensemble.joblib`](../amp_screen/data/trained_complex_ensemble.joblib) (15 MB)
- **Canonical Training Dataset:** [`amp_screen/data/pdbbind_peptide_affinity_dataset.csv`](../amp_screen/data/pdbbind_peptide_affinity_dataset.csv) (1.3 MB)

In this directory (`data/`), the corresponding files are symbolic links pointing directly to the canonical files in `amp_screen/data/`. This eliminates data duplication in git and prevents the 15 MB model binary from being packaged twice in source distribution archives (`sdist`).

## Dataset Provenance & Terms
For detailed provenance information regarding the 1,054 curated co-crystal structures, public wwPDB coordinate sources, experimental affinity literature citations, and third-party database acknowledgments (PDBbind, PepBDB), see [`DATA_PROVENANCE_AND_LICENSE.md`](../DATA_PROVENANCE_AND_LICENSE.md).

## Dataset Regeneration
To regenerate the 96-dimensional physical descriptor matrix from scratch using public RCSB PDB coordinate downloads:

```bash
python data/generate_expanded_dataset.py
```
This script reads the curated entry list, fetches required co-crystal structures from RCSB PDB, calculates all 96 feature parameters in silico, and outputs the result directly to `amp_screen/data/pdbbind_peptide_affinity_dataset.csv`.
