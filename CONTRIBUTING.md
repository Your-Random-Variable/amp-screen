# Contributing to AMP-Screen

Thank you for your interest in contributing to **AMP-Screen**! We welcome contributions ranging from bug fixes and feature enhancements to new benchmark complexes and documentation improvements.

---

## 🛠️ Development Setup

### 1. Clone the Repository
```bash
git clone https://github.com/Your-Random-Variable/amp-screen.git
cd amp-screen
```

### 2. Set Up a Virtual Environment

**Using Conda (Recommended):**
```bash
conda env create -f environment.yml
conda activate amp-screen
```

**Using Pip / venv:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
pip install -e ".[relaxation]"
```

---

## 🧪 Running the Test Suite

Before submitting any code changes, ensure all automated tests pass cleanly:

```bash
# Run the complete automated test suite
pytest tests/

# Run standardized co-crystal benchmarks
python run_benchmarks.py

# Run external blind validation benchmarks
python run_external_benchmark.py
```

All 60+ unit and integration tests must pass without regressions.

---

## 📐 Code Style & Conventions

- **Python Version:** Compatible with Python 3.9, 3.10, and 3.11.
- **Physical Coordinates:** All molecular coordinates must remain in 100% native Cartesian coordinates ($\text{\AA}$) using Gemmi and NumPy. Do not reintroduce empirical grid file formats or external binary calls to the core docking pipeline.
- **Thermodynamic Conventions:** Use standard reference temperature $T = 298.15\text{ K}$ ($25^\circ\text{C}$, $RT = 0.5925\text{ kcal/mol}$, $RT \ln(10) \approx 1.3642\text{ kcal/mol}$).
- **Documentation:** Preserve existing comments and docstrings. Add docstrings to all new functions describing input shapes, units, and expected return types.

---

## 🔀 Submitting a Pull Request (PR)

1. **Fork the Repository:** Create your feature branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. **Commit Your Changes:** Keep commits focused and provide clear, descriptive commit messages:
   ```bash
   git commit -m "feat: add support for cyclic stapled peptide perception"
   ```
3. **Verify Tests:** Ensure `pytest tests/` passes 100% green.
4. **Open a PR:** Push your branch to your fork and submit a Pull Request to `origin/main` on GitHub. Include a summary of changes, motivation, and any benchmark comparisons.

---

## 🐛 Reporting Issues

If you discover a bug, numerical instability, or unexpected behavior:
1. Check the [Issue Tracker](https://github.com/Your-Random-Variable/amp-screen/issues) to see if it has already been reported.
2. Open a new issue with:
   - Minimal reproducible example (input sequence, target PDB ID, CLI command).
   - Expected vs. observed output.
   - Operating system and Python environment details (`python --version`, dependency versions).

---

## 📄 License

By contributing to **AMP-Screen**, you agree that your contributions will be licensed under the project's [MIT License](LICENSE).
