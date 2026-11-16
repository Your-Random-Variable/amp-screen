import sys
import unittest
import tempfile
from pathlib import Path
if str(Path(__file__).parent.parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).parent.parent))
from amp_screen.benchmark import evaluate_docked_complex
try:
    from run_benchmarks import run_all_benchmarks
except ImportError:
    from amp_screen.benchmark import run_all_benchmarks

class TestBenchmark(unittest.TestCase):
    def test_capri_evaluation_benchmark(self):
        pdb_path = Path("benchmarks/1awr.pdb")
        if pdb_path.exists():
            metrics = evaluate_docked_complex(
                docked_pdb=pdb_path,
                reference_pdb=pdb_path,
                docked_rec_chain="A",
                docked_pep_chain="B",
                ref_rec_chain="A",
                ref_pep_chain="B"
            )
            self.assertLess(metrics["receptor_ca_rmsd"], 1e-3)
            self.assertLess(metrics["l_rmsd_backbone"], 1e-3)
            self.assertEqual(metrics["f_nat_percent"], 100.0)
            self.assertEqual(metrics["capri_quality"], "High (***)")

    def test_run_all_benchmarks(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_csv = Path(tmp_dir) / "benchmark_results.csv"
            df = run_all_benchmarks(benchmarks_dir=Path("benchmarks"), out_csv=out_csv)
            self.assertGreaterEqual(len(df), 6)
            self.assertEqual(len(df), len(list(Path("benchmarks").glob("*.pdb"))))
            self.assertTrue(out_csv.exists())
            self.assertIn("pdb_id", df.columns)
            self.assertIn("haddock_score", df.columns)

    def test_run_all_external_benchmarks(self):
        from run_external_benchmark import run_external_benchmark
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_csv = Path(tmp_dir) / "external_benchmark_results.csv"
            df = run_external_benchmark(benchmarks_dir=Path("benchmarks/external_test_set"), out_csv=out_csv)
            self.assertEqual(len(df), 52)
            self.assertTrue(out_csv.exists())
            self.assertIn("pdb_id", df.columns)
            self.assertIn("ml_delta_g", df.columns)
            self.assertIn("phys_delta_g", df.columns)
            self.assertIn("exp_delta_g", df.columns)

if __name__ == "__main__":
    unittest.main()
