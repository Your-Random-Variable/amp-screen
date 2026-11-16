"""
Verification test for model packaging, resource resolution, and ensemble persistence.
"""

from pathlib import Path
import tempfile
import unittest
import joblib

from amp_screen.utils import get_default_model_path, get_default_dataset_path
from amp_screen.ml_scorer import MLPeptideScorer


class TestModelPackaging(unittest.TestCase):
    def test_default_resource_resolution(self):
        """Verifies that default model and dataset paths resolve to real files."""
        model_path = get_default_model_path()
        self.assertIsNotNone(model_path, "Default model path must not be None")
        self.assertTrue(Path(model_path).exists(), f"Model file must exist at {model_path}")
        self.assertGreater(Path(model_path).stat().st_size, 1000000, "Model file must be > 1 MB")

        dataset_path = get_default_dataset_path()
        self.assertIsNotNone(dataset_path, "Default dataset path must not be None")
        self.assertTrue(Path(dataset_path).exists(), f"Dataset file must exist at {dataset_path}")
        self.assertGreater(Path(dataset_path).stat().st_size, 500000, "Dataset file must be > 500 KB")

    def test_default_scorer_auto_loads_calibrated(self):
        """Verifies that MLPeptideScorer() initializes in calibrated mode with zero arguments."""
        scorer = MLPeptideScorer()
        self.assertTrue(scorer.is_calibrated, "MLPeptideScorer() must auto-load default model and be calibrated")
        self.assertIsNotNone(scorer.rf_model, "Random Forest model must be present")
        self.assertIsNotNone(scorer.cb_model, "CatBoost model must be present in Pent-Ensemble")
        self.assertIsNotNone(scorer.lgb_model, "LightGBM model must be present in Pent-Ensemble")
        self.assertIsNotNone(scorer.ensemble_weights, "Ensemble weights must be present")

    def test_model_save_load_round_trip(self):
        """Verifies that save_model and load_model preserve all ensemble components without dropping."""
        scorer = MLPeptideScorer()
        self.assertTrue(scorer.is_calibrated)

        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_model_path = Path(tmp_dir) / "test_roundtrip_ensemble.joblib"
            scorer.save_model(temp_model_path)
            self.assertTrue(temp_model_path.exists())

            # Inspect saved dictionary directly
            data = joblib.load(str(temp_model_path))
            self.assertIn("rf", data)
            self.assertIn("gb", data)
            self.assertIn("et", data)
            self.assertIn("cb", data)
            self.assertIn("lgb", data)
            self.assertIn("ensemble_weights", data)
            self.assertIn("features", data)
            self.assertEqual(len(data["features"]), 96)

            # Reload into a fresh scorer
            reloaded_scorer = MLPeptideScorer(model_path=temp_model_path)
            self.assertTrue(reloaded_scorer.is_calibrated)
            self.assertIsNotNone(reloaded_scorer.cb_model)
            self.assertIsNotNone(reloaded_scorer.lgb_model)
            self.assertEqual(reloaded_scorer.ensemble_weights, scorer.ensemble_weights)


if __name__ == "__main__":
    unittest.main()
