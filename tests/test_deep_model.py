"""
Unit tests for PyTorch Deep Learning ResNet architecture and hybrid ensemble in AMP-Screen.
"""

import pytest
import numpy as np
import torch
from pathlib import Path

from amp_screen.deep_model import DeepAffinityResNet, DeepAffinityRegressor, TabularResBlock
from amp_screen.ml_scorer import MLPeptideScorer, FEATURE_NAMES


def test_tabular_resblock_forward():
    block = TabularResBlock(dim=64, dropout=0.1)
    x = torch.randn(8, 64)
    out = block(x)
    assert out.shape == (8, 64)
    assert not torch.isnan(out).any()


def test_deep_affinity_resnet_architecture():
    model = DeepAffinityResNet(in_features=96, hidden_dim=128, n_blocks=2, dropout=0.1)
    x = torch.randn(4, 96)
    out = model(x)
    assert out.shape == (4,)
    assert not torch.isnan(out).any()

    # Verify gradients backpropagate
    loss = out.sum()
    loss.backward()
    for param in model.parameters():
        if param.requires_grad:
            assert param.grad is not None


def test_deep_affinity_regressor_fit_predict():
    np.random.seed(42)
    X = np.random.randn(40, 96).astype(np.float32)
    y = (X[:, 0] * 1.5 - X[:, 1] * 0.8 + 6.2).astype(np.float32)

    reg = DeepAffinityRegressor(
        in_features=96,
        hidden_dim=64,
        n_blocks=1,
        epochs=10,
        batch_size=16,
        random_state=42
    )
    reg.fit(X, y)

    preds = reg.predict(X[:5])
    assert len(preds) == 5
    assert not np.isnan(preds).any()

    # Serialization test
    state = reg.to_dict()
    assert "state_dict" in state
    assert "mean" in state
    assert "std" in state

    reg_loaded = DeepAffinityRegressor.from_dict(state)
    preds_loaded = reg_loaded.predict(X[:5])
    np.testing.assert_allclose(preds, preds_loaded, atol=1e-5)


def test_ml_scorer_with_hybrid_ensemble():
    model_path = Path("data/trained_complex_ensemble.joblib")
    if not model_path.exists():
        pytest.skip("Trained ensemble not found on disk")

    scorer = MLPeptideScorer(model_path=model_path)
    assert scorer.is_calibrated is True
    assert scorer.rf_model is not None
    assert scorer.gb_model is not None
    assert scorer.et_model is not None
    assert scorer.dl_model is not None
    assert len(scorer.feature_names) == 96
