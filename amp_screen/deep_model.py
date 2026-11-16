"""
Deep Learning Architecture for Protein-Peptide Affinity Prediction.
Provides a tabular ResNet-MLP model (DeepAffinityResNet) with residual connections,
LayerNorm, GELU activations, and smooth Huber loss optimization.
"""

from pathlib import Path
from typing import Optional, Dict, Any, Union
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader

from amp_screen.utils import logger


class TabularResBlock(nn.Module):
    """Residual block for continuous tabular physical features with LayerNorm and GELU."""

    def __init__(self, dim: int, dropout: float = 0.15):
        super().__init__()
        self.fc1 = nn.Linear(dim, dim)
        self.ln1 = nn.LayerNorm(dim)
        self.act1 = nn.GELU()
        self.drop1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(dim, dim)
        self.ln2 = nn.LayerNorm(dim)
        self.act2 = nn.GELU()
        self.drop2 = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        out = self.drop1(self.act1(self.ln1(self.fc1(x))))
        out = self.drop2(self.ln2(self.fc2(out)))
        return self.act2(out + residual)


class DeepAffinityResNet(nn.Module):
    """
    Deep Tabular ResNet for 96-dimensional physical & structural interface features.
    Architecture:
      Input (96D) -> Linear(96, 256) -> LayerNorm -> GELU -> Dropout
      -> TabularResBlock(256) -> TabularResBlock(256)
      -> Head: Linear(256, 64) -> GELU -> Dropout -> Linear(64, 1)
    """

    def __init__(
        self,
        in_features: int = 96,
        hidden_dim: int = 256,
        n_blocks: int = 2,
        dropout: float = 0.15
    ):
        super().__init__()
        self.in_features = in_features
        self.hidden_dim = hidden_dim

        self.input_layer = nn.Sequential(
            nn.Linear(in_features, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout)
        )

        self.res_blocks = nn.ModuleList([
            TabularResBlock(hidden_dim, dropout=dropout)
            for _ in range(n_blocks)
        ])

        self.head = nn.Sequential(
            nn.Linear(hidden_dim, 64),
            nn.LayerNorm(64),
            nn.GELU(),
            nn.Dropout(dropout * 0.75),
            nn.Linear(64, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.input_layer(x)
        for block in self.res_blocks:
            out = block(out)
        out = self.head(out)
        return out.squeeze(-1)


class DeepAffinityRegressor:
    """
    Scikit-learn compatible wrapper around DeepAffinityResNet.
    Features:
      - Automatic feature standardization (mean / std)
      - AdamW optimizer with CosineAnnealingLR
      - Huber / Smooth L1 loss (robust against outlier affinities)
      - Early stopping with validation checkpointing
      - Serializable to dictionary for joblib storage
    """

    def __init__(
        self,
        in_features: int = 96,
        hidden_dim: int = 64,
        n_blocks: int = 1,
        dropout: float = 0.35,
        lr: float = 1e-3,
        weight_decay: float = 5e-3,
        epochs: int = 60,
        batch_size: int = 32,
        device: Optional[str] = None,
        random_state: int = 42
    ):
        self.in_features = in_features
        self.hidden_dim = hidden_dim
        self.n_blocks = n_blocks
        self.dropout = dropout
        self.lr = lr
        self.weight_decay = weight_decay
        self.epochs = epochs
        self.batch_size = batch_size
        self.random_state = random_state

        if device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = device

        self.mean_: Optional[np.ndarray] = None
        self.std_: Optional[np.ndarray] = None
        self.model_: Optional[DeepAffinityResNet] = None
        self.best_loss_: float = float("inf")

    def _init_model(self):
        torch.manual_seed(self.random_state)
        self.model_ = DeepAffinityResNet(
            in_features=self.in_features,
            hidden_dim=self.hidden_dim,
            n_blocks=self.n_blocks,
            dropout=self.dropout
        ).to(self.device)

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        val_data: Optional[tuple] = None,
        verbose: bool = False
    ) -> "DeepAffinityRegressor":
        """Fits the DeepAffinityResNet model on feature matrix X and targets y."""
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.float32)

        # Compute standardization parameters
        self.mean_ = np.mean(X, axis=0)
        self.std_ = np.std(X, axis=0)
        self.std_[self.std_ < 1e-6] = 1.0  # Prevent divide-by-zero for constant features

        X_norm = (X - self.mean_) / self.std_

        self._init_model()
        self.model_.train()

        dataset = TensorDataset(torch.from_numpy(X_norm), torch.from_numpy(y))
        loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True, drop_last=False)

        optimizer = torch.optim.AdamW(
            self.model_.parameters(),
            lr=self.lr,
            weight_decay=self.weight_decay
        )
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=self.epochs, eta_min=1e-5
        )

        has_val = val_data is not None
        if has_val:
            X_val, y_val = val_data
            X_val_norm = (np.asarray(X_val, dtype=np.float32) - self.mean_) / self.std_
            t_X_val = torch.from_numpy(X_val_norm).to(self.device)
            t_y_val = torch.from_numpy(np.asarray(y_val, dtype=np.float32)).to(self.device)

        best_weights = None
        best_val_loss = float("inf")

        for epoch in range(self.epochs):
            self.model_.train()
            train_loss = 0.0
            for batch_x, batch_y in loader:
                batch_x = batch_x.to(self.device)
                batch_y = batch_y.to(self.device)

                optimizer.zero_grad()
                pred = self.model_(batch_x)
                loss = F.smooth_l1_loss(pred, batch_y, beta=0.5)
                loss.backward()
                nn.utils.clip_grad_norm_(self.model_.parameters(), max_norm=1.0)
                optimizer.step()
                train_loss += loss.item() * len(batch_x)

            scheduler.step()
            train_loss /= len(X)

            if has_val:
                self.model_.eval()
                with torch.no_grad():
                    val_pred = self.model_(t_X_val)
                    val_loss = F.smooth_l1_loss(val_pred, t_y_val, beta=0.5).item()

                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_weights = {k: v.cpu().clone() for k, v in self.model_.state_dict().items()}
            else:
                if train_loss < best_val_loss:
                    best_val_loss = train_loss
                    best_weights = {k: v.cpu().clone() for k, v in self.model_.state_dict().items()}

        if best_weights is not None:
            self.model_.load_state_dict(best_weights)
        self.model_.eval()
        self.best_loss_ = best_val_loss
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predicts affinity (pKd) for feature matrix X."""
        if self.model_ is None or self.mean_ is None or self.std_ is None:
            raise RuntimeError("DeepAffinityRegressor must be fitted or loaded before predict.")

        X = np.asarray(X, dtype=np.float32)
        X_norm = (X - self.mean_) / self.std_

        self.model_.eval()
        with torch.no_grad():
            t_x = torch.from_numpy(X_norm).to(self.device)
            preds = self.model_(t_x)
            return preds.cpu().numpy()

    def to_dict(self) -> Dict[str, Any]:
        """Exports model parameters and state for serialization."""
        if self.model_ is None:
            raise RuntimeError("Cannot serialize unfitted DeepAffinityRegressor.")
        return {
            "in_features": self.in_features,
            "hidden_dim": self.hidden_dim,
            "n_blocks": self.n_blocks,
            "dropout": self.dropout,
            "mean": self.mean_,
            "std": self.std_,
            "state_dict": {k: v.cpu() for k, v in self.model_.state_dict().items()}
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any], device: Optional[str] = None) -> "DeepAffinityRegressor":
        """Reconstructs a DeepAffinityRegressor from a state dictionary."""
        instance = cls(
            in_features=data["in_features"],
            hidden_dim=data["hidden_dim"],
            n_blocks=data["n_blocks"],
            dropout=data.get("dropout", 0.15),
            device=device
        )
        instance.mean_ = np.asarray(data["mean"], dtype=np.float32)
        instance.std_ = np.asarray(data["std"], dtype=np.float32)
        instance._init_model()
        instance.model_.load_state_dict(data["state_dict"])
        instance.model_.eval()
        return instance
