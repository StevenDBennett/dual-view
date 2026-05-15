"""
training.py
-----------
Ghost-regularised PyTorch quantized training on MNIST.

Provides a QuantizedMLP with straight-through estimator (STE)
quantization and optional ghost regularization via GhostMap
or thermodynamics diagnostics via SeedThermodynamics.

PyTorch is an optional dependency — the module degrades gracefully
if torch is not available.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple
import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    import torch.optim as optim
    from torch.utils.data import DataLoader
    from torchvision import datasets, transforms
    _TORCH_AVAILABLE = True
except ImportError:
    _TORCH_AVAILABLE = False

from .core import _mask, modinv_newton
from .regularization import GhostMap, ghost_penalty


if _TORCH_AVAILABLE:
    class _STEQuantize(torch.autograd.Function):
        """Straight-Through Estimator rounding."""

        @staticmethod
        def forward(ctx, x, k):
            ctx.save_for_backward(x)
            scale = (1 << (k - 1)) - 1
            return torch.round(torch.clamp(x, -scale, scale))

        @staticmethod
        def backward(ctx, grad_output):
            return grad_output, None

    def ste_quantize(x: torch.Tensor, k: int) -> torch.Tensor:
        """Apply STE quantization to k bits."""
        return _STEQuantize.apply(x, k)

    class QuantizedMLP(nn.Module):
        """
        Two-layer MLP with STE-quantized weights.

        Architecture: 784 → 128 → 10 (MNIST).
        """

        def __init__(self, k: int = 8) -> None:
            super().__init__()
            self.k = k
            self.fc1 = nn.Linear(784, 128)
            self.fc2 = nn.Linear(128, 10)

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            x = x.view(-1, 784)
            h = F.relu(self.fc1(x))
            return self.fc2(h)

        def get_weights_numpy(self) -> Dict[str, np.ndarray]:
            """Extract quantized weights as numpy arrays."""
            weights = {}
            for name, param in self.named_parameters():
                if 'weight' in name:
                    w_int = (param.data * (1 << (self.k - 1))).round().to(torch.int32)
                    weights[name] = w_int.cpu().numpy()
            return weights

    def build_loaders(
        batch_size: int = 64, data_root: str = './data'
    ) -> Tuple[DataLoader, DataLoader]:
        """MNIST data loaders with standard normalisation."""
        transform = transforms.Compose([
            transforms.ToTensor(),
            transforms.Normalize((0.1307,), (0.3081,)),
        ])
        train_loader = DataLoader(
            datasets.MNIST(data_root, train=True, download=True, transform=transform),
            batch_size=batch_size, shuffle=True,
        )
        test_loader = DataLoader(
            datasets.MNIST(data_root, train=False, transform=transform),
            batch_size=batch_size, shuffle=False,
        )
        return train_loader, test_loader

    def train(
        model: QuantizedMLP,
        train_loader: DataLoader,
        test_loader: DataLoader,
        epochs: int = 5,
        lr: float = 0.001,
        ghost_map: Optional[GhostMap] = None,
        ghost_scale: float = 0.01,
        device: Optional[torch.device] = None,
        use_thermodynamics: bool = False,
        thermo_k: int = 8,
    ) -> Dict:
        """
        Training loop with optional ghost regularization.

        When use_thermodynamics=True, per-layer cliff_risk and
        alpha_fraction are tracked via SeedThermodynamics.

        Returns history dict with loss, accuracy, grad_norm,
        update_norm, ghost_penalty, cliff_risk, alpha_fraction.
        """
        if device is None:
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        model = model.to(device)
        optimizer = optim.Adam(model.parameters(), lr=lr)

        history: Dict[str, List] = {
            'loss': [], 'acc': [], 'grad_norm': [], 'update_norm': [],
            'ghost_penalty': [],
        }
        if use_thermodynamics:
            from .thermodynamics import SeedThermodynamics
            history['cliff_risk'] = []
            history['alpha_fraction'] = []

        for epoch in range(epochs):
            model.train()
            epoch_loss = 0.0
            batch_count = 0

            for data, target in train_loader:
                data, target = data.to(device), target.to(device)
                optimizer.zero_grad()
                output = model(data)
                loss = F.cross_entropy(output, target)

                # Ghost penalty
                penalty = 0.0
                if ghost_map is not None:
                    for name, param in model.named_parameters():
                        if 'weight' in name:
                            w_np = param.data.detach().cpu().numpy()
                            w_int = np.round(w_np * (1 << (model.k - 1))).astype(np.int32)
                            p, _ = ghost_penalty(w_int, ghost_map)
                            penalty += p

                    total_penalty = ghost_scale * penalty
                    loss = loss + total_penalty

                loss.backward()
                grad_norm = sum(
                    p.grad.norm().item() for p in model.parameters() if p.grad is not None
                )

                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)

                optimizer.step()
                update_norm = sum(
                    p.data.norm().item() for p in model.parameters()
                )

                epoch_loss += loss.item()
                batch_count += 1

            # Evaluation
            model.eval()
            correct = 0
            total = 0
            with torch.no_grad():
                for data, target in test_loader:
                    data, target = data.to(device), target.to(device)
                    output = model(data)
                    pred = output.argmax(dim=1)
                    correct += pred.eq(target).sum().item()
                    total += target.size(0)

            acc = correct / total
            avg_loss = epoch_loss / max(batch_count, 1)

            history['loss'].append(avg_loss)
            history['acc'].append(acc)
            history['grad_norm'].append(grad_norm)
            history['update_norm'].append(update_norm)
            history['ghost_penalty'].append(penalty)

            # Per-layer thermodynamics tracking
            thermo_str = ""
            if use_thermodynamics:
                st = SeedThermodynamics(k=thermo_k)
                cliffs = []
                alphas = []
                for name, param in model.named_parameters():
                    if 'weight' in name:
                        w_np = param.data.detach().cpu().numpy()
                        w_int = np.round(w_np * (1 << (model.k - 1))).astype(np.int64)
                        stats = st.analyse(w_int)
                        cliffs.append(stats["cliff_risk"])
                        alphas.append(stats["alpha_fraction"])
                epoch_cliff = float(np.mean(cliffs)) if cliffs else 0.0
                epoch_alpha = float(np.mean(alphas)) if alphas else 0.0
                history['cliff_risk'].append(epoch_cliff)
                history['alpha_fraction'].append(epoch_alpha)
                thermo_str = f", cliff={epoch_cliff:.3f}, α={epoch_alpha:.3f}"

            ghost_str = f", ghost={penalty:.4f}" if ghost_map else ""
            print(
                f"Epoch {epoch+1:2d}/{epochs}: "
                f"loss={avg_loss:.4f}, acc={acc:.4f}{ghost_str}{thermo_str}"
            )

        return history

else:
    class QuantizedMLP:
        """Placeholder when PyTorch is not available."""
        def __init__(self, k: int = 8):
            raise ImportError("PyTorch is required for QuantizedMLP.  Install with: pip install torch")

    def build_loaders(*args, **kwargs):
        raise ImportError("PyTorch is required for data loading.  Install with: pip install torch")

    def train(*args, **kwargs):
        raise ImportError("PyTorch is required for training.  Install with: pip install torch")
