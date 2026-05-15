"""
regularization.py
-----------------
Ghost-aware regularisation for quantized neural-network training.

The GhostMap pre-computes the Newton convergence fate for every odd
residue modulo 2^k and provides queries to find stable alternatives
and compute surrogate gradients for backpropagation.

k is limited to ≤ 10 because full enumeration over the 2^(k-2)
exponent domain becomes prohibitive at larger k.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from .basin import BasinExplorer
from .core import _mask, _valuation


class GhostMap:
    """
    Pre-compute Newton convergence ratios for all odd residues mod 2^k.

    Parameters
    ----------
    k : int
        Bit precision (3 ≤ k ≤ 10).  Limited to k ≤ 10 because
        full enumeration of 2^(k-2) seeds × odd values is O(2^k).
    g : int
        Generator; must satisfy g ≡ 5 (mod 8).
    max_iter : int
        Maximum Newton iterations per seed.  Default 64.
    """

    def __init__(self, k: int, g: int = 5, max_iter: int = 64) -> None:
        if k < 3:
            raise ValueError("k must be ≥ 3")
        if k > 10:
            raise ValueError(
                f"k={k} is too large for GhostMap.  "
                f"Full enumeration at k={k} requires analysing "
                f"2^{k-2} × 2^{k-2} ≈ {2**(2*k-4)} trajectories, "
                f"which is impractical.  Limit k ≤ 10."
            )
        self.k = k
        self.g = g
        self.max_iter = max_iter
        self.mask = _mask(k)
        self._ratio: Dict[int, float] = {}
        self._build()

    def _build(self) -> None:
        """Compute convergence ratios for all odd values mod 2^k."""
        for a in range(1, self.mask + 1, 2):  # odd values only
            try:
                explorer = BasinExplorer(self.k, self.g, a)
                portrait = explorer.portrait()
                n_converged = len(portrait['converged'])
                total = n_converged + len(portrait['cycle'])
                ratio = n_converged / total if total > 0 else 0.0
            except Exception:
                ratio = 0.0
            self._ratio[a] = ratio

    def ratio(self, a: int) -> float:
        """Return the convergence ratio for odd weight a."""
        a_int = int(a)  # coerce numpy.int64 etc.
        return self._ratio.get(a_int & self.mask, 0.0)

    def nearest_stable(self, a: int, search_radius: int = 4) -> Tuple[int, float]:
        """
        Find the nearest odd weight with a better convergence ratio.

        Returns (weight, ratio).
        """
        a_int = int(a) & self.mask
        best = (a_int, self.ratio(a_int))

        for delta in range(1, search_radius + 1):
            for candidate in (a_int - delta, a_int + delta):
                if candidate > 0 and (candidate & 1):
                    r = self.ratio(candidate)
                    if r > best[1]:
                        best = (candidate, r)
        return best

    def __repr__(self) -> str:
        return f"GhostMap(k={self.k}, n_odd={len(self._ratio)})"


def local_ratio_gradient(
    weight_int: int, ghost_map: GhostMap
) -> List[Tuple[int, float]]:
    """
    Compute local improvement candidates for a weight.

    Returns list of (delta, ratio) for delta in {-2, -1, 1, 2}
    that improve the convergence ratio.
    """
    current = ghost_map.ratio(weight_int)
    results: List[Tuple[int, float]] = []
    for delta in (-2, -1, 1, 2):
        candidate = int(weight_int) + delta
        if candidate > 0 and (candidate & 1):
            r = ghost_map.ratio(candidate)
            if r > current:
                results.append((delta, r))
    results.sort(key=lambda x: -x[1])
    return results


def ghost_penalty(
    weights: np.ndarray, ghost_map: GhostMap, step_scale: float = 1.0
) -> Tuple[float, np.ndarray]:
    """
    Ghost regularisation penalty and surrogate gradient.

    The penalty is P = mean(1 - ratio(w)) over the weight array.
    The surrogate gradient points toward the nearest odd integer
    with a better convergence ratio.

    Parameters
    ----------
    weights : np.ndarray
        Integer weight array.
    ghost_map : GhostMap
        Pre-computed convergence map.
    step_scale : float
        Scaling factor for the surrogate gradient.

    Returns
    -------
    (penalty, gradient) where gradient has the same shape as weights.
    """
    flat = weights.ravel()
    penalties = np.array([1.0 - ghost_map.ratio(int(w)) for w in flat])
    penalty = float(np.mean(penalties))

    grad = np.zeros_like(flat, dtype=np.float64)
    for i, w in enumerate(flat):
        w_int = int(w)
        if w_int & 1:
            improvers = local_ratio_gradient(w_int, ghost_map)
            if improvers:
                best_delta, _ = improvers[0]
                grad[i] = best_delta * step_scale

    return penalty, grad.reshape(weights.shape)
