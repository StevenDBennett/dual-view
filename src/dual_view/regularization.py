"""
regularization.py
-----------------
Ghost-aware regularisation for quantized neural-network training.

.. deprecated::
    The GhostMap convergence-ratio signal collapsed to uniform 1.0 after
    the α=1 sector fix (see bug_history.md).  The v₂(e_true) scores
    provided here are still useful as a graded stability diagnostic, but
    ``ghost_penalty`` and ``local_ratio_gradient`` return zero gradient
    for odd weights and should not be used for training.

    Use ``thermodynamics.SeedThermodynamics`` instead for genuine
    weight-stability diagnostics.

The GhostMap pre-computes the 2-adic exponent stability score v₂(e_true)
for every odd residue modulo 2^k and provides queries to find stable
alternatives and compute surrogate gradients for backpropagation.

Unlike the original convergence-ratio formulation (which collapsed to
uniform 1.0 after the α=1 fix), v₂(e_true) is a genuinely graded
stability measure: weights deeper in the 2-adic congruence filtration
(higher v₂(e_true)) are more stable under quantisation.

k is limited to ≤ 16 because the map enumerates all 2^{k-1} odd
residues — acceptable up to k=16 (32768 entries).
"""
from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Tuple
import numpy as np

from .core import _mask, _valuation, two_adic_dlog


class GhostMap:
    """
    Pre-compute 2-adic stability scores for all odd residues mod 2^k.

    The stability score for an odd weight a is based on
    v₂(e_true) — the valuation of the true discrete-log exponent:

        score = min(v₂(e_true), k-2) / (k-2)

    Higher scores mean the weight lies deeper in the 2-adic congruence
    filtration, hence more stable under quantisation.

    Parameters
    ----------
    k : int
        Bit precision (3 ≤ k ≤ 16).
    g : int
        Generator; must satisfy g ≡ 5 (mod 8).  Kept for API compatibility.
    """

    def __init__(self, k: int, g: int = 5) -> None:
        warnings.warn(
            "GhostMap is deprecated after the α=1 fix: convergence ratios "
            "are now uniform (1.0 for all odd weights). Use "
            "thermodynamics.SeedThermodynamics for graded v₂(e_true) "
            "stability diagnostics.",
            DeprecationWarning,
            stacklevel=2,
        )
        if k < 3:
            raise ValueError("k must be ≥ 3")
        if k > 16:
            raise ValueError(
                f"k={k} is too large for GhostMap.  "
                f"Enumerating all 2^{k-1} odd residues is impractical.  "
                f"Limit k ≤ 16."
            )
        self.k = k
        self.g = g
        self.mask = _mask(k)
        self._score: Dict[int, float] = {}
        self._build()

    def _build(self) -> None:
        """Compute stability scores for all odd values mod 2^k."""
        max_score = self.k - 2
        if max_score < 1:
            max_score = 1
        for a in range(1, self.mask + 1, 2):
            try:
                dlog = two_adic_dlog(a, self.k)
            except Exception:
                self._score[a] = 0.0
                continue
            if dlog is None:
                self._score[a] = 0.0
                continue
            _, e_true = dlog
            if e_true == 0:
                v2_e = self.k
            else:
                v2_e = _valuation(e_true)
            self._score[a] = min(v2_e, max_score) / max_score

    def ratio(self, a: int) -> float:
        """Return the stability score for weight a (handles even and zero)."""
        a_int = int(a) & self.mask
        if a_int == 0:
            return 1.0
        if a_int & 1:
            return self._score.get(a_int, 0.0)
        # even: use the 2-adic valuation of a directly
        v = _valuation(a_int)
        max_score = self.k - 2 if self.k > 2 else 1
        return min(v, max_score) / max_score

    def nearest_stable(self, a: int, search_radius: int = 4) -> Tuple[int, float]:
        """
        Find the nearest odd weight with a better stability score.

        Returns (weight, score).
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
        return f"GhostMap(k={self.k}, n_odd={len(self._score)})"


def local_ratio_gradient(
    weight_int: int, ghost_map: GhostMap
) -> List[Tuple[int, float]]:
    """
    Compute local improvement candidates for a weight.

    .. deprecated::
        This function is deprecated.  After the α=1 fix, GhostMap ratios
        are uniform and this returns no useful gradient.

    Returns list of (delta, score) for delta in {-2, -1, 1, 2}
    that improve the stability score.
    """
    warnings.warn(
        "local_ratio_gradient is deprecated: GhostMap ratios are now "
        "uniform after the α=1 fix. Use SeedThermodynamics instead.",
        DeprecationWarning,
        stacklevel=2,
    )
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

    .. deprecated::
        After the α=1 fix, this function returns zero gradient for all
        odd weights and is not useful for training.  Use
        ``thermodynamics.SeedThermodynamics`` for genuine stability
        diagnostics.

    The penalty is P = mean(1 - score(w)) over the weight array,
    where score is the v₂(e_true)-based stability from GhostMap.
    The surrogate gradient points toward the nearest odd integer
    with a better stability score.

    Parameters
    ----------
    weights : np.ndarray
        Integer weight array.
    ghost_map : GhostMap
        Pre-computed stability map.
    step_scale : float
        Scaling factor for the surrogate gradient.

    Returns
    -------
    (penalty, gradient) where gradient has the same shape as weights.
    """
    warnings.warn(
        "ghost_penalty is deprecated: after the α=1 fix, this returns "
        "zero gradient for odd weights. Use SeedThermodynamics instead.",
        DeprecationWarning,
        stacklevel=2,
    )
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
