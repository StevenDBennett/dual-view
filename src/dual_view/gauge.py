"""
gauge.py
--------
Gauge invariants for weighted cyclic operators.

Given a cycle of integer weights w₀, …, w_{n-1}, the gauge invariant
quantities are built from the cyclic product W = Π w_i (mod 2^k) and
the spectral determinant det(D_w) = 1 - W (mod 2^k).  The 2-adic
dual-view decomposition of det(D_w) yields a scalar invariant
H = v + α·e that encodes the rigidity of the weight cycle under
tidal (2-adic) perturbations.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from .core import _mask, _valuation, two_adic_dlog


def cycle_product(weights: Sequence[int], k: int) -> int:
    """
    Cyclic product of weights modulo 2^k.

    Returns Π_{i} w_i (mod 2^k).
    """
    mask = _mask(k)
    prod = 1
    for w in weights:
        prod = (prod * int(w)) & mask
    return prod


def spectral_det(weights: Sequence[int], k: int, mod: Optional[int] = None) -> int:
    """
    Spectral determinant det(D_w) = 1 - W (mod 2^k).

    D_w is the weighted shift operator on the cycle; its determinant
    measures how far the cycle deviates from being a unitary.
    """
    if mod is not None:
        k = mod.bit_length() - 1
        if k < 3:
            k = 3
    W = cycle_product(weights, k)
    return (1 - W) & _mask(k)


def det_coordinates(weights: Sequence[int], k: int) -> Optional[Tuple[int, int, int]]:
    """
    Dual-view decomposition of the spectral determinant.

    Returns (v, α, e) or None if the determinant is zero.
    """
    det_val = spectral_det(weights, k)
    if det_val == 0:
        return None
    v = _valuation(det_val)
    odd = (det_val >> v) & _mask(k)
    result = two_adic_dlog(odd, k)
    if result is None:
        return None
    alpha, e = result
    return (v, alpha, e)


def tidal_scalar(weights: Sequence[int], k: int) -> Optional[int]:
    """
    Tidal scalar H = v + α·e — a single integer encoding the 2-adic
    rigidity of the weight cycle.

    Larger values of H correspond to deeper spectral determinent
    in the 2-adic congruence filtration, indicating greater numeric
    stability.
    """
    coords = det_coordinates(weights, k)
    if coords is None:
        return None
    v, alpha, e = coords
    return v + alpha * e


class GaugeLayer:
    """
    High-level interface for gauge analysis of a weight cycle.

    Pre-computes the product, spectral determinant, dual coordinates,
    and tidal scalar for a list of weights modulo 2^k.
    """

    def __init__(self, weights: Sequence[int], k: int) -> None:
        self.weights = list(weights)
        self.k = k
        self.mask = _mask(k)
        self.product = cycle_product(weights, k)
        self.det_val = (1 - self.product) & self.mask
        self.coords = det_coordinates(weights, k)
        self.tidal = tidal_scalar(weights, k)

    def report(self) -> str:
        """Return a formatted report of all gauge invariants."""
        lines = [
            "GaugeLayer Report",
            f"  k             = {self.k}",
            f"  weights       = {self.weights}",
            f"  product (W)   = {self.product}",
            f"  det(D_w)      = {self.det_val}",
        ]
        if self.coords is not None:
            v, alpha, e = self.coords
            lines.append(f"  coords        = (v={v}, α={alpha}, e={e})")
        else:
            lines.append("  coords        = None (det is zero)")
        lines.append(f"  tidal scalar  = {self.tidal}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"GaugeLayer({len(self.weights)} weights, k={self.k}, "
            f"product={self.product}, tidal={self.tidal})"
        )
