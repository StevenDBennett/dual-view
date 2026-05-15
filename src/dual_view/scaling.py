"""
scaling.py
----------
Float-to-integer scaling for quantized neural network weights.

Provides three rounding modes (round, floor, stochastic) and an
ensure_odd option that nudges even values so every weight lands
in the multiplicative group (Z/2^k)^× for 2-adic analysis.
"""
from __future__ import annotations

from typing import Dict, Tuple
import numpy as np

from .core import _valuation


def scale_weights(
    W: np.ndarray,
    scale: float,
    mode: str = "round",
    ensure_odd: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """
    Scale float weights to integers for 2-adic analysis.

    Parameters
    ----------
    W : np.ndarray
        Float weight array.
    scale : float
        Multiplicative scale factor.
    mode : str
        'round' (default), 'floor', or 'stochastic'.
    ensure_odd : bool
        If True, add 1 to even values so all weights are odd.

    Returns
    -------
    (W_int, meta) where W_int is int32 array and meta is a dict
    with keys: scale, mode, ensure_odd, v2_hist (2-adic valuation
    histogram for quick assessment).
    """
    if mode == "round":
        W_int = np.round(W * scale).astype(np.int32)
    elif mode == "floor":
        W_int = (np.abs(W * scale) * np.sign(W)).astype(np.int32)
    elif mode == "stochastic":
        frac = W * scale
        floor = np.floor(frac).astype(np.int32)
        prob = frac - floor
        mask = np.random.random(W.shape) < prob
        W_int = floor + mask.astype(np.int32)
    else:
        raise ValueError(f"Unknown mode '{mode}'.  Use round, floor, or stochastic.")

    if ensure_odd:
        even_mask = (W_int & 1) == 0
        W_int[even_mask] += 1

    n_even = int(np.sum((W_int & 1) == 0))
    n_zero = int(np.sum(W_int == 0))
    v2_vals = [_valuation(int(w)) if w != 0 else -1 for w in W_int.flat]
    v2_hist = {}
    for v in v2_vals:
        v2_hist[v] = v2_hist.get(v, 0) + 1

    meta = {
        "scale": scale,
        "mode": mode,
        "ensure_odd": ensure_odd,
        "n_even": n_even,
        "n_zero": n_zero,
        "range": (int(W_int.min()), int(W_int.max())),
        "v2_hist": v2_hist,
    }
    return W_int, meta


def auto_scale(W: np.ndarray, target_bits: int = 8) -> float:
    """
    Automatically choose a scale factor from the 99th percentile.

    Scale is chosen so that abs(W) * scale fits in (target_bits-1)
    bits (one bit reserved for sign).
    """
    p99 = float(np.percentile(np.abs(W), 99))
    if p99 == 0.0:
        return 1.0
    return (2 ** (target_bits - 1) - 1) / p99


def common_scales() -> Dict[str, float]:
    """Return a dict of standard quantization bit depths.
    
    For k-bit signed quantization, scale = 2^(k-1) — the maximum
    representable integer magnitude.
    """
    return {
        "INT7": 64.0,
        "INT8": 128.0,
        "INT9": 256.0,
        "INT10": 512.0,
        "INT12": 2048.0,
        "INT16": 32768.0,
    }
