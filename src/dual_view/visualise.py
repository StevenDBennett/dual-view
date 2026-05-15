"""
visualise.py
------------
Visualisation primitives for 2-adic cliff diagnostics.

Renders SeedThermodynamics cliff scores back into the original
weight-tensor shape and provides ASCII heatmap output.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from .core import _mask, _valuation, two_adic_dlog


def cliff_matrix(
    st, original_shape: Tuple[int, ...]
) -> np.ndarray:
    """
    Reshape SeedThermodynamics cliff scores to original shape.

    Even/zero weights get NaN.
    """
    flat = np.full(np.prod(original_shape), np.nan, dtype=np.float64)
    cliffs = st.cliffs
    for idx in range(len(cliffs)):
        c = cliffs.get(idx)
        if c is not None:
            flat[idx] = float(c)
    return flat.reshape(original_shape)


def sector_matrix(weights_int: np.ndarray, k: int) -> np.ndarray:
    """
    Map each odd weight to its α-sector (0 or 1), NaN for evens.

    The α-sector is the sign of the odd part: 0 for +5^e, 1 for -5^e.
    """
    flat = weights_int.ravel()
    result = np.full(len(flat), np.nan, dtype=np.float64)
    for i, w in enumerate(flat):
        w_int = int(w)
        if w_int & 1:
            result_i = two_adic_dlog(w_int & _mask(k), k)
            if result_i is not None:
                result[i] = float(result_i[0])
    return result.reshape(weights_int.shape)


def valuation_matrix(weights_int: np.ndarray) -> np.ndarray:
    """
    Map each weight to its 2-adic valuation v2.  Zeros get -1.
    """
    flat = weights_int.ravel()
    result = np.full(len(flat), -1.0, dtype=np.float64)
    for i, w in enumerate(flat):
        w_int = int(w)
        v = _valuation(w_int)
        result[i] = -1.0 if (v == float("inf")) else float(v)
    return result.reshape(weights_int.shape)


def print_cliff_ascii(
    C: np.ndarray,
    title: str = "Cliff Matrix",
    max_rows: int = 40,
    max_cols: int = 80,
) -> None:
    """
    Render a 2D cliff matrix as an ASCII density heatmap.

    Down-samples large matrices.
    """
    if C.ndim == 1:
        C_2d = C.reshape(1, -1)
    elif C.ndim >= 3:
        C_2d = C.reshape(C.shape[0], -1)
    else:
        C_2d = C

    if C_2d.shape[0] > max_rows:
        row_step = C_2d.shape[0] // max_rows + 1
        C_2d = C_2d[::row_step, :]
    if C_2d.shape[1] > max_cols:
        col_step = C_2d.shape[1] // max_cols + 1
        C_2d = C_2d[:, ::col_step]

    chars = " .:-=+*#%@"
    c_min = float(np.nanmin(C_2d)) if not np.all(np.isnan(C_2d)) else 0.0
    c_max = float(np.nanmax(C_2d)) if not np.all(np.isnan(C_2d)) else 1.0
    c_range = c_max - c_min if c_max > c_min else 1.0

    print(f"\n{title}  ({C_2d.shape[0]}×{C_2d.shape[1]})")
    print("┌" + "─" * C_2d.shape[1] + "┐")
    for row in C_2d:
        line = "│"
        for val in row:
            if np.isnan(val):
                line += "·"
            else:
                idx = int((val - c_min) / c_range * (len(chars) - 1))
                line += chars[min(idx, len(chars) - 1)]
        line += "│"
        print(line)
    print("└" + "─" * C_2d.shape[1] + "┘")


def cliff_stats_by_layer(layers: Dict[str, np.ndarray]) -> str:
    """
    Per-layer summary table of cliff statistics.

    Parameters
    ----------
    layers : dict of name → cliff matrix (from cliff_matrix)

    Returns formatted string.
    """
    lines = ["Layer Cliff Statistics", f"{'Layer':<20} {'Mean Cliff':>12} {'Min':>6} {'Max':>6} {'NaN%':>8}"]
    lines.append("-" * 54)

    for name, C in layers.items():
        valid = C[~np.isnan(C)]
        if len(valid) > 0:
            mean_c = float(np.mean(valid))
            min_c = float(np.min(valid))
            max_c = float(np.max(valid))
        else:
            mean_c = min_c = max_c = 0.0
        nan_pct = float(np.isnan(C).sum()) / C.size * 100
        lines.append(
            f"{name:<20} {mean_c:>12.2f} {min_c:>6.0f} {max_c:>6.0f} {nan_pct:>7.1f}%"
        )

    return "\n".join(lines)
