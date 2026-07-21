"""
bridge.py
---------
Unification of dual-view tidal coordinates and butterfly_v2-style
spectral seed analysis for neural network weight matrices.

Three complementary seeds per layer
------------------------------------
S1  DEPTH HISTOGRAM — circulant companion of (h - null)
S2  MAP SEED — W as a linear operator  (butterfly_v2 native)
S3  SIGN SEED — C_2 butterfly factor

When S1 and S2 agree: character is unambiguous.
When they split: that disagreement IS the finding.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "src")

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
from numpy.linalg import eigvals

from dual_view.core import _mask, _valuation

try:
    from scipy.linalg import circulant
except ImportError:
    circulant = None  # optional; needed for S1 seed only


# ═════════════════════════════════════════════════════════════════════════════
# 1.  SpectralThermodynamics  (spectral analysis of seed matrices)
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class SpectralThermodynamics:
    spectral_radius: float
    min_eigenvalue_magnitude: float
    max_eigenvalue_magnitude: float
    entropy_rate: float
    is_unitary: bool
    is_conservative: bool
    is_contractive: bool
    is_expansive: bool
    is_nilpotent: bool
    lyapunov_exponent: float

    @classmethod
    def analyze(cls, S: np.ndarray, tol: float = 1e-10) -> "SpectralThermodynamics":
        S = np.asarray(S, dtype=complex)
        eigs = eigvals(S)
        mags = np.abs(eigs)
        rho = float(np.max(mags))
        is_unitary = bool(np.allclose(S @ S.conj().T, np.eye(S.shape[0]), atol=tol))
        char_poly = np.poly(S)
        is_nilpotent = bool(np.all(np.abs(char_poly[1:]) < tol))
        return cls(
            spectral_radius=rho,
            min_eigenvalue_magnitude=float(np.min(mags)),
            max_eigenvalue_magnitude=rho,
            entropy_rate=float(np.log(max(rho, 1e-300))),
            is_unitary=is_unitary,
            is_conservative=bool(np.allclose(mags, 1.0, atol=tol)),
            is_contractive=bool(np.all(mags < 1 - tol)),
            is_expansive=bool(np.any(mags > 1 + tol)),
            is_nilpotent=is_nilpotent,
            lyapunov_exponent=float(np.log(rho)) if rho > 0 else -np.inf,
        )

    def character(self) -> str:
        if self.is_nilpotent:    return "NILPOTENT"
        if self.is_unitary:      return "CONSERVATIVE/UNITARY"
        if self.is_conservative: return "CONSERVATIVE"
        if self.is_contractive:  return "CONTRACTIVE"
        if self.is_expansive:    return "EXPANSIVE"
        return "MIXED"

    def __str__(self) -> str:
        return (f"{self.character():<22} "
                f"\u03c1={self.spectral_radius:.6f}  "
                f"\u03bb={self.lyapunov_exponent:.4f}")


# ═════════════════════════════════════════════════════════════════════════════
# 2.  Quantisation
# ═════════════════════════════════════════════════════════════════════════════


def quantize(W: np.ndarray, k: int) -> np.ndarray:
    """Float tensor to Z/2^k Z via symmetric min-max quantisation."""
    levels = 1 << k
    lo, hi = W.min(), W.max()
    if hi == lo:
        return np.zeros_like(W, dtype=np.int64)
    scaled = (W - lo) / (hi - lo) * (levels - 1)
    return ((scaled - levels // 2).astype(np.int64)) % levels


def _depths(W_int: np.ndarray, k: int) -> np.ndarray:
    """2-adic valuation of every element, capped at k-1."""
    mask = _mask(k)
    flat = W_int.ravel().astype(np.int64) & mask
    return np.array(
        [min(k - 1, _valuation(int(w)) if w != 0 else k - 1) for w in flat],
        dtype=int,
    )


# ═════════════════════════════════════════════════════════════════════════════
# 3.  The three seeds
# ═════════════════════════════════════════════════════════════════════════════


def _geometric_null(k: int) -> np.ndarray:
    """P(v=j) = 1/2^(j+1) for pure-random integers mod 2^k."""
    h = np.array([0.5 ** (j + 1) for j in range(k)])
    h[-1] += 1.0 - h.sum()
    return h


def seed_S1(W_int: np.ndarray, k: int):
    """Depth-histogram circulant seed. Returns (h, dev, S1, thermo)."""
    if circulant is None:
        raise ImportError("scipy.linalg.circulant is required for seed_S1. "
                          "Install with: pip install dual-view[scipy]")
    d = _depths(W_int, k)
    h = np.bincount(d, minlength=k).astype(float)
    h /= h.sum() + 1e-12
    null = _geometric_null(k)
    dev = h - null
    S1 = circulant(dev)
    thermo = SpectralThermodynamics.analyze(S1)
    return h, dev, S1, thermo


def seed_S2(W_float: np.ndarray):
    """Map seed — W as a linear operator. Returns (S2, thermo)."""
    W = np.asarray(W_float, dtype=float)
    W_n = W / (np.abs(W).max() + 1e-12)
    S2 = W_n if W_n.shape[0] == W_n.shape[1] else (W_n @ W_n.T) / W_n.shape[1]
    thermo = SpectralThermodynamics.analyze(S2)
    return S2, thermo


def seed_S3(W_int: np.ndarray, k: int):
    """Sign-branch seed — C_2 butterfly factor. Returns (S3, thermo)."""
    mask = _mask(k)
    flat = W_int.ravel().astype(np.int64) & mask
    nz = [int(w) for w in flat if w != 0]
    p1 = sum(1 for w in nz if (w >> 1) & 1) / (len(nz) + 1e-12)
    p0 = 1.0 - p1
    S3 = np.array([[p0, p1], [p1, p0]])
    thermo = SpectralThermodynamics.analyze(S3)
    return S3, thermo


def depth_char(h: np.ndarray, dev: np.ndarray, k: int) -> str:
    """Classify depth deviation relative to geometric null."""
    if circulant is None:
        raise ImportError("scipy.linalg.circulant is required for depth_char. "
                          "Install with: pip install dual-view[scipy]")
    rho_c = float(np.max(np.abs(eigvals(circulant(dev).astype(complex)))))
    prefix = "STRUCTURED/" if rho_c > 0.3 else ""
    if dev[:2].mean() > 0.03:    return prefix + "EXPANSIVE"
    if dev[k // 2:].mean() > 0.03: return prefix + "CONTRACTIVE"
    return prefix + "NEUTRAL"


# ═════════════════════════════════════════════════════════════════════════════
# 4.  Layer and model reports
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class LayerReport:
    name:       str
    shape:      Tuple[int, ...]
    k:          int
    n:          int
    mean_v:     float
    depth_entropy: float
    zero_frac:  float
    alpha_frac: float
    depth_hist: np.ndarray
    depth_dev:  np.ndarray
    depth_char: str
    thermo_S1:  SpectralThermodynamics
    thermo_S2:  SpectralThermodynamics
    thermo_S3:  SpectralThermodynamics

    def consensus(self) -> str:
        s1c = self.depth_char.split("/")[-1]
        s2c = self.thermo_S2.character().split("/")[0]
        relevant = {c for c in (s1c, s2c) if c not in ("NEUTRAL", "NILPOTENT")}
        if not relevant:           return "NEUTRAL  (matches arithmetic null)"
        if len(relevant) == 1:     return f"CONSENSUS: {relevant.pop()}"
        return f"SPLIT   depth={s1c}  map={s2c}"

    def report(self) -> str:
        null4 = _geometric_null(self.k)[:4]
        lines = [
            f"  {self.name}",
            f"    shape={self.shape}  n={self.n:,}  k={self.k}",
            f"    S0  TIDAL: mean_v={self.mean_v:.3f}  depth_entropy={self.depth_entropy:.3f}  "
            f"zero_frac={self.zero_frac:.3f}  alpha_frac={self.alpha_frac:.3f}",
            f"    S1  DEPTH HISTOGRAM: {self.thermo_S1}",
            f"      char = {self.depth_char}",
            f"      h[:4]={self.depth_hist[:4].round(4)}  dev[:4]={self.depth_dev[:4].round(4)}",
            f"    S2  MAP SEED: {self.thermo_S2}",
            f"    S3  SIGN SEED: {self.thermo_S3}",
            f"    CONSENSUS: {self.consensus()}",
        ]
        return "\n".join(lines)


@dataclass
class ModelReport:
    layers: List[LayerReport] = field(default_factory=list)

    def add(self, r: LayerReport) -> None:
        self.layers.append(r)

    def _traj(self, fn) -> np.ndarray:
        return np.array([fn(r) for r in self.layers])

    def boundaries(self) -> List[str]:
        out = []
        for a, b in zip(self.layers, self.layers[1:]):
            ca = a.thermo_S2.character()
            cb = b.thermo_S2.character()
            if ca != cb:
                out.append(f"  {a.name}  ->  {b.name}  [{ca} -> {cb}]")
        return out

    def report(self) -> str:
        lines = [
            "UNIFIED BUTTERFLY / TIDAL ANALYSIS",
            "Three-seed 2-adic structure of neural network weights",
            "",
        ]
        for r in self.layers:
            lines += [r.report(), ""]

        lya1 = self._traj(lambda r: r.thermo_S1.lyapunov_exponent)
        lya2 = self._traj(lambda r: r.thermo_S2.lyapunov_exponent)
        mv = self._traj(lambda r: r.mean_v)

        lines += [
            "TRAJECTORY (per layer):",
            f"  {'Layer':<26} {'lambda(S1)':>8} {'lambda(S2)':>8} {'mean_v':>7}  consensus",
            "-" * 72,
        ]
        for r, l1, l2, d in zip(self.layers, lya1, lya2, mv):
            lines.append(
                f"  {r.name:<26} {l1:>8.3f} {l2:>8.3f} {d:>7.3f}  {r.consensus()}"
            )

        bounds = self.boundaries()
        lines.append("")
        if bounds:
            lines += ["Butterfly boundaries (S2 map transitions):"] + bounds
        else:
            lines.append("No S2 transitions detected across layers.")

        splits = [r.name for r in self.layers if "SPLIT" in r.consensus()]
        if splits:
            lines += ["", "Split layers (depth != map):"]
            lines += [f"  {s}" for s in splits]

        return "\n".join(lines)


# ═════════════════════════════════════════════════════════════════════════════
# 5.  Bridge — main entry point
# ═════════════════════════════════════════════════════════════════════════════


class ButterflyBridge:
    """
    Unified butterfly_v2 + tidal_coordinates weight analyser.

    Parameters
    ----------
    k : int   bit-width for Z/2^k Z  (default 8)
    """

    def __init__(self, k: int = 8):
        self.k = k

    def analyse_layer(self,
                      W_float: np.ndarray,
                      name: str = "unnamed") -> LayerReport:
        k = self.k
        W_i = quantize(W_float, k)
        d = _depths(W_i, k)
        mask = _mask(k)
        flat = W_i.ravel().astype(np.int64) & mask

        mean_v = float(d.mean())
        h_ent = np.bincount(d, minlength=k).astype(float)
        h_ent /= h_ent.sum() + 1e-12
        depth_entropy = float(-np.sum(h_ent * np.log2(h_ent + 1e-12)))
        zero_frac = float(np.sum(flat == 0)) / len(flat)
        nz = [int(w) for w in flat if w != 0]
        alpha_frac = (sum(1 for w in nz if (w >> 1) & 1) / (len(nz) + 1e-12))

        h, dev, S1, t1 = seed_S1(W_i, k)
        dc = depth_char(h, dev, k)
        S2, t2 = seed_S2(W_float)
        S3, t3 = seed_S3(W_i, k)

        return LayerReport(
            name=name, shape=W_float.shape, k=k,
            n=W_float.size, mean_v=mean_v, depth_entropy=depth_entropy,
            zero_frac=zero_frac, alpha_frac=alpha_frac,
            depth_hist=h, depth_dev=dev, depth_char=dc,
            thermo_S1=t1, thermo_S2=t2, thermo_S3=t3,
        )

    def analyse_model(self, layers: Dict[str, np.ndarray]) -> ModelReport:
        report = ModelReport()
        for name, W in layers.items():
            report.add(self.analyse_layer(W, name=name))
        return report
