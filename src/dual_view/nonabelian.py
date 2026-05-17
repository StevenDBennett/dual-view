"""
nonabelian.py
-------------
Non-Abelian CRT-dual system with GL(2) holonomy.

Extends the gauge-theoretic invariants to matrix-valued weights
acting on a cycle over Z/(2^k · p)Z.  The holonomy H = Π M_i
lives in GL(2), and its determinant (mod 2^k) and trace (mod p)
provide crossed invariants that combine via CRT.

Note on ramp_break_strength
---------------------------
The original eps_crit metric in ramp_break_strength is degenerate:
for a shear matrix [[1+ε, 1], [0, 1]], det = 1+ε, and the smallest
ε giving an even determinant is always 1.  The genuinely interesting
signal is whether the perturbation flips the α-sector of the
determinant (phase alignment).  Use phase_alignment_experiment()
for this purpose.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from .core import _mask, _valuation, modinv_newton, two_adic_dlog, _mat_mul, _mat_det
from .crt import CRTDualProcessor


def _perturbation_matrix(d_target: int, mod_full: int) -> List[List[int]]:
    """Create a shear matrix [[d, 1], [0, 1]] mod mod_full."""
    return [[d_target % mod_full, 1], [0, 1]]


class NonAbelianCRTDual:
    """
    Non-Abelian (matrix-valued) CRT-dual system.

    Computes holonomy for a cycle of 2×2 matrices over
    Z/(2^k · p)Z, extracts joint invariants, and tests
    their correlation with 2-adic ghost stability.
    """

    def __init__(self, k: int, p: int) -> None:
        self.k = k
        self.p = p
        self.mod2 = 1 << k
        self.mod_full = self.mod2 * p
        self.crt = CRTDualProcessor(k, p)

    def holonomy(self, mats: List[List[List[int]]]) -> List[List[int]]:
        """Product of all matrices in the cycle modulo mod_full."""
        H = [[1, 0], [0, 1]]
        for M in mats:
            H = _mat_mul(H, M, self.mod_full)
        return H

    def invariants(self, mats: List[List[List[int]]]) -> Dict[str, int]:
        """
        Extract invariants: det mod 2^k, det dual view, trace mod p,
        and CRT merge.
        """
        H = self.holonomy(mats)
        det_mod2k = _mat_det(H, self.mod2)
        trace_modp = (H[0][0] + H[1][1]) % self.p

        # Dual-view decomposition of determinant
        det_odd = det_mod2k if (det_mod2k & 1) else det_mod2k + 1
        dlog_result = two_adic_dlog(det_odd, self.k)
        if dlog_result is not None:
            alpha_det, e_det = dlog_result
        else:
            alpha_det, e_det = 0, 0

        # CRT merge
        inv_2k = self.crt.crt_reconstruct(det_mod2k, trace_modp)
        inv_p = self.crt.crt_reconstruct(trace_modp, det_mod2k)

        return {
            "det_mod2k": det_mod2k,
            "alpha_det": alpha_det,
            "e_det": e_det,
            "trace_modp": trace_modp,
            "crt_2k_view": inv_2k,
            "crt_p_view": inv_p,
        }

    def convergence_ratio_full(self, mats: List[List[List[int]]]) -> float:
        """Convergence ratio of the holonomy determinant 2-adic component."""
        from .basin import BasinExplorer

        H = self.holonomy(mats)
        det = _mat_det(H, self.mod2)
        if det & 1 == 0:
            return 0.0
        try:
            explorer = BasinExplorer(self.k, 5, det)
            portrait = explorer.portrait()
            n_converged = len(portrait['converged'])
            total = n_converged + len(portrait['cycle'])
            return n_converged / total if total > 0 else 0.0
        except Exception:
            return 0.0


def ramp_break_strength(
    k: int, p: int, cycle_length: int = 4, num_cycles: int = 30,
    bit_shift: int = 1, max_epsilon: int = 8, break_threshold: float = 0.5,
) -> Dict:
    """
    Measure correlation between Newton convergence ratio and
    perturbation tolerance.

    DEPRECATED: The eps_crit metric is degenerate — for the shear
    perturbation matrix [[d, 1], [0, 1]], det = d, and the smallest
    ε giving an even determinant is always 1.  The genuinely
    interesting quantity is the α-sector flip (phase alignment),
    tracked via 'phase_alignment' and 'alpha_correlation' in the
    output dict.

    Use phase_alignment_experiment() instead.
    """
    import warnings
    warnings.warn(
        "ramp_break_strength.eps_crit is degenerate (always 1). "
        "Use phase_alignment_experiment() for the α-sector flip metric.",
        DeprecationWarning, stacklevel=2,
    )

    nc = NonAbelianCRTDual(k, p)
    results: List[Dict] = []

    for _ in range(num_cycles):
        mats = [
            [[np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)],
             [np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)]]
            for _ in range(cycle_length)
        ]

        inv = nc.invariants(mats)
        alpha0 = inv["alpha_det"]
        conv0 = nc.convergence_ratio_full(mats)

        eps_crit = None
        alpha_flipped = False
        for eps in range(1, max_epsilon + 1):
            d_new = (inv["det_mod2k"] + eps) & (nc.mod2 - 1)
            if d_new & 1 == 0:
                eps_crit = eps
            # Check α-sector flip
            dlog_new = two_adic_dlog(d_new | 1, k)  # force odd
            if dlog_new is not None:
                alpha_new = dlog_new[0]
                if alpha_new != alpha0:
                    alpha_flipped = True
                    if eps_crit is None:
                        eps_crit = eps

        results.append({
            "alpha0": alpha0,
            "conv0": conv0,
            "eps_crit": eps_crit if eps_crit else max_epsilon,
            "alpha_flipped": alpha_flipped,
        })

    conv_vals = [r["conv0"] for r in results]
    eps_vals = [r["eps_crit"] for r in results]
    n_flipped = sum(1 for r in results if r["alpha_flipped"])

    return {
        "mean_conv": float(np.mean(conv_vals)),
        "mean_eps_crit": float(np.mean(eps_vals)),
        "phase_alignment": n_flipped / len(results) if results else 0.0,
        "n_trials": num_cycles,
        "n_flipped": n_flipped,
        "deprecated": True,
    }


def phase_alignment_experiment(
    k: int, p: int, N_cycle: int = 4, n_cycles: int = 30,
) -> Dict[str, float]:
    """
    Test whether a single-bit perturbation flips the α-sector of the
    holonomy determinant.  This is the replacement metric for the
    degenerate eps_crit in ramp_break_strength.

    Returns dict with 'alignment' (fraction of trials where a flip
    occurred) and 'n_trials'.
    """
    nc = NonAbelianCRTDual(k, p)
    flips = 0

    for _ in range(n_cycles):
        mats = [
            [[np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)],
             [np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)]]
            for _ in range(N_cycle)
        ]

        inv = nc.invariants(mats)
        alpha0 = inv["alpha_det"]

        eps_mat = _perturbation_matrix(inv["det_mod2k"] + 1, nc.mod_full)
        perturbed = mats[:]
        perturbed[0] = _mat_mul(perturbed[0], eps_mat, nc.mod_full)
        inv_p = nc.invariants(perturbed)

        if inv_p["alpha_det"] != alpha0:
            flips += 1

    return {
        "alignment": flips / n_cycles if n_cycles > 0 else 0.0,
        "n_trials": n_cycles,
    }
