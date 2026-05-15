"""
mersenne.py
-----------
Mersenne Ghost Theorem and bootstrap optimality analysis.

Mersenne Ghost Theorem
    For w = 2^n - 1 (a Mersenne number), the 2-adic dual coordinates
    are always α = 1 (ghost sector) and e_true = 2^(n-2) within a
    stable window.  Consequently v₂(e_true) = n-2 and the quantisation
    cliff is at k* = n+2.

Bootstrap Optimality
    In the Viglietta discrete-log algorithm, the bootstrap phase
    starting precision eprec₀ = k/2 is optimal (saves log₂(√k) - 1
    Newton steps compared to the √k heuristic currently used).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
from functools import lru_cache
import math

from .core import (
    _mask, _valuation, modinv_newton,
    two_adic_log5, two_adic_dlog, DualNumber,
)


# ── Mersenne coordinates ────────────────────────────────────────────────────

def mersenne_coordinates(n: int, k: int) -> Optional[Tuple[int, int, int]]:
    """
    Return (α, e_true, v₂(e_true)) for Mersenne weight w = 2^n - 1.

    Requires n ≥ 3 and k ≥ n + 2.
    """
    if n < 3:
        raise ValueError("n must be ≥ 3")
    w = (1 << n) - 1
    result = two_adic_dlog(w, k)
    if result is None:
        return None
    alpha, e_true = result
    v2_e = _valuation(e_true) if e_true != 0 else k
    return alpha, e_true, v2_e


def verify_core_identity(n_max: int = 12) -> Dict[int, bool]:
    """
    Verify the core identity: 5^(2^(n-2)) ≡ 1 - 2^n (mod 2^(n+1)).

    This identity holds for all n ≥ 3 and is the foundation of the
    Mersenne Ghost Theorem.
    """
    results: Dict[int, bool] = {}
    for n in range(3, n_max + 1):
        lhs = pow(5, 1 << (n - 2), 1 << (n + 1))
        rhs = (1 - (1 << n)) & _mask(n + 1)
        results[n] = lhs == rhs
    return results


def mersenne_cliff_table(n_max: int = 12) -> List[Dict]:
    """
    Find the cliff precision k* for each Mersenne weight 2^n - 1.

    Returns list of dicts with n, k*, α, e_true, v₂(e_true), and
    the formula prediction k_pred = n + 2.
    """
    rows = []
    for n in range(3, n_max + 1):
        stable_k = None
        for k in range(n + 1, n + 10):
            result = two_adic_dlog((1 << n) - 1, k)
            if result is not None:
                _, e_true = result
                if e_true == (1 << (n - 2)):
                    stable_k = k
                else:
                    break
        if stable_k is not None:
            _, e_true, v2_e = mersenne_coordinates(n, stable_k + 1) or (0, 0, 0)
            rows.append({
                "n": n,
                "k*": stable_k,
                "k_pred": n + 2,
                "alpha": 1,
                "e_true": 1 << (n - 2),
                "v2_e": v2_e,
            })
    return rows


# ── Bootstrap optimality ────────────────────────────────────────────────────

def bootstrap_cost(eprec0: int, k: int) -> int:
    """
    Total bits processed by the Viglietta dlog algorithm.

    Cost = eprec0 (bootstrap) + Σ_{i} 2·eprec_i (Newton steps)
    where eprec_i doubles each iteration.
    """
    cost = eprec0
    eprec = eprec0
    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        cost += new_eprec  # each Newton step processes ~new_eprec bits
        eprec = new_eprec
    return cost


def optimal_bootstrap(k_values: List[int] = None) -> Dict[int, int]:
    """
    Find the optimal eprec₀ for each k by minimising total cost.

    Returns dict mapping k → optimal eprec₀ (consistently near k/2).
    """
    if k_values is None:
        k_values = list(range(8, 65, 4))

    results: Dict[int, int] = {}
    for k in k_values:
        best_cost = float("inf")
        best_eprec = 0
        for eprec0 in range(2, k - 2):
            cost = bootstrap_cost(eprec0, k)
            if cost < best_cost:
                best_cost = cost
                best_eprec = eprec0
        results[k] = best_eprec
    return results


def compare_bootstrap_strategies(k_values: List[int] = None) -> Dict[int, Dict]:
    """
    Compare sqrt(k) heuristic vs k/2 optimal vs LUT b=8.

    Returns dict mapping k → {sqrt_cost, half_cost, lut_cost}.
    """
    if k_values is None:
        k_values = [16, 24, 32, 48, 64]

    results: Dict[int, Dict] = {}
    for k in k_values:
        sqrt_eprec = max(4, int(math.isqrt(k)) + 2)
        half_eprec = max(4, k // 2 + 2)
        results[k] = {
            "sqrt_eprec": sqrt_eprec,
            "sqrt_cost": bootstrap_cost(sqrt_eprec - 2, k),
            "half_eprec": half_eprec,
            "half_cost": bootstrap_cost(half_eprec - 2, k),
            "lut_cost": bootstrap_cost(6, k),  # 8-bit LUT gives 6 bits
        }
    return results


# ── LUT-based dlog ──────────────────────────────────────────────────────────

@lru_cache(maxsize=32)
def _build_lut(b: int = 8) -> Dict[int, int]:
    """
    Build a lookup table for dlog modulo 2^b.

    Returns dict mapping a → e such that 5^e ≡ a (mod 2^b),
    for all odd a ≡ 1 (mod 4).
    """
    mod = 1 << b
    lut: Dict[int, int] = {}
    for e in range(1 << (b - 2)):
        a = pow(5, e, mod)
        lut[a] = e
    return lut


def dlog_with_lut(a: int, k: int, b: int = 8) -> int:
    """
    Discrete logarithm using a pre-computed LUT for bootstrap.

    Eliminates the O(k) bootstrap phase — the LUT gives exactly
    6 bits of precision at b=8, and each Newton step doubles.
    """
    if k <= b - 2:
        return _dlog_bit_by_bit(a, k)
    if a & 3 != 1:
        raise ValueError("dlog_with_lut requires a ≡ 1 (mod 4)")

    lut = _build_lut(b)
    a_b = a & ((1 << b) - 1)
    e = lut.get(a_b, 0) & _mask(b - 2)
    eprec = b - 2

    L = two_adic_log5(k) >> 2
    mask_full = _mask(k)
    a &= mask_full

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = _mask(bits)
        emask = _mask(new_eprec)
        pow5e = pow(5, e, 1 << bits)
        f = (pow5e - a) & mask
        df_unit = (pow5e * L) & emask
        df_inv = modinv_newton(df_unit, new_eprec)
        delta = ((f >> 2) * df_inv) & emask
        e = (e - delta) & emask
        eprec = new_eprec

    return e


def _dlog_bit_by_bit(a: int, k: int) -> int:
    """Simple O(k) bit-by-bit discrete log (fallback)."""
    if k <= 2:
        return 0
    mask = _mask(k)
    a &= mask
    e, pow5 = (0, 1) if (a & 7) == 1 else (1, 5)
    mult = 25 & mask
    for n in range(3, k):
        if ((a >> n) & 1) != ((pow5 >> n) & 1):
            e |= (1 << (n - 2))
            pow5 = (pow5 * mult) & mask
        mult = (mult * mult) & mask
    return e


def verify_lut_dlog(k: int, b: int = 8, n_trials: int = 100) -> bool:
    """
    Verify that LUT-based dlog matches the standard dlog.
    """
    for _ in range(n_trials):
        import random
        a = random.randrange(1, 1 << k)
        if a & 1 == 0:
            continue
        a_mod = (a if (a & 3) == 1 else (-a) & _mask(k))
        e1 = dlog_with_lut(a_mod, k, b)
        e2 = _dlog_bit_by_bit(a_mod, k)
        if e1 != e2:
            return False
    return True
