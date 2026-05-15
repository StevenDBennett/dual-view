"""
crt.py
------
Chinese Remainder Theorem dual system:  Z/(2^k · p)Z  for odd prime p.

Combines the 2-adic dual-view decomposition with a modular (prime-field)
residue, allowing joint analysis of weights in the product ring.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import math
import numpy as np

from .core import _mask, _valuation, modinv_newton, two_adic_log5, two_adic_dlog, DualNumber


def _primitive_root(p: int) -> Optional[int]:
    """Find the smallest primitive root modulo prime p."""
    if p == 2:
        return 1
    phi = p - 1
    factors = set()
    n = phi
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.add(n)

    for g in range(2, p):
        ok = True
        for q in factors:
            if pow(g, phi // q, p) == 1:
                ok = False
                break
        if ok:
            return g
    return None


def _prime_dlog(a: int, p: int, g_p: int) -> int:
    """Brute-force discrete logarithm modulo prime p."""
    cur = 1
    for e in range(p):
        if cur == a:
            return e
        cur = (cur * g_p) % p
    return 0


class CRTDualNumber:
    """
    Element of Z/(2^k · p)Z with dual-coordinate representation.

    Attributes
    ----------
    component_2 : DualNumber — the 2-adic part
    residue_p  : int — residue modulo p
    dlog_p     : int — discrete-log of residue_p modulo p (base g_p)
    """

    def __init__(self, n: int, k: int, p: int, g_p: int) -> None:
        self.k = k
        self.p = p
        self.g_p = g_p
        self.mod_full = (1 << k) * p
        self.component_2 = DualNumber(n % (1 << k), k)
        self.residue_p = n % p
        self.dlog_p = _prime_dlog(self.residue_p, p, g_p) if self.residue_p != 0 else -1

    def verify(self) -> bool:
        """CRT round-trip check."""
        return self.component_2.verify()

    def __repr__(self) -> str:
        c2 = self.component_2
        if c2.is_zero:
            return f"CRTDualNumber(0, k={self.k}, p={self.p})"
        sign = "-" if c2.alpha else "+"
        return (
            f"CRTDualNumber(k={self.k}, p={self.p})"
            f"  =  2^{c2.v} · {sign}5^{c2.e}  (mod 2^{self.k})"
            f"  ×  g_p^{self.dlog_p}  (mod {self.p})"
        )


class CRTDualProcessor:
    """
    Arithmetic processor for CRTDualNumbers over Z/(2^k · p)Z.
    """

    def __init__(self, k: int, p: int, g_p: Optional[int] = None) -> None:
        self.k = k
        self.p = p
        self.mod2 = 1 << k
        self.mod_full = self.mod2 * p
        if g_p is None:
            g_p = _primitive_root(p)
            if g_p is None:
                raise ValueError(f"No primitive root found for p={p}")
        self.g_p = g_p

    def crt_reconstruct(self, r2: int, rp: int) -> int:
        """
        CRT: combine residue r2 (mod 2^k) and rp (mod p) into
        a single integer mod 2^k·p.
        """
        inv_p = modinv_newton(self.p, self.k)
        t = ((r2 - rp) * inv_p) & (self.mod2 - 1)
        return (t * self.p + rp) % self.mod_full

    def mul(self, A: CRTDualNumber, B: CRTDualNumber) -> CRTDualNumber:
        """Multiply two CRTDualNumbers directly."""
        n2 = (A.component_2.value * B.component_2.value) % self.mod2
        np_val = (A.residue_p * B.residue_p) % self.p
        n_full = self.crt_reconstruct(n2, np_val)
        return CRTDualNumber(n_full, self.k, self.p, self.g_p)

    def product(self, weights: List[int]) -> CRTDualNumber:
        """Product of raw integers in CRT space."""
        nums = [CRTDualNumber(w, self.k, self.p, self.g_p) for w in weights]
        return self.cycle_product(nums)

    def cycle_product(self, numbers: List[CRTDualNumber]) -> CRTDualNumber:
        """Product of a list of CRTDualNumbers."""
        if not numbers:
            return CRTDualNumber(1, self.k, self.p, self.g_p)
        prod = numbers[0]
        for num in numbers[1:]:
            prod = self.mul(prod, num)
        return prod

    def convergence_ratio_2adic(self, P: CRTDualNumber) -> float:
        """Convergence ratio of the 2-adic component of P."""
        from .basin import BasinExplorer

        a = P.component_2.value
        if a == 0 or (a & 1) == 0:
            return 0.0
        try:
            explorer = BasinExplorer(self.k, 5, a)
            portrait = explorer.portrait()
            n_converged = len(portrait['converged'])
            total = n_converged + len(portrait['cycle'])
            return n_converged / total if total > 0 else 0.0
        except Exception:
            return 0.0

    def __repr__(self) -> str:
        return f"CRTDualProcessor(k={self.k}, p={self.p}, g_p={self.g_p})"


def combined_stability(
    k: int, p: int, num_cycles: int = 50, cycle_length: int = 4
) -> Dict[str, float]:
    """
    Randomised test: correlation between 2-adic convergence ratio
    and v_2 of change under single-bit weight-flip perturbation.

    Returns dict with pearson_r, n_samples, and mean_ratio.
    """
    proc = CRTDualProcessor(k, p)
    ratios_orig: List[float] = []
    delta_v2s: List[float] = []

    for _ in range(num_cycles):
        weights = [np.random.randint(0, proc.mod_full) for _ in range(cycle_length)]
        P = proc.cycle_product([CRTDualNumber(w, k, p, proc.g_p) for w in weights])
        r_orig = proc.convergence_ratio_2adic(P)

        flip_idx = np.random.randint(0, cycle_length)
        w_flipped = weights[flip_idx] ^ (1 << np.random.randint(0, min(k, 16)))
        weights_pert = weights[:]
        weights_pert[flip_idx] = w_flipped
        P_pert = proc.cycle_product(
            [CRTDualNumber(w, k, p, proc.g_p) for w in weights_pert]
        )
        r_pert = proc.convergence_ratio_2adic(P_pert)

        delta = abs(P.component_2.value - P_pert.component_2.value)
        v2_delta = _valuation(delta) if delta > 0 else 0

        ratios_orig.append(r_orig)
        delta_v2s.append(float(v2_delta))

    ratios_arr = np.array(ratios_orig)
    deltas_arr = np.array(delta_v2s)
    std_r = ratios_arr.std()
    std_d = deltas_arr.std()

    if std_r > 0 and std_d > 0:
        pearson_r = float(np.corrcoef(ratios_arr, deltas_arr)[0, 1])
    else:
        pearson_r = 0.0

    return {
        "pearson_r": pearson_r,
        "n_samples": num_cycles,
        "mean_ratio": float(ratios_arr.mean()),
    }
