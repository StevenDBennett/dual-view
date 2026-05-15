"""
core.py
-------
Fast modular arithmetic in power-of-two rings via group-algebra duality.

Public names
~~~~~~~~~~~~
    modinv_newton   – Newton-lifted modular inverse mod 2^k
    two_adic_log5   – 2-adic logarithm of 5 (truncated to k bits)
    two_adic_dlog   – Full decomposition: n = 2^v · (-1)^alpha · 5^e
    DualNumber      – Coordinate triple (v, alpha, e) for n in Z/2^k
    TwoAdicProcessor– Arithmetic on DualNumbers
    run_all_tests   – Basic self-check
"""
from __future__ import annotations

from functools import lru_cache
from typing import Optional, Tuple
import math


# ── 8-bit discrete log LUT ─────────────────────────────────────────────────────
# Precomputed: for each odd a ≡ 1 (mod 4) in [0, 256), store e s.t. 5^e ≡ a (mod 256).
# 64 entries × 8 bytes = 512 bytes. Gives O(1) bootstrap to 6-bit precision.
_DLOG8_LUT: dict[int, int] = {
    pow(5, e, 256): e for e in range(64)
}


# ── private utilities ──────────────────────────────────────────────────────────

def _mask(k: int) -> int:
    return (1 << k) - 1


def _valuation(n: int) -> int:
    """2-adic valuation v_2(n)."""
    if n == 0:
        return float("inf")
    return (n & -n).bit_length() - 1


# ── Newton modular inverse ─────────────────────────────────────────────────────

def modinv_newton(a: int, k: int) -> int:
    """
    a^{-1} mod 2^k via quadratic Newton lifting.  Requires a odd.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    if a & 1 == 0:
        raise ValueError("a must be odd to be invertible mod 2^k")
    a &= _mask(k)
    x = a & 7
    bits = 3
    while bits < k:
        bits = min(bits * 2, k)
        x = (x * (2 - a * x)) & _mask(bits)
    return x & _mask(k)


# ── 2-adic log of 5 ────────────────────────────────────────────────────────────

@lru_cache(maxsize=256)
def two_adic_log5(k: int) -> int:
    """
    Compute the 2-adic logarithm of 5, truncated to k bits.

    This is the unique L in Z_2 satisfying exp_2(L) = 5 in the
    2-adic sense.  Concretely it is computed via the 2-adic log series
    applied to (5-1).

    Primary use: derivative of the map e ↦ 5^e is 5^e · (L >> 2).
    Cached per k.
    """
    mask = _mask(k)
    result = 0
    n = 1
    while True:
        v = _valuation(n)
        exp = 2 * n - v
        if exp >= k:
            break
        odd_part = n >> v
        inv_odd = modinv_newton(odd_part, k)
        term = ((1 << exp) * inv_odd) & mask
        if n % 2 == 0:
            term = (-term) & mask
        result = (result + term) & mask
        n += 1
    return result


# ── discrete-log helpers (private) ─────────────────────────────────────────────

def _dlog_bootstrap(a: int, k: int) -> int:
    """
    Bit-by-bit dlog for a ≡ 1 (mod 8).
    Returns e with 5^e ≡ a (mod 2^k).  O(k) steps; used to seed Newton.
    """
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


def _dlog_newton(a: int, k: int, L: Optional[int] = None) -> int:
    """
    Newton-lifted dlog: 5^e ≡ a (mod 2^k).  Requires a ≡ 1 (mod 4).

    Bootstrap uses an 8-bit LUT for k ≤ 34 (covers 6-bit precision),
    otherwise uses bit-by-bit to k/2.  Per the Mersenne boostrap
    optimality analysis, k/2 saves ⌈log₂(√k)⌉−1 Newton steps vs √k.
    """
    if k <= 2:
        return 0
    if a & 3 != 1:
        raise ValueError("_dlog_newton requires a ≡ 1 (mod 4)")
    mask_full = _mask(k)
    a &= mask_full
    if L is None:
        L = two_adic_log5(k)
    L_unit = L >> 2

    if k <= 34:
        a_b = a & 0xFF
        e_raw = _DLOG8_LUT.get(a_b, 0)
        eprec = min(6, k - 2)
        e = e_raw & _mask(eprec)
    else:
        bootstrap_k = max(4, k // 2 + 2)
        eprec = min(bootstrap_k - 2, k - 2)
        e = _dlog_bootstrap(a, bootstrap_k) & _mask(eprec)

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = _mask(bits)
        emask = _mask(new_eprec)

        pow5e = pow(5, e, 1 << bits)
        f = (pow5e - a) & mask
        df_unit = (pow5e * L_unit) & emask
        df_inv = modinv_newton(df_unit, new_eprec)
        delta = ((f >> 2) * df_inv) & emask
        e = (e - delta) & emask
        eprec = new_eprec

    return e


# ── public discrete-log API ────────────────────────────────────────────────────

def two_adic_dlog(
    a: int, k: int, L: Optional[int] = None
) -> Optional[Tuple[int, int]]:
    """
    Decompose the odd part of a as (-1)^alpha · 5^e (mod 2^k).

    Parameters
    ----------
    a : int   — integer to decompose (any residue mod 2^k)
    k : int   — bit precision (≥ 3)
    L : int   — precomputed two_adic_log5(k); computed if not supplied

    Returns
    -------
    (alpha, e)  with alpha ∈ {0,1} and e ∈ [0, 2^(k-2)), or
    None        if a is even.
    """
    if a & 1 == 0:
        return None
    mask = _mask(k)
    a &= mask
    alpha = (a >> 1) & 1
    if alpha:
        a = (-a) & mask
    e = _dlog_newton(a, k, L)
    return alpha, e


# ── DualNumber ─────────────────────────────────────────────────────────────────

class DualNumber:
    """
    Represent n ∈ Z/2^k in dual coordinates (v, alpha, e) where

        n  ≡  2^v · (-1)^alpha · 5^e   (mod 2^k)

    v     = v_2(n)  (2-adic valuation)
    alpha ∈ {0, 1}  (sign of the odd unit part)
    e     ∈ [0, 2^(k-2))  (discrete log base 5 of the odd unit part)

    Zero is a distinguished element with v = ∞.

    Construction
    ------------
    From an integer:       DualNumber(42, k=16)
    From coordinates:      DualNumber.from_coords(v, alpha, e, k)
    """

    __slots__ = ("k", "mask", "_ord", "v", "alpha", "e", "_n", "is_zero")

    def __init__(self, n: int, k: int = 64):
        if k < 3:
            raise ValueError("k must be ≥ 3")
        self.k = k
        self.mask = _mask(k)
        self._ord = 1 << (k - 2)

        n_mod = n & self.mask
        if n_mod == 0:
            self.is_zero = True
            self.v = float("inf")
            self.alpha = 0
            self.e = 0
            self._n = 0
            return

        self.is_zero = False
        self.v = _valuation(n_mod)
        odd = (n_mod >> self.v) & self.mask

        result = two_adic_dlog(odd, k)
        if result is None:
            raise ValueError(f"Cannot decompose {n} mod 2^{k}: odd part is even")
        self.alpha, self.e = result
        self._n = n_mod

    @classmethod
    def from_coords(cls, v: int, alpha: int, e: int, k: int = 64) -> "DualNumber":
        """Build a DualNumber directly from coordinates."""
        obj = cls.__new__(cls)
        obj.k = k
        obj.mask = _mask(k)
        obj._ord = 1 << (k - 2)
        if v >= k:
            obj.is_zero = True
            obj.v = float("inf")
            obj.alpha = 0
            obj.e = 0
            obj._n = 0
            return obj
        obj.is_zero = False
        obj.v = v
        obj.alpha = alpha & 1
        obj.e = e % obj._ord
        obj._n = obj._to_int()
        return obj

    def _to_int(self) -> int:
        if self.is_zero:
            return 0
        odd = pow(5, self.e, 1 << self.k)
        if self.alpha:
            odd = (-odd) & self.mask
        return (odd << self.v) & self.mask

    def verify(self) -> bool:
        """Round-trip check: coordinates → integer matches stored integer."""
        if self.is_zero:
            return self._n == 0
        return self._to_int() == self._n

    @property
    def value(self) -> int:
        """The integer n mod 2^k."""
        return self._n

    def coords(self) -> Tuple[int, int, int]:
        """Return (v, alpha, e)."""
        return (self.v, self.alpha, self.e)

    def __repr__(self) -> str:
        if self.is_zero:
            return f"DualNumber(0, k={self.k})"
        sign = "-" if self.alpha else "+"
        return (
            f"DualNumber({self._n}, k={self.k})"
            f"  =  2^{self.v} · {sign}5^{self.e}"
        )


# ── TwoAdicProcessor ──────────────────────────────────────────────────────────

class TwoAdicProcessor:
    """
    Arithmetic on DualNumbers in Z/2^k.

    Multiplication, inversion, and exponentiation operate directly in
    coordinate space, making the group structure explicit.
    """

    def __init__(self, k: int = 64):
        if k < 3:
            raise ValueError("k must be ≥ 3")
        self.k = k
        self.mask = _mask(k)
        self._ord = 1 << (k - 2)
        self.L = two_adic_log5(k)

    def _check(self, *args: DualNumber) -> None:
        for a in args:
            if a.k != self.k:
                raise ValueError(
                    f"DualNumber has k={a.k} but processor has k={self.k}"
                )

    def mul(self, a: DualNumber, b: DualNumber) -> DualNumber:
        """Multiply two elements: coordinates add componentwise."""
        self._check(a, b)
        if a.is_zero or b.is_zero:
            return DualNumber(0, self.k)
        v = a.v + b.v
        if v >= self.k:
            return DualNumber(0, self.k)
        result = DualNumber.from_coords(
            v=v,
            alpha=a.alpha ^ b.alpha,
            e=(a.e + b.e) % self._ord,
            k=self.k,
        )
        assert result.verify()
        return result

    def inv(self, a: DualNumber) -> DualNumber:
        """Invert a unit (v = 0).  Raises ValueError if a is not a unit."""
        self._check(a)
        if a.is_zero or a.v != 0:
            raise ValueError("Only units (v=0) are invertible")
        result = DualNumber.from_coords(
            v=0,
            alpha=a.alpha,
            e=(-a.e) % self._ord,
            k=self.k,
        )
        assert result.verify()
        return result

    def pow(self, a: DualNumber, n: int) -> DualNumber:
        """Raise a to an integer power.  Negative exponents require a to be a unit."""
        self._check(a)
        if a.is_zero:
            return DualNumber(0, self.k)
        if n < 0:
            return self.pow(self.inv(a), -n)
        v = a.v * n
        if v >= self.k:
            return DualNumber(0, self.k)
        result = DualNumber.from_coords(
            v=v,
            alpha=a.alpha if (n % 2 == 1) else 0,
            e=(a.e * n) % self._ord,
            k=self.k,
        )
        assert result.verify()
        return result

    def dlog(self, a: DualNumber) -> Tuple[int, int, int]:
        """Return full coordinate triple (v, alpha, e)."""
        self._check(a)
        return a.coords()


# ── self-test ─────────────────────────────────────────────────────────────────

def run_all_tests(k: int = 16, verbose: bool = True) -> None:
    """
    Basic sanity checks for the core arithmetic.

    Tests round-trip encoding, multiplication, inversion, and powering.
    Raises AssertionError on any failure.
    """
    proc = TwoAdicProcessor(k)
    mask = _mask(k)
    failures: list[str] = []

    samples = [1, 3, 5, 7, 9, 15, 255, (1 << (k // 2)) + 1, mask - 2]
    for n in samples:
        d = DualNumber(n, k)
        if not d.verify():
            failures.append(f"Round-trip failed: n={n}")

    z = DualNumber(0, k)
    if not z.is_zero:
        failures.append("Zero not recognised as is_zero")

    d3 = DualNumber(3, k)
    d5 = DualNumber(5, k)
    d15 = proc.mul(d3, d5)
    if d15.value != (15 & mask):
        failures.append(f"mul: 3×5 → {d15.value}, expected 15")

    inv3 = proc.inv(d3)
    if proc.mul(d3, inv3).value != 1:
        failures.append("inv: 3 × 3⁻¹ ≠ 1")

    if proc.pow(d3, 4).value != (81 & mask):
        failures.append(f"pow: 3^4 → {proc.pow(d3, 4).value}, expected {81 & mask}")

    if failures:
        for msg in failures:
            print(f"FAIL  {msg}")
        raise AssertionError(f"{len(failures)} test(s) failed at k={k}")

    if verbose:
        print(f"run_all_tests: all checks passed  (k={k})")
