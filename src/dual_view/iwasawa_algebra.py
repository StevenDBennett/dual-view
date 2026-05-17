"""
iwasawa_algebra.py
------------------
Iwasawa algebra Z_2[[G]] and profinite filtered modules.

Provides the Iwasawa algebra Z_2[[G]] as power series in the topological
generator (1 - gamma), together with the classification of shift-covariant
differential operators via the augmentation ideal.

Classes
-------
IwasawaElement   -- An element of Z_2[[G]] as a power series in (1-gamma).
IwasawaAlgebra   -- Factory and operations for the Iwasawa algebra.
ProModule        -- A profinite Z_2-module with valuation filtration.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


# ============================================================================
# IWASAWA ALGEBRA Z_2[[G]]
# ============================================================================

@dataclass
class IwasawaElement:
    """
    An element of the Iwasawa algebra Z_2[[G]] represented as a power series
    in the topological generator gamma:
        mu = sum_{n>=0} c_n (1-gamma)^n

    The augmentation ideal I = (1-gamma) is principal.
    """
    coeffs: List[int]
    precision: int

    def __post_init__(self):
        self.original_degree = len(self.coeffs)

        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()
        if not self.coeffs:
            self.coeffs = [0]

        self.truncation_degree = len(self.coeffs)
        self.truncation_error = 0

    @classmethod
    def from_generator(cls, precision: int = 32) -> IwasawaElement:
        """The canonical generator 1-gamma of the augmentation ideal."""
        return cls(coeffs=[0, 1], precision=precision)

    @classmethod
    def unit(cls, precision: int = 32) -> IwasawaElement:
        """The multiplicative identity."""
        return cls(coeffs=[1], precision=precision)

    @classmethod
    def zero(cls, precision: int = 32) -> IwasawaElement:
        return cls(coeffs=[0], precision=precision)

    def valuation(self) -> int | float:
        """2-adic valuation: the smallest n with c_n odd."""
        for i, c in enumerate(self.coeffs):
            if c % 2 != 0:
                return i
        return float('inf')

    def is_unit(self) -> bool:
        """An element is a unit iff its constant term is odd."""
        return len(self.coeffs) > 0 and self.coeffs[0] % 2 == 1

    def is_generator_of_aug_ideal(self) -> bool:
        """Check if this element generates the augmentation ideal."""
        if self.valuation() == 0:
            return False
        if len(self.coeffs) < 2:
            return False
        return self.coeffs[1] % 2 == 1

    def __add__(self, other: IwasawaElement) -> IwasawaElement:
        assert self.precision == other.precision
        max_len = max(len(self.coeffs), len(other.coeffs))
        mod = 1 << self.precision
        new_coeffs = []
        for i in range(max_len):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = other.coeffs[i] if i < len(other.coeffs) else 0
            new_coeffs.append((a + b) % mod)
        return IwasawaElement(new_coeffs, self.precision)

    def __mul__(self, other: IwasawaElement) -> IwasawaElement:
        """Multiplication in Z_2[[G]] via Cauchy product."""
        assert self.precision == other.precision
        mod = 1 << self.precision
        result_len = min(len(self.coeffs) + len(other.coeffs) - 1, self.precision)
        new_coeffs = [0] * result_len

        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                if i + j < result_len:
                    new_coeffs[i + j] = (new_coeffs[i + j] + a * b) % mod

        return IwasawaElement(new_coeffs, self.precision)

    def truncation_status(self) -> str:
        """Report truncation metadata."""
        return (f"original_degree={self.original_degree}, "
                f"truncation_degree={self.truncation_degree}, "
                f"precision=2^{self.precision}")

    def __repr__(self) -> str:
        terms = []
        for i, c in enumerate(self.coeffs):
            if c == 0:
                continue
            if i == 0:
                terms.append(f"{c}")
            elif i == 1:
                terms.append(f"{c}(1-gamma)")
            else:
                terms.append(f"{c}(1-gamma)^{i}")
        base = " + ".join(terms) if terms else "0"
        return f"{base}  [{self.truncation_status()}]"


class IwasawaAlgebra:
    """Factory and operations for the Iwasawa algebra."""

    @staticmethod
    def aug_ideal_generator(precision: int = 32) -> IwasawaElement:
        return IwasawaElement.from_generator(precision)

    @staticmethod
    def classify_dirac_operator(mu: IwasawaElement) -> Dict:
        """Classify an element as a shift-covariant differential operator."""
        val = mu.valuation()
        result = {
            'is_valid': False,
            'is_unit_multiple': False,
            'generator_form': None,
        }

        if val == 0:
            return result
        if not mu.is_generator_of_aug_ideal():
            return result

        result['is_valid'] = True
        result['is_unit_multiple'] = True
        result['generator_form'] = mu
        return result


# ============================================================================
# PROFINITE FILTERED MODULES
# ============================================================================

@dataclass
class ProModule:
    """A profinite Z_2-module with valuation filtration."""
    name: str
    dimension: int
    precision: int

    def __post_init__(self):
        self.filtration: Dict[int, List[int]] = {}
        self.graded: Dict[int, List[int]] = {}

    def truncate(self, k: int) -> ProModule:
        """Apply the truncation functor T_k."""
        return ProModule(f"{self.name}^({k})", self.dimension, k)

    def valuation(self, v: List[int]) -> int | float:
        """Compute the valuation of an element (largest n where all
        coefficients are divisible by 2^n, bounded by precision)."""
        if not v:
            return float('inf')
        n = 0
        while n < self.precision and all(c % (1 << (n + 1)) == 0 for c in v):
            n += 1
        return n
