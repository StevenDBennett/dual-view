"""
operators.py
------------
Symbolic operator algebra on the function space Z/N → Z/2^k.

Provides an environment of linear operators (identity, shift, forward
difference, multiplication, summation) as Python objects, together
with a Connes-style spectral triple and a Newton projector that lifts
discrete-log seeds to full precision.

Closed-form identities verified in isometry.py:
    avg^2 = N·avg
    D·avg = 0
    avg·D = 0

PERFORMANCE NOTE: The operator algebra is implemented via Python
closures and is O(N) per application (N = 2^{k-2}).  For k ≥ 9
(N ≥ 128) this becomes a pedagogical / decorative interface rather
than a practical computational tool.  The underlying theory is
independent of the implementation — for production use at large k,
the discrete-log Newton projector in core.py is the efficient path.

Bug fixes applied (from the original 2-Adic-Newton-Dynamics codebase):
    1. OperatorContext.g attribute was missing — added as a property.
    2. NewtonProjector._step: was inverting an even number (4*f*Lg).
       Fixed to invert df_unit = 5^e * L (always odd, since the
       derivative of the map e -> 5^e is 5^e * log(5) which is odd
       for all e).  The factor of 4 -> f>>2 shift is moved to the
       numerator.
    3. R operator (now called avg): was O(N^2) due to nested
       summation closures.  Fixed to O(N) direct-sum lambda, BUT
       still O(N) per application, so the operator algebra as a
       whole is impractical for k ≥ 10.
    4. Closures in __add__/__sub__: used late-binding — fixed with
       default-argument capture (s=self, o=other).
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple


class _Operator:
    """
    Symbolic linear operator on the space of functions Z/N → Z/2^k.

    Operators support addition, subtraction, composition (via *), and
    integer powers.  They are callable:  (op)(f, e) applies the
    operator to function f at exponent e.
    """

    def __init__(self, label: str, action: Callable, mask: int = 0) -> None:
        self._label = label
        self._action = action
        self._mask = mask

    def __call__(self, f: Callable, e: int) -> int:
        return self._action(f, e)

    def __add__(self, other: _Operator) -> _Operator:
        m = self._mask or other._mask or (1 << 64) - 1
        def action(f, e, s=self, o=other):
            return (s(f, e) + o(f, e)) & m
        return _Operator(f"({self._label}+{other._label})", action, m)

    def __sub__(self, other: _Operator) -> _Operator:
        m = self._mask or other._mask or (1 << 64) - 1
        def action(f, e, s=self, o=other):
            return (s(f, e) - o(f, e)) & m
        return _Operator(f"({self._label}-{other._label})", action, m)

    def __mul__(self, other: _Operator) -> _Operator:
        m = self._mask or other._mask or (1 << 64) - 1
        def action(f, e, s=self, o=other):
            return s(lambda e: o(f, e), e)
        return _Operator(f"({self._label}∘{other._label})", action, m)

    def __pow__(self, n: int) -> _Operator:
        result = self
        for _ in range(n - 1):
            result = result * self
        return result

    def __repr__(self) -> str:
        return f"Op[{self._label}]"


class Identity(_Operator):
    """Identity operator I: f ↦ f."""

    def __init__(self, mask: int = 0) -> None:
        super().__init__("I", lambda f, e: f(e), mask)


class Shift(_Operator):
    """Shift operator S: (Sf)(e) = f(e + 1)."""

    def __init__(self, N: int, mask: int = 0) -> None:
        self.N = N
        super().__init__("S", lambda f, e: f((e + 1) % self.N), mask)


class Multiplication(_Operator):
    """Multiplication operator M_h: (M_h f)(e) = h(e)·f(e)."""

    def __init__(self, h: Callable, label: str = "M", mask: int = 0) -> None:
        super().__init__(label, lambda f, e: (h(e) * f(e)) & mask, mask)


class OperatorContext:
    """
    Create a complete operator environment for precision k and generator g.

    Provides:
        I      – identity
        S      – shift
        diff   – forward difference I - S
        avg    – direct-sum average (O(N) per application)
        M(h)   – multiplication by h
        g      – the generator (fix for missing attribute)

    NOTE: All operators are implemented via Python closures and run in
    O(N) time per application (N = 2^{k-2}).  For k ≥ 9 (N ≥ 128) this
    module becomes a pedagogical / decorative interface.  For efficient
    computation at large k use the discrete-log tools in core.py.
    """

    def __init__(self, k: int, g: int) -> None:
        self.k = k
        self.mask = (1 << k) - 1
        self.N = 1 << (k - 2)
        self._g = g

        self.I = Identity(self.mask)
        self.S = Shift(self.N, self.mask)
        self.diff = self.I - self.S
        self._g_pow = [pow(g, e, 1 << k) for e in range(self.N)]

        def avg_action(f, e, ctx=self):
            total = 0
            for t in range(ctx.N):
                total += f(t)
            return total & ctx.mask
        self.avg = _Operator("avg", avg_action, self.mask)

    @property
    def g(self) -> int:
        return self._g

    def multiply_by(self, h: Callable, name: str = "M") -> _Operator:
        """Create a Multiplication operator by function h."""
        return Multiplication(h, name, self.mask)


class SpectralTriple:
    """
    Connes-style spectral triple (A, H, D) on the exponent domain.

    H is the Hilbert space of functions Z/N → Z/2^k (discrete,
    finite-dimensional).  D = diff is the Dirac operator.
    """

    def __init__(self, ctx: OperatorContext) -> None:
        self.ctx = ctx
        self.dirac = ctx.diff

    def one_form(self, h: Callable) -> _Operator:
        """ω = [D, M_h] — a Connes one-form."""
        return self.dirac * self.ctx.multiply_by(h) - self.ctx.multiply_by(h) * self.dirac

    def curvature(self, h1: Callable, h2: Callable) -> _Operator:
        """Ω = [ω_1, ω_2] — curvature of the one-form connection."""
        w1 = self.one_form(h1)
        w2 = self.one_form(h2)
        return w1 * w2 - w2 * w1

    def gauge_transform(self, potential: _Operator) -> _Operator:
        """D' = D + ω — gauge-transformed Dirac operator."""
        return self.dirac + potential

    def trace(self, f: Callable) -> int:
        """Tr(f) = Σ_e f(e).  Manually computed for finite dimension."""
        total = 0
        for e in range(self.ctx.N):
            total += f(e)
        return total & self.ctx.mask


class NewtonProjector:
    """
    Newton projection operator for solving 5^e ≡ a (mod 2^k).

    Starting from a seed e₀, repeated application of the Newton step
    converges quadratically to the true root e_true.
    """

    def __init__(self, ctx: OperatorContext, a: int, steps: int = 16) -> None:
        self.ctx = ctx
        self.a = a & ctx.mask
        self.L = _compute_log5(ctx.k) >> 2
        self.steps = steps

    def _step(self, e: int) -> int:
        mask = self.ctx.mask
        g = pow(5, e, 1 << self.ctx.k)
        f = (g - self.a) & mask
        df = (g * self.L) & mask
        e_mask = (1 << (self.ctx.k - 2)) - 1
        df_inv = _modinv_newton(df, self.ctx.k - 2)
        delta = ((f >> 2) * df_inv) & e_mask
        return (e - delta) & e_mask

    def project_point(self, e_guess: int) -> int:
        """Apply Newton refinement for self.steps iterations."""
        e = e_guess
        for _ in range(self.steps):
            e = self._step(e)
        return e


# Internal helpers (avoid circular imports with core.py on __init__)

def _compute_log5(k: int) -> int:
    from .core import two_adic_log5
    return two_adic_log5(k)


def _modinv_newton(a: int, k: int) -> int:
    from .core import modinv_newton
    return modinv_newton(a, k)
