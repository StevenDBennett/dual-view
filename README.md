# dual-view

**dual-view: a unified mathematical framework for 2-adic number systems, p-adic Newton dynamics, noncommutative spectral geometry, and quantum butterfly compilation.**

This package provides a unified mathematical framework spanning 2-adic number theory, p-adic dynamics, spectral triples, gauge theory, and quantum circuit synthesis. The core insight is that every odd
integer modulo `2^k` decomposes uniquely as:

```
n = 2^v · (-1)^α · 5^e   (mod 2^k)
```

where `v` is the 2-adic valuation, `α` is the sign sector, and `e` is the
discrete logarithm base 5 — the **dual-view coordinates**.

This coordinate system reveals **quantization cliffs**: specific bit-precisions
where certain weights become numerically unstable under Newton iteration
("ghost attractors"). The Mersenne Ghost Theorem shows that Mersenne numbers
`2^n - 1` are maximally fragile, with cliff at `k* = n + 2`.

## Mathematical Foundation

Everything follows from a single identity — Lifting the Exponent Lemma (LTE):

```
v₂(5^e₁ − 5^e₂) = v₂(e₁ − e₂) + 2
```

The exponential map `e ↦ 5^e` is a **scaled 2-adic isometry** with scale factor 4.
All core theorems are corollaries:

| # | Theorem | Proof |
|---|---------|-------|
| T1 | Quadratic convergence of Newton dlog map | LTE + linearisation |
| T2 | Trajectory separation: `n*(s) = ⌈log₂(s)⌉ − 1` | LTE + additive dynamics |
| T3 | Basin dichotomy: α=0 globally stable, no ghosts | Coset argument |
| T4 | Ghost formula: `e* = dlog(a+2, k)` for α=1 targets | T3 + LTE |
| T5 | Mersenne cliff: `k* = n+2`, `v₂(e_true) = n−2` | LTE at `e = 2^(n-2)` |
| T6 | Trivial Julia set — linearisable maps have no fractal structure | LTE + Berkovich |

### The Unified Cliff Constant

For generator `g ≡ 5 (mod 8)`:

```
c(g) = v₂(g − g₀) − 2     where g₀ = exp₂(−4) ≡ −123 (mod 2¹³)
```

Three regimes (`s = v₂(g − 5)`):
- `s < 7`:  `c(g) = s − 2`        (linear)
- `s > 7`:  `c(g) = 5`            (capped)
- `s = 7`:  `c(g) = 5 + v₂(m+1)` (boundary — the ε(n) correction)

## Features

- **Core arithmetic**: Fast modular inverse, Viglietta discrete log with
  10-bit LUT bootstrap (256-entry sweet spot), 2-adic logarithm
- **LTE dual addition**: `dual_add` — exact addition in `(v, α, e)` coordinates
  via LTE, no group representation round-trip
- **General 2-adic exp/log**: `padic_exp` and `padic_log` for arbitrary
  2-adic arguments (not just log(5))
- **DualNumber**: 2-adic dual-view decomposition `(v, α, e)` for any integer
- **Mahler calculus**: Binomial basis with Dirac/Volterra operators,
  proven boundary asymmetry `D∘T ≠ T∘D` at `e₁`
- **Operator calculus**: Shift/difference/average operators, Connes-style
  spectral triple, Newton projector
- **Newton basin analysis**: Full ghost detection, precision sweep
- **Thermodynamics**: Graded `v₂(e_true)` stability diagnostic, random
  comparison with z-scores
- **Gauge invariants**: Cyclic product invariants, tidal scalar `H = v + α·e`
- **CRT extension**: Combined `Z/(2^k · p)Z` analysis
- **Non-Abelian extension**: GL(2) holonomy with phase alignment
- **Kronecker cliff scoring**: Per-factor analysis for butterfly matrices
- **Trajectory Separation Theorem**: `n*(s) = ceil(log₂(s)) - 1` (proven,
  zero variance)
- **Fourier analysis**: DFT of step-count function, dyadic spectrum
- **p-adic root finding**: Newton, Halley, composed methods (order 2/3/4/8)
- **p-adic Newton dynamics**: Dynatomic polynomials, multiplier analysis,
  clean primes {7, 103, 181}, 4 proven theorems (special values, product
  formula, period-4 identity), 6-property landscape synthesis
  (deterministic, discrete, information-doubling, phase-transitional,
  universal, ultrametric)
- **Iwasawa decomposition**: GL(2) congruence filtration, LDU, commutator
  depth theorem
- **Iwasawa algebra**: Z₂[[G]] power series, augmentation ideal, Dirac
  operator classification
- **Mersenne Ghost Theorem**: Full proof, bootstrap optimality, cliff
  density theory (`Pr[c=0]=7/8`, `E[c]=1/4`)
- **Butterfly seed**: Dual-view Newton projector as butterfly-compilable seed,
  clean-prime analysis, OpenQASM circuit generation
- **PyTorch integration**: QuantizedMLP with ghost-regularized training
  (optional dependency)

## Installation

```bash
# Core (NumPy only)
pip install dual-view

# With PyTorch support
pip install dual-view[torch]

# With scipy support (bridge module — depth histogram seed S1)
pip install dual-view[scipy]

# Development (testing, linting, type-checking)
pip install dual-view[dev]
```

## Quick Start

```python
from dual_view import DualNumber, TwoAdicProcessor

# Decompose an integer into dual coordinates
d = DualNumber(42, k=16)
print(d)          # DualNumber(42, k=16) = 2^1 · +5^2
print(d.coords()) # (1, 0, 2)

# Arithmetic in coordinate space
proc = TwoAdicProcessor(16)
a = DualNumber(3, 16)
b = DualNumber(7, 16)
c = proc.mul(a, b)  # 21
```

```python
from dual_view import SeedThermodynamics
import numpy as np

# Analyse ghost cliffs in a weight matrix
W = np.random.randint(0, 256, size=(32, 32)).astype(np.int64)
st = SeedThermodynamics(k=8)
stats = st.analyse(W)
print(stats["mean_v2_e"], stats["cliff_risk"])
```

## Module Overview

| Module | Description |
|--------|-------------|
| `core` | DualNumber, modular inverse, 2-adic exp/log, cliff centre g₀ |
| `exponent` | Additive coordinate chart on Z/2^(k-2) |
| `mahler` | Mahler basis, Dirac/Volterra operators, boundary asymmetry |
| `operators` | Symbolic operator algebra (shift, difference, average) — O(N) per application, decorative for k ≥ 9 |
| `basin` | Newton basin analysis, ghost detection |
| `thermodynamics` | Graded `v₂(e_true)` weight stability diagnostic |
| `regularization` | Ghost-aware regularisation for NN training |
| `gauge` | Gauge invariants for weighted cyclic operators |
| `crt` | CRT extension to composite moduli |
| `nonabelian` | GL(2) gauge theory, holonomy invariants |
| `scaling` | Float-to-int quantization scaling |
| `visualise` | Cliff matrix rendering and ASCII heatmaps |
| `butterfly` | Kronecker factor cliff scoring |
| `separation` | Trajectory Separation Theorem |
| `fourier` | DFT of Newton step-count function |
| `padic_roots` | Multi-order p-adic root finding (Newton, Halley) |
| `newton_dynamics` | Dynatomic polynomials, multipliers, clean primes |
| `iwasawa` | GL(2) congruence filtration, LDU decomposition |
| `iwasawa_algebra` | Iwasawa algebra Z₂[[G]], profinite filtered modules |
| `mersenne` | Mersenne Ghost Theorem, cliff constant proofs, bootstrap optimality |
| `isometry` | Exponential isometry, operator algebra theorems |
| `butterfly_seed` | Dual-view Newton projector, clean-prime analysis, QASM |
| `bridge` | Three-seed 2-adic weight analysis (depth histogram, map, sign) |
| `training` | PyTorch QuantizedMLP with ghost regularisation |

```python
from dual_view.newton_dynamics import dynatomic_polynomial, compute_iterates

# Compute dynatomic polynomial for period-2
iters = compute_iterates(2)
phi2 = dynatomic_polynomial(2, iters)  # → [2, 5, 20] = 20u² + 5u + 2

# Period-4: degree 24, Φ(0) = 2¹², Φ(1) = 3³⁶
iters = compute_iterates(4)
phi4 = dynatomic_polynomial(4, iters)  # 25 coefficients
```

## Running Tests

```bash
pytest tests/ -v
```

## Demo

```bash
python -m dual_view.demo
python -m dual_view.demo --quick  # skip slow tests
```

## License

MIT
