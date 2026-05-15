# dual-view

**2-adic dual-view diagnostics for quantized neural network weights.**

This package provides a mathematical framework for analyzing quantized integer
weight matrices using 2-adic arithmetic. The core insight is that every odd
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

## Features

- **Core arithmetic**: Fast modular inverse, Viglietta discrete log with
  8-bit LUT bootstrap, 2-adic logarithm
- **DualNumber**: 2-adic dual-view decomposition `(v, α, e)` for any integer
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
- **Iwasawa decomposition**: GL(2) congruence filtration, LDU, commutator
  depth theorem
- **Mersenne Ghost Theorem**: Full proof and bootstrap optimality
- **PyTorch integration**: QuantizedMLP with ghost-regularized training
  (optional dependency)

## Installation

```bash
# Core (NumPy only)
pip install dual-view

# With PyTorch support
pip install dual-view[torch]

# Development
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
| `core` | DualNumber, modular inverse, Viglietta discrete log |
| `exponent` | Additive coordinate chart on Z/2^(k-2) |
| `operators` | Symbolic operator algebra (shift, difference, average) |
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
| `iwasawa` | GL(2) congruence filtration, LDU decomposition |
| `mersenne` | Mersenne Ghost Theorem, bootstrap optimality |
| `isometry` | Exponential isometry, operator algebra theorems |
| `training` | PyTorch QuantizedMLP with ghost regularisation |

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
