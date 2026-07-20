# dual-view

**A unified mathematical framework spanning 2-adic number theory, p-adic Newton dynamics, noncommutative spectral geometry, and quantum butterfly compilation.**

Every odd integer modulo `2^k` decomposes uniquely as a **dual-view coordinate triple**:

```
n = 2^v · (-1)^α · 5^e   (mod 2^k)
```

where `v` is the 2-adic valuation, `α` is the sign sector, and `e` is the discrete logarithm base 5. This decomposition reveals **quantization cliffs** — bit-precisions where weights become numerically unstable under Newton iteration. The Mersenne Ghost Theorem proves that Mersenne numbers `2^n - 1` are maximally fragile, with cliff at `k* = n + 2`.

## Highlights

### 15 Clean Primes (Verified to 30,000,000)

```
{7, 31, 41, 59, 103, 181, 359, 659, 811, 8111, 14159, 31741, 115679, 162251, 403549}
```

A prime `p` is *clean* for `N(x) = (2x^3+1)/(3x^2)` over `F_p` iff the Newton functional graph is a rooted forest with no cycles (ghost attractors). The original conjecture of {7, 103, 181} was refuted by exhaustive DFS search to 30M — 15 primes found, set still appears finite.

### Nilpotent Basin Structure (S^M = 0)

For every clean prime, ordering elements by basin depth makes the Newton shift operator `S` strictly upper-triangular, hence nilpotent: `S^M = 0`. The resolvent expands as an exact finite Neumann series:

```
(I - S)^{-1} = I + S + S^2 + ... + S^{M-1}
```

### Butterfly Routing Compiler

The Newton iteration `x_{n+1} = N(x_n)` parallelises to depth `⌈log₂(M)⌉` via a butterfly routing network. Each stage applies N `2^s` times, converging in `⌈log₂(M)⌉` stages instead of M serial steps. At p=403549 (M=1223): 11 butterfly stages vs 1223 serial Newton steps — a **111× speedup**.

### Basin-Optimised QASM Circuit

For clean primes, the quantum circuit uses `D = ⌈log₂(M)⌉` exponent qubits instead of the standard `k-2`, dropping circuit depth from `O(k²)` to `O(D²)`. At k=32, p=403549: D² = 121 vs k² = 1024 controlled-phase gates — **8.5× reduction**.

### Discriminant Theorem (Δ = 108(x³-1))

The identical early-depth structure shared by all 3-root clean primes (depths 0-5) is proven: the Newton preimage cubic discriminant `Δ = 108(x³-1)` vanishes at the roots, and the small-depth x-values are constrained to a fixed set of algebraic expressions, independent of p.

## Mathematical Foundation

Everything follows from Lifting the Exponent Lemma (LTE):

```
v₂(5^e₁ − 5^e₂) = v₂(e₁ − e₂) + 2
```

The exponential map `e ↦ 5^e` is a **scaled 2-adic isometry** with scale factor 4.

| # | Theorem | Proof |
|---|---------|-------|
| T1 | Quadratic convergence of Newton dlog map | LTE + linearisation |
| T2 | Trajectory separation: `n*(s) = ⌈log₂(s)⌉ − 1` | LTE + additive dynamics |
| T3 | Basin dichotomy: α=0 globally stable, no ghosts | Coset argument |
| T4 | Ghost formula: `e* = dlog(a+2, k)` for α=1 targets | T3 + LTE |
| T5 | Mersenne cliff: `k* = n+2`, `v₂(e_true) = n−2` | LTE at `e = 2^(n-2)` |
| T6 | Trivial Julia set — linearisable maps have no fractal structure | LTE + Berkovich |

## Module Overview

| Module | Description |
|--------|-------------|
| `core` | DualNumber, modular inverse, 2-adic exp/log, cliff centre g₀ |
| `exponent` | Additive coordinate chart on Z/2^(k-2) |
| `mahler` | Mahler basis, Dirac/Volterra operators, boundary asymmetry |
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
| `newton_dynamics` | Dynatomic polynomials, multipliers, clean primes |
| `iwasawa` | GL(2) congruence filtration, LDU decomposition |
| `iwasawa_algebra` | Iwasawa algebra Z₂[[G]], profinite filtered modules |
| `mersenne` | Mersenne Ghost Theorem, cliff constant proofs |
| `isometry` | Exponential isometry, operator algebra theorems |
| `butterfly_seed` | Dual-view Newton projector, clean-prime analysis |
| `butterfly_emitter` | Basin-optimised OpenQASM emission (nilpotent depth reduction) |
| `bridge` | Three-seed 2-adic weight analysis (depth, map, sign) |
| `training` | PyTorch QuantizedMLP with ghost regularisation |

## Research

| File | What it proves |
|------|----------------|
| `research/routing_simulator.py` | Classical butterfly routing convergence on all 15 primes |
| `research/butterfly_compiler.py` | Nilpotent shift operator S, Neumann series, routing stages |
| `research/expA-D` | Basin Newton operator, depth spectrum, forest isomorphism, Kronecker clean signals |
| `research/REPORT.md` | Full 9-experiment results |
| `research/BUTTERFLY_COMPILER.md` | Compiler prototype documentation |

## Installation

```bash
pip install dual-view
pip install dual-view[torch]   # with PyTorch support
pip install dual-view[dev]     # testing, linting, type-checking
```

## Quick Start

```python
from dual_view import DualNumber, TwoAdicProcessor

d = DualNumber(42, k=16)
print(d)          # DualNumber(42, k=16) = 2^1 · +5^2
print(d.coords()) # (1, 0, 2)

proc = TwoAdicProcessor(16)
c = proc.mul(DualNumber(3, 16), DualNumber(7, 16))  # 21
```

```python
from dual_view import analyze_prime
from dual_view.butterfly_emitter import basin_qasm_emitter

prof = analyze_prime(103)
print(f"Clean: {prof.is_clean}, M={prof.nilpotency_index}")  # Clean: True, M=15

qasm = basin_qasm_emitter(103, k=16, target=5)
print(qasm[:500])  # Basin-optimised circuit (4 exponent qubits vs 14)
```

## Running Tests

```bash
pytest tests/ -v
```

## License

MIT
