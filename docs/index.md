# Dual-View Documentation Index

## Project Overview

**dual-view** is a unified mathematical framework for 2-adic number systems, p-adic Newton dynamics, noncommutative spectral geometry, and quantum butterfly compilation. The core insight is that every odd integer modulo `2^k` decomposes uniquely as:

```
n = 2^v · (-1)^α · 5^e   (mod 2^k)
```

where:
- `v = v₂(n)` — the 2-adic valuation (power of 2 dividing n)
- `α ∈ {0, 1}` — sign sector (0 for `n ≡ 1 mod 4`, 1 for `n ≡ 3 mod 4`)
- `e ∈ [0, 2^(k-2))` — discrete logarithm base 5

This dual-view coordinate system reveals **quantization cliffs**: specific bit-precisions where weights become numerically unstable under Newton iteration ("ghost attractors").

## Module Reference

| Module | Lines | Description |
|--------|-------|-------------|
| `core.py` | 523 | DualNumber, modular inverse, 2-adic exp/log, cliff centre g₀ |
| `exponent.py` | 81 | ExponentSpace — additive coordinate chart on Z/2^(k-2) |
| `mahler.py` | 123 | Mahler basis, Dirac/Volterra operators, boundary asymmetry |
| `operators.py` | 215 | OperatorContext, SpectralTriple, NewtonProjector |
| `basin.py` | 346 | BasinExplorer, LayerGhostDiagnosticV2, GhostHunt |
| `thermodynamics.py` | 285 | SeedThermodynamics — graded 2-adic stability diagnostics |
| `regularization.py` | 150 | GhostMap, ghost_penalty, local_ratio_gradient |
| `gauge.py` | 120 | Gauge invariants for weighted cyclic operators |
| `crt.py` | 205 | CRT extension to Z/(2^k·p)Z |
| `nonabelian.py` | 214 | GL(2) gauge theory with phase alignment |
| `scaling.py` | 106 | Float-to-int quantization scaling |
| `visualise.py` | 135 | Cliff matrix rendering, ASCII heatmaps |
| `butterfly.py` | 145 | Kronecker factor cliff scoring |
| `separation.py` | 185 | Trajectory Separation Theorem |
| `fourier.py` | 185 | DFT of Newton step-count function |
| `padic_roots.py` | 311 | Multi-order p-adic root finding |
| `iwasawa.py` | 235 | GL(2) congruence filtration, LDU |
| `mersenne.py` | 585 | Mersenne Ghost Theorem, bootstrap optimality, cliff constant proofs |
| `isometry.py` | 257 | Exponential isometry, operator algebra theorems |
| `butterfly_seed.py` | 502 | Dual-view Newton projector, clean-prime analysis, QASM |
| `bridge.py` | 326 | Three-seed 2-adic weight analysis (depth histogram, map, sign) |
| `training.py` | 230 | PyTorch QuantizedMLP with ghost regularization |
| `newton_dynamics/` | 5 modules | p-adic Newton dynamics for N(x) = (2x³+1)/(3x²) |
| `demo.py` | 223 | Runnable demonstration suite |

## Theorem Reference

| ID | Theorem | Module | Status |
|----|---------|--------|--------|
| T1 | Gain law: `v(e_new - e_true) = 2j+1` (quadratic + LTE bonus) | `core.py` | Proven |
| T1b | General 2-adic exp/log: `exp(x)`, `log(g)` via exact arithmetic | `core.py` | Verified |
| T1c | Cliff density: `Pr[c=0]=7/8`, `E[c]=1/4` | `core.py` | Proven |
| T1d | Real vs 2-adic reconciliation: 10–12× speedup over squaring loop | — | Verified |
| — | Divisor optimality: `d=2` unique optimal Newton divisor | `core.py` | Proven |
| T2 | Trajectory separation: `n*(s) = ceil(log₂(s)) - 1` | `separation.py` | Proven, zero variance |
| T3 | Basin dichotomy: α=0 globally stable | `basin.py` | Proven |
| T4 | Ghost formula: `e* = dlog(a+2, k)` for α=1 | `basin.py` | Proven |
| T5 | Mersenne cliff: `k* = n+2` for `w = 2^n - 1` | `mersenne.py` | Verified n=3..11 |
| T6a | Exponential map isometry: `v₂(5^e-1) = v₂(e)+2` | `isometry.py` | Proven |
| T6b | Operator algebra: `avg² = N·avg`, `D·avg = avg·D = 0` | `isometry.py` | Proven |
| T6d | Mahler basis: `D∘T = id` on `ker(ε)`, `T∘D = id` on `n≥2` | `mahler.py` | Proven |
| T6c | Trace-mod-p independence (GL(2) holonomy) | `isometry.py` | Statistical |
| — | Commutator depth: `depth([M,N]) ≥ depth(M)+depth(N)` | `iwasawa.py` | Verified |
| — | Mersenne cliff constant: `c(g) = v₂(g - exp₂(-4)) - 2` | `mersenne.py` | Proven |
| — | p-adic convergence law: `v_p(x_n-x*) = m^n·v_p(x₀-x*)` | `padic_roots.py` | Proven, zero variance |
| — | Newton correction uniformity (first step) | `padic_roots.py` | Empirical |
| T7 | Dynatomic special value: `Φ_n^*(0) = 2^{μ₃(n)/6}` | `newton_dynamics/` | Proven |
| T8 | Dynatomic special value: `Φ_n^*(1) = 3^{μ₃(n)/2}` | `newton_dynamics/` | Proven |
| T9 | Universal multiplier product: `∏μ = 6^{μ₃(n)/6}` | `newton_dynamics/` | Proven |

## Fast Function Reference

| Function | Location | Purpose |
|----------|----------|---------|
| `DualNumber(n, k)` | `core.py` | Decompose n into (v, α, e) |
| `DualNumber.from_coords(v, α, e, k)` | `core.py` | Build from coordinates |
| `modinv_newton(a, k)` | `core.py` | `a⁻¹ mod 2^k` via Newton |
| `two_adic_dlog(a, k)` | `core.py` | `(α, e)` decomposition |
| `two_adic_log5(k)` | `core.py` | 2-adic log of 5, cached |
| `ExponentSpace(g, k)` | `exponent.py` | `e ↦ g^e` with difference calculus |
| `OperatorContext(k, g)` | `operators.py` | I, S, diff, avg, M(h) |
| `BasinExplorer(k, g, a)` | `basin.py` | Newton basin portrait |
| `SeedThermodynamics(k)` | `thermodynamics.py` | Weight stability analysis |
| `GhostMap(k)` | `regularization.py` | Precomputed convergence ratios |
| `ghost_penalty(W, gm)` | `regularization.py` | Penalty + surrogate gradient |
| `scale_weights(W, scale)` | `scaling.py` | Float-to-int scaling |
| `cliff_matrix(st, shape)` | `visualise.py` | Reshape cliff scores |
| `KroneckerCliffScorer(factors)` | `butterfly.py` | Factor cliff scoring |
| `mersenne_cliff_table(n_max)` | `mersenne.py` | Mersenne cliff thresholds |
| `cliff_constant(g, k)` | `mersenne.py` | Mersenne cliff constant `c = v₂(log₂(g)/4+1)` |
| `exp2_neg4(k)` | `mersenne.py` | 2-adic zero of `log₂(g)/4+1` |
| `mersenne_cliff_theorem()` | `mersenne.py` | Verify the full Mersenne cliff theorem |
| `lift_root(a, p, k)` | `padic_roots.py` | Hensel lift cube root to mod p^k |
| `verify_isometry(k)` | `isometry.py` | Verify T6a empirically |
| `compute_iterates(k)` | `newton_dynamics/` | Compute Newton iterates A_d, B_d |
| `dynatomic_polynomial(n, iters)` | `newton_dynamics/` | Dynatomic polynomial Φ_n^*(u) |
| `poly_mul(p, q)` | `newton_dynamics/` | Polynomial multiplication |
| `tonelli_shanks(n, p)` | `newton_dynamics/` | Modular square root |
| `is_cube(a, p)` | `newton_dynamics/` | Cube residue modulo p |
| `analyze_prime(p)` | `butterfly_seed.py` | Classify prime by Newton functional graph |
| `DualViewSeed(k, a)` | `butterfly_seed.py` | Position-dependent butterfly seed builder |
| `dual_view_qasm_emitter(k, a, p)` | `butterfly_seed.py` | OpenQASM 2.0 circuit generator |
| `SpectralThermodynamics.analyze(S)` | `bridge.py` | Spectral classification of seed matrices |
| `ButterflyBridge(k)` | `bridge.py` | Three-seed unified weight analysis |
| `quantize(W, k)` | `bridge.py` | Symmetric min-max quantisation to Z/2^k |

## Open Problems

See `research_opportunities.md` for the full list. Key open questions:

1. **Secondary correction**: `ε(n) = v₂(n) - 1` in Mersenne cliff formula at powers of 2
2. **Bootstrap optimization**: `eprec₀ = k/2` not deployed in core dlog
3. **CRT stability correlation**: Empirical correlation without theoretical prediction
4. **Popcount compression**: Unverified correlation from parallel report
5. **Trace-mod-p bridge**: Why does α(det H) correlate with Tr(H) mod p?

## Documents

| Document | Description |
|----------|-------------|
| `index.md` | This file — overview, file table, theorem reference |
| `mathematics.md` | Full mathematical background |
| `api.md` | Complete API reference |
| `mersenne_ghost_theorem.md` | Mersenne Ghost Theorem proof |
| `bug_history.md` | Bug fix narrative and audit history |
| `research_opportunities.md` | Open problems, experiments, future directions |
| `newton_dynamics/index.md` | p-adic Newton dynamics — overview, quick-start |
| `newton_dynamics/theorems.md` | Full proofs of all 4 theorems |
| `newton_dynamics/computation_data.md` | Period-by-period data tables |
| `newton_dynamics/clean_primes.md` | Clean prime analysis (7, 103, 181) |
| `newton_dynamics/status.md` | Research status and open problems |
| `unified_dual_view_butterfly_guidance.md` | Synthesis: dual-view × butterfly compiler research roadmap |
| `audit.md` | Package audit: bug register, scaling laws, boundary test theorems |
| `padic_newton_dynamics.md` | p-Adic Newton dynamics: higher-order methods, convergence law verification |
