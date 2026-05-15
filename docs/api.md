# API Reference

## Core Module (`dual_view.core`)

### Functions

#### `_mask(k: int) -> int`
Bit mask `(1 << k) - 1`.

#### `_valuation(n: int) -> int`
2-adic valuation `v₂(n)`. Returns `float('inf')` for `n = 0`.

#### `modinv_newton(a: int, k: int) -> int`
`a⁻¹ mod 2^k` via quadratic Newton lifting. Requires `a` odd.

#### `two_adic_log5(k: int) -> int`
2-adic logarithm of 5 truncated to `k` bits. The unique `L ∈ Z₂` satisfying `exp₂(L) = 5`. Cached per `k` (LRU, 256 entries).

Used as derivative scaling factor in Newton dlog: the derivative of `e ↦ 5^e` is `5^e · (L >> 2)`.

#### `two_adic_dlog(a: int, k: int, L: Optional[int] = None) -> Optional[Tuple[int, int]]`
Full 2-adic decomposition of odd part: `a = (-1)^α · 5^e (mod 2^k)`.
Returns `(alpha, e)` or `None` if `a` is even. Uses 8-bit LUT bootstrap for `k ≤ 34`, bit-by-bit to `k/2` otherwise.

#### `run_all_tests(k: int = 16, verbose: bool = True) -> None`
Self-check for core arithmetic: round-trip, multiplication, inversion, powering.

### Classes

#### `DualNumber(n: int, k: int = 64)`
2-adic dual-view decomposition of `n` modulo `2^k`.

**Attributes**: `v` (valuation, float('inf') for zero), `alpha` (0 or 1), `e` (exponent), `is_zero`, `value` (the integer `n mod 2^k`).

**Methods**:
- `verify()` — round-trip check: coordinates → integer matches stored integer
- `coords()` — returns `(v, alpha, e)` tuple
- `from_coords(v, alpha, e, k)` — classmethod, build from coordinates

#### `TwoAdicProcessor(k: int = 64)`
Arithmetic on DualNumbers in coordinate space.

**Methods**:
- `mul(a, b)` — multiply, coordinates add componentwise
- `inv(a)` — invert unit (v=0)
- `pow(a, n)` — integer power (negative → inverse)
- `dlog(a)` — returns `(v, alpha, e)`

## Exponent Module (`dual_view.exponent`)

#### `ExponentSpace(g: int, k: int)`
Model the multiplicative group `⟨g⟩ ≅ Z/N` where `N = 2^(k-2)`.

**Methods**:
- `lift(e)` — `g^e mod 2^k`
- `difference(f, e)` — forward difference `(Df)(e) = f(e+1) - f(e)`
- `integrate(f)` — discrete integral `Σ f(e) mod 2^k`
- `is_eigenfunction(f, e)` — check `D(g^e) = (g-1)·g^e`

## Operators Module (`dual_view.operators`)

#### `OperatorContext(k: int, g: int)`
Operator environment. Provides:
- `I` — identity
- `S` — shift: `(Sf)(e) = f(e+1 mod N)`
- `diff` — forward difference `I - S`
- `avg` — summation: `(avg f)(e) = Σ f(t)`
- `multiply_by(h, name)` — multiplication operator `M_h`
- `g` — the generator

#### `SpectralTriple(ctx: OperatorContext)`
Connes-style spectral triple `(A, H, D)`.

**Methods**:
- `one_form(h)` — `ω = [D, M_h]`
- `curvature(h1, h2)` — `Ω = [ω₁, ω₂]`
- `gauge_transform(potential)` — `D' = D + ω`
- `trace(f)` — `Tr(f) = Σ f(e)`

#### `NewtonProjector(ctx: OperatorContext, a: int, steps: int = 16)`
Newton projection operator for solving `5^e ≡ a`.

**Methods**:
- `project_point(e_guess)` — apply Newton refinement

## Basin Module (`dual_view.basin`)

#### `BasinExplorer(k: int, g: int, target_a: int)`
Newton basin landscape analysis.

**Methods**:
- `newton_step(e)` — single Newton iteration
- `classify(e0)` — returns `(fate, value, path)`, fate ∈ {'converged', 'cycle', 'diverged'}
- `portrait()` — full basin portrait dict over all seeds
- `fate_vector()` — compact encoding (0=converged, 1=cycle, 2=diverged)

#### `LayerGhostDiagnosticV2(k: int, g: int = 5, max_iter: int = 100)`
Per-layer diagnostic for weight matrices.

**Methods**:
- `diagnostic_matrix(W)` — returns `(fate, conv_ratio, ghost_ratio, mean_e, v2_e)`

#### `GhostHunt(g: int = 5, max_iter: int = 100)`
Precision threshold hunting.

**Methods**:
- `precision_threshold_sweep(k_min, k_max, target_e)` — print sweep table
- `quantization_cliff(W, k_min, k_max)` — ghost density vs k for a matrix

#### `precision_sweep(k_min, k_max, g=5, target_e=2) -> List[Tuple[int, float]]`
Standalone precision sweep.

## Thermodynamics Module (`dual_view.thermodynamics`)

#### `SeedThermodynamics(k: int = 16, g: int = 5)`
Graded 2-adic weight stability diagnostics.

**Coordinate analysis**:
- `weight_coordinates(w)` — return `(v, α, e)`
- `analyse(W)` — return dict with alpha_fraction, mean_v2_e, std_v2_e, cliff_risk
- `mersenne_cliff_score(w)` — returns `k* = v₂(e_true) + 2`
- `compare_to_random(W, n_samples)` — z-score comparison

**Precision-sweep analysis**:
- `__call__(weights, k_range)` — configure for sweep
- `compute()` — lazy computation of ghost profiles
- `profiles` — property: per-weight ghost density profiles
- `cliffs` — property: per-weight cliff thresholds
- `stable_weights(max_k=None)` — indices of stable weights
- `ghost_weights(max_k=None)` — indices of ghost-affected weights
- `cliff_histogram()` — cliff distribution
- `summary()` — aggregate statistics
- `report()` — formatted ASCII report with bar chart

## Regularization Module (`dual_view.regularization`)

#### `GhostMap(k: int, g: int = 5, max_iter: int = 64)`
Precomputed convergence ratios. `k` limited to ≤ 10.

**Methods**:
- `ratio(a)` — convergence ratio for weight `a`
- `nearest_stable(a, search_radius=4)` — nearest weight with better ratio

#### `local_ratio_gradient(weight_int, ghost_map) -> List[Tuple[int, float]]`
Local improvement candidates.

#### `ghost_penalty(weights, ghost_map, step_scale=1.0) -> Tuple[float, ndarray]`
Ghost regularization penalty and surrogate gradient.

## CRT Module (`dual_view.crt`)

#### `CRTDualNumber(n: int, k: int, p: int, g_p: int)`
Element of `Z/(2^k·p)Z` with dual coordinates.

#### `CRTDualProcessor(k: int, p: int, g_p: Optional[int] = None)`
CRT arithmetic.

**Methods**:
- `crt_reconstruct(r2, rp)` — CRT combine residues
- `mul(A, B)` — multiply two CRTDualNumbers
- `product(weights)` — product of raw integers
- `cycle_product(numbers)` — product of CRTDualNumber list
- `convergence_ratio_2adic(P)` — 2-adic convergence ratio

#### `combined_stability(k, p, num_cycles=50, cycle_length=4) -> Dict`
Randomised correlation test.

## Nonabelian Module (`dual_view.nonabelian`)

#### `NonAbelianCRTDual(k: int, p: int)`
GL(2) gauge theory on a cycle.

**Methods**:
- `holonomy(mats)` — product of all matrices
- `invariants(mats)` — dict with det_mod2k, alpha_det, trace_modp, crt views
- `convergence_ratio_full(mats)` — ghost ratio of determinant

#### `ramp_break_strength(...) -> Dict`
**Deprecated**: Use `phase_alignment_experiment()` instead.

#### `phase_alignment_experiment(k, p, N_cycle=4, n_cycles=30) -> Dict`
Test α-sector flip under perturbation.

## Scaling Module (`dual_view.scaling`)

#### `scale_weights(W, scale, mode='round', ensure_odd=False) -> Tuple[ndarray, Dict]`
Float-to-int quantization.

**Meta keys**: scale, mode, ensure_odd, n_even, n_zero, range, v2_hist

#### `auto_scale(W, target_bits=8) -> float`
Auto-compute scale from 99th percentile.

#### `common_scales() -> Dict[str, float]`
Standard bit depths: INT7=64, INT8=128, INT9=256, etc.

## Visualise Module (`dual_view.visualise`)

#### `cliff_matrix(st, original_shape) -> ndarray`
Reshape cliff scores to original shape.

#### `sector_matrix(weights_int, k) -> ndarray`
Map odd weights to α-sector (0/1), NaN for evens.

#### `valuation_matrix(weights_int) -> ndarray`
Map weights to v₂, -1 for zeros.

#### `print_cliff_ascii(C, title, max_rows=40, max_cols=80)`
ASCII density heatmap.

#### `cliff_stats_by_layer(layers) -> str`
Per-layer statistics table.

## Butterfly Module (`dual_view.butterfly`)

#### `KroneckerCliffScorer(factors, k_range, scale_mode='round', ensure_odd=True)`
Kronecker factor cliff scoring.

**Methods**:
- `score_factors()` — run analysis on all factors
- `composition_cliff()` — min of factor cliffs
- `fragile_factors(threshold=6.0)` — indices below threshold
- `print_report()` — formatted results

#### `semiring_cliff_score(factor_cliffs, semiring='standard') -> Optional[float]`
Aggregate under different semirings (standard/tropical/boolean).

## Advanced Theorem Modules

### Separation (`dual_view.separation`)
- `newton_trajectory(a, k, e_seed, steps=10)` — per-step Newton history
- `separation_step(a, a_prime, k, e_seed)` — first divergence step
- `predicted_separation(s, method_order=2)` — theoretical `n*(s)`
- `verify_separation(k, s_values, n_trials=50)` — zero-variance verification
- `ultrametric_ball_tree(k, e_true, depth=3)` — ASCII tree
- `step_count_profile(k, e_true)` — v₂ level counts

### Fourier (`dual_view.fourier`)
- `step_count_fn(k, e_true)` — numeric step-count function
- `analytic_step_count(k, e_true)` — closed-form O(N) construction
- `dft(f)` — numpy FFT wrapper
- `power_spectrum(f)` — |DFT|²
- `dyadic_coefficients(f)` — extract at N/2, N/4, ...
- `analytic_coefficients(k)` — closed-form Fourier coefficients
- `fourier_summary(k, e_true)` — complete analysis
- `ultrametric_uncertainty(k)` — uncertainty principle statement

### p-adic Roots (`dual_view.padic_roots`)
- `newton_step(x, a, pk)` — order 2
- `halley_step(x, a, pk)` — order 3
- `newton2_step(x, a, pk)` — order 4 (composed Newton)
- `newton3_step(x, a, pk)` — order 8 (triple Newton)
- `convergence_profile(x0, a, p, k, step_fn, x_true)` — track v_p
- `compare_methods(p, k, n_trials)` — rate comparison
- `verify_order(primes, k, n_trials)` — verify convergence ratio
- `newton_correction_uniformity(p, k, n_seeds)` — chi-square test
- `popcount_compression(k, n_trials)` — popcount correlation

### Iwasawa (`dual_view.iwasawa`)
- `congruence_depth(M, k)` — depth in GL(2) filtration
- `filtration_residue(M, depth, k)` — gl(2, F₂) direction
- `ldu_decompose(M, k)` — LDU factorisation
- `matrix_coordinates(M, k)` — full coordinate decomposition
- `holonomy_depth_profile(k, p, ...)` — depth under perturbation
- `filtration_portrait(k)` — GL(2) quotient sizes
- `matrix_commutator(M, N, k)` — `[M,N]` computation
- `verify_commutator_depth(k, ...)` — depth theorem verification
- `MatrixCoordinates` — dataclass

### Mersenne (`dual_view.mersenne`)
- `mersenne_coordinates(n, k)` — `(α, e_true, v₂(e_true))`
- `verify_core_identity(n_max)` — 5^(2^(n-2)) ≡ 1 - 2^n
- `mersenne_cliff_table(n_max)` — cliff thresholds
- `bootstrap_cost(eprec0, k)` — Viglietta bit-cost
- `optimal_bootstrap(k_values)` — minimiser search
- `compare_bootstrap_strategies(k_values)` — sqrt vs k/2 vs LUT
- `dlog_with_lut(a, k, b=8)` — LUT-based dlog
- `verify_lut_dlog(k, b=8, n_trials)` — correctness check

### Isometry (`dual_view.isometry`)
- `verify_isometry(k, n_trials)` — v₂(5^e-1) = v₂(e)+2
- `isometry_pair_test(k, n_trials)` — pair form
- `isometry_summary(k)` — full conditioning picture
- `verify_operator_algebra(k_values)` — avg², D·avg, avg·D
- `trace_alpha_independence(k, p, ...)` — chi-square test
- `trace_exponent_independence(k, p, ...)` — ANOVA F-test
- `exponent_valuation_profile(k, n_samples)` — v₂(e_true) distribution

## Training Module (`dual_view.training`)

Requires `torch` (optional dependency).

#### `QuantizedMLP(k=8)`
Two-layer MLP (784→128→10) for MNIST.

**Methods**:
- `forward(x)` — forward pass
- `get_weights_numpy()` — quantized weights as numpy

#### `build_loaders(batch_size=64, data_root='./data')`
MNIST data loaders.

#### `train(model, train_loader, test_loader, epochs=5, lr=0.001, ghost_map=None, ghost_scale=0.01, device=None, use_thermodynamics=False, thermo_k=8) -> Dict`
Training loop with optional ghost regularisation and thermodynamics tracking. Returns history with loss, acc, grad_norm, ghost_penalty, (optionally cliff_risk, alpha_fraction).
