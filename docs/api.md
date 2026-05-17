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

#### `padic_exp(x: int, k: int) -> int`
General 2-adic exponential `exp(x) mod 2^k`. Requires `v₂(x) ≥ 2`. Uses exact integer arithmetic with analytic valuation tracking for termination.

#### `padic_log(g: int, k: int) -> int`
General 2-adic logarithm `log(g) mod 2^k`. Requires `g ≡ 1 mod 4` (i.e. `v₂(g−1) ≥ 2`). Uses exact integer arithmetic.

#### `g0(k: int) -> int`
The cliff centre `g₀ = exp₂(−4) mod 2^k`. The unique 2-adic unit with `log(g₀) = −4`. The hardware approximation `−123` agrees with `g₀` to 13 bits.

#### `dlog_residual_tracking(a: int, k: int, L: Optional[int] = None) -> Tuple[int, List[Dict]]`
Viglietta discrete log with normalised residual tracking at each Newton step. Returns `(e, history)` where each entry in history contains `tau_before`/`tau_after` (the residual `(5^e · a⁻¹ − 1) / 4`) and its 2-adic valuation, confirming the quadratic convergence gain law. Requires `a ≡ 1 (mod 4)`.

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

## Mahler Module (`dual_view.mahler`)

#### `MahlerCalculus`
Operations on the Mahler (binomial) basis for integer-valued functions.

**Static methods**:
- `mahler_polynomial(n, x)` — binomial coefficient `C(x, n)`
- `to_mahler(f, max_degree)` — convert function values to Mahler coefficients via finite-difference table
- `from_mahler(coeffs, x)` — evaluate function from Mahler coefficients at point `x`
- `dirac_operator(coeffs)` — forward difference on coefficients: `D(a₀, a₁, …) = (−a₁, −a₂, …)`
- `volterra_operator(coeffs)` — right inverse of Dirac: `T(0, a₁, …) = (0, 0, −a₁, …)`; requires `a₀ = 0`
- `truncate(coeffs, k)` — reduce coefficients modulo `2^k`

**Boundary behaviour**:
```
D∘T = id  on  ker(ε) = span{e_n : n ≥ 1}
T∘D = id  on  span{e_n : n ≥ 2}   (NOT on all of ker(ε))
```

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
Randomised correlation test: correlates `v₂(prod)` with the change in `v₂` under 2-adic multiplicative perturbation `w → w·(1+2^t)`.

## Nonabelian Module (`dual_view.nonabelian`)

#### `NonAbelianCRTDual(k: int, p: int)`
GL(2) gauge theory on a cycle.

**Methods**:
- `holonomy(mats)` — product of all matrices
- `invariants(mats)` — dict with det_mod2k, alpha_det, trace_modp, crt views
- `convergence_ratio_full(mats)` — ghost ratio of determinant

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

#### `show_dual_bits(n, k=16, label="") -> str`
2-adic bit structure with dual-view annotation. Displays the bit pattern of `n` modulo `2^k` annotated with valuation (`v`), the fixed leading-1 bit (`1`), the α sign bit (`a`), and the discrete-log exponent bits (`e`). For odd `n` the full `(v, α, e)` decomposition and a reconstruction check are included.

## Butterfly Seed Module (`dual_view.butterfly_seed`)

Bridges the 2-adic Newton dynamics with the butterfly compiler's position-dependent operad framework.

#### `analyze_prime(p: int) -> CleanPrimeProfile`
Classify prime `p` by the thermodynamics of its Newton functional graph over `F_p^*`. Returns `CleanPrimeProfile` with fields: `is_clean`, `roots` (cube roots of 1), `nilpotency_index`, `basin_ordering`, `tree_depths`, `obstruction` ("clean", "ghost_cycle", "pole_chain", "mixed").

Clean primes (known: 7, 103, 181) have functional graphs that are rooted forests with exactly 3 trees.

#### `CleanPrimeProfile` (dataclass)
- `p: int` — the prime
- `is_clean: bool` — whether the functional graph is a 3-tree forest
- `roots: Tuple[int, ...]` — cube roots of 1 mod p
- `nilpotency_index: int` — max basin depth (0 if not clean)
- `basin_ordering: List[int]` — F_p^* ordered by depth (deepest first)
- `tree_depths: dict` — x → depth until root
- `obstruction: str` — classification string

#### `DualViewSeed(k: int, target_a: int, g: int = 5)`
2-adic Newton projector on the exponent space as a position-dependent butterfly seed.

**Methods**:
- `newton_step_e(e)` — single Newton update on exponent e (mod N)
- `build_position_dependent_seeds()` — list-of-lists of 2×2 complex seeds for each butterfly stage
- `thermodynamic_signature()` — dict with spectral_radius, is_unitary, etc.
- `solvability_report()` — Lie algebra solvability analysis (metabelian, depth 2)

#### `dual_view_qasm_emitter(k: int, target_a: int, p_clean: Optional[int] = None) -> str`
Generate OpenQASM 2.0 for the full dual-view Newton pipeline:
1. State preparation via Hensel bootstrap
2. QFT on exponent register
3. Newton diagonal phase accumulation
4. Inverse QFT
5. Valuation guard (cliff detector)
6. Measurement

If `p_clean` is provided, circuit is vacuum-optimised.

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
- `lift_root(a, p, k)` — Hensel lift cube root from mod p to mod p^k
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

## Iwasawa Algebra Module (`dual_view.iwasawa_algebra`)

The Iwasawa algebra `Z₂[[G]]` represented as power series in the topological generator `(1-γ)`. Provides classification of shift-covariant differential operators via the augmentation ideal.

### Classes

#### `IwasawaElement(coeffs: List[int], precision: int)`
An element of the Iwasawa algebra `Z₂[[G]]` as a power series `μ = Σ cₙ (1-γ)ⁿ`.

**Classmethods**:
- `from_generator(precision=32)` — canonical `1-γ` generator of the augmentation ideal
- `unit(precision=32)` — multiplicative identity
- `zero(precision=32)` — zero element

**Properties / methods**:
- `valuation()` — 2-adic valuation: smallest `n` with `cₙ` odd; returns `inf` if all even
- `is_unit()` — constant term is odd
- `is_generator_of_aug_ideal()` — checks `valuation()≥1` and `c₁` odd
- `truncation_status()` — metadata: `original_degree`, `truncation_degree`, `precision`

**Arithmetic**:
- `__add__` — coefficient-wise addition mod `2^precision`
- `__mul__` — Cauchy product (power series multiplication in `Z₂[[G]]`)

#### `IwasawaAlgebra`
Factory and classification operations.

**Static methods**:
- `aug_ideal_generator(precision=32)` — returns `IwasawaElement.from_generator(precision)`
- `classify_dirac_operator(mu)` — classify `μ` as a shift-covariant differential operator. Returns `{'is_valid', 'is_unit_multiple', 'generator_form'}`. Valid iff `μ` is a unit multiple of `(1-γ)`.

#### `ProModule(name: str, dimension: int, precision: int)`
A profinite `Z₂`-module with valuation filtration.

**Methods**:
- `truncate(k)` — apply truncation functor `Tₖ`, returns new `ProModule`
- `valuation(v)` — compute filtration level of element `v`

### Mersenne (`dual_view.mersenne`)
- `mersenne_coordinates(n, k)` — `(α, e_true, v₂(e_true))`
- `verify_core_identity(n_max)` — 5^(2^(n-2)) ≡ 1 - 2^n
- `mersenne_cliff_table(n_max)` — cliff thresholds
- `bootstrap_cost(eprec0, k)` — Viglietta bit-cost
- `optimal_bootstrap(k_values)` — minimiser search
- `compare_bootstrap_strategies(k_values)` — sqrt vs k/2 vs LUT
- `dlog_with_lut(a, k, b=8)` — LUT-based dlog
- `verify_lut_dlog(k, b=8, n_trials)` — correctness check
- `cliff_constant(g, k)` — compute `c = v₂(log₂(g)/4 + 1)`
- `cliff_formula(g)` — human-readable c(g) formula
- `mersenne_cliff_theorem(verbose)` — state and verify the full theorem
- `prove_cliff_constant(verbose)` — prove `c=5` from 4 log-series terms
- `prove_c_formula(verbose)` — prove `c(g) = v₂(g-5) - 2`
- `exp2_neg4(k)` — compute `exp₂(-4) mod 2^k`, the zero of `log₂(g)/4+1`
- `cliff_constant_unified(g, k)` — unified formula via Newton-Taylor lemma
- `verify_unified_formula(g_values, k)` — verify unified matches direct
- `proof_connection(verbose)` — show all proofs connected via `log₂(5) ≡ -4 (mod 128)`

### Isometry (`dual_view.isometry`)
- `verify_isometry(k, n_trials)` — v₂(5^e-1) = v₂(e)+2
- `isometry_pair_test(k, n_trials)` — pair form
- `isometry_summary(k)` — full conditioning picture
- `verify_operator_algebra(k_values)` — avg², D·avg, avg·D
- `trace_alpha_independence(k, p, ...)` — chi-square test
- `trace_exponent_independence(k, p, ...)` — ANOVA F-test
- `exponent_valuation_profile(k, n_samples)` — v₂(e_true) distribution

## Bridge Module (`dual_view.bridge`)

Unifies butterfly_v2-style spectral seed analysis with tidal-coordinate 2-adic decomposition for neural network weight matrices. Provides three complementary seeds per layer:

- **S1 (Depth Histogram)**: Circulant companion of the depth histogram deviation from a geometric null distribution
- **S2 (Map Seed)**: The weight matrix (or its symmetrised form) as a linear operator
- **S3 (Sign Seed)**: The C₂ butterfly factor capturing α-sector bias

When S1 and S2 agree, the thermodynamic character is unambiguous. When they disagree, the split itself is the finding.

### Functions

#### `quantize(W: np.ndarray, k: int) -> np.ndarray`
Symmetric min-max quantisation: float tensor to `Z/2^k Z`.

### Classes

#### `SpectralThermodynamics`
Spectral eigenvalue analysis of seed matrices. Classifies a matrix's thermodynamic character via its eigenvalue spectrum.

**Classmethod**:
- `analyze(S, tol=1e-10)` — analyse matrix `S`, returns `SpectralThermodynamics`

**Attributes**: `spectral_radius`, `min_eigenvalue_magnitude`, `max_eigenvalue_magnitude`, `entropy_rate`, `is_unitary`, `is_conservative`, `is_contractive`, `is_expansive`, `is_nilpotent`, `lyapunov_exponent`

**Methods**:
- `character()` — returns one of: `NILPOTENT`, `CONSERVATIVE/UNITARY`, `CONSERVATIVE`, `CONTRACTIVE`, `EXPANSIVE`, `MIXED`

#### `ButterflyBridge(k: int = 8)`
Main entry point for unified butterfly + tidal analysis.

**Methods**:
- `analyse_layer(W_float, name="unnamed")` — analyse a single weight matrix, returns `LayerReport`
- `analyse_model(layers: Dict[str, np.ndarray])` — analyse a dict of named weight matrices, returns `ModelReport`

#### `LayerReport` (dataclass)
Per-layer analysis report.

**Attributes**: `name`, `shape`, `k`, `n`, `mean_v`, `H_val`, `zero_frac`, `alpha_frac`, `depth_hist`, `depth_dev`, `depth_char`, `thermo_S1`, `thermo_S2`, `thermo_S3`

**Methods**:
- `consensus()` — compare S1 and S2 character: `NEUTRAL`, `CONSENSUS: ...`, or `SPLIT`
- `report()` — formatted ASCII report

#### `ModelReport` (dataclass)
Multi-layer report aggregating per-layer analyses.

**Methods**:
- `add(r: LayerReport)` — add a layer
- `boundaries()` — list S2 character transitions between consecutive layers
- `report()` — formatted multi-layer report with trajectory table, boundaries, and split layers

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

## Newton Dynamics Module (`dual_view.newton_dynamics`)

p-adic Newton dynamics for the rational map `N(x) = (2x³+1)/(3x²)`.

### Polynomial Arithmetic (`poly.py`)

#### `poly_mul(p: List[int], q: List[int]) -> List[int]`
Polynomial multiplication (coefficient-list convolution).

#### `poly_add(p: List[int], q: List[int]) -> List[int]`
Polynomial addition.

#### `poly_scalar_mul(p: List[int], c: int) -> List[int]`
Polynomial scalar multiplication.

#### `poly_pow(p: List[int], n: int) -> List[int]`
Polynomial power.

#### `poly_divmod(dividend: List[int], divisor: List[int]) -> Tuple[List[int], List[int]]`
Polynomial division. Returns `(quotient, remainder)`.

### Newton Iterates (`iterates.py`)

#### `mobius(n: int) -> int`
Möbius function μ(n).

#### `compute_iterates(max_k: int) -> List[Tuple[List[int], List[int]]]`
Compute (A_d, B_d) for d = 0, …, max_k for the Newton map recurrence:
- `A_{d+1} = 2·A_d³ + B_d³`
- `B_{d+1} = 3·A_d²·B_d`
- `A_0 = x`, `B_0 = 1`

### Dynatomic Polynomial (`dynatomic.py`)

#### `dynatomic_polynomial(period: int, iterates: List[Tuple]) -> List[int]`
Compute Φ_n^*(u) in u = x³ via Möbius inversion. Returns coefficient list from constant term to leading term.

### Clean Primes (`clean_primes.py`)

#### `is_cube(a: int, p: int) -> bool`
Check whether a is a cube modulo prime p.

#### `tonelli_shanks(n: int, p: int) -> Optional[int]`
Tonelli–Shanks modular square root. Returns x with x² ≡ n (mod p), or None.

#### `check_quadratic_cube_roots(a, b, c, p) -> bool`
Check whether ax² + bx + c ≡ 0 (mod p) has a root that is a cube modulo p.

### Precomputed Data (`data.py`)

#### `COEFFS_PERIOD4: List[int]`
Period-4 dynatomic coefficients (25 entries, degree 24 in u).

#### `COEFFS_PERIOD5: List[int]`
Period-5 dynatomic coefficients (81 entries, degree 80 in u).

#### `MULTIPLIERS_PERIOD5: List[Tuple[float, float]]`
16 period-5 multipliers as (real, imag) pairs.

#### `load_period6_coefficients(path: str = "") -> List[int]`
Load period-6 coefficients from file. Returns 233 entries, degree 232 in u.

#### `PERIOD6_PREDICTED: Dict`
Predicted invariants for period 6: mu3, deg_u, deg_x, phi_at_0, phi_at_1_power, prod_mu_power.
