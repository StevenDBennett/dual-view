# Bug Fix History and Audit

## Summary

Across the three predecessor projects, **10 bugs** were identified and fixed. The unified dual-view project incorporates all fixes.

## Critical Bugs

### 1. α=1 Sector Ghost Attractor (BasinExplorer.newton_step)

**Location**: `basin.py` / `BasinExplorer`

**Bug**: For targets with `α = 1` (`a ≡ 3 mod 4`), the Newton iteration solved `5^e ≡ a` instead of `5^e ≡ -a`. This caused convergence to a ghost fixed point at `e* = dlog(a+2, k)` for all α=1 targets.

**Fix**: Added check in `__init__`:
```python
if alpha:
    self.a = (-target) & self.mask
else:
    self.a = target
```

**Impact**: Affected v2 and v3 — both lacked this fix. Without it, ALL weights in the α=1 sector (roughly half of all odd weights) would be classified as ghosts, making the GhostMap binary rather than graded.

**Found**: Original v1 codebase contained the fix; v2/v3 regressed when code was rewritten.

### 2. NewtonProjector Modulus Mismatch

**Location**: `operators.py` / `NewtonProjector._step`

**Bug**: In v2, the modular inverse was computed modulo `k` (full ring) instead of `k-2` (exponent domain): `modinv_newton(df_unit, self.ctx.k)`. Since `df_unit` is masked to `k-2` bits, passing `k` as the modulus was incorrect.

**Fix**: Changed to `modinv_newton(df_unit, self.ctx.k - 2)`.

**Impact**: Only affected v2. V1 and v3 had the correct modulus.

### 3. OperatorContext.g Attribute Missing

**Location**: `operators.py` / `OperatorContext`

**Bug**: The `g` (generator) attribute was computed inside the initialization but not stored as `self.g`, causing `AttributeError` when `NewtonProjector` accessed `ctx.g`.

**Fix**: Added `self._g = g` and `@property def g(self): return self._g`.

**Impact**: v1 original had this bug; fixed in v1 and preserved in v2/v3.

## High-Severity Bugs

### 4. Average Operator O(N²) Performance

**Location**: `operators.py` / `OperatorContext.avg`

**Bug**: The average operator was implemented as `(1/N) Σ Sⁱ`, where the `_Operator.__pow__` chained `N` shift operators, each adding a closure layer. For `N = 2^(k-2)`, this created `O(N²)` closures at construction time, which was fatal for `k ≥ 8` (256 closures → 65,536 compositions).

**Fix**: Replaced with a direct O(N) summation lambda:
```python
def avg_action(f, e, ctx=self):
    total = 0
    for t in range(ctx.N):
        total += f(t)
    return total & ctx.mask
```

### 5. numpy.int64 Crash in _v2()

**Location**: `regularization.py` / `GhostMap`

**Bug**: `_v2(n)` calls `n.bit_length()`, which returns `0` for `numpy.int64` zero (instead of raising). This caused crashes in `GhostMap.ratio()` when called with numpy types.

**Fix**: Added `int()` coercion in all public GhostMap methods: `a_int = int(a)`.

### 6. Deprecated eps_crit Metric in ramp_break_strength

**Location**: `nonabelian.py` / `ramp_break_strength`

**Bug**: The function measured `eps_crit` — the smallest `ε` such that a perturbation `d → d + ε` makes `det` even. For the shear matrix `[[d, 1], [0, 1]]`, the determinant is always `d`, and the smallest `ε` giving an even result is always `1`. The metric is degenerate.

**Fix**: Added `warnings.warn(..., DeprecationWarning)` and replaced with `phase_alignment_experiment()` which measures whether the α-sector of the determinant flips under perturbation.

## Medium-Severity Bugs

### 7. Training History Overwritten Each Epoch

**Location**: `training.py` / `train()`

**Bug**: History lists were overwritten each epoch: `history['loss'] = ...` instead of `history['loss'].append(...)`.

**Fix**: Changed from assignment to `.append()`.

### 8. weight_product mod=1 Bug

**Location**: `gauge.py` / `weight_product` (v1)

**Bug**: For plain integer lists with `mod=1`, the product returned 0 instead of the modular reduction.

**Fix**: Added explicit mod check: `mod = max(mod, 1 << k)`.

### 9. Closures in _Operator.__add__ / __sub__

**Location**: `operators.py` / `_Operator.__add__`, `__sub__`

**Bug**: Closures captured loop variables by reference (late binding), so all composed operators used the final values.

**Fix**: Added default-argument capture: `def action(f, e, s=self, o=other):`.

### 10. Test 11 Was a Tautology

**Location**: Original test suite (pre-dating dual-view)

**Bug**: "Test 11" tested Newton correction uniformity by comparing the first-step correction distribution to a uniform distribution — but the test used the same data to fit the distribution parameters and test against them, making it a tautology.

**Fix**: Proper separation of training and test data.

## Audit Findings (from v1/v2 Audits)

### Boundary Testing (33 pushes)

The boundary testing campaign found:
- **2 new theorems**: global convergence for α=0, ghost fixed point formula
- **4 bug classes**: all fixed
- **3 fatal flaw analyses**: ghost regularization signal, ramp break metric, popcount compression correlation

### Fatal Flaw: Ghost Regularization Signal

The ghost regularization signal (from `GhostMap`) is **binary** after the α=1 fix — every odd weight has convergence ratio 1.0. The genuinely graded stability measure is `v₂(e_true)`, which is the quantity tracked by `SeedThermodynamics`. The original ghost penalty provides no useful gradient for odd weights.

This is why `ghost_penalty()` returns zero gradient for odd weights, and the regularization must use the thermodynamic `v₂(e)` metric instead.

### Phase Alternation in GL(2) Holonomy

The α-sector of the holonomy determinant flips under single-bit perturbations at a rate of ~68% (empirical, 30 trials). This phase alignment signal is the replacement for the degenerate `eps_crit` metric.

## Bug Fix Chronology

| # | Bug | Found | Fixed in | Type |
|---|-----|-------|----------|------|
| 1 | α=1 ghost attractor | v1 original | v1 → dual-view | Critical |
| 2 | NewtonProjector modulus | v2 | v1/v3 → dual-view | Critical |
| 3 | ctx.g missing | v1 | v1 → dual-view | Critical |
| 4 | avg O(N²) | v3 addendum | v3 → dual-view | High |
| 5 | numpy.int64 crash | v3 addendum | v3 → dual-view | High |
| 6 | eps_crit degenerate | v3 addendum | v3 → dual-view | High |
| 7 | History overwrite | v1 audit | v1 → dual-view | Medium |
| 8 | mod=1 bug | v1 | v1 → dual-view | Medium |
| 9 | Closure capture | v1 → v3 | v3 → dual-view | Medium |
| 10 | Test 11 tautology | v3 addendum | v3 audit | Low |

All 10 bugs are fixed in the unified dual-view codebase.
