"""
demo.py
-------
Runnable demonstration of the dual-view package.

Usage:
    python -m dual_view.demo          # full demo
    python -m dual_view.demo --quick   # skip slow tests
"""
from __future__ import annotations

import argparse
import time
import sys


def demo_core(quick: bool = False) -> None:
    """Core arithmetic: DualNumber roundtrip, modinv, processor."""
    print("\n=== Core Arithmetic ===")
    from .core import run_all_tests, DualNumber, TwoAdicProcessor, _mask

    run_all_tests(k=16)

    proc = TwoAdicProcessor(32)
    a = DualNumber(42, 32)
    b = DualNumber(7, 32)
    c = proc.mul(a, b)
    print(f"  42 × 7 = {c.value}  (expected {42 * 7 & _mask(32)})")

    inv7 = proc.inv(b)
    check = proc.mul(b, inv7)
    print(f"  7 × 7⁻¹ = {check.value}  (expected 1)")
    print("  PASS" if c.value == (42 * 7 & _mask(32)) else "  FAIL")


def demo_operators(quick: bool = False) -> None:
    """Operator algebra: eigenfunction identity."""
    print("\n=== Operator Algebra ===")
    from .exponent import ExponentSpace

    es = ExponentSpace(5, 8)
    f = lambda e: pow(5, e, 1 << 8)
    all_eigen = all(es.is_eigenfunction(f, e) for e in range(es.N))
    print(f"  D(g^e) = (g-1)·g^e  for all e: {all_eigen}")
    print("  PASS" if all_eigen else "  FAIL")


def demo_projector(quick: bool = False) -> None:
    """Newton projector convergence."""
    print("\n=== Newton Projector ===")
    from .operators import OperatorContext, NewtonProjector

    ctx = OperatorContext(8, 5)
    proj = NewtonProjector(ctx, pow(5, 7, 1 << 8))
    e = proj.project_point(0)
    print(f"  Project seed 0 → e = {e}  (expected 7)")
    print("  PASS" if e == 7 else "  FAIL")


def demo_precision_sweep(quick: bool = False) -> None:
    """Precision sweep: find ghost cliff."""
    print("\n=== Precision Sweep ===")
    from .basin import precision_sweep

    results = precision_sweep(4, 10, 5, 3)
    for k, frac in results:
        marker = " ← GHOST!" if frac > 0 else ""
        print(f"  k={k:2d}: ghost fraction = {frac:.3f}{marker}")
    print("  Done")


def demo_gauge(quick: bool = False) -> None:
    """Gauge invariants."""
    print("\n=== Gauge Invariants ===")
    from .gauge import GaugeLayer

    gl = GaugeLayer([3, 5, 7, 9], k=16)
    print(f"  {gl.report()}")
    print("  PASS" if gl.tidal is not None else "  FAIL")


def demo_crt_stability(quick: bool = False) -> None:
    """CRT combined stability."""
    print("\n=== CRT Stability ===")
    from .crt import combined_stability

    result = combined_stability(k=6, p=7, num_cycles=20)
    print(f"  Pearson r = {result['pearson_r']:.3f}  (n={result['n_samples']})")
    print("  Done")


def demo_ramp_break(quick: bool = False) -> None:
    """Non-Abelian ramp-break strength."""
    if quick:
        print("\n=== Non-Abelian Ramp Break (skipped: --quick) ===")
        return
    print("\n=== Non-Abelian Ramp Break ===")
    import warnings
    from .nonabelian import ramp_break_strength

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = ramp_break_strength(k=6, p=7, num_cycles=10)
    print(f"  Phase alignment: {result['phase_alignment']:.1%}")
    print(f"  Mean conv ratio: {result['mean_conv']:.3f}")
    print("  Done")


def demo_mersenne_cliff(quick: bool = False) -> None:
    """Mersenne cliff table."""
    print("\n=== Mersenne Cliff Table ===")
    from .mersenne import mersenne_cliff_table

    rows = mersenne_cliff_table(n_max=10)
    print(f"  {'n':>3} {'k*':>4} {'k_pred':>6} {'match':>6}")
    print("  " + "-" * 22)
    for r in rows:
        match = "✓" if r['k*'] == r['k_pred'] else "✗"
        print(f"  {r['n']:>3} {r['k*']:>4} {r['k_pred']:>6} {match:>6}")
    print("  Done")


def demo_thermodynamics(quick: bool = False) -> None:
    """SeedThermodynamics analysis."""
    print("\n=== Thermodynamics ===")
    from .thermodynamics import SeedThermodynamics
    import numpy as np

    W = np.random.randint(0, 256, size=(8, 8)).astype(np.int64)
    st = SeedThermodynamics(k=8, g=5)
    stats = st.analyse(W)
    for key, val in stats.items():
        print(f"  {key}: {val:.4f}")
    print("  Done")


def main() -> None:
    parser = argparse.ArgumentParser(description="dual-view demo suite")
    parser.add_argument("--quick", action="store_true", help="Skip slow tests")
    args = parser.parse_args()

    start = time.time()
    print("dual-view: 2-adic Dual-View Diagnostics for Quantized Weights")

    demo_core(args.quick)
    demo_operators(args.quick)
    demo_projector(args.quick)
    demo_precision_sweep(args.quick)
    demo_gauge(args.quick)
    demo_crt_stability(args.quick)
    demo_mersenne_cliff(args.quick)
    demo_thermodynamics(args.quick)
    demo_ramp_break(args.quick)

    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.2f}s")


if __name__ == "__main__":
    main()
