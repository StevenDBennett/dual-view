"""
Experiment 1: QASM circuit size comparison for all 15 clean primes.
Compare circuit gate counts with and without vacuum optimisation.
"""
import sys
sys.path.insert(0, "src")

from dual_view.butterfly_seed import dual_view_qasm_emitter, analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES

print("=" * 72)
print("Experiment 1: QASM Circuit Size Comparison")
print("=" * 72)

for k in (8, 12, 16):
    print(f"\n--- Precision k={k} ---")
    print(f"{'Prime':>8} {'Nroots':>6} {'Nilpot':>6} {'Std lines':>10} {'Vac lines':>10} {'Reduction':>10}")
    print("-" * 60)
    for p in KNOWN_CLEAN_PRIMES:
        prof = analyze_prime(p)
        qasm_std = dual_view_qasm_emitter(k, 17, p_clean=None)
        qasm_vac = dual_view_qasm_emitter(k, 17, p_clean=p)
        std_lines = len(qasm_std.splitlines())
        vac_lines = len(qasm_vac.splitlines())
        reduction = std_lines - vac_lines
        print(f"{p:>8} {len(prof.roots):>6} {prof.nilpotency_index:>6} {std_lines:>10} {vac_lines:>10} {reduction:>+10}")

print("\nDone.")
