"""
Experiment 2: Mersenne cliff at n=16 (Priority 1 from research_opportunities.md).
The secondary correction ε(n) = v2(n) - 1 predicts that at n=16 (a power of 2),
the cliff should be k* = n + 2 + (v2(16) - 1) = 16 + 2 + (4 - 1) = 21.
"""
import sys
sys.path.insert(0, "src")

from dual_view.mersenne import mersenne_cliff_table, cliff_constant

c = cliff_constant(g=5)
table = mersenne_cliff_table(n_max=16)

print("=" * 72)
print("Experiment 2: Mersenne Cliff at n=16")
print("=" * 72)

print(f"\nUnified cliff constant c(5) = {c}")
print("\nPrediction: k* = n + 2 + v2(n) - 1")
print("  For n=16: k* = 16 + 2 + (4 - 1) = 21\n")

print(f"{'n':>4} {'k*':>4} {'k_pred':>7} {'v2(e)':>6} {'Power2?':>8} {'Match':>10}")
print("-" * 45)
for row in table:
    n = row["n"]
    k_star = row["k*"]
    k_pred = row["k_pred"]
    v2_e = row["v2_e"]
    is_pow2 = n & (n - 1) == 0
    expected_base = n + 2
    if is_pow2:
        correction = n.bit_length() - 2
        expected = n + 2 + correction
    else:
        expected = expected_base
    match = "OK" if k_star == expected else "MISMATCH"
    print(f"{n:>4} {k_star:>4} {k_pred:>7} {v2_e:>6} {str(is_pow2):>8} {match:>10}")

print(f"\nExtended to n=16: {'PASS' if any(r['n']==16 for r in table) else 'FAIL - n=16 row missing'}")
for r in table:
    if r['n'] == 16:
        print(f"  n=16: k*={r['k*']}, expected=21, match={r['k*']==21}")

print("\nDone.")
