"""
Experiment 5: Nilpotency index vs prime size analysis.
"""
import sys
sys.path.insert(0, "src")

import math
from dual_view.butterfly_seed import analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES

print("=" * 72)
print("Experiment 5: Nilpotency Index Analysis")
print("=" * 72)

print(f"\n{'Prime':>8} {'p mod 3':>8} {'Nroots':>6} {'Nilpot':>8} {'BasinSz':>8} {'p/Basin':>8} {'log2(p)':>8} {'Nilt/log':>8}")
print("-" * 72)

results = []
for p in KNOWN_CLEAN_PRIMES:
    prof = analyze_prime(p)
    nroots = len(prof.roots)
    nilpot = prof.nilpotency_index
    basin_size = len(prof.basin_ordering)
    p_mod_3 = p % 3
    ratio = p / max(basin_size, 1)
    log2p = math.log2(p)
    nil_per_log = nilpot / max(log2p, 1)
    results.append((p, nroots, nilpot, basin_size, ratio, log2p, nil_per_log))
    print(f"{p:>8} {p_mod_3:>8} {nroots:>6} {nilpot:>8} {basin_size:>8} {ratio:>8.2f} {log2p:>8.2f} {nil_per_log:>8.4f}")

print(f"\nSummary stats:")
nilpots = [r[2] for r in results]
print(f"  Nilpotency: min={min(nilpots)}, max={max(nilpots)}, mean={sum(nilpots)/len(nilpots):.1f}")

# Correlation between p and nilpotency
import numpy as np
ps = [r[0] for r in results]
logps = [math.log2(p) for p in ps]
corr = np.corrcoef(logps, nilpots)[0, 1]
print(f"  Correlation log2(p) vs nilpotency: r={corr:.4f}")

print("\nDone.")
