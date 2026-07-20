"""
Experiment B: Basin depth spectrum for all clean primes.
Histogram depths, compare 1-root vs 3-root.
"""
import sys
sys.path.insert(0, "src")

from collections import Counter
from dual_view.butterfly_seed import analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES

print("=" * 72)
print("Experiment B: Basin Depth Spectrum")
print("=" * 72)

print(f"\n{'Prime':>8} {'Nroots':>6} {'Nilpot':>6} {'BasinSz':>8} {'DepthDist':>40}")
print("-" * 72)

for p in KNOWN_CLEAN_PRIMES:
    prof = analyze_prime(p)
    depths = [prof.tree_depths[x] for x in prof.basin_ordering if prof.tree_depths.get(x, -1) >= 0]
    if not depths:
        continue
    counter = Counter(depths)
    sorted_depths = sorted(counter.items())
    depth_str = " ".join(f"{d}:{c}" for d, c in sorted_depths[:8])
    if len(sorted_depths) > 8:
        depth_str += f" ... (+{len(sorted_depths)-8} more)"
    print(f"{p:>8} {len(prof.roots):>6} {prof.nilpotency_index:>6} {len(depths):>8} {depth_str:>40}")

print("\nDone.")
