"""
Experiment D: KroneckerCliffScorer with clean-prime structured factors.
Build synthetic matrices using basin depth distributions as weight templates.
"""
import sys
sys.path.insert(0, "src")

import numpy as np
from dual_view.butterfly_seed import analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES
from dual_view.butterfly import KroneckerCliffScorer

print("=" * 72)
print("Experiment D: Kronecker Cliff Scores for Clean-Prime Structured Factors")
print("=" * 72)

for p in KNOWN_CLEAN_PRIMES[:6]:  # first 6 primes only (runs fast)
    prof = analyze_prime(p)
    depths = [prof.tree_depths[x] for x in prof.basin_ordering if prof.tree_depths.get(x, -1) >= 0]
    if not depths:
        continue
    n = min(len(depths), 64)
    depths_arr = np.array(depths[:n], dtype=np.float64)
    depths_arr = depths_arr.reshape(1, -1)
    factors = [depths_arr, depths_arr.T]
    scorer = KroneckerCliffScorer(factors, k_range=range(4, 10))
    try:
        results = scorer.score_factors()
        report = scorer.print_report()
        print(f"\np={p}:")
        print(report)
    except Exception as e:
        print(f"\np={p}: ERROR {e}")

print("\nDone.")
