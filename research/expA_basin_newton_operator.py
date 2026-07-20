"""
Experiment A: Build the basin shift operator (non-root elements only).
The Newton map on non-root elements shifts each element to its parent,
forming a nilpotent strictly upper-triangular matrix in the basin ordering.
"""
import sys
sys.path.insert(0, "src")

import numpy as np
from dual_view.butterfly_seed import analyze_prime, _newton_fp
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES

print("=" * 72)
print("Experiment A: Basin Shift Operator — Nilpotent Structure")
print("=" * 72)

print(f"\n{'Prime':>8} {'Nroots':>6} {'NonRoot':>8} {'NilpotIdx':>10} {'UpperTri':>10} {'S^n=0':>10} {'Rank':>8}")
print("-" * 65)

for p in KNOWN_CLEAN_PRIMES:
    prof = analyze_prime(p)
    ordering = prof.basin_ordering
    roots_set = set(prof.roots)
    non_roots = [x for x in ordering if x not in roots_set]
    n = len(non_roots)
    if n > 10000:
        print(f"{p:>8} {len(prof.roots):>6} {n:>8} {prof.nilpotency_index:>10} {'SKIP (too large)':>28}")
        continue
    idx_of = {x: i for i, x in enumerate(non_roots)}
    S = np.zeros((n, n), dtype=np.float64)
    for x in non_roots:
        parent = _newton_fp(x, p)
        if parent is not None and parent in idx_of:
            j = idx_of[x]
            i = idx_of[parent]
            S[i, j] = 1.0
    is_upper = bool(np.allclose(S, np.triu(S)))
    nilpot_check = bool(np.allclose(np.linalg.matrix_power(S, prof.nilpotency_index), 0))
    rank = int(np.linalg.matrix_rank(S))
    print(f"{p:>8} {len(prof.roots):>6} {n:>8} {prof.nilpotency_index:>10} "
          f"{str(is_upper):>10} {str(nilpot_check):>10} {rank:>8}")

print("\nDone.")
