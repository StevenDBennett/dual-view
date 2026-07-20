"""
Experiment 3: Spectral thermodynamics on all 15 clean primes.
Build the Newton adjacency from the basin ordering and classify each.
Skip primes whose basin size would make the matrix too large (>10000).
"""
import sys
sys.path.insert(0, "src")

import numpy as np
from dual_view.butterfly_seed import analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES
from dual_view.bridge import SpectralThermodynamics

MAX_MATRIX_SIZE = 5000

print("=" * 72)
print("Experiment 3: Spectral Thermodynamics of Clean Primes")
print("=" * 72)

print(f"\n{'Prime':>8} {'Nroots':>6} {'Nilpot':>8} {'BasinSz':>8} {'Character':>22} {'IsNilpot':>10} {'Lyapunov':>10} {'SpectralR':>10}")
print("-" * 85)

for p in KNOWN_CLEAN_PRIMES:
    prof = analyze_prime(p)
    ordering = prof.basin_ordering
    if not ordering:
        print(f"{p:>8} {len(prof.roots):>6} {prof.nilpotency_index:>8} {0:>8} {'EMPTY':>22}")
        continue
    if len(ordering) > MAX_MATRIX_SIZE:
        print(f"{p:>8} {len(prof.roots):>6} {prof.nilpotency_index:>8} {len(ordering):>8} {'SKIP (matrix too large)':>22}")
        continue
    M = np.zeros((len(ordering), len(ordering)), dtype=np.float64)
    for i in range(len(ordering) - 1):
        M[i, i + 1] = 1.0
    sig = SpectralThermodynamics.analyze(M)
    char = sig.character()
    print(f"{p:>8} {len(prof.roots):>6} {prof.nilpotency_index:>8} {len(ordering):>8} "
          f"{char:>22} {str(sig.is_nilpotent):>10} {sig.lyapunov_exponent:>10.4f} "
          f"{sig.spectral_radius:>10.4f}")

print("\nDone.")
