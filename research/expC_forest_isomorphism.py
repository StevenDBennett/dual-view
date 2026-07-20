"""
Experiment C: Cross-prime forest isomorphism.
Compare basin depth profiles to detect structural similarity.
"""
import sys
sys.path.insert(0, "src")

from collections import Counter
from dual_view.butterfly_seed import analyze_prime
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES

print("=" * 72)
print("Experiment C: Forest Isomorphism — Structural Similarity")
print("=" * 72)

profiles = []
for p in KNOWN_CLEAN_PRIMES:
    prof = analyze_prime(p)
    depths = Counter()
    for x, d in prof.tree_depths.items():
        if d >= 0:
            depths[d] += 1
    max_depth = max(depths.keys()) if depths else 0
    n_roots = len(prof.roots)
    n_nodes = sum(depths.values())
    sig = tuple(sorted(depths.items()))
    profiles.append((p, n_roots, max_depth, n_nodes, sig))

print(f"\n{'Prime':>8} {'Nroots':>6} {'MaxDep':>8} {'Nodes':>8} {'DepthSig':>40}")
print("-" * 72)

for p, nr, md, nn, sig in profiles:
    sig_str = " ".join(f"{d}:{c}" for d, c in sig[:10])
    if len(sig) > 10:
        sig_str += f" ... (+{len(sig)-10})"
    print(f"{p:>8} {nr:>6} {md:>8} {nn:>8} {sig_str:>40}")

sig_map = {}
for p, nr, md, nn, sig in profiles:
    sig_map.setdefault(sig, []).append(p)

isomorphisms = {s: ps for s, ps in sig_map.items() if len(ps) > 1}
if isomorphisms:
    print(f"\nISOMORPHIC PAIRS (same depth signature):")
    for sig, primes in isomorphisms.items():
        print(f"  {primes}  sig={dict(sig)}")
else:
    print(f"\nNo isomorphic pairs — all 15 depth signatures are unique.")

print("\nDone.")
