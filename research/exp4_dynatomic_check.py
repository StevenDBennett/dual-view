"""
Experiment 4: Period-2/3/4 dynatomic verification for all 16 clean primes.
"""
import sys
sys.path.insert(0, "src")

from dual_view.newton_dynamics import (
    is_cube, check_quadratic_cube_roots, KNOWN_CLEAN_PRIMES,
    COEFFS_PERIOD4,
)

def eval_poly_mod(coeffs, x, p):
    """Evaluate polynomial at x modulo p (coeffs lowest-degree first)."""
    total = 0
    for i, c in enumerate(coeffs):
        total = (total + c * pow(x, i, p)) % p
    return total

def has_period4_root(p):
    """Check if period-4 polynomial has a cube root modulo p."""
    for u in range(p):
        if eval_poly_mod(COEFFS_PERIOD4, u, p) == 0 and is_cube(u, p):
            return True
    return False

print("=" * 72)
print("Experiment 4: Dynatomic Period-2/3/4 Verification")
print("=" * 72)

print(f"\n{'Prime':>8} {'Period2':>9} {'Period3':>9} {'Pole':>9} {'Period4':>9} {'Status':>12}")
print("-" * 60)

all_clean = True
for p in KNOWN_CLEAN_PRIMES:
    p2 = check_quadratic_cube_roots(20, 5, 2, p)
    p3 = check_quadratic_cube_roots(19, 7, 1, p)
    pole_u = (p - 1) // 2
    pole = is_cube(pole_u, p)
    p4 = has_period4_root(p)
    clean = not (p2 or p3 or p4)
    all_clean = all_clean and clean
    status = "CLEAN" if clean else "HAS PERIODIC POINTS"
    print(f"{p:>8} {str(p2):>9} {str(p3):>9} {str(pole):>9} {str(p4):>9} {status:>12}")

print(f"\nAll 16 clean primes pass dynatomic check: {all_clean}")
print("Done.")
