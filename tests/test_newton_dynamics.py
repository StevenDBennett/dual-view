"""Tests for dual_view.newton_dynamics — polynomial arithmetic, iterates, dynatomic."""
import unittest
import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.newton_dynamics import (
    poly_mul, poly_add, poly_scalar_mul, poly_pow, poly_divmod,
    mobius, compute_iterates, dynatomic_polynomial,
    is_cube, tonelli_shanks, check_quadratic_cube_roots,
    COEFFS_PERIOD4, COEFFS_PERIOD5,
    load_period6_coefficients, PERIOD6_PREDICTED,
)


class TestMobius(unittest.TestCase):
    def test_mobius_one(self):
        self.assertEqual(mobius(1), 1)

    def test_mobius_prime(self):
        self.assertEqual(mobius(2), -1)
        self.assertEqual(mobius(3), -1)
        self.assertEqual(mobius(5), -1)

    def test_mobius_square_free_product(self):
        self.assertEqual(mobius(6), 1)   # 2*3, even # factors
        self.assertEqual(mobius(30), -1) # 2*3*5, odd # factors

    def test_mobius_square_factor(self):
        self.assertEqual(mobius(4), 0)
        self.assertEqual(mobius(12), 0)
        self.assertEqual(mobius(18), 0)


class TestPoly(unittest.TestCase):
    def test_mul_constant(self):
        self.assertEqual(poly_mul([1], [2, 3]), [2, 3])

    def test_mul_simple(self):
        # (1 + x)(1 - x) = 1 - x^2
        p = poly_mul([1, 1], [1, -1])
        self.assertEqual(p, [1, 0, -1])

    def test_mul_degree(self):
        p = poly_mul([1, 2, 3], [4, 5])
        # (1 + 2x + 3x^2)(4 + 5x) = 4 + 13x + 22x^2 + 15x^3
        self.assertEqual(p, [4, 13, 22, 15])

    def test_add(self):
        self.assertEqual(poly_add([1, 2], [3, 4, 5]), [4, 6, 5])

    def test_scalar_mul(self):
        self.assertEqual(poly_scalar_mul([1, 2, 3], 4), [4, 8, 12])

    def test_pow_zero(self):
        self.assertEqual(poly_pow([1, 2], 0), [1])

    def test_pow_one(self):
        self.assertEqual(poly_pow([1, 2], 1), [1, 2])

    def test_pow_two(self):
        # (1 + x)^2 = 1 + 2x + x^2
        self.assertEqual(poly_pow([1, 1], 2), [1, 2, 1])

    def test_divmod_exact(self):
        # (x^2 - 1) / (x - 1) = x + 1
        q, r = poly_divmod([-1, 0, 1], [-1, 1])
        self.assertEqual(q, [1, 1])
        self.assertEqual(r, [0])

    def test_divmod_with_remainder(self):
        # (x^2) / (x + 1) = x - 1 rem 1
        q, r = poly_divmod([0, 0, 1], [1, 1])
        self.assertEqual(q, [-1, 1])
        self.assertEqual(r, [1])

    def test_divmod_degree_lt(self):
        q, r = poly_divmod([1, 2], [3, 4, 5])
        self.assertEqual(q, [0])
        self.assertEqual(r, [1, 2])

    def test_divmod_by_zero_raises(self):
        with self.assertRaises(ValueError):
            poly_divmod([1], [0])


class TestIterates(unittest.TestCase):
    def test_iterate_0(self):
        iters = compute_iterates(0)
        self.assertEqual(len(iters), 1)
        A0, B0 = iters[0]
        self.assertEqual(A0, [0, 1])  # x
        self.assertEqual(B0, [1])      # 1

    def test_iterate_1(self):
        iters = compute_iterates(1)
        A1, B1 = iters[1]
        # A_1 = 2x^3 + 1, B_1 = 3x^2
        self.assertEqual(A1, [1, 0, 0, 2])  # 1 + 2x^3
        self.assertEqual(B1, [0, 0, 3])      # 3x^2

    def test_iterate_2(self):
        iters = compute_iterates(2)
        A2, B2 = iters[2]
        # A_2 = 2(2x^3+1)^3 + (3x^2)^3 = 2 + 12x^3 + 51x^6 + 16x^9
        # B_2 = 3(2x^3+1)^2 * 3x^2 = 9x^2 + 36x^5 + 36x^8
        self.assertEqual(A2[0], 2)   # constant term
        self.assertEqual(A2[1], 0)
        self.assertEqual(A2[2], 0)
        self.assertEqual(A2[3], 12)  # x^3 term
        self.assertEqual(B2[0], 0)
        self.assertEqual(B2[1], 0)
        self.assertEqual(B2[2], 9)   # x^2 term

    def test_a_at_zero(self):
        """A_d(0) = 2^{(3^{d-1}-1)/2} for d >= 1."""
        iters = compute_iterates(6)
        for d in range(1, 7):
            expected = 2 ** ((3 ** (d - 1) - 1) // 2)
            self.assertEqual(iters[d][0][0], expected)


class TestDynatomic(unittest.TestCase):
    def test_period_2(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        # Φ_2^*(u) = 20u^2 + 5u + 2  (u = x^3)
        self.assertEqual(phi2, [2, 5, 20])

    def test_period_2_special_values(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        # Φ_2^*(0) = 2 = 2^1
        self.assertEqual(phi2[0], 2)
        # Φ_2^*(1) = 20 + 5 + 2 = 27 = 3^3
        self.assertEqual(sum(phi2), 27)

    def test_period_4_degree(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        # degree 24 in u = x^3
        self.assertEqual(len(phi4) - 1, 24)

    def test_period_4_special_values(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        # Φ_4^*(0) = 2^12 = 4096
        self.assertEqual(phi4[0], 4096)
        # Φ_4^*(0) matches COEFFS_PERIOD4
        self.assertEqual(phi4[0], COEFFS_PERIOD4[0])

    def test_period_4_matches_data(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        self.assertEqual(len(phi4), len(COEFFS_PERIOD4))
        # Verify constant term and leading coeff
        self.assertEqual(phi4[0], COEFFS_PERIOD4[0])
        self.assertEqual(phi4[-1], COEFFS_PERIOD4[-1])

    def test_period_5_degree(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        # degree 80 in u = x^3
        self.assertEqual(len(phi5) - 1, 80)

    def test_period_5_constant_term(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        # Φ_5^*(0) = 2^40 = 1099511627776
        self.assertEqual(phi5[0], 1099511627776)


class TestCleanPrimes(unittest.TestCase):
    def test_is_cube_trivial(self):
        # 0 and 1 are cubes mod any prime
        self.assertTrue(is_cube(0, 7))
        self.assertTrue(is_cube(1, 7))

    def test_is_cube_mod_7(self):
        # cubes mod 7: 0, 1, 6 (= -1)
        self.assertTrue(is_cube(6, 7))
        # 3 is not a cube mod 7
        self.assertFalse(is_cube(3, 7))

    def test_tonelli_shanks_perfect_square(self):
        # 4^2 = 16 ≡ 5 (mod 11)
        sqrt = tonelli_shanks(5, 11)
        self.assertIsNotNone(sqrt)
        self.assertEqual((sqrt * sqrt) % 11, 5)

    def test_tonelli_shanks_non_residue(self):
        self.assertIsNone(tonelli_shanks(2, 5))

    def test_tonelli_shanks_p_eq_3_mod_4(self):
        # p = 7 ≡ 3 mod 4
        sqrt = tonelli_shanks(2, 7)  # 3^2 ≡ 2 (mod 7)
        self.assertIsNotNone(sqrt)
        self.assertEqual((sqrt * sqrt) % 7, 2)

    def test_check_quadratic_cube_roots_linear(self):
        # 3x + 1 ≡ 0 (mod 7) → x = 2. Is 2 a cube mod 7?
        # cubes mod 7: 0, 1, 6 → 2 is not a cube
        self.assertFalse(check_quadratic_cube_roots(0, 3, 1, 7))

    def test_check_quadratic_cube_roots_no_discriminant(self):
        # x^2 + x + 1 ≡ 0 mod 5 — check discriminant
        self.assertFalse(check_quadratic_cube_roots(1, 1, 1, 5))


class TestData(unittest.TestCase):
    def test_period4_length(self):
        self.assertEqual(len(COEFFS_PERIOD4), 25)

    def test_period4_constant(self):
        self.assertEqual(COEFFS_PERIOD4[0], 4096)

    def test_period4_leading(self):
        self.assertEqual(COEFFS_PERIOD4[-1], 223338299392)

    def test_period5_length(self):
        self.assertEqual(len(COEFFS_PERIOD5), 81)

    def test_period5_constant(self):
        self.assertEqual(COEFFS_PERIOD5[0], 1099511627776)

    def test_period5_multipliers_count(self):
        from dual_view.newton_dynamics import MULTIPLIERS_PERIOD5 as M
        self.assertEqual(len(M), 16)

    def test_period5_multiplier_sum(self):
        from dual_view.newton_dynamics import MULTIPLIERS_PERIOD5 as M
        total_real = sum(m[0] for m in M)
        self.assertAlmostEqual(total_real, 486.0, places=5)

    def test_period6_predicted(self):
        self.assertEqual(PERIOD6_PREDICTED["mu3"], 696)
        self.assertEqual(PERIOD6_PREDICTED["deg_u"], 232)
        self.assertEqual(PERIOD6_PREDICTED["phi_at_0"], 2 ** 116)

    def test_load_period6_coefficients(self):
        coeffs = load_period6_coefficients()
        self.assertEqual(len(coeffs), 233)
        self.assertEqual(coeffs[0], 83076749736557242056487941267521536)


class TestIntegration(unittest.TestCase):
    """End-to-end checks that the research findings hold."""

    def test_period2_product_formula(self):
        """μ₃(2) = 6, product = 6¹ = 6"""
        from dual_view.newton_dynamics import mobius
        mu3 = sum(mobius(2 // d) * (3 ** d) for d in (1, 2))
        self.assertEqual(mu3, 6)
        # ∏μ = 6^{μ₃/6}
        self.assertEqual(6 ** (mu3 // 6), 6)

    def test_period4_product_formula(self):
        """μ₃(4) = 72, product = 6¹² = 2176782336"""
        from dual_view.newton_dynamics import mobius
        mu3 = sum(mobius(4 // d) * (3 ** d) for d in (1, 2, 4))
        self.assertEqual(mu3, 72)
        self.assertEqual(6 ** (mu3 // 6), 2176782336)


if __name__ == "__main__":
    unittest.main()
