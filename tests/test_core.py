"""Tests for dual_view.core — DualNumber, modinv, dlog, TwoAdicProcessor."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.core import (
    _mask, _valuation, modinv_newton, two_adic_log5,
    two_adic_dlog, DualNumber, TwoAdicProcessor, run_all_tests,
)


class TestMask(unittest.TestCase):
    def test_mask_identity(self):
        for k in [1, 2, 3, 8, 16, 32]:
            self.assertEqual(_mask(k), (1 << k) - 1)

    def test_mask_range(self):
        for k in range(1, 33):
            m = _mask(k)
            self.assertEqual(m.bit_length(), k)


class TestValuation(unittest.TestCase):
    def test_powers_of_two(self):
        for n in range(10):
            self.assertEqual(_valuation(1 << n), n)

    def test_zero(self):
        self.assertEqual(_valuation(0), float("inf"))

    def test_odd(self):
        for n in [1, 3, 5, 7, 9, 255]:
            self.assertEqual(_valuation(n), 0)

    def test_mixed(self):
        cases = {12: 2, 24: 3, 100: 2, 1024: 10, 1025: 0}
        for n, expected in cases.items():
            self.assertEqual(_valuation(n), expected)


class TestModinvNewton(unittest.TestCase):
    def test_small_k(self):
        for k in [8, 16, 32]:
            for a in [1, 3, 5, 7, 9]:
                inv = modinv_newton(a, k)
                self.assertEqual((a * inv) & _mask(k), 1)

    def test_random(self):
        import random
        for k in [8, 12, 16, 24]:
            for _ in range(20):
                a = random.randrange(1, 1 << k, 2)
                inv = modinv_newton(a, k)
                self.assertEqual((a * inv) & _mask(k), 1)

    def test_even_raises(self):
        with self.assertRaises(ValueError):
            modinv_newton(2, 8)

    def test_k_zero_raises(self):
        with self.assertRaises(ValueError):
            modinv_newton(3, 0)


class TestTwoAdicLog5(unittest.TestCase):
    def test_suffix_stability(self):
        L8 = two_adic_log5(8)
        L16 = two_adic_log5(16)
        self.assertEqual(L16 & _mask(8), L8)

    def test_caching(self):
        L1 = two_adic_log5(12)
        L2 = two_adic_log5(12)
        self.assertIs(L1, L2)


class TestTwoAdicDlog(unittest.TestCase):
    def test_known_values(self):
        for e_true in [0, 1, 2, 3, 7, 15, 31, 63]:
            k = 16
            a = pow(5, e_true, 1 << k)
            result = two_adic_dlog(a, k)
            self.assertIsNotNone(result)
            alpha, e = result
            # a ≡ 1 mod 4 → alpha = 0, e = e_true
            if a & 3 == 1:
                self.assertEqual(e, e_true)
                self.assertEqual(alpha, 0)

    def test_odd_returns_tuple(self):
        for a in [1, 3, 5, 7, 255]:
            result = two_adic_dlog(a, 16)
            self.assertIsNotNone(result)
            alpha, e = result
            self.assertIn(alpha, (0, 1))

    def test_even_returns_none(self):
        for a in [0, 2, 4, 6, 8]:
            self.assertIsNone(two_adic_dlog(a, 16))

    def test_verify_roundtrip(self):
        for a in [1, 3, 5, 7, 255, 1023, 32767]:
            k = 16
            result = two_adic_dlog(a, k)
            if result is not None:
                alpha, e = result
                expected = a
                recomputed = pow(5, e, 1 << k)
                if alpha:
                    recomputed = (-recomputed) & _mask(k)
                self.assertEqual(recomputed, expected)


class TestDualNumber(unittest.TestCase):
    def test_zero(self):
        d = DualNumber(0, 16)
        self.assertTrue(d.is_zero)
        self.assertEqual(d.value, 0)

    def test_one(self):
        d = DualNumber(1, 16)
        self.assertFalse(d.is_zero)
        self.assertEqual(d.value, 1)
        self.assertTrue(d.verify())

    def test_negative_one(self):
        d = DualNumber(-1 & _mask(16), 16)
        self.assertTrue(d.verify())

    def test_power_of_two(self):
        d = DualNumber(8, 16)
        self.assertTrue(d.verify())
        self.assertEqual(d.v, 3)

    def test_overflow(self):
        d = DualNumber(1 << 17, 16)
        self.assertTrue(d.is_zero)

    def test_from_coords_roundtrip(self):
        for v in [0, 1, 2]:
            for alpha in [0, 1]:
                for e in [0, 1, 5, 10]:
                    k = 12
                    d = DualNumber.from_coords(v, alpha, e, k)
                    self.assertTrue(d.verify())
                    self.assertEqual(d.v, v)
                    self.assertEqual(d.alpha, alpha)
                    if not d.is_zero:
                        self.assertEqual(d.e, e % (1 << (k - 2)))

    def test_k_too_small_raises(self):
        with self.assertRaises(ValueError):
            DualNumber(1, 2)


class TestTwoAdicProcessor(unittest.TestCase):
    def setUp(self):
        self.proc = TwoAdicProcessor(16)

    def test_mul(self):
        a = DualNumber(3, 16)
        b = DualNumber(5, 16)
        c = self.proc.mul(a, b)
        self.assertEqual(c.value, 15)

    def test_mul_commutative(self):
        a = DualNumber(7, 16)
        b = DualNumber(11, 16)
        ab = self.proc.mul(a, b)
        ba = self.proc.mul(b, a)
        self.assertEqual(ab.value, ba.value)

    def test_mul_associative(self):
        a = DualNumber(3, 16)
        b = DualNumber(5, 16)
        c = DualNumber(7, 16)
        ab_c = self.proc.mul(self.proc.mul(a, b), c)
        a_bc = self.proc.mul(a, self.proc.mul(b, c))
        self.assertEqual(ab_c.value, a_bc.value)

    def test_inv(self):
        a = DualNumber(3, 16)
        inv_a = self.proc.inv(a)
        prod = self.proc.mul(a, inv_a)
        self.assertEqual(prod.value, 1)

    def test_inv_non_unit_raises(self):
        a = DualNumber(2, 16)
        with self.assertRaises(ValueError):
            self.proc.inv(a)

    def test_pow(self):
        a = DualNumber(3, 16)
        a4 = self.proc.pow(a, 4)
        self.assertEqual(a4.value, 81)

    def test_pow_zero(self):
        a = DualNumber(3, 16)
        a0 = self.proc.pow(a, 0)
        self.assertEqual(a0.value, 1)

    def test_pow_negative(self):
        a = DualNumber(3, 16)
        a_neg = self.proc.pow(a, -1)
        prod = self.proc.mul(a, a_neg)
        self.assertEqual(prod.value, 1)

    def test_mul_zero(self):
        a = DualNumber(3, 16)
        z = DualNumber(0, 16)
        c = self.proc.mul(a, z)
        self.assertTrue(c.is_zero)

    def test_overflow(self):
        a = DualNumber(1 << 14, 16)
        b = DualNumber(4, 16)  # v=2 → v_total = 14 + 2 = 16 → overflow
        c = self.proc.mul(a, b)
        self.assertTrue(c.is_zero)


class TestSynthesisFindings(unittest.TestCase):
    """Synthesis-level correctness tests."""

    def test_exponential_isometry(self):
        k = 12
        for e in range(1, 32):
            diff = (pow(5, e, 1 << k) - 1) & _mask(k)
            v2_diff = _valuation(diff) if diff != 0 else k
            v2_e = _valuation(e)
            self.assertEqual(v2_diff, v2_e + 2)

    def test_lut_bootstrap_correctness(self):
        k = 16
        for e_true in range(64):
            a = pow(5, e_true, 1 << k)
            result = two_adic_dlog(a, k)
            if result is not None:
                alpha, e = result
                if a & 3 == 1:
                    self.assertEqual(e, e_true)


if __name__ == "__main__":
    unittest.main()
