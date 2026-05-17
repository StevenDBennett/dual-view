"""Tests for dual_view.iwasawa_algebra — IwasawaElement, IwasawaAlgebra, ProModule."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.iwasawa_algebra import IwasawaElement, IwasawaAlgebra, ProModule


# ============================================================================
# IwasawaElement — construction and basic properties
# ============================================================================

class TestIwasawaElementConstruction(unittest.TestCase):
    def test_zero(self):
        z = IwasawaElement.zero(precision=16)
        self.assertEqual(z.coeffs, [0])
        self.assertEqual(z.precision, 16)
        self.assertEqual(z.valuation(), float('inf'))
        self.assertFalse(z.is_unit())
        self.assertFalse(z.is_generator_of_aug_ideal())

    def test_unit_from_classmethod(self):
        u = IwasawaElement.unit(precision=32)
        self.assertEqual(u.coeffs, [1])
        self.assertTrue(u.is_unit())
        self.assertEqual(u.valuation(), 0)

    def test_from_generator(self):
        gen = IwasawaElement.from_generator(precision=16)
        self.assertEqual(gen.coeffs, [0, 1])
        self.assertEqual(gen.valuation(), 1)
        self.assertFalse(gen.is_unit())

    def test_from_coeffs_strips_trailing_zeros(self):
        a = IwasawaElement([1, 2, 0, 0], precision=16)
        self.assertEqual(a.coeffs, [1, 2])
        self.assertEqual(a.original_degree, 4)
        self.assertEqual(a.truncation_degree, 2)

    def test_all_zeros_becomes_single_zero(self):
        a = IwasawaElement([0, 0, 0], precision=8)
        self.assertEqual(a.coeffs, [0])

    def test_single_non_zero_preserved(self):
        a = IwasawaElement([5], precision=8)
        self.assertEqual(a.coeffs, [5])

    def test_repr_zero(self):
        z = IwasawaElement.zero(precision=8)
        self.assertIn("0", repr(z))

    def test_repr_generator(self):
        gen = IwasawaElement.from_generator(precision=16)
        r = repr(gen)
        self.assertIn("1(1-gamma)", r)
        self.assertIn("original_degree=2", r)

    def test_repr_with_terms(self):
        a = IwasawaElement([3, 0, 5], precision=16)
        r = repr(a)
        self.assertIn("3", r)
        self.assertIn("5(1-gamma)^2", r)


# ============================================================================
# IwasawaElement — valuation
# ============================================================================

class TestIwasawaElementValuation(unittest.TestCase):
    def test_valuation_zero(self):
        self.assertEqual(IwasawaElement.zero(16).valuation(), float('inf'))

    def test_valuation_constant_odd(self):
        self.assertEqual(IwasawaElement([3], 16).valuation(), 0)

    def test_valuation_constant_even(self):
        self.assertEqual(IwasawaElement([4], 16).valuation(), float('inf'))

    def test_valuation_first_odd_at_index(self):
        a = IwasawaElement([2, 4, 1], 16)
        self.assertEqual(a.valuation(), 2)

    def test_valuation_generator(self):
        gen = IwasawaElement.from_generator(16)
        self.assertEqual(gen.valuation(), 1)

    def test_valuation_all_even(self):
        a = IwasawaElement([2, 4, 6], 16)
        self.assertEqual(a.valuation(), float('inf'))


# ============================================================================
# IwasawaElement — unit detection
# ============================================================================

class TestIwasawaElementUnit(unittest.TestCase):
    def test_unit_from_odd_constant(self):
        self.assertTrue(IwasawaElement([1], 16).is_unit())
        self.assertTrue(IwasawaElement([3], 16).is_unit())
        self.assertTrue(IwasawaElement([255], 16).is_unit())

    def test_not_unit_from_even_constant(self):
        self.assertFalse(IwasawaElement([0], 16).is_unit())
        self.assertFalse(IwasawaElement([2], 16).is_unit())
        self.assertFalse(IwasawaElement([4, 1], 16).is_unit())

    def test_unit_with_higher_terms(self):
        self.assertTrue(IwasawaElement([1, 2, 3], 16).is_unit())

    def test_not_unit_from_generator(self):
        self.assertFalse(IwasawaElement.from_generator(16).is_unit())


# ============================================================================
# IwasawaElement — augmentation ideal generator detection
# ============================================================================

class TestIwasawaElementAugIdeal(unittest.TestCase):
    def test_generator_is_in_aug_ideal(self):
        gen = IwasawaElement.from_generator(16)
        self.assertTrue(gen.is_generator_of_aug_ideal())

    def test_constant_is_not_generator(self):
        u = IwasawaElement.unit(16)
        self.assertFalse(u.is_generator_of_aug_ideal())

    def test_zero_is_not_generator(self):
        z = IwasawaElement.zero(16)
        self.assertFalse(z.is_generator_of_aug_ideal())

    def test_val0_is_not_generator(self):
        a = IwasawaElement([1, 1], 16)
        self.assertEqual(a.valuation(), 0)
        self.assertFalse(a.is_generator_of_aug_ideal())

    def test_even_second_coeff_is_not_generator(self):
        a = IwasawaElement([0, 2], 16)
        self.assertEqual(a.valuation(), float('inf'))
        self.assertFalse(a.is_generator_of_aug_ideal())

    def test_truncated_coeffs_not_generator(self):
        a = IwasawaElement([0], 16)
        self.assertFalse(a.is_generator_of_aug_ideal())


# ============================================================================
# IwasawaElement — arithmetic
# ============================================================================

class TestIwasawaElementArithmetic(unittest.TestCase):
    def test_add_two_elements(self):
        a = IwasawaElement([1, 2], 16)
        b = IwasawaElement([3, 4], 16)
        c = a + b
        self.assertEqual(c.coeffs, [4, 6])

    def test_add_with_broadcasting(self):
        a = IwasawaElement([1], 16)
        b = IwasawaElement([2, 3], 16)
        c = a + b
        self.assertEqual(c.coeffs, [3, 3])

    def test_add_with_mod_reduction(self):
        mod = 1 << 4
        a = IwasawaElement([mod - 1, 1], 4)
        b = IwasawaElement([1, 0], 4)
        c = a + b
        self.assertEqual(c.coeffs, [0, 1])

    def test_mul_identity(self):
        u = IwasawaElement.unit(16)
        a = IwasawaElement([2, 3], 16)
        prod = a * u
        self.assertEqual(prod.coeffs, [2, 3])

    def test_mul_by_zero(self):
        z = IwasawaElement.zero(16)
        a = IwasawaElement([2, 3], 16)
        prod = a * z
        self.assertEqual(prod.coeffs, [0])

    def test_mul_generator_square(self):
        gen = IwasawaElement.from_generator(16)
        sq = gen * gen
        # (1-gamma)^2 = 1 - 2gamma + gamma^2 = 0 + 0(1-g) + 1(1-g)^2
        self.assertEqual(sq.coeffs, [0, 0, 1])

    def test_mul_generator_by_unit(self):
        gen = IwasawaElement.from_generator(16)
        u = IwasawaElement([3], 16)
        prod = u * gen
        self.assertEqual(prod.coeffs, [0, 3])

    def test_mul_cauchy_product(self):
        a = IwasawaElement([1, 1], 16)
        b = IwasawaElement([1, 1], 16)
        # (1 + (1-g))^2 = 1 + 2(1-g) + (1-g)^2
        prod = a * b
        self.assertEqual(prod.coeffs, [1, 2, 1])

    def test_mul_respects_precision_bound(self):
        a = IwasawaElement([1] * 20, 4)
        b = IwasawaElement([1] * 20, 4)
        prod = a * b
        self.assertLessEqual(len(prod.coeffs), 4)


# ============================================================================
# IwasawaElement — truncation metadata
# ============================================================================

class TestIwasawaElementTruncation(unittest.TestCase):
    def test_truncation_status_preserves_original_degree(self):
        a = IwasawaElement([1, 0, 0, 2], 16)
        status = a.truncation_status()
        self.assertIn("original_degree=4", status)
        self.assertIn("truncation_degree=4", status)
        self.assertIn("precision=2^16", status)

    def test_truncation_status_with_trailing_zeros(self):
        a = IwasawaElement([1, 2, 0, 0], 16)
        status = a.truncation_status()
        self.assertIn("original_degree=4", status)
        self.assertIn("truncation_degree=2", status)


# ============================================================================
# IwasawaAlgebra — factory and classification
# ============================================================================

class TestIwasawaAlgebra(unittest.TestCase):
    def test_aug_ideal_generator(self):
        gen = IwasawaAlgebra.aug_ideal_generator(precision=16)
        self.assertIsInstance(gen, IwasawaElement)
        self.assertTrue(gen.is_generator_of_aug_ideal())
        self.assertEqual(gen.valuation(), 1)
        self.assertIn("truncation_degree=2", gen.truncation_status())

    def test_classify_valid_generator(self):
        gen = IwasawaAlgebra.aug_ideal_generator(16)
        result = IwasawaAlgebra.classify_dirac_operator(gen)
        self.assertTrue(result['is_valid'])
        self.assertTrue(result['is_unit_multiple'])
        self.assertIsNotNone(result['generator_form'])

    def test_classify_unit_multiple_of_generator(self):
        gen = IwasawaElement.from_generator(16)
        u = IwasawaElement([3], 16)
        D = u * gen
        result = IwasawaAlgebra.classify_dirac_operator(D)
        self.assertTrue(result['is_valid'])
        self.assertTrue(result['is_unit_multiple'])

    def test_classify_constant_fails(self):
        u = IwasawaElement.unit(16)
        result = IwasawaAlgebra.classify_dirac_operator(u)
        self.assertFalse(result['is_valid'])

    def test_classify_zero_fails(self):
        z = IwasawaElement.zero(16)
        result = IwasawaAlgebra.classify_dirac_operator(z)
        self.assertFalse(result['is_valid'])

    def test_classify_ideal_element_not_generator(self):
        a = IwasawaElement([0, 2], 16)  # second coeff even → not a generator
        result = IwasawaAlgebra.classify_dirac_operator(a)
        self.assertFalse(result['is_valid'])


# ============================================================================
# ProModule
# ============================================================================

class TestProModule(unittest.TestCase):
    def test_construction(self):
        m = ProModule("V_reg", dimension=100, precision=16)
        self.assertEqual(m.name, "V_reg")
        self.assertEqual(m.dimension, 100)
        self.assertEqual(m.precision, 16)

    def test_truncate(self):
        m = ProModule("V", 10, 16)
        t = m.truncate(8)
        self.assertEqual(t.name, "V^(8)")
        self.assertEqual(t.precision, 8)
        self.assertEqual(t.dimension, 10)

    def test_truncate_preserves_dimension(self):
        m = ProModule("V", 10, 16)
        t = m.truncate(8)
        self.assertEqual(t.dimension, m.dimension)

    def test_valuation_all_zeros_returns_precision(self):
        m = ProModule("V", 5, 16)
        val = m.valuation([0, 0, 0, 0, 0])
        self.assertEqual(val, 16)

    def test_valuation_with_mixed_zeros(self):
        m = ProModule("V", 5, 16)
        val = m.valuation([16, 0, 0, 0, 0])
        self.assertEqual(val, 4)

    def test_graded_initialized_empty(self):
        m = ProModule("V", 10, 16)
        self.assertEqual(m.graded, {})


# ============================================================================
# Smoke test: end-to-end classification pipeline
# ============================================================================

class TestClassificationPipeline(unittest.TestCase):
    def test_classify_1_minus_gamma(self):
        """Canonical Dirac operator is valid."""
        D = IwasawaAlgebra.aug_ideal_generator(32)
        result = IwasawaAlgebra.classify_dirac_operator(D)
        self.assertTrue(result['is_valid'])
        self.assertEqual(D.coeffs, [0, 1])

    def test_classify_scaled_by_unit(self):
        """Any unit multiple of (1-gamma) is also a valid Dirac operator."""
        D = IwasawaAlgebra.aug_ideal_generator(32)
        for unit_val in [1, 3, 5, 7, 9]:
            U = IwasawaElement([unit_val], 32)
            D_scaled = U * D
            result = IwasawaAlgebra.classify_dirac_operator(D_scaled)
            self.assertTrue(result['is_valid'],
                            f"Unit {unit_val} * (1-gamma) should be valid")

    def test_bad_element_not_classified(self):
        """Elements not of the form unit * (1-gamma) are rejected."""
        bad_cases = [
            IwasawaElement([1], 32),           # constant unit
            IwasawaElement([0], 32),           # zero
            IwasawaElement([0, 2], 32),        # even second coeff
            IwasawaElement([2, 2], 32),        # valuation 0 (even constant)
            IwasawaElement([0, 1, 2], 32),     # has higher terms, still valid
        ]
        for bad in bad_cases:
            result = IwasawaAlgebra.classify_dirac_operator(bad)
            # The last one should actually be valid (it IS a generator)
            if bad.coeffs == [0, 1, 2]:
                self.assertTrue(result['is_valid'])
            else:
                self.assertFalse(result['is_valid'],
                                f"Expected invalid for {bad}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
