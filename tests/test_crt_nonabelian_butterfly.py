"""Tests for dual_view.crt, dual_view.nonabelian, dual_view.butterfly."""
import unittest

from dual_view.crt import (
    CRTDualNumber, CRTDualProcessor, combined_stability,
    _primitive_root, _prime_dlog,
)
from dual_view.nonabelian import (
    NonAbelianCRTDual, phase_alignment_experiment,
    _mat_mul, _mat_det,
)
from dual_view.butterfly import KroneckerCliffScorer, semiring_cliff_score
import numpy as np


class TestCRTHelpers(unittest.TestCase):
    def test_primitive_root_known(self):
        roots = {3: 2, 5: 2, 7: 3, 11: 2, 13: 2, 17: 3}
        for p, expected in roots.items():
            self.assertEqual(_primitive_root(p), expected)

    def test_prime_dlog_known(self):
        self.assertEqual(_prime_dlog(1, 7, 3), 0)
        self.assertEqual(_prime_dlog(3, 7, 3), 1)


class TestCRTDualNumber(unittest.TestCase):
    def test_initialization(self):
        n = CRTDualNumber(42, k=6, p=7, g_p=3)
        self.assertIsNotNone(n.component_2)
        self.assertEqual(n.residue_p, 42 % 7)

    def test_verify(self):
        n = CRTDualNumber(42, k=6, p=7, g_p=3)
        self.assertTrue(n.verify())

    def test_zero_component(self):
        n = CRTDualNumber(0, k=6, p=7, g_p=3)
        self.assertTrue(n.component_2.is_zero)


class TestCRTDualProcessor(unittest.TestCase):
    def setUp(self):
        self.proc = CRTDualProcessor(k=6, p=7)

    def test_crt_reconstruct(self):
        r = self.proc.crt_reconstruct(5, 3)
        self.assertEqual(r % (1 << 6), 5)
        self.assertEqual(r % 7, 3)

    def test_mul(self):
        a = CRTDualNumber(3, k=6, p=7, g_p=self.proc.g_p)
        b = CRTDualNumber(5, k=6, p=7, g_p=self.proc.g_p)
        c = self.proc.mul(a, b)
        self.assertTrue(c.verify())
        self.assertEqual(c.component_2.value, (3 * 5) & ((1 << 6) - 1))
        self.assertEqual(c.residue_p, (3 * 5) % 7)

    def test_product(self):
        P = self.proc.product([3, 5, 7])
        self.assertTrue(P.verify())
        expected_2 = (3 * 5 * 7) & ((1 << 6) - 1)
        expected_p = (3 * 5 * 7) % 7
        self.assertEqual(P.component_2.value, expected_2)
        self.assertEqual(P.residue_p, expected_p)

    def test_cycle_product(self):
        nums = [CRTDualNumber(i, k=6, p=7, g_p=self.proc.g_p) for i in [3, 5, 7]]
        P = self.proc.cycle_product(nums)
        self.assertTrue(P.verify())

    def test_convergence_ratio(self):
        n = CRTDualNumber(5, k=6, p=7, g_p=self.proc.g_p)
        r = self.proc.convergence_ratio_2adic(n)
        self.assertGreaterEqual(r, 0.0)
        self.assertLessEqual(r, 1.0)


class TestCombinedStability(unittest.TestCase):
    def test_returns_dict(self):
        result = combined_stability(k=6, p=7, num_cycles=10)
        self.assertIn("pearson_r", result)
        self.assertIn("n_samples", result)
        self.assertIn("mean_v2", result)


class TestNonAbelianMatrixOps(unittest.TestCase):
    def test_mat_mul(self):
        A = [[1, 2], [3, 4]]
        B = [[5, 6], [7, 8]]
        C = _mat_mul(A, B, 100)
        self.assertEqual(C[0][0], (1*5 + 2*7) % 100)
        self.assertEqual(C[0][1], (1*6 + 2*8) % 100)

    def test_mat_det(self):
        M = [[1, 2], [3, 4]]
        det = _mat_det(M, 100)
        self.assertEqual(det, (1*4 - 2*3) % 100)


class TestNonAbelianCRTDual(unittest.TestCase):
    def setUp(self):
        self.nc = NonAbelianCRTDual(k=6, p=7)

    def test_holonomy_identity(self):
        mats = [[[1, 0], [0, 1]] for _ in range(3)]
        H = self.nc.holonomy(mats)
        self.assertEqual(H, [[1, 0], [0, 1]])

    def test_invariants_keys(self):
        mats = [[[2, 1], [1, 2]]]
        inv = self.nc.invariants(mats)
        for key in ("det_mod2k", "alpha_det", "trace_modp"):
            self.assertIn(key, inv)

    def test_convergence_ratio(self):
        mats = [[[2, 1], [1, 2]]]
        r = self.nc.convergence_ratio_full(mats)
        self.assertGreaterEqual(r, 0.0)


class TestPhaseAlignmentExperiment(unittest.TestCase):
    def test_phase_alignment_experiment(self):
        result = phase_alignment_experiment(k=5, p=7, n_cycles=10)
        self.assertIn("alignment", result)
        self.assertIn("n_trials", result)


class TestKroneckerCliffScorer(unittest.TestCase):
    def test_score_factors(self):
        factors = [np.random.randn(4, 4).astype(np.float32) for _ in range(2)]
        scorer = KroneckerCliffScorer(factors, k_range=range(4, 8))
        results = scorer.score_factors()
        self.assertEqual(len(results), 2)
        for key in ("factor_0", "factor_1"):
            self.assertIn(key, results)
            self.assertIn("summary", results[key])

    def test_composition_cliff(self):
        factors = [np.random.randn(4, 4).astype(np.float32) for _ in range(2)]
        scorer = KroneckerCliffScorer(factors, k_range=range(4, 8))
        cc = scorer.composition_cliff()
        if cc is not None:
            self.assertGreater(cc, 0)

    def test_fragile_factors(self):
        factors = [np.random.randn(4, 4).astype(np.float32) for _ in range(2)]
        scorer = KroneckerCliffScorer(factors, k_range=range(4, 8))
        fragile = scorer.fragile_factors(threshold=10)
        self.assertIsInstance(fragile, list)

    def test_print_report(self):
        factors = [np.random.randn(4, 4).astype(np.float32) for _ in range(2)]
        scorer = KroneckerCliffScorer(factors, k_range=range(4, 8))
        report = scorer.print_report()
        self.assertIsInstance(report, str)


class TestSemiringCliffScore(unittest.TestCase):
    def test_standard_is_min(self):
        self.assertEqual(semiring_cliff_score([3.0, 7.0, 5.0], "standard"), 3.0)

    def test_tropical_is_max(self):
        self.assertEqual(semiring_cliff_score([3.0, 7.0, 5.0], "tropical"), 7.0)

    def test_all_none_returns_none(self):
        self.assertIsNone(semiring_cliff_score([None, None]))

    def test_partial_none(self):
        self.assertEqual(semiring_cliff_score([None, 5.0], "standard"), 5.0)


if __name__ == "__main__":
    unittest.main()
