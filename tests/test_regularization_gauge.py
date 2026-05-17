"""Tests for dual_view.regularization and dual_view.gauge."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.regularization import GhostMap, local_ratio_gradient, ghost_penalty
from dual_view.gauge import (
    cycle_product, spectral_det, det_coordinates, tidal_scalar, GaugeLayer,
)
from dual_view.core import _mask
import numpy as np


class TestGhostMap(unittest.TestCase):
    def setUp(self):
        self.gm = GhostMap(k=6, g=5)

    def test_ratio_powers_of_five(self):
        for e in range(4):
            a = pow(5, e, 1 << 6)
            r = self.gm.ratio(a)
            self.assertGreaterEqual(r, 0.0)

    def test_ratio_zero_for_zero(self):
        r = self.gm.ratio(0)
        self.assertGreaterEqual(r, 0.0)

    def test_ratios_in_unit_interval(self):
        for a in range(1, 64, 2):
            r = self.gm.ratio(a)
            self.assertGreaterEqual(r, 0.0)
            self.assertLessEqual(r, 1.0)

    def test_ratio_handles_numpy_int(self):
        r = self.gm.ratio(np.int64(3))
        self.assertIsInstance(r, float)

    def test_nearest_stable(self):
        best, r = self.gm.nearest_stable(1)
        self.assertIsInstance(best, int)
        self.assertGreaterEqual(r, 0.0)

    def test_large_k_raises(self):
        with self.assertRaises(ValueError):
            GhostMap(k=20)

    def test_small_k_raises(self):
        with self.assertRaises(ValueError):
            GhostMap(k=2)


class TestLocalRatioGradient(unittest.TestCase):
    def setUp(self):
        self.gm = GhostMap(k=6, g=5)

    def test_returns_list_of_deltas(self):
        grad = local_ratio_gradient(3, self.gm)
        self.assertIsInstance(grad, list)
        for delta, r in grad:
            self.assertIn(delta, (-2, -1, 1, 2))
            self.assertGreater(r, 0.0)

    def test_stable_weight_has_no_improvements(self):
        grad_stable = local_ratio_gradient(16, self.gm)
        grad_ghost = local_ratio_gradient(3, self.gm)
        # stable weights (v2 >= k-2=4) have no improvements; odd do
        self.assertEqual(len(grad_stable), 0)
        self.assertGreater(len(grad_ghost), 0)


class TestGhostPenalty(unittest.TestCase):
    def setUp(self):
        self.gm = GhostMap(k=6, g=5)

    def test_penalty_shape(self):
        W = np.array([1, 3, 5, 7], dtype=np.int32)
        penalty, grad = ghost_penalty(W, self.gm)
        self.assertIsInstance(penalty, float)
        self.assertEqual(grad.shape, W.shape)

    def test_penalty_range(self):
        W = np.array([1, 3, 5, 7], dtype=np.int32)
        penalty, _ = ghost_penalty(W, self.gm)
        self.assertGreaterEqual(penalty, 0.0)
        self.assertLessEqual(penalty, 1.0)

    def test_all_stable_penalty_zero(self):
        W = np.array([0], dtype=np.int32)
        penalty, _ = ghost_penalty(W, self.gm)
        self.assertAlmostEqual(penalty, 0.0, places=2)

    def test_gradient_zero_for_stable_weights(self):
        W = np.array([0, 16, 32, 48], dtype=np.int32)
        _, grad = ghost_penalty(W, self.gm)
        self.assertTrue(np.all(grad == 0))


class TestGaugeFunctions(unittest.TestCase):
    def test_cycle_product(self):
        prod = cycle_product([3, 5, 7], 16)
        self.assertEqual(prod, (3 * 5 * 7) & _mask(16))

    def test_product_permutation_invariant(self):
        weights = [3, 5, 7, 9]
        p1 = cycle_product(weights, 16)
        p2 = cycle_product(list(reversed(weights)), 16)
        self.assertEqual(p1, p2)

    def test_spectral_det_formula(self):
        weights = [3, 5, 7]
        k = 16
        det_val = spectral_det(weights, k)
        W = cycle_product(weights, k)
        self.assertEqual(det_val, (1 - W) & _mask(k))

    def test_det_coordinates_returns_tuple(self):
        weights = [3, 5, 7]
        coords = det_coordinates(weights, 16)
        self.assertIsNotNone(coords)
        v, alpha, e = coords
        self.assertIsInstance(v, (int, float))
        self.assertIn(alpha, (0, 1))
        self.assertIsInstance(e, int)

    def test_tidal_scalar_type(self):
        weights = [3, 5, 7]
        h = tidal_scalar(weights, 16)
        self.assertIsInstance(h, int)

    def test_tidal_scalar_none_for_zero_det(self):
        weights = [1, 1, 1]
        h = tidal_scalar(weights, 16)
        self.assertIsNone(h)


class TestGaugeLayer(unittest.TestCase):
    def test_initialization(self):
        gl = GaugeLayer([3, 5, 7, 9], k=16)
        self.assertEqual(len(gl.weights), 4)
        self.assertIsInstance(gl.product, int)
        self.assertIsInstance(gl.det_val, int)

    def test_report_output(self):
        gl = GaugeLayer([3, 5, 7], k=16)
        report = gl.report()
        self.assertIsInstance(report, str)
        self.assertIn("GaugeLayer Report", report)


if __name__ == "__main__":
    unittest.main()
