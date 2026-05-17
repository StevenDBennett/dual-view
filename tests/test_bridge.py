"""Tests for dual_view.bridge."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.bridge import (
    ButterflyBridge, LayerReport, ModelReport, SpectralThermodynamics,
    quantize, _depths, _geometric_null, seed_S1, seed_S2, seed_S3, depth_char,
)
from dual_view.core import _valuation, _mask
import numpy as np


class TestQuantize(unittest.TestCase):
    def test_quantize_identity(self):
        W = np.array([[-1.0, 0.0, 1.0]])
        Q = quantize(W, k=3)
        self.assertEqual(Q.shape, W.shape)
        self.assertEqual(Q.dtype, np.int64)
        self.assertTrue(np.all(Q >= 0))
        self.assertTrue(np.all(Q < 8))

    def test_quantize_uniform(self):
        W = np.ones((3, 3))
        Q = quantize(W, k=4)
        self.assertTrue(np.all(Q == 0))

    def test_quantize_range(self):
        W = np.linspace(-1, 1, 5).reshape(1, 5)
        Q = quantize(W, k=3)
        self.assertTrue(np.all(Q >= 0))
        self.assertTrue(np.all(Q < 8))


class TestDepths(unittest.TestCase):
    def test_depths_known_values(self):
        W = np.array([[0, 1, 2, 4, 8]], dtype=np.int64)
        d = _depths(W, k=6)
        # 0 -> capped to k-1 = 5
        self.assertEqual(d[0], 5)
        # 1 -> v2(1)=0
        self.assertEqual(d[1], 0)
        # 2 -> v2(2)=1
        self.assertEqual(d[2], 1)
        # 4 -> v2(4)=2
        self.assertEqual(d[3], 2)
        # 8 -> v2(8)=3
        self.assertEqual(d[4], 3)

    def test_depths_capped(self):
        W = np.array([[1 << 10]], dtype=np.int64)
        d = _depths(W, k=4)
        self.assertEqual(d[0], 3)


class TestGeometricNull(unittest.TestCase):
    def test_sums_to_one(self):
        for k in range(2, 10):
            h = _geometric_null(k)
            self.assertAlmostEqual(h.sum(), 1.0, places=10)

    def test_decreasing(self):
        h = _geometric_null(6)
        for i in range(len(h) - 2):
            self.assertGreater(h[i], h[i + 1])

    def test_length(self):
        self.assertEqual(len(_geometric_null(8)), 8)


class TestSpectralThermodynamics(unittest.TestCase):
    def test_identity_matrix(self):
        S = np.eye(4)
        t = SpectralThermodynamics.analyze(S)
        self.assertTrue(t.is_unitary)
        self.assertTrue(t.is_conservative)
        self.assertAlmostEqual(t.spectral_radius, 1.0)

    def test_nilpotent_matrix(self):
        S = np.array([[0, 1], [0, 0]], dtype=float)
        t = SpectralThermodynamics.analyze(S)
        self.assertTrue(t.is_nilpotent)
        self.assertEqual(t.character(), "NILPOTENT")

    def test_expansive_matrix(self):
        S = np.array([[2, 0], [0, 0.5]], dtype=float)
        t = SpectralThermodynamics.analyze(S)
        self.assertTrue(t.is_expansive)
        self.assertAlmostEqual(t.spectral_radius, 2.0)

    def test_contractive_matrix(self):
        S = np.array([[0.5, 0], [0, 0.25]], dtype=float)
        t = SpectralThermodynamics.analyze(S)
        self.assertTrue(t.is_contractive)
        self.assertAlmostEqual(t.spectral_radius, 0.5)

    def test_zero_matrix(self):
        S = np.zeros((3, 3))
        t = SpectralThermodynamics.analyze(S)
        self.assertAlmostEqual(t.spectral_radius, 0.0)
        self.assertEqual(t.lyapunov_exponent, -float("inf"))

    def test_string_representation(self):
        S = np.eye(2)
        t = SpectralThermodynamics.analyze(S)
        s = str(t)
        self.assertIn("CONSERVATIVE", s)
        self.assertIn("=", s)


class TestSeedS2(unittest.TestCase):
    def test_square_map(self):
        W = np.random.randn(8, 8)
        S2, t = seed_S2(W)
        self.assertEqual(S2.shape, (8, 8))

    def test_non_square_map(self):
        W = np.random.randn(8, 16)
        S2, t = seed_S2(W)
        self.assertEqual(S2.shape, (8, 8))

    def test_thermo_returned(self):
        W = np.random.randn(4, 4)
        _, t = seed_S2(W)
        self.assertIsInstance(t, SpectralThermodynamics)


class TestSeedS3(unittest.TestCase):
    def test_shape(self):
        W = np.array([[1, 3, 5, 7]], dtype=np.int64)
        S3, t = seed_S3(W, k=4)
        self.assertEqual(S3.shape, (2, 2))

    def test_symmetric(self):
        W = np.array([[1, 2, 3]], dtype=np.int64)
        S3, _ = seed_S3(W, k=4)
        self.assertEqual(S3[0, 1], S3[1, 0])

    def test_all_odd_values(self):
        W = np.array([[1, 1, 1, 1]], dtype=np.int64)
        S3, t = seed_S3(W, k=4)
        # all have alpha=0 (bit 1 is 0)
        self.assertAlmostEqual(S3[0, 0], 1.0)
        self.assertAlmostEqual(S3[0, 1], 0.0)

    def test_all_alpha_one(self):
        # numbers with bit 1 set: 3, 7, 11, 15
        W = np.array([[3, 7, 11, 15]], dtype=np.int64)
        S3, _ = seed_S3(W, k=4)
        self.assertAlmostEqual(S3[0, 1], 1.0)
        self.assertAlmostEqual(S3[1, 0], 1.0)


class TestSeedS1(unittest.TestCase):
    def test_returns_expected_shapes(self):
        W = np.array([[1, 3, 5, 7]], dtype=np.int64)
        h, dev, S1, t = seed_S1(W, k=4)
        self.assertEqual(len(h), 4)
        self.assertEqual(len(dev), 4)
        self.assertEqual(S1.shape, (4, 4))
        self.assertIsInstance(t, SpectralThermodynamics)

    def test_histogram_sums_to_one(self):
        W = np.random.randint(0, 256, size=(20, 20)).astype(np.int64)
        h, _, _, _ = seed_S1(W, k=8)
        self.assertAlmostEqual(h.sum(), 1.0, places=10)


class TestDepthChar(unittest.TestCase):
    def test_returns_string(self):
        W = np.random.randint(0, 256, size=(20, 20)).astype(np.int64)
        h, dev, _, _ = seed_S1(W, k=8)
        c = depth_char(h, dev, 8)
        self.assertIsInstance(c, str)
        self.assertTrue(c in ("NEUTRAL", "EXPANSIVE", "CONTRACTIVE") or
                        c.startswith("STRUCTURED/"))


class TestButterflyBridge(unittest.TestCase):
    def setUp(self):
        self.bridge = ButterflyBridge(k=8)
        self.rng = np.random.default_rng(42)

    def test_analyse_layer_returns_layer_report(self):
        W = self.rng.normal(0, 0.02, (16, 16))
        report = self.bridge.analyse_layer(W, name="test")
        self.assertIsInstance(report, LayerReport)
        self.assertEqual(report.name, "test")
        self.assertEqual(report.shape, (16, 16))
        self.assertEqual(report.k, 8)

    def test_layer_report_has_all_fields(self):
        W = self.rng.normal(0, 0.02, (4, 4))
        r = self.bridge.analyse_layer(W, name="layer1")
        self.assertIsInstance(r.mean_v, float)
        self.assertIsInstance(r.depth_entropy, float)
        self.assertIsInstance(r.zero_frac, float)
        self.assertIsInstance(r.alpha_frac, float)
        self.assertIsInstance(r.thermo_S1, SpectralThermodynamics)
        self.assertIsInstance(r.thermo_S2, SpectralThermodynamics)
        self.assertIsInstance(r.thermo_S3, SpectralThermodynamics)

    def test_analyse_model_returns_model_report(self):
        layers = {
            "embed": self.rng.normal(0, 0.1, (8, 8)),
            "proj": self.rng.normal(0, 0.01, (8, 16)),
        }
        report = self.bridge.analyse_model(layers)
        self.assertIsInstance(report, ModelReport)
        self.assertEqual(len(report.layers), 2)

    def test_model_report_trajectory_length(self):
        layers = {
            "a": self.rng.normal(0, 0.1, (4, 4)),
            "b": self.rng.normal(0, 0.01, (4, 4)),
            "c": self.rng.normal(0, 0.001, (4, 4)),
        }
        report = self.bridge.analyse_model(layers)
        lya1 = report._traj(lambda r: r.thermo_S1.lyapunov_exponent)
        self.assertEqual(len(lya1), 3)

    def test_model_report_report_string(self):
        layers = {
            "test": self.rng.normal(0, 0.02, (4, 4)),
        }
        report = self.bridge.analyse_model(layers)
        s = report.report()
        self.assertIn("test", s)

    def test_layer_report_report_string(self):
        W = self.rng.normal(0, 0.02, (4, 4))
        r = self.bridge.analyse_layer(W, name="layer_x")
        s = r.report()
        self.assertIn("layer_x", s)
        self.assertIn("CONSENSUS", s)

    def test_boundaries_empty_with_single_layer(self):
        layers = {"only": self.rng.normal(0, 0.02, (4, 4))}
        report = self.bridge.analyse_model(layers)
        self.assertEqual(report.boundaries(), [])

    def test_consensus_known_weights(self):
        W = self.rng.normal(0, 0.02, (4, 4))
        r = self.bridge.analyse_layer(W)
        c = r.consensus()
        self.assertIn(c, ("NEUTRAL  (matches arithmetic null)",
                          "CONSENSUS: EXPANSIVE",
                          "CONSENSUS: CONTRACTIVE",
                          "CONSENSUS: MIXED",
                          "SPLIT   depth=EXPANSIVE  map=EXPANSIVE",
                          "SPLIT   depth=CONTRACTIVE  map=EXPANSIVE"))
        # consensus is always a string
        self.assertIsInstance(c, str)

    def test_different_k_values(self):
        for k in (4, 6, 8):
            bridge = ButterflyBridge(k=k)
            W = self.rng.normal(0, 0.02, (4, 4))
            r = bridge.analyse_layer(W)
            self.assertEqual(r.k, k)


class TestEdgeCases(unittest.TestCase):
    def test_all_zero_weights(self):
        W = np.zeros((4, 4))
        bridge = ButterflyBridge(k=8)
        r = bridge.analyse_layer(W)
        self.assertEqual(r.n, 16)
        self.assertAlmostEqual(r.zero_frac, 1.0)

    def test_single_element(self):
        W = np.array([[42.0]])
        bridge = ButterflyBridge(k=8)
        r = bridge.analyse_layer(W)
        self.assertEqual(r.shape, (1, 1))

    def test_negative_weights(self):
        W = np.array([[-1.0, -2.0, -3.0]])
        bridge = ButterflyBridge(k=8)
        r = bridge.analyse_layer(W)
        self.assertEqual(r.shape, (1, 3))

    def test_large_k(self):
        W = np.random.randn(10, 10)
        bridge = ButterflyBridge(k=16)
        r = bridge.analyse_layer(W)
        self.assertEqual(r.k, 16)

    def test_small_k(self):
        W = np.array([[1.0, 2.0], [3.0, 4.0]])
        bridge = ButterflyBridge(k=3)
        r = bridge.analyse_layer(W)
        self.assertEqual(r.k, 3)

    def test_empty_model(self):
        report = ModelReport()
        self.assertEqual(len(report.layers), 0)
        self.assertEqual(report.boundaries(), [])

    def test_model_report_empty(self):
        report = ModelReport()
        s = report.report()
        self.assertIn("UNIFIED BUTTERFLY", s)


class TestSpectralThermodynamicsEdgeCases(unittest.TestCase):
    def test_character_nilpotent_string(self):
        S = np.array([[0, 0], [1, 0]], dtype=float)
        t = SpectralThermodynamics.analyze(S)
        self.assertEqual(t.character(), "NILPOTENT")

    def test_character_unitary_string(self):
        S = np.eye(3)
        t = SpectralThermodynamics.analyze(S)
        self.assertEqual(t.character(), "CONSERVATIVE/UNITARY")

    def test_character_mixed(self):
        S = np.array([[2, 0], [0, 0.5]], dtype=float)
        t = SpectralThermodynamics.analyze(S)
        # expansive takes precedence
        self.assertEqual(t.character(), "EXPANSIVE")

    def test_lyapunov_for_zero_radius(self):
        S = np.zeros((2, 2))
        t = SpectralThermodynamics.analyze(S)
        self.assertEqual(t.lyapunov_exponent, -float("inf"))

    def test_min_eigenvalue(self):
        S = np.diag([3.0, 1.0, 0.5])
        t = SpectralThermodynamics.analyze(S)
        self.assertAlmostEqual(t.min_eigenvalue_magnitude, 0.5)


class TestDepthCharEdgeCases(unittest.TestCase):
    def test_flat_weights(self):
        W = np.ones((10, 10), dtype=np.int64)
        h, dev, _, _ = seed_S1(W, k=8)
        c = depth_char(h, dev, 8)
        self.assertIsInstance(c, str)


class TestSymmetricMinMax(unittest.TestCase):
    def test_quantize_symmetry(self):
        W = np.array([[-1.0, 0.0, 1.0]])
        Q = quantize(W, k=4)
        self.assertTrue(np.all(Q >= 0))
        self.assertTrue(np.all(Q < 16))

    def test_quantize_constant_input(self):
        W = np.ones((5,)) * 42.0
        Q = quantize(W, k=8)
        self.assertTrue(np.all(Q == 0))


if __name__ == "__main__":
    unittest.main()
