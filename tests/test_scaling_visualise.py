"""Tests for dual_view.scaling and dual_view.visualise."""
import unittest
import io

from dual_view.scaling import scale_weights, auto_scale, common_scales
from dual_view.visualise import (
    cliff_matrix, sector_matrix, valuation_matrix,
    print_cliff_ascii, cliff_stats_by_layer, show_dual_bits,
)
from dual_view.thermodynamics import SeedThermodynamics
from dual_view.core import _valuation
import numpy as np


class TestScaleWeights(unittest.TestCase):
    def test_shape_preserved(self):
        W = np.random.randn(4, 4).astype(np.float32)
        W_int, meta = scale_weights(W, 128.0)
        self.assertEqual(W_int.shape, W.shape)

    def test_dtype_int32(self):
        W = np.random.randn(10).astype(np.float32)
        W_int, _ = scale_weights(W, 128.0)
        self.assertEqual(W_int.dtype, np.int32)

    def test_meta_keys(self):
        W = np.random.randn(10).astype(np.float32)
        _, meta = scale_weights(W, 128.0)
        for key in ('scale', 'mode', 'ensure_odd', 'n_even', 'n_zero', 'range'):
            self.assertIn(key, meta)

    def test_meta_n_even_n_zero(self):
        W = np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32)
        _, meta = scale_weights(W, 1.0)
        self.assertIn('n_even', meta)
        self.assertIn('n_zero', meta)
        self.assertIn('range', meta)
        self.assertIsInstance(meta['range'], tuple)
        self.assertEqual(len(meta['range']), 2)

    def test_meta_v2_hist_int_keys(self):
        W = np.array([1.0, 2.0, 4.0], dtype=np.float32)
        _, meta = scale_weights(W, 1.0)
        self.assertIn('v2_hist', meta)
        # keys should be ints
        for k in meta['v2_hist']:
            self.assertIsInstance(k, int)

    def test_scale_applied(self):
        W = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        W_int, _ = scale_weights(W, 2.0)
        self.assertTrue(np.allclose(W_int, [2, 4, 6], atol=1))

    def test_ensure_odd_produces_all_odd(self):
        W = np.array([2.0, 4.0, 6.0], dtype=np.float32)
        W_int, _ = scale_weights(W, 1.0, ensure_odd=True)
        self.assertTrue(np.all(W_int % 2 == 1))

    def test_floor_mode_handles_negatives(self):
        W = np.array([-1.5, 0.0, 1.5], dtype=np.float32)
        W_int, _ = scale_weights(W, 2.0, mode="floor")
        self.assertEqual(W_int.dtype, np.int32)

    def test_invalid_mode_raises(self):
        W = np.random.randn(5).astype(np.float32)
        with self.assertRaises(ValueError):
            scale_weights(W, 128.0, mode="invalid")


class TestAutoScale(unittest.TestCase):
    def test_scale_reasonable_range(self):
        W = np.random.randn(100).astype(np.float32)
        s = auto_scale(W, 8)
        self.assertGreater(s, 0)
        self.assertLess(s, 500)

    def test_zero_input_returns_one(self):
        s = auto_scale(np.zeros((10,), dtype=np.float32), 8)
        self.assertEqual(s, 1.0)


class TestCommonScales(unittest.TestCase):
    def test_contains_standard_depths(self):
        scales = common_scales()
        for name in ('INT7', 'INT8', 'INT9', 'INT16'):
            self.assertIn(name, scales)


class TestVisualise(unittest.TestCase):
    def setUp(self):
        self.W = np.array([[1, 2], [3, 4]], dtype=np.int64)
        self.st = SeedThermodynamics(k=8)
        self.st(self.W, range(4, 10))
        self.st.compute()

    def test_cliff_matrix_shape(self):
        C = cliff_matrix(self.st, self.W.shape)
        self.assertEqual(C.shape, self.W.shape)

    def test_cliff_matrix_nan_for_even(self):
        C = cliff_matrix(self.st, self.W.shape)
        # w=2 and w=4 are even -> NaN
        self.assertTrue(np.isnan(C[0, 1]))
        self.assertTrue(np.isnan(C[1, 1]))

    def test_sector_matrix_shape(self):
        S = sector_matrix(self.W, 8)
        self.assertEqual(S.shape, self.W.shape)

    def test_sector_matrix_values(self):
        S = sector_matrix(self.W, 8)
        # Odd weights have sector 0 or 1
        if not np.isnan(S[0, 0]):
            self.assertIn(S[0, 0], (0.0, 1.0))

    def test_valuation_matrix(self):
        V = valuation_matrix(self.W)
        self.assertEqual(V.shape, self.W.shape)
        # v2(1) = 0, v2(2) = 1, v2(3) = 0, v2(4) = 2
        self.assertEqual(V[0, 0], 0.0)  # 1
        self.assertEqual(V[0, 1], 1.0)  # 2
        self.assertEqual(V[1, 0], 0.0)  # 3
        self.assertEqual(V[1, 1], 2.0)  # 4

    def test_print_cliff_ascii(self):
        C = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        f = io.StringIO()
        # Capture stdout by redirecting
        import contextlib
        with contextlib.redirect_stdout(f):
            print_cliff_ascii(C, "Test")
        output = f.getvalue()
        self.assertIn("Test", output)

    def test_cliff_stats_by_layer(self):
        C = np.array([[1.0, np.nan], [3.0, 4.0]], dtype=np.float64)
        layers = {"layer1": C}
        stats = cliff_stats_by_layer(layers)
        self.assertIn("layer1", stats)
        self.assertIn("Mean Cliff", stats)


class TestShowDualBits(unittest.TestCase):
    """Tests for show_dual_bits — 2-adic bit annotation."""

    def test_contains_correct_n(self):
        s = show_dual_bits(42, k=16)
        self.assertIn("n = 42", s)

    def test_contains_annotation_key(self):
        s = show_dual_bits(7, k=8)
        self.assertIn("annotation", s)

    def test_contains_valuation(self):
        s = show_dual_bits(8, k=8)
        self.assertIn("v = 3", s)

    def test_zero_output(self):
        s = show_dual_bits(0, k=8)
        self.assertIn("n = 0", s)

    def test_odd_includes_decomposition(self):
        s = show_dual_bits(7, k=8)
        self.assertIn("Dual decomposition", s)

    def test_even_no_decomposition(self):
        s = show_dual_bits(8, k=8)
        self.assertNotIn("Dual decomposition", s)

    def test_label_appears(self):
        s = show_dual_bits(5, k=8, label="TEST")
        self.assertIn("TEST", s)

    def test_reconstruction_pass(self):
        s = show_dual_bits(5, k=8)
        self.assertIn("PASS", s)

    def test_k_too_small_raises(self):
        with self.assertRaises(ValueError):
            show_dual_bits(1, k=2)


if __name__ == "__main__":
    unittest.main()
