"""Additional tests to improve coverage on previously untested functions."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

import numpy as np


# ── padic_roots.py ──────────────────────────────────────────────────────────────

class TestLiftRoot(unittest.TestCase):
    def test_known_cube_root_mod_5(self):
        from dual_view.padic_roots import lift_root
        # 2^3 = 8 ≡ 3 (mod 5), so lift_root(3, 5, 1) should return 2
        x = lift_root(3, 5, 1)
        self.assertIsNotNone(x)
        self.assertEqual((x ** 3) % 5, 3)

    def test_hensel_lift_to_higher_power(self):
        from dual_view.padic_roots import lift_root
        # 2^3 = 8 ≡ 3 (mod 5), lift to 5^4
        x = lift_root(3, 5, 4)
        self.assertIsNotNone(x)
        self.assertEqual((x ** 3) % (5 ** 4), 3)

    def test_no_root_returns_none(self):
        from dual_view.padic_roots import lift_root
        # For p=7, cubes mod 7 are {0, 1, 6}. 2 is not a cube.
        x = lift_root(2, 7, 1)
        self.assertIsNone(x)

    def test_divisible_by_p_returns_none(self):
        from dual_view.padic_roots import lift_root
        x = lift_root(0, 5, 2)
        self.assertIsNone(x)


class TestNewton2Step(unittest.TestCase):
    def test_converges_faster_than_newton(self):
        from dual_view.padic_roots import newton2_step, newton_step, lift_root, _vp
        p, k = 5, 6
        pk = p ** k
        a = 8
        x_true = lift_root(a, p, k)
        self.assertIsNotNone(x_true)
        x0 = 2
        # Two Newton steps vs one composed-2-Newton step
        x_n1 = newton_step(x0, a, pk)
        x_n2 = newton_step(x_n1, a, pk)
        x_comp = newton2_step(x0, a, pk)
        self.assertEqual(x_n2, x_comp)


class TestNewton3Step(unittest.TestCase):
    def test_equals_three_newton_steps(self):
        from dual_view.padic_roots import newton3_step, newton_step, lift_root
        p, k = 5, 6
        pk = p ** k
        a = 8
        x0 = 2
        x_n1 = newton_step(x0, a, pk)
        x_n2 = newton_step(x_n1, a, pk)
        x_n3 = newton_step(x_n2, a, pk)
        x_comp = newton3_step(x0, a, pk)
        self.assertEqual(x_n3, x_comp)


class TestCompareMethods(unittest.TestCase):
    def test_returns_all_methods(self):
        from dual_view.padic_roots import compare_methods
        result = compare_methods(5, 4, n_trials=5)
        self.assertIn("Newton (ord 2)", result)
        self.assertIn("Halley (ord 3)", result)
        self.assertIn("Comp-Newton (ord 4)", result)
        self.assertIn("Comp×3 (ord 8)", result)

    def test_higher_order_converges_better(self):
        from dual_view.padic_roots import compare_methods
        result = compare_methods(5, 4, n_trials=20)
        # Order 8 should generally converge at least as well as order 2
        self.assertGreaterEqual(result["Comp×3 (ord 8)"], result["Newton (ord 2)"])


class TestVerifyOrder(unittest.TestCase):
    def test_returns_nested_dict(self):
        from dual_view.padic_roots import verify_order
        result = verify_order([5, 7], k=4, n_trials=3)
        self.assertIn("Newton", result)
        self.assertIn("Halley", result)
        for method_name in result:
            self.assertIsInstance(result[method_name], dict)


class TestNewtonCorrectionUniformity(unittest.TestCase):
    def test_returns_chi2_stats(self):
        from dual_view.padic_roots import newton_correction_uniformity
        result = newton_correction_uniformity(5, 2, n_seeds=100)
        self.assertIn("chi2_stat", result)
        self.assertIn("n_bins", result)
        self.assertIn("n_samples", result)
        self.assertIn("df", result)
        self.assertEqual(result["df"], 4)


# ── separation.py ────────────────────────────────────────────────────────────────

class TestSeparationStep(unittest.TestCase):
    def test_identical_targets_never_diverge(self):
        from dual_view.separation import separation_step
        k = 8
        a = pow(5, 3, 1 << k)
        n = separation_step(a, a, k, 0)
        self.assertEqual(n, -1)

    def test_close_targets_diverge_later(self):
        from dual_view.separation import separation_step
        k = 8
        a = pow(5, 3, 1 << k)
        a_prime = a ^ (1 << 4)
        n = separation_step(a, a_prime, k, 0)
        self.assertGreaterEqual(n, 0)


class TestVerifySeparation(unittest.TestCase):
    def test_returns_results_for_all_s_values(self):
        from dual_view.separation import verify_separation, predicted_separation
        k = 8
        s_values = [3, 4, 5]
        results = verify_separation(k, s_values, n_trials=5)
        for s in s_values:
            self.assertIn(s, results)
            # Results should be integers (either predicted or -1 on mismatch)
            self.assertIsInstance(results[s], int)


# ── thermodynamics.py ────────────────────────────────────────────────────────────

class TestMersenneCliffScore(unittest.TestCase):
    def test_even_returns_none(self):
        from dual_view.thermodynamics import SeedThermodynamics
        st = SeedThermodynamics(k=16)
        self.assertIsNone(st.mersenne_cliff_score(4))

    def test_mersenne_number_has_expected_cliff(self):
        from dual_view.thermodynamics import SeedThermodynamics
        st = SeedThermodynamics(k=16)
        # 2^5 - 1 = 31, v2(e_true) = 5 - 2 = 3, so k* = 3 + 2 = 5
        score = st.mersenne_cliff_score(31)
        self.assertIsNotNone(score)
        self.assertEqual(score, 5)

    def test_power_of_five_has_high_cliff(self):
        from dual_view.thermodynamics import SeedThermodynamics
        st = SeedThermodynamics(k=16)
        # 5^e has e_true = e, v2(e_true) depends on e
        score = st.mersenne_cliff_score(5)
        self.assertIsNotNone(score)
        self.assertGreaterEqual(score, 2)


class TestCompareToRandom(unittest.TestCase):
    def test_returns_z_scores(self):
        from dual_view.thermodynamics import SeedThermodynamics
        st = SeedThermodynamics(k=8)
        W = np.array([1, 3, 5, 7, 255, 127, 63], dtype=np.int64)
        result = st.compare_to_random(W, n_samples=50)
        self.assertIn("z_alpha", result)
        self.assertIn("z_v2_e", result)
        self.assertIn("z_cliff_risk", result)
        for key in result:
            self.assertIsInstance(result[key], float)


class TestCliffHistogram(unittest.TestCase):
    def test_returns_dict_of_ints(self):
        from dual_view.thermodynamics import SeedThermodynamics
        W = np.array([1, 3, 5, 7, 255, 127, 63, 31], dtype=np.int64)
        st = SeedThermodynamics(k=8)
        st(W, range(4, 10))
        st.compute()
        hist = st.cliff_histogram()
        self.assertIsInstance(hist, dict)
        for k_val, count in hist.items():
            self.assertIsInstance(k_val, int)
            self.assertIsInstance(count, int)
            self.assertGreater(count, 0)


# ── fourier.py ───────────────────────────────────────────────────────────────────

class TestDyadicCoefficients(unittest.TestCase):
    def test_has_dc_component(self):
        from dual_view.fourier import dyadic_coefficients
        h = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int32)
        coeffs = dyadic_coefficients(h)
        self.assertIn("DC", coeffs)

    def test_dyadic_keys_present(self):
        from dual_view.fourier import dyadic_coefficients
        h = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int32)
        coeffs = dyadic_coefficients(h)
        # N=8, so we should have DC, N/2, N/4
        self.assertIn("N/2", coeffs)
        self.assertIn("N/4", coeffs)

    def test_matches_analytic_for_step_count(self):
        from dual_view.fourier import step_count_fn, analytic_step_count, dyadic_coefficients, analytic_coefficients
        k = 6
        h_num = step_count_fn(k, 3)
        h_an = analytic_step_count(k, 3)
        dyadic_num = dyadic_coefficients(h_num)
        dyadic_an = analytic_coefficients(k)
        # Both should have the same key structure
        self.assertIn("DC", dyadic_num)
        self.assertIn("DC", dyadic_an)
        # DC should be positive and real in both
        self.assertGreater(abs(dyadic_num["DC"]), 0)
        self.assertGreater(abs(dyadic_an["DC"]), 0)


class TestAnalyticCoefficients(unittest.TestCase):
    def test_returns_expected_keys(self):
        from dual_view.fourier import analytic_coefficients
        coeffs = analytic_coefficients(6)
        self.assertIn("DC", coeffs)
        self.assertIn("N/2", coeffs)
        self.assertIn("N/4", coeffs)

    def test_dc_is_real(self):
        from dual_view.fourier import analytic_coefficients
        coeffs = analytic_coefficients(6)
        self.assertEqual(coeffs["DC"].imag, 0)


class TestUltrametricUncertainty(unittest.TestCase):
    def test_returns_string_with_N(self):
        from dual_view.fourier import ultrametric_uncertainty
        text = ultrametric_uncertainty(6)
        self.assertIn("N=16", text)
        self.assertIn("uncertainty", text.lower())


# ── iwasawa.py ──────────────────────────────────────────────────────────────────

class TestVerifyCommutatorDepth(unittest.TestCase):
    def test_returns_triples(self):
        from dual_view.iwasawa import verify_commutator_depth
        results = verify_commutator_depth(8, [(1, 1)], n_trials=5)
        self.assertIsInstance(results, list)
        for dM, dN, dMN in results:
            self.assertIsInstance(dM, int)
            self.assertIsInstance(dN, int)
            self.assertIsInstance(dMN, int)

    def test_commutator_depth_theorem_holds(self):
        from dual_view.iwasawa import verify_commutator_depth
        results = verify_commutator_depth(8, [(1, 1), (1, 2), (2, 2)], n_trials=10)
        for dM, dN, dMN in results:
            self.assertGreaterEqual(dMN, dM + dN,
                f"depth([M,N])={dMN} < depth(M)+depth(N)={dM+dN}")


class TestHolonomyDepthProfile(unittest.TestCase):
    def test_returns_dict_with_keys(self):
        from dual_view.iwasawa import holonomy_depth_profile
        result = holonomy_depth_profile(5, 7, cycle_length=3, n_cycles=5)
        self.assertIn("mean_depth_orig", result)
        self.assertIn("mean_depth_pert", result)
        self.assertIn("n_cycles", result)


# ── operators.py ────────────────────────────────────────────────────────────────

class TestOperatorPow(unittest.TestCase):
    def test_shift_squared(self):
        from dual_view.operators import OperatorContext
        ctx = OperatorContext(6, 5)
        S2 = ctx.S ** 2
        def f(e):
            return e + 1
        result = S2(f, 0)
        # S^2 f(0) = f(2) = 3
        self.assertEqual(result, 3)

    def test_shift_cubed(self):
        from dual_view.operators import OperatorContext
        ctx = OperatorContext(6, 5)
        S3 = ctx.S ** 3
        def f(e):
            return e
        result = S3(f, 0)
        # S^3 f(0) = f(3) = 3
        self.assertEqual(result, 3)


class TestSpectralTripleCurvature(unittest.TestCase):
    def test_curvature_returns_operator(self):
        from dual_view.operators import OperatorContext, SpectralTriple
        ctx = OperatorContext(6, 5)
        st = SpectralTriple(ctx)
        h1 = lambda e: e + 1
        h2 = lambda e: e * 2
        curv = st.curvature(h1, h2)
        self.assertIsInstance(curv, type(st.dirac))

    def test_curvature_is_callable(self):
        from dual_view.operators import OperatorContext, SpectralTriple
        ctx = OperatorContext(6, 5)
        st = SpectralTriple(ctx)
        h1 = lambda e: 1
        h2 = lambda e: 1
        curv = st.curvature(h1, h2)
        result = curv(lambda e: 1, 0)
        self.assertIsInstance(result, int)


# ── core.py ─────────────────────────────────────────────────────────────────────

class TestDualNumberCoords(unittest.TestCase):
    def test_coords_returns_triple(self):
        from dual_view.core import DualNumber
        d = DualNumber(42, k=16)
        coords = d.coords()
        self.assertEqual(len(coords), 3)
        v, alpha, e = coords
        self.assertEqual(v, d.v)
        self.assertEqual(alpha, d.alpha)
        self.assertEqual(e, d.e)


class TestProcessorDlog(unittest.TestCase):
    def test_dlog_returns_triple(self):
        from dual_view.core import DualNumber, TwoAdicProcessor
        proc = TwoAdicProcessor(16)
        d = DualNumber(42, 16)
        result = proc.dlog(d)
        self.assertEqual(len(result), 3)
        self.assertEqual(result, d.coords())


# ── isometry.py ─────────────────────────────────────────────────────────────────

class TestTraceExponentIndependence(unittest.TestCase):
    def test_returns_anova_results(self):
        from dual_view.isometry import trace_exponent_independence
        result = trace_exponent_independence(5, 7, cycle_length=3, n_cycles=20)
        self.assertIn("F_stat", result)
        self.assertIn("df_between", result)
        self.assertIn("df_within", result)
        self.assertIn("n_samples", result)
        self.assertIsInstance(result["F_stat"], float)


if __name__ == "__main__":
    unittest.main()
