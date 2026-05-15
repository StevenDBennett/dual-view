"""Tests for dual_view.exponent and dual_view.operators."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.exponent import ExponentSpace
from dual_view.operators import OperatorContext, SpectralTriple, NewtonProjector
from dual_view.core import _mask


class TestExponentSpace(unittest.TestCase):
    def test_lift_zero(self):
        es = ExponentSpace(5, 8)
        self.assertEqual(es.lift(0), 1)

    def test_lift_one(self):
        es = ExponentSpace(5, 8)
        self.assertEqual(es.lift(1), 5)

    def test_lift_matches_pow(self):
        es = ExponentSpace(5, 12)
        for e in range(8):
            self.assertEqual(es.lift(e), pow(5, e, 1 << 12))

    def test_periodicity(self):
        es = ExponentSpace(5, 8)
        N = 1 << 6
        self.assertEqual(es.lift(0), es.lift(N))

    def test_difference_operator(self):
        es = ExponentSpace(5, 8)
        f = lambda e: pow(5, e, 1 << 8)
        for e in range(4):
            diff = es.difference(f, e)
            expected = (f(e + 1) - f(e)) & _mask(8)
            self.assertEqual(diff, expected)

    def test_integral(self):
        es = ExponentSpace(5, 8)
        f = lambda e: 1
        integral = es.integrate(f)
        self.assertEqual(integral, es.N)


class TestOperatorContext(unittest.TestCase):
    def setUp(self):
        self.ctx = OperatorContext(8, 5)
        self.g_pow = [pow(5, e, 1 << 8) for e in range(self.ctx.N)]
        self.f = lambda e: self.g_pow[e]

    def test_identity(self):
        for e in range(4):
            self.assertEqual(self.ctx.I(self.f, e), self.f(e))

    def test_shift(self):
        for e in range(4):
            self.assertEqual(self.ctx.S(self.f, e), self.f((e + 1) % self.ctx.N))

    def test_difference_identity(self):
        for e in range(4):
            diff = (self.ctx.I - self.ctx.S)(self.f, e)
            expected = self.f(e) - self.f((e + 1) % self.ctx.N)
            self.assertEqual(diff, expected & self.ctx.mask)

    def test_multiplication(self):
        h = lambda e: e + 1
        M = self.ctx.multiply_by(h, "M_test")
        for e in range(4):
            self.assertEqual(M(self.f, e), (h(e) * self.f(e)) & self.ctx.mask)

    def test_average_operator(self):
        total = sum(self.f(e) for e in range(self.ctx.N)) & self.ctx.mask
        self.assertEqual(self.ctx.avg(self.f, 0), total)

    def test_operator_add(self):
        op = self.ctx.I + self.ctx.S
        for e in range(4):
            self.assertEqual(
                op(self.f, e),
                (self.f(e) + self.f((e + 1) % self.ctx.N)) & self.ctx.mask,
            )

    def test_operator_compose(self):
        SS = self.ctx.S * self.ctx.S
        for e in range(4):
            self.assertEqual(
                SS(self.f, e),
                self.f((e + 2) % self.ctx.N),
            )


class TestSpectralTriple(unittest.TestCase):
    def setUp(self):
        self.ctx = OperatorContext(8, 5)
        self.st = SpectralTriple(self.ctx)
        self.h1 = lambda e: 1
        self.h2 = lambda e: e + 1

    def test_one_form_non_trivial(self):
        w = self.st.one_form(self.h1)
        g_pow = [pow(5, e, 1 << 8) for e in range(self.ctx.N)]
        f = lambda e: g_pow[e]
        val = w(f, 0)
        self.assertIsInstance(val, int)

    def test_gauge_transform(self):
        f = lambda e: 1
        w = self.st.one_form(f)
        D_prime = self.st.gauge_transform(w)
        g_pow = [pow(5, e, 1 << 8) for e in range(self.ctx.N)]
        fn = lambda e: g_pow[e]
        val = D_prime(fn, 0)
        self.assertIsInstance(val, int)

    def test_trace_int(self):
        f = lambda e: 1
        tr = self.st.trace(f)
        self.assertEqual(tr, self.ctx.N & self.ctx.mask)


class TestNewtonProjector(unittest.TestCase):
    def test_project_converges(self):
        ctx = OperatorContext(8, 5)
        target_e = 7
        a = pow(5, target_e, 1 << 8)
        proj = NewtonProjector(ctx, a)
        result = proj.project_point(0)
        self.assertEqual(result, target_e)

    def test_negative_sector(self):
        ctx = OperatorContext(8, 5)
        target_e = 7
        a = (-pow(5, target_e, 1 << 8)) & _mask(8)
        proj = NewtonProjector(ctx, a)
        result = proj.project_point(0)
        self.assertIsInstance(result, int)


if __name__ == "__main__":
    unittest.main()
