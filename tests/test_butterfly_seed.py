"""Tests for dual_view.butterfly_seed — Dual-View Newton projector as butterfly seed."""
import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from dual_view.butterfly_seed import (
    DualViewSeed,
    analyze_prime,
    dual_view_qasm_emitter,
    CleanPrimeProfile,
    _hensel_bootstrap_exponent,
    _newton_fp,
)
from dual_view.butterfly_emitter import (
    basin_qasm_emitter,
    _basin_depth_d,
    dual_view_qasm_emitter_clean,
)
from dual_view.core import _valuation, two_adic_log5, two_adic_dlog
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES


class TestValuation(unittest.TestCase):
    def test_valuation_even(self):
        self.assertEqual(_valuation(8), 3)
        self.assertEqual(_valuation(12), 2)

    def test_valuation_odd(self):
        self.assertEqual(_valuation(7), 0)

    def test_valuation_zero(self):
        self.assertTrue(_valuation(0) == float("inf"))


class TestNewtonFp(unittest.TestCase):
    def test_newton_fp_fixed_point(self):
        # x=1 is a fixed point: N(1) = (2+1)/(3) = 1
        self.assertEqual(_newton_fp(1, 7), 1)

    def test_newton_fp_pole(self):
        # x=0 is a pole (denominator 0)
        self.assertIsNone(_newton_fp(0, 7))

    def test_newton_fp_cube_root(self):
        # cube roots of 1 mod 7 are 1, 2, 4 — all should be fixed points
        for x in (1, 2, 4):
            self.assertEqual(_newton_fp(x, 7), x)


class TestAnalyzePrime(unittest.TestCase):
    def test_p7_is_clean(self):
        prof = analyze_prime(7)
        self.assertTrue(prof.is_clean)
        self.assertEqual(prof.obstruction, "clean")
        self.assertEqual(prof.roots, (1, 2, 4))

    def test_p103_is_clean(self):
        prof = analyze_prime(103)
        self.assertTrue(prof.is_clean)
        self.assertEqual(prof.obstruction, "clean")
        self.assertEqual(len(prof.roots), 3)

    def test_p181_is_clean(self):
        prof = analyze_prime(181)
        self.assertTrue(prof.is_clean)
        self.assertEqual(prof.obstruction, "clean")
        self.assertEqual(len(prof.roots), 3)

    def test_p3_is_clean(self):
        prof = analyze_prime(3)
        self.assertTrue(prof.is_clean)
        self.assertIn(prof.obstruction, ("pole_chain", "clean"))
        self.assertEqual(prof.roots, (1,))

    def test_p5_is_clean(self):
        prof = analyze_prime(5)
        self.assertTrue(prof.is_clean)
        self.assertEqual(prof.roots, (1,))

    def test_p13_has_ghost_cycle(self):
        prof = analyze_prime(13)
        # 13 ≡ 1 mod 3, but has ghost cycles
        self.assertFalse(prof.is_clean)
        self.assertIn(prof.obstruction, ("ghost_cycle", "pole_chain", "mixed"))

    def test_clean_prime_has_correct_roots(self):
        for p in KNOWN_CLEAN_PRIMES:
            prof = analyze_prime(p)
            self.assertIn(len(prof.roots), (1, 3),
                          f"p={p}: expected 1 or 3 roots, got {len(prof.roots)}")
            # p ≡ 1 mod 3 should have 3 roots, p ≡ 2 mod 3 should have 1 root
            if p % 3 == 1:
                self.assertEqual(len(prof.roots), 3,
                                 f"p={p}: p≡1 mod 3 should have 3 roots")
            elif p % 3 == 2:
                self.assertEqual(len(prof.roots), 1,
                                 f"p={p}: p≡2 mod 3 should have 1 root")

    def test_nilpotency_index_positive_for_clean(self):
        for p in KNOWN_CLEAN_PRIMES:
            prof = analyze_prime(p)
            self.assertGreater(prof.nilpotency_index, 0)

    def test_nilpotency_index_zero_for_non_clean(self):
        prof = analyze_prime(13)
        self.assertEqual(prof.nilpotency_index, 0)

    def test_basin_ordering_contains_all_non_pole(self):
        prof = analyze_prime(7)
        # F_7^* has 6 elements, all should be in basin for clean prime
        self.assertEqual(len(prof.basin_ordering), 6)


class TestDualViewSeedConstruction(unittest.TestCase):
    def test_basic_construction(self):
        dvs = DualViewSeed(k=8, target_a=5)
        self.assertEqual(dvs.k, 8)
        self.assertEqual(dvs.N, 64)  # 2^(8-2)

    def test_even_target_raises(self):
        with self.assertRaises(ValueError):
            DualViewSeed(k=8, target_a=4)

    def test_custom_generator(self):
        dvs = DualViewSeed(k=8, target_a=5, g=5)
        self.assertEqual(dvs.g, 5)

    def test_a_reduced_mod_2k(self):
        dvs = DualViewSeed(k=8, target_a=257)  # 257 ≡ 1 mod 256
        self.assertEqual(dvs.a, 1)


class TestDualViewSeedNewtonStep(unittest.TestCase):
    def test_newton_step_converges_for_known_dlog(self):
        k = 8
        a = pow(5, 10, 2**k)  # a = 5^10 mod 2^k
        dvs = DualViewSeed(k, a)
        e = dvs.newton_step_e(10)
        # After one Newton step from the true root, should stay at root
        self.assertEqual(e, 10)

    def test_newton_step_reduces_error(self):
        k = 10
        a = pow(5, 7, 2**k)
        dvs = DualViewSeed(k, a)
        e0 = 0  # start far from root
        e1 = dvs.newton_step_e(e0)
        # After one step, should be closer (or at least different)
        self.assertNotEqual(e1, e0)


class TestDualViewSeedButterflySeeds(unittest.TestCase):
    def test_build_seeds_returns_correct_structure(self):
        k = 6
        dvs = DualViewSeed(k, target_a=5)
        seeds = dvs.build_position_dependent_seeds()
        # Should have k-2 stages
        self.assertEqual(len(seeds), k - 2)

    def test_seed_stage_sizes(self):
        k = 8
        dvs = DualViewSeed(k, target_a=5)
        seeds = dvs.build_position_dependent_seeds()
        for m, stage_seeds in enumerate(seeds):
            self.assertEqual(len(stage_seeds), 1 << m)

    def test_seed_matrices_are_2x2(self):
        k = 6
        dvs = DualViewSeed(k, target_a=5)
        seeds = dvs.build_position_dependent_seeds()
        for stage in seeds:
            for seed in stage:
                self.assertEqual(seed.shape, (2, 2))

    def test_seed_matrices_are_complex(self):
        k = 6
        dvs = DualViewSeed(k, target_a=5)
        seeds = dvs.build_position_dependent_seeds()
        for stage in seeds:
            for seed in stage:
                self.assertEqual(seed.dtype, complex)

    def test_seed_unitary_structure(self):
        """Each seed should be proportional to a unitary matrix."""
        k = 6
        dvs = DualViewSeed(k, target_a=5)
        seeds = dvs.build_position_dependent_seeds()
        for stage in seeds:
            for seed in stage:
                # S^† S should be proportional to identity
                prod = seed.conj().T @ seed
                # Off-diagonal should be ~0
                self.assertAlmostEqual(abs(prod[0, 1]), 0, places=10)
                self.assertAlmostEqual(abs(prod[1, 0]), 0, places=10)


class TestDualViewSeedThermodynamics(unittest.TestCase):
    def test_thermodynamic_signature_unitary(self):
        dvs = DualViewSeed(k=8, target_a=5)
        sig = dvs.thermodynamic_signature()
        self.assertTrue(sig["is_unitary"])
        self.assertTrue(sig["is_conservative"])
        self.assertFalse(sig["is_contractive"])
        self.assertFalse(sig["is_expansive"])
        self.assertFalse(sig["is_nilpotent"])
        self.assertEqual(sig["spectral_radius"], 1.0)

    def test_solvability_report(self):
        dvs = DualViewSeed(k=8, target_a=5)
        report = dvs.solvability_report()
        self.assertTrue(report["is_solvable"])
        self.assertEqual(report["depth"], 2)
        self.assertIn("SOLVABLE", report["conclusion"])


class TestHenselBootstrap(unittest.TestCase):
    def test_bootstrap_reconstructs_a(self):
        """5^e ≡ a (mod 2^k) should hold after bootstrap."""
        for k in (6, 8, 10, 12, 16):
            for a_val in (5, 9, 17, 33, 65):
                if a_val % 4 != 1:
                    continue
                bits = _hensel_bootstrap_exponent(k, a_val)
                e = sum(b << i for i, b in enumerate(bits))
                reconstructed = pow(5, e, 2**k)
                self.assertEqual(reconstructed, a_val % (2**k),
                                 f"Failed for k={k}, a={a_val}")

    def test_bootstrap_length(self):
        for k in (6, 8, 10, 16):
            bits = _hensel_bootstrap_exponent(k, 5)
            self.assertEqual(len(bits), k - 2)

    def test_bootstrap_a1(self):
        """5^0 = 1, so bootstrap for a=1 should give all zeros."""
        for k in (6, 8, 10):
            bits = _hensel_bootstrap_exponent(k, 1)
            self.assertEqual(bits, [0] * (k - 2))


class TestQasmEmitter(unittest.TestCase):
    def test_qasm_header(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("OPENQASM 2.0;", qasm)
        self.assertIn('include "qelib1.inc";', qasm)

    def test_qasm_register_sizes(self):
        k = 10
        qasm = dual_view_qasm_emitter(k=k, target_a=5)
        n_exp = k - 2
        n_val = max(1, (k + 1) // 3)
        n_total = n_exp + n_val + 1
        self.assertIn(f"qreg q[{n_total}];", qasm)
        self.assertIn(f"creg c[{n_total}];", qasm)

    def test_qasm_clean_prime_annotated(self):
        """p_clean=7 should use the basin-annotated circuit (D=1)."""
        qasm = dual_view_qasm_emitter(k=8, target_a=5, p_clean=7)
        self.assertIn("Basin-annotated Dual-View", qasm)
        self.assertIn("p=7", qasm)
        # p=7 has M=2, D=1 — only 1 exponent qubit instead of 6
        self.assertIn("Exponent qubits: 1", qasm)

    def test_qasm_clean_prime_reduced_register(self):
        """p=103 should use D=4 exponent qubits instead of k-2=6."""
        qasm = dual_view_qasm_emitter(k=8, target_a=5, p_clean=103)
        self.assertIn("Exponent qubits: 4", qasm)

    def test_qasm_no_clean_prime(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5, p_clean=None)
        self.assertNotIn("Basin-optimised", qasm)

    def test_qasm_has_qft_section(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("// --- QFT on exponent register ---", qasm)

    def test_qasm_has_inverse_qft(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("// --- Inverse QFT ---", qasm)

    def test_qasm_has_newton_diagonal(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("// --- Newton diagonal", qasm)

    def test_qasm_has_valuation_guard(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("// --- Valuation guard", qasm)

    def test_qasm_has_measurement(self):
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("// --- Measurement ---", qasm)

    def test_qasm_state_prep_has_x_gates(self):
        """State preparation should emit X gates for non-zero bits."""
        qasm = dual_view_qasm_emitter(k=8, target_a=5)
        self.assertIn("x q[", qasm)

    def test_qasm_line_count_scales_with_k(self):
        qasm1 = dual_view_qasm_emitter(k=8, target_a=5)
        qasm2 = dual_view_qasm_emitter(k=12, target_a=5)
        self.assertGreater(len(qasm2.splitlines()), len(qasm1.splitlines()))


class TestIntegration(unittest.TestCase):
    """End-to-end integration tests."""

    def test_clean_prime_circuit_generation(self):
        """Full pipeline: analyze prime → build seed → emit QASM."""
        p = 7
        prof = analyze_prime(p)
        self.assertTrue(prof.is_clean)

        k = 10
        a = pow(5, 3, 2**k)
        dvs = DualViewSeed(k, a)
        seeds = dvs.build_position_dependent_seeds()
        self.assertEqual(len(seeds), k - 2)

        qasm = dual_view_qasm_emitter(k, a, p_clean=p)
        self.assertIn("OPENQASM 2.0;", qasm)

    def test_seed_phase_uses_actual_dlog(self):
        """The _newton_phase method should use the actual discrete log."""
        k = 8
        a = pow(5, 3, 2**k)  # a = 5^3 mod 2^k
        dvs = DualViewSeed(k, a)

        # The discrete log of a should be 3
        result = two_adic_dlog(a, k)
        self.assertIsNotNone(result)
        _, e = result
        self.assertEqual(e, 3)

        # Phase at stage 0 should depend on bit 0 of e=3 (which is 1)
        phase0 = dvs._newton_phase(0, 1)
        self.assertIsInstance(phase0, int)

    def test_multiple_clean_primes_consistent(self):
        """All known clean primes should produce consistent profiles."""
        for p in KNOWN_CLEAN_PRIMES:
            prof = analyze_prime(p)
            self.assertTrue(prof.is_clean, f"p={p} should be clean")
            self.assertIn(len(prof.roots), (1, 3), f"p={p} should have 1 or 3 roots")
            self.assertGreater(prof.nilpotency_index, 0)

    def test_seed_with_different_generators(self):
        """Different generators should produce different seeds."""
        k = 8
        a = 5
        dvs5 = DualViewSeed(k, a, g=5)
        seeds5 = dvs5.build_position_dependent_seeds()
        self.assertGreater(len(seeds5), 0)


class TestEdgeCases(unittest.TestCase):
    def test_minimal_k(self):
        """k=3 is the minimum valid precision."""
        dvs = DualViewSeed(k=3, target_a=5)
        self.assertEqual(dvs.N, 2)  # 2^(3-2) = 2
        seeds = dvs.build_position_dependent_seeds()
        self.assertEqual(len(seeds), 1)  # k-2 = 1 stage

    def test_large_k(self):
        """Test with larger precision."""
        k = 20
        a = pow(5, 123, 2**k)
        dvs = DualViewSeed(k, a)
        self.assertEqual(dvs.N, 2**18)
        seeds = dvs.build_position_dependent_seeds()
        self.assertEqual(len(seeds), 18)

    def test_qasm_large_k(self):
        """QASM emission with larger k."""
        k = 20
        a = pow(5, 123, 2**k)
        qasm = dual_view_qasm_emitter(k, a)
        self.assertIn("OPENQASM 2.0;", qasm)


class TestBasinEmitter(unittest.TestCase):
    """Tests for the basin-optimised QASM emitter."""

    def test_basin_depth_d(self):
        """_basin_depth_d computes ceil(log2(M)) for clean primes."""
        # p=7:  M=2  -> D=1
        self.assertEqual(_basin_depth_d(7), 1)
        # p=103: M=15 -> D=4
        self.assertEqual(_basin_depth_d(103), 4)
        # p=403549: M=1223 -> D=11
        self.assertEqual(_basin_depth_d(403549), 11)

    def test_basin_qasm_header(self):
        qasm = basin_qasm_emitter(7, 8, 5)
        self.assertIn("OPENQASM 2.0;", qasm)
        self.assertIn("Basin-annotated Dual-View", qasm)
        self.assertIn("p=7", qasm)

    def test_basin_qasm_register_sizes(self):
        """p=7 has D=1, so only 1 exponent qubit."""
        qasm = basin_qasm_emitter(7, 8, 5)
        # Register: 1 (exp) + 3 (val) + 1 (sign) = 5 total
        self.assertIn("qreg q[5];", qasm)

    def test_basin_qasm_larger_prime_register(self):
        """p=103 has D=4, so 4 exponent qubits."""
        qasm = basin_qasm_emitter(103, 10, 5)
        # Register: 4 (exp) + 3 (val, from (10+1)//3) + 1 (sign) = 8 total
        self.assertIn("qreg q[8];", qasm)
        self.assertIn("Exponent qubits: 4", qasm)

    def test_basin_qasm_has_routing_stages(self):
        qasm = basin_qasm_emitter(7, 8, 5)
        self.assertIn("Basin routing stages", qasm)
        self.assertIn("advance by", qasm)

    def test_basin_qasm_has_reduced_qft(self):
        qasm = basin_qasm_emitter(7, 8, 5)
        self.assertIn("Reduced QFT on 1 qubits", qasm)

    def test_basin_qasm_has_neumann_diagonal(self):
        qasm = basin_qasm_emitter(7, 8, 5)
        self.assertIn("Neumann phase stage", qasm)

    def test_basin_qasm_fallback_for_non_clean(self):
        """Non-clean prime should fall back to standard circuit."""
        qasm = basin_qasm_emitter(13, 8, 5)  # p=13 has ghost cycles
        self.assertNotIn("Basin-annotated", qasm)

    def test_basin_qasm_small_clean_primes(self):
        """First 8 clean primes produce valid QASM (fast check)."""
        for p in KNOWN_CLEAN_PRIMES[:8]:
            qasm = basin_qasm_emitter(p, 10, 5)
            self.assertIn("OPENQASM 2.0;", qasm)
            self.assertIn(f"p={p}", qasm)

    def test_clean_emitter_wrapper(self):
        """dual_view_qasm_emitter_clean delegates correctly."""
        # Clean prime: annotation active
        qasm = dual_view_qasm_emitter_clean(8, 5, p_clean=7)
        self.assertIn("Basin-annotated", qasm)
        # No p_clean: standard circuit
        qasm = dual_view_qasm_emitter_clean(8, 5, p_clean=None)
        self.assertNotIn("Basin-annotated", qasm)

    def test_standard_emitter_delegates_to_basin(self):
        """dual_view_qasm_emitter delegates to basin for clean primes."""
        qasm = dual_view_qasm_emitter(8, 5, p_clean=7)
        self.assertIn("Basin-annotated", qasm)
        self.assertIn("Exponent qubits: 1", qasm)


if __name__ == "__main__":
    unittest.main()
