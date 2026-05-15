"""
dual-view (v1.0.0)
===================
2-adic dual-view diagnostics for quantized neural network weights.

Every odd integer modulo 2^k decomposes uniquely as a dual-view
coordinate triple (v, α, e) where:

    n = 2^v · (-1)^α · 5^e   (mod 2^k)

This package provides the complete mathematical framework for
analysing the 2-adic structure of weight matrices, detecting
quantization cliffs, and regularising training.

Submodules
----------
core             — DualNumber, modular inverse, discrete log
exponent         — Additive coordinate chart on Z/2^(k-2)
operators        — Symbolic operator algebra (shift, difference)
gauge            — Gauge invariants for weighted cyclic operators
basin            — Newton basin analysis and ghost detection
thermodynamics   — Graded weight stability diagnostics
regularization   — Ghost-aware regularisation for NN training
crt              — CRT extension to composite moduli
nonabelian       — GL(2) gauge theory for matrix-valued weights
scaling          — Float-to-int quantization scaling
visualise        — Cliff matrix rendering and ASCII heatmaps
butterfly        — Kronecker factor cliff scoring
separation       — Trajectory Separation Theorem
fourier          — Discrete Fourier analysis on exponent domain
padic_roots      — Multi-order p-adic root finding
iwasawa          — GL(2) congruence filtration and LDU decomposition
mersenne         — Mersenne Ghost Theorem and bootstrap optimality
isometry         — Exponential isometry and operator algebra theorems
training         — PyTorch quantized MLP with ghost reg. (if torch avail)
demo             — Runnable demonstration suite
"""

from ._version import __version__

from .core import (
    modinv_newton,
    two_adic_log5,
    two_adic_dlog,
    DualNumber,
    TwoAdicProcessor,
    run_all_tests,
)

from .exponent import ExponentSpace
from .operators import OperatorContext, SpectralTriple, NewtonProjector
from .gauge import cycle_product, spectral_det, det_coordinates, tidal_scalar, GaugeLayer
from .basin import BasinExplorer, precision_sweep
from .thermodynamics import SeedThermodynamics
from .regularization import GhostMap, local_ratio_gradient, ghost_penalty

from .crt import CRTDualNumber, CRTDualProcessor, combined_stability
from .nonabelian import NonAbelianCRTDual, ramp_break_strength, phase_alignment_experiment
from .scaling import scale_weights, auto_scale, common_scales
from .visualise import (
    cliff_matrix, sector_matrix, valuation_matrix,
    print_cliff_ascii, cliff_stats_by_layer,
)
from .butterfly import KroneckerCliffScorer, semiring_cliff_score
from .separation import (
    newton_trajectory, separation_step, predicted_separation,
    verify_separation, ultrametric_ball_tree, step_count_profile,
)
from .fourier import (
    step_count_fn, analytic_step_count, dft, power_spectrum,
    dyadic_coefficients, analytic_coefficients, fourier_summary,
    ultrametric_uncertainty,
)
from .padic_roots import (
    newton_step, halley_step, newton2_step, newton3_step,
    convergence_profile, compare_methods, verify_order,
    newton_correction_uniformity, popcount_compression,
)
from .iwasawa import (
    congruence_depth, filtration_residue, ldu_decompose,
    matrix_coordinates, holonomy_depth_profile, filtration_portrait,
    matrix_commutator, verify_commutator_depth, MatrixCoordinates,
)
from .mersenne import (
    mersenne_coordinates, verify_core_identity, mersenne_cliff_table,
    bootstrap_cost, optimal_bootstrap, compare_bootstrap_strategies,
    dlog_with_lut, verify_lut_dlog,
)
from .isometry import (
    verify_isometry, isometry_pair_test, isometry_summary,
    verify_operator_algebra, trace_alpha_independence,
    trace_exponent_independence, exponent_valuation_profile,
)

__all__ = [
    # core
    "modinv_newton", "two_adic_log5", "two_adic_dlog",
    "DualNumber", "TwoAdicProcessor", "run_all_tests",
    # exponent
    "ExponentSpace",
    # operators
    "OperatorContext", "SpectralTriple", "NewtonProjector",
    # gauge
    "cycle_product", "spectral_det", "det_coordinates", "tidal_scalar",
    "GaugeLayer",
    # basin
    "BasinExplorer", "precision_sweep",
    # thermodynamics
    "SeedThermodynamics",
    # regularization
    "GhostMap", "local_ratio_gradient", "ghost_penalty",
    # crt
    "CRTDualNumber", "CRTDualProcessor", "combined_stability",
    # nonabelian
    "NonAbelianCRTDual", "ramp_break_strength", "phase_alignment_experiment",
    # scaling
    "scale_weights", "auto_scale", "common_scales",
    # visualise
    "cliff_matrix", "sector_matrix", "valuation_matrix",
    "print_cliff_ascii", "cliff_stats_by_layer",
    # butterfly
    "KroneckerCliffScorer", "semiring_cliff_score",
    # separation
    "newton_trajectory", "separation_step", "predicted_separation",
    "verify_separation", "ultrametric_ball_tree", "step_count_profile",
    # fourier
    "step_count_fn", "analytic_step_count", "dft", "power_spectrum",
    "dyadic_coefficients", "analytic_coefficients", "fourier_summary",
    "ultrametric_uncertainty",
    # padic_roots
    "newton_step", "halley_step", "newton2_step", "newton3_step",
    "convergence_profile", "compare_methods", "verify_order",
    "newton_correction_uniformity", "popcount_compression",
    # iwasawa
    "congruence_depth", "filtration_residue", "ldu_decompose",
    "matrix_coordinates", "holonomy_depth_profile", "filtration_portrait",
    "matrix_commutator", "verify_commutator_depth", "MatrixCoordinates",
    # mersenne
    "mersenne_coordinates", "verify_core_identity", "mersenne_cliff_table",
    "bootstrap_cost", "optimal_bootstrap", "compare_bootstrap_strategies",
    "dlog_with_lut", "verify_lut_dlog",
    # isometry
    "verify_isometry", "isometry_pair_test", "isometry_summary",
    "verify_operator_algebra", "trace_alpha_independence",
    "trace_exponent_independence", "exponent_valuation_profile",
]
