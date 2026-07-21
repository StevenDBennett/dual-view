"""
dual-view (v1.0.0)
===================
dual-view: a mathematical framework for 2-adic number systems, with diagnostics for quantized neural network weights.

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
regularization   — GhostMap stability scores (deprecated for training; use thermodynamics)
crt              — CRT extension to composite moduli
nonabelian       — GL(2) gauge theory for matrix-valued weights
scaling          — Float-to-int quantization scaling
visualise        — Cliff matrix rendering and ASCII heatmaps
butterfly        — Kronecker factor cliff scoring
separation       — Trajectory Separation Theorem
fourier          — Discrete Fourier analysis on exponent domain
padic_roots      — Multi-order p-adic root finding
iwasawa          — GL(2) congruence filtration and LDU decomposition
iwasawa_algebra  — Iwasawa algebra Z_2[[G]] and profinite filtered modules
mersenne         — Mersenne Ghost Theorem and bootstrap optimality
isometry         — Exponential isometry and operator algebra theorems
newton_dynamics  — p-adic Newton dynamics for N(x) = (2x³+1)/(3x²)
butterfly_seed   — Dual-view Newton projector as butterfly-compilable seed
training         — PyTorch quantized MLP with ghost reg. (if torch avail)
demo             — Runnable demonstration suite
"""

from ._version import __version__

from .core import (
    modinv_newton,
    two_adic_log5,
    two_adic_dlog,
    dual_add,
    DualNumber,
    TwoAdicProcessor,
    padic_exp,
    padic_log,
)

from .exponent import ExponentSpace
from .mahler import MahlerCalculus
from .operators import OperatorContext, SpectralTriple, NewtonProjector
from .gauge import cycle_product, spectral_det, det_coordinates, tidal_scalar, GaugeLayer
from .basin import BasinExplorer, precision_sweep
from .thermodynamics import SeedThermodynamics
from .regularization import GhostMap, local_ratio_gradient, ghost_penalty

from .crt import CRTDualNumber, CRTDualProcessor, combined_stability
from .nonabelian import NonAbelianCRTDual, phase_alignment_experiment
from .scaling import scale_weights, auto_scale, common_scales
from .visualise import (
    cliff_matrix, sector_matrix, valuation_matrix,
    print_cliff_ascii, cliff_stats_by_layer, show_dual_bits,
)
from .butterfly import KroneckerCliffScorer, semiring_cliff_score
from .separation import newton_trajectory, separation_step, predicted_separation, step_count_profile
from .fourier import step_count_fn, analytic_step_count, dft, power_spectrum, dyadic_coefficients, fourier_summary
from .padic_roots import lift_root, newton_step, halley_step, newton2_step, newton3_step, convergence_profile
from .iwasawa import congruence_depth, filtration_residue, ldu_decompose, matrix_coordinates, matrix_commutator, MatrixCoordinates

from .newton_dynamics import (
    poly_mul, poly_add, poly_scalar_mul, poly_pow, poly_divmod,
    mobius, compute_iterates, dynatomic_polynomial,
    is_cube, tonelli_shanks, check_quadratic_cube_roots,
    COEFFS_PERIOD4, COEFFS_PERIOD5,
    MULTIPLIERS_PERIOD4, MULTIPLIERS_PERIOD5,
    load_period6_coefficients, PERIOD6_PREDICTED,
)

from .mersenne import (
    mersenne_coordinates, mersenne_cliff_table,
    cliff_constant, cliff_formula, mersenne_cliff_theorem,
    cliff_constant_unified, dlog_with_lut,
)

from .isometry import verify_isometry, isometry_pair_test, isometry_summary, verify_operator_algebra
from .butterfly_seed import DualViewSeed, analyze_prime, CleanPrimeProfile
from .training import QuantizedMLP
from . import demo

from .iwasawa_algebra import IwasawaElement, IwasawaAlgebra, ProModule

__all__ = [
    "modinv_newton", "two_adic_log5", "two_adic_dlog", "dual_add",
    "DualNumber", "TwoAdicProcessor", "padic_exp", "padic_log",
    "ExponentSpace", "MahlerCalculus",
    "OperatorContext", "SpectralTriple", "NewtonProjector",
    "cycle_product", "spectral_det", "det_coordinates", "tidal_scalar", "GaugeLayer",
    "BasinExplorer", "precision_sweep",
    "SeedThermodynamics",
    "GhostMap", "local_ratio_gradient", "ghost_penalty",
    "CRTDualNumber", "CRTDualProcessor", "combined_stability",
    "NonAbelianCRTDual", "phase_alignment_experiment",
    "scale_weights", "auto_scale", "common_scales",
    "cliff_matrix", "sector_matrix", "valuation_matrix",
    "print_cliff_ascii", "cliff_stats_by_layer", "show_dual_bits",
    "KroneckerCliffScorer", "semiring_cliff_score",
    "newton_trajectory", "separation_step", "predicted_separation", "step_count_profile",
    "step_count_fn", "analytic_step_count", "dft", "power_spectrum", "dyadic_coefficients", "fourier_summary",
    "lift_root", "newton_step", "halley_step", "newton2_step", "newton3_step", "convergence_profile",
    "congruence_depth", "filtration_residue", "ldu_decompose", "matrix_coordinates", "matrix_commutator", "MatrixCoordinates",
    "mersenne_coordinates", "mersenne_cliff_table", "cliff_constant", "cliff_formula",
    "mersenne_cliff_theorem", "cliff_constant_unified", "dlog_with_lut",
    "poly_mul", "poly_add", "poly_scalar_mul", "poly_pow", "poly_divmod",
    "mobius", "compute_iterates", "dynatomic_polynomial",
    "is_cube", "tonelli_shanks", "check_quadratic_cube_roots",
    "COEFFS_PERIOD4", "COEFFS_PERIOD5",
    "MULTIPLIERS_PERIOD4", "MULTIPLIERS_PERIOD5",
    "load_period6_coefficients", "PERIOD6_PREDICTED",
    "verify_isometry", "isometry_pair_test", "isometry_summary", "verify_operator_algebra",
    "DualViewSeed", "analyze_prime", "CleanPrimeProfile",
    "IwasawaElement", "IwasawaAlgebra", "ProModule",
    "QuantizedMLP",
    "demo",
]
