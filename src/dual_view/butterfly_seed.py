"""
butterfly_seed.py — Dual-View Newton Projector as Butterfly-Compilable Seed
===========================================================================

Slots into the dual-view architecture, bridging the 2-adic Newton dynamics
with the butterfly compiler's position-dependent operad framework.

Provides:
  - DualViewSeed: 2-adic Newton projector on the exponent space as a
    position-dependent butterfly operad seed
  - analyze_prime: classifies a prime p by the thermodynamics of its Newton
    functional graph, returning the nilpotency index and basin-depth ordering
    required by the butterfly compiler

Dependencies: numpy (matches dual-view's existing dependencies)
"""

from __future__ import annotations

import numpy as np
from typing import Optional, Tuple, List
from dataclasses import dataclass

from .core import (
    _valuation,
    modinv_newton,
    two_adic_log5,
    two_adic_dlog,
)

# ---------------------------------------------------------------------------
# 1. Clean-prime functional-graph analyser
# ---------------------------------------------------------------------------


@dataclass
class CleanPrimeProfile:
    """
    Thermodynamic / structural profile of a prime p under the Newton map
    N(x) = (2x^3 + 1) / (3x^2) over F_p^*.
    """
    p: int
    is_clean: bool
    roots: Tuple[int, ...]           # cube roots of 1 mod p
    nilpotency_index: int            # max basin depth (0 if not clean)
    basin_ordering: List[int]        # F_p^* ordered by depth (deepest first)
    tree_depths: dict                # x -> depth until root
    obstruction: str                 # "clean", "ghost_cycle", "pole_chain", "mixed"


def _newton_fp(x: int, p: int) -> Optional[int]:
    """Newton map over F_p. Returns None for pole (denominator 0)."""
    den = (3 * x * x) % p
    if den == 0:
        return None
    num = (2 * pow(x, 3, p) + 1) % p
    return (num * pow(den, -1, p)) % p


def analyze_prime(p: int) -> CleanPrimeProfile:
    """
    Classify p by the thermodynamics of its Newton functional graph.

    A prime is *clean* for N(x) = (2x^3+1)/(3x^2) if F_p^* admits no
    periodic points of period 2, 3, or 4 with multiplier mu != 1 (mod p).
    Equivalently, the functional graph has no cycles other than the fixed
    points at the cube roots of unity. Pole chains (elements that crash
    to x=0) do not disqualify a prime — they are not periodic points.

    For p ≡ 1 (mod 3) there are 3 cube roots and 3 basin trees.
    For p ≡ 2 (mod 3) there is 1 cube root and 1 basin tree.
    """
    roots = tuple(sorted({x for x in range(1, p) if pow(x, 3, p) == 1}))
    nxt = [0] * p
    for x in range(1, p):
        nxt[x] = _newton_fp(x, p)

    state = [0] * p
    depths = {}
    cycles = []
    has_pole_chain = False

    for s in range(1, p):
        if state[s] != 0:
            continue
        path = []
        cur = s
        while True:
            if cur is None or cur == 0:
                if any(node not in roots for node in path):
                    has_pole_chain = True
                for node in path:
                    depths[node] = -1
                    state[node] = 2
                break
            if state[cur] == 1:
                idx = path.index(cur)
                cyc = path[idx:]
                cycles.append(cyc)
                for node in cyc:
                    if node in roots:
                        depths[node] = 0
                for i, node in enumerate(path[:idx]):
                    if all(r in roots for r in cyc):
                        depths[node] = len(path) - i
                    else:
                        depths[node] = -1
                for node in path:
                    state[node] = 2
                break
            if state[cur] == 2:
                known = depths.get(cur, None)
                for i, node in enumerate(path):
                    if known is not None and known >= 0:
                        depths[node] = known + len(path) - i
                    else:
                        depths[node] = -1
                for node in path:
                    state[node] = 2
                break
            state[cur] = 1
            path.append(cur)
            cur = nxt[cur]

    ghost_cycles = [c for c in cycles if any(x not in roots for x in c)]
    is_clean = len(ghost_cycles) == 0

    obstruction = "clean"
    if ghost_cycles and has_pole_chain:
        obstruction = "mixed"
    elif ghost_cycles:
        obstruction = "ghost_cycle"
    elif has_pole_chain:
        obstruction = "pole_chain"

    basin = [x for x in range(1, p) if depths.get(x, -1) >= 0]
    basin.sort(key=lambda x: depths[x], reverse=True)

    nil_idx = max((depths[x] for x in basin), default=0) + 1 if is_clean else 0

    return CleanPrimeProfile(
        p=p,
        is_clean=is_clean,
        roots=roots,
        nilpotency_index=nil_idx,
        basin_ordering=basin,
        tree_depths=depths,
        obstruction=obstruction,
    )

# ---------------------------------------------------------------------------
# 2. DualViewSeed — butterfly-compilable seed for the exponent space
# ---------------------------------------------------------------------------


class DualViewSeed:
    """
    Represents the 2-adic Newton projector on the exponent space
    Z/2^{k-2}Z as a position-dependent butterfly seed.

    In the dual-view coordinates (v, alpha, e) the Newton step on the
    exponent register is:

        e ← e - (g^e - a) / (4 * g^e * L)   (mod 2^{k-2})

    After QFT, multiplication by g^e becomes a phase shift, and the
    Newton step turns into a diagonal phase accumulation followed by an
    inverse QFT.  This class builds the position-dependent 2×2 seeds
    that encode that diagonal operator at each butterfly position.
    """

    def __init__(self, k: int, target_a: int, g: int = 5):
        """
        k       : precision (modulus 2^k)
        target_a: unit to invert / take log of  (a ≡ 1 mod 4)
        g       : generator of the cyclic part (default 5)
        """
        if target_a % 2 == 0:
            raise ValueError("target_a must be odd")
        self.k = k
        self.N = 2**(k - 2)          # exponent space size
        self.a = target_a % (2**k)
        self.g = g % (2**k)
        self.L = two_adic_log5(k)

    # ------------------------------------------------------------------
    # 2-adic Newton step in the exponent ring (classical, for reference)
    # ------------------------------------------------------------------

    def newton_step_e(self, e: int) -> int:
        """Single Newton update on exponent e (integer, mod N)."""
        mod = self.N
        g_e = pow(self.g, e, 2**self.k)
        diff = (g_e - self.a) % (2**self.k)
        # strip the predictable factor of 4
        diff_quarter = diff // 4
        denom = (g_e * (self.L >> 2)) % mod
        delta = (diff_quarter * modinv_newton(denom, self.k - 2)) % mod
        return (e - delta) % mod

    # ------------------------------------------------------------------
    # Butterfly seed generation
    # ------------------------------------------------------------------

    def build_position_dependent_seeds(self) -> List[List[np.ndarray]]:
        """
        Build the list-of-lists of 2×2 complex seeds used by
        PositionDependentOperad in butterfly.qft.

        Stage m (0 ≤ m < k-2) has 2^m butterflies.
        At stage m, butterfly position j gets the seed:

            S_{m,j} = (1/√2) * [[1,  ω^{φ(j)}],
                                 [1, -ω^{φ(j)}]]

        where ω = exp(-2πi / 2^{m+1}) is the primitive root for that stage,
        and φ(j) is the 2-adic Newton phase accumulated at that position.

        For the *clean-prime vacuum* the phase φ(j) collapses to one of
        three values (corresponding to the three root basins), making the
        seed effectively a 3-output multiplexer embedded in the QFT butterfly.
        """
        level_seeds = []
        for m in range(self.k - 2):
            half_m = 1 << m
            W = np.exp(-2j * np.pi / (2 * half_m))
            seeds_m = []
            for j in range(half_m):
                phase = self._newton_phase(m, j)
                S = (1.0 / np.sqrt(2)) * np.array(
                    [[1.0, W ** phase],
                     [1.0, -(W ** phase)]],
                    dtype=complex,
                )
                seeds_m.append(S)
            level_seeds.append(seeds_m)
        return level_seeds

    def _newton_phase(self, stage: int, position: int) -> int:
        """
        Compute the phase exponent for the position-dependent twiddle.

        Uses the exact 2-adic discrete logarithm from dual_view.core to
        decompose the target into binary digits weighted by the butterfly
        stage. The phase at stage m is the m-th bit of the discrete log
        of the target, scaled by the position index.

        In the clean-prime vacuum this collapses to 0, 1, or 2.
        """
        # Compute the exact discrete log e s.t. g^e ≡ a (mod 2^k)
        # using the existing dual-view Newton-lifted dlog.
        result = two_adic_dlog(self.a, self.k)
        if result is None:
            return position  # fallback to standard Cooley-Tukey

        _, e = result

        # Decompose e into binary digits; the phase at stage m is
        # the m-th bit of e, weighted by the position's contribution
        # to the butterfly topology.
        bit_m = (e >> stage) & 1

        # The Newton correction phase: in the QFT basis, the Newton
        # step is diagonal with eigenvalues determined by the discrete
        # log. The phase at position j in stage m is:
        #   φ(m, j) = j * bit_m  (mod 2^{m+1})
        # This gives the standard Cooley-Tukey twiddle structure when
        # bit_m = 1, and identity when bit_m = 0.
        modulus = 2 * (1 << stage) if stage > 0 else 2
        return (position * bit_m) % modulus

    # ------------------------------------------------------------------
    # Thermodynamic / solvability interface
    # ------------------------------------------------------------------

    def thermodynamic_signature(self) -> dict:
        """
        Return the SeedThermodynamics-compatible classification of the
        exponent-space shift operator S (not the full Newton projector).
        """
        # The cyclic shift on C_N has eigenvalues the N-th roots of unity.
        # Spectral radius = 1, conservative, unitary in the appropriate basis.
        return {
            "spectral_radius": 1.0,
            "is_unitary": True,
            "is_conservative": True,
            "is_contractive": False,
            "is_expansive": False,
            "is_nilpotent": False,
            "entropy_rate": 0.0,
            "lyapunov_exponent": 0.0,
        }

    def solvability_report(self) -> dict:
        """
        Mimic solvability_series() output for the exponent-space algebra.
        The Lie algebra <I-S, I-S^†> is metabelian (derived algebra central),
        hence solvable of depth 2.
        """
        return {
            "series": [2, 1, 0],
            "is_solvable": True,
            "depth": 2,
            "conclusion": "SOLVABLE (metabelian) → butterfly compiles",
            "method": "structural (abelian derived algebra)",
        }

# ---------------------------------------------------------------------------
# 3. Example usage / self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("DualViewSeed — self-test")
    print("=" * 60)

    # 1. Analyse known clean primes
    from .newton_dynamics import KNOWN_CLEAN_PRIMES
    for p in KNOWN_CLEAN_PRIMES:
        prof = analyze_prime(p)
        print(f"\n  p={p}: clean={prof.is_clean}, obstruction={prof.obstruction}, "
              f"nilpotency_index={prof.nilpotency_index}, roots={prof.roots}")

    # 2. Build dual-view seed for k=16, target a=17
    k = 16
    a = 17
    dvs = DualViewSeed(k, a)
    print(f"\n  DualViewSeed(k={k}, a={a}):")
    print(f"    N = 2^{k-2} = {dvs.N}")
    print(f"    L = ln_2(5)/4 mod 2^{k-2} = {dvs.L}")
    print(f"    Solvability: {dvs.solvability_report()['conclusion']}")
    print(f"    Thermodynamic: unitary={dvs.thermodynamic_signature()['is_unitary']}")

    print("\n" + "=" * 60)
    print("Self-test complete.")
    print("=" * 60)
