"""
butterfly_emitter.py — QASM emitter with basin-routing annotations
===================================================================

Emits OpenQASM 2.0 circuits annotated with classical basin routing
information for clean primes.  Uses D = ceil(log2(M)) exponent qubits
to reduce routing-table size, but the QFT structure is unchanged.

A true quantum depth reduction exploiting nilpotency (3-way multiplexer,
swap-routing on qubits) remains an open problem — see research/.
"""

from __future__ import annotations

import math
import numpy as np
from typing import Optional

from .butterfly_seed import (
    analyze_prime,
    CleanPrimeProfile,
    _newton_fp,
    _valuation,
    dual_view_qasm_emitter,
)
from .core import two_adic_log5


_LARGE_BASIN_LIMIT = 5000


def _basin_depth_d(p: int) -> int:
    """Compute the butterfly depth D = ceil(log2(M)) for a clean prime."""
    prof = analyze_prime(p)
    if not prof.is_clean:
        return 0
    M = prof.nilpotency_index
    return max(1, (M - 1).bit_length()) if M > 0 else 0


def _swap_count_by_depth(prof: CleanPrimeProfile, D: int) -> list[dict]:
    """
    Compute swap counts per stage efficiently using depth data.

    For large basins, avoids enumerating every element. Instead counts
    how many elements at each depth need to advance.
    """
    basin = prof.basin_ordering
    roots_set = set(prof.roots)
    depths = prof.tree_depths

    stages = []
    for s in range(D):
        step = 1 << s
        count = 0
        for x in basin:
            if x in roots_set:
                continue
            d = depths.get(x, -1)
            if d >= 0 and d > 0:  # has room to advance
                count += 1
        stages.append({
            "stage": s,
            "step": step,
            "swaps": [],  # no individual swap data for large basins
            "total_swaps": count,
        })
    return stages


def _basin_swap_network(prof: CleanPrimeProfile, D: int) -> list[dict]:
    """
    Build the multi-stage swap network from the basin forest.

    For basins <= _LARGE_BASIN_LIMIT, enumerates individual swaps.
    For larger basins, uses depth-based count only.
    """
    basin = prof.basin_ordering
    roots_set = set(prof.roots)

    if len(basin) > _LARGE_BASIN_LIMIT:
        return _swap_count_by_depth(prof, D)

    pos = {x: i for i, x in enumerate(basin)}

    stages = []
    for s in range(D):
        step = 1 << s
        swaps = []
        for x in basin:
            if x in roots_set:
                continue
            cur = x
            for _ in range(step):
                parent = _newton_fp(cur, prof.p)
                if parent is None or parent in roots_set:
                    cur = parent
                    break
                cur = parent
            if cur is not None and cur != x and cur in pos:
                i, j = pos[x], pos[cur]
                if i != j:
                    swaps.append((i, j))
        deduped = set((min(i, j), max(i, j)) for i, j in swaps)
        deduped_list = sorted(deduped)[:20]
        stages.append({
            "stage": s,
            "step": step,
            "swaps": deduped_list,
            "total_swaps": len(deduped),
        })
    return stages


def basin_qasm_emitter(p: int, k: int, target: int = 5,
                       _prof: Optional[CleanPrimeProfile] = None) -> str:
    """
    Generate OpenQASM 2.0 annotated with basin routing information.

    The circuit uses D = ceil(log2(M)) exponent qubits and a D-stage
    QFT instead of the full (k-2)-stage QFT — this reduces the classical
    routing-table size but does NOT change the quantum gate count
    (the QFT structure is identical, just smaller).

    The basin routing swap network is documented as QASM comments.
    A true quantum depth reduction using nilpotency (3-way multiplexer,
    actual swap-routing on qubits) remains an open problem.

    Falls back to the standard emitter if p is not clean.
    """
    prof = _prof if _prof is not None else analyze_prime(p)
    if not prof.is_clean:
        return dual_view_qasm_emitter(k, target, p_clean=None)

    roots = prof.roots
    M = prof.nilpotency_index
    D = max(1, (M - 1).bit_length()) if M > 0 else 1
    n_basin = len(prof.basin_ordering)
    basin = prof.basin_ordering

    # Register sizes
    n_exp = D
    n_val = max(1, (k + 1) // 3)
    n_sign = 1
    n_flag = 1
    n_total = n_exp + n_val + n_sign

    # Classical precomputation for state prep
    mod = 2**k
    a = target % mod
    if a % 4 != 1:
        a = (-a) % mod
    from .core import _dlog_newton
    e_classical = _dlog_newton(a, k)

    # Build swap network documentation
    stages = _basin_swap_network(prof, D)

    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"// Basin-annotated Dual-View  (p={p}, k={k}, a={target})",
        f"// Nilpotency M={M},  Butterfly depth D={D},  Basin size={n_basin}",
        f"// Roots: {roots}",
        f"// Exponent qubits: {n_exp} (routing table size, not quantum advantage)",
        "",
        f"qreg q[{n_total}];",
        f"qreg flag[{n_flag}];",
        f"creg c[{n_total}];",
        f"creg cflag[{n_flag}];",
        "",
        "// --- 1. State preparation (Hensel bootstrap) ---",
    ]

    # State prep: encode basin depth from classical exponent
    depth_bits = []
    for i in range(n_exp):
        bit = (e_classical >> i) & 1
        depth_bits.append(bit)
        if bit:
            lines.append(f"x q[{i}];  // e[{i}] = 1")
    lines.append("")

    # --- 2. Reduced QFT on D qubits ---
    lines.append(f"// --- Reduced QFT on {D} qubits (basin-ordered) ---")
    for m in range(n_exp):
        lines.append(f"h q[{m}];")
        for j in range(m):
            denom = 1 << (m - j)
            lines.append(f"cp(pi/{denom}) q[{j}], q[{m}];")
    lines.append("")

    # --- 3. Neumann-series phase accumulation ---
    lines.append("// --- Newton diagonal (Neumann series finite sum) ---")
    L = two_adic_log5(k)
    for m in range(n_exp):
        # Phase computed from basin-path length, truncated to M terms
        phase = (target * L) % mod
        # Scale phase by M / 2^k to match the reduced register
        phase_reduced = (phase * M) // max(1, mod)
        angle = 2 * np.pi * phase_reduced / max(1, 2**n_exp)
        lines.append(f"p({angle:.10f}) q[{m}];  // Neumann phase stage {m} (M={M})")
    lines.append("")

    # --- 4. Reduced inverse QFT ---
    lines.append(f"// --- Inverse QFT on {D} qubits ---")
    for m in range(n_exp - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            denom = 1 << (m - j)
            lines.append(f"cp(-pi/{denom}) q[{j}], q[{m}];")
        lines.append(f"h q[{m}];")
    lines.append("")

    # --- 5. Basin routing network (documented as comments) ---
    lines.append("// --- Basin routing stages (conceptual element swaps) ---")
    lines.append(f"//  {len(stages)} stages,  {sum(s['total_swaps'] for s in stages)} total swaps")
    lines.append("//  (element indices in F_p^*, not qubit indices)")
    for s in stages[:6]:  # cap display at 6 stages
        lines.append(f"//  Stage {s['stage']}: advance by {s['step']} step(s)  "
                      f"({s['total_swaps']} swaps)")
        for i, j in s['swaps'][:6]:
            i_elem = basin[i] if i < len(basin) else i
            j_elem = basin[j] if j < len(basin) else j
            lines.append(f"//    basin[{i}]={i_elem} <-> basin[{j}]={j_elem}")
        if s['total_swaps'] > 6:
            lines.append(f"//    ... and {s['total_swaps'] - 6} more")
    if len(stages) > 6:
        lines.append(f"//  ... and {len(stages) - 6} more stages (structural)")
    lines.append("")

    # --- 6. Valuation guard ---
    lines.append("// --- Valuation guard ---")
    val_start = n_exp
    g_val = 5
    c_cliff = max(0, _valuation(g_val + 123) - 2)
    threshold = c_cliff
    lines.append(f"//   Valuation register: q[{val_start}..{val_start + n_val - 1}]")
    lines.append(f"//   Cliff constant: c(g) = max(0, v2(g+123)-2) = {threshold}")
    for i in range(n_val):
        lines.append(f"measure q[{val_start + i}] -> c[{val_start + i}];")
    lines.append(f"if (c == 0) x flag[0];")
    lines.append("barrier q, flag;")
    lines.append("")

    # --- 7. Measurement ---
    lines.append("// --- Measurement ---")
    for i in range(n_total):
        lines.append(f"measure q[{i}] -> c[{i}];")
    lines.append(f"measure flag[0] -> cflag[0];")

    return "\n".join(lines)


def dual_view_qasm_emitter_clean(k: int, target_a: int,
                                  p_clean: Optional[int] = None) -> str:
    """
    Wrapper: delegates to the basin-annotated emitter for clean primes,
    falls back to the standard emitter otherwise.

    This is the public entry point.  Import as:
        from dual_view import dual_view_qasm_emitter_clean
    """
    if p_clean is not None:
        prof = analyze_prime(p_clean)
        if prof.is_clean:
            return basin_qasm_emitter(p_clean, k, target_a, _prof=prof)
    return dual_view_qasm_emitter(k, target_a, p_clean=None)


# ---------------------------------------------------------------------------
# Self-test / demonstration
# ---------------------------------------------------------------------------

def print_comparison(p: int, k: int = 16):
    """Compare standard vs optimised circuit sizes for a clean prime."""
    import sys
    prof = analyze_prime(p)
    M = prof.nilpotency_index
    D = max(1, (M - 1).bit_length()) if M > 0 else 1
    n_std = k - 2
    std_gates = n_std * (n_std + 1) // 2  # QFT cp pairs
    opt_gates = D * (D + 1) // 2

    print(f"  p={p:>6}  M={M:>5}  D={D:>3}  "
          f"std={n_std}q {std_gates:>4}cp  "
          f"opt={D}q {opt_gates:>4}cp  "
          f"ratio={std_gates/max(opt_gates,1):.1f}x")


if __name__ == "__main__":
    from .newton_dynamics import KNOWN_CLEAN_PRIMES

    print("=" * 70)
    print("Butterfly Emitter — Depth comparison (k=16)")
    print("=" * 70)
    print(f"  {'Prime':>6}  {'M':>5}  {'D':>3}  Standard          Optimised         Ratio")
    print("  " + "-" * 58)
    for p in KNOWN_CLEAN_PRIMES:
        print_comparison(p, k=16)

    print("\n" + "=" * 70)
    demo = basin_qasm_emitter(7, 8, 5)
    print(f"Demo circuit (p=7, k=8):\n{demo}")
