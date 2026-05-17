"""
butterfly.py
------------
Kronecker/operadic factor cliff scoring.

For Kronecker-factored weight matrices (e.g. in butterfly or low-rank
factorisations), each factor can be analysed independently for its
2-adic ghost cliff.  The composition cliff of the Kronecker product
is the minimum of the factor cliffs (since precision lost in any
single factor propagates to the product).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from .scaling import scale_weights, auto_scale
from .thermodynamics import SeedThermodynamics


class KroneckerCliffScorer:
    """
    Score each Kronecker factor independently for ghost cliffs.

    Parameters
    ----------
    factors : list of np.ndarray
        Float weight matrices (one per Kronecker factor).
    k_range : range
        Precision range for SeedThermodynamics analysis.
    scale_mode : str
        Scaling mode for scale_weights ('round', 'floor', 'stochastic').
    ensure_odd : bool
        Whether to nudge even integers to odd.
    """

    def __init__(
        self,
        factors: List[np.ndarray],
        k_range: range = range(4, 13),
        scale_mode: str = "round",
        ensure_odd: bool = True,
    ) -> None:
        self.factors = factors
        self.k_range = k_range
        self.scale_mode = scale_mode
        self.ensure_odd = ensure_odd
        self._results: Optional[Dict] = None

    def score_factors(self) -> Dict:
        """Run SeedThermodynamics on each factor and collect results."""
        results = {}
        for idx, F in enumerate(self.factors):
            s = auto_scale(F)
            F_int, meta = scale_weights(F, s, self.scale_mode, self.ensure_odd)
            st = SeedThermodynamics(k=max(self.k_range), g=5)
            st(F_int, self.k_range)
            st.compute()
            summary = st.summary()
            results[f"factor_{idx}"] = {
                "shape": F.shape,
                "scale": s,
                "summary": summary,
                "meta": meta,
            }
        self._results = results
        return results

    def composition_cliff(self) -> Optional[float]:
        """
        Lower bound on the precision at which the Kronecker product
        breaks: min of all factor cliff minima (or None if no data).
        """
        if self._results is None:
            self.score_factors()
        cliffs = []
        for res in self._results.values():
            c = res["summary"].get("min_cliff")
            if c is not None and c > 0:
                cliffs.append(c)
        return min(cliffs) if cliffs else None

    def fragile_factors(self, threshold: float = 6.0) -> List[int]:
        """
        Indices of factors with mean_cliff below threshold.
        """
        if self._results is None:
            self.score_factors()
        fragile = []
        for key, res in self._results.items():
            idx = int(key.split("_")[1])
            mean_c = res["summary"].get("mean_cliff", 0.0)
            if mean_c < threshold:
                fragile.append(idx)
        return fragile

    def print_report(self) -> str:
        """Formatted report of all factor scores."""
        if self._results is None:
            self.score_factors()
        lines = ["Kronecker Factor Cliff Report", ""]
        for key, res in self._results.items():
            s = res["summary"]
            flag = " [!]" if s.get("ghost_fraction", 0) > 0.5 else " [*]"
            lines.append(
                f"  {key}: shape={res['shape']}, "
                f"scale={res['scale']:.1f}, "
                f"ghost={s.get('ghost_fraction', 0):.0%}, "
                f"cliff={s.get('mean_cliff', 0):.1f}{flag}"
            )
        comp = self.composition_cliff()
        lines.append(f"\n  Composition cliff: {comp}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"KroneckerCliffScorer({len(self.factors)} factors)"


def semiring_cliff_score(
    factor_cliffs: List[Optional[float]],
    semiring: str = "standard",
) -> Optional[float]:
    """
    Aggregate factor cliff scores under different algebraic regimes.

    Parameters
    ----------
    factor_cliffs : list of float or None
        Per-factor minimum cliffs.
    semiring : str
        'standard' (min, for standard arithmetic),
        'tropical' (max, for tropical semiring),
        'boolean' (min, same as standard).

    Returns
    -------
    Aggregate cliff, or None if all are None.
    """
    valid = [c for c in factor_cliffs if c is not None]
    if not valid:
        return None
    if semiring == "tropical":
        return float(max(valid))
    else:
        return float(min(valid))
