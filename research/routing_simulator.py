"""
Butterfly Routing Simulator
===========================
Simulates the butterfly routing network for a clean prime p.
Shows elements converging to their roots through ceil(log2(M)) stages.

The routing network implements parallel Newton iteration: each stage
advances elements by 2^k steps in the basin tree (binary exponentiation).
"""
import sys, math, random
sys.path.insert(0, "src")

from dual_view.butterfly_seed import analyze_prime, _newton_fp
from dual_view.newton_dynamics import KNOWN_CLEAN_PRIMES


class RoutingSimulator:
    """
    Simulate butterfly routing for a clean prime.

    Given a starting element x, the routing network advances it
    through ceil(log2(M)) stages, each stage applying N 2^k times.
    """

    def __init__(self, p: int):
        self.p = p
        self.prof = analyze_prime(p)
        self.roots = set(self.prof.roots)
        self.M = self.prof.nilpotency_index
        self.depth = max(1, (self.M - 1).bit_length()) if self.M > 0 else 0

    def _newton_k(self, x: int, k: int) -> int:
        """Apply N k times (Newton iteration)."""
        cur = x
        for _ in range(k):
            cur2 = _newton_fp(cur, self.p)
            if cur2 is None:
                return None
            cur = cur2
        return cur

    def simulate(self, start: int, verbose: bool = True) -> dict:
        """
        Simulate butterfly routing for a starting element.

        At stage s (0-indexed), the element advances by 2^s steps.
        After depth stages, it should reach a root.
        """
        if start % self.p == 0:
            return {"start": start, "error": "pole"}

        path = [start]
        cur = start
        while cur not in self.roots:
            cur = _newton_fp(cur, self.p)
            if cur is None:
                return {"start": start, "error": "pole_path", "path": path}
            path.append(cur)

        root = path[-1]
        path_len = len(path) - 1

        if verbose:
            print(f"\n  Routing path: {' → '.join(str(x) for x in path[:8])}"
                  f"{' … → ' + str(root) if path_len > 7 else ''}")
            print(f"  Path length: {path_len}, Butterfly depth: {self.depth}")

        # Butterfly routing: stage s applies N 2^s times
        state = start
        route = [state]
        for s in range(self.depth):
            step = 1 << s
            state = self._newton_k(state, step)
            if state is None:
                route.append(None)
                if verbose:
                    print(f"  Stage {s}: advance {step} steps → POLE (crashed)")
                break
            route.append(state)
            if verbose:
                arrow = "→" + "─" * min(step * 2, 20) + "→"
                print(f"  Stage {s}: N^{step:<3} {start:>6} {arrow} {state:>6}"
                      f"{'  ✓ ROOT' if state in self.roots else ''}")

        converged = state in self.roots if state is not None else False
        return {
            "start": start,
            "root": root,
            "path": path,
            "path_length": path_len,
            "route": route,
            "butterfly_stages": self.depth,
            "converged": converged,
            "final_state": state,
        }


def run_trials(p: int, n_trials: int = 5):
    """Run random routing trials for a clean prime."""
    rs = RoutingSimulator(p)
    print(f"\n{'='*60}")
    print(f"Routing Simulator — p={p}")
    print(f"  Roots: {rs.prof.roots}")
    print(f"  Nilpotency M: {rs.M}")
    print(f"  Butterfly depth: {rs.depth}")
    print(f"{'='*60}")

    successes = 0
    poles = 0
    for _ in range(n_trials):
        start = random.randint(1, p - 1)
        result = rs.simulate(start, verbose=True)
        if result.get("converged"):
            successes += 1
        elif result.get("error") in ("pole", "pole_path"):
            poles += 1

    print(f"\n  Trials: {n_trials}, Converged: {successes}, Pole crashes: {poles}")
    return successes


if __name__ == "__main__":
    random.seed(42)

    # Demo: specific elements for small primes
    for p in (7, 31, 103):
        rs = RoutingSimulator(p)
        for start in [3, 5, 6]:  # depth-1 elements for p=7
            if start < p:
                rs.simulate(start, verbose=True)

    # Random trials for key primes
    for p in (7, 103, 181, 14159, 403549):
        run_trials(p, n_trials=3)

    print("\nDone.")
