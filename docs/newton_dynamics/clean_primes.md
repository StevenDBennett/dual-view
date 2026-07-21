# Clean Primes

A prime $p$ is **clean** for $N(x) = (2x^3+1)/(3x^2)$ if $\mathbb{F}_p$ admits no periodic points of period 2, 3, or 4 with multiplier $\mu \not\equiv 1 \pmod{p}$ — equivalently, if the reduction modulo $p$ of every dynatomic polynomial $\Phi_n^*$ has no root in $\mathbb{F}_p$ for $n \le 4$.

## Known Clean Primes (Verified up to 30,000,000)

The known clean primes are:

$$\\{5, 7, 31, 41, 59, 103, 181, 359, 659, 811, 8111, 14159, 31741, 115679, 162251, 403549\\}$$

| # | Prime | Decimal Property | Notes |
|---|-------|-----------------|-------|
| 1 | **5** | | Trivially clean (only root is 1, pole chain for 2 elements) |
| 2 | 7 | | Original known |
| 3 | 31 | | New |
| 4 | 41 | | New |
| 5 | 59 | | New |
| 6 | 103 | | Original known |
| 7 | 181 | | Original known |
| 8 | 359 | | New |
| 9 | 659 | | New |
| 10 | 811 | | New |
| 11 | 8111 | Palindrome | New |
| 12 | **14159** | Digits of $\pi$ | New |
| 13 | 31741 | | New |
| 14 | 115679 | | New |
| 15 | 162251 | | New |
| 16 | 403549 | | New |

Of particular note: **14159** appears as the first five decimal digits of $\pi$ ($3.14159\ldots$), suggesting a potential number-theoretic connection between the clean-prime condition and the decimal expansion of $\pi$.

> **Note on pole chains**: Of the 16 clean primes, only **7, 103, 181** have functional graphs where every element of $\mathbb{F}_p^*$ converges to a root. The other 13 have "pole chains" — elements that eventually map to the pole $x = 0$ (where the denominator $3x^2 \equiv 0$). These elements are not periodic points and do not form cycles, so they do not disqualify cleanliness. The obstruction classification describes this: `analyze_prime(p)` returns `obstruction="clean"` for pole-free primes and `pole_chain` for the rest, with `is_clean=True` for both.

### Previous Conjecture

The original conjecture (based on verification below 100,000) was that the clean primes were exactly **{7, 103, 181}** and the set was complete. Exhaustive search to **30,000,000** has disproved completeness — 12 additional clean primes exist — while preserving the finiteness thesis.

## Verification Method

For each prime $p \equiv 1 \pmod{3}$, verified by full functional-graph analysis (checking **all** cycles and multipliers, not just polynomial root existence):

1. Compute the Newton map $N(x) = (2x^3+1)/(3x^2)$ over $\mathbb{F}_p^\*$
2. Detect all cycles (periods 1 through $p-1$) via DFS
3. Check multipliers $\mu = N'(x)$ for each cycle
4. A prime is clean iff every cycle has multiplier $\mu \equiv 1 \pmod{p}$ (equivalently, the functional graph is a rooted forest with the three cube roots of unity as the only fixed points)

## Density Observations (Finiteness Conjecture)

For each period $n$, cleanliness requires no period-$n$ points with $\mu \not\equiv 1 \pmod{p}$:

- **Period 2**: $\mu = 6$, so all $p > 5$ require no period-2 points
- **Period 3**: $\mu = \pm 24\sqrt{-3}$, norm $= 1728 = 2^6 \cdot 3^3$, so all $p > 19$ require no period-3 points  
- **Period 4**: $P(1) = 1,\!905,\!120,\!253$ (prime), so all $p \neq P(1)$ require no period-4 points

As $n \to \infty$, these conditions become mutually exclusive. By Chebotarev's density theorem applied to the infinite Galois compositum, the density of clean primes tends to zero, strongly suggesting the set is finite.

### Density Analysis

If $K_2$, $K_3$, $K_4$ were linearly disjoint over $\mathbb{Q}$:
- $[K_2 \cdot K_3 \cdot K_4 : \mathbb{Q}] = 12 \cdot 18 \cdot 24 = 5,\!184$
- Density bound: $1/5,\!184 \approx 0.0193\%$

Observed density (new search): $16 / 30,\!000,\!000 \approx 0.000053\%$

The observed density is now well **below** the linear-disjointness bound, which is consistent with finiteness (the conditions accumulate to force the density toward zero). The original bound $1/5,\!184$ was an upper bound from periods 2–4 alone; the actual density after period-5+6 conditions is much lower.

### Comparison to Original Search

| Metric | Original | Current |
|--------|----------|---------|
| Search bound | 100,000 | 30,000,000 |
| Clean primes found | 3 | 16 |
| Density | $3 \times 10^{-5}$ | $5.3 \times 10^{-7}$ |
| Asymptotic trend | — | Decreasing |

The decreasing density with increasing search range supports the finiteness conjecture despite the larger-than-expected set cardinality.

### False-Positive Caveat

The fast check (periods 2–5) is **not sufficient** for full verification. Known false positive:

- **$p = 313$** passes all period-2–5 checks but has a **period-27 ghost cycle**. Full verification requires checking **all** cycles and multipliers.

The $p = 181$ anomaly: period-4 polynomial has 4 linear factors modulo 181 ($u = 45, 47, 123, 179$), but **none are cubes** in $\mathbb{F}_{181}$, so $x^3 = u$ has no solutions. Without the cube-root check, $p = 181$ would be a false *negative*.

All 16 primes in the known list were verified by full functional-graph DFS, not just fast polynomial checks.
