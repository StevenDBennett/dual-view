# Research Status

*Research session completed 2026-05-11. Ported to dual-view 2026-05-16.*

---

## Completed

### Theorems Proven

| Theorem | Result | Proof |
|---------|--------|-------|
| $Φ_n^*(0)$ | $2^{\mu_3(n)/6}$ | $A_d(0) = 2^{(3^{d-1}-1)/2}$, Möbius inversion |
| $Φ_n^*(1)$ | $3^{\mu_3(n)/2}$ | $B_d(1) = 3^{(3^d-1)/2}$, L'Hôpital |
| $\prod\mu$ | $6^{\mu_3(n)/6}$ | Ratio of $Φ_n^*(1)/Φ_n^*(0)$ |
| $\mu_u = \mu_x$ | Period-4 identity | Cycle equation $\prod(2u_i+1) = 81\prod u_i$ |

### Computations

- **Period-5 dynatomic**: Degree 80 in $u$, irreducible, discriminant 5,605 digits
- **Period-6 dynatomic**: Degree 232 in $u$, coefficients computed (233 entries)
- **Clean primes**: Confirmed {7, 103, 181} below 100,000
- **Galois group (period-4)**: Subgroup of $S_6 \wr V_4$

| Period | deg(u) | $\Phi(0)$ | $\Phi(1)$ | $\prod\mu$ | $\sum\mu$ | Cycles |
|--------|--------|-----------|-----------|------------|-----------|--------|
| 2 | 2 | $2^1$ | $3^3$ | $6^1$ | 6 | 1 |
| 4 | 24 | $2^{12}$ | $3^{36}$ | $6^{12}$ | 90 | 6 |
| 5 | 80 | $2^{40}$ | $3^{120}$ | $6^{40}$ | 486 | 16 |
| 6 | 232 | $2^{116}$ | $3^{348}$ | $6^{116}$ | — | — |

## In Progress

### Galois Group Identification
- **Evidence**: Subgroup of $S_6 \wr V_4$, order 2,949,120
- **Next**: Check parity of block permutations; feed constraints into GAP/Magma for exact identification

### Period-6 Dynatomic
- **Status**: Coefficients computed (degree 232 in $u$)
- **Strategy**: Evaluation-interpolation to avoid coefficient explosion
- **Next**: Multiplier computation from coefficients

## Open Problems

### Critical
1. **Identify exact period-4 Galois group** — feed constraints into GAP/Magma
2. **Prove clean primes = {7, 103, 181}** — use obstruction accumulation
3. **Prove $Φ_n^*(1) = 3^{\mu_3(n)/2}$** — full theoretical proof (currently has complete proof — see `theorems.md`)

### High
4. **Compute period-6 multipliers** — numerical computation from coefficients
5. **Find closed form for individual multipliers** — not just product

### Medium
6. **Relate $K_2$ to ray class fields of $\mathbb{Q}(\zeta_3)$**
7. **Generalize to degree $d$ maps**
