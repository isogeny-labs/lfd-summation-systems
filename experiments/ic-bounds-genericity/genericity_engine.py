"""Group-algebra closure engine for the genericity study (N6).

Model: A_k^{(c)} = F_p[x1, x2] / (x1^{2^k} - c, x2^{2^k} - c).  Monomials are exponent pairs
in [0, 2^k)^2, the degree is the Hamming weight, and the chain variable u_{i,j} acts as the
translation by 2^j in block i; when the addition overflows 2^k the monomial wraps and picks
up the factor c.  c = 1 is the ordinary chain algebra A_k (subgroup H of order 2^k), c = 0 is
the truncated polynomial ring (no wrap), and a general c is the coset mu^{-1} H with
c = mu^{-2^k}: by the weight-preserving linear map phi_mu (x^e -> mu^{-|e|} x^e) the closure
of g(mu^2 a, mu^3 b, mu t) in A_k equals phi_mu of the closure of g(a, b, t) in A_k^{(c)}
(docs/proofs/N6-genericity.md, Proposition 3.2).

The closure V_D is the least subspace of M_{<=D} containing the generator rows tau_mu g
(wt mu <= D - 2) and closed under tau_{2^j e_i} applied to its elements of weight < D.
`closure` is the once-per-pivot shift engine of csd/ic_bounds_fall.py with the twisted
shift and a final verification pass; it returns the weight-sorted RREF basis.
"""
from __future__ import annotations

import numpy as np

from csd import elim
from csd.ic_bounds_fall import Grid, make_grid, rref


def s3_pairs(a: int, b: int, t: int, p: int) -> dict[tuple[int, int], int]:
    """Coefficients of S_3(x1, x2, t) for y^2 = x^3 + a x + b (textbook formula)."""
    c = {(2, 0): t * t, (0, 2): t * t, (1, 1): -2 * t * t - 2 * a, (2, 1): -2 * t, (1, 2): -2 * t,
         (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b, (0, 0): -4 * b * t + a * a, (2, 2): 1}
    return {e: v % p for e, v in c.items() if v % p}


def generator(a: int, b: int, t: int, p: int, k: int) -> dict[int, int]:
    N = 1 << k
    g: dict[int, int] = {}
    for (e1, e2), v in s3_pairs(a, b, t, p).items():
        m = (e1 % N) | ((e2 % N) << k)
        g[m] = (g.get(m, 0) + v) % p
    return {m: v for m, v in g.items() if v}


class TwistGrid:
    """Grid plus, for every variable, the wrap flag of every column."""

    def __init__(self, k: int, D: int):
        self.grid: Grid = make_grid(k, D)
        self.k, self.D, self.N = k, D, 1 << k
        self.wraps = []
        for v in range(2 * k):
            i, j = divmod(v, k)
            flag = np.zeros(len(self.grid.masks), dtype=bool)
            for col, mask in enumerate(self.grid.masks):
                a = [mask & (self.N - 1), mask >> k]
                flag[col] = a[i] + (1 << j) >= self.N
            self.wraps.append(flag)

    @property
    def masks(self):
        return self.grid.masks

    @property
    def weights(self):
        return self.grid.weights

    @property
    def const_column(self) -> int:
        return len(self.grid.masks) - 1


def generator_rows(g: dict[int, int], tg: TwistGrid, p: int) -> np.ndarray:
    """tau_mu g for wt mu <= D - deg g.  These never wrap (exponents <= 2^{k-1} + 2^{k-2} + 2)."""
    grid, N, k, D = tg.grid, tg.N, tg.k, tg.D
    deg = max(m.bit_count() for m in g)
    rows = []
    for mu in grid.masks:
        if mu.bit_count() + deg > D:
            continue
        row = np.zeros(len(grid.masks), dtype=np.int64)
        m1, m2 = mu & (N - 1), mu >> k
        for mask, c in g.items():
            e1, e2 = m1 + (mask & (N - 1)), m2 + (mask >> k)
            assert e1 < N and e2 < N
            col = grid.column(e1, e2)
            row[col] = (row[col] + c) % p
        rows.append(row)
    return np.asarray(rows, dtype=np.int64)


def twisted_shift(rows: np.ndarray, tg: TwistGrid, v: int, c: int, p: int) -> np.ndarray:
    perm = tg.grid.perms[v]
    assert np.all(rows[:, perm < 0] == 0)
    keep = perm >= 0
    factor = np.where(tg.wraps[v], c % p, 1)[keep]
    out = np.zeros_like(rows)
    out[:, perm[keep]] = rows[:, keep] * factor % p
    return out


def closure(g: dict[int, int], tg: TwistGrid, p: int, c: int = 1, *, verify: bool = True,
            product_cap: int = 20000) -> dict:
    """V_D in A_k^{(c)}: weight-sorted RREF basis, pivots, rank history, closure flag."""
    k, D = tg.k, tg.D
    basis = generator_rows(g, tg, p)
    processed: set[int] = set()
    history, rounds = [], 0
    budget = max(1, product_cap // (2 * k))
    while True:
        out = []
        rank, pivots, _ = elim.fp_echelon(basis, p, out)
        basis = out[0]
        history.append(int(rank))
        if tg.const_column in pivots:
            break
        pending = [i for i, col in enumerate(pivots) if tg.weights[col] < D and col not in processed]
        if not pending:
            break
        chosen = pending[:budget]
        processed.update(int(pivots[i]) for i in chosen)
        basis = np.vstack([basis] + [twisted_shift(basis[chosen], tg, v, c, p) for v in range(2 * k)])
        rounds += 1
    closed = None
    if verify and tg.const_column not in pivots:
        low = [i for i, col in enumerate(pivots) if tg.weights[col] < D]
        closed = True
        for start in range(0, len(low), budget):
            chunk = low[start:start + budget]
            stacked = np.vstack([basis] + [twisted_shift(basis[chunk], tg, v, c, p) for v in range(2 * k)])
            if elim.fp_echelon(stacked, p)[0] != rank:
                closed = False
                break
    basis = rref(basis, list(pivots), p)
    pivots = [int(x) for x in pivots]
    return dict(basis=basis, pivots=pivots, rank=int(rank), history=history, rounds=rounds, closed=closed,
                contains_one=tg.const_column in pivots,
                profile=[int(sum(1 for col in pivots if tg.weights[col] <= d)) for d in range(D + 1)])


def row_space_relation(A: np.ndarray, B: np.ndarray, p: int) -> dict:
    ra = elim.fp_echelon(A, p)[0] if A.shape[0] else 0
    rb = elim.fp_echelon(B, p)[0] if B.shape[0] else 0
    rab = elim.fp_echelon(np.vstack([A, B]), p)[0]
    return dict(rank_a=int(ra), rank_b=int(rb), rank_sum=int(rab), equal=bool(ra == rb == rab),
                a_in_b=bool(rab == rb), b_in_a=bool(rab == ra))


def subgroup_generator(p: int, k: int) -> int:
    """An element of exact order 2^k in F_p^* (2^k | p - 1)."""
    N = 1 << k
    assert (p - 1) % N == 0
    for x in range(2, p):
        z = pow(x, (p - 1) // N, p)
        if pow(z, N // 2, p) != 1:
            return z
    raise ValueError


def coset_representatives(p: int, k: int) -> list[int]:
    """One representative of every coset of H = {x : x^{2^k} = 1} in F_p^*: the powers of a
    primitive root below (p-1)/2^k."""
    N = 1 << k
    m = (p - 1) // N
    # primitive root
    phi = p - 1
    fac = []
    n = phi
    d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        fac.append(n)
    for r in range(2, p):
        if all(pow(r, phi // q, p) != 1 for q in fac):
            break
    return [pow(r, i, p) for i in range(m)]
