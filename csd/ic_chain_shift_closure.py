"""Shift-structured closure engine for two-summand chain systems at depth k >= 8.

The chain quotient F_p[u]/(chains) is the group algebra of G = (Z/2^k)^2:
a squarefree chain monomial is the pair of block exponents (a1, a2) and
multiplication by a chain variable u_{i,j} adds 2^j to a_i modulo 2^k, so
the product of a row by a variable is a permutation of its coefficients and
is never expanded through a Macaulay product matrix.  The ordinary degree of
a monomial is the Hamming weight of (a1, a2).

For a degree D the engine establishes

* an exact lower profile through degree D0 < D with the dense engine
  (`ic_fall_certificates.chain_degree_closure`, exclusion-certified), which
  fixes the strict falls through D0 and the dimension of V_D0;
* a *sound* upper certificate at D: rows that lie in V_D by construction
  (normalized multiples of the generator, and shifts of rows of weight < D of
  the running row space), eliminated over F_p; if their rank reaches
  dim(I ∩ M_{<=D}) = |M_{<=D}| - |Z| and the evaluation of M_{<D} is onto the
  solution set Z, then V_D = I ∩ M_{<=D} and the full ideal has LFD <= D
  (local-to-global theorem).  Every low pivot row is processed once, so when
  no low pivot is left the row space is closed under all admissible products
  and is exactly V_D (confirmed by a verification pass); the exact dimension
  then decides the strict fall at D and serves as the reference for D+1.

Exact value: with the certificate at D and the exact profile through D0 = D-1,
V_D ∩ M_{<D} = I ∩ M_{<D} has dimension |M_{<D}| - |Z|; if dim V_{D-1} is
smaller there is a strict fall at D and the exact LFD is D, otherwise the exact
LFD is the last fall of the dense profile.  Unit ideals (Z empty) use the same
rank target |M_{<=D}| and the dense profile's absence of 1 through D-1.
"""
from __future__ import annotations

import math
import random
from itertools import combinations

import numpy as np

from . import elim
from . import ic_fall_certificates as C
from .ic_curve_sweep import prime_evaluation_rank


def monomials(k: int, D: int) -> list[int]:
    r = 2 * k
    masks = [sum(1 << j for j in choice) for d in range(min(r, D) + 1) for choice in combinations(range(r), d)]
    return sorted(masks, key=lambda mask: (mask.bit_count(), mask), reverse=True)


def shift_permutation(k: int, masks: list[int], index: dict[int, int], variable: int) -> np.ndarray:
    """Column map of multiplication by chain variable `variable` (0..2k-1): block i = variable // k,
    depth j = variable % k.  Columns whose image leaves the degree window map to -1."""
    N = 1 << k
    i, j = divmod(variable, k)
    out = np.full(len(masks), -1, dtype=np.int64)
    for col, mask in enumerate(masks):
        a = [mask & (N - 1), mask >> k]
        a[i] = (a[i] + (1 << j)) % N
        new = a[0] | (a[1] << k)
        if new in index:
            out[col] = index[new]
    return out


def generator_rows(generator: dict[int, int], k: int, D: int, masks: list[int], index: dict[int, int], p: int) -> np.ndarray:
    """NF(mu * g) for every squarefree mu with wt(mu) + deg(g) <= D."""
    N = 1 << k
    deg = max(mask.bit_count() for mask in generator)
    rows = []
    for mu in masks:
        if mu.bit_count() + deg > D:
            continue
        row = np.zeros(len(masks), dtype=np.int64)
        m1, m2 = mu & (N - 1), mu >> k
        for mask, c in generator.items():
            a1, a2 = (m1 + (mask & (N - 1))) % N, (m2 + (mask >> k)) % N
            col = index[a1 | (a2 << k)]
            row[col] = (row[col] + c) % p
        rows.append(row)
    return np.asarray(rows, dtype=np.int64).reshape((-1, len(masks)))


def shift_rows(rows: np.ndarray, perm: np.ndarray) -> np.ndarray:
    assert np.all(rows[:, perm < 0] == 0), 'A shifted row left the degree window.'
    out = np.zeros_like(rows)
    keep = perm >= 0
    out[:, perm[keep]] = rows[:, keep]
    return out


def shift_closure(generator: dict[int, int], k: int, p: int, D: int, target_rank: int, *,
                  product_cap: int = 9000, max_cells: int = 200_000_000, verify: bool = True, log=None):
    """V_D modulo the chains by processing every low pivot row once.

    Rows: NF(mu*g), then the 2k shifts of every echelon row whose pivot has
    weight < D, each processed once by pivot column (later tail reductions of
    a processed row are combinations of rows with smaller pivots, which are
    processed in turn).  Every row lies in V_D by construction, so reaching
    `target_rank` is a sound upper certificate.  When every low pivot has been
    processed and `verify` is set, the shifts of all final low rows are
    re-reduced to confirm closure, which makes the final row space exactly V_D.
    Returns a dict with rank, reached, closed, history and low_dimension (= dim V_D ∩ M_{<D} when complete)."""
    masks = monomials(k, D)
    index = {mask: col for col, mask in enumerate(masks)}
    perms = [shift_permutation(k, masks, index, v) for v in range(2 * k)]
    weights = np.asarray([mask.bit_count() for mask in masks])
    basis = generator_rows(generator, k, D, masks, index, p)
    processed = set()
    history = []
    budget = max(1, product_cap // (2 * k))
    rounds = 0
    while True:
        rounds += 1
        if basis.size > max_cells:
            raise ValueError('Shift-closure matrix budget exceeded.')
        out = []
        rank, pivots, _ = elim.fp_echelon(basis, p, out)
        basis = out[0]
        history.append(rank)
        if log:
            log(dict(round=rounds, rank=rank, target=target_rank, processed=len(processed)))
        low_count = sum(1 for col in pivots if weights[col] < D)
        if rank >= target_rank:
            return dict(rank=rank, reached=True, closed=False, history=history, low_dimension=low_count)
        pending = [i for i, col in enumerate(pivots) if weights[col] < D and col not in processed]
        if not pending:
            break
        chosen = pending[:budget]
        processed.update(pivots[i] for i in chosen)
        basis = np.vstack([basis] + [shift_rows(basis[chosen], perm) for perm in perms])
    closed = False
    if verify:
        low = [i for i, col in enumerate(pivots) if weights[col] < D]
        closed = True
        for start in range(0, len(low), budget):
            chunk = low[start:start + budget]
            stacked = np.vstack([basis] + [shift_rows(basis[chunk], perm) for perm in perms])
            if elim.fp_echelon(stacked, p)[0] != rank:
                closed = False
                break
        if log:
            log(dict(verification='closed' if closed else 'not closed', rank=rank))
    return dict(rank=rank, reached=False, closed=closed, history=history, low_dimension=low_count)


def shift_profile(generator: dict[int, int], k: int, p: int, points: list[int], generator_value: int, *,
                  dense_degree: int = 4, certificate_degrees=(5, 6), product_cap: int = 9000,
                  dense_max_columns: int = 16384, dense_max_cells: int = 200_000_000,
                  shift_max_cells: int = 200_000_000, log=None):
    """Exact LFD of {g} + chain (two summands) at large depth; see the module docstring."""
    r = 2 * k
    N = 1 << k
    points = list(points)
    assert len(points) == len(set(points))
    for bits in points:
        E1, E2 = bits & (N - 1), bits >> k
        assert sum(c * pow(generator_value, ((mask & (N - 1)) * E1 + (mask >> k) * E2) % N, p)
                   for mask, c in generator.items()) % p == 0
    count = lambda D: sum(math.comb(r, d) for d in range(min(r, D) + 1))
    # exact dense profile through dense_degree
    spaces = []
    falls = []
    for D in range(dense_degree + 1):
        space = C.chain_degree_closure([generator], 2, k, p, D, max_columns=dense_max_columns, max_cells=dense_max_cells)
        if spaces and space.low_degree_dimension(D - 1) > spaces[-1].rank:
            falls.append(D)
        spaces.append(space)
        if space.contains_one:
            return dict(status='exact-unit', exact_lfd=D, lfd_lower_bound=D, lfd_upper_bound=D, strict_falls=falls,
                        scanned_through=D, dense_degree=D, point_count=0, engine='dense-through-unit')
        assert C.check_closed_exclusion_certificate(space, [generator])
        if points and D > 0:
            rank = prime_evaluation_rank(points, [m for m in space.masks if m.bit_count() < D], k, generator_value, p,
                                         dense_max_cells)
            codim = len(space.masks) - space.quotient_rank
            if rank == len(points) == codim:
                last = max(falls, default=0)
                return dict(status='exact-proper-local', exact_lfd=last, lfd_lower_bound=last, lfd_upper_bound=last,
                            strict_falls=falls, scanned_through=D, local_upper_bound=D, dense_degree=D,
                            point_count=len(points), engine='dense-local')
    prev_dim = spaces[-1].quotient_rank  # exact dim of V_{dense_degree} modulo chains
    prev_degree = dense_degree
    checks = []
    for D in certificate_degrees:
        if D != prev_degree + 1:
            break
        if points:
            onto = prime_evaluation_rank(points, [m for m in monomials(k, D) if m.bit_count() < D], k, generator_value,
                                         p, dense_max_cells) == len(points)
        else:
            onto = True
        target = count(D) - len(points)
        closure = shift_closure(generator, k, p, D, target, product_cap=product_cap, max_cells=shift_max_cells, log=log)
        checks.append(dict(degree=D, target_rank=target, rank=closure['rank'], rounds=len(closure['history']),
                           reached=closure['reached'], closed_when_complete=closure['closed'], evaluation_onto=onto,
                           low_dimension=closure['low_dimension'], previous_dimension=prev_dim,
                           rank_history=closure['history']))
        if closure['reached'] and onto:
            # V_D = I ∩ M_{<=D}, so V_D ∩ M_{<D} = I ∩ M_{<D} has dimension count(D-1) - |Z|
            fall_at_D = prev_dim < count(D - 1) - len(points)
            exact = D if fall_at_D else max(falls, default=0)
            status = 'exact-proper-shift' if points else 'exact-unit-shift'
            return dict(status=status, exact_lfd=exact, lfd_lower_bound=exact, lfd_upper_bound=exact,
                        strict_falls=falls + ([D] if fall_at_D else []), scanned_through=D, local_upper_bound=D,
                        dense_degree=dense_degree, point_count=len(points), certificate_checks=checks,
                        engine='shift-closure')
        if closure['closed']:
            # exact V_D without the certificate: record the fall at D and continue one degree up
            if closure['low_dimension'] > prev_dim:
                falls.append(D)
            prev_dim = closure['rank']
            prev_degree = D
            continue
        break
    # every certificate degree failed; if the closure at the last degree was verified complete, LFD > that degree
    lower = max(falls, default=0)
    return dict(status='undecided-shift', exact_lfd=None, lfd_lower_bound=lower,
                lfd_upper_bound=2 * k + 1, strict_falls=falls, scanned_through=prev_degree, dense_degree=dense_degree,
                point_count=len(points), certificate_checks=checks, engine='shift-closure')
