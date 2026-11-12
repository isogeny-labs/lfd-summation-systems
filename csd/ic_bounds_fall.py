"""Fall-space structure of the two-summand prime chain closure (bounds program, item 5).

Setting (see paper/curve-lfd-bounds.tex, Section 7).  Modulo the chains the
two-summand prime presentation lives in the group algebra A = F_p[G],
G = (Z/2^k)^2: a squarefree chain monomial is a pair of block exponents
(a1, a2), its degree is the Hamming weight wt(a1) + wt(a2), and
multiplication by the chain variable u_{i,j} is the translation
tau_{2^j e_i}.  The normalized generator g = NF(S_3(x1, x2, t)) is supported
on {0,1,2}^2 (levels 0 and 1 only) and the ideal is I = gA.

The degree-D closure V_D is the smallest subspace of M_{<=D} containing every
translate tau_mu g with wt(mu) <= D-2 and closed under tau_{2^j e_i} applied to
its elements of weight < D.  This module

* recomputes V_D in the group algebra (the once-per-pivot shift engine of
  csd/ic_chain_shift_closure.py, with the round in which each pivot first
  appears recorded and the final row space verified closed), and returns a
  reduced row echelon basis with columns ordered by weight (heaviest first), so
  that the rows whose pivot has weight <= d span exactly V_D ∩ M_{<=d};
* intersects V_D with monomial subspaces indexed by level sets or by exponent
  patterns (dimension and basis), which is the level-profile decomposition
  asked for in item 5;
* for unit ideals (Z empty) computes the unique multiplier q = v / g of every
  basis element through the 2-adic Fourier transform over F_p, so a closure
  element is described canonically as the combination sum_h q_h tau_h g of
  translates of g (translation weights of supp q are the "which shifts, which
  carries" data);
* tests candidate containing spaces: a subspace W of M_{<=D} (given by a
  spanning set) is checked to contain the generator rows and to be closed under
  every admissible variable product, and 1 is tested for membership.  This is
  the exclusion certificate of ic_fall_certificates in the group algebra; it is
  a sound lower-bound certificate for the instance only (not a proof for all k).

Evidence labels follow the task: everything here is `exact` for the instances
computed, never a bound for other k.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from . import elim
from . import ic_dreg as I
from . import ic_curve_sweep as S


# ----------------------------------------------------------------------------
# monomial grid
# ----------------------------------------------------------------------------


@dataclass
class Grid:
    k: int
    D: int
    masks: list[int]           # squarefree chain masks (block 1 bits 0..k-1, block 2 bits k..2k-1)
    index: dict[int, int]
    weights: np.ndarray
    perms: list[np.ndarray]    # column image of multiplication by variable v (-1 = leaves the window)

    @property
    def N(self):
        return 1 << self.k

    def pair(self, mask: int) -> tuple[int, int]:
        return mask & (self.N - 1), mask >> self.k

    def levels(self, mask: int) -> int:
        """Bit set of levels touched by either block."""
        a1, a2 = self.pair(mask)
        return a1 | a2

    def column(self, a1: int, a2: int) -> int:
        return self.index[(a1 % self.N) | ((a2 % self.N) << self.k)]


def monomials(k: int, D: int) -> list[int]:
    r = 2 * k
    masks = [sum(1 << j for j in c) for d in range(min(r, D) + 1) for c in combinations(range(r), d)]
    return sorted(masks, key=lambda m: (m.bit_count(), m), reverse=True)


def make_grid(k: int, D: int) -> Grid:
    masks = monomials(k, D)
    index = {m: i for i, m in enumerate(masks)}
    N = 1 << k
    weights = np.asarray([m.bit_count() for m in masks], dtype=np.int64)
    perms = []
    for v in range(2 * k):
        i, j = divmod(v, k)
        out = np.full(len(masks), -1, dtype=np.int64)
        for col, mask in enumerate(masks):
            a = [mask & (N - 1), mask >> k]
            a[i] = (a[i] + (1 << j)) % N
            new = a[0] | (a[1] << k)
            if new in index:
                out[col] = index[new]
        perms.append(out)
    return Grid(k, D, masks, index, weights, perms)


# ----------------------------------------------------------------------------
# instances
# ----------------------------------------------------------------------------


def random_instance(p: int, k: int, seed: int, *, min_order_bits: int = 8, max_tries: int = 500) -> dict:
    """Random curve y^2 = x^3 + a x + b over F_p (a, b nonzero, nonsingular, j not in {0, 1728})
    and a random rational target of large order (order > 2^min_order_bits, and never 2).
    Returns the generator g (chain masks at depth k, and as depth-independent exponent pairs), the
    scalar solution set Z at depth k (generator_value_at_k is the depth-k subgroup generator; both
    are depth-specific and must be recomputed for another depth) and the point order."""
    rng = random.Random(seed)
    N = 1 << k
    assert (p - 1) % N == 0
    for _ in range(max_tries):
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        E = I.PrimeCurve(p, a, b)
        try:
            order = I.group_order(E, random.Random(rng.randrange(1 << 30)))
        except ValueError:
            continue
        R = S.random_rational_target(E, rng, p)
        if R[1] == 0:
            continue
        o = I._point_order(E, R, order)
        if o <= (1 << min_order_bits):
            continue
        g0 = S.subgroup_generator_value(p, k)
        Z = S.prime_points(E, k, g0, 1, 1, R[0])
        g = I.pkm_generator(E, R[0], k)
        return dict(prime=p, a=a, b=b, k=k, seed=seed, target=[int(R[0]), int(R[1])], curve_order=int(order),
                    target_order=int(o), generator_value_at_k=int(g0), points_at_k=[int(z) for z in Z],
                    generator={int(m): int(c) for m, c in g.items()},
                    generator_pairs=sorted(((m & (N - 1), m >> k, int(c)) for m, c in g.items())))
    raise RuntimeError('no instance found')


def encode_generator(pairs, k: int) -> dict[int, int]:
    """Chain-mask encoding of g at depth k from exponent pairs (e1, e2, c): block 2 is shifted by k,
    so the mask dictionary must be rebuilt for every depth (the pairs are depth-independent)."""
    N = 1 << k
    g: dict[int, int] = {}
    for e1, e2, c in pairs:
        m = (e1 % N) | ((e2 % N) << k)
        g[m] = g.get(m, 0) + int(c)
    return {m: c for m, c in g.items() if c}


def instance_hash(inst: dict) -> str:
    key = dict(prime=inst['prime'], a=inst['a'], b=inst['b'], k=inst['k'], target=inst['target'],
               generator=sorted((int(m), int(c)) for m, c in inst['generator'].items()))
    return hashlib.sha256(json.dumps(key, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


# ----------------------------------------------------------------------------
# closure
# ----------------------------------------------------------------------------


def generator_rows(g: dict[int, int], grid: Grid, p: int, max_shift_weight: int | None = None) -> tuple[np.ndarray, list[int]]:
    """Rows NF(u^mu g) = tau_mu g for wt(mu) + deg g <= D (admissible generator multiples).
    Returns the rows and the list of mu."""
    N = grid.N
    deg = max(m.bit_count() for m in g)
    D = grid.D
    rows, mus = [], []
    for mu in grid.masks:
        if mu.bit_count() + deg > D:
            continue
        if max_shift_weight is not None and mu.bit_count() > max_shift_weight:
            continue
        row = np.zeros(len(grid.masks), dtype=np.int64)
        m1, m2 = mu & (N - 1), mu >> grid.k
        for mask, c in g.items():
            col = grid.column(m1 + (mask & (N - 1)), m2 + (mask >> grid.k))
            row[col] = (row[col] + c) % p
        rows.append(row)
        mus.append(mu)
    return np.asarray(rows, dtype=np.int64).reshape((-1, len(grid.masks))), mus


def shift_rows(rows: np.ndarray, perm: np.ndarray) -> np.ndarray:
    assert np.all(rows[:, perm < 0] == 0), 'a shifted row left the degree window'
    out = np.zeros_like(rows)
    keep = perm >= 0
    out[:, perm[keep]] = rows[:, keep]
    return out


def rref(basis: np.ndarray, pivots: list[int], p: int) -> np.ndarray:
    """Reduced row echelon form from an echelon basis (pivots normalised to 1)."""
    B = basis.copy() % p
    for r in range(len(pivots) - 1, 0, -1):
        c = pivots[r]
        above = np.flatnonzero(B[:r, c])
        if above.size:
            B[above] = (B[above] - np.outer(B[above, c], B[r])) % p
    return B


def closure(g: dict[int, int], grid: Grid, p: int, *, product_cap: int = 20000, log=None,
            verify: bool = True, reduced: bool = True) -> dict:
    """Exact V_D modulo the chains, with the birth round of every pivot column.

    Round 0: the generator rows.  Round r >= 1: all 2k shifts of every echelon row whose pivot
    has weight < D and has not been processed yet.  Stops when no unprocessed low pivot is left;
    then (verify) the shifts of every final low row are re-reduced to confirm closure, which makes
    the row space exactly V_D (ic_chain_shift_closure docstring).  The result holds the RREF basis
    with pivots, their weights and birth rounds, the rank history and the closure flag."""
    k, D = grid.k, grid.D
    basis, mus = generator_rows(g, grid, p)
    processed: set[int] = set()
    birth: dict[int, int] = {}
    history = []
    rounds = 0
    budget = max(1, product_cap // (2 * k))
    while True:
        out = []
        rank, pivots, _ = elim.fp_echelon(basis, p, out)
        basis = out[0]
        for c in pivots:
            birth.setdefault(int(c), rounds)
        history.append(int(rank))
        if log:
            log(dict(round=rounds, rank=int(rank), processed=len(processed)))
        if len(grid.masks) - 1 in pivots:
            break
        pending = [i for i, c in enumerate(pivots) if grid.weights[c] < D and c not in processed]
        if not pending:
            break
        chosen = pending[:budget]
        processed.update(int(pivots[i]) for i in chosen)
        basis = np.vstack([basis] + [shift_rows(basis[chosen], perm) for perm in grid.perms])
        rounds += 1
    closed = None
    if verify and len(grid.masks) - 1 not in pivots:
        low = [i for i, c in enumerate(pivots) if grid.weights[c] < D]
        closed = True
        for start in range(0, len(low), budget):
            chunk = low[start:start + budget]
            stacked = np.vstack([basis] + [shift_rows(basis[chunk], perm) for perm in grid.perms])
            if elim.fp_echelon(stacked, p)[0] != rank:
                closed = False
                break
    if reduced:
        basis = rref(basis, list(pivots), p)
    pivots = [int(c) for c in pivots]
    return dict(basis=basis, pivots=pivots, rank=int(rank), history=history, rounds=rounds, closed=closed,
                contains_one=(len(grid.masks) - 1) in pivots,
                pivot_weights=[int(grid.weights[c]) for c in pivots],
                birth=[birth[c] for c in pivots], generator_mus=mus)


def low_part(result: dict, grid: Grid, d: int) -> np.ndarray:
    """Basis of V_D ∩ M_{<=d}: the echelon rows whose pivot has weight <= d (columns are sorted by
    weight descending, so such a row lives in M_{<=d}); exact because the basis is echelon."""
    rows = [i for i, c in enumerate(result['pivots']) if grid.weights[c] <= d]
    return result['basis'][rows]


def low_profile(result: dict, grid: Grid) -> list[int]:
    return [int(sum(1 for w in result['pivot_weights'] if w <= d)) for d in range(grid.D + 1)]


# ----------------------------------------------------------------------------
# intersections with monomial subspaces
# ----------------------------------------------------------------------------


def intersect_columns(basis: np.ndarray, keep: np.ndarray, p: int, want_basis: bool = False):
    """dim (row space ∩ span of the columns `keep`) = rank - rank(basis restricted to the other columns).
    With want_basis the rows of an echelon form (other columns first) whose pivot lies in `keep`
    are returned: they span the intersection."""
    other = np.flatnonzero(~keep)
    if basis.shape[0] == 0:
        return (0, basis) if want_basis else 0
    rank_all = basis.shape[0]
    if other.size == 0:
        return (rank_all, basis) if want_basis else rank_all
    sub_rank = elim.fp_echelon(basis[:, other], p)[0] if other.size else 0
    dim = rank_all - sub_rank
    if not want_basis:
        return dim
    order = np.concatenate([other, np.flatnonzero(keep)])
    out = []
    rank, piv, _ = elim.fp_echelon(basis[:, order], p, out)
    E = out[0]
    rows = [i for i, c in enumerate(piv) if c >= other.size]
    back = np.zeros_like(E[rows])
    back[:, order] = E[rows]
    assert len(rows) == dim
    return dim, back


def level_mask_columns(grid: Grid, levels: int) -> np.ndarray:
    """Columns whose monomial touches only the levels in the bit set `levels`."""
    return np.asarray([(grid.levels(m) & ~levels) == 0 for m in grid.masks])


def level_profile(basis: np.ndarray, grid: Grid, p: int, level_sets: list[int]) -> dict[int, int]:
    return {L: intersect_columns(basis, level_mask_columns(grid, L), p) for L in level_sets}


def exact_level_counts(dims: dict[int, int]) -> dict[int, int]:
    """Moebius inversion over the Boolean lattice of level sets present in `dims` (which must be
    downward closed): e(L) = sum_{L' ⊆ L} (-1)^{|L|-|L'|} dim(L').  When the space is a direct sum
    of its intersections with the M_L these are the dimensions of the exact-level-set pieces."""
    out = {}
    for L in dims:
        total = 0
        sub = L
        while True:
            if sub in dims:
                total += (-1) ** ((L ^ sub).bit_count()) * dims[sub]
            if sub == 0:
                break
            sub = (sub - 1) & L
        out[L] = total
    return out


# ----------------------------------------------------------------------------
# Fourier transform over F_p on (Z/2^k)^2 and the multiplier q = v / g (unit ideals)
# ----------------------------------------------------------------------------


def primitive_root_of_unity(p: int, k: int) -> int:
    """A primitive 2^k-th root of unity in F_p (2^k | p-1)."""
    N = 1 << k
    for x in range(2, p):
        z = pow(x, (p - 1) // N, p)
        if pow(z, N // 2, p) != 1:
            return z
    raise ValueError


def _fft1(a: np.ndarray, zeta: int, k: int, p: int) -> np.ndarray:
    """Length-2^k DFT mod p along the last axis: A[j] = sum_i a[i] zeta^{ij}  (iterative radix-2)."""
    N = 1 << k
    a = a.copy() % p
    # bit reversal
    rev = np.zeros(N, dtype=np.int64)
    for i in range(N):
        r = 0
        x = i
        for _ in range(k):
            r = (r << 1) | (x & 1)
            x >>= 1
        rev[i] = r
    a = a[..., rev]
    length = 2
    while length <= N:
        w = pow(zeta, N // length, p)
        half = length // 2
        tw = np.array([pow(w, j, p) for j in range(half)], dtype=np.int64)
        a = a.reshape(a.shape[:-1] + (N // length, length))
        even = a[..., :half]
        odd = (a[..., half:] * tw) % p
        a = np.concatenate([(even + odd) % p, (even - odd) % p], axis=-1)
        a = a.reshape(a.shape[:-2] + (N,))
        length *= 2
    return a


def fft2(a: np.ndarray, zeta: int, k: int, p: int, inverse: bool = False) -> np.ndarray:
    """2-D DFT of an N x N array (index [a1, a2]) over F_p; inverse divides by N^2."""
    N = 1 << k
    z = pow(zeta, -1, p) if inverse else zeta
    b = _fft1(a, z, k, p)
    b = _fft1(b.T, z, k, p).T
    if inverse:
        b = (b * pow(N * N, -1, p)) % p
    return b % p


def to_array(row: np.ndarray, grid: Grid, p: int) -> np.ndarray:
    N = grid.N
    arr = np.zeros((N, N), dtype=np.int64)
    nz = np.flatnonzero(row)
    for col in nz:
        a1, a2 = grid.pair(grid.masks[col])
        arr[a1, a2] = row[col] % p
    return arr


def multipliers(basis: np.ndarray, g: dict[int, int], grid: Grid, p: int) -> np.ndarray:
    """q_v = v / g in F_p[G] for every row v (unit ideal only: g must be invertible, i.e. its
    Fourier transform never vanishes, which is checked).  Returns an array (rows, N, N)."""
    k, N = grid.k, grid.N
    zeta = primitive_root_of_unity(p, k)
    G = np.zeros((N, N), dtype=np.int64)
    for mask, c in g.items():
        a1, a2 = mask & (N - 1), mask >> k
        G[a1, a2] = c % p
    Gh = fft2(G, zeta, k, p)
    if np.any(Gh == 0):
        raise ValueError('g is a zero divisor (proper ideal): the multiplier is not unique')
    Ginv = np.array([[pow(int(x), -1, p) for x in row] for row in Gh], dtype=np.int64)
    out = np.zeros((basis.shape[0], N, N), dtype=np.int64)
    for i in range(basis.shape[0]):
        V = to_array(basis[i], grid, p)
        Q = fft2((fft2(V, zeta, k, p) * Ginv) % p, zeta, k, p, inverse=True)
        out[i] = Q
    return out


def weight_table(k: int) -> np.ndarray:
    N = 1 << k
    w = np.array([int(x).bit_count() for x in range(N)], dtype=np.int64)
    return w[:, None] + w[None, :]


# ----------------------------------------------------------------------------
# containing-space certificate in the group algebra
# ----------------------------------------------------------------------------


def check_containing_space(rows: np.ndarray, g: dict[int, int], grid: Grid, p: int, *, product_cap: int = 20000) -> dict:
    """Given a spanning set `rows` of a subspace W of M_{<=D}: W contains all admissible generator
    rows, is closed under every admissible variable product, and 1 ∉ W  =>  1 ∉ V_D (exclusion
    certificate, sound without trusting how W was built).  Returns the three facts and dim W."""
    out = []
    rank, piv, _ = elim.fp_echelon(rows, p, out)
    W = out[0]
    gen, _ = generator_rows(g, grid, p)
    contains_gen = elim.fp_echelon(np.vstack([W, gen]), p)[0] == rank
    low = [i for i, c in enumerate(piv) if grid.weights[c] < grid.D]
    budget = max(1, product_cap // (2 * grid.k))
    closed = True
    for start in range(0, len(low), budget):
        chunk = low[start:start + budget]
        stacked = np.vstack([W] + [shift_rows(W[chunk], perm) for perm in grid.perms])
        if elim.fp_echelon(stacked, p)[0] != rank:
            closed = False
            break
    one = (len(grid.masks) - 1) in piv
    return dict(dim=int(rank), contains_generator_rows=bool(contains_gen), closed=bool(closed), contains_one=bool(one),
                certifies_one_not_in_VD=bool(contains_gen and closed and not one))


def row_digest(arr: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


# ----------------------------------------------------------------------------
# free high extension: the k-uniform candidate containing space
# ----------------------------------------------------------------------------


def embed_rows(rows: np.ndarray, grid_from: Grid, grid_to: Grid, shift: tuple[int, int] = (0, 0)) -> np.ndarray:
    """Read rows of A_{k0} (exponents in [0, 2^k0)^2) as elements of A_k, k >= k0, translated by the
    exponent pair `shift` (added blockwise modulo 2^k).  Every image monomial must lie in the target
    window (asserted)."""
    N0, N = grid_from.N, grid_to.N
    assert N <= 1 << 60 and N0 <= N
    out = np.zeros((rows.shape[0], len(grid_to.masks)), dtype=np.int64)
    cols = np.flatnonzero(rows.any(axis=0))
    for c in cols:
        a1, a2 = grid_from.pair(grid_from.masks[c])
        b1, b2 = (a1 + shift[0]) % N, (a2 + shift[1]) % N
        out[:, grid_to.column(b1, b2)] = rows[:, c]
    return out


def high_patterns(k: int, k0: int, max_weight: int) -> list[tuple[int, int]]:
    """Exponent pairs supported on levels k0..k-1 (both blocks) of weight <= max_weight."""
    positions = [(i, j) for i in range(2) for j in range(k0, k)]
    out = []
    for w in range(max_weight + 1):
        for choice in combinations(positions, w):
            h = [0, 0]
            for i, j in choice:
                h[i] |= 1 << j
            out.append((h[0], h[1]))
    return out


def free_extension(result: dict, grid0: Grid, grid: Grid, p: int) -> tuple[np.ndarray, dict]:
    """Free_k(V) = sum over high patterns h (levels >= k0, wt h <= D) of tau_h (V ∩ M_{<= D - wt h}),
    for V the row space of `result` (an echelon basis in A_{k0}, columns weight-sorted) and D = grid.D.
    Returns the stacked rows (a spanning set; the sum is direct because distinct high patterns give
    disjoint monomial supports) and the predicted dimension sum_h dim(V ∩ M_{<= D - wt h})."""
    k0, k, D = grid0.k, grid.k, grid.D
    assert grid0.D == D and k >= k0
    parts = []
    predicted = 0
    pieces = {}
    for h in high_patterns(k, k0, D):
        w = h[0].bit_count() + h[1].bit_count()
        low = low_part(result, grid0, D - w)
        if low.shape[0] == 0:
            continue
        predicted += low.shape[0]
        pieces[w] = pieces.get(w, 0) + low.shape[0]
        parts.append(embed_rows(low, grid0, grid, h))
    rows = np.vstack(parts) if parts else np.zeros((0, len(grid.masks)), dtype=np.int64)
    return rows, dict(predicted_dim=int(predicted), rows_by_high_weight=pieces)


def same_row_space(A: np.ndarray, B: np.ndarray, p: int) -> dict:
    ra = elim.fp_echelon(A, p)[0] if A.shape[0] else 0
    rb = elim.fp_echelon(B, p)[0] if B.shape[0] else 0
    rab = elim.fp_echelon(np.vstack([A, B]), p)[0]
    return dict(rank_a=int(ra), rank_b=int(rb), rank_sum=int(rab), equal=bool(ra == rb == rab),
                a_in_b=bool(rab == rb), b_in_a=bool(rab == ra))
