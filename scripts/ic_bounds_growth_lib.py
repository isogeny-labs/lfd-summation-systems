"""Small exact engine for docs/proofs/N8-growth.md (self-contained; numpy only).

Model: A_k = F_p[(Z/2^k)^2]; monomial = exponent pair (a1, a2), weight = popcount(a1) + popcount(a2).
V_D(k) = least subspace of span(M_<=D) containing g = S_3(x1, x2, t) and closed under
f -> x_i^(2^j) f for f of weight < D.  Everything is dense linear algebra mod p (int64), for small k.
"""
from __future__ import annotations
import itertools
import numpy as np


def s3_coeffs(p, a, b, t):
    c = {}
    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p
    add((2, 0), t * t); add((1, 1), -2 * t * t); add((0, 2), t * t)
    add((2, 1), -2 * t); add((1, 2), -2 * t); add((1, 0), -2 * t * a); add((0, 1), -2 * t * a)
    add((0, 0), -4 * t * b)
    add((2, 2), 1); add((1, 1), -2 * a); add((0, 0), a * a)
    add((1, 0), -4 * b); add((0, 1), -4 * b)
    return {e: v for e, v in c.items() if v}


def popcount(x):
    return bin(x).count('1')


class Layout:
    """Monomials of weight <= D at depth k, ordered by weight DESCENDING (ties by (a1, a2))."""
    def __init__(self, k, D):
        self.k, self.D, self.N = k, D, 1 << k
        mons = [(a1, a2) for a1 in range(self.N) for a2 in range(self.N)
                if popcount(a1) + popcount(a2) <= D]
        mons.sort(key=lambda m: (-(popcount(m[0]) + popcount(m[1])), m))
        self.mons = mons
        self.idx = {m: i for i, m in enumerate(mons)}
        self.wt = np.array([popcount(m[0]) + popcount(m[1]) for m in mons])
        self.n = len(mons)
        self.a1 = np.array([m[0] for m in mons]); self.a2 = np.array([m[1] for m in mons])

    def shift_perm(self, i, j):
        """Column map of multiplication by x_i^(2^j) in A_k (with wrap); -1 where the image has weight > D."""
        N = self.N
        out = np.full(self.n, -1, dtype=np.int64)
        for c, (a1, a2) in enumerate(self.mons):
            m = ((a1 + (1 << j)) % N, a2) if i == 0 else (a1, (a2 + (1 << j)) % N)
            out[c] = self.idx.get(m, -1)
        return out

    def vec(self, d, p):
        v = np.zeros(self.n, dtype=np.int64)
        for e, c in d.items():
            v[self.idx[e]] = c % p
        return v


def rref(M, p):
    """Reduced row echelon form mod p (leading = leftmost column). Returns (R, pivots)."""
    M = np.array(M, dtype=np.int64) % p
    r = 0
    piv = []
    rows, cols = M.shape
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0:
            continue
        s = r + nz[0]
        if s != r:
            M[[r, s]] = M[[s, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        piv.append(c)
        r += 1
    return M[:r], piv


def reduce_against(R, piv, Y, p):
    """Reduce rows Y against RREF (R, piv)."""
    Y = np.array(Y, dtype=np.int64) % p
    if len(piv) == 0 or Y.size == 0:
        return Y
    C = Y[:, piv]
    return (Y - (C @ R) % p) % p if C.any() else Y


def nullspace(E, p):
    """Basis (rows) of {v : E v = 0} mod p, in RREF w.r.t. the given column order."""
    R, piv = rref(E, p)
    n = E.shape[1]
    free = [c for c in range(n) if c not in set(piv)]
    Nb = np.zeros((len(free), n), dtype=np.int64)
    for r_, f in enumerate(free):
        Nb[r_, f] = 1
        if piv:
            Nb[r_, piv] = (-R[:, f]) % p
    return Nb


def merge(B, piv, Y, p):
    """Add rows Y to the RREF basis (B, piv); returns the new RREF basis, rows sorted by pivot column."""
    res = reduce_against(B, piv, Y, p)
    res = res[np.any(res != 0, axis=1)]
    if res.shape[0] == 0:
        return B, piv, []
    C, pc = rref(res, p)
    # a row of B has zeros left of its pivot, so only rows of C with a larger pivot column are subtracted
    if len(piv):
        B = (B - (B[:, pc] @ C) % p) % p
    allp = list(piv) + list(pc)
    order = np.argsort(allp)
    B = np.vstack([B, C])[order]
    return B, [allp[i] for i in order], pc


def closure(L, g, p, shifts=None, verbose=False, chunk=256):
    """V_D(k): RREF basis (rows) of the closure of g. shifts: list of (i, j); default all 2k.

    Columns are in weight-descending order, so the rows with a pivot of weight < D span V ∩ M_<D.  Every such
    row is processed (all its shifts added) once, in its state at that time; a low row is afterwards modified
    only by subtracting later low rows (which are processed too), so the final low rows lie in the span of the
    processed states and the final space is closed.  ic_bounds_growth_ideal.py re-checks closedness of the
    result by a full sweep (`verify_closed`)."""
    k, D = L.k, L.D
    if shifts is None:
        shifts = [(i, j) for i in range(2) for j in range(k)]
    perms = {s: L.shift_perm(*s) for s in shifts}
    B, piv = rref(L.vec(g, p)[None, :], p)
    processed = set()
    while True:
        todo = [c for c in piv if L.wt[c] < D and c not in processed][:chunk]
        if not todo:
            break
        pos = {c: r for r, c in enumerate(piv)}
        new = B[[pos[c] for c in todo]]
        processed.update(todo)
        prods = []
        for s in shifts:
            pm = perms[s]
            Y = np.zeros_like(new)
            src = np.nonzero(pm >= 0)[0]
            Y[:, pm[src]] = new[:, src]
            assert not new[:, pm < 0].any()
            prods.append(Y)
        B, piv, _ = merge(B, piv, np.vstack(prods), p)
        if verbose:
            print('  closure rank', len(piv), 'low', sum(1 for c in piv if L.wt[c] < D), flush=True)
    return B, piv


def verify_closed(L, B, piv, p, shifts=None):
    """Full sweep: number of (low row, shift) products that do not reduce to zero against (B, piv)."""
    k, D = L.k, L.D
    if shifts is None:
        shifts = [(i, j) for i in range(2) for j in range(k)]
    low = B[[r for r, c in enumerate(piv) if L.wt[c] < D]]
    bad = 0
    for s in shifts:
        pm = L.shift_perm(*s)
        Y = np.zeros_like(low); src = np.nonzero(pm >= 0)[0]
        Y[:, pm[src]] = low[:, src]
        bad += int(np.any(reduce_against(B, piv, Y, p) != 0, axis=1).sum())
    return bad


def profile(L, piv):
    return [sum(1 for c in piv if L.wt[c] <= j) for j in range(L.D + 1)]


def curve_points(p, a, b, t, limit=None):
    """Affine F_p-points (x1, x2) of S_3(x1, x2, t) = 0."""
    g = s3_coeffs(p, a, b, t)
    sq = {}
    for y in range(p):
        sq.setdefault(y * y % p, y)
    pts = []
    for x1 in range(p):
        co = [0, 0, 0]
        for (e1, e2), c in g.items():
            co[e2] = (co[e2] + c * pow(x1, e1, p)) % p
        c0, c1, c2 = co
        if c2 == 0:
            if c1:
                pts.append((x1, (-c0 * pow(c1, p - 2, p)) % p))
            continue
        disc = (c1 * c1 - 4 * c2 * c0) % p
        if disc not in sq:
            continue
        s = sq[disc]; inv = pow(2 * c2, p - 2, p)
        for r in {s, (-s) % p}:
            pts.append((x1, ((-c1 + r) * inv) % p))
        if limit and len(pts) >= limit:
            break
    return pts


def eval_matrix(L, pts, p):
    """Rows = points, columns = monomials of L evaluated as integer-exponent monomials x1^a1 x2^a2."""
    N = L.N
    P1 = np.ones((len(pts), N), dtype=np.int64); P2 = np.ones((len(pts), N), dtype=np.int64)
    x1 = np.array([q[0] for q in pts], dtype=np.int64); x2 = np.array([q[1] for q in pts], dtype=np.int64)
    for e in range(1, N):
        P1[:, e] = P1[:, e - 1] * x1 % p; P2[:, e] = P2[:, e - 1] * x2 % p
    return P1[:, L.a1] * P2[:, L.a2] % p
