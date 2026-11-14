"""Small exact engine for two-summand prime chain systems modulo the chains.

Independent of csd/: the chain quotient is the group algebra A = F_p[(Z/N1) x (Z/N2)],
N_i = 2^k_i; a monomial is an exponent pair (e1, e2), its degree is wt(e1) + wt(e2) (binary
Hamming weight), the chain variable u_{i,j} acts by e_i -> e_i + 2^j (mod N_i).  The generator is
g = S_3(x1, x2, t) (exponents 0..2).  V_D is the smallest subspace of M_<=D containing g and
closed under u_{i,j} applied to its elements of degree < D (the quotient closure of
formal/ICFallDegrees/ChainClosureBridge.lean).

All arithmetic is exact modulo p; float64 BLAS products are used only where p^2 * rank < 2^53.
"""
from __future__ import annotations

import numpy as np


def s3_coeffs(p, a, b, t):
    """S_3(x1, x2, t) for y^2 = x^3 + a x + b as {(e1, e2): c}."""
    c = {}

    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p

    add((2, 0), t * t); add((1, 1), -2 * t * t); add((0, 2), t * t)
    add((2, 1), -2 * t); add((1, 2), -2 * t); add((1, 0), -2 * t * a); add((0, 1), -2 * t * a)
    add((0, 0), -4 * t * b)
    add((2, 2), 1); add((1, 1), -2 * a); add((0, 0), a * a)
    add((1, 0), -4 * b); add((0, 1), -4 * b)
    return {e: v for e, v in c.items() if v}


def s3_value(p, a, b, t, x1, x2):
    return ((x1 - x2) ** 2 * t * t - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t
            + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p


def subgroup(p, k):
    N = 1 << k
    assert (p - 1) % N == 0
    if k == 0:
        return 1, [1]
    for z in range(2, p):
        h = pow(z, (p - 1) // N, p)
        if pow(h, N // 2, p) != 1:
            return h, [pow(h, e, p) for e in range(N)]
    raise ValueError


def solutions(p, a, b, t, k1, k2=None):
    """Z as exponent pairs (E1, E2) w.r.t. the subgroup generators, and as values."""
    k2 = k1 if k2 is None else k2
    _, H1 = subgroup(p, k1)
    _, H2 = subgroup(p, k2)
    return [(E1, E2) for E1, x1 in enumerate(H1) for E2, x2 in enumerate(H2)
            if s3_value(p, a, b, t, x1, x2) == 0]


def wt(e):
    return bin(e).count('1')


class Window:
    """Monomials of degree <= D, heaviest first (so a row's leading column carries its degree)."""

    def __init__(self, k1, k2, D):
        self.k1, self.k2, self.D = k1, k2, D
        self.N1, self.N2 = 1 << k1, 1 << k2
        mons = [(e1, e2) for e1 in range(self.N1) for e2 in range(self.N2) if wt(e1) + wt(e2) <= D]
        mons.sort(key=lambda m: (-(wt(m[0]) + wt(m[1])), m))
        self.mons = mons
        self.n = len(mons)
        self.index = {m: i for i, m in enumerate(mons)}
        self.weights = np.array([wt(m[0]) + wt(m[1]) for m in mons], dtype=np.int64)
        self.vars = [(0, j) for j in range(k1)] + [(1, j) for j in range(k2)]
        self.perms = []
        for (i, j) in self.vars:
            perm = np.full(self.n, -1, dtype=np.int64)
            for c, (e1, e2) in enumerate(mons):
                m = ((e1 + (1 << j)) % self.N1, e2) if i == 0 else (e1, (e2 + (1 << j)) % self.N2)
                perm[c] = self.index.get(m, -1)
            self.perms.append(perm)

    def row(self, poly, p):
        r = np.zeros(self.n, dtype=np.int64)
        for m, c in poly.items():
            r[self.index[m]] = (r[self.index[m]] + c) % p
        return r

    def shift(self, rows, vi):
        perm = self.perms[vi]
        assert not np.any(rows[:, perm < 0]), 'product leaves the degree window'
        out = np.zeros_like(rows)
        keep = perm >= 0
        out[:, perm[keep]] = rows[:, keep]
        return out

    def eval_matrix(self, pts, p):
        """Rows = points (E1, E2), columns = monomials; entry h1^(E1 e1) h2^(E2 e2)."""
        _, H1 = subgroup(p, self.k1)
        _, H2 = subgroup(p, self.k2)
        M = np.zeros((len(pts), self.n), dtype=np.int64)
        for r, (E1, E2) in enumerate(pts):
            for c, (e1, e2) in enumerate(self.mons):
                M[r, c] = H1[(E1 * e1) % self.N1] * H2[(E2 * e2) % self.N2] % p
        return M


def _inv(x, p):
    return pow(int(x), p - 2, p)


class Echelon:
    """Incremental reduced row echelon basis modulo p (leftmost pivots)."""

    def __init__(self, n, p):
        self.n, self.p = n, p
        self.B = np.zeros((0, n), dtype=np.int64)
        self.piv = []

    @property
    def rank(self):
        return len(self.piv)

    def reduce(self, rows):
        rows = np.asarray(rows, dtype=np.int64).reshape((-1, self.n)) % self.p
        if self.rank and len(rows):
            assert self.p * self.p * max(1, self.rank) < 2 ** 52
            coef = rows[:, self.piv].astype(np.float64)
            rows = (rows - np.rint(coef @ self.B.astype(np.float64)).astype(np.int64) % self.p) % self.p
        return rows

    def add(self, rows, chunk=256):
        added = 0
        rows = np.asarray(rows, dtype=np.int64).reshape((-1, self.n))
        for s in range(0, len(rows), chunk):
            added += self._add(rows[s:s + chunk])
        return added

    def _add(self, rows):
        p = self.p
        R = self.reduce(rows)
        new, newpiv = [], []
        for i in range(len(R)):
            r = R[i]
            nz = np.flatnonzero(r)
            if not nz.size:
                continue
            c = int(nz[0])
            r = r * _inv(r[c], p) % p
            R[i] = r
            # clear column c in the later rows and in the rows already accepted
            f = R[i + 1:, c].copy()
            if f.any():
                R[i + 1:] = (R[i + 1:] - np.outer(f, r)) % p
            for q in range(len(new)):
                if new[q][c]:
                    new[q] = (new[q] - new[q][c] * r) % p
            new.append(r)
            newpiv.append(c)
        if not new:
            return 0
        Nw = np.array(new, dtype=np.int64)
        if self.rank:
            coef = self.B[:, newpiv].astype(np.float64)
            self.B = (self.B - np.rint(coef @ Nw.astype(np.float64)).astype(np.int64) % p) % p
        self.B = np.vstack([self.B, Nw])
        self.piv = self.piv + newpiv
        order = np.argsort(self.piv)
        self.B = self.B[order]
        self.piv = [self.piv[i] for i in order]
        return len(new)

    def contains(self, rows):
        return not np.any(self.reduce(rows))


def closure(win, gens, p, log=None):
    """V_D in the window: generators (rows) plus variable shifts of every element of degree < D.

    Loop: shift every basis row whose leading (= heaviest) column has degree < D and that has
    not been shifted in its present form; stop when a full pass over the final basis adds nothing
    (so the returned space is verified closed).  Returns the Echelon."""
    E = Echelon(win.n, p)
    E.add(gens)
    done = set()
    D = win.D
    while True:
        low = [i for i, c in enumerate(E.piv) if win.weights[c] < D]
        todo = [i for i in low if E.B[i].tobytes() not in done]
        if not todo:
            break
        rows = E.B[todo].copy()
        for r in rows:
            done.add(r.tobytes())
        grew = 0
        for vi in range(len(win.vars)):
            grew += E.add(win.shift(rows, vi))
        if log:
            log(dict(rank=E.rank, shifted=len(todo), grew=grew))
        if E.rank == win.n:
            break
    return E


def profile(E, win):
    """dim(V ∩ M_<=j), j = 0..D."""
    w = win.weights[E.piv] if E.rank else np.array([], dtype=np.int64)
    return [int(np.sum(w <= j)) for j in range(win.D + 1)]


def rank_mod_p(M, p):
    E = Echelon(M.shape[1], p)
    E.add(M)
    return E.rank


def kernel_basis(M, p):
    """Basis (rows) of {v : M v = 0} modulo p."""
    n = M.shape[1]
    E = Echelon(n, p)
    E.add(M)
    piv = E.piv
    pset = set(piv)
    free = [c for c in range(n) if c not in pset]
    K = np.zeros((len(free), n), dtype=np.int64)
    for i, f in enumerate(free):
        K[i, f] = 1
        if piv:
            K[i, piv] = (-E.B[:, f]) % p
    return K


def generator_row(win, p, a, b, t):
    g = {}
    for (e1, e2), c in s3_coeffs(p, a, b, t).items():
        m = (e1 % win.N1, e2 % win.N2)
        g[m] = (g.get(m, 0) + c) % p
    g = {m: c for m, c in g.items() if c}
    return win.row(g, p), g
