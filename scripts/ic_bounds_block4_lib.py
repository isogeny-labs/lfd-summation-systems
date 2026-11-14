"""Shared engine for docs/proofs/N9-block-four.md (conjectural block-of-four lemmas).

Model: A_{k1,k2} = F_p[Z/2^k1 x Z/2^k2]; monomial = exponent pair (a1, a2), weight = popcount(a1) + popcount(a2);
the shift u_{i,j} multiplies by x_i^(2^j) (exponent addition mod 2^{k_i}).  V_D = least subspace of M_<=D
containing the generator g and closed under the shifts on elements of weight < D.  Equal depths k1 = k2 = k give
the repository's V_D(k).

The closure loop, RREF and full-sweep check are those of scripts/ic_bounds_growth_lib.py (imported, unchanged);
only the layout is generalized to unequal depths.  Dense int64 arithmetic mod p (p < 2^17 keeps every product
sum below 2^63 for the sizes used here).
"""
from __future__ import annotations

import hashlib
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ic_bounds_growth_lib as GL  # noqa: E402

popcount = GL.popcount
rref = GL.rref
reduce_against = GL.reduce_against
merge = GL.merge
s3_coeffs = GL.s3_coeffs

# Unit instance of experiments/ic-bounds-k14 (y^2 = x^3 + x over F_65537, t = x(R) = 4271; empty solution set at
# every depth <= 15, re-checked by `unit_check` for the depths used).
UNIT = dict(p=65537, a=1, b=0, t=4271)
# Second unit instance of the same file (y^2 = x^3 + 3x over F_163841, t = 73730; empty for depths <= 14).
UNIT2 = dict(p=163841, a=3, b=0, t=73730)


class Layout2:
    """Monomials of weight <= D at depths (k1, k2), weight-descending (ties by (a1, a2))."""

    def __init__(self, k1, k2, D):
        self.k1, self.k2, self.D = k1, k2, D
        self.N1, self.N2 = 1 << k1, 1 << k2
        # attributes used by growth_lib helpers
        self.k = max(k1, k2)
        mons = []
        for a1 in range(self.N1):
            w1 = popcount(a1)
            if w1 > D:
                continue
            for a2 in range(self.N2):
                if w1 + popcount(a2) <= D:
                    mons.append((a1, a2))
        mons.sort(key=lambda m: (-(popcount(m[0]) + popcount(m[1])), m))
        self.mons = mons
        self.idx = {m: i for i, m in enumerate(mons)}
        self.wt = np.array([popcount(m[0]) + popcount(m[1]) for m in mons])
        self.n = len(mons)
        self.a1 = np.array([m[0] for m in mons]); self.a2 = np.array([m[1] for m in mons])
        self.shifts = [(0, j) for j in range(k1)] + [(1, j) for j in range(k2)]

    def shift_perm(self, i, j):
        out = np.full(self.n, -1, dtype=np.int64)
        for c, (a1, a2) in enumerate(self.mons):
            m = ((a1 + (1 << j)) % self.N1, a2) if i == 0 else (a1, (a2 + (1 << j)) % self.N2)
            out[c] = self.idx.get(m, -1)
        return out

    def vec(self, d, p):
        v = np.zeros(self.n, dtype=np.int64)
        for (e1, e2), c in d.items():
            key = (e1 % self.N1, e2 % self.N2)
            v[self.idx[key]] = (v[self.idx[key]] + c) % p
        return v


def generator(inst):
    return s3_coeffs(inst['p'], inst['a'], inst['b'], inst['t'])


def closure(L, g, p, verbose=False, chunk=256, extra=None):
    """V_D at depths (L.k1, L.k2); `extra` = further input vectors (weight <= D) added to the generator."""
    B, piv = GL.rref(L.vec(g, p)[None, :], p)
    if extra is not None and len(extra):
        B, piv, _ = GL.merge(B, piv, np.array(extra, dtype=np.int64), p)
    D = L.D
    perms = {s: L.shift_perm(*s) for s in L.shifts}
    processed = set()
    while True:
        todo = [c for c in piv if L.wt[c] < D and c not in processed][:chunk]
        if not todo:
            break
        pos = {c: r for r, c in enumerate(piv)}
        new = B[[pos[c] for c in todo]]
        processed.update(todo)
        prods = []
        for s in L.shifts:
            pm = perms[s]
            Y = np.zeros_like(new)
            src = np.nonzero(pm >= 0)[0]
            Y[:, pm[src]] = new[:, src]
            assert not new[:, pm < 0].any()
            prods.append(Y)
        B, piv, _ = GL.merge(B, piv, np.vstack(prods), p)
        if verbose:
            print('  closure rank', len(piv), flush=True)
    return B, piv


def verify_closed(L, B, piv, p):
    """Full sweep: number of (low row, shift) products not reducing to zero (0 = closed)."""
    low = B[[r for r, c in enumerate(piv) if L.wt[c] < L.D]]
    bad = 0
    for s in L.shifts:
        pm = L.shift_perm(*s)
        Y = np.zeros_like(low); src = np.nonzero(pm >= 0)[0]
        Y[:, pm[src]] = low[:, src]
        if Y.shape[0]:
            bad += int(np.any(GL.reduce_against(B, piv, Y, p) != 0, axis=1).sum())
    return bad


def profile(L, piv):
    return [sum(1 for c in piv if L.wt[c] <= j) for j in range(L.D + 1)]


def member(L, B, piv, d, p):
    """Is the element with exponent->coefficient dict d in the row space?"""
    v = L.vec(d, p)  # raises KeyError if d has a monomial of weight > D (then it is not in V_D)
    return not np.any(GL.reduce_against(B, piv, v[None, :], p))


def has_one(L, piv):
    return L.idx[(0, 0)] in set(piv)


def basis_hash(B):
    return hashlib.sha256(np.ascontiguousarray(B % (1 << 32)).astype(np.uint32).tobytes()).hexdigest()[:16]


def root_of_unity(p, k):
    N = 1 << k
    assert (p - 1) % N == 0
    for t in range(2, p):
        h = pow(t, (p - 1) // N, p)
        if pow(h, N // 2, p) != 1:
            return h
    raise ValueError


def unit_check(inst, k1, k2):
    """Number of (alpha, beta) in mu_{2^k1} x mu_{2^k2} with S_3(alpha, beta, t) = 0."""
    p = inst['p']
    co = generator(inst)
    h1 = root_of_unity(p, k1); h2 = root_of_unity(p, k2)
    H2 = np.array([pow(h2, e, p) for e in range(1 << k2)], dtype=np.int64)
    pw2 = [np.ones_like(H2), H2, H2 * H2 % p]
    cnt = 0
    for e in range(1 << k1):
        x1 = pow(h1, e, p)
        pw1 = [1, x1, x1 * x1 % p]
        v = np.zeros_like(H2)
        for (e1, e2), c in co.items():
            v = (v + (c * pw1[e1] % p) * pw2[e2]) % p
        cnt += int(np.count_nonzero(v == 0))
    return cnt


# --------------------------------------------------------------------------- fast engine (unequal depths)
# Uses the block-RREF basis class `Part` of csd/ic_bounds_d5_verify.py (imported, unchanged; float64 BLAS with
# exact integer arithmetic for p < 2^24).  No symmetry split (unequal depths break the swap x1 <-> x2).

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
from csd.ic_bounds_d5_verify import Part  # noqa: E402

F64 = np.float64


class FastLayout:
    def __init__(self, k1, k2, D):
        self.k1, self.k2, self.D = k1, k2, D
        N1, N2 = 1 << k1, 1 << k2
        self.N1, self.N2 = N1, N2
        pc1 = np.array([popcount(x) for x in range(N1)]); pc2 = np.array([popcount(x) for x in range(N2)])
        s1 = np.flatnonzero(pc1 <= D); s2 = np.flatnonzero(pc2 <= D)
        A1, A2 = np.meshgrid(s1, s2, indexing='ij'); A1, A2 = A1.ravel(), A2.ravel()
        W = pc1[A1] + pc2[A2]
        keep = W <= D
        A1, A2, W = A1[keep], A2[keep], W[keep]
        key = A1 | (A2 << k1)
        order = np.lexsort((key, -W))
        self.a1, self.a2, self.w, self.key = A1[order], A2[order], W[order], key[order]
        self.c = len(self.key)
        self._sorted = np.argsort(self.key); self._skeys = self.key[self._sorted]
        self.low = np.flatnonzero(self.w < D); self.lowc = len(self.low)
        assert np.all(self.low == np.arange(self.c - self.lowc, self.c))
        self.shift_list = [(0, j) for j in range(k1)] + [(1, j) for j in range(k2)]
        self.shift = {}
        for (i, j) in self.shift_list:
            b1 = self.a1[self.low].copy(); b2 = self.a2[self.low].copy()
            if i == 0:
                b1 = (b1 + (1 << j)) % N1
            else:
                b2 = (b2 + (1 << j)) % N2
            self.shift[(i, j)] = self.index(b1, b2)

    def index(self, b1, b2):
        kk = np.asarray(b1) | (np.asarray(b2) << self.k1)
        pos = np.searchsorted(self._skeys, kk)
        if not np.all(self._skeys[np.minimum(pos, self.c - 1)] == kk):
            raise KeyError('monomial of weight > D')
        return self._sorted[pos]

    def vec(self, d, p):
        v = np.zeros(self.c, dtype=F64)
        for (e1, e2), c in d.items():
            col = int(self.index(np.array([e1 % self.N1]), np.array([e2 % self.N2]))[0])
            v[col] = (v[col] + c) % p
        return v

    def shifted(self, F, ij):
        out = np.zeros((F.shape[0], self.c), dtype=F64)
        out[:, self.shift[ij]] = F[:, self.c - self.lowc:]
        return out


class FastClosure:
    def __init__(self, k1, k2, D, p, inputs):
        """inputs: list of exponent->coefficient dicts (weight <= D)."""
        self.L = FastLayout(k1, k2, D)
        self.p, self.D = p, D
        self.P = Part(self.L.w, p)
        self.processed = set()
        H = np.array([self.L.vec(d, p) for d in inputs])
        self.P.add(H)

    def low_rows(self):
        return np.flatnonzero(self.P.w[self.P.piv] < self.D)

    def _shifts(self, idx):
        F = self.P.rows(idx)
        assert not np.any(F[:, :self.L.c - self.L.lowc])
        return np.vstack([self.L.shifted(F, ij) for ij in self.L.shift_list])

    def run(self, batch=64):
        while True:
            todo = [i for i in self.low_rows() if int(self.P.piv[i]) not in self.processed][:batch]
            if not todo:
                return
            for i in todo:
                self.processed.add(int(self.P.piv[i]))
            self.P.add(self._shifts(todo))

    def sweep_full(self, chunk=64):
        bad = 0
        low = self.low_rows()
        for a in range(0, len(low), chunk):
            Y = self._shifts(low[a:a + chunk])
            bad += int(np.count_nonzero(np.any(self.P.residual(Y) != 0, axis=1)))
        return bad

    def contains(self, d):
        try:
            v = self.L.vec(d, self.p)
        except KeyError:
            return False
        return not np.any(self.P.residual(v[None, :]))

    def contains_one(self):
        return self.contains({(0, 0): 1})

    def profile(self):
        return self.P.profile(self.D)

    def dim(self):
        return self.P.rank
