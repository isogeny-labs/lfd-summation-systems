"""Independent re-implementation of the degree-D chain closure, written for the verification of
item 3 ("degree-5 closure growth").  It does NOT import csd/ic_bounds_d5.py, the producer's scripts,
csd/ic_fall_certificates.py, csd/ic_chain_shift_closure.py or csd/ic_chain_block_closure.py.

Setting.  Two summands, depth k, prime p with 2^k | p - 1.  Modulo the chains the ring is the group
algebra F_p[(Z/2^k)^2]; a squarefree chain monomial is an exponent pair (a1, a2), its degree is the
Hamming weight popcount(a1) + popcount(a2), and multiplication by u_{i,j} adds 2^j to a_i mod 2^k.
V_D is the least subspace of span(M_<=D) containing NF(g) and closed under f -> u_{i,j} f for every
f of weight < D.

Engine.
  * columns ordered by weight DESCENDING (ties by the integer a1 | a2 << k), so the leading column of
    a row is a column of maximal weight in its support;
  * optional summand-swap split (sigma(a1, a2) = (a2, a1)); V_D = V^+ (+) V^- when g is symmetric;
    parts are kept in orbit coordinates; with the split only the k block-1 shifts are applied and
    each shifted row is projected to both parts (h + sigma h, h - sigma h);
  * each part is a reduced row echelon basis stored as (pivot columns, coefficients on the free
    columns), float32 storage, exact float64 products (every partial sum < 2^53);
  * closure: every row whose pivot has weight < D is processed (its shifts added) once, at its
    current state; the loop ends when every low pivot has been processed.  Low rows are only ever
    modified by later low rows, so the final low rows lie in the span of the processed states and
    the final space is closed (the invariant argument in VERIFY.md);
  * certificates computed independently of that argument: a randomized sweep (t random
    combinations of the low rows of each part, all shifts, all projections, must reduce to zero;
    a non-closed space passes with probability <= 4k p^-t) and, when asked, the full
    deterministic sweep (every shift of every final low row must reduce to zero).
"""
from __future__ import annotations

import hashlib
import time

import numpy as np

F64 = np.float64


# --------------------------------------------------------------------------- instance (own code)

def s3_coeffs(p: int, a: int, b: int, X: int) -> dict:
    """S_3(x1, x2, X) for y^2 = x^3 + a x + b, expanded here from
    (x1 - x2)^2 X^2 - 2((x1 + x2)(x1 x2 + a) + 2b) X + (x1 x2 - a)^2 - 4b(x1 + x2)."""
    c = {}

    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p

    # (x1 - x2)^2 X^2 = X^2 x1^2 - 2 X^2 x1 x2 + X^2 x2^2
    add((2, 0), X * X); add((1, 1), -2 * X * X); add((0, 2), X * X)
    # -2X[(x1 + x2)(x1 x2 + a) + 2b] = -2X[x1^2 x2 + x1 x2^2 + a x1 + a x2 + 2b]
    add((2, 1), -2 * X); add((1, 2), -2 * X); add((1, 0), -2 * X * a); add((0, 1), -2 * X * a)
    add((0, 0), -4 * X * b)
    # (x1 x2 - a)^2 = x1^2 x2^2 - 2a x1 x2 + a^2
    add((2, 2), 1); add((1, 1), -2 * a); add((0, 0), a * a)
    # -4b(x1 + x2)
    add((1, 0), -4 * b); add((0, 1), -4 * b)
    return {e: v for e, v in c.items() if v}


def chain_generator(p: int, a: int, b: int, X: int, k: int) -> dict:
    """NF of S_3 modulo the chains: x_i^e -> exponent e mod 2^k in block i.  {(a1, a2): coeff}."""
    N = 1 << k
    g = {}
    for (e1, e2), v in s3_coeffs(p, a, b, X).items():
        key = (e1 % N, e2 % N)
        g[key] = (g.get(key, 0) + v) % p
    return {e: v for e, v in g.items() if v}


def subgroup(p: int, k: int):
    """A generator h of the subgroup of order 2^k of F_p^* and the list h^E, E = 0..2^k-1."""
    N = 1 << k
    assert (p - 1) % N == 0
    for t in range(2, p):
        h = pow(t, (p - 1) // N, p)
        if pow(h, N // 2, p) != 1:
            return h, [pow(h, e, p) for e in range(N)]
    raise ValueError


def solutions(p: int, a: int, b: int, X: int, k: int) -> list:
    """All (E1, E2) with S_3(h^E1, h^E2, X) = 0 (vectorized over E2)."""
    h, H = subgroup(p, k)
    Hv = np.array(H, dtype=np.int64)
    co = s3_coeffs(p, a, b, X)
    out = []
    pw2 = [np.ones_like(Hv), Hv, Hv * Hv % p]
    for E1, x1 in enumerate(H):
        pw1 = [1, x1, x1 * x1 % p]
        v = np.zeros_like(Hv)
        for (e1, e2), c in co.items():
            v = (v + (c * pw1[e1] % p) * pw2[e2]) % p
        for E2 in np.flatnonzero(v == 0):
            out.append((E1, int(E2)))
    return out


# --------------------------------------------------------------------------- layout

class Layout:
    def __init__(self, k: int, D: int, split: bool):
        self.k, self.D, self.split = k, D, split
        N = 1 << k
        pc = np.array([bin(x).count('1') for x in range(N)], dtype=np.int64)
        small = np.flatnonzero(pc <= D)
        A1, A2 = np.meshgrid(small, small, indexing='ij')
        A1, A2 = A1.ravel(), A2.ravel()
        W = pc[A1] + pc[A2]
        keep = W <= D
        A1, A2, W = A1[keep], A2[keep], W[keep]
        key = A1 | (A2 << k)
        order = np.lexsort((key, -W))       # weight descending, then key ascending
        self.a1, self.a2, self.w, self.key = A1[order], A2[order], W[order], key[order]
        self.c = len(self.key)
        self._sorted = np.argsort(self.key)
        self._skeys = self.key[self._sorted]
        self.low = np.flatnonzero(self.w < D)            # full columns of weight < D
        self.lowc = len(self.low)
        assert np.all(self.low == np.arange(self.c - self.lowc, self.c))
        # shifts u_{i,j} on low columns
        self.shift = {}
        for i in (0, 1):
            for j in range(k):
                b1 = self.a1[self.low]
                b2 = self.a2[self.low]
                if i == 0:
                    b1 = (b1 + (1 << j)) % N
                else:
                    b2 = (b2 + (1 << j)) % N
                self.shift[(i, j)] = self.index(b1, b2)
        if split:
            sig = self.index(self.a2, self.a1)
            self.sigma = sig
            symrep = np.flatnonzero(self.a1 <= self.a2)
            antirep = np.flatnonzero(self.a1 < self.a2)
            # part columns are full columns in the inherited order (weight descending)
            self.parts = []
            for reps, sign in ((symrep, 1), (antirep, -1)):
                pos = -np.ones(self.c, dtype=np.int64)
                pos[reps] = np.arange(len(reps))
                lift_idx = np.where(pos >= 0, pos, pos[sig])
                lift_sign = np.where(pos >= 0, 1, sign).astype(F64)
                if sign == -1:
                    lift_sign[self.a1 == self.a2] = 0.0
                    lift_idx[self.a1 == self.a2] = 0
                self.parts.append(dict(reps=reps, rep_sig=sig[reps], sign=sign, w=self.w[reps],
                                       lift_idx=lift_idx, lift_sign=lift_sign, c=len(reps)))
            self.shift_list = [(0, j) for j in range(k)]
        else:
            self.parts = [dict(reps=np.arange(self.c), w=self.w, c=self.c)]
            self.shift_list = [(i, j) for i in (0, 1) for j in range(k)]

    def index(self, b1, b2):
        kk = np.asarray(b1) | (np.asarray(b2) << self.k)
        pos = np.searchsorted(self._skeys, kk)
        assert np.all(self._skeys[np.minimum(pos, self.c - 1)] == kk)
        return self._sorted[pos]

    # part coordinates <-> full coordinates
    def lift(self, part: int, Y: np.ndarray) -> np.ndarray:
        if not self.split:
            return Y
        P = self.parts[part]
        return Y[:, P['lift_idx']] * P['lift_sign']

    def project(self, H: np.ndarray, p: int) -> list:
        if not self.split:
            return [np.mod(H, p)]
        out = []
        for P in self.parts:
            # entries in (-p, 2p): exact in every later float64 product; residual() reduces mod p
            Y = H[:, P['reps']]
            if P['sign'] == 1:
                Y += H[:, P['rep_sig']]
            else:
                Y -= H[:, P['rep_sig']]
            out.append(Y)
        return out

    def shifted(self, Hlow_full: np.ndarray, ij) -> np.ndarray:
        """u_{i,j} times rows supported on weight < D (full coordinates)."""
        out = np.zeros((Hlow_full.shape[0], self.c), dtype=F64)
        out[:, self.shift[ij]] = Hlow_full[:, self.c - self.lowc:]
        return out


# --------------------------------------------------------------------------- echelon

def _rref_rows(S: np.ndarray, p: int):
    """Row-by-row RREF of a small block (rows already reduced against everything else)."""
    f = S.shape[1]
    L = np.zeros((S.shape[0], f), dtype=F64)
    pl = []
    cnt = 0
    for v in S:
        v = v.copy()
        if cnt:
            v = np.mod(v - v[pl] @ L[:cnt], p)
        nz = np.flatnonzero(v)
        if nz.size == 0:
            continue
        lead = int(nz[0])
        v = np.mod(v * pow(int(v[lead]), -1, p), p)
        if cnt:
            L[:cnt] = np.mod(L[:cnt] - np.outer(L[:cnt, lead], v), p)
        L[cnt] = v
        pl.append(lead)
        cnt += 1
    return L[:cnt], pl


def rref_block(R: np.ndarray, p: int, sub: int = 32):
    """RREF of R (values in [0,p)); returns (L, pivots) with pivot = leading column of each row."""
    f = R.shape[1]
    L = np.zeros((0, f), dtype=F64)
    pl: list = []
    for a in range(0, R.shape[0], sub):
        S = R[a:a + sub]
        if pl:
            S = np.mod(S - S[:, pl] @ L, p)
        S = S[np.any(S != 0, axis=1)]
        if S.shape[0] == 0:
            continue
        Ln, pn = _rref_rows(S, p)
        if not pn:
            continue
        if pl:
            L = np.mod(L - L[:, pn] @ Ln, p)
        L = np.vstack([L, Ln])
        pl = pl + pn
    return L, pl


class Part:
    """RREF basis: row i = e_{piv[i]} + B[i] on the columns free (sorted)."""

    def __init__(self, w: np.ndarray, p: int, tile_entries: int = 1 << 23):
        self.w, self.p, self.c = w, p, len(w)
        self.piv = np.zeros(0, dtype=np.int64)
        self.free = np.arange(self.c)
        self.B = np.zeros((0, self.c), dtype=np.float32)
        self.tile_entries = tile_entries
        self.madds = 0
        self.peak_entries = 0
        assert p < (1 << 24)

    @property
    def rank(self):
        return len(self.piv)

    def _tile(self):
        return max(64, self.tile_entries // max(1, self.B.shape[1]))

    def residual(self, Y: np.ndarray) -> np.ndarray:
        p = self.p
        R = Y[:, self.free].astype(F64, copy=True)
        r = self.rank
        if r:
            P = Y[:, self.piv]
            t = self._tile()
            for t0 in range(0, r, t):
                R -= P[:, t0:t0 + t] @ self.B[t0:t0 + t].astype(F64)
            self.madds += Y.shape[0] * r * R.shape[1]
        np.mod(R, p, out=R)
        return R

    def add(self, Y: np.ndarray) -> int:
        R = self.residual(Y)
        R = R[np.any(R != 0, axis=1)]
        if R.shape[0] == 0:
            return 0
        L, q = rref_block(R, self.p)
        if not q:
            return 0
        q = np.array(q, dtype=np.int64)
        keep = np.ones(len(self.free), dtype=bool)
        keep[q] = False
        r, t = self.rank, len(q)
        NB = np.empty((r + t, int(keep.sum())), dtype=np.float32)
        Lk = L[:, keep]
        if r:
            Bq = self.B[:, q].astype(F64)
            tt = self._tile()
            for t0 in range(0, r, tt):
                t1 = min(t0 + tt, r)
                blk = self.B[t0:t1][:, keep].astype(F64) - Bq[t0:t1] @ Lk
                NB[t0:t1] = np.mod(blk, self.p)
            self.madds += r * t * Lk.shape[1]
        NB[r:] = Lk
        self.piv = np.concatenate([self.piv, self.free[q]])
        self.free = self.free[keep]
        self.B = NB
        self.peak_entries = max(self.peak_entries, NB.size)
        return t

    def rows(self, idx) -> np.ndarray:
        idx = np.asarray(idx, dtype=np.int64)
        M = np.zeros((len(idx), self.c), dtype=F64)
        M[np.arange(len(idx)), self.piv[idx]] = 1.0
        M[:, self.free] = self.B[idx]
        return M

    def check_rref(self, sample: int = 64, rng=None) -> bool:
        """pivot columns of the basis form an identity; each row's leading column is its pivot."""
        if self.rank == 0:
            return True
        rng = rng or np.random.default_rng(0)
        idx = rng.choice(self.rank, size=min(sample, self.rank), replace=False)
        M = self.rows(idx)
        if not np.array_equal(M[:, self.piv], np.eye(self.rank)[idx]):
            return False
        lead = np.argmax(M != 0, axis=1)
        return bool(np.all(lead == self.piv[idx]))

    def profile(self, D: int) -> list:
        wp = self.w[self.piv]
        return [int(np.sum(wp <= d)) for d in range(D + 1)]

    def sha(self) -> str:
        """Hash of the canonical RREF (rows sorted by pivot column)."""
        o = np.argsort(self.piv)
        h = hashlib.sha256()
        h.update(self.piv[o].astype(np.int64).tobytes())
        h.update(self.free.astype(np.int64).tobytes())
        h.update(np.ascontiguousarray(self.B[o]).astype(np.uint16 if self.p < 65536 else np.uint32).tobytes())
        return h.hexdigest()


# --------------------------------------------------------------------------- closure

class Closure:
    def __init__(self, g: dict, k: int, p: int, D: int, split: bool = True, gens: str = 'all', log=None):
        N = 1 << k
        if split:
            assert all(g.get((e2, e1), 0) == v for (e1, e2), v in g.items()), 'generator not swap-symmetric'
        self.g, self.k, self.p, self.D = g, k, p, D
        self.L = Layout(k, D, split)
        self.parts = [Part(P['w'], p) for P in self.L.parts]
        self.processed = [set() for _ in self.parts]
        self.log = log or (lambda **kw: None)
        self.t0 = time.time()
        self.rounds = 0
        gdeg = max(bin(e1).count('1') + bin(e2).count('1') for e1, e2 in g)
        # generator rows: NF(m g) for squarefree m = (b1, b2) of weight <= D - deg g (all in V_D)
        if gens == 'all':
            ms = [(b1, b2) for b1, b2 in zip(self.L.a1, self.L.a2) if bin(b1).count('1') + bin(b2).count('1') <= D - gdeg]
        else:
            ms = [(0, 0)]
        self.generator_rows = len(ms)
        terms = list(g.items())
        for a in range(0, len(ms), 256):
            chunk = ms[a:a + 256]
            H = np.zeros((len(chunk), self.L.c), dtype=F64)
            for r, (b1, b2) in enumerate(chunk):
                b1s = np.array([(b1 + e1) % N for (e1, e2), _ in terms])
                b2s = np.array([(b2 + e2) % N for (e1, e2), _ in terms])
                cols = self.L.index(b1s, b2s)
                for col, (_, v) in zip(cols, terms):
                    H[r, col] = (H[r, col] + v) % p
            for part, Y in zip(self.parts, self.L.project(H, p)):
                part.add(Y)

    def _shift_rows_of(self, pi: int, idx) -> list:
        """all shifts of rows idx of part pi, projected; returns per target part a matrix."""
        F = self.L.lift(pi, self.parts[pi].rows(idx))
        assert not np.any(F[:, :self.L.c - self.L.lowc]), 'low row with high support'
        outs = [[] for _ in self.parts]
        for ij in self.L.shift_list:
            for q, Y in enumerate(self.L.project(self.L.shifted(F, ij), self.p)):
                outs[q].append(Y)
        return [np.vstack(o) for o in outs]

    def low_rows(self, pi: int) -> np.ndarray:
        P = self.parts[pi]
        return np.flatnonzero(P.w[P.piv] < self.D)

    def run(self, batch: int = 48, deadline: float | None = None) -> bool:
        while True:
            pending = None
            for pi, P in enumerate(self.parts):
                low = self.low_rows(pi)
                todo = [i for i in low if int(P.piv[i]) not in self.processed[pi]]
                if todo:
                    pending = (pi, todo[:batch])
                    break
            if pending is None:
                return True
            pi, idx = pending
            for i in idx:
                self.processed[pi].add(int(self.parts[pi].piv[i]))
            for q, Y in enumerate(self._shift_rows_of(pi, idx)):
                self.parts[q].add(Y)
            self.rounds += 1
            if self.rounds % 20 == 0:
                self.log(stage='round', round=self.rounds, ranks=[P.rank for P in self.parts],
                         low=[len(self.low_rows(i)) for i in range(len(self.parts))],
                         processed=[len(s) for s in self.processed], sec=round(time.time() - self.t0, 1))
            if deadline and time.time() > deadline:
                return False

    def sweep_random(self, t: int = 6, seed: int = 1) -> dict:
        """t random combinations of the low rows of each part; every shift and projection must reduce
        to zero.  A non-closed space passes with probability <= (#maps) p^-t."""
        rng = np.random.default_rng(seed)
        bad = 0
        for pi, P in enumerate(self.parts):
            low = self.low_rows(pi)
            if len(low) == 0:
                continue
            C = rng.integers(0, self.p, size=(t, len(low))).astype(F64)
            M = np.zeros((t, P.c), dtype=F64)
            for a in range(0, len(low), 2048):
                M += C[:, a:a + 2048] @ P.rows(low[a:a + 2048])
                M = np.mod(M, self.p)
            F = self.L.lift(pi, M)
            for ij in self.L.shift_list:
                for q, Y in enumerate(self.L.project(self.L.shifted(F, ij), self.p)):
                    bad += int(np.count_nonzero(np.any(self.parts[q].residual(Y) != 0, axis=1)))
        return dict(t=t, seed=seed, nonzero_residuals=bad, maps=len(self.parts) ** 2 * len(self.L.shift_list),
                    passed=bad == 0)

    def sweep_full(self, chunk: int = 48, deadline: float | None = None) -> dict:
        """every shift of every final low row must reduce to zero (deterministic)."""
        bad = done = 0
        for pi in range(len(self.parts)):
            low = self.low_rows(pi)
            for a in range(0, len(low), chunk):
                for q, Y in enumerate(self._shift_rows_of(pi, low[a:a + chunk])):
                    bad += int(np.count_nonzero(np.any(self.parts[q].residual(Y) != 0, axis=1)))
                done += len(low[a:a + chunk])
                if deadline and time.time() > deadline:
                    return dict(complete=False, rows_checked=done, nonzero_residuals=bad)
        return dict(complete=True, rows_checked=done, nonzero_residuals=bad, passed=bad == 0)

    # results
    def rank(self):
        return sum(P.rank for P in self.parts)

    def profile(self):
        pr = [P.profile(self.D) for P in self.parts]
        return [sum(x) for x in zip(*pr)], pr

    def contains_one(self) -> bool:
        P = self.parts[0]
        zero_col = int(np.flatnonzero(P.w == 0)[0])
        return zero_col in set(P.piv.tolist())

    def random_elements_full(self, t: int, seed: int) -> np.ndarray:
        """t random elements of the final space in full coordinates (for soundness checks)."""
        rng = np.random.default_rng(seed)
        out = np.zeros((t, self.L.c), dtype=F64)
        for pi, P in enumerate(self.parts):
            C = rng.integers(0, self.p, size=(t, P.rank)).astype(F64)
            M = np.zeros((t, P.c), dtype=F64)
            for a in range(0, P.rank, 2048):
                M = np.mod(M + C[:, a:a + 2048] @ P.rows(np.arange(a, min(a + 2048, P.rank))), self.p)
            out = np.mod(out + self.L.lift(pi, M), self.p)
        return out


def evaluate_full(Lay: Layout, rows: np.ndarray, pts: list, p: int, k: int) -> np.ndarray:
    """values of full-coordinate rows at points (E1, E2): monomial (a1, a2) -> h^(a1 E1 + a2 E2)."""
    h, H = subgroup(p, k)
    Hv = np.array(H, dtype=np.int64)
    N = 1 << k
    out = np.zeros((rows.shape[0], len(pts)), dtype=np.int64)
    R = rows.astype(np.int64)
    for z, (E1, E2) in enumerate(pts):
        ev = Hv[(Lay.a1 * E1 + Lay.a2 * E2) % N]
        out[:, z] = (R % p * ev % p).sum(axis=1) % p
    return out


def eval_rank(Lay: Layout, pts: list, p: int, k: int, maxw: int) -> int:
    """rank of evaluation of the monomials of weight <= maxw at the points."""
    if not pts:
        return 0
    h, H = subgroup(p, k)
    Hv = np.array(H, dtype=np.int64)
    N = 1 << k
    cols = np.flatnonzero(Lay.w <= maxw)
    M = np.zeros((len(pts), len(cols)), dtype=F64)
    for z, (E1, E2) in enumerate(pts):
        M[z] = Hv[(Lay.a1[cols] * E1 + Lay.a2[cols] * E2) % N]
    P = Part(np.zeros(len(cols), dtype=np.int64), p)
    P.add(M)
    return P.rank
