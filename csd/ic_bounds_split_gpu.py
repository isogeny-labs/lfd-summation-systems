"""Exact closures V_D of the split presentation with m summands on a device (numpy / CuPy).

Same object as csd.ic_bounds_split_m.closure (system `split_generators_m`, monomial space
`SplitSpace`, V_D modulo the chains), computed with a different engine:

* **compact reduced echelon form.**  A basis row is e_pivot + (entries on the free columns); only
  the free columns are stored, so the basis is a (rank x (n - rank)) matrix instead of the
  (rank x n) matrix of `InPlaceClosure`: at most n^2/4 cells instead of n^2, and one third of the
  multiply-adds.  The update that keeps the form reduced is done in place on the rows that meet a
  new pivot; columns that became pivots are dropped tile by tile once they exceed `compact_frac`
  of the stored columns, so there is never a second full copy (one tile of temporary).
* **summand-swap split.**  sigma exchanges the first two summands (blocks 0 and 1 of the mask).
  g_1 = S_3(x_1, x_2, z_1) is symmetric in x_1, x_2 and no other input contains them, so V_D is
  sigma-stable and V_D = V_D^+ (+) V_D^-.  The engine works in orbit coordinates (+ part: one
  column per orbit; - part: one per orbit of size two), which halves the dimension.  This is the
  whole symmetry group of the system among permutations of summands and auxiliaries (measured,
  experiments/gpu-split-m/README.md section 3); it is checked on the tables and on the inputs at
  run time and refused (one part, no split) when a check fails.
* **device.**  All matrices live behind the `xp` interface of csd.ic_bounds_big (`cpu` = numpy,
  `gpu` = CuPy, `fake-gpu` = numpy arrays that raise when a host array is mixed in).  Arithmetic is
  exact float64: residues in [0, p), p < 2^24, inner dimensions cut at `max_inner(p)`.
* **batch echelon with host pivoting.**  The rows of a batch are echelonised in blocks of `sub`
  rows; the device reports the first nonzero column of every row, the host eliminates the small
  (rows x leading columns) matrix and returns the transformation, the device applies it with one
  product.  About ten device calls per block, not per row.

Pivot rule.  Columns are sorted by degree, descending.  A pivot of a row lies in the highest
degree block where the row is nonzero (not necessarily its first nonzero column).  That is all the
closure argument uses: a row with a pivot of degree d is supported on degree <= d, and since the
form is reduced, an element of the span of degree <= d is a combination of the rows with pivot
degree <= d.  So `low_profile` (pivots of degree <= d) equals dim(V_D n M_<=d) for any such choice.

Products.  Every row of pivot degree < D is multiplied by every operator once, in its state at
that time.  That is enough: if b(t) was the row when it was processed and b the final row, then
b(t) = b + sum c_j b_j over rows b_j inserted later (reduced form), all of degree < D, so the
products of b follow from those of b(t) and of later rows by induction from the last row inserted.
`verify_closed` does not rely on this: it multiplies every final low row by every operator and
requires the reduction to vanish, and requires every input in the span.

Operators in the split.  For a variable v fixed by sigma (blocks >= 2, auxiliaries), v maps each
part to itself.  For u = u_{0,j}: the images of a part-p element w under u_{0,j} and u_{1,j} span
the same space as pi_+ (u w) and pi_- (u w), so the operators are (j, p -> q) for the four pairs;
u_{1,j} is not used.  All operators are gathers with coefficients +-1 in orbit coordinates.

Nothing here changes csd.ic_bounds_split_m or csd.ic_bounds_big.
"""
from __future__ import annotations

import hashlib
import json
import os
import time

import numpy as np

from . import ic_bounds_big as B
from .ic_bounds_split_m import SplitSpace, generator_degree

ENGINE = 'ic_bounds_split_gpu/1'
F64 = np.float64


def _sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


ENGINE_SHA256 = _sha(__file__)


# =============================================================================== symmetry tables

def sigma_permutation(sp: SplitSpace) -> np.ndarray:
    """Index permutation of the monomials under the exchange of blocks 0 and 1."""
    k = sp.k
    N = (1 << k) - 1
    out = np.empty(len(sp), dtype=np.int64)
    for i, (mask, es) in enumerate(sp.mons):
        b0, b1 = mask & N, (mask >> k) & N
        out[i] = sp.index[((mask & ~(N | (N << k))) | b1 | (b0 << k), es)]
    return out


def _swap_generator(g, k):
    N = (1 << k) - 1
    out = {}
    for (mask, es), c in g.items():
        b0, b1 = mask & N, (mask >> k) & N
        out[((mask & ~(N | (N << k))) | b1 | (b0 << k), es)] = c
    return out


def check_symmetry(sp: SplitSpace, gens, sigma: np.ndarray):
    """Exact run-time checks that sigma is a symmetry of the closure problem.  Returns (ok, reason)."""
    n, k = len(sp), sp.k
    if not np.array_equal(sigma[sigma], np.arange(n)):
        return False, 'sigma is not an involution'
    if not np.array_equal(sp.deg[sigma], sp.deg):
        return False, 'sigma changes degrees'
    cols = sp.cols
    for v in range(len(cols)):
        w = v + k if v < k else (v - k if v < 2 * k else v)       # sigma u_{0,j} sigma = u_{1,j}
        src = np.flatnonzero(cols[v] >= 0)
        if not np.array_equal(np.flatnonzero(cols[w] >= 0), np.sort(sigma[src])):
            return False, f'shift {v}: domain not carried to the domain of shift {w}'
        if not np.array_equal(sigma[cols[v, src]], cols[w, sigma[src]]):
            return False, f'shift {v} does not intertwine with shift {w}'
    for i, g in enumerate(gens):
        if _swap_generator(g, k) != g:
            return False, f'input {i} is not invariant under the exchange of the first two summands'
    return True, 'checked: involution, degrees, all shifts, all inputs'


class Coordinates:
    """Parts and operators.  split=False: one part (all monomials), one operator per variable.
    split=True: orbit coordinates of sigma, parts + and -, operators as in the module docstring."""

    def __init__(self, sp: SplitSpace, gens, split: bool = True):
        self.sp = sp
        n = len(sp)
        self.symmetry = 'none (not requested)'
        sigma = None
        if split and sp.m >= 2:
            sigma = sigma_permutation(sp)
            ok, why = check_symmetry(sp, gens, sigma)
            self.symmetry = ('split: ' if ok else 'refused: ') + why
            if not ok:
                sigma = None
        self.split = sigma is not None
        self.sigma = sigma
        nvar = len(sp.cols)
        inv = np.full((nvar, n), -1, dtype=np.int64)
        for v in range(nvar):
            src = np.flatnonzero(sp.cols[v] >= 0)
            dst = sp.cols[v, src]
            assert len(np.unique(dst)) == len(dst), 'a shift is not injective'
            inv[v, dst] = src
        if not self.split:
            self.parts = [dict(n=n, deg=sp.deg.copy(), rep=np.arange(n), sign=1)]
            self.ops = [dict(var=v, src=0, dst=0, terms=[(inv[v], np.where(inv[v] >= 0, 1.0, 0.0))])
                        for v in range(nvar)]
            self.orb = [np.arange(n)]
            self.sg = [np.ones(n)]
            return
        rep = np.minimum(np.arange(n), sigma)
        is_rep = rep == np.arange(n)
        selfp = sigma == np.arange(n)
        reps_plus = np.flatnonzero(is_rep)
        reps_minus = np.flatnonzero(is_rep & ~selfp)
        self.parts, self.orb, self.sg = [], [], []
        for s, reps in ((1, reps_plus), (-1, reps_minus)):
            col = np.full(n, -1, dtype=np.int64)
            col[reps] = np.arange(len(reps))
            orb = col[rep]                      # column of the orbit of x, -1 if the part has none
            sg = np.where(is_rep, 1.0, float(s))
            sg[orb < 0] = 0.0
            self.parts.append(dict(n=len(reps), deg=sp.deg[reps].copy(), rep=reps, sign=s))
            self.orb.append(orb)
            self.sg.append(sg)
        assert self.parts[0]['rep'][-1] == n - 1, 'the constant monomial must be the last + column'
        k = sp.k

        def term(p, y, extra=1.0):
            """gather term: coordinates of part p read at the full monomials y (-1: nothing)."""
            ok = y >= 0
            yy = np.where(ok, y, 0)
            idx = np.where(ok, self.orb[p][yy], -1)
            return idx, np.where(idx >= 0, extra * self.sg[p][yy], 0.0)

        self.ops = []
        for v in range(nvar):
            if k <= v < 2 * k:
                continue                         # u_{1,j}: covered by the two projections of u_{0,j}
            for p in (0, 1):
                if v < k:
                    for q in (0, 1):
                        rq = self.parts[q]['rep']
                        t1 = term(p, inv[v][rq])
                        t2 = term(p, inv[v][sigma[rq]], float(self.parts[q]['sign']))
                        self.ops.append(dict(var=v, src=p, dst=q, terms=[t1, t2]))
                else:
                    self.ops.append(dict(var=v, src=p, dst=p, terms=[term(p, inv[v][self.parts[p]['rep']])]))
        self.ops = [o for o in self.ops if any(np.any(idx >= 0) for idx, _ in o['terms'])]

    def project(self, v: np.ndarray, q: int, p: int) -> np.ndarray:
        """Coordinates in part q of (1 + s_q sigma) v for a full vector v (host, integers mod p)."""
        if not self.split:
            return v % p
        rq = self.parts[q]['rep']
        return (v[rq] + self.parts[q]['sign'] * v[self.sigma[rq]]) % p

    def lift(self, rows: np.ndarray, q: int, p: int) -> np.ndarray:
        """Full vectors of part-q coordinate rows (host)."""
        if not self.split:
            return rows % p
        ok = self.orb[q] >= 0
        out = np.zeros((rows.shape[0], len(self.sp)), dtype=np.int64)
        out[:, ok] = rows[:, self.orb[q][ok]] * self.sg[q][ok].astype(np.int64)
        return out % p


# =============================================================================== small host kernel

def host_rref_transform(A: np.ndarray, p: int):
    """A: small host matrix (m x c), integers.  Returns (T, r, pc): T is m x m invertible mod p,
    the first r rows of T A are the reduced echelon form of A with pivot columns pc (indices into
    the columns of A) and the remaining rows of T A are zero."""
    m, c = A.shape
    M = np.concatenate([np.asarray(A, dtype=np.int64) % p, np.eye(m, dtype=np.int64)], axis=1)
    r, pc = 0, []
    for j in range(c):
        if r == m:
            break
        nz = np.flatnonzero(M[r:, j])
        if not len(nz):
            continue
        i = r + int(nz[0])
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = M[r] * pow(int(M[r, j]), -1, p) % p
        col = M[:, j].copy()
        col[r] = 0
        nzr = np.flatnonzero(col)
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        pc.append(j)
        r += 1
    return M[:, c:], r, pc


def host_inverse(M: np.ndarray, p: int) -> np.ndarray:
    T, r, pc = host_rref_transform(M, p)
    assert r == M.shape[0] and pc == list(range(r)), 'singular block in the batch echelon'
    return T


def rref_batch(R, fdeg: np.ndarray, p: int, mi: int, sub: int, stats: dict, width: int = 4):
    """Reduced echelon form of the rows of R (b x f device array, entries in [0, p)) under the pivot
    rule of the module docstring.  fdeg (host): degree of each column, non-increasing.
    Returns (L, lead): L is t x f on the device with the identity on the columns `lead` (host)."""
    xp, DV = B.xp, B.DV
    b, f = R.shape
    assert sub <= mi
    L = None
    t = 0
    lead = np.zeros(0, dtype=np.int64)
    for a in range(0, b, sub):
        S = xp.array(R[a:a + sub])
        if t:
            B.matmul_sub(S, S[:, DV.to_dev(lead)], L[:t], p, mi)
        rows, pn = [], []
        while True:
            nzm = S != 0
            alive = np.flatnonzero(DV.to_host(xp.any(nzm, axis=1)))
            if not len(alive):
                break
            if len(alive) < S.shape[0]:
                ad = DV.to_dev(alive)
                S, nzm = S[ad], nzm[ad]
            first = DV.to_host(xp.argmax(nzm, axis=1)).astype(np.int64)
            d = fdeg[first]
            top = d.max()
            sel = np.flatnonzero(d == top)
            sd = DV.to_dev(sel)
            Ssel = S[sd]
            # candidate pivot columns: the leading columns of the selected rows and, up to `width` per
            # row, the first columns of the same degree block where a selected row is nonzero
            lo = int(np.searchsorted(-fdeg, -top, side='left'))
            hi = int(np.searchsorted(-fdeg, -top, side='right'))
            occupied = lo + np.flatnonzero(DV.to_host(xp.any(nzm[sd][:, lo:hi], axis=0)))
            J = np.union1d(first[sel], occupied[:width * len(sel)])
            A = DV.to_host(Ssel[:, DV.to_dev(J)]).astype(np.int64)
            T, r, pc = host_rref_transform(A, p)
            stats['rounds'] += 1
            stats['host_cells'] += A.size
            Ssel = DV.to_dev(T.astype(F64)) @ Ssel          # inner dimension <= sub <= mi
            B.modp_(Ssel, p)
            rows.append(Ssel[:r])
            pn += [int(J[c]) for c in pc]
            others = np.setdiff1d(np.arange(S.shape[0]), sel)
            rest = [Ssel[r:]] if r < len(sel) else []
            if len(others):
                rest.append(S[DV.to_dev(others)])
            if not rest:
                break
            S = xp.concatenate(rest) if len(rest) > 1 else rest[0]
        if not pn:
            continue
        Ln = xp.concatenate(rows) if len(rows) > 1 else rows[0]
        pnh = np.asarray(pn, dtype=np.int64)
        tn = len(pn)
        if tn > 1:                                           # reduce the new rows among themselves
            M = DV.to_host(Ln[:, DV.to_dev(pnh)]).astype(np.int64)
            if not np.array_equal(M, np.eye(tn, dtype=np.int64)):
                Ln = DV.to_dev(host_inverse(M, p).astype(F64)) @ Ln
                B.modp_(Ln, p)
        if L is None:
            L = xp.zeros((b, f), dtype=F64)
        if t:
            X = L[:t][:, DV.to_dev(pnh)]
            Lt = L[:t]
            B.matmul_sub(Lt, X, Ln, p, mi)
        L[t:t + tn] = Ln
        t += tn
        lead = np.concatenate([lead, pnh])
        stats['madds_rref'] += 2 * tn * t * f
    if L is None:
        return None, lead
    return L[:t], lead


# =============================================================================== one part

class PartBasis:
    """Reduced echelon basis of one part: tiles of (rows x free columns) on the device."""

    def __init__(self, n: int, deg: np.ndarray, p: int, D: int, eng):
        self.n, self.deg, self.p, self.D, self.E = n, deg, p, D, eng
        self.mi = B.max_inner(p)
        self.tiles = []                                 # device arrays, tile_rows x nphys, rows filled in order
        self.starts = [0]                               # starts[i]: first row of tile i; starts[-1] = rank
        self.piv = np.zeros(0, dtype=np.int64)
        self.phys = np.arange(n, dtype=np.int64)        # column of the part behind each stored column
        self.alive = np.ones(n, dtype=bool)             # False: the column has become a pivot (dead)
        self.compactions = 0
        self.low = np.zeros(0, dtype=np.int64)          # rows with pivot degree < D, insertion order
        self._dev = {}

    # -- bookkeeping
    @property
    def rank(self):
        return len(self.piv)

    @property
    def nphys(self):
        return len(self.phys)

    @property
    def nfree(self):
        return self.n - self.rank

    def dev(self, name):
        if name not in self._dev:
            self._dev[name] = B.DV.to_dev(getattr(self, name))
        return self._dev[name]

    def cells(self):
        return self.rank * self.nphys

    def tile(self, ti):
        """The filled rows of tile ti (a view)."""
        return self.tiles[ti][:self.starts[ti + 1] - self.starts[ti]]

    def filled(self):
        return [self.tile(ti) for ti in range(len(self.tiles))]

    def _gather(self, rows: np.ndarray):
        """rows (host, ascending) of the stored matrix: len(rows) x nfree device array (a copy)."""
        xp, DV = B.xp, B.DV
        out = []
        starts = np.asarray(self.starts)
        ti = np.searchsorted(starts, rows, side='right') - 1
        for t in np.unique(ti):
            loc = rows[ti == t] - starts[t]
            out.append(self.tiles[t][DV.to_dev(loc)])       # loc < filled rows
        return xp.concatenate(out) if len(out) > 1 else out[0]

    def expand(self, rows: np.ndarray):
        """Basis rows as dense vectors of the part: len(rows) x n device array."""
        xp, DV = B.xp, B.DV
        X = xp.zeros((len(rows), self.n), dtype=F64)
        if self.nphys:
            X[:, self.dev('phys')] = self._gather(rows)
        X[DV.to_dev(np.arange(len(rows))), DV.to_dev(self.piv[rows])] = 1.0
        return X

    # -- arithmetic
    def reduce(self, P):
        """P: b x n device rows with entries in [0, p).  Returns the residues on the stored columns
        (exactly zero on the dead ones: a stored row has the identity on its own dead pivot column)."""
        xp, DV = B.xp, B.DV
        st = self.E.stats
        R = P[:, self.dev('phys')]
        if not self.rank:
            return R
        C = P[:, self.dev('piv')]
        b = P.shape[0]
        pending = 0
        for ti, tile in enumerate(self.filled()):
            r0, r1 = self.starts[ti], self.starts[ti + 1]
            Ct = C[:, r0:r1]
            hit = np.flatnonzero(DV.to_host(xp.any(Ct != 0, axis=0)))
            if not len(hit):
                continue
            if len(hit) < self.E.hit_frac * (r1 - r0):
                hd = DV.to_dev(hit)
                Ct, Tt = Ct[:, hd], tile[hd]
            else:
                Tt = tile
            h = Ct.shape[1]
            for a in range(0, h, self.mi):
                c = min(self.mi, h - a)
                if pending + c > self.mi:
                    B.modp_(R, self.p)
                    pending = 0
                R -= Ct[:, a:a + c] @ Tt[a:a + c]
                pending += c
            st['madds_reduce'] += b * h * self.nphys
            st['madds_reduce_dense'] += b * (r1 - r0) * self.nphys
        if pending:
            B.modp_(R, self.p)
        return R

    def absorb(self, R) -> int:
        """R: rows reduced against the basis (on the free columns).  Echelonise, update, append."""
        xp, DV = B.xp, B.DV
        st = self.E.stats
        t0 = time.perf_counter()
        L, lead = rref_batch(R, self.deg[self.phys], self.p, self.mi, self.E.sub, st)
        DV.sync()
        st['s_rref'] += time.perf_counter() - t0
        t = len(lead)
        if not t:
            return 0
        assert self.alive[lead].all(), 'a dead column became a pivot'
        t0 = time.perf_counter()
        ld = DV.to_dev(lead)
        for ti in range(len(self.tiles)):                 # old rows <- row - row[lead] @ L, in place
            tile = self.tile(ti)
            X = tile[:, ld]
            touched = np.flatnonzero(DV.to_host(xp.any(X != 0, axis=1)))
            if len(touched):
                if len(touched) < self.E.hit_frac * X.shape[0]:
                    td = DV.to_dev(touched)
                    part = tile[td]
                    B.matmul_sub(part, X[td], L, self.p, self.mi)
                    tile[td] = part
                else:
                    B.matmul_sub(tile, X, L, self.p, self.mi)
                st['madds_update'] += len(touched) * t * self.nphys
            st['madds_update_dense'] += X.shape[0] * t * self.nphys
        a = 0
        while a < t:                                      # append into the last tile, new tiles as needed
            used = self.starts[-1] - self.starts[-2] if self.tiles else self.E.tile_rows
            if used == self.E.tile_rows:
                self.tiles.append(xp.zeros((self.E.tile_rows, self.nphys), dtype=F64))
                self.starts.append(self.starts[-1])
                used = 0
            c = min(t - a, self.E.tile_rows - used)
            self.tiles[-1][used:used + c] = L[a:a + c]
            self.starts[-1] += c
            a += c
        del L
        r0 = self.rank
        newpiv = self.phys[lead]
        self.piv = np.concatenate([self.piv, newpiv])
        self.alive[lead] = False
        lowrows = r0 + np.flatnonzero(self.deg[newpiv] < self.D)
        self.low = np.concatenate([self.low, lowrows])
        self._dev.pop('piv', None)
        self.E.note_memory()
        if self.nphys - int(self.alive.sum()) > self.E.compact_frac * self.nphys:
            self.compact()
        DV.sync()
        st['s_update'] += time.perf_counter() - t0
        return t

    def compact(self):
        """Drop the dead columns, tile by tile (one tile of temporary)."""
        keep = np.flatnonzero(self.alive)
        kd = B.DV.to_dev(keep)
        for ti in range(len(self.tiles)):
            fill = self.starts[ti + 1] - self.starts[ti]
            if fill == self.E.tile_rows:
                self.tiles[ti] = self.tiles[ti][:, kd]
            else:
                new = B.xp.zeros((self.E.tile_rows, len(keep)), dtype=F64)
                new[:fill] = self.tiles[ti][:fill][:, kd]
                self.tiles[ti] = new
        self.phys = self.phys[keep]
        self.alive = np.ones(len(keep), dtype=bool)
        self._dev = {}
        self.compactions += 1
        B.DV.free_pool()

    # -- queries
    def profile(self, D: int):
        dp = self.deg[self.piv]
        return [int(np.sum(dp <= d)) for d in range(D + 1)]

    def host_rows(self) -> np.ndarray:
        """All basis rows as dense host integer vectors of the part (small cases only)."""
        out = np.zeros((self.rank, self.n), dtype=np.int64)
        if self.rank:
            if self.nphys:
                M = np.concatenate([B.DV.to_host(t) for t in self.filled()])
                out[:, self.phys] = M.astype(np.int64)
            out[np.arange(self.rank), self.piv] = 1
        return out


# =============================================================================== engine

class Stop(Exception):
    pass


class SplitClosure:
    """V_D of the split presentation.  gens: output of split_generators_m (normal forms)."""

    def __init__(self, gens, m, k, p, D, *, device='cpu', split=True, batch=1024, tile_rows=4096, sub=64,
                 hit_frac=0.7, compact_frac=0.125, fault=None, log=None, budget_s=None, mem_budget_gb=None, checkpoint=None,
                 checkpoint_every_s=None):
        if p >= 1 << 24 or p % 2 == 0:
            raise ValueError('p must be an odd prime below 2^24 (exact float64 reduction)')
        self.m, self.k, self.p, self.D = m, k, p, D
        self.device = device
        B.use_device(device)
        if device == 'fake-gpu':
            B._Live.peak = B._Live.cur
        self.batch, self.tile_rows, self.hit_frac, self.compact_frac = batch, tile_rows, hit_frac, compact_frac
        self.sub = min(sub, B.max_inner(p))
        self.fault, self.log = fault, log
        self.budget_s, self.mem_budget = budget_s, (mem_budget_gb * 2 ** 30 if mem_budget_gb else None)
        self.checkpoint, self.checkpoint_every_s = checkpoint, checkpoint_every_s
        t0 = time.perf_counter()
        self.sp = SplitSpace(m, k, D, m - 2)
        self.n = len(self.sp)
        self.gens = [g for g in gens]
        self.co = Coordinates(self.sp, gens, split)
        self.parts = [PartBasis(pt['n'], pt['deg'], p, D, self) for pt in self.co.parts]
        self.ops = self.co.ops
        self.all_ops = self.ops
        if fault == 'drop-operator':      # test hook: the closure loop skips one operator
            self.ops = self.ops[1:]
        self.cursor = [0] * len(self.ops)
        self._opdev = {}
        self.contains_one = False
        self.stats = dict(products=0, zero_rows=0, pivots=0, batches=0, rounds=0, host_cells=0,
                          madds_reduce=0, madds_reduce_dense=0, madds_update=0, madds_update_dense=0,
                          madds_rref=0, s_setup=time.perf_counter() - t0, s_products=0.0, s_reduce=0.0,
                          s_rref=0.0, s_update=0.0, s_verify=0.0, s_closure=0.0, peak_cells=0,
                          peak_device_bytes=0)
        self.verified_closed = None
        self.status = 'new'
        self._t_start = time.perf_counter()
        self._t_ckpt = time.perf_counter()
        self._elapsed_before = 0.0

    # -- helpers
    def say(self, msg):
        if self.log:
            self.log(msg)

    def note_memory(self):
        st = self.stats
        st['peak_cells'] = max(st['peak_cells'], sum(pt.cells() for pt in self.parts))
        used = (B.DV.pool_bytes() if self.device != 'cpu'
                else 8 * sum(sum(t.size for t in pt.tiles) for pt in self.parts))
        st['peak_device_bytes'] = max(st['peak_device_bytes'], int(used))
        if self.mem_budget and 8 * sum(pt.cells() for pt in self.parts) > self.mem_budget:
            raise Stop(f'memory budget: basis {8 * sum(pt.cells() for pt in self.parts) / 2 ** 30:.2f} GiB')

    def seconds(self):
        return self._elapsed_before + time.perf_counter() - self._t_start

    def _op_terms(self, oi):
        if oi not in self._opdev:
            out = []
            for ti, (idx, coef) in enumerate(self.ops[oi]['terms']):
                if self.fault == 'wrong-gather' and oi == 0 and ti == 0:   # test hook
                    idx = np.roll(idx, 1)
                out.append((B.DV.to_dev(np.where(idx >= 0, idx, 0)), B.DV.to_dev(coef.astype(F64))))
            self._opdev[oi] = out
        return self._opdev[oi]

    def image(self, oi: int, rows: np.ndarray):
        """Operator oi applied to basis rows of its source part: len(rows) x n_dst, entries in [0, p)."""
        op = self.ops[oi]
        X = self.parts[op['src']].expand(rows)
        P = None
        for idx, coef in self._op_terms(oi):
            T = X[:, idx]
            T *= coef
            P = T if P is None else P + T
        B.modp_(P, self.p)
        return P

    def _insert(self, q: int, P):
        """Reduce the rows P (part q, entries in [0,p)) and absorb what is left."""
        xp, DV = B.xp, B.DV
        st = self.stats
        part = self.parts[q]
        t0 = time.perf_counter()
        R = part.reduce(P)
        nz = np.flatnonzero(DV.to_host(xp.any(R != 0, axis=1)))
        DV.sync()
        st['s_reduce'] += time.perf_counter() - t0
        st['zero_rows'] += P.shape[0] - len(nz)
        if not len(nz):
            return 0
        if len(nz) < R.shape[0]:
            R = R[DV.to_dev(nz)]
        t = part.absorb(R)
        st['pivots'] += t
        if q == 0 and t and part.piv[-t:].max() == part.n - 1:
            self.contains_one = True
        return t

    def add_inputs(self):
        self.input_rows = []
        for g in self.gens:
            if generator_degree(g) <= self.D:
                self.input_rows.append(self.sp.row(g, self.p))
        if not self.input_rows:
            return
        full = np.asarray(self.input_rows, dtype=np.int64)
        for q in range(len(self.parts)):
            proj = np.stack([self.co.project(v, q, self.p) for v in full])
            if q == 1:
                assert not proj.any(), 'an input has a - part although the split was accepted'
                continue
            self._insert(q, B.DV.to_dev(proj.astype(F64)))

    # -- closure
    def pending(self):
        return [len(self.parts[op['src']].low) - self.cursor[oi] for oi, op in enumerate(self.ops)]

    def saturate(self):
        st = self.stats
        t_run = time.perf_counter()
        self.status = 'running'
        try:
            while not self.contains_one:
                pend = self.pending()
                if not any(pend):
                    break
                for oi, op in enumerate(self.ops):
                    while not self.contains_one:
                        src = self.parts[op['src']]
                        todo = len(src.low) - self.cursor[oi]
                        if todo <= 0:
                            break
                        take = min(todo, self.batch)
                        rows = src.low[self.cursor[oi]:self.cursor[oi] + take]
                        t0 = time.perf_counter()
                        P = self.image(oi, rows)
                        B.DV.sync()
                        st['s_products'] += time.perf_counter() - t0
                        self.cursor[oi] += take
                        st['products'] += take
                        st['batches'] += 1
                        self._insert(op['dst'], P)
                        del P
                        if st['batches'] % 50 == 0:
                            self.say(f'  D={self.D} batches={st["batches"]} rank={[pt.rank for pt in self.parts]} '
                                     f'pending={sum(self.pending())} t={self.seconds():.1f}s '
                                     f'dev={st["peak_device_bytes"] / 2 ** 30:.2f}GiB')
                        if self.budget_s and self.seconds() > self.budget_s:
                            raise Stop(f'time budget {self.budget_s}s')
                        if (self.checkpoint and self.checkpoint_every_s
                                and time.perf_counter() - self._t_ckpt > self.checkpoint_every_s):
                            self.save_checkpoint()
                        if self.one_per_op:
                            break
                    if self.contains_one:
                        break
            self.status = 'one' if self.contains_one else 'closed'
        except Stop as e:
            self.status = f'stopped: {e}'
            if self.checkpoint:
                self.save_checkpoint()
        st['s_closure'] += time.perf_counter() - t_run
        return self.status

    one_per_op = True      # serve one batch per operator in turn (False: drain an operator first)

    def verify_closed(self) -> bool:
        """Every operator image of every final low row reduces to zero, and every input is in the span."""
        xp, DV = B.xp, B.DV
        if self.contains_one:
            return True
        t0 = time.perf_counter()
        ok = True
        saved = self.ops, self._opdev, self.fault
        self.ops, self._opdev, self.fault = self.all_ops, {}, None      # never the faulted operators
        for oi, op in enumerate(self.ops):
            src, dst = self.parts[op['src']], self.parts[op['dst']]
            for a in range(0, len(src.low), self.batch):
                R = dst.reduce(self.image(oi, src.low[a:a + self.batch]))
                if bool(DV.to_host(xp.any(R != 0))):
                    ok = False
        self.ops, self._opdev, self.fault = saved
        if self.input_rows:
            full = np.asarray(self.input_rows, dtype=np.int64)
            for q in range(len(self.parts)):
                proj = np.stack([self.co.project(v, q, self.p) for v in full])
                R = self.parts[q].reduce(DV.to_dev(proj.astype(F64)))
                if bool(DV.to_host(xp.any(R != 0))):
                    ok = False
        self.stats['s_verify'] += time.perf_counter() - t0
        return ok

    # -- results
    @property
    def rank(self):
        return self.n if self.contains_one else sum(pt.rank for pt in self.parts)

    def low_profile(self):
        if self.contains_one:
            return [int(np.sum(self.sp.deg <= d)) for d in range(self.D + 1)]
        prof = [pt.profile(self.D) for pt in self.parts]
        return [sum(x[d] for x in prof) for d in range(self.D + 1)]

    def ops_record(self):
        st = self.stats
        return dict(schema='optlevel.ops/v1', kernel=ENGINE, unit='float64 multiply-add (exact, mod p)',
                    counted=dict(reduce=st['madds_reduce'], update=st['madds_update'], rref=st['madds_rref']),
                    total=st['madds_reduce'] + st['madds_update'] + st['madds_rref'],
                    dense_equivalent=st['madds_reduce_dense'] + st['madds_update_dense'] + st['madds_rref'],
                    seconds=dict(reduce=round(st['s_reduce'], 3), update=round(st['s_update'], 3),
                                 rref=round(st['s_rref'], 3)),
                    note='reduce and update counts are the GEMM sizes actually issued (rows without a '
                         'nonzero coefficient are skipped); rref is an upper estimate')

    def record(self):
        st = dict(self.stats)
        for key in list(st):
            if key.startswith('s_'):
                st[key] = round(st[key], 3)
        return dict(engine=ENGINE, engine_sha256=ENGINE_SHA256, device=self.device, m=self.m, k=self.k, p=self.p,
                    D=self.D, monomials_le_D=self.n, symmetry=self.co.symmetry,
                    part_columns=[pt.n for pt in self.parts], part_rank=[pt.rank for pt in self.parts],
                    operators=len(self.ops), status=self.status, dim_V_D=int(self.rank),
                    low_profile=self.low_profile(), contains_one=bool(self.contains_one),
                    verified_closed=self.verified_closed, batch=self.batch, tile_rows=self.tile_rows, sub=self.sub,
                    guard_hits=int(B.GUARD_HITS[0]), seconds=round(self.seconds(), 3), stats=st,
                    products_per_pivot=round(st['products'] / max(1, st['pivots']), 3),
                    peak_basis_gib=round(8 * st['peak_cells'] / 2 ** 30, 4),
                    peak_device_gib=round(st['peak_device_bytes'] / 2 ** 30, 4),
                    maxrss_mb=B.maxrss_mb(), loadavg=os.getloadavg()[0], ops=self.ops_record())

    # -- checkpoints
    def save_checkpoint(self):
        path = self.checkpoint
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        data = dict(meta=np.array(json.dumps(dict(
            engine=ENGINE, m=self.m, k=self.k, p=self.p, D=self.D, split=self.co.split, cursor=self.cursor,
            contains_one=self.contains_one, stats=self.stats, elapsed=self.seconds(),
            ops=len(self.ops), gens=hashlib.sha256(repr(sorted(map(lambda g: sorted(g.items()), self.gens))).encode()).hexdigest()))))
        for q, pt in enumerate(self.parts):
            data[f'piv{q}'], data[f'phys{q}'], data[f'alive{q}'], data[f'low{q}'] = pt.piv, pt.phys, pt.alive, pt.low
            data[f'M{q}'] = (np.concatenate([B.DV.to_host(t) for t in pt.filled()]) if pt.tiles
                             else np.zeros((0, pt.nphys)))
        tmp = path + '.tmp.npz'
        np.savez(tmp, **data)
        os.replace(tmp, path)
        self._t_ckpt = time.perf_counter()
        self.say(f'  checkpoint {path} rank={[pt.rank for pt in self.parts]} t={self.seconds():.1f}s')

    def load_checkpoint(self) -> bool:
        path = self.checkpoint
        if not path or not os.path.exists(path):
            return False
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z['meta']))
        mine = hashlib.sha256(repr(sorted(map(lambda g: sorted(g.items()), self.gens))).encode()).hexdigest()
        want = dict(engine=ENGINE, m=self.m, k=self.k, p=self.p, D=self.D, split=self.co.split, ops=len(self.ops), gens=mine)
        for key, v in want.items():
            if meta[key] != v:
                raise ValueError(f'checkpoint {path} is for another run: {key} = {meta[key]!r}, expected {v!r}')
        for q, pt in enumerate(self.parts):
            pt.piv, pt.phys, pt.alive, pt.low = z[f'piv{q}'], z[f'phys{q}'], z[f'alive{q}'], z[f'low{q}']
            M = z[f'M{q}']
            pt.tiles, pt.starts, pt._dev = [], [0], {}
            for a in range(0, M.shape[0], self.tile_rows):
                blk = np.zeros((self.tile_rows, M.shape[1]))
                rows = M[a:a + self.tile_rows]
                blk[:len(rows)] = rows
                pt.tiles.append(B.DV.to_dev(blk))
                pt.starts.append(pt.starts[-1] + len(rows))
        self.cursor = list(meta['cursor'])
        self.contains_one = meta['contains_one']
        self.stats.update(meta['stats'])
        self._elapsed_before = meta['elapsed']
        self._t_start = time.perf_counter()
        self.resumed = True
        return True


def closure(gens, m, k, p, D, *, device='cpu', split=True, verify=True, resume=False, **kw):
    """Exact V_D; returns (record, engine).  The record has the fields of
    csd.ic_bounds_split_m.closure (D, monomials_le_D, dim_V_D, low_profile, contains_one,
    verified_closed) and the engine's own."""
    eng = SplitClosure(gens, m, k, p, D, device=device, split=split, **kw)
    if not (resume and eng.load_checkpoint()):
        eng.add_inputs()
    else:
        eng.input_rows = [eng.sp.row(g, p) for g in eng.gens if generator_degree(g) <= D]
    eng.saturate()
    if eng.status == 'closed' and verify:
        eng.verified_closed = bool(eng.verify_closed())
        if not eng.verified_closed:
            eng.status = 'verification-failed'
    rec = eng.record()
    rec['inputs_used'] = len(eng.input_rows)
    return rec, eng


def predicted_cells(n_parts):
    """Upper bound of the stored cells: sum over parts of max_r r (n - r) = n^2 / 4."""
    return sum(nq * nq // 4 for nq in n_parts)
