"""Degree lift through the quotient: from a finished degree-D chain closure V_D (two summands, depth k,
swap split) decide whether 1 lies in V_{D+1}, and compute dim V_{D+1}, its profile and the canonical
reduced echelon basis of V_{D+1} n M_<=D, with linear algebra in the small space
T = F^Sp (+) (M_<=D / V_D) instead of a closure in M_<=D+1.

Needs ic_bounds_big.py (the closure engine) next to this file or in the same directory (csd/).  numpy /
cupy through the engine's device abstraction (``cpu``, ``gpu``, ``fake-gpu``).

MATHEMATICS (statement and proof: experiments/ic-bounds-k14/LIFT.md).
  L = span(M_<=D), W = V_{D+1} n L.  V_{D+1} = W + sum_u u W, and W is the limit of
      W_0 = V_D,   X_i = W_i + sum_u u W_i,   W_{i+1} = X_i n L.
  A weight-D monomial is REGULAR if its swap orbit is a pivot of both parts of V_D (self-paired: of the
  + part), else FREE-TYPE.  A weight-(D+1) monomial is SPECIAL if one of its weight-D divisors is
  free-type, else REGULAR.  For a regular weight-D monomial m there is l_m = e_m + phi(m) in V_D with
  phi(m) supported on free-type monomials, and u l_m = e_{um} + (special top part) + (part in L).
  Modulo V_D (the residual against the V_D basis) and for one target part q, every generator of X_i^q
  becomes a vector of F^{Reg_q} (+) T_q, T_q = F^{Sp_q} (+) F^{Free_q}, with at most one nonzero regular
  coordinate.  For a regular column r all generators with top r (one per divisor variable) are
  normalised and chained by differences; the differences and the generators without a regular coordinate
  span K_q in T_q.  Then X^q = V_D^q (+) span(one reductor per regular column) (+) lift(K_q), the falls
  are the elements of K_q with zero Sp part (rows of the echelon form of K_q whose pivot is a Free
  column, Sp columns ordered first), their shifts are again vectors of T, and at the fixed point
      V_{D+1}^q n L = V_D^q (+) falls_q,    dim V_{D+1}^q = dim V_D^q + #Reg_q + dim K_q.

ALGORITHM (target part q; everything is reduced modulo the current K_q as it is produced)
  * the image of a basis row b = e_c + B[c] (B on the free columns of its part p) under
    C = project_q . u_{1,j} . lift_p is  red(C e_c) + B[c] @ R_{p,j,q},  where the rows of R_{p,j,q} are
    red(C e_f), f free in part p: one matrix (#free_p x dim T_q) per (p, j), built from gathers of rows of
    the V_D^q basis, kept reduced modulo K_q on K_q's physical columns ("follower" of K_q, synchronised
    lazily, float32, LRU cache under a byte budget).  All heavy work is the product B[rows] @ follower.
  * a regular weight-D orbit o gives, per shift j, the unit (o, j): rows y+ (part 0) and y- (part 1),
    g1 = (y+ + y-)/2 = project_q(u_{1,j} l_m), g2 = (y+ - y-)/2 = project_q(u_{1,j} l_{sigma m}), each with
    a single top term.  Units are grouped by the index set J = bits(a1 | a2) | {j}: every generator of a
    regular column r has the same J, so a batch made of whole J-groups contains all D+1 generators of each
    of its regular columns (checked: the count must be D+1, or popcount(a) for a self-paired column).
  * rows whose orbit is a pivot in exactly one part (and, with ``include_low``, the low rows) have a
    special top part only and are processed directly.
  * every candidate (difference or direct vector) is reduced and added to K_q: nothing is sampled.
  * fall phase: every fall row of K_{q0} (processed once per target part, at its state then; the engine's
    closure invariant) is shifted with the same followers into K_0 and K_1; rounds until none is pending.
  * finally the falls are added to the V_D parts (engine ``Part.add_residuals``), which become the
    canonical reduced echelon bases of V_{D+1} n M_<=D (degree-D column order): hash and profile.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from collections import OrderedDict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import ic_bounds_big as B  # noqa: E402

F64, F32 = np.float64, np.float32
LIFT_VERSION = 'ic_bounds_lift/1'
LIFT_SHA256 = B.file_sha256(__file__)
EXIT_SIGTERM = 75


def _d(h):
    return B.DV.to_dev(h)


def _h(a):
    return B.DV.to_host(a)


def addmul(acc, A, Bm, p: int, mi: int, sign: int = 1):
    """acc <- (acc + sign * A @ Bm) mod p exactly; entries in [0, p); the inner dimension is split into
    chunks of <= mi (= max_inner(p)) with a reduction after each chunk (same bound as the engine)."""
    n = A.shape[1]
    for a in range(0, n, mi):
        prod = A[:, a:a + mi] @ Bm[a:a + mi]
        if sign < 0:
            acc -= prod
        else:
            acc += prod
        del prod
        B.modp_(acc, p)
    return acc


def rows_dev(P, rows):
    """basis rows of an engine Part on its physical columns: n x ncols device float64 (pivot not included)."""
    rows = np.asarray(rows, dtype=np.int64)
    out = B.xp.zeros((len(rows), P.ncols), dtype=F64)
    for ti, idx, off in P._groups(rows):
        out[_d(idx)] = P._rows(ti, off, slice(None)).astype(F64)
    return out


# =============================================================================== K part and followers

class _KEng:
    """minimal engine stand-in for the K parts (no spilling)."""
    tile_budget = None

    def __init__(self, tile_bytes, work_bytes, compact_frac, col_block_opt):
        self.tile_bytes, self.work_bytes = int(tile_bytes), int(work_bytes)
        self.compact_frac, self.col_block_opt = compact_frac, col_block_opt

    def activate(self, part):
        pass

    def ensure_room(self, part, nbytes):
        pass

    def col_block(self, n_rows, tr):
        if self.col_block_opt:
            return int(self.col_block_opt)
        return max(64, int(0.25 * self.work_bytes // (8 * (tr + n_rows))))


class Follower:
    """A matrix (n x K.ncols, float32, device) whose rows are kept reduced modulo K on K's physical columns
    (dead columns exactly zero).  ``sync`` reduces against the rows K gained since the last sync: K is in
    reduced echelon form, so a vector that vanishes on the old pivots stays so, and subtracting
    v[c] * row(c) for the new pivots c gives the reduced vector."""

    def __init__(self, K, M):
        self.K, self.M, self.rank_synced = K, M, K.rank

    @property
    def nbytes(self):
        return int(self.M.shape[0]) * int(self.M.shape[1]) * 4

    def sync(self, chunk: int = 4096):
        K = self.K
        if self.rank_synced == K.rank:
            return
        xp = B.xp
        p, mi = K.p, K.mi
        n = self.M.shape[0]
        new = np.arange(self.rank_synced, K.rank)
        for u in range(0, len(new), chunk):
            rows = new[u:u + chunk]
            lead = np.searchsorted(K.pc2col, K.piv[rows])
            assert np.all(K.pc2col[lead] == K.piv[rows]), 'follower: pivot column already compacted away'
            ld = _d(lead)
            X = self.M[:, ld].astype(F64)
            if not bool(_h(xp.any(X != 0))):
                continue
            c0 = int(lead.min())
            cw = K.E.col_block(n, K.TR)
            groups = list(K._groups(rows))
            Xs = [X[:, _d(idx)] for _, idx, _ in groups]
            for a in range(c0, K.ncols, cw):
                b = min(K.ncols, a + cw)
                blk = self.M[:, a:b].astype(F64)
                for (ti, idx, off), Xg in zip(groups, Xs):
                    addmul(blk, Xg, K._rows(ti, off, slice(a, b)).astype(F64), p, mi, -1)
                self.M[:, a:b] = blk.astype(F32)
            self.M[:, ld] = 0.0
            K.madds_follow += n * len(rows) * (K.ncols - c0)
        self.rank_synced = K.rank


class KPart(B.Part):
    """engine Part (reduced echelon basis on T columns) plus followers and a row-preserving reduction."""

    def __init__(self, c: int, p: int, eng, tile_rows=None):
        super().__init__(np.ones(c, dtype=np.int64), c, p, 0, eng, tile_rows)
        self.followers: list = []
        self.madds_follow = 0

    def reduce_keep(self, T: np.ndarray, Yt):
        """rows of Yt (n x len(T), T sorted unique part columns, entries in [0,p), device) reduced against
        the basis: n x ncols device float64 on the physical columns (all rows kept, zero rows included).
        Same computation as ``Part.residual``."""
        xp = B.xp
        n = Yt.shape[0]
        R = xp.zeros((n, self.ncols), dtype=F64)
        if Yt.shape[1] == 0 or n == 0:
            return R
        nzc = _h(xp.any(Yt != 0, axis=0))
        if not nzc.any():
            return R
        T = T[nzc]
        Yt = Yt[:, _d(np.flatnonzero(nzc))]
        rowid = self.col2row[T]
        isp = rowid >= 0
        fcol = np.flatnonzero(~isp)
        if len(fcol):
            fp = self.col2phys[T[fcol]]
            assert np.all(fp >= 0)
            R[:, _d(fp)] = Yt[:, _d(fcol)]
        sel = rowid[isp]
        if len(sel):
            Ysel = Yt[:, _d(np.flatnonzero(isp))]
            c0 = int(np.searchsorted(self.pc2col, int(self.piv[sel].min())))
            self.accumulate(sel, Ysel, R, c0, -1)
            self.madds['reduce'] += n * len(sel) * (self.ncols - c0)
        return R

    def compact(self, force=False):
        dead = self.ncols - int(self.alive.sum())
        if dead == 0:
            return
        for F in self.followers:
            F.sync()
        keep = np.flatnonzero(self.alive)
        super().compact(force)
        kd = _d(keep)
        for F in self.followers:
            F.M = F.M[:, kd]
        B.DV.free_pool()


class Stop(Exception):
    pass


# =============================================================================== the lift

class Lift:
    def __init__(self, eng, work_gb: float = 2.0, cache_gb=None, tile_mb: float = 256.0, tile_rows=None,
                 col_block=None, batch_units=None, direct_rows=None, compact_frac: float = 0.1,
                 stop_on_one: bool = False, include_low: bool = False, fall_every: int = 8, status=None, log=True,
                 log_seconds: float = 20.0, ckpt=None, ckpt_minutes: float = 0.0, driver_sha: str = '',
                 test_drop_shift=None, test_no_fall_shifts: bool = False, vd_info=None):
        """``eng``: an ``ic_bounds_big.Engine`` holding the finished degree-D closure (split)."""
        assert eng.split, 'the lift needs the swap split'
        self.eng = eng
        self.k, self.D, self.p = eng.k, eng.D, eng.p
        assert self.p % 2 == 1
        self.mi = B.max_inner(self.p)
        self.inv2 = float(pow(2, -1, self.p))
        self.work_bytes = int(work_gb * (1 << 30))
        self.cache_bytes = None if cache_gb is None else int(cache_gb * (1 << 30))
        self.batch_units_opt, self.direct_rows_opt = batch_units, direct_rows
        self.stop_on_one, self.include_low = stop_on_one, include_low
        self.fall_every = max(1, int(fall_every))
        self.status_path, self.logging, self.log_seconds = status, log, log_seconds
        self.ckpt, self.ckpt_minutes, self.driver_sha = ckpt, ckpt_minutes, driver_sha
        self.test_drop_shift, self.test_no_fall_shifts = test_drop_shift, test_no_fall_shifts
        self.vd_info = vd_info or {}
        self.t_start, self.elapsed_prev, self.last_log, self.last_ckpt = time.time(), 0.0, 0.0, time.time()
        self.term = False
        self.st = dict(images=0, eterm_rows=0, madds_gemm=0, candidates=0, candidates_nonzero=0, batches=0,
                       follower_builds=0, follower_evictions=0, fall_images=0)
        self.cache: OrderedDict = OrderedDict()
        self.target = None
        self.contains_one = False
        self.fixed_point = False
        self.finished = False
        self.resumes = 0
        self.rounds = 0
        self.falls_new: list = []
        self.stage = dict(name='main', pos=0)
        self._static()
        ke = _KEng(tile_mb * (1 << 20), self.work_bytes, compact_frac, col_block)
        self.K = [KPart(self.Tdim[q], self.p, ke, tile_rows) for q in (0, 1)]
        self.fall_done = [[np.zeros(self.Tdim[q0], dtype=bool) for _ in (0, 1)] for q0 in (0, 1)]
        self.reg_count = [np.zeros(self.nReg[q], dtype=np.int64) for q in (0, 1)]

    # ----------------------------------------------------------------------------- static bookkeeping (host)
    def _static(self):
        k, D = self.k, self.D
        eng = self.eng
        for P in eng.parts:
            P.compact()
        self.rankD = [P.rank for P in eng.parts]
        L1 = self.L1 = B.Layout(k, D + 1, True)
        LD = eng.L
        self.s0 = [L1.parts[q]['s0'] for q in (0, 1)]
        for q in (0, 1):                      # the low columns of the degree-(D+1) layout ARE the degree-D columns
            P1, PD = L1.parts[q], LD.parts[q]
            assert P1['c'] - self.s0[q] == PD['c'] == eng.parts[q].c
            r1 = P1['reps'][self.s0[q]:]
            assert np.array_equal(L1.a1[r1], LD.a1[PD['reps']]) and np.array_equal(L1.a2[r1], LD.a2[PD['reps']])
            assert np.array_equal(P1['w'][self.s0[q]:], PD['w'])
        P0, P1 = eng.parts
        self.fr = [P.free_cols().copy() for P in eng.parts]
        self.nfr = [len(f) for f in self.fr]
        self.frpos = []
        for q in (0, 1):
            pos = -np.ones(eng.parts[q].c, dtype=np.int64)
            pos[self.fr[q]] = np.arange(self.nfr[q])
            self.frpos.append(pos)
            assert np.array_equal(eng.parts[q].pc2col, self.fr[q])
        self.rowD = [P.col2row.copy() for P in eng.parts]
        # regular weight-D orbits: pivot in both parts (self-paired: pivot of the + part)
        c0 = np.flatnonzero((P0.w == D) & (self.rowD[0] >= 0))
        h = L1.parts[0]['reps'][self.s0[0] + c0]
        diag = L1.a1[h] == L1.a2[h]
        c1 = np.where(diag, -1, L1.parts[1]['lift_idx'][h] - self.s0[1])
        ok = diag | (self.rowD[1][np.maximum(c1, 0)] >= 0)
        c0, c1, h, diag = c0[ok], c1[ok], h[ok], diag[ok]
        self.o_c0, self.o_c1, self.o_diag = c0, c1, diag
        self.o_row0 = self.rowD[0][c0]
        self.o_row1 = np.where(diag, -1, self.rowD[1][np.maximum(c1, 0)])
        self.o_J = L1.a1[h] | L1.a2[h]
        self.nO = len(c0)
        self.regcol = [np.zeros(P.c, dtype=bool) for P in eng.parts]
        self.regcol[0][c0] = True
        self.regcol[1][c1[~diag]] = True
        # free-type weight-D monomials (full coordinates) and special weight-(D+1) monomials
        regfull = np.zeros(L1.c, dtype=bool)
        regfull[h] = True
        regfull[L1.index(L1.a2[h], L1.a1[h])] = True
        ft = np.flatnonzero((L1.w == D) & ~regfull)
        special = np.zeros(L1.c, dtype=bool)
        fa1, fa2 = L1.a1[ft], L1.a2[ft]
        for i in range(k):
            bit = 1 << i
            s = (fa1 & bit) == 0
            special[L1.index(fa1[s] | bit, fa2[s])] = True
            s = (fa2 & bit) == 0
            special[L1.index(fa1[s], fa2[s] | bit)] = True
        assert np.all(L1.w[special] == D + 1)
        self.n_freetype = len(ft)
        self.tcol_top, self.regid, self.nSp, self.nReg, self.Tdim, self.reg_expect = [], [], [], [], [], []
        for q in (0, 1):
            hq = L1.parts[q]['reps'][:self.s0[q]]
            sp = special[hq]
            tt = -np.ones(self.s0[q], dtype=np.int64)
            tt[sp] = np.arange(int(sp.sum()))
            rg = -np.ones(self.s0[q], dtype=np.int64)
            rg[~sp] = np.arange(int((~sp).sum()))
            self.tcol_top.append(tt)
            self.regid.append(rg)
            self.nSp.append(int(sp.sum()))
            self.nReg.append(int((~sp).sum()))
            self.Tdim.append(self.nSp[q] + self.nfr[q])
            hr = hq[~sp]
            pc = np.array([bin(int(x)).count('1') for x in L1.a1[hr]], dtype=np.int64) if len(hr) else np.zeros(0, np.int64)
            self.reg_expect.append(np.where(L1.a1[hr] == L1.a2[hr], pc, D + 1))
        # rows processed directly: weight-D pivots whose orbit is free-type (and the low rows on request)
        self.direct = []
        for q, P in enumerate(eng.parts):
            wp = P.w[P.piv]
            sel = (wp == D) & ~self.regcol[q][P.piv]
            if self.include_low:
                sel |= wp < D
            self.direct.append(np.flatnonzero(sel))
        # units (o, j) sorted by J-group
        Uo = np.repeat(np.arange(self.nO), k)
        Uj = np.tile(np.arange(k), self.nO)
        if self.test_drop_shift is not None:
            kp = Uj != self.test_drop_shift
            Uo, Uj = Uo[kp], Uj[kp]
        UJ = self.o_J[Uo] | (1 << Uj)
        o = np.lexsort((Uo, Uj, UJ))
        self.U_o, self.U_j, UJ = Uo[o], Uj[o], UJ[o]
        self.group_start = np.concatenate([[0], np.flatnonzero(np.diff(UJ)) + 1, [len(UJ)]]) if len(UJ) else np.array([0])
        self.zero_T = self.nSp[0] + int(self.frpos[0][P0.zero_col]) if self.frpos[0][P0.zero_col] >= 0 else -1

    def plan_batches(self, q):
        """static batch plan (whole J-groups, at most ``nu`` units; depends only on sizes fixed at the start)."""
        nu = self.batch_units_opt or max(64, int(0.7 * self.work_bytes // (80 * max(1, self.Tdim[q]) + 32 * max(self.nfr))))
        gs = self.group_start
        out, a = [], 0
        while a < len(gs) - 1:
            b = a + 1
            while b < len(gs) - 1 and gs[b + 1] - gs[a] <= nu:
                b += 1
            out.append((int(gs[a]), int(gs[b])))
            a = b
        return out

    def direct_chunk(self, q):
        return self.direct_rows_opt or max(16, int(0.7 * self.work_bytes // (40 * max(1, self.Tdim[q]) + 16 * max(self.nfr))))

    # ----------------------------------------------------------------------------- unit images
    def slots(self, p, j, q, dcols):
        """C = project_q . u_{1,j} . lift_p applied to the unit vectors of the degree-D columns ``dcols`` of
        part p: two terms (mA, cA) (image of the orbit representative) and (mB, cB) (image of its swap),
        m = column of part q in the degree-(D+1) layout, c = coefficient (0: no term)."""
        L1 = self.L1
        Pp, Q = L1.parts[p], L1.parts[q]
        s = np.asarray(dcols) + self.s0[p]
        hrep, hsig = Pp['reps'][s], Pp['rep_sig'][s]
        sh = L1.shift[(0, j)]
        hA = sh[hrep - L1.low0]
        hB = sh[hsig - L1.low0]
        mA, cA = Q['lift_idx'][hA], Q['proj_coef'][hA].copy()
        mB, cB = Q['lift_idx'][hB], Pp['sign'] * Q['proj_coef'][hB]
        cB[hsig == hrep] = 0
        return mA, cA, mB, cB

    def eterm(self, q, n, terms):
        """terms: list of (rowidx, m, coef) with unique rowidx inside each triple.  Splits the terms into
        regular top terms (returned as host arrays rowidx, regular id, coef) and the rest, which is reduced
        modulo V_D^q (low pivot column -> minus its basis row) and modulo K_q.
        Returns (n x K.ncols device float64 or None, (rowidx, regid, coef))."""
        xp = B.xp
        p = self.p
        K, P, s0, nSp = self.K[q], self.eng.parts[q], self.s0[q], self.nSp[q]
        regs, units, pivs = [], [], []
        skip = K.ncols == 0              # K_q = T_q: everything outside the regular columns reduces to zero
        for ri, m, cf in terms:
            cf = np.asarray(cf, dtype=np.int64) % p
            nz = cf != 0
            ri, m, cf = np.asarray(ri)[nz], np.asarray(m)[nz], cf[nz]
            top = m < s0
            mt = np.minimum(m, max(s0 - 1, 0))
            rg = np.where(top, self.regid[q][mt] if s0 else -1, -1)
            isreg = top & (rg >= 0)
            regs.append((ri[isreg], rg[isreg], cf[isreg]))
            sp = top & ~isreg
            if sp.any():
                units.append((ri[sp], self.tcol_top[q][m[sp]], cf[sp]))
            low = ~top
            d = np.where(low, m - s0, 0)
            row = self.rowD[q][d]
            fl = low & (row < 0)
            if fl.any():
                units.append((ri[fl], nSp + self.frpos[q][d[fl]], cf[fl]))
            pv = low & (row >= 0)
            if pv.any():
                pivs.append((ri[pv], row[pv], cf[pv]))
        reg = tuple(np.concatenate([r[i] for r in regs]) for i in range(3))
        if skip or (not units and not pivs):
            return None, reg
        cols = np.unique(np.concatenate([u[1] for u in units])) if units else np.zeros(0, dtype=np.int64)
        if pivs:
            cols = np.union1d(cols, np.arange(nSp, nSp + self.nfr[q]))
        Yt = xp.zeros((n, len(cols)), dtype=F64)
        for ri, c, cf in units:
            Yt[_d(ri), _d(np.searchsorted(cols, c))] += _d(cf.astype(F64))
        if pivs:
            f0 = int(np.searchsorted(cols, nSp))
            for ri, row, cf in pivs:
                G = rows_dev(P, row)
                G *= _d(((p - cf) % p).astype(F64))[:, None]
                B.modp_(G, p)
                Yt[_d(ri), f0:] += G
                del G
        B.modp_(Yt, p)
        self.st['eterm_rows'] += n
        return K.reduce_keep(cols, Yt), reg

    # ----------------------------------------------------------------------------- followers
    def set_target(self, q):
        """the followers of both targets stay cached (least recently used first out, under ``cache_bytes``)."""
        self.target = q

    def drop_cache(self):
        for K in self.K:
            K.followers = []
        self.cache.clear()
        B.DV.free_pool()

    def follower(self, p, j, keep=()):
        q = self.target
        K = self.K[q]
        key = (q, p, j)
        F = self.cache.get(key)
        if F is None:
            need = self.nfr[p] * K.ncols * 4
            if self.cache_bytes is not None:
                while sum(f.nbytes for f in self.cache.values()) + need > self.cache_bytes:
                    vic = next((x for x in self.cache if x != key and x not in keep), None)
                    if vic is None:
                        break
                    Fv = self.cache.pop(vic)
                    Fv.K.followers.remove(Fv)
                    del Fv
                    self.st['follower_evictions'] += 1
                    B.DV.free_pool()
            F = self._build_follower(p, j, q)
            self.cache[key] = F
            K.followers.append(F)
        else:
            F.sync()
        self.cache.move_to_end(key)
        return F

    def _build_follower(self, p, j, q):
        """rows red_q(C e_f) mod K_q for the free columns f of part p."""
        xp = B.xp
        K = self.K[q]
        n = self.nfr[p]
        M = xp.zeros((n, K.ncols), dtype=F32)
        ch = max(16, int(0.5 * self.work_bytes // (8 * (2 * self.Tdim[q] + 2 * self.nfr[q] + 1))))
        for a in range(0, n, ch):
            d = self.fr[p][a:a + ch]
            mA, cA, mB, cB = self.slots(p, j, q, d)
            ar = np.arange(len(d))
            E, reg = self.eterm(q, len(d), [(ar, mA, cA), (ar, mB, cB)])
            assert len(reg[0]) == 0, 'image of a free column has a regular top term'
            if E is not None:
                M[a:a + len(d)] = E.astype(F32)
            del E
        self.st['follower_builds'] += 1
        return Follower(K, M)

    def gemm(self, A, F):
        """A (n x nfr_p, device float64 in [0,p)) @ follower, mod p: n x K.ncols float64."""
        xp = B.xp
        K = F.K
        n, nc = A.shape[0], K.ncols
        out = xp.zeros((n, nc), dtype=F64)
        if n == 0 or nc == 0 or A.shape[1] == 0:
            return out
        cw = max(64, int(0.2 * self.work_bytes // (8 * (A.shape[1] + n))))
        for a in range(0, nc, cw):
            b = min(nc, a + cw)
            addmul(out[:, a:b], A, F.M[:, a:b].astype(F64), self.p, self.mi)
        self.st['madds_gemm'] += n * A.shape[1] * nc
        self.st['images'] += n
        return out

    # ----------------------------------------------------------------------------- adding candidates
    def add(self, q, C):
        """candidate rows (device, on K_q's physical columns, reduced modulo the current K_q) -> K_q."""
        xp = B.xp
        K = self.K[q]
        n = C.shape[0]
        self.st['candidates'] += n
        if n == 0 or K.ncols == 0:
            return 0
        nz = np.flatnonzero(_h(xp.any(C != 0, axis=1)))
        if len(nz) == 0:
            return 0
        self.st['candidates_nonzero'] += len(nz)
        R = C[_d(nz)] if len(nz) < n else C
        t = K.add_residuals(R)
        if K.nfree == 0 and K.ncols:        # K_q = T_q: drop the dead columns so that later work is skipped
            K.compact()
        if q == 0 and self.stop_on_one and self.one_now():
            raise Stop('one')
        return t

    def one_now(self) -> bool:
        P0 = self.eng.parts[0]
        if P0.zero_col >= 0 and self.rowD[0][P0.zero_col] >= 0:
            return True
        return self.zero_T >= 0 and bool(self.K[0].col2row[self.zero_T] >= 0)

    # ----------------------------------------------------------------------------- main pass
    def do_direct(self, q, p, j, rows):
        P = self.eng.parts[p]
        F = self.follower(p, j)
        Y = self.gemm(rows_dev(P, rows), F)
        mA, cA, mB, cB = self.slots(p, j, q, P.piv[rows])
        ar = np.arange(len(rows))
        E, reg = self.eterm(q, len(rows), [(ar, mA, cA), (ar, mB, cB)])
        assert len(reg[0]) == 0, 'direct row with a regular top term'
        if E is not None:
            Y += E
            B.modp_(Y, self.p)
        self.add(q, Y)

    def do_units(self, q, a, b):
        xp = B.xp
        p, D = self.p, self.D
        K = self.K[q]
        Uo, Uj = self.U_o[a:b], self.U_j[a:b]
        off = ~self.o_diag[Uo]
        ng = len(Uo) + int(off.sum())
        G = xp.zeros((ng, K.ncols), dtype=F64)
        gm = np.zeros(ng, dtype=np.int64)
        gc = np.zeros(ng, dtype=np.int64)
        pos = 0
        for j in np.unique(Uj):
            j = int(j)
            s = Uj == j
            os_, offj = Uo[s], off[s]
            n1, n2 = len(os_), int(offj.sum())
            mA, cA, mB, cB = self.slots(0, j, q, self.o_c0[os_])
            if K.ncols == 0:             # K_q = T_q: only the regular top terms matter (coverage check)
                gm[pos:pos + n1], gc[pos:pos + n1] = mA, cA
                gm[pos + n1:pos + n1 + n2], gc[pos + n1:pos + n1 + n2] = mB[offj], cB[offj]
                pos += n1 + n2
                continue
            F0 = self.follower(0, j)
            Y0 = self.gemm(rows_dev(self.eng.parts[0], self.o_row0[os_]), F0)
            if n2:
                F1 = self.follower(1, j, keep=((q, 0, j),))
                Y1 = self.gemm(rows_dev(self.eng.parts[1], self.o_row1[os_[offj]]), F1)
                od = _d(np.flatnonzero(offj))
                Yo = Y0[od]
                g2 = Yo - Y1
                g2 *= self.inv2
                B.modp_(g2, p)
                Yo += Y1
                Yo *= self.inv2
                B.modp_(Yo, p)
                Y0[od] = Yo
                G[pos + n1:pos + n1 + n2] = g2
                gm[pos + n1:pos + n1 + n2], gc[pos + n1:pos + n1 + n2] = mB[offj], cB[offj]
                del Y1, Yo, g2
            G[pos:pos + n1] = Y0
            gm[pos:pos + n1], gc[pos:pos + n1] = mA, cA
            del Y0
            pos += n1 + n2
        assert pos == ng
        E, (ri, rg, cf) = self.eterm(q, ng, [(np.arange(ng), gm, gc)])
        if E is not None:
            G += E
            B.modp_(G, p)
            del E
        isreg = np.zeros(ng, dtype=bool)
        isreg[ri] = True
        o = np.argsort(rg, kind='stable')
        ri, rg, cf = ri[o], rg[o], cf[o]
        if len(ri):
            inv = np.array([pow(int(c), -1, p) for c in np.unique(cf)], dtype=np.int64)
            invs = inv[np.searchsorted(np.unique(cf), cf)]
            ne = np.flatnonzero(invs != 1)
            if len(ne):
                X = G[_d(ri[ne])]
                X *= _d(invs[ne].astype(F64))[:, None]
                B.modp_(X, p)
                G[_d(ri[ne])] = X
                del X
            cnt = np.bincount(rg, minlength=self.nReg[q])
            ids = np.flatnonzero(cnt)
            if self.test_drop_shift is None:
                assert np.all(self.reg_count[q][ids] == 0), 'regular column met in two batches'
                assert np.array_equal(cnt[ids], self.reg_expect[q][ids]), 'wrong number of generators of a regular column'
            self.reg_count[q] += cnt
        same = np.flatnonzero(rg[1:] == rg[:-1]) if len(ri) else np.zeros(0, dtype=np.int64)
        dr = np.flatnonzero(~isreg)
        ns = len(same)
        if ns + len(dr) == 0:
            return
        C = xp.zeros((ns + len(dr), K.ncols), dtype=F64)
        if ns:
            Cd = G[_d(ri[same + 1])]
            Cd -= G[_d(ri[same])]
            B.modp_(Cd, p)
            C[:ns] = Cd
            del Cd
        if len(dr):
            C[ns:] = G[_d(dr)]
        del G
        self.add(q, C)

    def schedule(self):
        """static list of work items: the direct rows of both targets, then the unit batches, the two targets
        alternating and the batches of each target in a stride order (a spread of J-groups early makes K_q grow
        early, and everything later is reduced modulo it).  The order does not affect the result: every item is
        processed exactly once."""
        shifts = [j for j in range(self.k) if j != self.test_drop_shift]
        items = []
        for q in (0, 1):
            ch = self.direct_chunk(q)
            items += [('d', q, p, j, a, ch) for p in (0, 1) for j in shifts for a in range(0, len(self.direct[p]), ch)]
        bq = [self.plan_batches(q) for q in (0, 1)]
        orders = []
        for q in (0, 1):
            n = len(bq[q])
            st = max(1, int(round(n * 0.6180339887)))
            while n and np.gcd(st, n) != 1:
                st += 1
            orders.append([(i * st) % n for i in range(n)])
        for t in range(max(len(bq[0]), len(bq[1]))):
            for q in (0, 1):
                if t < len(bq[q]):
                    items.append(('u', q) + tuple(bq[q][orders[q][t]]))
        return items

    def main_phase(self):
        items = self.schedule()
        stg = self.stage
        since = 0
        while stg['pos'] < len(items):
            it = items[stg['pos']]
            self.set_target(it[1])
            if it[0] == 'd':
                _, q, p, j, a, ch = it
                self.do_direct(q, p, j, self.direct[p][a:a + ch])
            else:
                _, q, a, b = it
                self.do_units(q, a, b)
                self.st['batches'] += 1
                since += 1
            stg['pos'] += 1
            if since >= self.fall_every:
                since = 0
                self.fall_round()
            self.tick(total=len(items))
        if self.test_drop_shift is None:
            for q in (0, 1):
                assert np.array_equal(self.reg_count[q], self.reg_expect[q]), 'a regular column was not covered'

    def fall_rows_idx(self, q0):
        K = self.K[q0]
        return np.flatnonzero(K.piv >= self.nSp[q0])

    def fall_block(self, q0, rows):
        """fall rows of K_{q0} as coordinates on the free columns of V_D^{q0}: n x nfr device float64."""
        xp = B.xp
        K, nSp, n = self.K[q0], self.nSp[q0], len(rows)
        ps0 = int(np.searchsorted(K.pc2col, nSp))
        Aphys = xp.zeros((n, K.ncols - ps0), dtype=F64)
        for ti, idx, off in K._groups(rows):
            Aphys[_d(idx)] = K._rows(ti, off, slice(ps0, None)).astype(F64)
        A = xp.zeros((n, self.nfr[q0]), dtype=F64)
        A[:, _d(K.pc2col[ps0:] - nSp)] = Aphys
        del Aphys
        A[_d(np.arange(n)), _d(K.piv[rows] - nSp)] = 1.0
        return A

    def fall_round(self) -> bool:
        """one round: every fall row of K_+ and K_- that has not been shifted into target q yet (q = +, -) is
        read at its present state, marked and shifted.  Returns False if nothing was pending."""
        pend = {}
        for q0 in (0, 1):
            pv = self.K[q0].piv[self.fall_rows_idx(q0)]
            for q in (0, 1):
                pend[(q0, q)] = pv[~self.fall_done[q0][q][pv]]
        if all(len(v) == 0 for v in pend.values()):
            return False
        self.falls_new.append([int(max(len(pend[(q0, 0)]), len(pend[(q0, 1)]))) for q0 in (0, 1)])
        for q in (0, 1):
            if self.K[q].ncols == 0:          # K_q = T_q already: every image reduces to zero
                for q0 in (0, 1):
                    self.fall_done[q0][q][pend[(q0, q)]] = True
                continue
            self.set_target(q)
            for q0 in (0, 1):
                pv = pend[(q0, q)]
                ch = self.direct_chunk(q)
                for a in range(0, len(pv), ch):
                    cols = pv[a:a + ch]
                    rows = self.K[q0].col2row[cols]
                    A = self.fall_block(q0, rows)          # the state that is processed
                    self.fall_done[q0][q][cols] = True
                    for j in range(self.k):
                        Y = self.gemm(A, self.follower(q0, j))
                        self.st['fall_images'] += len(cols)
                        self.add(q, Y)
                        del Y
                    del A
        self.rounds += 1
        return True

    def absorb_falls(self):
        """add the falls to the V_D parts: they become the reduced echelon bases of V_{D+1} n M_<=D (or of
        the part of it found so far)."""
        if self.finished:
            return
        self.drop_cache()
        if self.eng.tile_budget is not None:
            self.eng.tile_budget = 1 << 62
        for q, P in enumerate(self.eng.parts):
            rows = self.fall_rows_idx(q)
            ch = max(16, int(0.3 * self.work_bytes // (8 * 3 * max(1, self.nfr[q]))))
            for a in range(0, len(rows), ch):
                A = self.fall_block(q, rows[a:a + ch])
                if P.ncols != self.nfr[q]:        # the part compacted: drop the columns that became pivots
                    A = A[:, _d(self.frpos[q][P.pc2col])]
                P.add_residuals(A)
                del A
        self.finished = True

    # ----------------------------------------------------------------------------- run
    def run(self) -> str:
        """'done' (fixed point), 'one' (stopped on 1 in V_{D+1}), 'sigterm', 'nofall' (test hook)."""
        stg = self.stage
        try:
            if stg['name'] == 'main':
                self.main_phase()
                stg.update(name='falls')
                self.maybe_checkpoint(force=bool(self.ckpt))
            if self.test_no_fall_shifts:
                self.absorb_falls()
                self.contains_one = self.eng.zero_one()
                return 'nofall'
            while self.fall_round():
                self.tick(total=0)
            self.fixed_point = True
        except Stop as e:
            if e.args and e.args[0] == 'one':
                self.contains_one = True
                self.emit(force=True, event='contains_one')
                self.absorb_falls()
                return 'one'
            return 'sigterm'
        self.absorb_falls()
        self.contains_one = self.eng.zero_one()
        stg.update(name='done')
        return 'done'

    # ----------------------------------------------------------------------------- logging, checkpoints
    def seconds(self):
        return self.elapsed_prev + time.time() - self.t_start

    def madds(self):
        km = sum(K.madds[x] for K in self.K for x in ('reduce', 'update', 'rref'))
        fo = sum(K.madds_follow for K in self.K)
        return dict(gemm=self.st['madds_gemm'], k_reduce_update_rref=km, follower_sync=fo, total=self.st['madds_gemm'] + km + fo)

    def snapshot(self, **extra):
        used, total = B.DV.mem_info()
        d = dict(t=round(self.seconds(), 1), stage=dict(self.stage), maxrss_mb=B.maxrss_mb(), dev_used_mb=used >> 20,
                 dev_total_mb=total >> 20, dev_pool_mb=B.DV.pool_bytes() >> 20, rank_K=[K.rank for K in self.K],
                 ncols_K=[K.ncols for K in self.K], dim_T=self.Tdim, cache_mb=sum(f.nbytes for f in self.cache.values()) >> 20,
                 cache_n=len(self.cache), madds=self.madds()['total'], rounds=self.rounds, **self.st)
        d.update(extra)
        return d

    def emit(self, force=False, **extra):
        now = time.time()
        if not force and now - self.last_log < self.log_seconds:
            return
        self.last_log = now
        d = self.snapshot(**extra)
        if self.logging:
            print(json.dumps(d), flush=True)
        if self.status_path:
            tmp = self.status_path + '.tmp'
            with open(tmp, 'w') as f:
                json.dump(d, f)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, self.status_path)

    def tick(self, total):
        """between batches: progress line with an estimate of the remaining time of the current stage,
        SIGTERM handling, timed checkpoint."""
        stg = self.stage
        extra = {}
        if total:
            if getattr(self, '_eta_t0', None) is None:
                self._eta_t0, self._eta_p0 = time.time(), stg['pos'] - 1
            done = stg['pos'] - self._eta_p0
            # remaining time of the main phase at the average rate so far (an upper estimate once K has grown:
            # later batches are reduced modulo a larger K and are cheaper)
            extra = dict(of=total, eta_main_s=round((time.time() - self._eta_t0) / max(1, done) * (total - stg['pos']), 1))
        self.emit(**extra)
        if self.term:
            self.emit(force=True, event='sigterm')
            if self.ckpt:
                self.save_checkpoint('sigterm')
            raise Stop('sigterm')
        self.maybe_checkpoint()

    def maybe_checkpoint(self, force=False):
        if not self.ckpt:
            return
        if force or (self.ckpt_minutes and time.time() - self.last_ckpt >= 60 * self.ckpt_minutes):
            self.save_checkpoint('timer')

    def _ident(self):
        return dict(k=self.k, D=self.D, p=self.p, rankD=self.rankD, Tdim=self.Tdim, nReg=self.nReg, lift_sha256=LIFT_SHA256,
                    work_bytes=self.work_bytes, batch_units=self.batch_units_opt, direct_rows=self.direct_rows_opt,
                    include_low=self.include_low, vd=self.vd_info.get('basis_id'), drop=self.test_drop_shift)

    def save_checkpoint(self, reason):
        """K_0, K_1, the fall rows already shifted and the position in the static schedule (item boundaries).
        One slot: written to a temporary directory, renamed, then the old one is removed."""
        import shutil
        t0 = time.time()
        root = self.ckpt
        os.makedirs(root, exist_ok=True)
        tmp, final = os.path.join(root, '.tmp-lift'), os.path.join(root, 'lift-ck')
        if os.path.exists(tmp):
            shutil.rmtree(tmp)
        os.makedirs(tmp)
        for q, K in enumerate(self.K):
            K.compact()
            rows = [_h(K.tiles[ti][:K.tile_rows_filled(ti)]) for ti in range(len(K.tiles)) if K.tile_rows_filled(ti) > 0]
            M = np.concatenate(rows) if rows else np.zeros((0, K.ncols), dtype=F32)
            np.savez(os.path.join(tmp, f'K{q}.npz'), piv=K.piv, M=M, reg_count=self.reg_count[q],
                     fall_done0=self.fall_done[q][0], fall_done1=self.fall_done[q][1])
        man = dict(ident=self._ident(), stage=self.stage, st=self.st, elapsed=self.seconds(), resumes=self.resumes,
                   rounds=self.rounds, falls_new=self.falls_new,
                   kmadds=[[K.madds, K.madds_follow] for K in self.K], reason=reason,
                   files={f'K{q}.npz': B.file_sha256(os.path.join(tmp, f'K{q}.npz')) for q in (0, 1)})
        with open(os.path.join(tmp, 'manifest.json'), 'w') as f:
            json.dump(man, f)
            f.flush()
            os.fsync(f.fileno())
        old = final + '.old'
        if os.path.exists(old):
            shutil.rmtree(old)
        if os.path.exists(final):
            os.rename(final, old)
        os.rename(tmp, final)
        if os.path.exists(old):
            shutil.rmtree(old)
        self.last_ckpt = time.time()
        self.emit(force=True, event='checkpoint', reason=reason, ckpt_seconds=round(time.time() - t0, 1))

    def load_checkpoint(self) -> bool:
        xp = B.xp
        root = self.ckpt
        for name in ('lift-ck', 'lift-ck.old'):
            path = os.path.join(root, name)
            mp = os.path.join(path, 'manifest.json')
            if not os.path.exists(mp):
                continue
            man = json.load(open(mp))
            if any(not os.path.exists(os.path.join(path, f)) or B.file_sha256(os.path.join(path, f)) != s for f, s in man['files'].items()):
                continue
            if man['ident'] != json.loads(json.dumps(self._ident())):
                raise SystemExit('lift checkpoint belongs to another run (basis, code or batch settings differ)')
            for q, K in enumerate(self.K):
                z = np.load(os.path.join(path, f'K{q}.npz'))
                K.piv = z['piv'].astype(np.int64)
                ispiv = np.zeros(K.c, dtype=bool)
                ispiv[K.piv] = True
                K.pc2col = np.flatnonzero(~ispiv)
                K.alive = np.ones(len(K.pc2col), dtype=bool)
                K.col2phys = -np.ones(K.c, dtype=np.int64)
                K.col2phys[K.pc2col] = np.arange(len(K.pc2col))
                K.col2row = -np.ones(K.c, dtype=np.int64)
                K.col2row[K.piv] = np.arange(len(K.piv))
                K.tiles, K.loc = [], []
                M = z['M']
                for a in range(0, M.shape[0], K.TR):
                    tile = xp.zeros((K.TR, len(K.pc2col)), dtype=F32)
                    blk = M[a:a + K.TR]
                    tile[:blk.shape[0]] = _d(blk)
                    K.tiles.append(tile)
                    K.loc.append('d')
                K.madds, K.madds_follow = man['kmadds'][q]
                self.reg_count[q] = z['reg_count'].astype(np.int64)
                self.fall_done[q] = [z['fall_done0'].astype(bool), z['fall_done1'].astype(bool)]
            self.stage, self.st, self.elapsed_prev = man['stage'], man['st'], man['elapsed']
            self.resumes = man['resumes'] + 1
            self.rounds, self.falls_new = man['rounds'], man['falls_new']
            self.t_start = time.time()
            self.emit(force=True, event='resumed', from_ckpt=name)
            return True
        return False

    def install_sigterm(self):
        import signal

        def h(signum, frame):
            self.term = True
        signal.signal(signal.SIGTERM, h)

    # ----------------------------------------------------------------------------- record
    def record(self, status: str) -> dict:
        D = self.D
        eng = self.eng
        complete = status == 'done' and self.fixed_point and self.test_drop_shift is None and not self.test_no_fall_shifts
        dimK = [K.rank for K in self.K]
        nfalls = [int(len(self.fall_rows_idx(q))) for q in (0, 1)]
        reg_done = [int(np.count_nonzero(c)) for c in self.reg_count]
        dim = sum(self.rankD) + sum(dimK) + sum(reg_done)
        pr = [P.profile(D) for P in eng.parts]           # after absorb_falls: V_{D+1} n M_<=D
        prof = [a + b for a, b in zip(*pr)] + [dim]
        rec = dict(k=self.k, D=D, D_target=D + 1, p=self.p, dim_VD=sum(self.rankD), part_ranks_VD=self.rankD,
                   free_columns_VD=self.nfr, regular_orbits_weight_D=self.nO, free_type_monomials_weight_D=self.n_freetype,
                   regular_top_orbits=self.nReg, special_top_orbits=self.nSp, dim_T=self.Tdim, dim_K=dimK,
                   falls=nfalls, falls_new_per_round=self.falls_new, rounds=self.rounds,
                   contains_one=bool(self.contains_one), fixed_point=bool(complete), status=status,
                   dim_VD1=dim, part_dims_VD1=[self.rankD[q] + dimK[q] + reg_done[q] for q in (0, 1)],
                   profile=prof, part_profiles_low=pr,
                   low_basis_sha256=[P.sha() for P in eng.parts], low_part_ranks=[P.rank for P in eng.parts],
                   bounds_only=not complete,
                   regular_columns_covered=reg_done,
                   reductions=self.st['images'] + self.st['eterm_rows'], counts=dict(self.st), madds=self.madds(),
                   k_madds_detail=[dict(K.madds) for K in self.K], seconds=round(self.seconds(), 1),
                   device=eng.device, maxrss_mb=B.maxrss_mb(), resumes=self.resumes,
                   settings=dict(work_gb=self.work_bytes / (1 << 30), cache_gb=None if self.cache_bytes is None else self.cache_bytes / (1 << 30),
                                 batch_units=self.batch_units_opt, direct_rows=self.direct_rows_opt, include_low=self.include_low,
                                 fall_every=self.fall_every, stop_on_one=self.stop_on_one, tile_rows_K=[K.TR for K in self.K]),
                   test_hooks=dict(drop_shift=self.test_drop_shift, no_fall_shifts=self.test_no_fall_shifts),
                   vd=self.vd_info, guard_hits=B.GUARD_HITS[0],
                   source_sha256=dict(lift=LIFT_SHA256, engine=B.ENGINE_SHA256, driver=self.driver_sha),
                   lift_version=LIFT_VERSION, engine_version=B.ENGINE_VERSION)
        if not complete:
            rec['note'] = ('NOT a fixed point: dim_VD1, profile and the basis hashes describe a subspace of V_{D+1} '
                           '(lower bounds only); contains_one=true is exact, contains_one=false is NOT a statement about V_{D+1}')
        used, total = B.DV.mem_info()
        rec.update(dev_used_mb=used >> 20, dev_total_mb=total >> 20,
                   dev_peak_live_mb=(B._Live.peak >> 20) if eng.device == 'fake-gpu' else None)
        return rec


# =============================================================================== V_D sources

def compute_vd(g: dict, k: int, p: int, D: int, device: str = 'cpu', gens: str = 'all', **kw):
    """run the engine's closure for V_D (small cases); returns the engine (closed, full sweep clean)."""
    eng = B.Engine(g, k, p, D, device=device, gens=gens, log=kw.pop('log', False), **kw)
    st = eng.run()
    if st != 'done' or not eng.closed:
        raise SystemExit(f'closure for V_D did not finish closed: status {st}')
    return eng


def load_vd(ckpt: str, device: str = 'cpu', basis_on: str = 'device', work_gb: float = 2.0, force: bool = False):
    """engine holding the basis of a finished degree-D closure read from an engine checkpoint directory.
    ``basis_on='host'``: the tiles stay in host memory and batch rows are transferred when needed."""
    path, man = B.Engine.find_checkpoint(ckpt)
    if path is None:
        raise SystemExit(f'no complete engine checkpoint in {ckpt}')
    cfg, st = man['config'], man['state']
    if not cfg['split']:
        raise SystemExit('the lift needs a split closure')
    if not st['closure_complete'] and not force:
        raise SystemExit('checkpoint is not a finished closure (closure_complete false); --force-ckpt overrides')
    kw = {}
    if basis_on == 'host' and device != 'cpu':
        kw['tile_budget_mb'] = 1 << 40           # budget set: every loaded tile starts on the host
    eng = B.Engine({(0, 0): 1}, cfg['k'], cfg['p'], cfg['D'], device=device, gens=cfg['gens'], log=False,
                   work_gb=work_gb, ckpt=ckpt, source_id=cfg['source_id'], **kw)
    B._Live.limit = None
    if not eng.load_checkpoint(force=True):
        raise SystemExit('checkpoint could not be loaded')
    info = dict(source='checkpoint', ckpt=os.path.basename(path), ckpt_reason=man.get('reason'), ckpt_phase=st['phase'],
                closure_complete=st['closure_complete'], ckpt_engine_sha256=man['engine_sha256'], ckpt_config=cfg,
                swept_closed='not known from the checkpoint: take it from the closure record of the same run',
                basis_id=hashlib.sha256(json.dumps(man['files'], sort_keys=True).encode()).hexdigest(), basis_on=basis_on)
    return eng, info


# =============================================================================== direct closure helpers (oracle)

def low_sha(P, s0: int):
    """From a part of a finished degree-(D+1) closure: (sha, rank) of the reduced echelon basis of
    V_{D+1} n M_<=D in the degree-D column order, i.e. the rows with pivot weight <= D restricted to the
    low columns (they vanish on the top columns; asserted).  Same hash definition as ``Part.sha``."""
    P.compact()
    low = P.low_rows()
    piv = P.piv[low] - s0
    assert np.all(piv >= 0)
    o = np.argsort(piv)
    free = P.free_cols()
    fl = free[free >= s0] - s0
    ps0 = int(np.searchsorted(P.pc2col, s0))
    h = hashlib.sha256()
    h.update(piv[o].astype(np.int64).tobytes())
    h.update(fl.astype(np.int64).tobytes())
    dt = np.uint16 if P.p < 65536 else np.uint32
    ch = max(1, (1 << 24) // max(1, P.ncols))
    for a in range(0, len(o), ch):
        rows = low[o[a:a + ch]]
        blk = np.zeros((len(rows), P.ncols), dtype=np.float32)
        for ti, idx, off in P._groups(rows):
            blk[idx] = P._rows_host(ti, off, slice(None))
        assert not np.any(blk[:, :ps0] != 0), 'low row with support on top columns'
        h.update(np.ascontiguousarray(blk[:, ps0:].astype(dt)).tobytes())
    return h.hexdigest(), int(len(low))


def full_space_sha(c: int) -> str:
    """``Part.sha`` of the full space on c columns (all columns pivots, no free column)."""
    h = hashlib.sha256()
    h.update(np.arange(c, dtype=np.int64).tobytes())
    h.update(np.zeros(0, dtype=np.int64).tobytes())
    return h.hexdigest()


def direct_record(g: dict, k: int, p: int, D1: int, device: str = 'cpu', **kw) -> dict:
    """direct closure V_{D1} with the engine and the quantities the lift is compared with."""
    eng = B.Engine(g, k, p, D1, device=device, log=False, **kw)
    st = eng.run()
    assert st == 'done' and eng.closed, f'direct closure not closed: {st}'
    prof, _ = eng.profile()
    lows = [low_sha(P, P.s0) for P in eng.parts]
    return dict(dim=eng.rank(), profile=prof, contains_one=bool(eng.contains_one), low_basis_sha256=[x[0] for x in lows],
                low_part_ranks=[x[1] for x in lows], part_ranks=[P.rank for P in eng.parts],
                seconds=round(eng.seconds(), 1), madds=eng.total_madds(), closed=bool(eng.closed))


def compare(lift_rec: dict, direct: dict) -> dict:
    keys = dict(dim=lift_rec['dim_VD1'] == direct['dim'], profile=lift_rec['profile'] == direct['profile'],
                contains_one=lift_rec['contains_one'] == direct['contains_one'],
                low_sha=lift_rec['low_basis_sha256'] == direct['low_basis_sha256'])
    if 'part_ranks' in direct:
        keys['part_dims'] = lift_rec['part_dims_VD1'] == direct['part_ranks']
    keys['all'] = all(keys.values())
    return keys


# =============================================================================== plan

def plan(k: int, D: int, free, free_top, work_gb: float, dev_gb=None, basis_on: str = 'host') -> dict:
    """predicted sizes for the lift D -> D+1 from the free-column counts of V_D per part (``free``) and the
    number of those of weight D (``free_top``).  Upper bounds unless stated."""
    from math import comb
    L = B.Layout(k, D, True)
    cs = [P['c'] for P in L.parts]
    wD = [int(np.sum(P['w'] == D)) for P in L.parts]
    rank = [cs[q] - free[q] for q in (0, 1)]
    top_full = comb(2 * k, D + 1)
    ndiag = comb(k, (D + 1) // 2) if (D + 1) % 2 == 0 else 0
    top = [(top_full + ndiag) // 2, (top_full - ndiag) // 2]
    ft_orbits = free_top[0] + free_top[1]
    sp = [min(top[q], ft_orbits * (2 * k - D)) for q in (0, 1)]
    T = [sp[q] + free[q] for q in (0, 1)]
    nO = min(wD[0] - free_top[0], wD[1] - free_top[1] + (comb(k, D // 2) if D % 2 == 0 else 0))
    units = nO * k
    images = 2 * (2 * units + (free_top[0] + free_top[1]) * k)
    madds_worst = sum(units * (free[0] + free[1]) * T[q] for q in (0, 1))
    fol = [[free[p] * T[q] * 4 for p in (0, 1)] for q in (0, 1)]
    cache_all = [k * sum(fol[q]) for q in (0, 1)]
    cache_group = [(D + 1) * sum(fol[q]) for q in (0, 1)]
    basis = [rank[q] * free[q] * 4 for q in (0, 1)]
    kworst = [T[q] * T[q] for q in (0, 1)]                      # max_r r (T - r) 4 bytes
    work = work_gb * (1 << 30)
    nu = max(64, int(0.7 * work // (80 * max(T) + 32 * max(free))))
    gb = 1 << 30
    out = dict(k=k, D=D, part_columns_D=cs, rank_D=rank, free=list(free), free_top=list(free_top),
               top_orbits=top, special_top_orbits_max=sp, regular_top_orbits_min=[top[q] - sp[q] for q in (0, 1)],
               dim_T_max=T, regular_weight_D_orbits=nO, units_per_target=units, images_total=images,
               madds_gemm_worst=madds_worst, madds_note='rows x free_p x dim T_q with K empty; the real count is rows x free_p x '
               '(physical columns of K_q), which falls as K_q grows and is 0 once K_q = T_q',
               follower_gb_each_max=[[round(x / gb, 2) for x in fq] for fq in fol],
               cache_all_followers_gb_max=[round(x / gb, 1) for x in cache_all],
               cache_one_J_group_gb_max=[round(x / gb, 1) for x in cache_group],
               basis_gb=[round(x / gb, 1) for x in basis], K_worst_gb=[round(x / gb, 2) for x in kworst],
               work_gb=work_gb, batch_units=nu, host_gb_basis_on_host=round(sum(basis) / gb + 6, 1))
    if dev_gb:
        dev = dev_gb * gb
        b = sum(basis) if basis_on == 'device' else max(basis)      # host mode: one part comes to the device at the end only
        cache = dev - work - max(kworst) - 1.5 * gb - (sum(basis) if basis_on == 'device' else 0)
        out.update(dev_gb=dev_gb, basis_on=basis_on, suggested_cache_gb=round(cache / gb, 1),
                   cache_holds_one_J_group=bool(cache >= max(cache_group)), cache_holds_all=bool(cache >= max(cache_all)),
                   final_step_device_gb=round((b + work) / gb, 1),
                   predicted_device_peak_gb=round((min(cache, max(cache_all)) + work + max(kworst) + 1.5 * gb
                                                   + (sum(basis) if basis_on == 'device' else 0)) / gb, 1))
    return out
