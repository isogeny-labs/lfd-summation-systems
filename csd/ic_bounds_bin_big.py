"""Closure engine for the Boolean (characteristic 2) degree-D closure V_D at large size: numpy (cpu) or
cupy (gpu) behind one code path, with the ``fake-gpu`` stand-in of ``ic_bounds_big`` for local validation.

DEFINITION (identical to csd/ic_boolean_closure.py, ``boolean_degree_closure``; checked in
experiments/ic-bounds-binary-deep/validation).  Variables x_0..x_{r-1} (r = 2k), Boolean ring x_i^2 = x_i.
Columns = squarefree monomials (bit masks) of degree <= D ordered degree DESCENDING, ties by mask DESCENDING
(``monomial_masks``: sorted by (popcount, mask), reverse=True).  V_D = the least subspace of span(columns)
that contains every input equation of degree <= D and is closed under f -> x_i f for every variable i and
every f in V_D of degree < D (the multiplication S -> S u {i} is not injective on monomials: two source
monomials landing on one target add mod 2).  With the order above the pivot (leading column) of a row is its
highest degree monomial, so "rows of degree < D" = rows whose pivot column has degree < D = "low" rows.
Fall at D (D >= 1): dim(V_D cap M_<=D-1) > dim V_(D-1), i.e. the number of rows of the reduced basis with
pivot degree <= D-1 exceeds dim V_(D-1) (this is ``low_degree_dimension(D-1) > spaces[-1].rank`` of
ic_local_certificate.local_profile with the ordinary rank written out).  1 in V_D (a row with pivot = the
constant column) makes V_D everything (unit ideal): LFD = that D.

STORAGE.  The reduced row echelon basis (RREF) for the fixed column order: row i = e_{piv[i]} + B[i] on the
free (non-pivot) columns.  B is stored BIT-PACKED (uint64 words, 1 bit per entry) in row tiles: a float32 or
uint8 row would cost 32 or 8 times more memory, and the worst case r x (N - r) is N^2/4 bits = 41 GB at
N = 1,149,017 columns (2k = 32, D = 6), against 1.3 TB as float32.  Pivot columns are deleted lazily (kept
exactly zero) and compacted tile by tile.  The matrix products are done on 0/1 float32 matrices that are
unpacked from the words block by block: every partial sum is an integer below 2^24 (asserted: the inner
dimension is at most the number of columns), hence exact in float32 in any summation order; the parity of the
result is packed back and XORed in.  The canonical hash ``rref_sha256`` depends only on the subspace and the
column order, not on the physical layout, tiles, batches or device.

ALGORITHM (mirror of csd/ic_bounds_big.py with the swap split removed: there is none in characteristic 2).
  * ``rref_block``: exact RREF of a block of rows by panels of 64 columns (Gauss-Jordan on the 64-bit panel
    value on the host, then one float32 product to apply it).
  * ``absorb``: new rows are added to the basis by one product Cold @ New per tile (old rows that meet the new
    pivots are updated in place, the new pivot columns become zero in every old row).
  * closure: low rows are processed once each, at their state when their batch starts; for every variable i
    the shifted block Y_i (m x (columns with bit i)) is reduced against the basis,
    R = Y_free + Y_piv @ B using only the pivot columns T on which the batch is nonzero (sparse-aware), and the
    nonzero residual rows are absorbed.
  * certificate: deterministic full sweep over the final basis: every variable times every row of pivot degree
    < D reduces to zero (and every input reduces to zero).  ``closed`` is true only after one complete clean
    sweep over an unchanged basis; a nonzero residual is absorbed (repair), the closure continues and a new
    sweep starts (``--no-repair`` stops instead: the negative control).
"""
from __future__ import annotations

import functools
import hashlib
import json
import math
import os
import resource
import sys
import time
from collections import deque

import numpy as np

try:
    from csd import ic_bounds_big as BG
except Exception:  # pragma: no cover - remote layout: both files in one directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import ic_bounds_big as BG

ENGINE_VERSION = 'ic_bounds_bin_big/1'
F32 = np.float32
MAXINNER = 1 << 24
WORK_BYTES = 1 << 30      # budget of the transient (rows x column block) matrices of one product; set by Engine
DV = BG.Device('cpu')
xp = np


def use_device(name: str):
    global DV, xp
    DV = BG.use_device(name)
    xp = DV.xp
    return DV


def file_sha256(path) -> str:
    return BG.file_sha256(path)


ENGINE_SHA256 = file_sha256(__file__)


def maxrss_mb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


# =============================================================================== bit helpers (device generic)

def _sh8():
    return DV.to_dev(np.arange(8, dtype=np.uint8))


def unpack(words):
    """(m, W) uint64 -> (m, 64 W) uint8 of 0/1; bit j of word w is column 64 w + j (little endian)."""
    m, W = words.shape
    words = xp.ascontiguousarray(words)
    b = words.view(np.uint8).reshape(m, W * 8)
    bits = (b[:, :, None] >> _sh8()[None, None, :]) & 1
    return bits.reshape(m, W * 64)


def pack(bits):
    """(m, 64 W) uint8 of 0/1 -> (m, W) uint64 (inverse of ``unpack``)."""
    m, n = bits.shape
    assert n % 64 == 0
    x = bits.reshape(m, n // 8, 8) << _sh8()[None, None, :]
    by = x.sum(axis=2, dtype=np.uint8)
    return xp.ascontiguousarray(by).view(np.uint64).reshape(m, n // 64)


def parity(counts):
    """float32 array of exact non-negative integers < 2^24 -> uint8 parity."""
    return (counts.astype(np.int32) & 1).astype(np.uint8)


class Flops:
    upd = 0.0
    res = 0.0
    gen = 0.0
    rounds = 0
    products = 0

    @classmethod
    def reset(cls):
        cls.upd = cls.res = cls.gen = 0.0
        cls.rounds = cls.products = 0


def xor_apply(rows, Cf, P, w_start=0, cb=1 << 16, kind='upd'):
    """rows ^= (Cf @ P) mod 2, in place, column block by column block.  rows (m, W) uint64, Cf (m, p) float32
    of 0/1, P (p, W) uint64.  Exact: the inner dimension p is < 2^24 (asserted)."""
    m, W = rows.shape
    p = Cf.shape[1]
    assert p < MAXINNER
    if m == 0 or p == 0:
        return rows
    step = max(1, min(cb, max(64, int(WORK_BYTES // (13 * m)) // 64 * 64)) // 64)
    for w0 in range(w_start, W, step):
        w1 = min(W, w0 + step)
        Pu = unpack(P[:, w0:w1]).astype(F32)
        acc = Cf @ Pu
        rows[:, w0:w1] ^= pack(parity(acc))
        Flops.products += 1
    f = 2.0 * m * p * 64 * (W - w_start)
    if kind == 'upd':
        Flops.upd += f
    else:
        Flops.res += f
    return rows


def panel_gauss(panel):
    """Gauss-Jordan on the 64-bit panel values (host).  Returns (picks, pivot bits, combos): picks = row indices
    forming a maximal independent subset (greedy by ascending column), pivot bit of each, and for each pick
    the bit mask (over pick order) of the original pick rows whose XOR is its reduced panel value."""
    work = panel.copy()
    cm = np.zeros(len(panel), dtype=np.uint64)
    taken = np.zeros(len(panel), dtype=bool)
    picks, pcs = [], []
    one = np.uint64(1)
    for c in range(64):
        has = ((work >> np.uint64(c)) & one).astype(bool)
        idx = np.flatnonzero(has & ~taken)
        if idx.size == 0:
            continue
        j = int(idx[0])
        taken[j] = True
        cm[j] ^= one << np.uint64(len(picks))
        picks.append(j)
        pcs.append(c)
        has[j] = False
        work[has] ^= work[j]
        cm[has] ^= cm[j]
    return np.array(picks, dtype=np.int64), np.array(pcs, dtype=np.int64), cm[np.array(picks, dtype=np.int64)]


def rref_block(R, cb=1 << 16):
    """Exact RREF of the rows of R (device, (m, W) uint64; consumed).  Returns (rows (t, W), pivot bit
    positions (t,) host int64 ascending).  Each returned row has its pivot as first set bit and zero at every
    other pivot position."""
    m, W = R.shape
    Acc = xp.zeros((max(m, 1), W), dtype=np.uint64)
    t = 0
    piv = []
    rem = R
    while rem.shape[0]:
        keep = xp.flatnonzero(rem.any(axis=1))
        if keep.shape[0] == 0:
            break
        if keep.shape[0] != rem.shape[0]:
            rem = rem[keep]
        first = DV.to_host(xp.flatnonzero(rem.any(axis=0))[:1])
        w0 = int(first[0])
        panel = DV.to_host(rem[:, w0]).astype(np.uint64)
        picks, pcs, cm = panel_gauss(panel)
        p = len(picks)
        Rp = rem[DV.to_dev(picks)]
        Tm = ((cm[:, None] >> np.arange(p, dtype=np.uint64)[None, :]) & np.uint64(1)).astype(F32)
        Tm_d = DV.to_dev(Tm)
        Pk = xp.zeros((p, W), dtype=np.uint64)
        xor_apply(Pk, Tm_d, Rp, cb=cb)
        pcs_d = DV.to_dev(pcs.astype(np.uint64))
        others = np.ones(rem.shape[0], dtype=bool)
        others[picks] = False
        oi = np.flatnonzero(others)
        if oi.size:
            rem = rem[DV.to_dev(oi)]
            C = ((rem[:, w0][:, None] >> pcs_d[None, :]) & 1).astype(F32)
            xor_apply(rem, C, Pk, w_start=w0, cb=cb)
        else:
            rem = xp.zeros((0, W), dtype=np.uint64)
        if t:
            C2 = ((Acc[:t, w0][:, None] >> pcs_d[None, :]) & 1).astype(F32)
            xor_apply(Acc[:t], C2, Pk, w_start=w0, cb=cb)
        Acc[t:t + p] = Pk
        t += p
        piv.extend((64 * w0 + pcs).tolist())
        Flops.rounds += 1
    return Acc[:t], np.array(piv, dtype=np.int64)


# =============================================================================== columns

@functools.lru_cache(maxsize=None)
def _combos(r: int, d: int):
    """ascending uint64 array of the masks of popcount d in r bits."""
    if d < 0 or d > r:
        return np.zeros(0, dtype=np.uint64)
    if d == 0:
        return np.zeros(1, dtype=np.uint64)
    if d == r:
        return np.array([(1 << r) - 1], dtype=np.uint64)
    a = _combos(r - 1, d)
    b = _combos(r - 1, d - 1) | np.uint64(1 << (r - 1))
    return np.concatenate([a, b])


def popcount(a):
    a = np.asarray(a, dtype=np.uint64)
    if hasattr(np, 'bitwise_count'):
        return np.bitwise_count(a).astype(np.int64)
    tab = np.array([bin(i).count('1') for i in range(256)], dtype=np.int64)
    b = a.reshape(-1, 1).view(np.uint8).reshape(-1, 8)
    return tab[b].sum(axis=1).reshape(a.shape)


def count_columns(r: int, D: int) -> int:
    return sum(math.comb(r, i) for i in range(0, min(r, D) + 1))


class Cols:
    """Squarefree monomials of degree <= D in r variables, degree descending then mask descending."""

    def __init__(self, r: int, D: int):
        assert 1 <= r <= 60 and D >= 0
        self.r, self.D = r, D
        D = min(D, r)
        self.dstart, self.cnt = {}, {}
        parts, pos = [], 0
        for d in range(D, -1, -1):
            a = _combos(r, d)
            self.dstart[d], self.cnt[d] = pos, len(a)
            pos += len(a)
            parts.append(a[::-1])
        self.masks = np.concatenate(parts)
        self.N = pos
        self.ntop = len(_combos(r, D))
        self.Nlow = self.N - self.ntop
        self.deg = popcount(self.masks)
        self._maps = {}
        self._maps_d = {}

    def index(self, m):
        m = np.asarray(m, dtype=np.uint64)
        dg = popcount(m)
        out = np.empty(len(m), dtype=np.int64)
        for d in np.unique(dg):
            d = int(d)
            sel = dg == d
            ps = np.searchsorted(_combos(self.r, d), m[sel])
            out[sel] = self.dstart[d] + self.cnt[d] - 1 - ps
        return out

    def shiftmap(self, i: int):
        if i not in self._maps:
            bit = np.uint64(1 << i)
            tgt = np.flatnonzero((self.masks & bit) != 0)
            a = np.where(tgt >= self.ntop, tgt - self.ntop, self.Nlow)
            b = self.index(self.masks[tgt] ^ bit) - self.ntop
            self._maps[i] = (tgt.astype(np.int64), a.astype(np.int64), b.astype(np.int64))
        return self._maps[i]

    def shiftmap_dev(self, i: int):
        if i not in self._maps_d:
            self._maps_d[i] = tuple(DV.to_dev(x) for x in self.shiftmap(i))
        return self._maps_d[i]

    def drop_dev(self):
        self._maps_d = {}


# =============================================================================== the engine

class Engine:
    def __init__(self, equations, r: int, D: int, points=None, device='cpu', batch=256, tile_rows=None, tile_mb=256.0,
                 work_gb=2.0, col_block=None, compact_frac=0.25, tiles_on='device', withhold=None, repair=True,
                 sweep=True, status=None, log_seconds=30.0, quiet=True):
        use_device(device)
        Flops.reset()
        global WORK_BYTES
        WORK_BYTES = int(work_gb * (1 << 30))
        self.device, self.r, self.D = device, r, D
        self.eqs = _clean(equations, r)
        self.points = points
        self.cols = Cols(r, D)
        c = self.cols
        self.N, self.ntop, self.Nlow = c.N, c.ntop, c.Nlow
        assert self.N < MAXINNER, 'float32 exactness needs fewer than 2^24 columns'
        self.batch = int(batch)
        self.compact_frac = compact_frac
        self.tiles_on = tiles_on
        self.withhold = withhold
        self.repair = repair
        self.do_sweep = sweep
        self.status_path = status
        self.log_seconds = log_seconds
        self.quiet = quiet
        self.phys = np.arange(self.N, dtype=np.int64)
        self.W = max(1, -(-self.N // 64))
        if tile_rows is None:
            tile_rows = max(8, int(tile_mb * (1 << 20) // (self.W * 8)))
        self.TR = max(8, min(int(tile_rows), self.N))
        if col_block is None:
            col_block = max(64, min(self.W * 64, int(work_gb * (1 << 30) / (max(1, self.batch) * 24)) // 64 * 64))
        self.cb = int(col_block)
        self.pos = np.arange(self.N, dtype=np.int64)
        self.pos_d = DV.to_dev(self.pos)
        self.rc = np.full(self.N, -1, dtype=np.int64)
        self.rc_d = DV.to_dev(self.rc)
        self.row_piv = np.zeros(self.N + 1, dtype=np.int64)
        self.nrows = 0
        self.tiles = []
        self.pending = deque()
        self.unit = False
        self.closed = False
        self.sweeps = 0
        self.repairs = 0
        self.batches = 0
        self.shifted_rows = 0
        self.compactions = 0
        self.t0 = time.time()
        self._last_log = 0.0
        self.sweep_clean_rows = 0

    # ------------------------------------------------------------------ tiles
    def _tile(self, ti):
        t = self.tiles[ti]
        return t['a'] if t['dev'] else DV.to_dev(t['a'])

    def _store(self, ti, arr):
        t = self.tiles[ti]
        t['a'] = arr if t['dev'] else DV.to_host(arr)

    def _new_tile(self):
        dev = self.tiles_on == 'device'
        a = xp.zeros((self.TR, self.W), dtype=np.uint64) if dev else np.zeros((self.TR, self.W), dtype=np.uint64)
        self.tiles.append(dict(a=a, n=0, dev=dev))

    def get_rows(self, idx):
        idx = np.asarray(idx, dtype=np.int64)
        out = xp.zeros((len(idx), self.W), dtype=np.uint64)
        tis = idx // self.TR
        for ti in np.unique(tis):
            sel = np.flatnonzero(tis == ti)
            a = self._tile(int(ti))
            out[DV.to_dev(sel)] = a[DV.to_dev(idx[sel] % self.TR)]
        return out

    def tile_bytes(self):
        return sum(t['a'].nbytes for t in self.tiles)

    # ------------------------------------------------------------------ absorbing rows
    def _update_old(self, Acc, piv):
        wpos = piv >> 6
        bpos = (piv & 63).astype(np.uint64)
        wd, bd = DV.to_dev(wpos), DV.to_dev(bpos)
        ws = int(wpos.min())
        for ti, t in enumerate(self.tiles):
            n = t['n']
            if n == 0:
                continue
            a = self._tile(ti)
            cols = a[:n][:, wd]
            Cb = ((cols >> bd[None, :]) & 1).astype(F32)
            rows_hit = xp.flatnonzero(Cb.any(axis=1))
            nh = int(rows_hit.shape[0])
            if nh == 0:
                continue
            if nh < 0.7 * n:
                sub = a[rows_hit]
                xor_apply(sub, Cb[rows_hit], Acc, ws, self.cb)
                a[rows_hit] = sub
            else:
                xor_apply(a[:n], Cb, Acc, ws, self.cb)
            if not t['dev']:
                self._store(ti, a)

    def absorb(self, R):
        """Add the rows of R (device (m, W), already reduced against the basis: zero on every pivot column) to
        the basis.  Returns the number of new rows."""
        Acc, piv = rref_block(R, self.cb)
        t = len(piv)
        if t == 0:
            return 0
        self._update_old(Acc, piv)
        # B-form: clear the own pivot bit
        ri = DV.to_dev(np.arange(t, dtype=np.int64))
        wi = DV.to_dev(piv >> 6)
        bit = DV.to_dev(np.uint64(1) << (piv & 63).astype(np.uint64))
        Acc[ri, wi] = Acc[ri, wi] ^ bit
        glob = self.phys[piv]
        first = self.nrows
        done = 0
        while done < t:
            if not self.tiles or self.tiles[-1]['n'] == self.TR:
                self._new_tile()
            ti = len(self.tiles) - 1
            tl = self.tiles[ti]
            take = min(t - done, self.TR - tl['n'])
            a = self._tile(ti)
            a[tl['n']:tl['n'] + take] = Acc[done:done + take]
            if not tl['dev']:
                self._store(ti, a)
            tl['n'] += take
            done += take
        self.row_piv[first:first + t] = glob
        self.rc[glob] = np.arange(first, first + t)
        self.rc_d[DV.to_dev(glob)] = DV.to_dev(np.arange(first, first + t, dtype=np.int64))
        self.nrows += t
        low = [first + j for j in range(t) if glob[j] >= self.ntop]
        self.pending.extend(low)
        if np.any(glob == self.N - 1):
            self.unit = True
        self.maybe_compact()
        return t

    # ------------------------------------------------------------------ compaction
    def maybe_compact(self, force=False):
        dead = len(self.phys) - (self.N - self.nrows)
        if dead <= 0:
            return
        if force or (dead >= self.compact_frac * len(self.phys) and dead >= 64):
            self.compact()

    def compact(self):
        keep = np.ones(len(self.phys), dtype=bool)
        pp = self.pos[self.row_piv[:self.nrows]]
        keep[pp[pp >= 0]] = False
        newphys = self.phys[keep]
        P2 = len(newphys)
        W2 = max(1, -(-P2 // 64))
        kd = DV.to_dev(np.flatnonzero(keep).astype(np.int64))
        oldW = self.W
        rc = max(1, int(((256 << 20) // max(1, oldW * 64))))
        for ti, t in enumerate(self.tiles):
            a = self._tile(ti)
            new = xp.zeros((self.TR, W2), dtype=np.uint64)
            for s in range(0, t['n'], rc):
                e = min(t['n'], s + rc)
                bits = unpack(a[s:e])
                padb = xp.zeros((e - s, W2 * 64), dtype=np.uint8)
                padb[:, :P2] = bits[:, kd]
                new[s:e] = pack(padb)
            del a
            t['a'] = new if t['dev'] else DV.to_host(new)
        self.phys, self.W = newphys, W2
        self.pos = np.full(self.N, -1, dtype=np.int64)
        self.pos[newphys] = np.arange(P2)
        self.pos_d = DV.to_dev(self.pos)
        self.compactions += 1

    # ------------------------------------------------------------------ residuals
    def residual(self, Y, tgt_d):
        """Reduce the rows of Y (device (m, nt) uint8, columns = global columns tgt) against the basis."""
        m, nt = Y.shape
        rows_of = self.rc_d[tgt_d]
        ispiv = rows_of >= 0
        fidx = xp.flatnonzero(~ispiv)
        Rd = xp.zeros((m, self.W * 64), dtype=np.uint8)
        if fidx.shape[0]:
            Rd[:, self.pos_d[tgt_d[fidx]]] = Y[:, fidx]
        R = pack(Rd)
        del Rd
        T = xp.flatnonzero(Y.any(axis=0) & ispiv)
        nT = int(T.shape[0])
        if nT:
            assert nT < MAXINNER
            Yp = Y[:, T].astype(F32)
            rows_T = rows_of[T]
            tile_T = DV.to_host(rows_T // self.TR)
            offs_all = rows_T % self.TR
            for ti in np.unique(tile_T):
                sel = DV.to_dev(np.flatnonzero(tile_T == ti))
                off = offs_all[sel]
                Yt = Yp[:, sel]
                a = self._tile(int(ti))
                nsel = int(sel.shape[0])
                step = max(1, min(self.cb, max(64, int(WORK_BYTES // (13 * max(m, nsel))) // 64 * 64)) // 64)
                for w0 in range(0, self.W, step):
                    w1 = min(self.W, w0 + step)
                    Bu = unpack(a[off, w0:w1]).astype(F32)
                    R[:, w0:w1] ^= pack(parity(Yt @ Bu))
                    Flops.products += 1
                Flops.res += 2.0 * m * int(sel.shape[0]) * 64 * self.W
        return R

    def dense_low(self, rows):
        R = self.get_rows(rows)
        bits = unpack(R)
        q0 = int(np.searchsorted(self.phys, self.ntop))
        P = len(self.phys)
        m = len(rows)
        Sg = xp.zeros((m, self.Nlow + 1), dtype=np.uint8)
        if P > q0:
            Sg[:, DV.to_dev(self.phys[q0:] - self.ntop)] = bits[:, q0:P]
        Sg[DV.to_dev(np.arange(m, dtype=np.int64)), DV.to_dev(self.row_piv[rows] - self.ntop)] = 1
        return Sg

    def process_batch(self, rows, sweep=False):
        """Shifts of the given low rows by every variable, reduced against the current basis.  Returns the
        number of variables whose block gave a nonzero residual."""
        rows = np.asarray(rows, dtype=np.int64)
        Sg = self.dense_low(rows)
        m = len(rows)
        wh = None
        if self.withhold is not None and not sweep:
            hit = np.flatnonzero(rows == self.withhold)
            wh = int(hit[0]) if hit.size else None
        bad = 0
        for i in range(self.r):
            tgt_d, a_d, b_d = self.cols.shiftmap_dev(i)
            Y = Sg[:, a_d] ^ Sg[:, b_d]
            if wh is not None:
                Y[wh] = 0
            if int(tgt_d.shape[0]) == 0:
                continue
            R = self.residual(Y, tgt_d)
            nz = bool(R.any())
            if nz:
                bad += 1
                if not sweep or self.repair:
                    self.absorb(R)
                    if sweep:
                        self.repairs += 1
                    if self.unit:
                        break
            self.shifted_rows += m
            self._tick()
        return bad

    # ------------------------------------------------------------------ driver
    def init_inputs(self):
        eq = [e for e in self.eqs if e and max(int(m).bit_count() for m in e) <= self.D]
        rows = np.zeros((len(eq), self.W * 64), dtype=np.uint8)
        for j, e in enumerate(eq):
            rows[j, self.cols.index(np.array(e, dtype=np.uint64))] ^= 1
        # an equation that is the empty list (0 = 0) carries no information; constant 1 is the mask 0
        if len(eq):
            self.absorb(pack(DV.to_dev(rows)))
        self.n_inputs = len(eq)

    def closure_loop(self):
        while self.pending and not self.unit:
            take = [self.pending.popleft() for _ in range(min(self.batch, len(self.pending)))]
            self.process_batch(take, sweep=False)
            self.batches += 1

    def check_inputs(self):
        eq = [e for e in self.eqs if e and max(int(m).bit_count() for m in e) <= self.D]
        if not eq:
            return True
        Y = np.zeros((len(eq), self.N), dtype=np.uint8)
        for j, e in enumerate(eq):
            Y[j, self.cols.index(np.array(e, dtype=np.uint64))] ^= 1
        tgt_d = DV.to_dev(np.arange(self.N, dtype=np.int64))
        R = self.residual(DV.to_dev(Y), tgt_d)
        return not bool(R.any())

    def sweep(self):
        """One deterministic pass over the final low rows.  Returns True iff no residual was nonzero."""
        self.sweeps += 1
        before = self.nrows
        ok = self.check_inputs()
        low = np.flatnonzero(self.row_piv[:self.nrows] >= self.ntop)
        bad = 0
        for s in range(0, len(low), self.batch):
            bad += self.process_batch(low[s:s + self.batch], sweep=True)
            if self.unit:
                break
        self.sweep_clean_rows = len(low)
        return ok and bad == 0 and self.nrows == before

    def run(self):
        self.init_inputs()
        while True:
            self.closure_loop()
            if self.unit:
                break
            if not self.do_sweep:
                break
            if self.sweep():
                self.closed = True
                break
            if not self.repair:
                break
            self.closure_loop()
            if self.unit:
                break
            # loop: another sweep over the repaired basis
        if self.unit:
            self.closed = True
        return self.record()

    # ------------------------------------------------------------------ checkpoint (one, at the end of a closure)
    def save_checkpoint(self, path, extra=None):
        """Write the basis: tiles as .npy, then the manifest last, in a temporary directory renamed atomically."""
        import shutil
        tmp = path + '.tmp'
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        for ti, t in enumerate(self.tiles):
            a = DV.to_host(self._tile(ti)) if t['dev'] else t['a']
            np.save(os.path.join(tmp, 'tile%05d.npy' % ti), a[:t['n']])
        np.save(os.path.join(tmp, 'phys.npy'), self.phys)
        np.save(os.path.join(tmp, 'row_piv.npy'), self.row_piv[:self.nrows])
        man = dict(engine=ENGINE_VERSION, r=self.r, D=self.D, N=self.N, nrows=int(self.nrows), TR=self.TR,
                   tiles=[t['n'] for t in self.tiles], unit=bool(self.unit), closed=bool(self.closed),
                   rref_sha256=None if self.unit else self.rref_sha256(), extra=extra or {})
        with open(os.path.join(tmp, 'manifest.json'), 'w') as f:
            json.dump(man, f)
            f.flush()
            os.fsync(f.fileno())
        shutil.rmtree(path, ignore_errors=True)
        os.replace(tmp, path)
        return man

    def load_checkpoint(self, path):
        with open(os.path.join(path, 'manifest.json')) as f:
            man = json.load(f)
        assert (man['r'], man['D'], man['N']) == (self.r, self.D, self.N), 'checkpoint is for another system size'
        self.TR = man['TR']
        self.phys = np.load(os.path.join(path, 'phys.npy'))
        self.W = max(1, -(-len(self.phys) // 64))
        self.pos = np.full(self.N, -1, dtype=np.int64)
        self.pos[self.phys] = np.arange(len(self.phys))
        self.pos_d = DV.to_dev(self.pos)
        self.nrows = man['nrows']
        self.row_piv[:self.nrows] = np.load(os.path.join(path, 'row_piv.npy'))
        self.rc[:] = -1
        self.rc[self.row_piv[:self.nrows]] = np.arange(self.nrows)
        self.rc_d = DV.to_dev(self.rc)
        self.tiles = []
        for ti, n in enumerate(man['tiles']):
            self._new_tile()
            a = np.load(os.path.join(path, 'tile%05d.npy' % ti))
            full = np.zeros((self.TR, self.W), dtype=np.uint64)
            full[:n] = a
            self.tiles[-1]['a'] = DV.to_dev(full) if self.tiles[-1]['dev'] else full
            self.tiles[-1]['n'] = n
        self.unit = man['unit']
        self.pending = deque()
        return man

    def resume_certificate(self, path):
        """Load a closure checkpoint and run the certificate sweep over it (no closure)."""
        self.load_checkpoint(path)
        self.t0 = time.time()
        if not self.unit:
            self.closed = self.sweep()
            if not self.closed and self.repair:
                self.closure_loop()
                self.closed = self.sweep()
        else:
            self.closed = True
        return self.record()

    # ------------------------------------------------------------------ outputs
    def pivots(self):
        return np.sort(self.row_piv[:self.nrows])

    def profile(self):
        """p_j = dim(V_D cap M_<=j) for j = 0..D."""
        dg = self.cols.deg[self.pivots()]
        return [int((dg <= j).sum()) for j in range(self.D + 1)]

    def rref_sha256(self):
        return canonical_hash(self.r, self.D, self.N, self.pivots(), self._free_rows_iter)

    def _free_rows_iter(self, free_cols):
        """yield the packed (little endian bits) free-column part of each row in pivot order."""
        order = np.argsort(self.row_piv[:self.nrows], kind='stable')
        fpos = DV.to_dev(self.pos[free_cols])
        nf = len(free_cols)
        chunk = max(1, int((64 << 20) // max(1, self.W * 64)))
        for s in range(0, len(order), chunk):
            idx = order[s:s + chunk]
            R = self.get_rows(idx)
            bits = DV.to_host(unpack(R)[:, fpos]) if nf else np.zeros((len(idx), 0), dtype=np.uint8)
            for row in bits:
                yield np.packbits(row, bitorder='little').tobytes()

    def record(self):
        piv = self.pivots()
        if self.unit:
            prof = [count_columns(self.r, j) for j in range(self.D + 1)]
            rank = self.N
            sha = canonical_hash_unit(self.r, self.D, self.N)
        else:
            prof = self.profile()
            rank = int(self.nrows)
            sha = self.rref_sha256()
        return dict(engine=ENGINE_VERSION, engine_sha256=ENGINE_SHA256, device=self.device, r=self.r, D=self.D,
                    columns=self.N, rank=rank, codim=self.N - rank, profile=prof, contains_one=bool(self.unit),
                    closed=bool(self.closed), unit=bool(self.unit), rref_sha256=sha, inputs=getattr(self, 'n_inputs', 0),
                    sweeps=self.sweeps, repairs=self.repairs, batches=self.batches, shifted_rows=self.shifted_rows,
                    compactions=self.compactions, low_rows=int((piv >= self.ntop).sum()) if not self.unit else None,
                    flops_update=Flops.upd, flops_residual=Flops.res, rref_rounds=Flops.rounds, products=Flops.products,
                    tile_rows=self.TR, col_block=self.cb, batch=self.batch, tiles=len(self.tiles),
                    tile_bytes=self.tile_bytes(), seconds=round(time.time() - self.t0, 3), peak_rss_mb=maxrss_mb(),
                    device_peak_bytes=int(DV.pool_bytes()) if DV.is_dev else 0)

    def _tick(self):
        now = time.time()
        if self.status_path and now - self._last_log >= self.log_seconds:
            self._last_log = now
            snap = dict(time=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), D=self.D, rank=int(self.nrows),
                        columns=self.N, pending=len(self.pending), batches=self.batches, sweeps=self.sweeps,
                        repairs=self.repairs, seconds=round(now - self.t0, 1), flops_update=Flops.upd,
                        flops_residual=Flops.res, tile_gb=round(self.tile_bytes() / 2**30, 2),
                        rss_mb=maxrss_mb(), device_used_gb=round(DV.mem_info()[0] / 2**30, 2))
            tmp = self.status_path + '.tmp'
            with open(tmp, 'w') as f:
                json.dump(snap, f)
            os.replace(tmp, self.status_path)
            if not self.quiet:
                print(json.dumps(snap), flush=True)


def _clean(equations, r):
    out = []
    for eq in equations:
        par = set()
        for m in eq:
            m = int(m)
            if not 0 <= m < 1 << r:
                raise ValueError('monomial outside the variable set')
            par ^= {m}
        out.append(tuple(sorted(par)))
    return out


# =============================================================================== canonical hash

def canonical_hash(r, D, N, pivots, rows_iter_fn):
    """sha256 of (r, D, N, rank, pivot columns ascending, free columns ascending, each row in pivot order as the
    little-endian packed bits of its entries on the free columns).  The row echelon basis of a subspace in
    reduced form is unique, so this identifies the subspace for the fixed column order."""
    h = hashlib.sha256()
    h.update(json.dumps(dict(fmt='bin-rref/1', r=r, D=D, N=N, rank=len(pivots))).encode())
    h.update(np.asarray(pivots, dtype='<u4').tobytes())
    mask = np.ones(N, dtype=bool)
    mask[pivots] = False
    free = np.flatnonzero(mask)
    h.update(np.asarray(free, dtype='<u4').tobytes())
    for b in rows_iter_fn(free):
        h.update(b)
    return h.hexdigest()


def canonical_hash_unit(r, D, N):
    return hashlib.sha256(json.dumps(dict(fmt='bin-rref/1', r=r, D=D, N=N, unit=True)).encode()).hexdigest()


def hash_from_dense(r, D, M):
    """canonical hash of the row space of a dense 0/1 matrix (host, small sizes) over the columns of Cols(r, D):
    independent reduction (used to compare with the existing engine's basis)."""
    M = np.array(M, dtype=np.uint8) & 1
    N = M.shape[1]
    # plain Gauss-Jordan on columns left to right
    A = M.copy()
    rk = 0
    pivs = []
    for c in range(N):
        if rk == A.shape[0]:
            break
        nz = np.flatnonzero(A[rk:, c])
        if nz.size == 0:
            continue
        j = rk + int(nz[0])
        if j != rk:
            A[[rk, j]] = A[[j, rk]]
        others = np.flatnonzero(A[:, c])
        others = others[others != rk]
        A[others] ^= A[rk]
        pivs.append(c)
        rk += 1
    A = A[:rk]
    pivs = np.array(pivs, dtype=np.int64)

    def rows_iter(free):
        for i in range(rk):
            yield np.packbits(A[i, free], bitorder='little').tobytes()
    return canonical_hash(r, D, N, pivs, rows_iter), rk, pivs


# =============================================================================== evaluation rank, local certificate

def eval_rank(points, r, D, device_batch=4096, cb=1 << 16, limit=None):
    """rank over F_2 of the evaluation of the squarefree monomials of degree < D on the points (a lower bound is
    enough to equal len(points); the loop stops as soon as the rank reaches the number of points)."""
    pts = np.array(sorted(points), dtype=np.uint64)
    Z = len(pts)
    if Z == 0:
        return 0
    Wz = -(-Z // 64)
    Acc = xp.zeros((0, Wz), dtype=np.uint64)
    rank = 0
    for d in range(0, min(D - 1, r) + 1):
        ms = _combos(r, d)
        for s in range(0, len(ms), device_batch):
            m = ms[s:s + device_batch]
            ev = ((pts[None, :] & m[:, None]) == m[:, None]).astype(np.uint8)
            padb = np.zeros((len(m), Wz * 64), dtype=np.uint8)
            padb[:, :Z] = ev
            blk = pack(DV.to_dev(padb))
            allr = xp.concatenate([Acc, blk], axis=0) if rank else blk
            Acc, piv = rref_block(allr, cb)
            rank = len(piv)
            if rank == Z:
                return rank
    return rank


def local_profile_bin(equations, r, points, max_degree, engine_kw=None, log=None):
    """Ladder D = 0, 1, ...: closure V_D, falls, unit test, local criterion (evaluation rank of M_<D on Z equals
    |Z| equals codim V_D).  Mirrors csd.ic_local_certificate.local_profile.  Returns the profile dict."""
    engine_kw = dict(engine_kw or {})
    points = list(points) if points is not None else []
    dims, falls, checks, recs = [], [], [], []
    status, lfd, cert = 'capped', None, None
    for D in range(0, max_degree + 1):
        eng = Engine(equations, r, D, points=points, **engine_kw)
        rec = eng.run()
        recs.append(rec)
        if log:
            log(rec)
        if not rec['closed']:
            status = 'not-closed'
            break
        N = rec['columns']
        dim = rec['rank']
        if D >= 1:
            p_low = rec['profile'][D - 1]
            if p_low > dims[D - 1]:
                falls.append(D)
        dims.append(dim)
        if rec['unit']:
            assert not points, 'unit ideal with a nonempty zero set'
            status, lfd, cert = 'exact-unit', D, D
            break
        if D > 0 and points:
            er = eval_rank(points, r, D, cb=eng.cb)
            codim = N - dim
            checks.append(dict(degree=D, evaluation_rank=er, codimension=codim))
            if er == len(points) == codim:
                status, lfd, cert = 'exact-proper-local', max(falls, default=0), D
                break
    return dict(status=status, exact_lfd=lfd, certificate_degree=cert, strict_falls=falls, quotient_ranks=dims,
                local_checks=checks, point_count=len(points), per_degree=recs)


# =============================================================================== plan

def plan(r: int, D: int, device='gpu', batch=1024, work_gb=20.0, tile_mb=1024.0, compact_frac=0.25,
         dev_gb=141.0, host_gb=250.0, nz=0, sgemm_tflops=30.0):
    """Predicted peak memory and an operation estimate for one closure at (r = 2k variables, degree D).
    Memory is exact bookkeeping of the buffers named below; the work estimate is an upper-bound model, see the
    notes (calibrated numbers are in experiments/ic-bounds-binary-deep/PLAN.md)."""
    N = count_columns(r, D)
    ntop = math.comb(r, D)
    Nlow = N - ntop
    # basis: rank x (free + lazily deleted) bits, worst case over the rank
    best, br = 0.0, 0
    step = max(1, N // 4000)
    for rk in range(0, N + 1, step):
        phys = min(N, (N - rk) * (1 + compact_frac) + 64) if rk else N
        by = rk * (-(-phys // 64)) * 8
        if by > best:
            best, br = by, rk
    nt = sum(math.comb(r - 1, d - 1) for d in range(1, min(D, r) + 1))     # columns containing a given variable
    W = -(-N // 64)
    TR = max(8, int(tile_mb * (1 << 20) // (W * 8)))
    m = batch
    cb = max(64, min(W * 64, int(work_gb * (1 << 30) / (m * 24)) // 64 * 64))
    bufs = dict(
        Sg=m * (Nlow + 1), Y=m * nt, Yp_f32=4 * m * min(nt, N), Rd=m * W * 64, R=m * W * 8,
        block_f32=4 * m * cb * 2, Bu_f32=4 * min(TR, N) * cb, tile_slack=TR * W * 8,
        shiftmaps=3 * 8 * nt * r // 4, panel_rows=3 * m * W * 8)
    block_terms = bufs['block_f32'] + bufs['Bu_f32']
    fixed = sum(v for kk, v in bufs.items() if kk not in ('block_f32', 'Bu_f32'))
    # the engine caps the transient (rows x column block) matrices of one product at WORK_BYTES = work_gb
    work = fixed + min(block_terms, work_gb * (1 << 30))
    compaction = (256 << 20) * 2
    ctx = 1.5 * (1 << 30)
    total = best + work + compaction + (ctx if device != 'cpu' else 0)
    # work: update products ~ N^3/3 (every new row meets every old row once), residual products: shifted rows x
    # nonzero pivot columns x free columns, bounded by (low rows x r) x nt x N; reported as a bound
    upd = N ** 3 / 3.0
    res_bound = Nlow * r * min(nt, N) * N * 2.0 / 8.0
    return dict(r=r, D=D, columns=N, top_columns=ntop, low_columns=Nlow, device=device,
                basis_worst_bytes=int(best), basis_worst_rank=br, basis_worst_gb=round(best / 2**30, 2),
                buffers_gb={k: round(v / 2**30, 3) for k, v in bufs.items()}, work_gb=round(work / 2**30, 2),
                predicted_peak_gb=round(total / 2**30, 1), fits_device=bool(total <= dev_gb * 2**30),
                fits_host_if_tiles_on_host=bool(best <= host_gb * 2**30),
                flops_update_model=upd, flops_residual_bound=res_bound,
                seconds_update_at_sgemm=round(upd / (sgemm_tflops * 1e12), 1),
                seconds_total_estimate=round(0.28 * N ** 3 / (sgemm_tflops * 1e12), 1), tile_rows=TR, col_block=cb,
                eval_points=nz,
                note='basis stored bit-packed (1 bit per entry); the worst case is rank N/2; float32 only in '
                     'transient blocks; fp32 matmul exact below 2^24 inner dimension')


# =============================================================================== selftest

def _rand_system(r, seed, n_eq=None, density=0.5):
    rng = np.random.RandomState(seed)
    n_eq = n_eq or max(2, r // 2)
    eqs = []
    for _ in range(n_eq):
        masks = set()
        for i in range(r):
            for j in range(i + 1, r):
                if rng.rand() < density / 3:
                    masks ^= {(1 << i) | (1 << j)}
        for i in range(r):
            if rng.rand() < 0.3:
                masks ^= {1 << i}
        if rng.rand() < 0.5:
            masks ^= {0}
        eqs.append(sorted(masks))
    return eqs


def selftest_kernels(device_name, seed=7):
    """exactness of the F_2 kernels against integer arithmetic; returns a dict of named booleans."""
    use_device(device_name)
    rng = np.random.RandomState(seed)
    out = {}
    # pack/unpack roundtrip
    w = rng.randint(0, 2**63, size=(5, 7), dtype=np.int64).astype(np.uint64) ^ (rng.randint(0, 2**63, size=(5, 7)).astype(np.uint64) << np.uint64(1))
    Wd = DV.to_dev(w)
    bits = DV.to_host(unpack(Wd))
    ref = np.array([[(int(w[i, j // 64]) >> (j % 64)) & 1 for j in range(7 * 64)] for i in range(5)], dtype=np.uint8)
    out['unpack_equals_integer_bits'] = bool((bits == ref).all())
    out['pack_inverse'] = bool((DV.to_host(pack(unpack(Wd))) == w).all())
    # xor_apply vs python integers, several column-block sizes
    m, p, W = 9, 6, 5
    A = rng.randint(0, 2, size=(m, p))
    Pm = rng.randint(0, 2**63, size=(p, W), dtype=np.int64).astype(np.uint64)
    R0 = rng.randint(0, 2**63, size=(m, W), dtype=np.int64).astype(np.uint64)
    exp = R0.copy()
    for i in range(m):
        for j in range(p):
            if A[i, j]:
                exp[i] ^= Pm[j]
    ok = True
    for cbv in (64, 128, 1 << 16):
        R = DV.to_dev(R0.copy())
        xor_apply(R, DV.to_dev(A.astype(F32)), DV.to_dev(Pm), 0, cbv)
        ok &= bool((DV.to_host(R) == exp).all())
    out['xor_apply_equals_integers'] = ok
    # exactness with a large inner dimension of ones: counts up to the inner dimension must stay exact
    inner = 200000
    A1 = np.ones((2, inner), dtype=F32)
    B1 = np.ones((inner, 64), dtype=F32)
    s = DV.to_host(DV.to_dev(A1) @ DV.to_dev(B1))
    out['float32_count_exact_200000'] = bool((s == inner).all())
    # rref_block vs host Gauss-Jordan
    ok = True
    for trial in range(6):
        mm, cc = int(rng.randint(1, 40)), int(rng.randint(1, 300))
        M = (rng.rand(mm, cc) < rng.choice([0.05, 0.5])).astype(np.uint8)
        if trial == 5 and mm > 3:
            M[2] = M[0] ^ M[1]
        Wc = -(-cc // 64)
        padm = np.zeros((mm, Wc * 64), dtype=np.uint8)
        padm[:, :cc] = M
        Rb = pack(DV.to_dev(padm))
        Acc, piv = rref_block(Rb, 64)
        got = DV.to_host(unpack(Acc))[:, :cc]
        h1, rk, pv = hash_from_dense(cc, 0, M)
        ok &= (rk == len(piv)) and (sorted(pv.tolist()) == sorted(piv.tolist()))
        # reduced: identity on pivots
        ok &= bool((got[:, piv] == np.eye(len(piv), dtype=np.uint8)).all())
    out['rref_block_equals_host_gauss'] = bool(ok)
    return out


def selftest_systems(device_name, cases, tiny=False):
    """run the embedded small systems; every case carries reference values from the existing engine."""
    res = []
    for c in cases:
        kw = dict(device=device_name)
        if tiny:
            kw.update(batch=3, tile_rows=8, col_block=64, compact_frac=0.05)
        prof = local_profile_bin(c['equations'], c['variables'], c['points'], c['max_degree'], kw)
        row = dict(name=c['name'], status=prof['status'], lfd=prof['exact_lfd'], cert=prof['certificate_degree'],
                   dims=prof['quotient_ranks'], falls=prof['strict_falls'],
                   hashes=[r['rref_sha256'] for r in prof['per_degree']])
        row['ok'] = (row['lfd'] == c['lfd'] and row['cert'] == c['cert'] and row['dims'] == c['dims']
                     and row['falls'] == c['falls'] and (c.get('hashes') is None or row['hashes'] == c['hashes']))
        res.append(row)
    return res


SELFTEST_CASES = []     # replaced by the generated literal at the end of the file
# ---- generated selftest cases (scripts/../validation/gen_selftest_cases.py) ----

SELFTEST_CASES = [{"name":"v12","system_hash":"80151f01b44c6f4e385b3dec319f63af441855ad1fa91d55197b0a63a2f613b5","variables":10,"equations":[[0,1,32,40,48,66,68,72,80,130,136,257,258,260,264,513,514,528],[1,2,4,8,32,33,34,36,40,64,65,68,72,128,129,130,132,136,256,257,258,260,264,272,520,528],[1,2,4,16,32,34,36,48,64,65,66,72,128,129,132,144,258,264,272,512,513,516,520],[1,2,16,32,33,36,40,48,64,72,80,129,144,257,258,512,513,514,516,528],[4,16,33,36,40,48,72,128,129,136,144,257,258,260,272,512,513,516,520],[16,34,36,40,65,68,129,130,132,144,257,512,516]],"points":[3,96,110,182,201,238,294,308,349,451,455,649,661,692,709,766,861,938,954,983],"max_degree":5,"lfd":4,"cert":4,"dims":[0,0,15,123,366],"falls":[2,3,4],"hashes":["18f98f4d96352a3931e6bc3e6e341b7176ab3ddccfc1b0012f123e2fc3083db8","7ac980f3dd425e16270979b4185f0859ea547406c9d29876987c9a90cc7498e6","5a9022f2d98d81b3e28969c1f244594aa6039640a4597892aa9a0af812533215","4e9ce0586d20a552f41e11532fc68f78ab49db55875ff931e3a76842dca84f88","1518969a0b5f10c0c513fc7f23127e1ee1c56bb5ee10daab3fe7cb217de76bea"]},{"name":"v24","system_hash":"2a5126dba54c057dda295e4bde1889e0735095c79cd3ff184c7228b0f352baec","variables":12,"equations":[[0,4,8,16,66,72,80,129,144,160,256,260,288,512,513,520,528,1024,1025,1026,1032,2050,2052,2080],[1,2,8,64,65,66,68,72,80,96,128,129,130,132,160,257,258,264,288,512,513,516,520,528,1025,1032,1056,2049,2050,2052,2064,2080],[1,2,64,66,80,96,128,129,130,136,160,260,288,514,1025,1056,2049,2050,2052,2064],[1,4,8,16,64,72,130,132,136,144,256,258,260,272,288,512,513,514,520,528,544,1024,1026,1028,1032,1040,2052,2056],[2,4,8,16,32,65,66,128,129,160,256,264,272,288,512,516,520,544,1024,1028,1040,2048,2050,2052,2056],[2,4,16,66,72,80,96,128,129,130,136,144,256,260,288,513,514,520,544,1024,1025,1026,2049,2052,2056,2080],[2,4,66,72,128,129,132,136,256,258,264,513,514,516,528,1032,1056,2064,2080],[2,8,32,80,96,128,144,260,272,512,544,1025,1026,1028,1056,2048,2049,2056,2064,2080],[4,8,66,68,72,80,129,130,132,136,144,160,256,257,258,272,288,512,513,514,520,1025,1026,1028,1040,2050,2052,2080],[4,32,65,68,72,80,136,144,256,257,260,264,288,513,514,516,520,528,1025,1026,1032,1040,2048,2052,2080]],"points":[1070,2040,2554,2960,3615,3751],"max_degree":4,"lfd":3,"cert":3,"dims":[0,0,21,293],"falls":[2,3],"hashes":["dfaca609fcf755d09dd92159cbb7a15988260fb4b787a898c07290450ee738cc","2de571b44229474d92842ce29d25fff79c39221cc55ed0d61f1cdaa7937d94c1","7a43a7d5b55ee7e6b9157e9beeba5202c0c4d435f963591cbbfe3fb37b0a0e9d","8f54a022919a5e0c2223758f380663cd3b71b489c64417f794cde69699f1c4de"]},{"name":"v13","system_hash":"c5f6226cb67217e58994a383d6ac9fe5b179cca9a1a020e92d95261365ea37d7","variables":10,"equations":[[0,1,8,32,33,34,40,65,72,136,144,256,257,258,260,272,516,520,528],[1,2,4,32,33,34,40,48,64,65,80,128,144,257,513,514,516,528],[1,2,8,16,32,34,48,64,65,136,144,256,260,512,513,516],[1,2,8,32,33,40,64,66,68,130,136,256,257,260,272,520,528],[1,2,16,32,34,40,48,64,65,80,136,144,257,260,272,512,513,514,516,520],[1,2,32,33,36,64,68,72,80,129,130,132,136,258,260,264,514],[1,4,8,16,32,34,36,48,65,66,68,128,129,130,132,136,256,260,264,272,512,513,520],[2,4,16,33,34,36,40,64,65,66,72,80,128,129,132,257,258,264,512,514,528],[2,4,33,34,36,64,65,66,72,128,129,132,258,528],[2,8,16,34,36,48,64,65,129,144,256,264,272,512,513,516,520],[8,34,40,48,65,66,68,72,80,130,132,256,257,258,272,513,514,520],[16,33,36,40,48,68,72,129,130,136,144,257,258,260,512,513,516]],"points":[],"max_degree":4,"lfd":3,"cert":3,"dims":[0,0,21,176],"falls":[2,3],"hashes":["18f98f4d96352a3931e6bc3e6e341b7176ab3ddccfc1b0012f123e2fc3083db8","7ac980f3dd425e16270979b4185f0859ea547406c9d29876987c9a90cc7498e6","19328638b6cc261a96a33aea17f31b924e21e27606348f85c41c68df359e5af2","f5462d9fd3381b8bbc38c949d5b0b6d3f5cd0a4c613ec9160950e79820d600bb"]},{"name":"v19","system_hash":"09014b53f74140a9461ebf3b431ef5737cb736da920cd53db9de02bef9bbbded","variables":10,"equations":[[0,1,16,32,33,34,40,48,65,68,72,130,257,258,272,512,513,520],[1,2,4,8,16,32,40,64,72,80,128,132,144,256,257,258,264,512,514,516],[1,2,4,16,32,33,48,64,68,72,80,128,130,132,136,144,258,260,264,512,513,514,516,528],[1,2,8,32,33,36,48,64,66,68,80,129,130,132,136,256,260,272,513,514,520],[1,32,33,36,66,72,80,129,132,136,144,258,260,264,272,514,516,520],[2,4,36,64,66,72,128,129,144,258,264,272,516,520,528],[2,16,33,48,64,68,72,80,130,132,136,144,258,260,264,512,513,514,516,528],[4,8,16,33,34,48,65,68,80,128,130,132,144,256,264,272,512,513,514,516,520],[8,33,34,36,65,66,80,129,132,144,256,264,514,516,528],[34,36,65,66,129,136,144,260,264,516]],"points":[],"max_degree":3,"lfd":2,"cert":2,"dims":[0,0,56],"falls":[2],"hashes":["18f98f4d96352a3931e6bc3e6e341b7176ab3ddccfc1b0012f123e2fc3083db8","7ac980f3dd425e16270979b4185f0859ea547406c9d29876987c9a90cc7498e6","c906e75dd27c25c2ce102b6f803e95769656d304d6cba5f0a0f567340ea5cebc"]},{"name":"v26","system_hash":"0898244a2d55639c7f81e3e21dc81abd5c0c6595a04b57f08d8c45a44ef4ae17","variables":12,"equations":[[0,65,66,68,80,96,129,130,132,257,258,260,264,272,516,1025,1028,1040,2049],[65,66,68,80,129,160,257,272,520,1025,1028,2050],[65,66,80,129,130,132,136,160,258,264,272,514,516,520,528,1025,1028,1032,1056,2050,2064],[65,72,136,260,272,288,513,514,528,1028,1032,2052,2080],[65,132,144,258,288,528,1026,1032,2052],[66,80,129,130,132,136,144,160,258,260,264,272,288,514,516,544,1025,1026,1028,1040,1056,2050,2052,2056,2064],[66,129,130,132,258,272,544,1028,1040,2056]],"points":[96,234,358,456,519,717,802,843,946,1271,1321,1371,1628,1749,1817,2029,2049,2188,2286,2437,2555,2644,2691,2911,2979,3214,3320,3539,3635,3709,3815,3961],"max_degree":4,"lfd":3,"cert":3,"dims":[0,0,7,267],"falls":[3],"hashes":["dfaca609fcf755d09dd92159cbb7a15988260fb4b787a898c07290450ee738cc","2de571b44229474d92842ce29d25fff79c39221cc55ed0d61f1cdaa7937d94c1","6fd09f5b3c0db5991ea46d41748d93ecfea7e507d0db06b3aa2aa990b9bafa23","b92e64b84805b13576ae70b9a64820d4c4db103e4859e94a38cb390ec020fa3a"]},{"name":"v37","system_hash":"59195bcafc58457f9a10a39bfc15c062131a86f0d031f3c5dbac6210839d62a0","variables":14,"equations":[[0,1,2,16,32,64,128,136,256,258,264,520,528,544,576,1025,1026,1028,1032,1088,2048,2052,2064,4096,4100,4128,4160,8192,8196,8200,8224,8256],[1,2,8,64,128,129,136,144,256,260,264,272,320,514,516,544,576,1024,1025,1026,1032,2049,2050,4100,8192,8194,8196],[1,2,8,64,128,132,144,192,256,260,264,288,513,514,528,576,1024,1026,1040,1088,2049,2052,2056,2064,2080,4098,4112,8192,8193,8196,8200,8256],[1,2,16,32,64,128,130,144,256,257,260,264,272,320,514,520,544,1026,1028,1032,1056,1088,2048,2049,2050,2064,2112,4096,4100,4104,4128,4160,8192,8194,8200,8208,8224,8256],[1,4,16,32,64,128,130,132,136,144,257,260,272,320,512,513,514,528,544,576,1025,1056,2048,2049,2050,2052,4096,4100,4104,4160,8192,8194,8196,8224,8256],[1,4,16,32,128,130,136,144,160,257,264,320,512,516,520,528,576,1025,1026,1028,1056,2048,2049,2052,4096,4097,4104,8194,8196,8256],[1,32,128,130,136,144,160,257,264,288,516,520,544,1025,1026,1028,1032,1040,2049,2056,2064,4096,4097,4098,4100],[1,64,128,129,130,132,144,160,257,258,260,288,513,514,520,576,1028,1032,1056,1088,2049,2112,4097,4098,4104,4128,8192,8196,8200,8208],[2,4,16,32,64,132,136,192,256,258,260,272,512,513,514,520,528,544,576,1025,1028,1032,1040,1088,2048,2050,2052,2056,2080,4096,4100,4112,8192,8193,8196,8200],[2,4,16,129,136,144,256,320,512,516,528,1025,1056,1088,2048,2049,2052,2064,2080,4104,4112,8194,8200]],"points":[886,1515,1694,3853,5192,5727,9256,11127,12159,12204,12400,13707,14432,15110,15318,16350],"max_degree":5,"lfd":4,"cert":4,"dims":[0,0,23,309,1455],"falls":[2,3,4],"hashes":["a4a9138490254945e7017befb660f805e8d95fab87e21eb0d4f1584aa2f0f4cb","66d6e11ff5a1b1f47e0f37993df34ea970cb255a328b1a4c99aa2df019e210fa","fc05fca8037975ca98cec579c064f35e21d2c05da611b5e334a82de0fb58e6ad","4e70bef1c70c74bd5818d559ff03268e524d1c8047664c30b64f59ded0e4a184","f26cbb30d20e006a9f877107fa1ff4a6e8bb4380d81bd9887ad83df2a1581c83"]},{"name":"sf12_6_a0","system_hash":"25f03d5a984438083e7eff6490cb527a5680107b847af68fba387fba55a93778","variables":12,"equations":[[0,68,144,257,544,1026,2056,2080],[65,80,136,264,288,514,516,1025,2052],[65,160,264,516,1040,2050],[66,129,144,272,520,544,1026,1028,2056],[72,96,130,132,258,513,528,1032,2049],[96,130,260,528,1032,1056,2049,2064]],"points":[68,145,239,257,374,400,479,552,626,666,722,797,880,911,974,1030,1090,1163,1253,1338,1372,1459,1533,1589,1639,1674,1761,1813,1868,1969,1991,2080,2139,2238,2286,2363,2387,2477,2521,2568,2667,2740,2793,2871,2918,2979,3011,3085,3166,3209,3286,3370,3416,3461,3564,3644,3711,3732,3812,3896,3927,4002,4089],"max_degree":5,"lfd":4,"cert":4,"dims":[0,0,6,234,731],"falls":[3,4],"hashes":["dfaca609fcf755d09dd92159cbb7a15988260fb4b787a898c07290450ee738cc","2de571b44229474d92842ce29d25fff79c39221cc55ed0d61f1cdaa7937d94c1","da40a4b17f527a7f3e5520da10322f7c338b3f985c96afef953c33ab162759bd","01e7d852076db06a09b93a74ee3623e6d9d51c8134d54dd909ec341f0853f74d","5dd7e0953b4776a5ac9a84e50b1d17ce70bfaf6300dc9413d87f340404e6d73a"]}]
