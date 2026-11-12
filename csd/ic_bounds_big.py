"""Closure engine for the degree-D chain closure V_D (two summands, depth k) at large size, with a CPU
(numpy) and a GPU (cupy) backend that share one code path.  Built for k = 14, D = 6, p = 65537
(499,178 columns).  Single self-contained file: Python 3.10+, numpy; cupy only for ``--device gpu``;
threadpoolctl only if importable.  It imports nothing from the repository.

MATHEMATICS (identical to csd/ic_bounds_d5_verify.py, the reference engine).  Modulo the chains the ring
is F_p[(Z/2^k)^2]; a squarefree chain monomial is an exponent pair (a1, a2) of weight popcount(a1) +
popcount(a2); u_{i,j} adds 2^j to a_i mod 2^k.  V_D is the least subspace of span(M_<=D) containing
NF(g) (and NF(m g), weight(m) <= D - deg g, with ``gens='all'``) and closed under f -> u_{i,j} f for every
f supported on weight < D.  Columns are ordered by weight descending, ties by key a1 | a2 << k.  With the
summand swap split (g symmetric) V_D = V+ (+) V-; parts are kept in orbit coordinates, only the k block-1
shifts are applied and every shifted row is projected to both parts.  Each part is a reduced row echelon
basis (RREF) for that column order: row i = e_{piv[i]} + B[i] on the free columns.  The RREF of a
subspace is unique, so the basis hash ``Part.sha()`` (same definition as the reference) is independent of
the order in which the closure is carried out; that is how this engine is validated.

DESIGN
  * ``xp`` is numpy or cupy.  Index bookkeeping (pivots, column maps, layouts, processed flags) stays on
    the host in numpy; only B (float32 tiles), batch blocks and residual blocks live on the device.
    Device specific code is limited to ``Device`` (transfer, free pool, memory info).  ``fake-gpu`` runs
    the device code path with a numpy stand-in that REFUSES to mix host and device arrays, to catch API
    slips without a GPU.
  * Exactness.  p < 2^24 so float32 storage is exact.  Every product is a float64 product of entries in
    [0, p); the inner dimension accumulated between two reductions is at most ``max_inner(p)`` =
    floor((2^50 - 4p) / (p-1)^2), so every partial sum is an integer of magnitude <= 2^50 < 2^53 and is
    exact whatever the summation order, FMA use or blocking.  The reduction kernel ``modp_`` is exact for
    |X| <= 2^50 (proof at the function) and is tested against Python integers in ``selftest_arith``.
  * Memory.  B is stored in row tiles; a rank-t update touches only the rows that meet the new pivots, in
    place, tile by tile and column block by column block.  New pivot columns are deleted lazily (kept
    exactly zero) and compacted tile by tile.  No operation copies a part's whole matrix.  Working
    buffers obey ``work_bytes`` (see ``Engine.batch_rows``, ``Engine.col_block``).
  * Sparse-aware reduction.  All rows of one batch (same source part, same shift, same target part) share
    their support T; R = Y[:, free] - Y[:, piv] @ B uses only the pivot columns of T on which the batch is
    nonzero, and only the columns of B at or after the smallest such pivot (B[r] vanishes before piv[r]).
    ``madds`` counts what is done, ``madds_dense`` what the dense method (n x rank x free) would do.
  * Closure invariant (as in the reference).  A low row (pivot weight < D) is processed once, at its state
    when its batch starts.  A low row is modified only by adding multiples of rows whose pivot is later in
    the column order; those rows are low too (a low row's support has weight < D), and are processed in
    their turn.  Hence the final low rows lie in the span of the processed states, whose shifts were all
    reduced into the basis, and the final space is closed.  The certificate (sweeps) does not rely on this.
  * Certificate.  Randomized sweep (t random combinations of the low rows of each part) then the
    deterministic full sweep: every shift of every final low row, both projections, must reduce to zero
    against the final, unchanged basis.  A nonzero residual is added (repair), the closure continues and
    a new sweep starts; ``closed`` is true only after one complete clean sweep over an unchanged basis.
  * Checkpoint/resume: single slot by default; temporary directory, fsync, manifest last, atomic rename,
    then the old slot is deleted.  SIGTERM: finish the batch, checkpoint, exit code 75.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import shutil
import signal
import sys
import time
import weakref

import numpy as np

try:  # optional: BLAS threads on the CPU path
    from threadpoolctl import threadpool_limits as _tpl
except Exception:  # pragma: no cover
    _tpl = None

F64 = np.float64
F32 = np.float32
MAXMAG = 1 << 50
ENGINE_VERSION = 'ic_bounds_big/2'
EXIT_SIGTERM = 75


def file_sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 22), b''):
            h.update(blk)
    return h.hexdigest()


ENGINE_SHA256 = file_sha256(__file__)


def cgroup_cpu_limit():
    """CPU limit of the cgroup (cpu.max of cgroup v2, else cfs quota of v1) rounded up, or None."""
    try:
        with open('/sys/fs/cgroup/cpu.max') as f:
            q, per = f.read().split()[:2]
        if q != 'max':
            return max(1, -(-int(q) // int(per)))
    except Exception:
        pass
    try:
        q = int(open('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read())
        per = int(open('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read())
        if q > 0:
            return max(1, -(-q // per))
    except Exception:
        pass
    return None


def default_threads() -> int:
    """CPUs this process may use: the minimum of the affinity count and the cgroup CPU quota (the affinity
    count alone can be the host's cores in a quota-limited pod; os.cpu_count() always is)."""
    try:
        n = max(1, len(os.sched_getaffinity(0)))
    except Exception:  # pragma: no cover
        n = 1
    lim = cgroup_cpu_limit()
    return max(1, min(n, lim)) if lim else n


def maxrss_mb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


# =============================================================================== device layer

class _Live:
    """bytes held by live fake-device arrays (owning arrays only): the stand-in for device memory use."""
    cur = 0
    peak = 0
    limit = None        # fake-gpu budget: exceeding it raises (catches accounting errors)
    reg: dict = {}      # id(root buffer) -> bytes, for every live root buffer behind a device array

    @classmethod
    def track(cls, a):
        """account the root buffer behind ``a`` once (numpy results of fancy indexing are views of a temporary)."""
        root = a
        while isinstance(root.base, np.ndarray):
            root = root.base
        k = id(root)
        if k in cls.reg or root.nbytes == 0:
            return
        cls.reg[k] = root.nbytes
        weakref.finalize(root, cls._free, k)
        cls.add(root.nbytes)

    @classmethod
    def _free(cls, k):
        cls.cur -= cls.reg.pop(k)

    @classmethod
    def add(cls, n):
        cls.cur += n
        cls.peak = max(cls.peak, cls.cur)
        if cls.limit is not None and cls.cur > cls.limit:
            raise MemoryError(f'fake-gpu: live device bytes {cls.cur} exceed the budget {cls.limit}')


class _Dev(np.ndarray):
    """numpy stand-in for a device array (``--device fake-gpu``): behaves like numpy but raises when a
    plain host ndarray is mixed into an operation, an index or an assignment, as cupy would."""

    @staticmethod
    def _chk(x, what):
        if isinstance(x, np.ndarray) and not isinstance(x, _Dev):
            raise TypeError(f'fake-gpu: host ndarray used in {what}')

    @classmethod
    def _wrap(cls, r):
        if isinstance(r, _Dev):
            return r
        if isinstance(r, np.ndarray):
            v = r.view(_Dev)
            _Live.track(v)
            return v
        if isinstance(r, np.generic):          # cupy reductions return 0-d device arrays
            return np.asarray(r).view(_Dev)
        if isinstance(r, tuple):
            return tuple(cls._wrap(x) for x in r)
        if isinstance(r, list):
            return [cls._wrap(x) for x in r]
        return r

    def __array_finalize__(self, obj):
        _Live.track(self)       # every new device array (owner or view) registers its root buffer once

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kw):
        for x in inputs:
            self._chk(x, ufunc.__name__)
        if out is not None:
            for x in out:
                self._chk(x, ufunc.__name__)
            kw['out'] = tuple(o.view(np.ndarray) if isinstance(o, _Dev) else o for o in out)
        ins = tuple(x.view(np.ndarray) if isinstance(x, _Dev) else x for x in inputs)
        return self._wrap(getattr(ufunc, method)(*ins, **kw))

    def __array_function__(self, func, types, args, kwargs):
        def flat(o):
            if isinstance(o, (list, tuple)):
                for x in o:
                    yield from flat(x)
            elif isinstance(o, dict):
                for x in o.values():
                    yield from flat(x)
            else:
                yield o
        for x in flat(args):
            self._chk(x, func.__name__)
        for x in flat(kwargs):
            self._chk(x, func.__name__)
        return self._wrap(super().__array_function__(func, types, args, kwargs))

    def __getitem__(self, idx):
        for x in (idx if isinstance(idx, tuple) else (idx,)):
            self._chk(x, 'index')
        return self._wrap(np.ndarray.__getitem__(self, idx))

    def __setitem__(self, idx, val):
        for x in (idx if isinstance(idx, tuple) else (idx,)):
            self._chk(x, 'index')
        self._chk(val, 'assignment')
        np.ndarray.__setitem__(self, idx, val)


class _FakeXP:
    """array-module facade for fake-gpu: results are _Dev, arguments must not be plain host arrays."""
    float64 = np.float64
    float32 = np.float32
    int64 = np.int64
    int8 = np.int8

    def __getattr__(self, name):
        f = getattr(np, name)
        if not callable(f):
            return f

        def call(*a, **k):
            for x in list(a) + list(k.values()):
                for y in (x if isinstance(x, (list, tuple)) else (x,)):
                    _Dev._chk(y, name)
            k = {key: (v.view(np.ndarray) if isinstance(v, _Dev) else v) for key, v in k.items()}
            return _Dev._wrap(f(*[(x.view(np.ndarray) if isinstance(x, _Dev) else x) for x in a], **k))
        return call


class Device:
    """The few device specific operations."""

    def __init__(self, name: str):
        self.name = name
        self.is_dev = name != 'cpu'
        if name == 'cpu':
            self.xp = np
        elif name == 'fake-gpu':
            self.xp = _FakeXP()
        elif name == 'gpu':
            import cupy as cp  # noqa: only when asked for
            self.cp = cp
            self.xp = cp
        else:
            raise ValueError(f'unknown device {name}')

    def to_dev(self, h):
        h = np.asarray(h)
        if self.name == 'cpu':
            return h
        if self.name == 'fake-gpu':
            return _Dev._wrap(np.array(h))
        return self.cp.asarray(h)

    def to_host(self, a):
        if self.name == 'cpu':
            return np.asarray(a)
        if self.name == 'fake-gpu':
            assert isinstance(a, _Dev), 'to_host of a host array'
            return np.array(a.view(np.ndarray))
        return self.cp.asnumpy(a)

    def sync(self):
        if self.name == 'gpu':
            self.cp.cuda.Stream.null.synchronize()

    def free_pool(self):
        if self.name == 'gpu':
            self.cp.get_default_memory_pool().free_all_blocks()
            self.cp.get_default_pinned_memory_pool().free_all_blocks()

    def mem_info(self):
        """(used bytes, total bytes) of the device; (0, 0) when there is none."""
        if self.name == 'gpu':
            free, total = self.cp.cuda.runtime.memGetInfo()
            return total - free, total
        if self.name == 'fake-gpu':
            return _Live.cur, 0
        return 0, 0

    def pool_bytes(self):
        if self.name == 'gpu':
            return int(self.cp.get_default_memory_pool().total_bytes())
        if self.name == 'fake-gpu':
            return _Live.peak          # peak live bytes of the stand-in
        return 0


DV = Device('cpu')
xp = np


def use_device(name: str) -> Device:
    global DV, xp
    DV = Device(name)
    xp = DV.xp
    return DV


# =============================================================================== exact arithmetic

def max_inner(p: int) -> int:
    """largest n with n (p-1)^2 + 4p <= 2^50: the inner dimension that may be accumulated between two
    reductions when all entries are in [0, p) and the accumulator starts in [0, p)."""
    return (MAXMAG - 4 * p) // ((p - 1) ** 2)


GUARD_HITS = [0]     # number of times the correction step below was needed (the proof says never)


def modp_(X, p: int, chunk: int = 1 << 22):
    """In place X <- X mod p (result in [0, p)) for a float64 array of integers with |X| <= 2^50, p < 2^24.

    Computes q = floor(fl(fl(X + 0.5) * fl(1/p))), r = X - q p, then a correction step.

    Proof.  X + 0.5 is exact (|X| < 2^52).  Put u = (X + 0.5)/p.  Since p u = X + 1/2 is not a multiple of
    p, u is at distance at least 1/(2p) from every integer.  fl(1/p) = (1/p)(1 + e1) and the product rounds
    with (1 + e2), |e1|, |e2| <= 2^-53, so the computed v satisfies |v - u| <= |u| 2.3 * 2^-53
    <= (2^50.01 / p) 2^-51.8 < 1/(2p): hence floor(v) = floor(u) = floor(X / p) (adding 1/2 to an
    integer X cannot cross a multiple of p).  q p is an integer of magnitude <= |X| + p < 2^53, so q p and
    X - q p are exactly representable and exact (also if a fused multiply-add is used), and
    r = X - q p lies in [0, p).  So the correction step (r < 0: r += p; r >= p: r -= p) never fires; it
    is kept as a guard that makes the result correct even if the rounding analysis were off by one
    quotient (it then repairs any q that is off by at most 1) and counts its use in GUARD_HITS.
    The result is tested against Python integers in ``selftest_arith`` (random data and the extremes
    0, +-1, +-p, +-2^50, multiples of p +- 1)."""
    inv = 1.0 / p
    if X.ndim == 1:
        X = X.reshape(1, -1)
    rc = max(1, chunk // max(1, X.shape[1]))
    for i in range(0, X.shape[0], rc):
        x = X[i:i + rc]
        if x.size == 0:
            continue
        q = x + 0.5
        q *= inv
        xp.floor(q, out=q)
        q *= p
        x -= q
        del q
        if float(x.min()) < 0 or float(x.max()) >= p:
            GUARD_HITS[0] += 1
            neg = (x < 0).astype(F64)
            neg *= p
            x += neg
            big = (x >= p).astype(F64)
            big *= p
            x -= big
    return X


def matmul_sub(acc, A, Bm, p: int, mi: int):
    """acc <- (acc - A @ Bm) mod p exactly.  acc, A, Bm have entries in [0, p); the inner dimension is
    split into chunks of <= mi with a reduction after each chunk."""
    n = A.shape[1]
    for a in range(0, n, mi):
        acc -= A[:, a:a + mi] @ Bm[a:a + mi]
        modp_(acc, p)
    return acc


# =============================================================================== instance (copied)

def s3_coeffs(p: int, a: int, b: int, X: int) -> dict:
    """S_3(x1, x2, X) for y^2 = x^3 + a x + b (copied from csd/ic_bounds_d5_verify.py)."""
    c = {}

    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p

    add((2, 0), X * X); add((1, 1), -2 * X * X); add((0, 2), X * X)
    add((2, 1), -2 * X); add((1, 2), -2 * X); add((1, 0), -2 * X * a); add((0, 1), -2 * X * a)
    add((0, 0), -4 * X * b)
    add((2, 2), 1); add((1, 1), -2 * a); add((0, 0), a * a)
    add((1, 0), -4 * b); add((0, 1), -4 * b)
    return {e: v for e, v in c.items() if v}


def chain_generator(p: int, a: int, b: int, X: int, k: int) -> dict:
    N = 1 << k
    g = {}
    for (e1, e2), v in s3_coeffs(p, a, b, X).items():
        key = (e1 % N, e2 % N)
        g[key] = (g.get(key, 0) + v) % p
    return {e: v for e, v in g.items() if v}


def subgroup(p: int, k: int):
    N = 1 << k
    assert (p - 1) % N == 0
    for t in range(2, p):
        h = pow(t, (p - 1) // N, p)
        if pow(h, N // 2, p) != 1:
            return h, [pow(h, e, p) for e in range(N)]
    raise ValueError


def solutions(p: int, a: int, b: int, X: int, k: int) -> list:
    """All (E1, E2) with S_3(h^E1, h^E2, X) = 0 (vectorized over E2; copied)."""
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


def count_columns(k: int, D: int) -> int:
    from math import comb
    return sum(comb(2 * k, w) for w in range(D + 1))


# =============================================================================== layout (host only)

class Layout:
    """Column order, shifts and swap split of the reference engine (built per weight pair, so that no
    (#small)^2 meshgrid is formed).  Pure numpy on the host."""

    def __init__(self, k: int, D: int, split: bool):
        self.k, self.D, self.split = k, D, split
        N = 1 << k
        pc = np.array([bin(x).count('1') for x in range(N)], dtype=np.int64)
        byw = [np.flatnonzero(pc == w) for w in range(D + 1)]
        A1s, A2s, Ws = [], [], []
        for w1 in range(D + 1):
            for w2 in range(D + 1 - w1):
                x1, x2 = np.meshgrid(byw[w1], byw[w2], indexing='ij')
                A1s.append(x1.ravel()); A2s.append(x2.ravel())
                Ws.append(np.full(x1.size, w1 + w2, dtype=np.int64))
        A1, A2, W = np.concatenate(A1s), np.concatenate(A2s), np.concatenate(Ws)
        key = A1 | (A2 << k)
        order = np.lexsort((key, -W))        # weight descending, then key ascending
        self.a1, self.a2, self.w, self.key = A1[order], A2[order], W[order], key[order]
        self.c = len(self.key)
        assert self.c == count_columns(k, D)
        self._sorted = np.argsort(self.key)
        self._skeys = self.key[self._sorted]
        self.low = np.flatnonzero(self.w < D)
        self.lowc = len(self.low)
        self.low0 = self.c - self.lowc
        assert np.all(self.low == np.arange(self.low0, self.c))
        blocks = (0,) if split else (0, 1)
        self.shift_list = [(i, j) for i in blocks for j in range(k)]
        self.shift = {}
        for (i, j) in self.shift_list:
            b1, b2 = self.a1[self.low], self.a2[self.low]
            if i == 0:
                b1 = (b1 + (1 << j)) % N
            else:
                b2 = (b2 + (1 << j)) % N
            self.shift[(i, j)] = self.index(b1, b2)
        self.parts = []
        if split:
            sig = self.index(self.a2, self.a1)
            diag = self.a1 == self.a2
            for reps, sign in ((np.flatnonzero(self.a1 <= self.a2), 1), (np.flatnonzero(self.a1 < self.a2), -1)):
                pos = -np.ones(self.c, dtype=np.int64)
                pos[reps] = np.arange(len(reps))
                lift_idx = np.where(pos >= 0, pos, pos[sig])
                lift_sign = np.where(pos >= 0, 1, sign).astype(np.int64)
                if sign == -1:
                    lift_sign[diag] = 0
                    lift_idx[diag] = 0
                # projection Y[m] = H[reps[m]] + sign H[sig reps[m]]: the full column h contributes
                # proj_coef[h] * H[h] to part column lift_idx[h] (2 on the symmetric diagonal)
                proj_coef = lift_sign.copy()
                if sign == 1:
                    proj_coef[diag] = 2
                self.parts.append(dict(reps=reps, rep_sig=sig[reps], sign=sign, w=self.w[reps], c=len(reps),
                                       lift_idx=lift_idx, lift_sign=lift_sign, proj_coef=proj_coef))
        else:
            ar = np.arange(self.c)
            one = np.ones(self.c, dtype=np.int64)
            self.parts.append(dict(reps=ar, rep_sig=ar, sign=1, w=self.w, c=self.c, lift_idx=ar,
                                   lift_sign=one, proj_coef=one))
        for P in self.parts:   # low columns of a part are a suffix (weight descending)
            P['s0'] = int(np.searchsorted(-P['w'], -(self.D - 1), side='left'))
            assert np.all(P['w'][P['s0']:] < self.D) and np.all(P['w'][:P['s0']] >= self.D)
        self._maps = {}

    def index(self, b1, b2):
        kk = np.asarray(b1) | (np.asarray(b2) << self.k)
        pos = np.searchsorted(self._skeys, kk)
        assert np.all(self._skeys[np.minimum(pos, self.c - 1)] == kk)
        return self._sorted[pos]

    def shift_map(self, pi: int, ij, q: int):
        """Sparse map (source part low column s) -> (target part column m, coef, layer) of
        lift, u_ij, project.  Within one layer the m are distinct (layer 0: H[reps[m]], layer 1:
        H[sigma reps[m]])."""
        key = (pi, ij, q)
        if key in self._maps:
            return self._maps[key]
        P, Q = self.parts[pi], self.parts[q]
        s = np.flatnonzero(P['w'] < self.D)
        cols = [P['reps'][s]]
        coefs = [np.ones(len(s), dtype=np.int64)]
        srcs = [s]
        if self.split:
            rs = P['rep_sig'][s]
            ok = rs != P['reps'][s]
            cols.append(rs[ok]); coefs.append(np.full(int(ok.sum()), P['sign'], dtype=np.int64)); srcs.append(s[ok])
        col, cf, src = np.concatenate(cols), np.concatenate(coefs), np.concatenate(srcs)
        h = self.shift[ij][col - self.low0]
        m = Q['lift_idx'][h]
        coef = cf * Q['proj_coef'][h]
        layer = np.where(Q['reps'][m] == h, 0, 1).astype(np.int8)
        keep = coef != 0
        out = (src[keep], m[keep], coef[keep].astype(F64), layer[keep])
        for ly in (0, 1):
            mm = out[1][out[3] == ly]
            assert len(np.unique(mm)) == len(mm), 'shift map layer not injective'
        self._maps[key] = out
        return out

    def lift_full(self, pi: int, M: np.ndarray) -> np.ndarray:
        """host rows in part coordinates -> full coordinates (as the reference ``lift``)."""
        if not self.split:
            return M
        P = self.parts[pi]
        return M[:, P['lift_idx']] * P['lift_sign']


# =============================================================================== echelon of a block

def _rref_rows(S, p: int, mi: int):
    """Row-by-row RREF of a small block S (m x f, entries in [0,p)); returns (L, pivots)."""
    m, f = S.shape
    L = xp.zeros((m, f), dtype=F64)
    pl = []
    cnt = 0
    for i in range(m):
        v = S[i:i + 1]
        if cnt:
            matmul_sub(v, v[:, DV.to_dev(np.array(pl, dtype=np.int64))], L[:cnt], p, mi)
        lead = int(xp.argmax((v[0] != 0).astype(xp.int8)))      # first nonzero (0 if the row is zero)
        val = int(v[0, lead])
        if val == 0:
            continue
        v *= float(pow(val, -1, p))      # entries < p, factor < p: product < 2^48, exact
        modp_(v, p)
        if cnt:
            L[:cnt] -= xp.outer(L[:cnt, lead], v[0])
            modp_(L[:cnt], p)
        L[cnt] = v[0]
        pl.append(lead)
        cnt += 1
    return L[:cnt], pl


def rref_block(R, p: int, mi: int, sub: int = 32):
    """RREF of the rows of R (n x f, entries in [0,p)), computed in blocks of ``sub`` rows.  R is used
    as storage: on return R[:t] is the RREF (identity on the returned leading columns) and the list
    ``lead`` holds the leading column of each row.  Returns (t, lead)."""
    n = R.shape[0]
    t = 0
    lead: list = []
    for a in range(0, n, sub):
        S = xp.array(R[a:a + sub])
        if t:
            matmul_sub(S, S[:, DV.to_dev(np.array(lead, dtype=np.int64))], R[:t], p, mi)
        keep = np.flatnonzero(DV.to_host(xp.any(S != 0, axis=1)))
        if keep.size == 0:
            continue
        S = S[DV.to_dev(keep)]
        Ln, pn = _rref_rows(S, p, mi)
        tn = len(pn)
        if tn == 0:
            continue
        if t:
            cs = min(pn)
            Lt = R[:t]
            X = Lt[:, DV.to_dev(np.array(pn, dtype=np.int64))]     # copy
            matmul_sub(Lt[:, cs:], X, Ln[:, cs:], p, mi)
        R[t:t + tn] = Ln      # rows [t, t+tn) lie below a+sub: already consumed (t <= a, tn <= sub)
        t += tn
        lead += pn
    return t, lead


# =============================================================================== part

class Part:
    """RREF basis of one part: row i = e_{piv[i]} + B[i] on the alive (free) columns.

    Physical columns: pc2col[x] = part column of physical column x (increasing), alive[x] = still free;
    dead physical columns (pivot columns not yet compacted away) are exactly zero in every row.  Rows live
    in float32 tiles of ``TR`` rows; row r is tiles[r // TR][r % TR]."""

    def __init__(self, w: np.ndarray, s0: int, p: int, D: int, eng, tile_rows=None):
        assert p < (1 << 24)
        self.w, self.s0, self.p, self.D, self.c = w, s0, p, D, len(w)
        self.E = eng
        self.mi = max_inner(p)
        self.pc2col = np.arange(self.c, dtype=np.int64)
        self.alive = np.ones(self.c, dtype=bool)
        self.col2phys = np.arange(self.c, dtype=np.int64)       # -1 once a pivot
        self.col2row = -np.ones(self.c, dtype=np.int64)
        self.piv = np.zeros(0, dtype=np.int64)
        self.TR = int(tile_rows) if tile_rows else min(max(16, eng.tile_bytes // max(1, 4 * self.c)), max(16, self.c))
        self.tiles: list = []          # xp arrays (loc 'd') or numpy arrays (loc 'h', spilled)
        self.loc: list = []
        self.madds = dict(reduce=0, reduce_dense=0, update=0, rref=0)
        self.secs = dict(reduce=0.0, update=0.0, rref=0.0)
        self.peak_entries = 0           # max rank x free (the reference's definition)
        self.peak_phys_entries = 0
        self.compactions = 0
        self.zero_col = int(np.flatnonzero(w == 0)[0]) if np.any(w == 0) else -1

    # -- properties
    @property
    def rank(self):
        return len(self.piv)

    @property
    def ncols(self):
        return len(self.pc2col)

    @property
    def nfree(self):
        return self.c - self.rank

    def free_cols(self):
        return self.pc2col[self.alive]

    def tile_rows_filled(self, ti):
        return min(self.TR, self.rank - ti * self.TR)

    def low_rows(self) -> np.ndarray:
        return np.flatnonzero(self.w[self.piv] < self.D)

    def tile_bytes_now(self):
        return self.TR * self.ncols * 4

    # -- tile location (device / host spill)
    def dev_bytes(self) -> int:
        return sum(t.nbytes for t, l in zip(self.tiles, self.loc) if l == 'd')

    def host_bytes(self) -> int:
        return sum(t.nbytes for t, l in zip(self.tiles, self.loc) if l == 'h')

    def _rows(self, ti, off, sl):
        """device copy of tile[off, sl] (off: host index array), wherever the tile is."""
        t = self.tiles[ti]
        if self.loc[ti] == 'd':
            return t[DV.to_dev(off), sl]
        return DV.to_dev(t[off, sl])

    def _rows_host(self, ti, off, sl):
        """host copy of tile[off, sl]."""
        t = self.tiles[ti]
        if self.loc[ti] == 'd':
            return DV.to_host(t[DV.to_dev(off), sl])
        return np.array(t[off, sl])

    # -- row groups
    def _groups(self, rows):
        rows = np.asarray(rows, dtype=np.int64)
        t = rows // self.TR
        o = np.argsort(t, kind='stable')
        ts = t[o]
        cuts = np.flatnonzero(np.diff(ts)) + 1
        for seg in np.split(np.arange(len(rows)), cuts):
            if len(seg) == 0:
                continue
            idx = o[seg]
            yield int(ts[seg[0]]), idx, rows[idx] - ts[seg[0]] * self.TR

    # -- reduction
    def accumulate(self, rows, A, out, c0: int, sign: int):
        """out[:, c0:] <- (out[:, c0:] + sign * A @ B[rows, c0:]) mod p, exactly.  A: n x len(rows) device
        float64 in [0,p) (columns aligned with ``rows``); out: n x ncols float64 in [0,p).  Column blocks
        bound the working set; the accumulated inner dimension is capped at mi between reductions."""
        p, mi = self.p, self.mi
        n = A.shape[0]
        groups = list(self._groups(rows))
        As = [A[:, DV.to_dev(idx)] for _, idx, _ in groups]
        cw = self.E.col_block(n, self.TR)
        for a in range(c0, self.ncols, cw):
            b = min(self.ncols, a + cw)
            blk = out[:, a:b]
            cnt = 0
            for (ti, idx, off), Ag in zip(groups, As):
                for u in range(0, len(off), mi):
                    o = off[u:u + mi]
                    if cnt + len(o) > mi:
                        modp_(blk, p)
                        cnt = 0
                    G = self._rows(ti, o, slice(a, b)).astype(F64)
                    prod = Ag[:, u:u + mi] @ G
                    if sign < 0:
                        blk -= prod
                    else:
                        blk += prod
                    cnt += len(o)
            modp_(blk, p)
        del As

    def residual(self, T: np.ndarray, Yt, count_dense=True):
        """Rows of Yt (n x len(T), part columns T sorted unique, entries in [0,p), device) reduced
        against the basis, on the physical columns.  Returns the NONZERO residual rows (n' x ncols,
        float64 in [0,p)) or None if every row reduces to zero."""
        p, n = self.p, Yt.shape[0]
        self.E.activate(self)
        t0 = time.perf_counter()
        nzc = DV.to_host(xp.any(Yt != 0, axis=0))
        if not nzc.any():
            return None
        T = T[nzc]
        Yt = Yt[:, DV.to_dev(np.flatnonzero(nzc))]
        rowid = self.col2row[T]
        isp = rowid >= 0
        R = xp.zeros((n, self.ncols), dtype=F64)
        fcol = np.flatnonzero(~isp)
        if len(fcol):
            fp = self.col2phys[T[fcol]]
            assert np.all(fp >= 0)
            R[:, DV.to_dev(fp)] = Yt[:, DV.to_dev(fcol)]
        if count_dense:
            self.madds['reduce_dense'] += n * self.rank * self.nfree
        sel = rowid[isp]
        if len(sel):
            Ysel = Yt[:, DV.to_dev(np.flatnonzero(isp))]
            c0 = int(np.searchsorted(self.pc2col, int(self.piv[sel].min())))
            self.accumulate(sel, Ysel, R, c0, -1)
            self.madds['reduce'] += n * len(sel) * (self.ncols - c0)
            del Ysel
        del Yt
        nz = np.flatnonzero(DV.to_host(xp.any(R != 0, axis=1)))
        DV.sync()
        self.secs['reduce'] += time.perf_counter() - t0
        if len(nz) == 0:
            return None
        if len(nz) < n:
            R = R[DV.to_dev(nz)]
        return R

    # -- adding rows
    def add_residuals(self, R) -> int:
        """R: nonzero rows already reduced against the basis (entries in [0,p), dead columns zero).  RREF
        them, update the basis in place, append.  Returns the number of new rows."""
        self.E.activate(self)
        t0 = time.perf_counter()
        n = R.shape[0]
        t, lead = rref_block(R, self.p, self.mi)
        self.madds['rref'] += n * n * self.ncols // 2
        self.secs['rref'] += time.perf_counter() - t0
        if t == 0:
            return 0
        L = R[:t]
        lead = np.asarray(lead, dtype=np.int64)
        t0 = time.perf_counter()
        self._update(L, lead)
        L[:, DV.to_dev(lead)] = 0.0
        self._append(L)
        DV.sync()
        self.secs['update'] += time.perf_counter() - t0
        newcols = self.pc2col[lead]
        self.piv = np.concatenate([self.piv, newcols])
        self.col2row[newcols] = np.arange(self.rank - t, self.rank)
        self.alive[lead] = False
        self.col2phys[newcols] = -1
        self.peak_entries = max(self.peak_entries, self.rank * self.nfree)
        self.peak_phys_entries = max(self.peak_phys_entries, self.rank * self.ncols)
        dead = self.ncols - int(self.alive.sum())
        if dead > self.E.compact_frac * self.ncols:
            self.compact()
        return t

    def _update(self, L, lead: np.ndarray):
        """old rows <- (row - row[lead] @ L) mod p, in place, only rows meeting the new pivots and only
        columns in [min lead, last nonzero column of L).  L has the identity on the lead columns, so
        the lead columns become exactly zero."""
        if self.rank == 0:
            return
        p, mi = self.p, self.mi
        t = L.shape[0]
        cs = int(lead.min())
        colnz = np.flatnonzero(DV.to_host(xp.any(L != 0, axis=0)))
        ce = int(colnz.max()) + 1
        lead_d = DV.to_dev(lead)
        cw = self.E.col_block(max(t, self.TR), self.TR)
        for ti, tile in enumerate(self.tiles):
            nr = self.tile_rows_filled(ti)
            if nr <= 0:
                continue
            X = tile[:nr][:, lead_d].astype(F64)           # copy: nr x t
            touched = np.flatnonzero(DV.to_host(xp.any(X != 0, axis=1)))
            if len(touched) == 0:
                continue
            td = DV.to_dev(touched)
            Xt = X[td]
            del X
            for a in range(cs, ce, cw):
                b = min(ce, a + cw)
                sub = tile[td, a:b].astype(F64)
                matmul_sub(sub, Xt, L[:, a:b], p, mi)
                tile[td, a:b] = sub.astype(F32)
            self.madds['update'] += len(touched) * t * (ce - cs)

    def _append(self, L):
        t = L.shape[0]
        r0 = self.rank
        i = 0
        while i < t:
            r = r0 + i
            ti, o = divmod(r, self.TR)
            if ti == len(self.tiles):
                self.E.ensure_room(self, self.TR * self.ncols * 4)
                self.tiles.append(xp.zeros((self.TR, self.ncols), dtype=F32))
                self.loc.append('d')
            m = min(t - i, self.TR - o)
            self.tiles[ti][o:o + m] = L[i:i + m].astype(F32)
            i += m

    def compact(self, force=False):
        """drop dead physical columns, tile by tile (never a second full copy)."""
        dead = self.ncols - int(self.alive.sum())
        if dead == 0:
            return
        keep = np.flatnonzero(self.alive)
        kd = DV.to_dev(keep)
        for ti in range(len(self.tiles)):
            old = self.tiles[ti]
            if self.loc[ti] == 'd':
                self.E.ensure_room(self, self.TR * len(keep) * 4)
                new = old[:, kd]
            else:
                new = np.ascontiguousarray(old[:, keep])
            self.tiles[ti] = new
            del old
        self.pc2col = self.pc2col[keep]
        self.alive = np.ones(len(keep), dtype=bool)
        self.col2phys[:] = -1
        self.col2phys[self.pc2col] = np.arange(len(keep))
        self.compactions += 1
        DV.free_pool()

    # -- reading
    def source(self, rows):
        """batch rows (all with pivot weight < D) on the low part columns: n x (c - s0) device float64,
        entries in [0,p)."""
        n = len(rows)
        s0 = self.s0
        Lp = self.c - s0
        ps0 = int(np.searchsorted(self.pc2col, s0))
        Sphys = xp.zeros((n, self.ncols - ps0), dtype=F64)
        for ti, idx, off in self._groups(rows):
            Sphys[DV.to_dev(idx)] = self._rows(ti, off, slice(ps0, None)).astype(F64)
        S = xp.zeros((n, Lp), dtype=F64)
        S[:, DV.to_dev(self.pc2col[ps0:] - s0)] = Sphys
        del Sphys
        pv = self.piv[rows]
        assert np.all(pv >= s0), 'batch row is not low'
        S[DV.to_dev(np.arange(n)), DV.to_dev(pv - s0)] = 1.0
        return S

    def combine(self, rows, C):
        """C (t x len(rows), device float64 in [0,p)) times the basis rows ``rows``: t x c device array
        in part coordinates, entries in [0,p)."""
        t = C.shape[0]
        self.E.activate(self)
        M = xp.zeros((t, self.ncols), dtype=F64)
        self.accumulate(rows, C, M, 0, +1)
        Mf = xp.zeros((t, self.c), dtype=F64)
        Mf[:, DV.to_dev(self.pc2col)] = M
        del M
        Mf[:, DV.to_dev(self.piv[rows])] = C
        return Mf

    def rows_host(self, rows) -> np.ndarray:
        """basis rows as a host float64 matrix in part coordinates (small row sets only)."""
        rows = np.asarray(rows, dtype=np.int64)
        M = np.zeros((len(rows), self.c), dtype=F64)
        M[np.arange(len(rows)), self.piv[rows]] = 1.0
        for ti, idx, off in self._groups(rows):
            blk = self._rows_host(ti, off, slice(None)).astype(F64)
            M[np.ix_(idx, self.pc2col)] += blk
        return M

    # -- results
    def check_rref(self, sample: int = 64, seed: int = 0) -> bool:
        """sampled rows: dead columns exactly zero, entries before the pivot zero, entries in [0,p),
        every pivot column free of other rows' pivots (checked through the dead-column invariant)."""
        if self.rank == 0:
            return True
        rng = np.random.default_rng(seed)
        idx = rng.choice(self.rank, size=min(sample, self.rank), replace=False)
        M = np.zeros((len(idx), self.ncols))
        for ti, ii, off in self._groups(idx):
            M[ii] = self._rows_host(ti, off, slice(None)).astype(F64)
        if np.any(M[:, ~self.alive] != 0) or np.any(M < 0) or np.any(M >= self.p) or np.any(M != np.floor(M)):
            return False
        for r, row in enumerate(idx):
            before = self.pc2col < self.piv[row]
            if np.any(M[r, before] != 0):
                return False
        return True

    def profile(self, D: int) -> list:
        wp = self.w[self.piv]
        return [int(np.sum(wp <= d)) for d in range(D + 1)]

    def sha(self) -> str:
        """Hash of the canonical RREF (rows sorted by pivot column), same definition as the reference
        Part.sha(), computed from host copies streamed in row chunks."""
        self.compact()
        o = np.argsort(self.piv)
        h = hashlib.sha256()
        h.update(self.piv[o].astype(np.int64).tobytes())
        h.update(self.free_cols().astype(np.int64).tobytes())
        dt = np.uint16 if self.p < 65536 else np.uint32
        ch = max(1, (1 << 24) // max(1, self.ncols))
        for a in range(0, len(o), ch):
            rows = o[a:a + ch]
            blk = np.zeros((len(rows), self.ncols), dtype=np.float32)
            for ti, idx, off in self._groups(rows):
                blk[idx] = self._rows_host(ti, off, slice(None))
            h.update(np.ascontiguousarray(blk.astype(dt)).tobytes())
        return h.hexdigest()


# =============================================================================== engine

class Stop(Exception):
    pass


class Engine:
    def __init__(self, g: dict, k: int, p: int, D: int, device: str = 'cpu', split: bool = True,
                 gens: str = 'all', work_gb: float = 4.0, tile_mb: float = 256.0, tile_rows=None,
                 col_block=None, batch: int = 512, threads=None, compact_frac: float = 0.05,
                 stop_on_one: bool = False, t_random: int = 6, sweep_chunk: int = 256, ckpt=None,
                 ckpt_minutes: float = 0.0, ckpt_slots: int = 1, status=None, log=True, no_repair=False,
                 withhold=None, source_id: str = '', driver_sha: str = '', log_seconds: float = 20.0,
                 sub: int = 32, dev_gb=None, tile_budget_mb=None):
        use_device(device)
        assert p < (1 << 24), 'p must be < 2^24'
        if split:
            assert all(g.get((e2, e1), 0) == v for (e1, e2), v in g.items()), 'generator not swap-symmetric'
        self.g, self.k, self.p, self.D = dict(g), k, p, D
        self.device, self.split, self.gens = device, split, gens
        self.work_bytes = int(work_gb * (1 << 30))
        self.tile_bytes = int(tile_mb * (1 << 20))
        self.tile_rows_opt, self.col_block_opt = tile_rows, col_block
        self.batch_opt, self.compact_frac = batch, compact_frac
        self.stop_on_one, self.t_random, self.sweep_chunk = stop_on_one, t_random, sweep_chunk
        self.ckpt, self.ckpt_minutes, self.ckpt_slots = ckpt, ckpt_minutes, max(1, ckpt_slots)
        self.status_path, self.logging, self.no_repair = status, log, no_repair
        self.withhold = withhold          # test hook (part, k-th low row in processing order)
        self.source_id, self.driver_sha = source_id, driver_sha
        self.log_seconds, self.sub = log_seconds, sub
        self.threads = int(threads) if threads else default_threads()
        if _tpl is not None and device == 'cpu':
            try:
                _tpl(limits=self.threads)
            except Exception:  # pragma: no cover
                pass
        # device budget for spilling: tiles of the inactive part move to host numpy arrays when the tiles
        # (plus work buffers, one compaction tile and the context slack) would exceed dev_gb
        self.tile_budget = None
        self.spill_stats = dict(spills=0, fetches=0, d2h_bytes=0, h2d_bytes=0, d2h_s=0.0, h2d_s=0.0, switches=0)
        self.active_part = None
        self.dev_limit = None
        self.t_start = time.time()
        self.elapsed_prev = 0.0
        self.last_log = 0.0
        self.last_ckpt = time.time()
        self.term = False
        self.L = Layout(k, D, split)
        self.parts = [Part(P['w'], P['s0'], p, D, self, tile_rows) for P in self.L.parts]
        if DV.is_dev and (dev_gb or tile_budget_mb):
            comp = max(P.TR * P.c * 4 for P in self.parts)       # one tile of the largest width: compaction transient
            slack = 0 if device == 'fake-gpu' else int(1.5 * (1 << 30))
            if tile_budget_mb:
                self.tile_budget = int(tile_budget_mb * (1 << 20))
                self.dev_limit = self.tile_budget + self.work_bytes + comp
            else:
                self.dev_limit = int(dev_gb * (1 << 30))
                self.tile_budget = self.dev_limit - self.work_bytes - comp - slack
            if self.tile_budget < max(P.TR * P.c * 4 for P in self.parts):
                raise SystemExit(f'device budget too small: tile budget {self.tile_budget} bytes after work buffers, '
                                 f'compaction tile and slack')
            if device == 'fake-gpu':
                _Live.limit = self.dev_limit
        self.processed = [np.zeros(P.c, dtype=bool) for P in self.parts]
        self.phase = 'gens'
        self.gen_next = 0
        self.rounds = 0
        self.repairs = 0
        self.repair_rows = 0
        self.resumes = 0
        self.withhold_done = False
        self.contains_one = False
        self.closure_complete = False
        self.closed = False
        self.sweep = dict(stage='random', q=0, pi=0, pos=0, rows_checked=0, bad=0, random=None, sweeps=0)
        self.sweep_full_result = None
        self.sweep_random_result = None
        self.seq = 0
        self.ckpt_count = 0
        self.madds_sweep = 0
        self.madds_closure = None      # total madds when the closure phase ended (sweeps counted after)
        gdeg = max(bin(e1).count('1') + bin(e2).count('1') for e1, e2 in g)
        if gens == 'all':
            sel = np.flatnonzero(self.L.w <= D - gdeg)
        else:
            sel = np.array([int(np.flatnonzero(self.L.w == 0)[0])])
        self.gen_cols = sel
        self.generator_rows = len(sel)

    # ----- budget
    def batch_rows(self, pi: int) -> int:
        """rows per batch from the work budget: R, filtered R and the in-place RREF (3 copies of n x
        ncols float64 plus a boolean), S, the image block and the gathered inner blocks (n x low
        columns, float64, about 6 copies), and 25% reserved for the column-block temporaries."""
        # physical columns <= free / (1 - compact_frac) (compaction fires when dead > frac * ncols); using
        # this logical bound instead of ncols makes the batch sizes independent of the lazy layout, so a
        # resumed run (which restarts compacted) follows the same sequence of states
        fr = min(self.compact_frac, 0.5)
        nc = max(int(np.ceil((P.c - P.rank) / (1.0 - fr))) for P in self.parts)
        nc = max(1, min(nc, max(P.c for P in self.parts)))
        low = max(P.c - P.s0 for P in self.parts)
        per = 8 * (3 * nc + 10 * low) + nc
        n = int((0.75 * self.work_bytes - (128 << 20)) // max(1, per))
        return max(1, min(self.batch_opt, n))

    def col_block(self, n_rows: int, tr: int) -> int:
        """column block width: temporaries (tr + n_rows) x cw float64 (gathered tile rows and the
        product) stay within 25% of the work budget."""
        if self.col_block_opt:
            return int(self.col_block_opt)
        return max(64, int(0.25 * self.work_bytes // (8 * (tr + n_rows))))

    # ----- spilling of tiles to the host
    def dev_tile_bytes(self) -> int:
        return sum(P.dev_bytes() for P in self.parts)

    def _move(self, part, ti, to):
        t0 = time.perf_counter()
        t = part.tiles[ti]
        n = t.nbytes
        if to == 'h':
            h = DV.to_host(t)
            DV.sync()
            part.tiles[ti], part.loc[ti] = h, 'h'
            del t
            DV.free_pool()
            st = self.spill_stats
            st['spills'] += 1
            st['d2h_bytes'] += n
            st['d2h_s'] += time.perf_counter() - t0
        else:
            d = DV.to_dev(t)
            DV.sync()
            part.tiles[ti], part.loc[ti] = d, 'd'
            st = self.spill_stats
            st['fetches'] += 1
            st['h2d_bytes'] += n
            st['h2d_s'] += time.perf_counter() - t0

    def ensure_room(self, part, nbytes: int):
        """before an allocation of nbytes of tile memory for ``part``: spill tiles of the OTHER parts
        (highest tile index first) until the tile budget holds."""
        if self.tile_budget is None:
            return
        while self.dev_tile_bytes() + nbytes > self.tile_budget:
            cand = [(Q, ti) for Q in self.parts if Q is not part for ti in range(len(Q.tiles) - 1, -1, -1) if Q.loc[ti] == 'd']
            if not cand:
                raise MemoryError(f'device tile budget {self.tile_budget} bytes cannot hold one part '
                                  f'({self.dev_tile_bytes()} on the device + {nbytes} requested)')
            Q, ti = cand[0]
            self._move(Q, ti, 'h')

    def activate(self, part):
        """make ``part`` the target: all its tiles on the device, spilling just enough tiles of the other
        part(s) first.  No-op without a budget or when the part is already complete on the device."""
        if self.tile_budget is None:
            return
        host = [ti for ti, l in enumerate(part.loc) if l == 'h']
        if self.active_part is not part:
            self.active_part = part
            self.spill_stats['switches'] += 1
        if not host:
            return
        need = sum(part.tiles[ti].nbytes for ti in host)
        self.ensure_room(part, need)
        for ti in host:
            self._move(part, ti, 'd')

    # ----- logging
    def seconds(self):
        return self.elapsed_prev + time.time() - self.t_start

    def snapshot(self, **extra):
        used, total = DV.mem_info()
        ml = {}
        for key in ('reduce', 'reduce_dense', 'update', 'rref'):
            ml[key] = sum(P.madds[key] for P in self.parts)
        d = dict(t=round(self.seconds(), 1), phase=self.phase, maxrss_mb=maxrss_mb(),
                 dev_used_mb=used >> 20, dev_total_mb=total >> 20, dev_pool_mb=DV.pool_bytes() >> 20,
                 ranks=[P.rank for P in self.parts], low=[int(len(P.low_rows())) for P in self.parts],
                 processed=[int(s.sum()) for s in self.processed], free=[P.nfree for P in self.parts],
                 ncols_phys=[P.ncols for P in self.parts],
                 madds=ml['reduce'] + ml['update'] + ml['rref'] + self.madds_sweep, madds_detail=ml,
                 rounds=self.rounds, repairs=self.repairs)
        if self.tile_budget is not None:
            d['spill'] = {k: (round(v, 2) if isinstance(v, float) else v) for k, v in self.spill_stats.items()}
            d['dev_tile_mb'] = self.dev_tile_bytes() >> 20
            d['host_tile_mb'] = sum(P.host_bytes() for P in self.parts) >> 20
        sec = self.seconds()
        d['rate_madds_per_s'] = round(d['madds'] / sec) if sec > 0 else 0
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

    # ----- helpers
    def total_madds(self) -> int:
        return sum(P.madds[key] for P in self.parts for key in ('reduce', 'update', 'rref')) + self.madds_sweep

    def zero_one(self) -> bool:
        P = self.parts[0]
        return P.zero_col >= 0 and bool(P.col2row[P.zero_col] >= 0)

    def check_term(self):
        """called between batches: checkpoint and exit on SIGTERM."""
        if self.tile_budget is not None:      # accounting invariant, checked between batches
            assert self.dev_tile_bytes() <= self.tile_budget, 'tiles on the device exceed the tile budget'
            if self.device == 'fake-gpu':
                assert _Live.cur >= self.dev_tile_bytes(), 'fake device accounting lost a tile'
        if self.term:
            self.emit(force=True, event='sigterm')
            if self.ckpt:
                self.save_checkpoint('sigterm')
            raise Stop()

    def maybe_checkpoint(self):
        if self.ckpt and self.ckpt_minutes and time.time() - self.last_ckpt >= 60 * self.ckpt_minutes:
            self.save_checkpoint('timer')

    # ----- generators
    def _gen_chunk(self, cols_idx):
        N = 1 << self.k
        L = self.L
        terms = list(self.g.items())
        e1 = np.array([t[0][0] for t in terms], dtype=np.int64)
        e2 = np.array([t[0][1] for t in terms], dtype=np.int64)
        vals = np.array([t[1] for t in terms], dtype=np.int64)
        b1 = L.a1[cols_idx][:, None]
        b2 = L.a2[cols_idx][:, None]
        full = L.index((b1 + e1[None, :]) % N, (b2 + e2[None, :]) % N)          # chunk x nterms
        n = len(cols_idx)
        for q, Q in enumerate(L.parts):
            coef = Q['proj_coef'][full] * vals[None, :]
            m = Q['lift_idx'][full]
            ok = Q['proj_coef'][full] != 0
            if not ok.any():
                continue
            T = np.unique(m[ok])
            Y = np.zeros((n, len(T)), dtype=np.int64)
            rr = np.broadcast_to(np.arange(n)[:, None], full.shape)
            np.add.at(Y, (rr[ok], np.searchsorted(T, m[ok])), coef[ok])
            Y %= self.p
            yield q, T, DV.to_dev(Y.astype(F64))

    def do_generators(self):
        n = len(self.gen_cols)
        while self.gen_next < n:
            a = self.gen_next
            for q, T, Y in self._gen_chunk(self.gen_cols[a:a + 256]):
                R = self.parts[q].residual(T, Y)
                del Y
                if R is not None:
                    self.parts[q].add_residuals(R)
                    del R
            self.gen_next = min(n, a + 256)
            self.emit(stage='gens', gen_next=self.gen_next, of=n)
            if self.stop_on_one and self.zero_one():
                raise Stop('one')
            self.check_term()
            self.maybe_checkpoint()
        self.phase = 'closure'

    # ----- images of a source block
    def images(self, pi: int, S, targets=None):
        """yield (q, T, Yt) for every target part and shift, TARGET-MAJOR (all images into part 0, then all
        into part 1), so that the target part changes once per batch: the shifted, projected rows of S (n x
        low columns of part pi), as dense blocks on their common support T (entries in [0,p))."""
        n = S.shape[0]
        s0 = self.parts[pi].s0
        nzs = DV.to_host(xp.any(S != 0, axis=0))
        for q in (range(len(self.parts)) if targets is None else targets):
            self.activate(self.parts[q])
            for ij in self.L.shift_list:
                src, m, coef, layer = self.L.shift_map(pi, ij, q)
                kp = nzs[src - s0]
                src, m, coef, layer = src[kp], m[kp], coef[kp], layer[kp]
                if len(m) == 0:
                    continue
                T = np.unique(m)
                Yt = xp.zeros((n, len(T)), dtype=F64)
                for ly in (0, 1):
                    s = layer == ly
                    if not s.any():
                        continue
                    pos = np.searchsorted(T, m[s])
                    Yt[:, DV.to_dev(pos)] += S[:, DV.to_dev(src[s] - s0)] * DV.to_dev(coef[s])[None, :]
                modp_(Yt, self.p)
                yield q, T, Yt

    # ----- closure
    def next_pending(self):
        for pi, P in enumerate(self.parts):
            low = P.low_rows()
            if len(low) == 0:
                continue
            pend = low[~self.processed[pi][P.piv[low]]]
            if len(pend):
                return pi, pend
        return None

    def process_batch(self, pi: int, rows):
        P = self.parts[pi]
        for r in rows:
            self.processed[pi][P.piv[r]] = True
        S = P.source(rows)
        for q, T, Yt in self.images(pi, S):
            R = self.parts[q].residual(T, Yt)
            del Yt
            if R is not None:
                self.parts[q].add_residuals(R)
                del R
            if self.stop_on_one and self.zero_one():
                raise Stop('one')
        del S
        self.rounds += 1

    def run_closure(self):
        while True:
            self.check_term()
            pend = self.next_pending()
            if pend is None:
                self.closure_complete = True
                return
            pi, rows = pend
            if self.withhold and not self.withhold_done and pi == self.withhold[0] and len(rows) > self.withhold[1]:
                r = rows[self.withhold[1]]          # test hook: mark one low row processed without shifting
                self.processed[pi][self.parts[pi].piv[r]] = True
                self.withhold_done = True
                self.emit(force=True, event='withheld', part=pi, pivot=int(self.parts[pi].piv[r]))
                rows = np.delete(rows, self.withhold[1])
                if len(rows) == 0:
                    continue
            rows = rows[:self.batch_rows(pi)]
            self.process_batch(pi, rows)
            self.emit(stage='closure')
            if self.stop_on_one and self.zero_one():
                raise Stop('one')
            self.maybe_checkpoint()

    # ----- certificate
    def _check_images(self, pi, S, repair: bool, targets=None):
        """reduce the images of S (into the target parts ``targets``, default all) against the current
        basis; nonzero residual rows are counted and, if ``repair``, added.  Returns their number."""
        bad = 0
        for q, T, Yt in self.images(pi, S, targets):
            R = self.parts[q].residual(T, Yt, count_dense=False)
            del Yt
            if R is not None:
                bad += R.shape[0]
                if repair:
                    self.parts[q].add_residuals(R)
                del R
        return bad

    def sweep_random_run(self):
        """t random combinations of the low rows of each part, all shifts and projections."""
        rng = np.random.default_rng(1)
        bad = 0
        if self.t_random <= 0:        # test setting: skip the randomized check (the full sweep still runs)
            return dict(t=0, seed=1, nonzero_residuals=0, passed=True, maps=0)
        srcs = []
        for pi, P in enumerate(self.parts):
            low = P.low_rows()
            if len(low) == 0:
                continue
            C = DV.to_dev(rng.integers(0, self.p, size=(self.t_random, len(low))).astype(F64))
            M = P.combine(low, C)
            assert not DV.to_host(xp.any(M[:, :P.s0] != 0)), 'low combination with high support'
            srcs.append((pi, xp.array(M[:, P.s0:])))
            del M
        for q in range(len(self.parts)):          # target-major: the target part changes once
            for pi, S in srcs:
                bad += self._check_images(pi, S, repair=not self.no_repair, targets=[q])
        return dict(t=self.t_random, seed=1, nonzero_residuals=bad, passed=bad == 0,
                    maps=len(self.parts) ** 2 * len(self.L.shift_list))

    def run_sweeps(self):
        sw = self.sweep
        while True:
            self.check_term()
            if sw['stage'] == 'random':
                before = [P.rank for P in self.parts]
                res = self.sweep_random_run()
                self.sweep_random_result = res
                self.emit(force=True, stage='sweep_random', result=res)
                if not res['passed']:
                    if self.no_repair:
                        self.closed = False
                        self.sweep_full_result = dict(complete=False, rows_checked=0, nonzero_residuals=res['nonzero_residuals'],
                                                      passed=False, note='randomized sweep failed (no repair)')
                        return
                    self.repairs += 1
                    self.repair_rows += sum(P.rank for P in self.parts) - sum(before)
                    self.phase = 'closure'
                    self.closure_complete = False
                    self.run_closure()
                    self.phase = 'sweep'
                    continue
                sw.update(stage='full', q=0, pi=0, pos=0, rows_checked=0, bad=0, sweeps=sw['sweeps'] + 1)
            if sw['stage'] == 'full':
                dirty = False
                # TARGET-MAJOR order: all source rows of both parts against target part 0, then against part 1
                # (the basis does not change during a clean sweep; a nonzero residual ends the pass)
                while sw['q'] < len(self.parts):
                  while sw['pi'] < len(self.parts):
                    pi = sw['pi']
                    P = self.parts[pi]
                    low = P.low_rows()
                    while sw['pos'] < len(low):
                        rows = low[sw['pos']:sw['pos'] + min(self.sweep_chunk, self.batch_rows(pi))]
                        S = P.source(rows)
                        rk = sum(Q.rank for Q in self.parts)
                        bad = self._check_images(pi, S, repair=not self.no_repair, targets=[sw['q']])
                        del S
                        sw['bad'] += bad
                        sw['rows_checked'] += len(rows)
                        sw['pos'] += len(rows)
                        if bad:
                            dirty = True
                            self.repair_rows += sum(Q.rank for Q in self.parts) - rk
                            break
                        self.emit(stage='sweep_full', sweep_q=sw['q'], sweep_pi=pi, sweep_pos=sw['pos'], of=len(low))
                        self.maybe_checkpoint()
                        self.check_term()
                    if dirty:
                        break
                    sw['pi'] += 1
                    sw['pos'] = 0
                  if dirty:
                    break
                  sw['q'] += 1
                  sw['pi'] = 0
                  sw['pos'] = 0
                if not dirty:
                    self.sweep_full_result = dict(complete=True, rows_checked=sw['rows_checked'],
                                                  nonzero_residuals=0, passed=True, sweeps_run=sw['sweeps'])
                    self.closed = True
                    return
                if self.no_repair:
                    self.sweep_full_result = dict(complete=False, rows_checked=sw['rows_checked'],
                                                  nonzero_residuals=sw['bad'], passed=False, sweeps_run=sw['sweeps'])
                    self.closed = False
                    return
                self.repairs += 1
                self.emit(force=True, event='repair', nonzero_rows=sw['bad'])
                # reset the sweep BEFORE the closure continues, so a checkpoint taken inside it resumes
                # with a fresh sweep over the repaired basis
                sw.update(stage='full', q=0, pi=0, pos=0, rows_checked=0, bad=0, sweeps=sw['sweeps'] + 1)
                self.phase = 'closure'
                self.closure_complete = False
                self.run_closure()
                self.phase = 'sweep'

    # ----- main
    def run(self) -> str:
        """returns 'done', 'one' (stopped on 1 in V_D), 'sigterm', or 'failed_sweep' (no-repair)."""
        try:
            if self.phase == 'gens':
                self.do_generators()
            if self.phase == 'closure':
                self.run_closure()
                if self.madds_closure is None:
                    self.madds_closure = self.total_madds()
                if self.ckpt:
                    self.phase = 'sweep'
                    self.save_checkpoint('closure_done')
                self.phase = 'sweep'
            if self.phase == 'sweep':
                self.run_sweeps()
            self.phase = 'done'
        except Stop as e:
            if e.args and e.args[0] == 'one':
                self.contains_one = True
                self.closure_complete = False
                self.phase = 'stopped_on_one'
                self.emit(force=True, event='contains_one')
                return 'one'
            return 'sigterm'
        self.contains_one = self.zero_one()
        if self.sweep_full_result and not self.sweep_full_result.get('passed'):
            return 'failed_sweep'
        return 'done'

    # ----- results
    def profile(self):
        pr = [P.profile(self.D) for P in self.parts]
        return [sum(x) for x in zip(*pr)], pr

    def rank(self):
        return sum(P.rank for P in self.parts)

    def record(self) -> dict:
        prof, parts = self.profile()
        ml = {key: sum(P.madds[key] for P in self.parts) for key in ('reduce', 'reduce_dense', 'update', 'rref')}
        rec = dict(k=self.k, D=self.D, p=self.p, split=self.split, gens=self.gens, columns=self.L.c,
                   generator_rows=self.generator_rows, rank=self.rank(), profile=prof, part_profiles=parts,
                   part_ranks=[P.rank for P in self.parts], contains_one=self.contains_one,
                   closure_complete=self.closure_complete, closed=self.closed, rounds=self.rounds,
                   seconds=round(self.seconds(), 1),
                   madds=ml['reduce'] + ml['update'] + ml['rref'] + self.madds_sweep, madds_detail=ml,
                   madds_dense=ml['reduce_dense'] + ml['update'] + ml['rref'],
                   madds_closure_phase=self.madds_closure,
                   seconds_detail={key: round(sum(P.secs[key] for P in self.parts), 2) for key in ('reduce', 'update', 'rref')},
                   peak_entries=[P.peak_entries for P in self.parts],
                   peak_phys_entries=[P.peak_phys_entries for P in self.parts],
                   compactions=[P.compactions for P in self.parts],
                   rref_ok=all(P.check_rref() for P in self.parts),
                   basis_sha256=[P.sha() for P in self.parts],
                   sweep_random=self.sweep_random_result, sweep_full=self.sweep_full_result,
                   repairs=self.repairs, repair_rows=self.repair_rows, resumes=self.resumes,
                   device=self.device, threads=self.threads, blas_threads_env=os.environ.get('OPENBLAS_NUM_THREADS'), maxrss_mb=maxrss_mb(),
                   engine_sha256=ENGINE_SHA256, engine_version=ENGINE_VERSION,
                   work_gb=self.work_bytes / (1 << 30), tile_mb=self.tile_bytes / (1 << 20),
                   tile_rows=[P.TR for P in self.parts], batch=self.batch_opt, compact_frac=self.compact_frac,
                   t_random=self.t_random, sweep_chunk=self.sweep_chunk, stop_on_one=self.stop_on_one,
                   no_repair=self.no_repair, compress=False,
                   spill=dict(self.spill_stats, tile_budget_bytes=self.tile_budget, dev_limit_bytes=self.dev_limit,
                              host_tile_bytes=sum(P.host_bytes() for P in self.parts)))
        used, total = DV.mem_info()
        rec.update(dev_used_mb=used >> 20, dev_total_mb=total >> 20, dev_peak_live_mb=(_Live.peak >> 20) if self.device == 'fake-gpu' else None)
        return rec

    # ----- soundness / evaluation checks (mirror the reference driver)
    def random_elements_full(self, t: int, seed: int) -> np.ndarray:
        rng = np.random.default_rng(seed)
        out = np.zeros((t, self.L.c), dtype=F64)
        for pi, P in enumerate(self.parts):
            C = DV.to_dev(rng.integers(0, self.p, size=(t, P.rank)).astype(F64))
            M = DV.to_host(P.combine(np.arange(P.rank), C))
            out = np.mod(out + self.L.lift_full(pi, M), self.p)
        return out

    # ----- checkpoint / resume
    def _state_json(self):
        return dict(phase=self.phase, gen_next=self.gen_next, rounds=self.rounds, repairs=self.repairs,
                    repair_rows=self.repair_rows, resumes=self.resumes, withhold_done=self.withhold_done,
                    closure_complete=self.closure_complete, sweep=self.sweep,
                    sweep_random_result=self.sweep_random_result, elapsed=self.seconds(),
                    madds=[P.madds for P in self.parts], secs=[P.secs for P in self.parts],
                    peak_entries=[P.peak_entries for P in self.parts],
                    peak_phys_entries=[P.peak_phys_entries for P in self.parts],
                    compactions=[P.compactions for P in self.parts], tile_rows=[P.TR for P in self.parts],
                    madds_sweep=self.madds_sweep, madds_closure=self.madds_closure, seq=self.seq,
                    spill_stats=self.spill_stats)

    def _config_json(self):
        return dict(k=self.k, p=self.p, D=self.D, split=self.split, gens=self.gens, source_id=self.source_id)

    def save_checkpoint(self, reason: str):
        t0 = time.time()
        root = self.ckpt
        os.makedirs(root, exist_ok=True)
        for P in self.parts:
            P.compact()
        self.seq += 1
        tmp = os.path.join(root, f'.tmp-{self.seq:08d}')
        final = os.path.join(root, f'ck-{self.seq:08d}')
        if os.path.exists(tmp):
            shutil.rmtree(tmp)
        os.makedirs(tmp)
        files = {}

        def fsync_file(path):
            fd = os.open(path, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)

        def put(name):
            path = os.path.join(tmp, name)
            fsync_file(path)
            files[name] = [os.path.getsize(path), file_sha256(path)]

        sleep = float(os.environ.get('ICBIG_TEST_CKPT_SLEEP', '0') or 0)     # test hook (kill -9 mid-write)
        if self.seq < int(os.environ.get('ICBIG_TEST_CKPT_SLEEP_FROM', '1')):
            sleep = 0.0
        for q, P in enumerate(self.parts):
            for ti, tile in enumerate(P.tiles):
                nr = P.tile_rows_filled(ti)
                np.save(os.path.join(tmp, f'p{q}_t{ti:05d}.npy'), DV.to_host(tile[:nr]) if P.loc[ti] == 'd' else np.asarray(tile[:nr]))
                put(f'p{q}_t{ti:05d}.npy')
                if sleep and ti == 0:
                    time.sleep(sleep)
            np.savez(os.path.join(tmp, f'p{q}_host.npz'), piv=P.piv, processed=self.processed[q])
            put(f'p{q}_host.npz')
        man = dict(version=1, engine_sha256=ENGINE_SHA256, driver_sha256=self.driver_sha, reason=reason,
                   config=self._config_json(), state=self._state_json(), files=files,
                   ntiles=[len(P.tiles) for P in self.parts])
        mp = os.path.join(tmp, 'manifest.json.part')
        with open(mp, 'w') as f:
            json.dump(man, f)
            f.flush()
            os.fsync(f.fileno())
        os.replace(mp, os.path.join(tmp, 'manifest.json'))
        fsync_file(tmp)
        os.rename(tmp, final)
        fsync_file(root)
        olds = sorted(d for d in os.listdir(root) if d.startswith('ck-'))
        for d in olds[:-self.ckpt_slots]:
            shutil.rmtree(os.path.join(root, d))
        fsync_file(root)
        self.last_ckpt = time.time()
        self.ckpt_count += 1
        self.emit(force=True, event='checkpoint', reason=reason, seq=self.seq, ckpt_seconds=round(time.time() - t0, 1))

    @staticmethod
    def find_checkpoint(root, verify=True):
        """newest complete checkpoint directory (manifest present, every file's size and sha256 match)."""
        if not os.path.isdir(root):
            return None, None
        for d in os.listdir(root):
            if d.startswith('.tmp-'):
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)
        for d in sorted((d for d in os.listdir(root) if d.startswith('ck-')), reverse=True):
            mp = os.path.join(root, d, 'manifest.json')
            if not os.path.exists(mp):
                continue
            try:
                man = json.load(open(mp))
                ok = True
                for name, (size, sha) in man['files'].items():
                    path = os.path.join(root, d, name)
                    if not os.path.exists(path) or os.path.getsize(path) != size or (verify and file_sha256(path) != sha):
                        ok = False
                        break
            except Exception:
                ok = False
            if ok:
                return os.path.join(root, d), man
        return None, None

    def load_checkpoint(self, force=False) -> bool:
        path, man = self.find_checkpoint(self.ckpt)
        if path is None:
            return False
        if man['config'] != self._config_json() and not force:
            raise SystemExit(f'checkpoint source mismatch: {man["config"]} != {self._config_json()} (use --force-resume)')
        if (man['engine_sha256'] != ENGINE_SHA256 or man['driver_sha256'] != self.driver_sha) and not force:
            raise SystemExit('checkpoint engine/driver sha256 differs (use --force-resume)')
        st = man['state']
        self.phase, self.gen_next, self.rounds = st['phase'], st['gen_next'], st['rounds']
        self.repairs, self.repair_rows, self.resumes = st['repairs'], st['repair_rows'], st['resumes'] + 1
        self.withhold_done, self.closure_complete = st['withhold_done'], st['closure_complete']
        self.sweep, self.sweep_random_result = st['sweep'], st['sweep_random_result']
        self.elapsed_prev, self.madds_sweep, self.seq = st['elapsed'], st['madds_sweep'], st['seq']
        self.madds_closure = st['madds_closure']
        self.spill_stats = st['spill_stats']
        self.t_start = time.time()
        for q, P in enumerate(self.parts):
            hz = np.load(os.path.join(path, f'p{q}_host.npz'))
            P.piv = hz['piv'].astype(np.int64)
            self.processed[q] = hz['processed'].astype(bool)
            P.TR = st['tile_rows'][q]
            ispiv = np.zeros(P.c, dtype=bool)
            ispiv[P.piv] = True
            P.pc2col = np.flatnonzero(~ispiv)
            P.alive = np.ones(len(P.pc2col), dtype=bool)
            P.col2phys = -np.ones(P.c, dtype=np.int64)
            P.col2phys[P.pc2col] = np.arange(len(P.pc2col))
            P.col2row = -np.ones(P.c, dtype=np.int64)
            P.col2row[P.piv] = np.arange(len(P.piv))
            P.tiles, P.loc = [], []
            for ti in range(man['ntiles'][q]):
                a = np.load(os.path.join(path, f'p{q}_t{ti:05d}.npy'))
                if self.tile_budget is not None:      # with a device budget every tile starts on the host
                    tile = np.zeros((P.TR, len(P.pc2col)), dtype=F32)
                    tile[:a.shape[0]] = a
                    P.tiles.append(tile)
                    P.loc.append('h')
                else:
                    tile = xp.zeros((P.TR, len(P.pc2col)), dtype=F32)
                    tile[:a.shape[0]] = DV.to_dev(a)
                    P.tiles.append(tile)
                    P.loc.append('d')
            P.madds, P.secs = st['madds'][q], st['secs'][q]
            P.peak_entries, P.peak_phys_entries = st['peak_entries'][q], st['peak_phys_entries'][q]
            P.compactions = st['compactions'][q]
        self.emit(force=True, event='resumed', from_ckpt=os.path.basename(path))
        return True

    def install_sigterm(self):
        def h(signum, frame):
            self.term = True
        signal.signal(signal.SIGTERM, h)


# =============================================================================== evaluation (host)

def evaluate_full(L: Layout, rows: np.ndarray, pts: list, p: int, k: int) -> np.ndarray:
    """values of full-coordinate rows at points (E1, E2): monomial (a1, a2) -> h^(a1 E1 + a2 E2)."""
    h, H = subgroup(p, k)
    Hv = np.array(H, dtype=np.int64)
    N = 1 << k
    out = np.zeros((rows.shape[0], len(pts)), dtype=np.int64)
    R = rows.astype(np.int64) % p
    for z, (E1, E2) in enumerate(pts):
        ev = Hv[(L.a1 * E1 + L.a2 * E2) % N]
        out[:, z] = (R * ev % p).sum(axis=1) % p
    return out


def eval_rank(L: Layout, pts: list, p: int, k: int, maxw: int, eng_kwargs: dict, batch: int = 256) -> int:
    """rank of the evaluation of the monomials of weight <= maxw at the points (via the Part machinery)."""
    if not pts:
        return 0
    h, H = subgroup(p, k)
    Hv = np.array(H, dtype=np.int64)
    N = 1 << k
    cols = np.flatnonzero(L.w <= maxw)
    class _E:      # minimal engine stand-in for Part
        tile_bytes = eng_kwargs.get('tile_bytes', 64 << 20)
        compact_frac = 0.05
        tile_budget = None
        work_bytes = eng_kwargs.get('work_bytes', 1 << 30)
        def activate(self, part):
            pass
        def ensure_room(self, part, nbytes):
            pass
        def col_block(self, n_rows, tr):
            return max(64, int(0.25 * self.work_bytes // (8 * (tr + n_rows))))
    P = Part(np.ones(len(cols), dtype=np.int64), len(cols), p, 0, _E(), eng_kwargs.get('tile_rows'))
    T = np.arange(len(cols))
    for a in range(0, len(pts), batch):
        M = np.zeros((min(batch, len(pts) - a), len(cols)), dtype=F64)
        for z in range(M.shape[0]):
            E1, E2 = pts[a + z]
            M[z] = Hv[(L.a1[cols] * E1 + L.a2[cols] * E2) % N]
        R = P.residual(T, DV.to_dev(M), count_dense=False)
        if R is not None:
            P.add_residuals(R)
    return P.rank


# =============================================================================== plan

def plan(k: int, D: int, device: str, work_gb: float, tile_mb: float, split: bool = True,
         compact_frac: float = 0.05, batch: int = 512, dev_gb=None) -> dict:
    """predicted peak memory: worst case of the two matrices (every rank, including lazily deleted
    columns and one partly filled tile per part) plus working buffers."""
    L = Layout(k, D, split)
    cs = [P['c'] for P in L.parts]
    lows = [P['c'] - P['s0'] for P in L.parts]
    worst_each, peak_rank = [], []
    for c in cs:
        best, br = 0, 0
        step = max(1, c // 4000)
        for r in range(0, c + 1, step):
            ncols_phys = (c - r) * (1 + compact_frac) if r else c        # lazy columns up to compact_frac
            ncols_phys = min(c, ncols_phys + 1)
            ent = r * ncols_phys * 4
            if ent > best:
                best, br = ent, r
        worst_each.append(best)
        peak_rank.append(br)
    tile_slack = len(cs) * tile_mb * (1 << 20)
    matrices = sum(worst_each) + tile_slack
    work = work_gb * (1 << 30)
    nc = max(cs)
    low = max(lows)
    per = 8 * (3 * nc + 10 * low) + nc
    nrows = max(1, min(batch, int((0.75 * work - (128 << 20)) // per)))
    ctx = 1.5 * (1 << 30)       # CUDA context, cuBLAS workspace, index arrays: allowance
    compaction = tile_mb * (1 << 20)      # compaction holds an old and a new tile while R is alive
    total = matrices + work + compaction + (ctx if device != 'cpu' else 0)
    spill = None
    if dev_gb:
        dev = dev_gb * (1 << 30)
        ctx_ = ctx if device == 'gpu' or device == 'cpu' else 0
        tile_budget = dev - work - compaction - ctx_
        target_need = max(worst_each) + tile_mb * (1 << 20)      # the target part wholly on the device, one tile of slack
        host_spill = max(0.0, matrices - tile_budget)
        spill = dict(dev_gb=dev_gb, tile_budget_bytes=tile_budget, tile_budget_gb=round(tile_budget / (1 << 30), 1),
                     target_part_worst_gb=round(target_need / (1 << 30), 1), feasible=bool(target_need <= tile_budget),
                     predicted_device_peak_gb=round((min(matrices, tile_budget) + work + compaction + ctx_) / (1 << 30), 1),
                     predicted_host_spill_worst_gb=round(host_spill / (1 << 30), 1),
                     note='both parts at their worst rank at the same time; the host holds the part of the matrices that does '
                          'not fit the tile budget (the inactive part, tiles from the highest index), a switch of the target part '
                          'moves at most that amount in each direction')
    return dict(k=k, D=D, device=device, columns=L.c, spill=spill, part_columns=cs, low_columns=lows,
                matrix_worst_bytes_per_part=worst_each, matrix_worst_rank_per_part=peak_rank,
                tile_slack_bytes=tile_slack, matrices_bytes=matrices, work_bytes=work, compaction_transient_bytes=compaction,
                batch_rows_at_start=nrows, context_allowance_bytes=ctx if device != 'cpu' else 0,
                predicted_peak_gb=round(total / (1 << 30), 1),
                matrices_gb=round(matrices / (1 << 30), 1), work_gb=work_gb,
                note='matrices: max over rank r of r x (free + lazily deleted) float32 per part (+ one tile per '
                     'part); the final rank is below this worst case, which is reached at rank about half the columns')


# =============================================================================== arithmetic self-tests

def selftest_arith(device_name=None, seed: int = 12345) -> dict:
    """(a) modp_ against Python integers (random data and extremes), several primes;
    (b) float64 products with exact result near 2^53 against int64 host arithmetic;
    (c) matmul_sub with the inner-dimension split."""
    if device_name:
        use_device(device_name)
    rng = np.random.default_rng(seed)
    out = dict(modp_cases=0, modp_failures=0, product_cases=0, product_failures=0, matmul_sub_failures=0)
    for p in (3, 12289, 40961, 65537, 16777213, 16777199):
        ext = [0, 1, -1, p, -p, p - 1, -(p - 1), p + 1, 1 << 50, -(1 << 50), (1 << 50) - 1, -(1 << 50) + 1,
               ((1 << 50) // p) * p, ((1 << 50) // p) * p - 1, ((1 << 50) // p) * p + 1,
               -((1 << 50) // p) * p, -((1 << 50) // p) * p + 1, -((1 << 50) // p) * p - 1]
        rnd = rng.integers(-(1 << 50), (1 << 50) + 1, size=200000, dtype=np.int64)
        mult = rng.integers(-(1 << 50) // p, (1 << 50) // p, size=100000, dtype=np.int64) * p
        near = np.concatenate([mult - 1, mult, mult + 1])
        near = near[np.abs(near) <= (1 << 50)]
        small = rng.integers(-3 * p, 3 * p, size=50000, dtype=np.int64)
        vals = np.concatenate([np.array(ext, dtype=np.int64), rnd, near, small])
        X = DV.to_dev(vals.astype(F64).reshape(-1, 1))
        modp_(X, p)
        got = DV.to_host(X).reshape(-1).astype(np.int64)
        want = np.array([int(v) % p for v in vals], dtype=np.int64)
        out['modp_cases'] += len(vals)
        out['modp_failures'] += int(np.count_nonzero(got != want))
    # (b) products near 2^53
    for p in (65537, 40961, 16777213):
        m = 3
        hi = p - 1
        n = ((1 << 53) - 1) // (hi * hi) - 1
        A = hi - rng.integers(0, 3, size=(m, n))
        B = hi - rng.integers(0, 3, size=(n, m))
        exact = A.astype(np.int64) @ B.astype(np.int64)
        assert int(exact.max()) < (1 << 53)
        got = DV.to_host(DV.to_dev(A.astype(F64)) @ DV.to_dev(B.astype(F64)))
        out['product_cases'] += 1
        out['product_failures'] += int(np.any(got != exact.astype(F64)))
        out.setdefault('product_max_exact', []).append(int(exact.max()))
    # (c) matmul_sub with the split of the inner dimension
    for p in (65537, 16777213):
        mi = max_inner(p)
        n = 3 * mi + 5 if mi < 100000 else 70000
        A = rng.integers(0, p, size=(2, n)).astype(F64)
        B = rng.integers(0, p, size=(n, 3)).astype(F64)
        acc = rng.integers(0, p, size=(2, 3)).astype(F64)
        want = (acc.astype(object) - A.astype(np.int64).astype(object) @ B.astype(np.int64).astype(object)) % p
        got = matmul_sub(DV.to_dev(acc.copy()), DV.to_dev(A), DV.to_dev(B), p, mi)
        out['matmul_sub_failures'] += int(np.any(DV.to_host(got).astype(object) != want))
    out['guard_hits'] = GUARD_HITS[0]
    out['passed'] = out['modp_failures'] == 0 and out['product_failures'] == 0 and out['matmul_sub_failures'] == 0
    return out
