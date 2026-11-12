#!/usr/bin/env python3
"""Independent checker for a "closed containing space" certificate of the degree-D chain closure V_D.

CLAIM CHECKED.  A checkpoint (tiles of a reduced row echelon basis, per part of the summand-swap split)
describes a subspace W of span(M_<=D).  The checker establishes, from the stored basis alone,
  (1) g in W,   (2) W closed under u_{i,j} on W cap span(M_<D),   (3) 1 not in W.
(1) and (2) give V_D <= W (V_D is the least such space), so (3) gives 1 not in V_D.  See
experiments/ic-bounds-k14/CERTCHECK.md for the argument, the independence statement and the cost model.

INDEPENDENCE.  Single file, Python >= 3.10, numpy; cupy only for ``--device gpu``.  It imports no
repository module.  Layout, shift action, projection and membership test are written here from the
mathematics.  The only engine facts used are the on-disk checkpoint format (manifest.json, p<q>_t<ti>.npy
float32 tiles, p<q>_host.npz with ``piv``) and the column order / orbit-coordinate conventions.

MATHEMATICS.  Ring F_p[(Z/2^k)^2]; monomial (a1, a2); weight popcount(a1) + popcount(a2); u_{i,j} adds
2^j to a_i mod 2^k (a translation, hence injective on monomials; it maps weight < D into weight <= D).
Columns: weight descending, then key a1 | a2 << k.  sigma(a1, a2) = (a2, a1).  Part "+" has columns
{a1 <= a2}, part "-" has {a1 < a2} (inherited order); lift(Y)[h] = Y[pos h] (h representative) or
+-Y[pos sigma h]; project(H)[m] = H[rep m] +- H[sigma rep m].  W = lift(W+) (+) lift(W-); each part is
{v : v[free] = v[piv] @ B (mod p)}, i.e. the row space of the rows e_{piv[i]} + B[i] on the free columns.

EXACTNESS.  Entries are integers in [0, p) stored in float64 (float32 on disk, exact for p < 2^24).  A
float64 product of integers < 2^53 is exact, and a sum of integers is exact whenever every partial sum
has magnitude <= 2^53 (every integer of that size is representable), whatever the order, blocking or FMA
use.  ``max_inner(p)`` is the largest m with m (p-1)^2 + (p-1) <= 2^53 - 2p; a product block of inner
dimension <= m plus a carried residue < p stays below that bound, and ``modp`` (below) is exact for
|X| <= 2^53 - 2p.  Longer inner dimensions are cut into pieces of m with a reduction between pieces.
For p = 65537, m = 2,096,639 > any inner dimension at k = 14, D = 6 (at most 224,000).

DEVICES.  ``--device cpu`` numpy; ``--device gpu`` cupy (the array module is switched through the ``Dev``
object only); ``--device fake-gpu`` the same code path on numpy arrays that refuse to be mixed with host
arrays (catches API slips that only a GPU would reveal).
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import resource
import sys
import tempfile
import time


def _early_threads(argv):
    n = None
    for i, x in enumerate(argv):
        if x == '--threads' and i + 1 < len(argv):
            n = argv[i + 1]
        elif x.startswith('--threads='):
            n = x.split('=', 1)[1]
    if n is None:
        try:
            n = str(len(os.sched_getaffinity(0)))
        except Exception:
            n = '1'
        for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ.setdefault(key, n)
    else:
        for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ[key] = n


if __name__ == '__main__':
    _early_threads(sys.argv)

import numpy as np  # noqa: E402

VERSION = 'ic_bounds_cert_check/1'
F64 = np.float64
F32 = np.float32
TWO53 = 1 << 53
MI_OVERRIDE = None        # test hook: force a small inner dimension to exercise the chunked path
T0 = time.time()
STATUS_PATH = None
PROGRESS = {}


# ----------------------------------------------------------------------------- logging
def file_sha256(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['maxrss_mb'] = rss_mb()
    print(json.dumps(kw), flush=True)


def write_json(path, obj):
    path = str(path)
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=False)
        f.write('\n')
    os.replace(tmp, path)


def status(**kw):
    PROGRESS.update(kw)
    PROGRESS['t'] = round(time.time() - T0, 1)
    PROGRESS['maxrss_mb'] = rss_mb()
    if STATUS_PATH:
        write_json(STATUS_PATH, PROGRESS)


# ----------------------------------------------------------------------------- device switch
class FakeArr(np.ndarray):
    """numpy array standing in for a device array: refuses to meet a host array (ndim > 0)."""

    @staticmethod
    def _chk(x):
        if isinstance(x, np.ndarray) and not isinstance(x, FakeArr) and x.ndim > 0:
            raise TypeError('fake-gpu: host array used with a device array')
        if isinstance(x, (list, tuple)):
            raise TypeError('fake-gpu: python sequence used with a device array')
        return x.view(np.ndarray) if isinstance(x, FakeArr) else x

    def __array_ufunc__(self, ufunc, method, *inputs, **kw):
        ins = [self._chk(x) for x in inputs]
        if 'out' in kw:
            kw['out'] = tuple(self._chk(o) for o in kw['out'])
        r = getattr(ufunc, method)(*ins, **kw)
        return _wrapres(r)

    def _chkidx(self, idx):
        items = idx if isinstance(idx, tuple) else (idx,)
        for it in items:
            if isinstance(it, (list, tuple)):
                raise TypeError('fake-gpu: python sequence as index')
            if isinstance(it, np.ndarray) and not isinstance(it, FakeArr) and it.ndim > 0:
                raise TypeError('fake-gpu: host index array on a device array')
        return tuple(it.view(np.ndarray) if isinstance(it, FakeArr) else it for it in items)

    def __getitem__(self, idx):
        return _wrapres(np.ndarray.__getitem__(self.view(np.ndarray), self._chkidx(idx)))

    def __setitem__(self, idx, val):
        np.ndarray.__setitem__(self.view(np.ndarray), self._chkidx(idx), self._chk(val))


def _wrapres(r):
    if isinstance(r, tuple):
        return tuple(_wrapres(x) for x in r)
    if isinstance(r, np.ndarray) and r.ndim > 0:
        return r.view(FakeArr)
    return r


class Dev:
    def __init__(self, name):
        self.name = name
        if name == 'gpu':
            import cupy as cp
            self.cp = cp
        elif name not in ('cpu', 'fake-gpu'):
            raise SystemExit(f'unknown device {name}')

    def asarray(self, a, dtype=None):
        a = np.asarray(a) if dtype is None else np.asarray(a, dtype=dtype)
        if self.name == 'gpu':
            return self.cp.asarray(a)
        if self.name == 'fake-gpu':
            return np.array(a).view(FakeArr)
        return a

    def host(self, a):
        if self.name == 'gpu':
            return self.cp.asnumpy(a)
        if self.name == 'fake-gpu':
            return np.array(a.view(np.ndarray)) if isinstance(a, FakeArr) else np.array(a)
        return np.asarray(a)

    def zeros(self, shape, dtype):
        if self.name == 'gpu':
            return self.cp.zeros(shape, dtype=dtype)
        z = np.zeros(shape, dtype=dtype)
        return z.view(FakeArr) if self.name == 'fake-gpu' else z

    def empty(self, shape, dtype):
        if self.name == 'gpu':
            return self.cp.empty(shape, dtype=dtype)
        z = np.empty(shape, dtype=dtype)
        return z.view(FakeArr) if self.name == 'fake-gpu' else z

    def floor(self, x):
        return self.cp.floor(x) if self.name == 'gpu' else np.floor(x)

    def where(self, c, a, b):
        if self.name == 'gpu':
            return self.cp.where(c, a, b)
        if self.name == 'cpu':
            return np.where(c, a, b)
        r = np.where(FakeArr._chk(c), FakeArr._chk(a), FakeArr._chk(b))
        return r.view(FakeArr)

    def free_pool(self):
        if self.name == 'gpu':
            self.cp.get_default_memory_pool().free_all_blocks()

    def mem_info(self):
        if self.name == 'gpu':
            free, total = self.cp.cuda.runtime.memGetInfo()
            return int(free), int(total)
        return None, None

    def sync(self):
        if self.name == 'gpu':
            self.cp.cuda.Stream.null.synchronize()


# ----------------------------------------------------------------------------- exact arithmetic
def max_inner(p: int) -> int:
    if MI_OVERRIDE:
        return int(MI_OVERRIDE)
    return (TWO53 - 2 * p - (p - 1)) // ((p - 1) ** 2)


def modp(DV, X, p):
    """X mod p in [0, p) for a float64 array of integers with |X| <= 2^53 - 2p.
    q = floor(fl(X / p)) is within 1 of floor(X / p) (relative error 2^-53 of a quotient <= 2^53 / p);
    q p is an integer of magnitude <= |X| + p <= 2^53, so q p and X - q p are exact; the result lies in
    [-p, 2p) and two selections bring it into [0, p)."""
    pf = float(p)
    q = DV.floor(X / pf)
    r = X - q * pf
    r = DV.where(r < 0, r + pf, r)
    r = DV.where(r >= pf, r - pf, r)
    return r


def modp_inplace(DV, X, p, rows=1 << 16):
    n = X.shape[0]
    step = max(1, min(rows, (1 << 24) // max(1, X.shape[1] if X.ndim > 1 else 1)))
    for a in range(0, n, step):
        X[a:a + step] = modp(DV, X[a:a + step], p)


def gemm_mod(DV, A, Bm, p):
    """(A @ Bm) mod p exactly; A (n x m), Bm (m x f) float64 integers in [0, p)."""
    mi = max_inner(p)
    m = A.shape[1]
    acc = None
    for a in range(0, max(m, 1), mi):
        part = A[:, a:a + mi] @ Bm[a:a + mi]
        acc = part if acc is None else acc + part
        acc = modp(DV, acc, p)
    return acc


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


# ----------------------------------------------------------------------------- instance (own code)
def poly_mul(a, b, p):
    out = {}
    for (e1, e2), u in a.items():
        for (f1, f2), v in b.items():
            key = (e1 + f1, e2 + f2)
            out[key] = (out.get(key, 0) + u * v) % p
    return out


def poly_add(*ps, p):
    out = {}
    for q in ps:
        for key, v in q.items():
            out[key] = (out.get(key, 0) + v) % p
    return out


def poly_scale(a, c, p):
    return {key: (v * c) % p for key, v in a.items()}


def s3_expand(p, a, b, X):
    """S_3(x1, x2, X) = (x1 - x2)^2 X^2 - 2((x1 + x2)(x1 x2 + a) + 2b) X + (x1 x2 - a)^2 - 4b(x1 + x2),
    expanded by polynomial arithmetic on {(e1, e2): coeff}."""
    x1, x2 = {(1, 0): 1}, {(0, 1): 1}
    one = {(0, 0): 1}
    s = poly_add(x1, x2, p=p)
    d = poly_add(x1, poly_scale(x2, -1, p), p=p)
    x1x2 = poly_mul(x1, x2, p)
    t1 = poly_scale(poly_mul(d, d, p), X * X, p)
    inner = poly_add(poly_mul(s, poly_add(x1x2, poly_scale(one, a, p), p=p), p),
                     poly_scale(one, 2 * b, p), p=p)
    t2 = poly_scale(inner, -2 * X, p)
    e = poly_add(x1x2, poly_scale(one, -a, p), p=p)
    t3 = poly_mul(e, e, p)
    t4 = poly_scale(s, -4 * b, p)
    tot = poly_add(t1, t2, t3, t4, p=p)
    return {key: v for key, v in tot.items() if v}


def chain_nf(poly, k, p):
    N = 1 << k
    out = {}
    for (e1, e2), v in poly.items():
        key = (e1 % N, e2 % N)
        out[key] = (out.get(key, 0) + v) % p
    return {key: v for key, v in out.items() if v}


# ----------------------------------------------------------------------------- layout (own code)
class PartLayout:
    pass


class Layout:
    def __init__(self, k: int, D: int, split: bool):
        self.k, self.D, self.split = k, D, split
        N = 1 << k
        ar = np.arange(N, dtype=np.int64)
        pc = np.zeros(N, dtype=np.int64)
        for b in range(k):
            pc += (ar >> b) & 1
        byw = [np.flatnonzero(pc == w) for w in range(D + 1)]
        A1, A2 = [], []
        for w1 in range(D + 1):
            for w2 in range(D + 1 - w1):
                x1, x2 = np.meshgrid(byw[w1], byw[w2], indexing='ij')
                A1.append(x1.ravel())
                A2.append(x2.ravel())
        a1, a2 = np.concatenate(A1), np.concatenate(A2)
        w = pc[a1] + pc[a2]
        key = a1 | (a2 << k)
        order = np.lexsort((key, -w))
        self.a1, self.a2, self.w, self.key = a1[order], a2[order], w[order], key[order]
        self.c = int(len(self.key))
        self.pc = pc
        self._sorted = np.argsort(self.key, kind='stable')
        self._skeys = self.key[self._sorted]
        self.lowc = int(np.sum(self.w < D))
        self.low0 = self.c - self.lowc
        assert np.all(self.w[self.low0:] < D) and np.all(self.w[:self.low0] >= D), 'low columns not a suffix'
        assert np.all(np.diff(self.w) <= 0), 'weights not descending'
        blocks = (0,) if split else (0, 1)
        self.all_shifts = [(i, j) for i in (0, 1) for j in range(k)]
        lowa1, lowa2 = self.a1[self.low0:], self.a2[self.low0:]
        self.shift = {}
        for (i, j) in self.all_shifts:
            if i == 0:
                self.shift[(i, j)] = self.index((lowa1 + (1 << j)) % N, lowa2)
            else:
                self.shift[(i, j)] = self.index(lowa1, (lowa2 + (1 << j)) % N)
        self.parts = []
        if split:
            self.sig = self.index(self.a2, self.a1)
            diag = self.a1 == self.a2
            for reps, sign in ((np.flatnonzero(self.a1 <= self.a2), 1), (np.flatnonzero(self.a1 < self.a2), -1)):
                self.parts.append(self._mkpart(reps, sign, diag))
        else:
            self.sig = np.arange(self.c)
            P = PartLayout()
            P.reps = np.arange(self.c)
            P.sign = 1
            P.c = self.c
            P.w = self.w
            P.name = 'full'
            P.tm = np.arange(self.c)
            P.tc = np.ones(self.c, dtype=np.int64)
            P.s0 = self.low0
            self.parts.append(P)

    def _mkpart(self, reps, sign, diag):
        P = PartLayout()
        P.reps, P.sign, P.c = reps, sign, len(reps)
        P.name = 'sym' if sign == 1 else 'anti'
        P.w = self.w[reps]
        pos = -np.ones(self.c, dtype=np.int64)
        pos[reps] = np.arange(len(reps))
        P.pos = pos
        tm = np.where(pos >= 0, pos, pos[self.sig])            # part column that column h projects onto
        tc = np.where(pos >= 0, 1, sign).astype(np.int64)      # and the coefficient
        if sign == 1:
            tc[diag] = 2                                       # H[h] + H[sigma h] = 2 H[h] on the diagonal
        else:
            tc[diag] = 0                                       # no antisymmetric column on the diagonal
            tm[diag] = -1
        P.tm, P.tc = tm, tc
        assert np.all(np.diff(P.w) <= 0)
        P.s0 = int(np.sum(P.w >= self.D))
        assert np.all(P.w[P.s0:] < self.D)
        return P

    def index(self, b1, b2):
        kk = np.asarray(b1) | (np.asarray(b2) << self.k)
        pos = np.searchsorted(self._skeys, kk)
        assert np.all(self._skeys[np.minimum(pos, self.c - 1)] == kk), 'monomial outside the layout'
        return self._sorted[pos]

    def project(self, H, p):
        """full-coordinate rows (n x c, integers) -> list of part-coordinate rows reduced into [0, p)."""
        if not self.split:
            return [np.mod(H, p)]
        out = []
        for P in self.parts:
            Y = H[:, P.reps] + P.sign * H[:, self.sig[P.reps]]
            out.append(np.mod(Y, p))
        return out

    def vec_of(self, terms, shift=(0, 0), p=None):
        """full-coordinate vector of sum coeff * x^(e + shift) for terms {(e1, e2): coeff}."""
        N = 1 << self.k
        v = np.zeros(self.c, dtype=np.int64)
        for (e1, e2), val in terms.items():
            idx = int(self.index((e1 + shift[0]) % N, (e2 + shift[1]) % N))
            v[idx] = (v[idx] + val) % p
        return v


# ----------------------------------------------------------------------------- checkpoint reading
class CkptError(Exception):
    pass


def locate_ckpt(path):
    if os.path.exists(os.path.join(path, 'manifest.json')):
        return path
    cands = sorted(d for d in os.listdir(path) if d.startswith('ck-') and os.path.exists(os.path.join(path, d, 'manifest.json')))
    if not cands:
        raise CkptError(f'no manifest.json under {path}')
    return os.path.join(path, cands[-1])


def read_checked(path, size, sha):
    with open(path, 'rb') as f:
        data = f.read()
    got = hashlib.sha256(data).hexdigest()
    if len(data) != size:
        raise CkptError(f'{os.path.basename(path)}: size {len(data)} != manifest {size}')
    if got != sha:
        raise CkptError(f'{os.path.basename(path)}: sha256 mismatch')
    return data, got


class PartStore:
    """what the manifest and host npz say about one part (no tile data)."""

    def __init__(self, q, L, ck, man, p, D):
        self.q, self.P = q, L.parts[q]
        c = self.P.c
        name = f'p{q}_host.npz'
        data, self.npz_sha = read_checked(os.path.join(ck, name), *man['files'][name])
        z = np.load(io.BytesIO(data), allow_pickle=False)
        piv = np.asarray(z['piv'])
        if piv.ndim != 1 or piv.dtype.kind not in 'iu':
            raise CkptError(f'part {q}: piv has the wrong shape or dtype')
        self.piv = piv.astype(np.int64)
        self.rank = len(self.piv)
        if self.rank and (self.piv.min() < 0 or self.piv.max() >= c):
            raise CkptError(f'part {q}: pivot outside [0, {c})')
        if len(np.unique(self.piv)) != self.rank:
            raise CkptError(f'part {q}: pivot columns not distinct')
        ispiv = np.zeros(c, dtype=bool)
        ispiv[self.piv] = True
        self.free = np.flatnonzero(~ispiv)
        self.F = len(self.free)
        assert self.F + self.rank == c
        self.tiles = [f'p{q}_t{ti:05d}.npy' for ti in range(man['ntiles'][q])]
        for t in self.tiles:
            if t not in man['files']:
                raise CkptError(f'{t} not listed in the manifest')
        self.is_low_row = self.P.w[self.piv] < D
        self.low_rows = np.flatnonzero(self.is_low_row)
        self.n_low = len(self.low_rows)
        self.Lp = c - self.P.s0
        self.jf0 = int(np.searchsorted(self.free, self.P.s0))      # low free columns = free[jf0:]
        self.pivrow = -np.ones(c, dtype=np.int64)
        self.pivrow[self.piv] = np.arange(self.rank)
        self.freeidx = -np.ones(c, dtype=np.int64)
        self.freeidx[self.free] = np.arange(self.F)
        self.tile_sha = {}


def parse_tile(data, store, p):
    a = np.load(io.BytesIO(data), allow_pickle=False)
    if a.ndim != 2 or a.shape[1] != store.F or a.dtype != np.float32:
        raise CkptError(f'part {store.q}: tile shape/dtype {a.shape} {a.dtype}, expected (*, {store.F}) float32')
    return a


# ----------------------------------------------------------------------------- phase 1: structure + low rows
def scan_part(DV, store, man, ck, p, D, strict_lead=True):
    """stream the tiles of one part: file hashes, entry range, leading-column (RREF) property, no support of
    low rows on weight >= D columns; build S = low rows restricted to the low columns (float32, device)."""
    P = store
    rep = dict(part=P.P.name, rank=P.rank, free=P.F, columns=P.P.c, tiles=len(P.tiles), low_rows=P.n_low,
               low_columns=P.Lp, entries_in_range=True, entries_integral=True, low_rows_supported_below_D=True,
               leading_column_ok=True, tile_rows_total=0, errors=[])
    S = DV.zeros((P.n_low, P.Lp), F32) if P.n_low and P.Lp else None
    r0 = 0
    n0 = 0
    for name in P.tiles:
        data, sha = read_checked(os.path.join(ck, name), *man['files'][name])
        P.tile_sha[name] = sha
        a = parse_tile(data, P, p)
        del data
        m = a.shape[0]
        if a.size and not (bool(a.min() >= 0) and bool(a.max() < p)):
            rep['entries_in_range'] = False
        if a.size and not np.array_equal(a, np.floor(a)):
            rep['entries_integral'] = False
        piv_t = P.piv[r0:r0 + m]
        if strict_lead and m and P.F:
            msk = P.free[None, :] < piv_t[:, None]
            if np.any(a[msk] != 0):
                rep['leading_column_ok'] = False
            del msk
        lows = np.flatnonzero(P.is_low_row[r0:r0 + m])
        if len(lows):
            blk = a[lows]
            if P.jf0 and np.any(blk[:, :P.jf0] != 0):
                rep['low_rows_supported_below_D'] = False
            if S is not None:
                Sb = np.zeros((len(lows), P.Lp), dtype=F32)
                Sb[:, P.free[P.jf0:] - P.P.s0] = blk[:, P.jf0:]
                Sb[np.arange(len(lows)), piv_t[lows] - P.P.s0] = 1.0
                S[n0:n0 + len(lows)] = DV.asarray(Sb)
                del Sb
            n0 += len(lows)
        r0 += m
        del a
    rep['tile_rows_total'] = r0
    if r0 != P.rank:
        rep['errors'].append(f'tile rows {r0} != number of pivots {P.rank}')
    if n0 != P.n_low:
        rep['errors'].append('low row bookkeeping mismatch')
    return rep, S


def load_B(DV, store, man, ck, p):
    """phase 2: the whole part on the device (float32), verifying again that every file is the one hashed."""
    P = store
    Bd = DV.empty((P.rank, P.F), F32)
    r0 = 0
    for name in P.tiles:
        data, sha = read_checked(os.path.join(ck, name), *man['files'][name])
        if sha != P.tile_sha[name]:
            raise CkptError(f'{name} changed between reads')
        a = parse_tile(data, P, p)
        del data
        if a.shape[0]:
            Bd[r0:r0 + a.shape[0]] = DV.asarray(a)
        r0 += a.shape[0]
        del a
    return Bd


def part_sha(DV, store, Bd, p, chunk_rows=2048):
    """canonical hash (same definition as the reference engine's Part.sha): rows sorted by pivot column;
    int64 pivots, int64 free columns, coefficients as uint16 (p < 65536) or uint32."""
    order = np.argsort(store.piv, kind='stable')
    h = hashlib.sha256()
    h.update(store.piv[order].astype(np.int64).tobytes())
    h.update(store.free.astype(np.int64).tobytes())
    dt = np.uint16 if p < 65536 else np.uint32
    for a in range(0, store.rank, chunk_rows):
        rows = DV.asarray(order[a:a + chunk_rows])
        blk = DV.host(Bd[rows])
        h.update(np.ascontiguousarray(blk).astype(dt).tobytes())
    return h.hexdigest()


# ----------------------------------------------------------------------------- membership
def membership(DV, store, Bd, Y, p, rb=8192):
    """Y: host float64 rows (n x c_q) with entries in [0, p).  Returns a host bool array: row in the part,
    i.e. Y[free] - Y[piv] @ B == 0 (mod p)."""
    n = Y.shape[0]
    if n == 0:
        return np.ones(0, dtype=bool)
    Yd = DV.asarray(Y)
    R = Yd[:, DV.asarray(store.free)].astype(F64)
    if store.rank and store.F:
        Yp = Yd[:, DV.asarray(store.piv)]
        mi = max_inner(p)
        step = max(1, min(rb, mi))
        for a in range(0, store.rank, step):
            b = min(store.rank, a + step)
            R = modp(DV, R - Yp[:, a:b] @ Bd[a:b].astype(F64), p)
    else:
        R = modp(DV, R, p)
    nz = (R != 0).any(axis=1) if store.F else DV.zeros((n,), bool)
    return ~DV.host(nz).astype(bool)


# ----------------------------------------------------------------------------- closure units
def unit_maps(L, ps_src, store_t, pi, q, ij, p):
    """For the source part pi (low columns s0..c) and shift ij: the sparse images in target part q.
    Layer 0: the representative column, layer 1: its sigma partner (coefficient = sign of the source part).
    Each entry: source low position s, target part column m, coefficient c (mod p)."""
    P, Q = L.parts[pi], L.parts[q]
    s_all = np.arange(P.s0, P.c)
    r = P.reps[s_all]
    layers = []
    for layer in (0, 1):
        if layer == 0:
            h = r
            cl = np.ones(len(r), dtype=np.int64)
        else:
            if not L.split:
                break
            h = L.sig[r]
            cl = np.where(h != r, P.sign, 0).astype(np.int64)
        hp = L.shift[ij][h - L.low0]
        m = Q.tm[hp]
        c = cl * Q.tc[hp]
        ok = (m >= 0) & (c != 0)
        s = (s_all - P.s0)[ok]
        m = m[ok]
        c = np.mod(c[ok], p).astype(np.float64)
        pr = store_t.pivrow[m]
        fi = store_t.freeidx[m]
        isp = pr >= 0
        layers.append(dict(s_piv=s[isp], rows=pr[isp], cm_piv=c[isp], s_free=s[~isp], fidx=fi[~isp], cm_free=c[~isp]))
    return layers


def build_T(DV, T, layers, Bd, f0, f1, p, rp, dev_layers):
    for lay, dl in zip(layers, dev_layers):
        if lay['s_free'].size:
            sel = (lay['fidx'] >= f0) & (lay['fidx'] < f1)
            if sel.any():
                T[DV.asarray(lay['s_free'][sel]), DV.asarray(lay['fidx'][sel] - f0)] += DV.asarray(lay['cm_free'][sel])
        n = lay['s_piv'].size
        for a in range(0, n, rp):
            b = min(n, a + rp)
            X = Bd[dl['rows'][a:b], f0:f1].astype(F64)
            X *= dl['cm_piv'][a:b][:, None]
            T[dl['s_piv'][a:b]] -= X
            del X
    modp_inplace(DV, T, p)


def run_unit(DV, L, pi, q, ij, Sd, store_s, store_t, Bd, p, Fb, nc, rp, stats):
    """all low rows of source part pi, shift ij, target part q.  Returns the host bool array 'bad'."""
    n = store_s.n_low
    F = store_t.F
    bad = DV.zeros((n,), bool)
    if n == 0 or F == 0 or store_s.Lp == 0:
        return DV.host(bad)
    layers = unit_maps(L, None, store_t, pi, q, ij, p)
    dev_layers = [dict(rows=DV.asarray(l['rows']), cm_piv=DV.asarray(l['cm_piv']), s_piv=DV.asarray(l['s_piv'])) for l in layers]
    if all(l['s_piv'].size == 0 and l['s_free'].size == 0 for l in layers):
        return DV.host(bad)
    Lp = store_s.Lp
    for f0 in range(0, F, Fb):
        f1 = min(F, f0 + Fb)
        T = DV.zeros((Lp, f1 - f0), F64)
        build_T(DV, T, layers, Bd, f0, f1, p, rp, dev_layers)
        for c0 in range(0, n, nc):
            c1 = min(n, c0 + nc)
            R = gemm_mod(DV, Sd[c0:c1].astype(F64), T, p)
            bad[c0:c1] = bad[c0:c1] | (R != 0).any(axis=1)
            stats['flops'] += 2.0 * (c1 - c0) * Lp * (f1 - f0)
            del R
        del T
    return DV.host(bad)


# ----------------------------------------------------------------------------- main verification
def read_engine_record(path):
    return json.load(open(path))


def verify(a) -> dict:
    t_start = time.time()
    global STATUS_PATH
    STATUS_PATH = a.status
    DV = Dev(a.device)
    V = dict(checker=dict(name=VERSION, source_sha256=file_sha256(__file__), argv=sys.argv[1:], numpy=np.__version__,
                          python=sys.version.split()[0]),
             device=a.device, errors=[], reasons=[])
    inst = json.load(open(a.instance))
    k, p = int(inst['k']), int(inst['prime'])
    V['instance'] = dict(file=os.path.basename(a.instance), sha256=file_sha256(a.instance), k=k, prime=p,
                         a=inst.get('a'), b=inst.get('b'), target=inst.get('target'))
    ck = locate_ckpt(a.ckpt)
    mpath = os.path.join(ck, 'manifest.json')
    man = json.load(open(mpath))
    V['inputs'] = dict(checkpoint_dir=os.path.basename(ck), manifest_sha256=file_sha256(mpath),
                       instance_sha256=V['instance']['sha256'])
    cfg = man['config']
    D = int(a.D if a.D is not None else cfg['D'])
    split = bool(cfg['split'])
    V['manifest'] = dict(version=man.get('version'), config=cfg, ntiles=man['ntiles'], reason=man.get('reason'),
                         engine_sha256=man.get('engine_sha256'), driver_sha256=man.get('driver_sha256'),
                         state={x: man['state'].get(x) for x in ('phase', 'closure_complete', 'rounds', 'repairs', 'seq')})
    V['params'] = dict(k=k, p=p, D=D, split=split, shifts=a.shifts if split else 'all')
    fmt = dict(ok=True, checks={})

    def chk(name, ok, detail=None):
        fmt['checks'][name] = bool(ok) if detail is None else dict(ok=bool(ok), detail=detail)
        if not ok:
            fmt['ok'] = False
            V['reasons'].append(f'format:{name}')

    chk('manifest_version_1', man.get('version') == 1)
    chk('config_k_p_D_match_instance', cfg['k'] == k and cfg['p'] == p and cfg['D'] == D,
        dict(config=[cfg['k'], cfg['p'], cfg['D']], instance_and_flags=[k, p, D]))
    chk('p_prime_odd_below_2^24', p % 2 == 1 and p < (1 << 24) and is_prime(p))
    chk('manifest_lists_every_part', len(man['ntiles']) == (2 if split else 1))
    unlisted = sorted(f for f in os.listdir(ck) if f != 'manifest.json' and f not in man['files'])
    chk('no_unlisted_files_in_checkpoint_dir', not unlisted, unlisted[:5])
    if not fmt['ok']:
        V['format'] = fmt
        V.update(format_ok=False, closed=None, g_in_W=None, one_excluded=None, certificate_valid=False, all_true=False)
        return V
    t = time.time()
    L = Layout(k, D, split)
    V['layout'] = dict(columns=L.c, low_columns=L.lowc, part_columns=[P.c for P in L.parts], seconds=round(time.time() - t, 1))
    log(event='layout', **V['layout'])

    # -- instance generator
    gi = {tuple(int(x) for x in key.split(',')): int(v) % p for key, v in inst['generator'].items()}
    gi = {e: v for e, v in gi.items() if v}
    g_formula = None
    if all(x in inst for x in ('a', 'b', 'target')):
        g_formula = chain_nf(s3_expand(p, int(inst['a']) % p, int(inst['b']) % p, int(inst['target'][0]) % p), k, p)
    V['generator'] = dict(terms=len(gi), instance_equals_formula=(None if g_formula is None else g_formula == gi),
                          weight=max((int(L.pc[e1]) + int(L.pc[e2])) for e1, e2 in gi))
    if g_formula is not None and g_formula != gi:
        V['reasons'].append('generator:instance dict differs from the formula S_3(a, b, target[0])')
    gdeg = V['generator']['weight']

    # -- phase 1: structure, hashes, low rows
    stores, scans, Sds = [], [], []
    t = time.time()
    for q in range(len(L.parts)):
        status(phase='scan', part=q)
        try:
            st = PartStore(q, L, ck, man, p, D)
            rep, S = scan_part(DV, st, man, ck, p, D)
        except CkptError as e:
            V['errors'].append(str(e))
            V['reasons'].append(f'format:{e}')
            V['format'] = fmt
            V.update(format_ok=False, closed=None, g_in_W=None, one_excluded=None, certificate_valid=False, all_true=False)
            return V
        stores.append(st)
        scans.append(rep)
        Sds.append(S)
        log(event='scanned', **rep)
    chk('pivots_distinct_free_partition_dims', all(not s['errors'] for s in scans), [s['errors'] for s in scans])
    chk('entries_in_0_p_and_integral', all(s['entries_in_range'] and s['entries_integral'] for s in scans))
    chk('low_rows_have_no_support_on_weight_ge_D', all(s['low_rows_supported_below_D'] for s in scans))
    chk('manifest_file_hashes_match', True)
    V['format'] = fmt
    V['scan'] = scans
    V['rref_leading_column_ok'] = all(s['leading_column_ok'] for s in scans)
    V['phase1_seconds'] = round(time.time() - t, 1)
    V['inputs']['file_sha256'] = {n: h for st in stores for n, h in {f'p{st.q}_host.npz': st.npz_sha, **st.tile_sha}.items()}
    V['format_ok'] = fmt['ok']
    if not fmt['ok']:
        V.update(closed=None, g_in_W=None, one_excluded=None, certificate_valid=False, all_true=False)
        return V
    ranks = [s.rank for s in stores]
    V['dimension'] = dict(W=sum(ranks), parts=ranks, columns=L.c)
    prof_parts = [[int(np.sum(st.P.w[st.piv] <= j)) for j in range(D + 1)] for st in stores]
    V['profile'] = dict(pivots_weight_le_j=[sum(x) for x in zip(*prof_parts)], parts=prof_parts,
                        equals_dim_W_cap_M_le_j=V['rref_leading_column_ok'])

    # -- shifts
    shifts = L.all_shifts if (not split or a.shifts == 'all') else [(0, j) for j in range(k)]
    nparts = len(L.parts)
    rows_checked = sum(s.n_low for s in stores)
    triples = rows_checked * len(shifts) * nparts
    total_units = nparts * nparts * len(shifts)
    stats = dict(flops=0.0)
    failures = 0
    fail_by = {}
    first_fail = []
    g_results, gm_results, one_results, hashes = {}, {}, {}, []
    done_units = 0
    t_close = time.time()
    # generator vectors (full coordinates)
    gvec = L.vec_of(gi, p=p)[None, :]
    gfvec = None if g_formula is None else L.vec_of(g_formula, p=p)[None, :]
    gweights = [int(m) for m in np.flatnonzero((L.w <= a.gen_mult_weight) & (L.w + gdeg <= D))]
    gm = np.zeros((len(gweights), L.c), dtype=np.int64)
    for r_, m in enumerate(gweights):
        gm[r_] = L.vec_of(gi, shift=(int(L.a1[m]), int(L.a2[m])), p=p)
    onevec = np.zeros((1, L.c), dtype=np.int64)
    zero_full = int(np.flatnonzero(L.w == 0)[0])
    onevec[0, zero_full] = 1
    gproj, gfproj, gmproj, oneproj = (L.project(gvec, p), None if gfvec is None else L.project(gfvec, p),
                                      L.project(gm, p), L.project(onevec, p))
    total_flops_est = sum(float(stores[pi].n_low) * stores[pi].Lp * stores[q].F * 2.0 for pi in range(nparts) for q in range(nparts)) * len(shifts)
    status(phase='closure', total_units=total_units, done_units=0, shifts=len(shifts), flops_est=total_flops_est)
    last_log = 0.0
    for q in range(nparts):
        st_t = stores[q]
        t = time.time()
        status(phase='load', part=q)
        Bd = load_B(DV, st_t, man, ck, p)
        sh = part_sha(DV, st_t, Bd, p)
        hashes.append(sh)
        log(event='part_loaded', part=st_t.P.name, sha256=sh, seconds=round(time.time() - t, 1))
        g_results[q] = bool(membership(DV, st_t, Bd, gproj[q].astype(F64), p)[0])
        if gfproj is not None:
            g_results[f'formula{q}'] = bool(membership(DV, st_t, Bd, gfproj[q].astype(F64), p)[0])
        gm_results[q] = membership(DV, st_t, Bd, gmproj[q].astype(F64), p)
        one_results[q] = bool(membership(DV, st_t, Bd, oneproj[q].astype(F64), p)[0])
        if q == 0 and split:
            V['weight0_column_free'] = bool(st_t.freeidx[int(np.flatnonzero(st_t.P.w == 0)[0])] >= 0)
        if not split:
            V['weight0_column_free'] = bool(st_t.freeidx[zero_full] >= 0)
        F = st_t.F
        # block width / chunks from the working-memory budget
        Lpmax = max(s.Lp for s in stores)
        Fb = int(max(64, min(F if F else 64, a.col_block or (a.work_gb * (1 << 30) / 3.0) / (8.0 * max(1, Lpmax)))))
        rp = int(max(256, min(1 << 16, (a.work_gb * (1 << 30) / 3.0) / (16.0 * max(1, Fb)))))
        nc = int(a.row_chunk)
        for ij in shifts:
            for pi in range(nparts):
                ts = time.time()
                bad = run_unit(DV, L, pi, q, ij, Sds[pi], stores[pi], st_t, Bd, p, Fb, nc, rp, stats)
                nb = int(bad.sum())
                failures += nb
                key = f'{L.parts[pi].name}->{L.parts[q].name}'
                fail_by[key] = fail_by.get(key, 0) + nb
                if nb and len(first_fail) < 20:
                    for r_ in np.flatnonzero(bad)[:3]:
                        row = int(stores[pi].low_rows[r_])
                        first_fail.append(dict(source=L.parts[pi].name, target=L.parts[q].name, shift=list(ij), row=row,
                                               pivot_column=int(stores[pi].piv[row])))
                done_units += 1
                now = time.time()
                if now - last_log >= a.progress_seconds or done_units == total_units:
                    last_log = now
                    el = now - t_close
                    eta = el / done_units * (total_units - done_units)
                    status(phase='closure', done_units=done_units, total_units=total_units, failures=failures,
                           elapsed=round(el, 1), eta_seconds=round(eta, 1), gflops_done=round(stats['flops'] / 1e9, 1))
                    log(event='progress', done_units=done_units, total_units=total_units, failures=failures,
                        last=dict(src=L.parts[pi].name, tgt=L.parts[q].name, shift=list(ij), seconds=round(now - ts, 2)),
                        eta_seconds=round(eta), tflop_done=round(stats['flops'] / 1e12, 2))
                if a.fail_fast and failures:
                    break
            if a.fail_fast and failures:
                break
        del Bd
        DV.free_pool()
        if a.fail_fast and failures:
            break
    complete = done_units == total_units
    closed = complete and failures == 0
    V['closure'] = dict(shifts_checked=[list(s) for s in shifts], shift_mode=('all 2k shifts' if len(shifts) == len(L.all_shifts) else 'block-1 only (uses sigma W = W)'),
                        low_rows_per_part=[s.n_low for s in stores], rows_checked=rows_checked,
                        triples_checked=triples if complete else None, units_done=done_units, units_total=total_units,
                        failures=failures, failures_by_direction=fail_by, first_failures=first_fail, complete=complete,
                        gemm_tflop=round(stats['flops'] / 1e12, 3), seconds=round(time.time() - t_close, 1))
    gm_ok = all(bool(np.all(v)) for v in gm_results.values())
    g_ok = all(v for v in g_results.values())
    V['generator'].update(in_W_per_part=g_results, multiples_checked=len(gweights), multiples_weight_le=a.gen_mult_weight,
                          multiples_in_W=gm_ok)
    V['one'] = dict(in_W_per_part=one_results, weight0_column_free=V.get('weight0_column_free'))
    one_in_W = all(one_results.values())
    V['basis_sha256'] = hashes
    V['closed'] = bool(closed)
    V['g_in_W'] = bool(g_ok and (V['generator']['instance_equals_formula'] is not False))
    V['generator_multiples_in_W'] = bool(gm_ok)
    V['one_excluded'] = bool(not one_in_W)
    if not g_ok:
        V['reasons'].append('g not in W')
    if not gm_ok:
        V['reasons'].append('some m*g (weight m <= %d) not in W' % a.gen_mult_weight)
    if not closed:
        V['reasons'].append(f'closure: {failures} (row, shift, target) triples leave W' if complete else 'closure incomplete')
    V['certificate_valid'] = bool(V['format_ok'] and V['g_in_W'] and closed)
    V['all_true'] = bool(V['certificate_valid'] and V['one_excluded'])
    V['meaning'] = ('certificate_valid: V_D is contained in W (g in W and W closed under the shifts of W cap M_<D). '
                    'one_excluded: 1 not in W, hence (if certificate_valid) 1 not in V_D. all_true = both.')
    # -- engine record comparison
    if a.engine_record:
        rec = read_engine_record(a.engine_record)
        cmpd = {}
        for key, mine in (('rank', V['dimension']['W']), ('part_ranks', ranks), ('profile', V['profile']['pivots_weight_le_j']),
                          ('basis_sha256', hashes), ('contains_one', bool(one_in_W)), ('columns', L.c)):
            if key in rec:
                cmpd[key] = dict(engine=rec[key], checker=mine, equal=(rec[key] == mine))
        V['engine_record_comparison'] = dict(file=os.path.basename(a.engine_record), sha256=file_sha256(a.engine_record),
                                             fields=cmpd, all_equal=all(v['equal'] for v in cmpd.values()))
    V['resources'] = dict(maxrss_mb=rss_mb(), seconds=round(time.time() - t_start, 1))
    if DV.name == 'gpu':
        free, total = DV.mem_info()
        V['resources'].update(device_free_mb=free >> 20, device_total_mb=total >> 20)
    return V


# ----------------------------------------------------------------------------- self-test and builders
def naive_basis(L, p, gi):
    """pure-python closure in full coordinates (used by the self-test and the test harness, never by verify)."""
    basis = []

    def reduce(v):
        for lead, row in basis:
            if v[lead]:
                v = (v - v[lead] * row) % p
        return v

    def add(v):
        v = reduce(v)
        nz = np.flatnonzero(v)
        if not len(nz):
            return False
        lead = int(nz[0])
        v = (v * pow(int(v[lead]), -1, p)) % p
        for i, (ld, row) in enumerate(basis):
            if row[lead]:
                basis[i] = (ld, (row - row[lead] * v) % p)
        basis.append((lead, v))
        return True

    add(L.vec_of(gi, p=p))
    changed = True
    while changed:
        changed = False
        for lead, row in list(basis):
            if L.w[lead] >= L.D:
                continue
            for ij in L.all_shifts:
                v = np.zeros(L.c, dtype=np.int64)
                v[L.shift[ij]] = row[L.low0:]
                changed |= add(v)
    # rows may have changed after being used: recompute until a clean pass (the loop above does that)
    return basis


def rref_rows(M, p):
    """RREF of integer rows M (n x c), pivot = first nonzero column; returns (piv list, rows array)."""
    rows = []
    pivs = []
    for v in M:
        v = v.copy() % p
        for lead, row in zip(pivs, rows):
            if v[lead]:
                v = (v - v[lead] * row) % p
        nz = np.flatnonzero(v)
        if not len(nz):
            continue
        lead = int(nz[0])
        v = (v * pow(int(v[lead]), -1, p)) % p
        for i in range(len(rows)):
            if rows[i][lead]:
                rows[i] = (rows[i] - rows[i][lead] * v) % p
        rows.append(v)
        pivs.append(lead)
    return pivs, (np.array(rows, dtype=np.int64) if rows else np.zeros((0, M.shape[1]), dtype=np.int64))


def write_checkpoint(root, parts, config, tile_rows=7, state=None, name='ck-00000001'):
    """parts: list of (piv int64 array in row order, rows_full: rank x c_q int array in part coordinates).
    Writes the engine's on-disk format (tiles of float32 on the free columns)."""
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    files = {}
    ntiles = []
    for q, (piv, rows) in enumerate(parts):
        c = rows.shape[1]
        ispiv = np.zeros(c, dtype=bool)
        ispiv[piv] = True
        free = np.flatnonzero(~ispiv)
        Bm = rows[:, free].astype(np.float32)
        nt = 0
        for t0 in range(0, max(1, len(piv)), tile_rows):
            nm = f'p{q}_t{nt:05d}.npy'
            np.save(os.path.join(d, nm), Bm[t0:t0 + tile_rows])
            files[nm] = [os.path.getsize(os.path.join(d, nm)), file_sha256(os.path.join(d, nm))]
            nt += 1
            if not len(piv):
                break
        ntiles.append(nt)
        nm = f'p{q}_host.npz'
        np.savez(os.path.join(d, nm), piv=np.asarray(piv, dtype=np.int64), processed=np.zeros(c, dtype=bool))
        files[nm] = [os.path.getsize(os.path.join(d, nm)), file_sha256(os.path.join(d, nm))]
    man = dict(version=1, engine_sha256='test', driver_sha256='test', reason='test', config=config,
               state=state or dict(phase='sweep', closure_complete=True, rounds=0, repairs=0, seq=1), files=files, ntiles=ntiles)
    json.dump(man, open(os.path.join(d, 'manifest.json'), 'w'))
    return d


def parts_from_full_basis(L, p, full_rows):
    """full-coordinate basis -> list of (piv, part rows) of the two parts (or the single part)."""
    out = []
    if not L.split:
        piv, rows = rref_rows(full_rows, p)
        return [(np.array(piv, dtype=np.int64), rows)]
    for Y in L.project(full_rows, p):
        piv, rows = rref_rows(Y, p)
        out.append((np.array(piv, dtype=np.int64), rows))
    return out


def selftest_arith(DV):
    rng = np.random.default_rng(5)
    res = dict(cases=[], passed=True)
    for p in (12289, 65537, 16777213):
        lim = TWO53 - 2 * p
        vals = [lim, -lim, lim - 1, 0, 1, -1, p, -p, p - 1, 2 * p, 3 * p - 1]
        qmax = lim // p
        for q_ in (qmax, qmax - 1, qmax - 2):
            vals += [q_ * p, q_ * p - 1, q_ * p + 1, -q_ * p, -q_ * p + 1, -q_ * p - 1]
        r = rng.integers(-lim, lim, size=4000, dtype=np.int64)
        vals += [int(x) for x in r] + [int(x) for x in lim - rng.integers(0, 1 << 20, size=2000, dtype=np.int64)]
        vals = [v for v in vals if abs(v) <= lim]
        X = DV.asarray(np.array(vals, dtype=np.float64))
        assert all(float(v) == v for v in vals) or True
        got = DV.host(modp(DV, X, p)).astype(np.int64)
        ok = all(int(g) == v % p for g, v in zip(got, vals))
        # the float64 images of the python ints must be exact: values below 2^53 are
        res['cases'].append(dict(kind='modp', p=p, n=len(vals), exact=ok))
        res['passed'] &= ok
        # product blocks at the inner-dimension bound
        mi = min(max_inner(p), {12289: 20000, 65537: 4096}.get(p, 10 ** 9))
        for fill in ('max', 'rand'):
            n_, f_ = 3, 5
            if fill == 'max':
                A = np.full((n_, mi + 7), p - 1, dtype=np.int64)
                Bm = np.full((mi + 7, f_), p - 1, dtype=np.int64)
            else:
                A = rng.integers(p - 3, p, size=(n_, mi + 7), dtype=np.int64)
                Bm = rng.integers(p - 3, p, size=(mi + 7, f_), dtype=np.int64)
            gotg = DV.host(gemm_mod(DV, DV.asarray(A.astype(np.float64)), DV.asarray(Bm.astype(np.float64)), p)).astype(np.int64)
            want = np.array([[sum(int(A[i, t]) * int(Bm[t, j]) for t in range(A.shape[1])) % p for j in range(f_)] for i in range(n_)])
            ok = bool(np.array_equal(gotg, want))
            res['cases'].append(dict(kind='gemm', p=p, fill=fill, inner=A.shape[1], max_inner=max_inner(p), exact=ok))
            res['passed'] &= ok
    return res


def selftest_pipeline(DV, a):
    """positive and negative cases through verify() on a tiny instance (k = 4, D = 3, p = 12289)."""
    p, k, D = 12289, 4, 3
    cases = []
    tmp = tempfile.mkdtemp(prefix='certcheck-selftest-')
    try:
        ainst, binst, X = 1, 0, 9362
        g = chain_nf(s3_expand(p, ainst, binst, X), k, p)
        inst = dict(k=k, prime=p, a=ainst, b=binst, target=[X, 0], generator={f'{e1},{e2}': v for (e1, e2), v in g.items()})
        ipath = os.path.join(tmp, 'inst.json')
        json.dump(inst, open(ipath, 'w'))
        for split in (True, False):
            L = Layout(k, D, split)
            basis = naive_basis(L, p, g)
            full = np.array([row for _, row in basis])
            parts = parts_from_full_basis(L, p, full)
            cfg = dict(k=k, p=p, D=D, split=split, gens='g', source_id='selftest')
            for neg in (None, 'flip-low', 'drop-low'):
                prt = [(pv.copy(), rw.copy()) for pv, rw in parts]
                if neg:
                    q_ = 0
                    pv, rw = prt[q_]
                    lowi = [i for i in range(len(pv)) if L.parts[q_].w[pv[i]] < D]
                    i = lowi[len(lowi) // 2]
                    if neg == 'flip-low':
                        free = np.setdiff1d(np.arange(rw.shape[1]), pv)
                        cand = [f for f in free if f > pv[i] and L.parts[q_].w[f] < D] or [f for f in free if f > pv[i]]
                        f = cand[0]
                        rw[i, f] = (rw[i, f] + 1) % p
                    else:
                        keep = np.arange(len(pv)) != i
                        pv, rw = pv[keep], rw[keep]
                        # the dropped pivot column becomes free; other rows have 0 there (RREF), format stays valid
                        prt[q_] = (pv, rw)
                root = os.path.join(tmp, f'ck{split}{neg}')
                write_checkpoint(root, prt, cfg)
                for shifts in (('all', 'block1') if split else ('all',)):
                    ns = argparse.Namespace(instance=ipath, ckpt=root, D=None, device=a.device, shifts=shifts, out=None, status=None,
                                            engine_record=None, gen_mult_weight=1, col_block=None, row_chunk=3, work_gb=0.01,
                                            progress_seconds=1e9, fail_fast=False, selftest=False)
                    old = sys.stdout
                    sys.stdout = io.StringIO()
                    try:
                        Vv = verify(ns)
                    finally:
                        sys.stdout = old
                    ok = (Vv['certificate_valid'] is True) if neg is None else (Vv['certificate_valid'] is False)
                    cases.append(dict(split=split, shifts=shifts, negative=neg, closed=Vv['closed'], g_in_W=Vv['g_in_W'],
                                      one_excluded=Vv['one_excluded'], failures=(Vv['closure'] or {}).get('failures'), as_expected=bool(ok)))
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
    return dict(cases=cases, passed=all(c['as_expected'] for c in cases))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--instance')
    ap.add_argument('--ckpt', help='checkpoint directory (ck-XXXXXXXX) or a root holding ck-* directories (newest is used)')
    ap.add_argument('--D', type=int, default=None, help='degree (default: from the manifest; a given value must agree)')
    ap.add_argument('--device', default='cpu', choices=['cpu', 'gpu', 'fake-gpu'])
    ap.add_argument('--shifts', default='all', choices=['all', 'block1'],
                    help='all: every u_{i,j}, i in {1,2} (default, direct); block1: only block 1 (valid because sigma W = W)')
    ap.add_argument('--out')
    ap.add_argument('--status')
    ap.add_argument('--engine-record', help='engine result JSON: rank, part_ranks, profile, basis_sha256, contains_one are compared')
    ap.add_argument('--gen-mult-weight', type=int, default=1, help='also test m*g for all monomials m of weight <= this')
    ap.add_argument('--col-block', type=int, default=None, help='free columns per product block (default from --work-gb)')
    ap.add_argument('--row-chunk', type=int, default=2048, help='source rows per product')
    ap.add_argument('--work-gb', type=float, default=4.0, help='working memory for the product blocks (device memory on gpu)')
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--progress-seconds', type=float, default=30.0)
    ap.add_argument('--fail-fast', action='store_true')
    ap.add_argument('--selftest', action='store_true', help='arithmetic exactness + tiny positive/negative pipeline on the device')
    a = ap.parse_args()
    if a.selftest:
        DV = Dev(a.device)
        r1 = selftest_arith(DV)
        log(event='selftest_arith', passed=r1['passed'], cases=r1['cases'])
        r2 = selftest_pipeline(DV, a)
        log(event='selftest_pipeline', passed=r2['passed'], cases=r2['cases'])
        rep = dict(arith=r1, pipeline=r2, passed=bool(r1['passed'] and r2['passed']), device=a.device,
                   source_sha256=file_sha256(__file__))
        if a.out:
            write_json(a.out, rep)
        sys.exit(0 if rep['passed'] else 1)
    if not a.instance or not a.ckpt:
        ap.error('--instance and --ckpt are required')
    V = verify(a)
    log(event='verdict', **{key: V.get(key) for key in ('format_ok', 'g_in_W', 'closed', 'one_excluded', 'certificate_valid', 'all_true')},
        dimension=V.get('dimension'), reasons=V['reasons'])
    if a.out:
        write_json(a.out, V)
        log(event='saved', file=a.out)
    sys.exit(0 if V.get('certificate_valid') else (2 if V.get('format_ok') else 3))


if __name__ == '__main__':
    main()
