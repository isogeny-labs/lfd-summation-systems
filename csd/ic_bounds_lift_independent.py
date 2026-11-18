"""Independent degree lift V_D -> V_{D+1} n span(M_<=D) for two-summand chain systems.

Second implementation of the lift of experiments/ic-bounds-k14/LIFT.md.  It imports no repository
module (numpy only; cupy when the caller passes it as ``xp``).  From the first implementation
(csd/ic_bounds_lift.py) and the engine (csd/ic_bounds_big.py) it shares the SPECIFICATION only: the
on-disk checkpoint format of the engine (read here by ``read_checkpoint``), its column order and
orbit convention (needed to read that file), and the definition of the canonical basis hash.

Setting.  F = F_p, monomials (a1, a2) of two k-bit masks, weight = popcount(a1) + popcount(a2); the
shift u_{i,j} adds 2^j to a_i modulo 2^k (weight goes up by at most 1).  L = span(M_<=D).  V_n is the
least subspace of span(M_<=n) containing g and closed under f -> u f for f in V_n of weight < n.

What is different from the first implementation.  Everything after the file is read happens in
ORDINARY monomial coordinates: no sigma-parts, no orbit columns, no projection coefficients.

  * From the two stored parts (row = e_pivot + B[row] on the free orbit columns; a stored vector Y
    of the +/- part stands for sum_o Y[o] (e_m +/- e_{sigma m}), a self-paired orbit for Y[o] e_m)
    we form, for every monomial x whose orbit is a pivot of both parts (self-paired: of the + part),
        l_x = (lift b+ +/- lift b-) / 2 = e_x + phi(x),      phi(x) supported on N,
    where Pd is the set of those x and N = M_<=D \\ Pd.  The remaining stored rows ("single" rows,
    orbit a pivot of exactly one part) are supported on N and span V_N.  So
        V_D = V_P (+) V_N,   V_P = span{l_x : x in Pd},   L = V_P (+) F^N.
  * A weight-(D+1) monomial ("top") is *special* if it is u n for some n in N, else *regular*.
    T = F^Sp (+) F^N.  red : span(M_<=D+1) -> F^Reg (+) T replaces e_x (x in Pd) by -phi(x); it is
    onto with kernel V_P, and an element of span(M_<=D+1) with red in F^N lies in L.
  * Generators gamma(x, u) = red(u l_x) = red(e_{ux}) + phi(x) R_u, where R_u : F^N -> T is
    e_n -> red(e_{un}).  If u x is a regular top t, gamma has regular coordinate e_t with coefficient
    exactly 1 and no other regular coordinate (n -> u n is injective and u N misses regular tops).
    A regular top t has exactly D + 1 such generators, x = t / u for the D + 1 variables dividing t.

Algorithm (Lemma 1 of LIFT.md: W_0 = V_D, W_{i+1} = (W_i + sum_u u W_i) n L, and W_{i+1} = W_i
implies W_i = V_{D+1} n L).  Write W_i = V_P (+) Phi_i, Phi_i c F^N, Phi_0 = V_N.  Then
    red(W_i + sum_u u W_i) = span{gamma(x, u)} + Phi_i + sum_u Phi_i R_u,
and its intersection with T is (triangular elimination of the regular coordinates, each with
pivot 1)
    K_i = span{gamma(t/u, u) - gamma(t/u0, u0) : t regular, u | t}
        + span{gamma(x, u) : u x not a regular top} + Phi_i + sum_u Phi_i R_u,
so Phi_{i+1} = K_i n F^N.  K is kept in reduced echelon form with the Sp columns first, hence the
rows of K with pivot in an N column (the "falls") are a basis of K n F^N.  The code grows one space
K: (1) V_N; (2) the generator differences and direct generators, taken in chunks of a fixed random
order of the regular tops and reduced modulo the K at the start of the chunk; (3) after each chunk,
every fall row not yet shifted is multiplied by every R_u, until none is pending (a fall row is
shifted once, at its state then: afterwards it changes only by multiples of later fall rows, which
are shifted in their turn, so the span of the shifted states plus the pending rows is always the
whole fall space); (4) at the end ALL fall rows are shifted again and must reduce to zero, an
explicit check of Phi R_u c K that does not use the argument in (3).
Only x of weight D are shifted (Lemma 2 of LIFT.md needs V_D closed); the generators skipped for
that reason are tested by ``check_closed`` (random combinations must already lie in V_N), or
included with ``include_low``.

Soundness and completeness.  Every vector ever added to K is red of an element of V_{D+1} with zero
regular coordinates, PROVIDED the input space is contained in V_D (this cannot be tested here; see
the README).  So "1 in K" proves 1 in V_{D+1} whatever part of the schedule has run.  If K becomes
all of T ("saturated") then V_{D+1} n L = L and V_{D+1} = span(M_<=D+1), again by soundness alone, and
the run stops.  Otherwise the result is exact when every chunk has been processed and (4) passed.

Arithmetic: float64 holding integers; a product block has inner dimension <= max_inner(p), so every
partial sum is below 2^53 (exact in any summation order); modp corrects the rounded quotient.
The stored basis is kept as float32 (exact for p < 2^24) in the engine's orbit layout, on the host or
(``basis_on='device'``) on the device; rows are converted to ordinary coordinates when they are used.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import time

import numpy as np

VERSION = 'ic_bounds_lift_independent/1'
I64 = np.int64


class Refuse(Exception):
    """the input fails a check; no result is produced."""


class Stop(Exception):
    """deadline reached or 1 found with stop_on_one: the record holds what is proved so far."""


# ------------------------------------------------------------------------------- exact arithmetic

def max_inner(p: int) -> int:
    """largest m with m (p-1)^2 + (p-1) <= 2^53 - 2p."""
    return (2 ** 53 - 2 * p - (p - 1)) // ((p - 1) ** 2)


def modp(X, p: int, xp):
    """in place: integers |X| <= 2^53 - 2p (float64) -> X mod p in [0, p)."""
    q = xp.floor(X / p)
    q *= p
    X -= q                      # in [-p, 2p): the rounded quotient is within 1 of the true one
    X += (X < 0) * float(p)
    X -= (X >= p) * float(p)
    return X


def mm(A, B, p: int, xp, mi: int | None = None):
    """(A @ B) mod p; requested positive chunks are capped at the exactness bound."""
    safe_inner = max_inner(p)
    if mi is not None and (not isinstance(mi, int) or mi <= 0):
        raise Refuse('mi must be a positive integer')
    mi = safe_inner if mi is None else min(mi, safe_inner)
    n = A.shape[1]
    if A.shape[0] == 0 or n == 0 or B.shape[1] == 0:
        return xp.zeros((A.shape[0], B.shape[1]), dtype=xp.float64)
    if n <= mi:
        return modp(A @ B, p, xp)
    acc = xp.zeros((A.shape[0], B.shape[1]), dtype=xp.float64)
    for a in range(0, n, mi):
        acc += A[:, a:a + mi] @ B[a:a + mi]
        modp(acc, p, xp)
    return acc


def to_host(a) -> np.ndarray:
    return a.get() if hasattr(a, 'get') else np.asarray(a)


def neg(X, p: int):
    """-X mod p for X in [0, p), in place."""
    X *= -1.0
    X += (X < 0) * float(p)
    return X


def rref(R, p: int, xp, block: int = 64, mi: int | None = None):
    """Reduced row echelon form of the row space of R (entries in [0, p)); the pivot of a row is its
    first nonzero column.  Returns (L, lead): L[i] has leading column lead[i]."""
    f = R.shape[1]
    L = xp.zeros((0, f), dtype=xp.float64)
    lead: list = []
    for a in range(0, R.shape[0], block):
        S = R[a:a + block].copy()
        if lead:
            S -= mm(S[:, xp.asarray(lead)], L, p, xp, mi)
            modp(S, p, xp)
        S = S[xp.any(S != 0, axis=1)]
        new: list = []
        for i in range(S.shape[0]):
            v = S[i]
            if new:                                      # clear the pivots found in this block
                v -= mm(v[xp.asarray(new)][None, :], S[:len(new)], p, xp, mi)[0]
                modp(v, p, xp)
            c = int(xp.argmax(v != 0))
            if float(v[c]) == 0.0:
                continue
            v *= float(pow(int(v[c]), -1, p))
            modp(v, p, xp)
            if new:
                S[:len(new)] -= S[:len(new), c][:, None] * v[None, :]
                modp(S[:len(new)], p, xp)
            S[len(new)] = v
            new.append(c)
        if not new:
            continue
        Ln = S[:len(new)]
        if lead:
            L -= mm(L[:, xp.asarray(new)], Ln, p, xp, mi)
            modp(L, p, xp)
        L = xp.concatenate([L, Ln])
        lead += new
    return L, np.array(lead, dtype=I64)


class Echelon:
    """Subspace of F^n in reduced echelon form: row i = e_{piv[i]} + B[i] on the columns ``free``
    (sorted).  ``add`` never changes an array in place, so a shallow copy is a frozen snapshot."""

    def __init__(self, n: int, p: int, xp, mi: int | None = None):
        self.n, self.p, self.xp, self.mi = n, p, xp, mi
        self.piv = np.zeros(0, dtype=I64)
        self.free = np.arange(n, dtype=I64)
        self.B = xp.zeros((0, n), dtype=xp.float64)
        self.madds = 0
        self._index()

    def _index(self):
        self.rowof = -np.ones(self.n, dtype=I64)
        self.rowof[self.piv] = np.arange(len(self.piv))
        self.freepos = -np.ones(self.n, dtype=I64)
        self.freepos[self.free] = np.arange(len(self.free))
        self.cache = {}

    @property
    def rank(self):
        return len(self.piv)

    def snapshot(self):
        return copy.copy(self)

    def reduce(self, V, r0: int = 0):
        """residual on ``free`` of full-width rows V (n columns, entries in [0, p)) that vanish on
        the pivots of the rows < r0."""
        xp = self.xp
        R = V[:, xp.asarray(self.free)]
        if self.rank > r0 and R.shape[1]:
            R = R - mm(V[:, xp.asarray(self.piv[r0:])], self.B[r0:], self.p, xp, self.mi)
            modp(R, self.p, xp)
            self.madds += V.shape[0] * (self.rank - r0) * R.shape[1]
        return R

    def add(self, V, r0: int = 0) -> int:
        xp, p = self.xp, self.p
        R = self.reduce(V, r0)
        R = R[xp.any(R != 0, axis=1)]
        if R.shape[0] == 0:
            return 0
        L, q = rref(R, p, xp, mi=self.mi)
        if len(q) == 0:
            return 0
        keep = np.ones(len(self.free), dtype=bool)
        keep[q] = False
        kd = xp.asarray(np.flatnonzero(keep))
        Lk = L[:, kd]
        old = self.B[:, kd] - mm(self.B[:, xp.asarray(q)], Lk, p, xp, self.mi)
        modp(old, p, xp)
        self.madds += self.rank * len(q) * Lk.shape[1]
        self.B = xp.concatenate([old, Lk])
        self.piv = np.concatenate([self.piv, self.free[q]])
        self.free = self.free[keep]
        self._index()
        return len(q)


# ------------------------------------------------------------------------------- instance

def s3_generator(p: int, a: int, b: int, X: int, k: int) -> dict:
    """NF of Semaev's S_3(x1, x2, X) for y^2 = x^3 + a x + b modulo x_i^(2^k) = 1, built here by
    polynomial arithmetic from
        (x1 - x2)^2 X^2 - 2 ((x1 + x2)(x1 x2 + a) + 2b) X + (x1 x2 - a)^2 - 4b (x1 + x2)."""
    def mul(f, g):
        h = {}
        for (i, j), c in f.items():
            for (s, t), d in g.items():
                h[(i + s, j + t)] = (h.get((i + s, j + t), 0) + c * d) % p
        return h

    def add(*fs):
        h = {}
        for f in fs:
            for e, c in f.items():
                h[e] = (h.get(e, 0) + c) % p
        return h

    def sc(c, f):
        return {e: c * v % p for e, v in f.items()}

    x1, x2, one = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}
    d = add(x1, sc(-1, x2))
    s, q = add(x1, x2), mul(x1, x2)
    qa = add(q, sc(-a, one))
    f = add(sc(X * X, mul(d, d)), sc(-2 * X, add(mul(s, add(q, sc(a, one))), sc(2 * b, one))),
            mul(qa, qa), sc(-4 * b, s))
    N = 1 << k
    g = {}
    for (i, j), c in f.items():
        g[(i % N, j % N)] = (g.get((i % N, j % N), 0) + c) % p
    return {e: c for e, c in g.items() if c}


# ------------------------------------------------------------------------------- monomial tables

def _popcounts(k: int) -> np.ndarray:
    x = np.arange(1 << k, dtype=I64)
    c = np.zeros(1 << k, dtype=I64)
    for j in range(k):
        c += (x >> j) & 1
    return c


def monomials(k: int, weights) -> tuple:
    """all (a1, a2) with popcount(a1) + popcount(a2) in ``weights``, ordered by weight descending,
    then by a1 | a2 << k (the order the engine's parts inherit)."""
    pc = _popcounts(k)
    byw = [np.flatnonzero(pc == w) for w in range(max(weights) + 1)]
    A1, A2, W = [], [], []
    for w in weights:
        for w1 in range(w + 1):
            A1.append(np.repeat(byw[w1], len(byw[w - w1])))
            A2.append(np.tile(byw[w - w1], len(byw[w1])))
            W.append(np.full(len(A1[-1]), w, dtype=I64))
    a1, a2, w = np.concatenate(A1), np.concatenate(A2), np.concatenate(W)
    o = np.lexsort((a1 | (a2 << k), -w))
    return a1[o], a2[o], w[o]


class Table:
    def __init__(self, a1, a2, w, k):
        self.a1, self.a2, self.w, self.k, self.n = a1, a2, w, k, len(a1)
        key = a1 | (a2 << k)
        self._o = np.argsort(key)
        self._s = key[self._o]

    def find(self, b1, b2) -> np.ndarray:
        """index of each (b1, b2), or -1."""
        key = np.asarray(b1) | (np.asarray(b2) << self.k)
        if self.n == 0:
            return -np.ones(len(key), dtype=I64)
        pos = np.minimum(np.searchsorted(self._s, key), self.n - 1)
        return np.where(self._s[pos] == key, self._o[pos], -1)


# ------------------------------------------------------------------------------- checkpoint reader

def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 24), b''):
            h.update(blk)
    return h.hexdigest()


def read_checkpoint(path: str, verify_files: bool = True) -> dict:
    """Engine checkpoint (directory ck-* or its parent): manifest.json, p{q}_host.npz (piv) and
    float32 tiles p{q}_t#####.npy holding the rows on the free columns of part q."""
    if not os.path.exists(os.path.join(path, 'manifest.json')):
        slots = sorted(d for d in os.listdir(path) if d.startswith('ck-'))
        if not slots:
            raise Refuse(f'no checkpoint under {path}')
        path = os.path.join(path, slots[-1])
    man = json.load(open(os.path.join(path, 'manifest.json')))
    cfg = man['config']
    if not cfg['split'] or len(man['ntiles']) != 2:
        raise Refuse('checkpoint is not a sigma-split basis')
    if not man['state'].get('closure_complete'):
        raise Refuse('checkpoint was written before the closure finished')
    files = {}
    for name, (size, sha) in man['files'].items():
        fp = os.path.join(path, name)
        if not os.path.exists(fp) or os.path.getsize(fp) != size:
            raise Refuse(f'{name}: missing or wrong size')
        if verify_files:
            if _sha256_file(fp) != sha:
                raise Refuse(f'{name}: sha256 differs from the manifest')
        files[name] = sha
    p = cfg['p']
    parts = []
    for q in (0, 1):
        piv = np.load(os.path.join(path, f'p{q}_host.npz'))['piv'].astype(I64)
        B, r = None, 0
        for ti in range(man['ntiles'][q]):
            t = np.load(os.path.join(path, f'p{q}_t{ti:05d}.npy'))
            if B is None:
                B = np.empty((len(piv), t.shape[1]), dtype=np.float32)
            if t.size and not (t.min() >= 0 and t.max() < p and np.array_equal(t, np.floor(t))):
                raise Refuse(f'part {q} tile {ti}: entries are not integers in [0, p)')
            B[r:r + len(t)] = t
            r += len(t)
        if B is None:
            B = np.zeros((0, 0), dtype=np.float32)
        if r != len(piv):
            raise Refuse(f'part {q}: {r} rows in the tiles, {len(piv)} pivots')
        parts.append((piv, B))
    basis_id = hashlib.sha256(json.dumps(sorted(files.items())).encode()).hexdigest()
    return dict(k=cfg['k'], p=p, D=cfg['D'], parts=parts, path=path, files_verified=verify_files,
                manifest_files_id=basis_id, engine_sha256=man.get('engine_sha256'), config=cfg)


def canonical_sha(piv: np.ndarray, B: np.ndarray, ncols: int, p: int) -> str:
    """the engine's Part.sha(): pivots ascending, free columns, rows in pivot order as uint16/uint32."""
    o = np.argsort(piv)
    free = np.setdiff1d(np.arange(ncols, dtype=I64), piv)
    h = hashlib.sha256()
    h.update(piv[o].astype(I64).tobytes())
    h.update(free.astype(I64).tobytes())
    dt = np.uint16 if p < 65536 else np.uint32
    step = max(1, (1 << 24) // max(1, B.shape[1]))
    for a in range(0, len(o), step):
        h.update(np.ascontiguousarray(B[o[a:a + step]].astype(dt)).tobytes())
    return h.hexdigest()


# ------------------------------------------------------------------------------- the lift

class Lift:
    def __init__(self, k: int, D: int, p: int, parts, g: dict, xp=np, include_low: bool = False,
                 batch: int = 2048, table_gb: float = 2.0, first_chunk: int = 2048, seed: int = 1,
                 mi: int | None = None, log=None, max_seconds: float | None = None, stop_on_one: bool = False,
                 basis_on: str = 'host'):
        self.k, self.D, self.p, self.g, self.xp = k, D, p, g, xp
        self.include_low, self.batch, self.table_gb, self.first_chunk = include_low, batch, table_gb, first_chunk
        self.seed, self.mi = seed, mi
        self.log = log or (lambda **kw: None)
        self.t0 = time.time()
        self.max_seconds, self.stop_on_one, self.stopped = max_seconds, stop_on_one, None
        self.cnt = dict(generator_rows=0, fall_images=0, followers=0, chunks=0, fall_rounds=0, madds_gemm=0, madds_red=0)
        if p >= 1 << 24 or p % 2 == 0:
            raise Refuse('p must be odd and below 2^24')
        self.inv2 = float(pow(2, -1, p))
        N = 1 << k
        self.shifts = [(i, j) for i in (0, 1) for j in range(k)]
        # --- tables: L = monomials of weight <= D, tops = weight D + 1
        self.L = L = Table(*monomials(k, range(D, -1, -1)), k)
        self.Tp = Table(*monomials(k, [D + 1]), k)
        sig = L.find(L.a2, L.a1)
        selfp = L.a1 == L.a2
        # --- the two stored parts (orbit columns in the inherited order; + has a1 <= a2, - has a1 < a2)
        self.Bq, rows, fpos, self.ncol = [], [], [], []
        for q, (piv, B) in enumerate(parts):
            reps = np.flatnonzero(L.a1 <= L.a2 if q == 0 else L.a1 < L.a2)
            nc = len(reps)
            if len(piv) == 0:                                   # an empty part has no tile file
                B = np.zeros((0, nc), dtype=np.float32)
            if len(np.unique(piv)) != len(piv) or (len(piv) and (piv.min() < 0 or piv.max() >= nc)) \
                    or B.shape != (len(piv), nc - len(piv)):
                raise Refuse(f'part {q}: shape of the basis does not fit {nc} columns')
            col = -np.ones(L.n, dtype=I64)
            col[reps] = np.arange(nc)
            col = np.where(col >= 0, col, col[sig])            # column of the orbit (-1: no column)
            rowof = -np.ones(nc + 1, dtype=I64)
            rowof[piv] = np.arange(len(piv))
            fp = -np.ones(nc + 1, dtype=I64)
            fp[np.setdiff1d(np.arange(nc), piv)] = np.arange(nc - len(piv))
            rows.append(rowof[col])                             # col = -1 reads the padding entry -1
            fpos.append(fp[col])
            self.Bq.append(xp.asarray(B) if basis_on == 'device' else B)   # float32 either way
            self.ncol.append(nc)
        self.rP, self.rM = rows
        self.selfp = selfp
        self.sgnL = np.where(L.a1 <= L.a2, 1.0, -1.0)
        self.Pd = (self.rP >= 0) & (selfp | (self.rM >= 0))
        # --- N and its gather indices into the free columns of each part (-1 -> zero column)
        self.Nmono = np.flatnonzero(~self.Pd)
        self.nN = len(self.Nmono)
        self.nidx = -np.ones(L.n, dtype=I64)
        self.nidx[self.Nmono] = np.arange(self.nN)
        self.gP, self.gM = fpos[0][self.Nmono], fpos[1][self.Nmono]
        self.sgnN = self.sgnL[self.Nmono]
        self.sigN = self.nidx[sig[self.Nmono]]
        if np.any(self.sigN < 0):
            raise Refuse('N is not stable under the swap')
        # --- single rows (orbit a pivot of exactly one part): representatives
        rep = L.a1 <= L.a2
        self.singleP = np.flatnonzero(rep & ~self.Pd & (self.rP >= 0))
        self.singleM = np.flatnonzero(rep & ~self.Pd & (self.rM >= 0))
        self.sig = sig
        # --- tops: special = u n for n in N
        special = np.zeros(self.Tp.n, dtype=bool)
        for s in range(2 * k):
            b1, b2 = self._shift(L.a1[self.Nmono], L.a2[self.Nmono], s)
            t = self.Tp.find(b1, b2)
            special[t[t >= 0]] = True
        self.spidx = -np.ones(self.Tp.n, dtype=I64)
        self.spidx[special] = np.arange(int(special.sum()))
        self.nSp = int(special.sum())
        self.regular = np.flatnonzero(~special)
        self.nT = self.nSp + self.nN
        # --- sources of generators and the direct ones (image not a regular top)
        src = np.flatnonzero(self.Pd & ((L.w == D) | include_low))
        hits = np.zeros(self.Tp.n, dtype=I64)
        dx, du = [], []
        for s in range(2 * k):
            _, _, reg = self.image(src, s)
            np.add.at(hits, reg[reg >= 0], 1)
            dx.append(src[reg < 0])
            du.append(np.full(int((reg < 0).sum()), s, dtype=I64))
        self.dir_x, self.dir_u = np.concatenate(dx), np.concatenate(du)
        self.n_sources = len(src)
        if not np.all(hits[self.regular] == D + 1) or np.any(hits[special]):
            raise Refuse('a regular top does not have exactly D + 1 generators')
        self.one = int(L.find(np.array([0]), np.array([0]))[0])
        self.one_in_VD = bool(self.Pd[self.one])
        self.K = Echelon(self.nT, p, xp, mi)
        self.done = np.zeros(self.nT, dtype=bool)       # fall pivots already shifted
        self.saturated = False

    # -- monomial maps
    def _shift(self, a1, a2, s):
        i, j = self.shifts[s]
        m = (1 << self.k) - 1
        return ((a1 + (1 << j)) & m, a2) if i == 0 else (a1, (a2 + (1 << j)) & m)

    def image(self, mono, s):
        """u_s times the monomials ``mono`` of L: (T column or -1, Pd monomial or -1, regular top or -1);
        exactly one of the three is >= 0."""
        b1, b2 = self._shift(self.L.a1[mono], self.L.a2[mono], s)
        li, ti = self.L.find(b1, b2), self.Tp.find(b1, b2)
        if not np.all((li >= 0) ^ (ti >= 0)):
            raise Refuse('a shifted monomial is outside M_<=D+1')
        inL = li >= 0
        nn = np.where(inL, self.nidx[np.maximum(li, 0)], -1)
        # At D >= 2k there is no top layer. np.where eagerly indexes both
        # branches, so avoid indexing the empty special-top array.
        sp = (-np.ones(len(li), dtype=I64) if self.Tp.n == 0 else
              np.where(inL, -1, self.spidx[np.maximum(ti, 0)]))
        tcol = np.where(nn >= 0, self.nSp + nn, np.where(sp >= 0, sp, -1))
        pd = np.where(inL & (nn < 0), li, -1)
        reg = np.where(~inL & (sp < 0), ti, -1)
        return tcol, pd, reg

    # -- reading the stored basis in ordinary coordinates
    def _rows(self, q, idx):
        B = self.Bq[q]
        return B[idx] if isinstance(B, np.ndarray) else B[self.xp.asarray(idx)]

    def _gatherN(self, A, idx, sgn=None):
        """rows A (float32) on the free columns of a part -> device float64 rows on N."""
        xp = self.xp
        E = xp.zeros((A.shape[0], A.shape[1] + 1), dtype=xp.float64)
        E[:, :A.shape[1]] = xp.asarray(A)
        E = E[:, xp.asarray(idx)]                     # idx = -1 reads the zero column
        if sgn is not None:
            E *= xp.asarray(sgn)[None, :]
        return E

    def phi(self, mono):
        """phi(x) on N for monomials x in Pd (device, entries in [0, p))."""
        xp, p = self.xp, self.p
        rp, rm = self.rP[mono], self.rM[mono]
        if np.any(rp < 0) or np.any((rm < 0) & ~self.selfp[mono]):
            raise Refuse('phi asked for a monomial outside Pd')
        out = self._gatherN(self._rows(0, rp), self.gP)
        pr = np.flatnonzero(rm >= 0)
        if len(pr):
            M = self._gatherN(self._rows(1, rm[pr]), self.gM, self.sgnN)
            M *= xp.asarray(self.sgnL[mono[pr]])[:, None]
            prd = xp.asarray(pr)
            X = out[prd] + M                           # in (-p, 2p)
            X *= self.inv2
            out[prd] = modp(X, p, xp)
        return out

    def red_L(self, C):
        """red of rows C given on all monomials of L (host, t x |L|, entries in [0, p)): rows on N.
        Uses products with whole parts (a second route to phi, compared with ``phi`` in the tests)."""
        xp, p, L = self.xp, self.p, self.L
        C = np.asarray(C, dtype=np.float64)
        rep = np.flatnonzero((L.a1 <= L.a2) & self.Pd)
        cs = C[:, self.sig[rep]]
        pair = (~self.selfp[rep]).astype(np.float64)
        Cp = np.zeros((C.shape[0], len(self.Bq[0])))
        Cm = np.zeros((C.shape[0], len(self.Bq[1])))
        Cp[:, self.rP[rep]] = np.mod((C[:, rep] + cs * pair) * np.where(pair > 0, self.inv2, 1.0), p)
        pr = rep[pair > 0]
        Cm[:, self.rM[pr]] = np.mod((C[:, pr] - C[:, self.sig[pr]]) * self.inv2, p)
        out = xp.asarray(C[:, self.Nmono])
        for q, (Cq, idx, sgn) in enumerate(((Cp, self.gP, None), (Cm, self.gM, self.sgnN))):
            B = self.Bq[q]
            acc = xp.zeros((C.shape[0], B.shape[1]), dtype=xp.float64)
            step = max(1, (1 << 27) // max(1, B.shape[1]))
            for a in range(0, len(B), step):
                acc += mm(xp.asarray(Cq[:, a:a + step]), xp.asarray(B[a:a + step], dtype=xp.float64), p, xp, self.mi)
                modp(acc, p, xp)
            E = xp.zeros((C.shape[0], B.shape[1] + 1), dtype=xp.float64)
            E[:, :B.shape[1]] = acc
            E = E[:, xp.asarray(idx)]
            if sgn is not None:
                E *= xp.asarray(sgn)[None, :]
            out -= E
            modp(out, p, xp)
        return out

    def single_rows(self):
        """the stored rows that are supported on N (a basis of V_N), as rows on N."""
        xp, p = self.xp, self.p
        out = []
        for q, reps, rw, idx, sgn in ((0, self.singleP, self.rP, self.gP, None), (1, self.singleM, self.rM, self.gM, self.sgnN)):
            if len(reps) == 0:
                continue
            E = self._gatherN(self._rows(q, rw[reps]), idx, sgn)
            ar = xp.arange(len(reps))
            E[ar, xp.asarray(self.nidx[reps])] += 1.0
            E[ar, xp.asarray(self.nidx[self.sig[reps]])] += 1.0 if q == 0 else -1.0
            out.append(modp(E, p, xp))
        return xp.concatenate(out) if out else xp.zeros((0, self.nN), dtype=xp.float64)

    # -- T coordinates and reduction modulo a snapshot K0 of K
    def embedN(self, VN):
        xp = self.xp
        E = xp.zeros((VN.shape[0], self.nT), dtype=xp.float64)
        E[:, self.nSp:] = VN
        return E

    def embed(self, V, K0):
        xp = self.xp
        E = xp.zeros((V.shape[0], self.nT), dtype=xp.float64)
        E[:, xp.asarray(K0.free)] = V
        return E

    def red_units(self, K0, tcol):
        """residual modulo K0 of the unit vectors e_c, c in tcol."""
        xp, p = self.xp, self.p
        out = xp.zeros((len(tcol), len(K0.free)), dtype=xp.float64)
        fr = np.flatnonzero(K0.freepos[tcol] >= 0)
        out[xp.asarray(fr), xp.asarray(K0.freepos[tcol[fr]])] = 1.0
        pv = np.flatnonzero(K0.rowof[tcol] >= 0)
        if len(pv) and out.shape[1]:
            out[xp.asarray(pv)] = neg(K0.B[xp.asarray(K0.rowof[tcol[pv]])], p)
        return out

    def red_N(self, K0, VN):
        """residual modulo K0 of rows VN supported on N.  Only the fall rows of K0 (pivot in N)
        can act: the other pivots are Sp columns, where VN vanishes."""
        xp, p = self.xp, self.p
        c = K0.cache
        if 'fall' not in c:
            fr = np.flatnonzero(K0.piv >= self.nSp)
            nf = np.flatnonzero(K0.free >= self.nSp)
            c['fall'] = (xp.asarray(K0.piv[fr] - self.nSp), K0.B[xp.asarray(fr)] if len(fr) < K0.rank else K0.B,
                         xp.asarray(nf), xp.asarray(K0.free[nf] - self.nSp), len(fr))
        pivN, Bf, pos, loc, nf_rows = c['fall']
        out = xp.zeros((VN.shape[0], len(K0.free)), dtype=xp.float64)
        out[:, pos] = VN[:, loc]
        if nf_rows and out.shape[1]:
            out -= mm(VN[:, pivN], Bf, p, xp, self.mi)
            modp(out, p, xp)
            self.cnt['madds_red'] += VN.shape[0] * nf_rows * out.shape[1]
        return out

    def follower(self, s, K0, tile: int | None = None):
        """R_s modulo K0: row n = residual of red(e_{u_s n}), n in N."""
        xp, p = self.xp, self.p
        tile = tile or self.batch
        tcol, pd, reg = self.image(self.Nmono, s)
        if np.any(reg >= 0):
            raise Refuse('u N meets a regular top')
        F = xp.zeros((self.nN, len(K0.free)), dtype=xp.float64)
        u = np.flatnonzero(tcol >= 0)
        for a in range(0, len(u), 8 * tile):
            F[xp.asarray(u[a:a + 8 * tile])] = self.red_units(K0, tcol[u[a:a + 8 * tile]])
        d = np.flatnonzero(pd >= 0)
        for a in range(0, len(d), tile):
            F[xp.asarray(d[a:a + tile])] = neg(self.red_N(K0, self.phi(pd[d[a:a + tile]])), p)
        self.cnt['followers'] += 1
        return F

    def _add(self, V, K0) -> int:
        """add rows given on K0.free (reduced modulo the snapshot K0) to K."""
        V = V[self.xp.any(V != 0, axis=1)]
        if V.shape[0] == 0:
            return 0
        n = self.K.add(self.embed(V, K0), r0=K0.rank)
        self.saturated = self.K.rank == self.nT
        self._check_stop()
        return n

    def _check_stop(self):
        if self.saturated:
            return
        if self.stop_on_one and not self.one_in_VD and self.K.rowof[self.nSp + int(self.nidx[self.one])] >= 0:
            raise Stop('one')
        if self.max_seconds and time.time() - self.t0 > self.max_seconds:
            raise Stop('deadline')

    # -- steps
    def start(self, closed_check: int = 2):
        """K = V_N; g must lie in V_D; the generators skipped by Lemma 2 must lie in V_N."""
        xp, p, L = self.xp, self.p, self.L
        VN = self.single_rows()
        for a in range(0, VN.shape[0], self.batch):
            self.K.add(self.embedN(VN[a:a + self.batch]))
        if self.K.rank != VN.shape[0]:
            raise Refuse('the single rows are dependent')
        self.n_single = int(VN.shape[0])
        gv = np.zeros((1, L.n))
        e = np.array(list(self.g.keys()), dtype=I64).reshape(-1, 2)
        gi = L.find(e[:, 0], e[:, 1])
        if np.any(gi < 0):
            raise Refuse('generator has a monomial of weight > D')
        gv[0, gi] = np.array(list(self.g.values()), dtype=np.float64) % p
        if bool(xp.any(self.K.reduce(self.embedN(self.red_L(gv))) != 0)):
            raise Refuse('the generator g is not in the stored space')
        self.closed_check = None
        self.log(event='start', n_single=self.n_single, generator_in_VD=True, sec=round(time.time() - self.t0, 1))
        if self.include_low or closed_check <= 0:
            return
        # low rows l_x (x in Pd, weight < D) must be supported on weight < D ...
        lowP = np.flatnonzero(self.Pd & (L.w < self.D))
        topN = xp.asarray(np.flatnonzero(L.w[self.Nmono] == self.D))
        for a in range(0, len(lowP), self.batch):
            if bool(xp.any(self.phi(lowP[a:a + self.batch])[:, topN] != 0)):
                raise Refuse('a stored row with pivot of weight < D has support in weight D')
        # ... and u (random combination of them) must lie in V_D for every shift u
        rng = np.random.RandomState(self.seed + 1000)
        t = closed_check
        w = np.zeros((t, L.n))
        w[:, lowP] = rng.randint(0, p, size=(t, len(lowP)))
        w[:, self.Nmono] = to_host(neg(self.red_L(w), p))      # w = sum c_x (e_x + phi(x))
        low = np.flatnonzero(L.w < self.D)
        sh = np.zeros((2 * self.k * t, L.n))
        for s in range(2 * self.k):
            b1, b2 = self._shift(L.a1[low], L.a2[low], s)
            sh[s * t:(s + 1) * t, L.find(b1, b2)] = w[:, low]          # a shift is injective
        bad = int(xp.count_nonzero(xp.any(self.K.reduce(self.embedN(self.red_L(sh))) != 0, axis=1)))
        self.closed_check = dict(random_combinations=t, shifts=2 * self.k, nonzero=bad)
        if bad:
            raise Refuse(f'the stored space is not closed: {bad} of {2 * self.k * t} shifted random low elements leave it')

    def chunk(self, tops, direct):
        """generator differences of the regular tops ``tops`` and the direct generators ``direct``."""
        xp, p, L, D = self.xp, self.p, self.L, self.D
        K0 = self.K.snapshot()
        f0 = len(K0.free)
        Ref = xp.zeros((len(tops), f0), dtype=xp.float64)
        have = np.zeros(len(tops) + 1, dtype=bool)                # one spare entry: slot -1 reads entry 0 under a mask
        ta = (self.Tp.a1[tops], self.Tp.a2[tops])
        for s, (i, j) in enumerate(self.shifts):
            slot = np.flatnonzero((ta[i] >> j) & 1)                 # tops divisible by u_s
            b = [ta[0][slot], ta[1][slot]]
            b[i] = b[i] ^ (1 << j)
            xs = L.find(b[0], b[1])
            if np.any(xs < 0) or not np.all(self.Pd[xs]):
                raise Refuse('a divisor of a regular top is not in Pd')
            dsel = direct[self.dir_u[direct] == s]
            xs = np.concatenate([xs, self.dir_x[dsel]])
            slot = np.concatenate([slot, -np.ones(len(dsel), dtype=I64)])
            if len(xs) == 0:
                continue
            F = self.follower(s, K0)
            for a in range(0, len(xs), self.batch):
                x, sl = xs[a:a + self.batch], slot[a:a + self.batch]
                G = mm(self.phi(x), F, p, xp, self.mi)
                self.cnt['generator_rows'] += len(x)
                self.cnt['madds_gemm'] += len(x) * self.nN * f0
                tcol, pd, reg = self.image(x, s)
                if not np.array_equal(reg >= 0, sl >= 0):
                    raise Refuse('generator type mismatch')
                un, dn = np.flatnonzero(tcol >= 0), np.flatnonzero(pd >= 0)
                if len(un):
                    G[xp.asarray(un)] += self.red_units(K0, tcol[un])
                if len(dn):
                    G[xp.asarray(dn)] += neg(self.red_N(K0, self.phi(pd[dn])), p)
                modp(G, p, xp)
                isreg = sl >= 0
                first = np.flatnonzero(isreg & ~have[np.maximum(sl, 0)])
                other = np.flatnonzero(isreg & have[np.maximum(sl, 0)])
                cand = [G[xp.asarray(np.flatnonzero(~isreg))]]
                if len(other):
                    C = G[xp.asarray(other)] - Ref[xp.asarray(sl[other])]
                    cand.append(modp(C, p, xp))
                if len(first):
                    Ref[xp.asarray(sl[first])] = G[xp.asarray(first)]
                    have[sl[first]] = True
                self._add(xp.concatenate(cand), K0)
                if self.saturated:
                    return
            self.log(event='chunk_shift', shift=s, of=2 * self.k, rows=len(xs), rank_K=self.K.rank, nT=self.nT,
                     sec=round(time.time() - self.t0, 1))
        self.cnt['chunks'] += 1

    def fall_rows(self, cols):
        """rows of K with the given pivot columns (all in N), as rows on N."""
        xp, K = self.xp, self.K
        out = xp.zeros((len(cols), self.nN), dtype=xp.float64)
        out[xp.arange(len(cols)), xp.asarray(cols - self.nSp)] = 1.0
        nf = np.flatnonzero(K.free >= self.nSp)
        if len(nf):
            out[:, xp.asarray(K.free[nf] - self.nSp)] = K.B[xp.asarray(K.rowof[cols])][:, xp.asarray(nf)]
        return out

    def falls(self, verify: bool = False) -> int:
        """shift every fall row not yet shifted (verify: all of them, and nothing may be new)."""
        xp, p = self.xp, self.p
        new = 0
        while not self.saturated:
            cols = self.K.piv[self.K.piv >= self.nSp]
            if not verify:
                cols = cols[~self.done[cols]]
            if len(cols) == 0:
                break
            self.done[cols] = True
            self.cnt['fall_rounds'] += 1
            for s in range(2 * self.k):
                K0 = self.K.snapshot()
                F = self.follower(s, K0)
                for a in range(0, len(cols), self.batch):
                    img = mm(self.fall_rows(cols[a:a + self.batch]), F, p, xp, self.mi)
                    self.cnt['fall_images'] += img.shape[0]
                    self.cnt['madds_gemm'] += img.shape[0] * self.nN * img.shape[1]
                    new += self._add(img, K0)
                    if self.saturated:
                        return new
                self.log(event='fall_shift', shift=s, of=2 * self.k, rows=len(cols), rank_K=self.K.rank, nT=self.nT,
                         verify=verify, sec=round(time.time() - self.t0, 1))
            self.log(event='fall_round', rows=len(cols), rank_K=self.K.rank, nT=self.nT, verify=verify,
                     sec=round(time.time() - self.t0, 1))
            if verify:
                break
        return new

    def run(self, closed_check: int = 2) -> dict:
        xp = self.xp
        self.start(closed_check)
        self.saturated = self.K.rank == self.nT
        rng = np.random.RandomState(self.seed)
        tops = self.regular[rng.permutation(len(self.regular))]
        direct = rng.permutation(len(self.dir_x))
        done, total, n = 0.0, max(1, len(tops)), self.first_chunk
        at = ad = 0
        complete, verify_new = False, None
        try:
            self._check_stop()
            self.falls()
            while not self.saturated and (at < len(tops) or ad < len(direct)):
                cap = int(self.table_gb * 2 ** 30 / (8 * max(1, self.nT - self.K.rank)))
                n = max(1, min(n, cap))
                done = min(1.0, done + n / total)
                bt = len(tops) if done >= 1.0 else int(done * len(tops))
                bd = len(direct) if done >= 1.0 else int(done * len(direct))
                self.chunk(tops[at:bt], direct[ad:bd])
                at, ad = bt, bd
                self.falls()
                self.log(event='chunk', tops_done=at, tops=len(tops), direct_done=ad, direct=len(direct),
                         rank_K=self.K.rank, nT=self.nT, sec=round(time.time() - self.t0, 1))
                n *= 2
            complete = at == len(tops) and ad == len(direct)
            if not self.saturated:
                self.max_seconds = None                       # the final check is not interrupted
                self.stop_on_one = False
                verify_new = self.falls(verify=True)
                if verify_new:
                    raise Refuse('fixed-point verification failed: a shifted fall row left K')
        except Stop as e:
            self.stopped, complete = str(e), False
        return self.result(complete, verify_new)

    def result(self, complete: bool, verify_new) -> dict:
        xp, p, L = self.xp, self.p, self.L
        K = self.K
        fcols = K.piv[K.piv >= self.nSp]
        nf = len(fcols)
        # sigma-parts of the fall space: ranks of Phi (1 + sigma) and Phi (1 - sigma)
        ranks = []
        for sign in (1.0, -1.0):
            E = Echelon(self.nN, p, xp, self.mi)
            for a in range(0, nf, self.batch):
                Fl = self.fall_rows(fcols[a:a + self.batch])
                E.add(modp(Fl + sign * Fl[:, xp.asarray(self.sigN)], p, xp))
            ranks.append(E.rank)
        rep = (L.a1 <= L.a2) & self.Pd
        pd_parts = [int(rep.sum()), int((rep & ~self.selfp).sum())]
        onecol = self.nSp + int(self.nidx[self.one])
        one_new = (not self.one_in_VD) and K.rowof[onecol] >= 0
        fixed = bool(self.saturated or complete)
        npd = int(self.Pd.sum())
        return dict(
            version=VERSION, k=self.k, D=self.D, D_target=self.D + 1, p=p,
            contains_one=bool(self.one_in_VD or one_new), one_in_VD=self.one_in_VD,
            fixed_point=fixed, saturated=bool(self.saturated), schedule_complete=bool(complete),
            stopped=self.stopped, bounds_only=not fixed,       # bounds_only: only contains_one = true is a statement
            verify_new_rows=verify_new,
            sigma_stable=bool(ranks[0] + ranks[1] == nf),
            dim_VD=npd + self.n_single,
            dim_low=npd + nf,                                   # dim V_{D+1} n M_<=D
            dim_low_parts=[pd_parts[0] + ranks[0], pd_parts[1] + ranks[1]],
            dim_VD1=npd + len(self.regular) + K.rank,
            columns=dict(L=L.n, top=self.Tp.n, total=L.n + self.Tp.n),
            n_Pd=npd, n_N=self.nN, n_N_weight_D=int(np.sum(L.w[self.Nmono] == self.D)),
            n_single=self.n_single, n_special_tops=self.nSp,
            n_regular_tops=len(self.regular), n_T=self.nT, rank_K=K.rank, falls=nf,
            falls_new=nf - self.n_single, n_sources=self.n_sources, n_direct=len(self.dir_x),
            include_low=self.include_low, closed_check=self.closed_check,
            counts=dict(self.cnt, madds_K=K.madds), seed=self.seed,
            seconds=round(time.time() - self.t0, 1))


def plan(nN: int, nSp: int, n_regular: int, rank_parts, free_parts, batch: int, table_gb: float,
         basis_on: str = 'host') -> dict:
    """Memory estimate (GB) from the counts, worst case over the rank r of K (f = nT - r free columns).
    Device arrays are float64: one follower nN x f; K as r x f (three copies while it is updated, one
    more for the snapshot a chunk works against); the fall rows of the snapshot nf x f; the reference
    table of a chunk; per batch about three nN-wide and four nT-wide blocks; at the end one echelon
    form on the nN columns of N for the sigma-parts."""
    nT = nN + nSp
    G = 2.0 ** 30
    basis = sum(r * f * 4 for r, f in zip(rank_parts, free_parts)) / G
    dev = dict(follower=nN * nT * 8 / G, K_worst=4 * (nT / 2) ** 2 * 8 / G, fall_rows_worst=(nN / 2) * (nT / 2) * 8 / G,
               table=table_gb, batch=batch * (3 * nN + 4 * nT + 2 * max(free_parts) + 2) * 8 / G,
               sigma_parts_at_end=3 * (nN / 2) ** 2 * 8 / G, basis=basis if basis_on == 'device' else 0.0)
    dev['total'] = sum(dev.values())
    return dict(n_N=nN, n_T=nT, basis_gb=round(basis, 1), basis_on=basis_on,
                host_gb=round((0.0 if basis_on == 'device' else basis) + 3.0, 1),
                device_gb={k: round(v, 1) for k, v in dev.items()},
                first_chunk_cap_tops=int(table_gb * G / (8 * nT)), regular_tops=n_regular,
                note='sizes are upper bounds; K_worst and fall_rows_worst are at r = f')
