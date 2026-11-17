"""Point decomposition on binary curves by linear algebra and one quadratic: no Groebner basis.

Problem (relation search of index calculus, csd.ic_dreg): E: y^2 + xy = x^3 + a x^2 + b over F_(2^n), V an
F_2-subspace of dimension k, factor base FB = {P : x(P) in V}.  Find every R = P_1 + ... + P_m, P_i in FB.

Two summands are linear.  P_1 + P_2 = T iff P_1, P_2, -T lie on a line y = lam (x + xT) + yT + xT.  Substituting
into the curve gives x^3 + (lam^2 + lam + a) x^2 + mu x + mu^2 + b with mu = lam xT + yT + xT, so

    e1 = x_1 + x_2 = lam^2 + lam + a + xT,          e2 = x_1 x_2 = xT lam^2 + ((yT + xT)^2 + b) / xT.

Both are F_2-AFFINE in lam.  x_1, x_2 in V forces e1 in V and e2 in V2 = span(V V).  Write e1 = sum t_i v_i,
c1 = a + xT, c2 = ((yT + xT)^2 + b)/xT.  Then lam = Q(e1 + c1) + eps with Q a fixed linear right inverse of
z -> z^2 + z (valid iff Tr(e1 + c1) = 0), lam^2 = lam + e1 + c1, and with W = Q + id

    e2 = xT W(e1) + eps xT + xT W(c1) + c2   in V2,       sum t_i Tr(v_i) = Tr(c1).

That is a linear system over F_2 in the k + 1 unknowns (t, eps) with n - dim V2 + 1 equations.  Each solution gives
one quadratic x^2 + e1 x + e2; it is kept when it splits with a root in V (the other root is then in V too).
The work per target is O(k) field multiplications and O(k^2) word operations when k + dim V2 <= n, and
2^(k + dim V2 - n) candidates otherwise; it never exceeds the 2^k of enumerating x_1.  dim V2 is 2k - 1 for
V = <1, w, ..., w^(k-1)> and min(n, k(k+1)/2) for a random V.

m summands: enumerate m - 2 factor-base points and solve the remaining two-summand problem for
T = R - P_3 - ... - P_m.  Cost |FB|^(m-2) linear systems, against |FB|^(m-1) curve additions for enumeration.

Complete: every line through -T except the vertical one is some lam, and the vertical line holds only T, -T, O.
The targets T = O are skipped (they are P + (-P), not decompositions) and x(T) = 0 falls back to enumeration.

The hot loops are in `_rootfind.c` (n <= 63).  `d2_reference` is the same algorithm in plain Python, used to
check the C kernel; `decompose_reference` is plain enumeration with the group law of csd.ic_dreg.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import ic_dreg as I

_SRC = Path(__file__).with_name("_rootfind.c")
u64 = ctypes.c_uint64


class Stats(ctypes.Structure):
    _fields_ = [(f, u64) for f in ("targets", "lin_solvable", "lin_candidates", "quad_split", "found",
                                   "degenerate", "max_kernel_dim")]

    def as_dict(self):
        return {f: int(getattr(self, f)) for f, _ in self._fields_}


@lru_cache(maxsize=1)
def _lib():
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        raise RuntimeError("csd.ic_rootfind needs a C compiler (and a CPU with PCLMULQDQ)")
    src = _SRC.read_bytes()
    tag = hashlib.sha256(src).hexdigest()[:16]
    cache = Path(os.environ.get("CSD_CACHE", Path.home() / ".cache" / "csd"))
    cache.mkdir(parents=True, exist_ok=True)
    so = cache / f"rootfind_{tag}.so"
    if not so.exists():
        tmp = so.with_suffix(f".{os.getpid()}.tmp")
        subprocess.run([cc, "-O3", "-march=native", "-shared", "-fPIC", "-o", str(tmp), str(_SRC)],
                       check=True, capture_output=True)
        os.replace(tmp, so)
    lib = ctypes.CDLL(str(so))
    p, L, i = ctypes.c_void_p, ctypes.c_long, ctypes.c_int
    lib.rf_ctx_size.restype = ctypes.c_size_t
    lib.rf_init.argtypes = [p, i, u64, u64, u64, i, p, i, p, p, p, p, u64]
    lib.rf_init.restype = None
    lib.rf_mul.argtypes = [p, u64, u64]
    lib.rf_mul.restype = u64
    lib.rf_inv.argtypes = [p, u64]
    lib.rf_inv.restype = u64
    lib.rf_add.argtypes = [p, p, p, p]
    lib.rf_add.restype = None
    lib.rf_fb.argtypes = [p, p, p, L]
    lib.rf_fb.restype = L
    lib.rf_lrf.argtypes = [p, i, p, p, p, L, L, L, p, L, p]
    lib.rf_lrf.restype = L
    lib.rf_brute.argtypes = [p, i, p, p, p, L, L, L, p, L]
    lib.rf_brute.restype = L
    lib.rf_mitm4.argtypes = [p, p, p, p, L, p, L, p]
    lib.rf_mitm4.restype = L
    return lib


# ---------------------------------------------------------------------------
# F_2-linear algebra on ints
# ---------------------------------------------------------------------------


def _echelon(vectors):
    """Reduced xor basis {pivot bit: vector}."""
    piv: dict[int, int] = {}
    for v in vectors:
        for pb, pv in piv.items():
            if v >> pb & 1:
                v ^= pv
        if v:
            pb = v.bit_length() - 1
            for qb in list(piv):
                if piv[qb] >> pb & 1:
                    piv[qb] ^= v
            piv[pb] = v
    return piv


def quotient_map(basis, n):
    """Images of the n unit vectors under a surjection F_2^n -> F_2^(n - dim) with kernel span(basis)."""
    piv = _echelon(basis)
    free = [b for b in range(n) if b not in piv]
    pos = {b: j for j, b in enumerate(free)}
    imgs = []
    for b in range(n):
        if b in pos:
            imgs.append(1 << pos[b])
        else:
            # unit vector e_b = piv[b] + (non-pivot bits of piv[b]); piv[b] is in the kernel
            r = piv[b] ^ (1 << b)
            imgs.append(sum(1 << pos[c] for c in range(n) if r >> c & 1))
    return imgs, len(free)


def product_space(F, basis):
    return list(_echelon(F.mul(u, v) for i, u in enumerate(basis) for v in basis[i:]).values())


def artin_schreier_inverse(F):
    """Images of the unit vectors under a linear Q with Q(c)^2 + Q(c) = c whenever Tr(c) = 0."""
    n = F.N
    # echelonize the images z^2 + z of the unit vectors, tracking preimages in the high half
    rows = [((1 << i) << n) | (F.sq(1 << i) ^ (1 << i)) for i in range(n)]
    low = (1 << n) - 1
    piv: dict[int, int] = {}
    for r in rows:
        for pb, pr in piv.items():
            if r >> pb & 1:
                r ^= pr
        if r & low:
            pb = (r & low).bit_length() - 1
            for qb in list(piv):
                if piv[qb] >> pb & 1:
                    piv[qb] ^= r
            piv[pb] = r
    imgs = []
    for b in range(n):
        v, pre = 1 << b, 0
        for pb in sorted(piv, reverse=True):
            if v >> pb & 1:
                v ^= piv[pb] & low
                pre ^= piv[pb] >> n
        imgs.append(pre)          # exact when e_b is in the image; otherwise the preimage of e_b + (a fixed coset rep)
    return imgs


def _apply(imgs, x):
    r = 0
    while x:
        b = (x & -x).bit_length() - 1
        r ^= imgs[b]
        x &= x - 1
    return r


class Solver:
    """Tables for one (curve, V); the methods return decompositions as sorted tuples of points."""

    def __init__(self, E: I.BinaryCurve, basis: list[int]):
        F = E.F
        n, k = F.N, len(basis)
        assert n <= 63 and I._f2_rank(basis) == k
        self.E, self.F, self.n, self.k, self.basis = E, F, n, k, list(basis)
        self.V2 = product_space(F, basis)
        self.dimV2 = len(self.V2)
        self.pi, wq = quotient_map(self.V2, n)
        self.chk, _ = quotient_map(basis, n)
        self.Q = artin_schreier_inverse(F)
        self.sqrt = [F.sqrt(1 << i) for i in range(n)]
        self.trmask = sum(1 << i for i in range(n) if F.trace(1 << i))
        self.w = wq + 1
        lib = _lib()
        self._ctx = ctypes.create_string_buffer(lib.rf_ctx_size())
        arr = lambda xs: (u64 * len(xs))(*xs)
        lib.rf_init(self._ctx, n, F.modulus, E.a, E.b, k, arr(basis), self.w,
                    arr(self.Q), arr(self.sqrt), arr(self.chk), arr(self.pi), self.trmask)
        cap = 2 << k
        xs, ys = np.zeros(cap, dtype=np.uint64), np.zeros(cap, dtype=np.uint64)
        cnt = lib.rf_fb(self._ctx, xs.ctypes.data, ys.ctypes.data, cap)
        assert cnt <= cap
        self.fx, self.fy = xs[:cnt].copy(), ys[:cnt].copy()
        self.nfb = int(cnt)

    # -- C kernels ----------------------------------------------------------

    def _rows(self, out, cnt, cap, m):
        if cnt > cap:
            raise OverflowError(f"{cnt} rows, capacity {cap}")
        rows = out[:cnt * 2 * m].reshape(cnt, m, 2)
        return sorted({tuple(sorted((int(x), int(y)) for x, y in row)) for row in rows})

    def _ranges(self, threads):
        if threads <= 1:
            return [(0, self.nfb)]
        step = max(1, self.nfb // (16 * threads))
        return [(lo, min(lo + step, self.nfb)) for lo in range(0, self.nfb, step)]

    def lrf(self, m: int, R, cap: int = 1 << 20, threads: int = 1):
        """Linearized root finding.  Returns (decompositions, stats, rows found before deduplication).
        With threads > 1 the first enumerated point is split into ranges (m >= 3); cap is per range."""
        lib, Rc = _lib(), (u64 * 2)(*R)

        def part(rng):
            out = np.zeros(cap * 2 * m, dtype=np.uint64)
            st = Stats()
            cnt = lib.rf_lrf(self._ctx, m, Rc, self.fx.ctypes.data, self.fy.ctypes.data, self.nfb,
                             rng[0], rng[1], out.ctypes.data, cap, ctypes.byref(st))
            if cnt < 0:
                raise ValueError("m is too large for the kernel")
            return self._rows(out, cnt, cap, m), st.as_dict(), int(cnt)

        ranges = self._ranges(threads if m >= 3 else 1)
        if len(ranges) == 1:
            return part(ranges[0])
        with ThreadPoolExecutor(threads) as ex:
            parts = list(ex.map(part, ranges))
        stats = {f: (max if f == "max_kernel_dim" else sum)(p[1][f] for p in parts) for f in parts[0][1]}
        return sorted(set().union(*(p[0] for p in parts))), stats, sum(p[2] for p in parts)

    def brute(self, m: int, R, cap: int = 1 << 20, threads: int = 1):
        lib, Rc = _lib(), (u64 * 2)(*R)

        def part(rng):
            out = np.zeros(cap * 2 * m, dtype=np.uint64)
            cnt = lib.rf_brute(self._ctx, m, Rc, self.fx.ctypes.data, self.fy.ctypes.data, self.nfb,
                               rng[0], rng[1], out.ctypes.data, cap)
            return self._rows(out, cnt, cap, m), int(cnt)

        ranges = self._ranges(threads)
        if len(ranges) == 1:
            return part(ranges[0])
        with ThreadPoolExecutor(threads) as ex:
            parts = list(ex.map(part, ranges))
        return sorted(set().union(*(p[0] for p in parts))), sum(p[1] for p in parts)

    def mitm4(self, R, cap: int = 1 << 20):
        out = np.zeros(cap * 8, dtype=np.uint64)
        tb = u64()
        cnt = _lib().rf_mitm4(self._ctx, (u64 * 2)(*R), self.fx.ctypes.data, self.fy.ctypes.data, self.nfb,
                              out.ctypes.data, cap, ctypes.byref(tb))
        if cnt < 0:
            raise MemoryError
        return self._rows(out, cnt, cap, 4), int(cnt), int(tb.value)

    # -- plain Python: the same algorithm, and plain enumeration -------------

    def d2_reference(self, T):
        """Two-summand decompositions of T by the linear system, in Python on csd.ic_dreg's field."""
        F, E = self.F, self.E
        xT, yT = T
        assert xT != 0
        ym = yT ^ xT
        c1 = E.a ^ xT
        c2 = F.mul(F.sq(ym) ^ E.b, F.inv(xT))
        Qf = lambda x: _apply(self.Q, x)
        pif = lambda x: _apply(self.pi, x)
        top = self.w - 1
        cols = [pif(F.mul(xT, Qf(v) ^ v)) | (F.trace(v) << top) for v in self.basis] + [pif(xT)]
        rhs = pif(F.mul(xT, Qf(c1) ^ c1) ^ c2) | (F.trace(c1) << top)
        sol = I._f2_solve(cols, rhs, self.w)
        out = set()
        if sol is None:
            return []
        for s in I._affine_span(*sol):
            e1 = 0
            for i, v in enumerate(self.basis):
                if s >> i & 1:
                    e1 ^= v
            eps = s >> self.k & 1
            lam = Qf(e1 ^ c1) ^ eps
            assert F.sq(lam) ^ lam == e1 ^ c1
            e2 = F.mul(xT, F.sq(lam)) ^ c2
            assert pif(e2) == 0
            if e1 == 0:
                roots = [F.sqrt(e2)] * 2
            else:
                d = F.mul(e2, F.inv(F.sq(e1)))
                if F.trace(d):
                    continue
                x1 = F.mul(e1, Qf(d))
                roots = [x1, x1 ^ e1]
            if _apply(self.chk, roots[0]):
                continue
            mu = F.mul(lam, xT) ^ ym
            P1, P2 = [(x, F.mul(lam, x) ^ mu) for x in roots]
            assert E.on_curve(P1) and E.on_curve(P2) and E.add(P1, P2) == T
            out.add(tuple(sorted((P1, P2))))
        return sorted(out)

    def decompose_reference(self, m: int, R):
        """Every multiset of m factor-base points with sum R, by enumeration with csd.ic_dreg's group law."""
        E = self.E
        fb = [(int(x), int(y)) for x, y in zip(self.fx, self.fy)]
        out = set()

        def rec(chosen, start, S):
            if len(chosen) == m - 1:
                if S is not None and _apply(self.chk, S[0]) == 0:
                    out.add(tuple(sorted(chosen + [S])))
                return
            for i in range(start, len(fb)):
                rec(chosen + [fb[i]], i, E.add(S, E.neg(fb[i])))
        rec([], 0, R)
        return sorted(out)
