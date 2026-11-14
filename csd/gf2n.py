"""GF(2^n) for n <= 32 as Python ints and numpy uint64 arrays, F_2[x]
factoring, and the F_2-subspaces of GF(2^n) the index-calculus sweep uses
as factor bases.

Elements are bit vectors in the polynomial basis modulo a fixed
irreducible (the lowest-weight trinomial, else pentanomial, chosen
deterministically). Squaring sigma(x) = x^2 is F_2-linear; a subspace V is
sigma-stable iff it is the kernel of L_g(x) = sum g_i x^(2^i) for a divisor
g of x^n - 1, which is also iff its subspace polynomial prod_{v in V}(X - v)
has coefficients in F_2. For prime n the stable subspaces are exactly the
sums of the kernels of the irreducible factors of x^n - 1 over F_2: one of
degree 1 and (n - 1) / ord_n(2) of degree ord_n(2).
"""
from __future__ import annotations

import itertools
import random

import numpy as np


# ------------------------------------------------------------------ F_2[x]
def pmul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def pdivmod(a: int, b: int) -> tuple[int, int]:
    assert b
    q, db = 0, b.bit_length()
    while a.bit_length() >= db:
        s = a.bit_length() - db
        q ^= 1 << s
        a ^= b << s
    return q, a


def pmod(a: int, b: int) -> int:
    return pdivmod(a, b)[1]


def pgcd(a: int, b: int) -> int:
    while b:
        a, b = b, pmod(a, b)
    return a


def pdeg(a: int) -> int:
    return a.bit_length() - 1


def ppowmod(a: int, e: int, m: int) -> int:
    r, a = 1, pmod(a, m)
    while e:
        if e & 1:
            r = pmod(pmul(r, a), m)
        a = pmod(pmul(a, a), m)
        e >>= 1
    return r


def is_irreducible(f: int) -> bool:
    """Ben-Or: f of degree n is irreducible iff gcd(x^(2^i) - x, f) = 1 for i <= n/2."""
    n = pdeg(f)
    if n <= 0:
        return False
    x = 2
    h = x
    for _ in range(n // 2):
        h = pmod(pmul(h, h), f)
        if pgcd(h ^ x, f) != 1:
            return False
    return True


def default_modulus(n: int) -> int:
    for k in range(1, n):
        f = (1 << n) | (1 << k) | 1
        if is_irreducible(f):
            return f
    for k1, k2, k3 in itertools.combinations(range(1, n), 3):
        f = (1 << n) | (1 << k3) | (1 << k2) | (1 << k1) | 1
        if is_irreducible(f):
            return f
    raise ValueError(n)


def _f2_nullspace(rows: list[int], ncols: int) -> list[int]:
    """Null space {v : M v = 0} of the matrix whose i-th row is the bit vector rows[i]."""
    piv: dict[int, int] = {}
    for r in rows:
        for c, pr in piv.items():
            if (r >> c) & 1:
                r ^= pr
        if r:
            c = r.bit_length() - 1
            for c2 in list(piv):
                if (piv[c2] >> c) & 1:
                    piv[c2] ^= r
            piv[c] = r
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for fc in free:
        v = 1 << fc
        for c, pr in piv.items():
            if (pr >> fc) & 1:
                v |= 1 << c
        out.append(v)
    return out


def factor_squarefree(f: int) -> list[int]:
    """Berlekamp over F_2 for a squarefree f; returns the irreducible factors."""
    n = pdeg(f)
    # Q - I as columns: column i is x^(2i) mod f - x^i.
    cols = [pmod(ppowmod(2, 2 * i, f), f) ^ (1 << i) for i in range(n)]
    rows = [sum(((cols[i] >> r) & 1) << i for i in range(n)) for r in range(n)]
    basis = _f2_nullspace(rows, n)
    k = len(basis)
    factors = [f]
    for v in basis:
        if len(factors) == k:
            break
        new = []
        for g in factors:
            if pdeg(g) == 1:
                new.append(g)
                continue
            a = pgcd(g, v)
            if 0 < pdeg(a) < pdeg(g):
                new += [a, pdivmod(g, a)[0]]
            else:
                new.append(g)
        factors = new
    assert len(factors) == k, (f, factors)
    return sorted(factors)


def multiplicative_order_2(n: int) -> int:
    return next(k for k in range(1, n + 1) if pow(2, k, n) == 1)


# ------------------------------------------------------------------ GF(2^n)
class GF2n:
    def __init__(self, n: int, modulus: int | None = None):
        assert 2 <= n <= 32
        self.n = n
        self.mod = modulus or default_modulus(n)
        assert self.mod >> n == 1 and is_irreducible(self.mod)
        self.order = 1 << n
        self.mask = self.order - 1
        # x^n, x^(n+1), ..., x^(2n-2) reduced: reduction table for numpy.
        self._red = [pmod(1 << (n + i), self.mod) for i in range(n - 1)]
        self.tr_mask = sum(1 << i for i in range(n) if self.trace(1 << i))

    # scalar
    def mul(self, a: int, b: int) -> int:
        return pmod(pmul(a, b), self.mod)

    def sqr(self, a: int) -> int:
        return self.mul(a, a)

    def pow(self, a: int, e: int) -> int:
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def inv(self, a: int) -> int:
        assert a
        return self.pow(a, self.order - 2)

    def frob(self, a: int, k: int = 1) -> int:
        for _ in range(k % self.n):
            a = self.sqr(a)
        return a

    def sqrt(self, a: int) -> int:
        return self.frob(a, self.n - 1)

    def trace(self, a: int) -> int:
        t, x = 0, a
        for _ in range(self.n):
            t ^= x
            x = self.sqr(x)
        assert t in (0, 1)
        return t

    def half_trace(self, c: int) -> int:
        assert self.n % 2 == 1
        h, x = 0, c
        for _ in range((self.n + 1) // 2):
            h ^= x
            x = self.sqr(self.sqr(x))
        return h

    def solve_z2_z(self, c: int) -> list[int]:
        """Roots of z^2 + z = c in GF(2^n), n odd."""
        if bin(c & self.tr_mask).count("1") & 1:
            return []
        z = self.half_trace(c)
        return [z, z ^ 1]

    # vectorised (uint64 arrays; n <= 32 keeps the product below 2^63)
    def vmul(self, a, b):
        a = np.asarray(a, dtype=np.uint64)
        b = np.asarray(b, dtype=np.uint64)
        a, b = np.broadcast_arrays(a, b)
        r = np.zeros(a.shape, dtype=np.uint64)
        one = np.uint64(1)
        for i in range(self.n):
            bit = (b >> np.uint64(i)) & one
            r ^= (a << np.uint64(i)) * bit
        return self._vreduce(r)

    def _vreduce(self, r):
        n = self.n
        one = np.uint64(1)
        for i in range(n - 2, -1, -1):
            bit = (r >> np.uint64(n + i)) & one
            r ^= np.uint64(self._red[i]) * bit
        return r & np.uint64(self.mask)

    def vsqr(self, a):
        return self.vmul(a, a)

    def vpow(self, a, e: int):
        a = np.asarray(a, dtype=np.uint64)
        r = np.ones_like(a)
        while e:
            if e & 1:
                r = self.vmul(r, a)
            a = self.vsqr(a)
            e >>= 1
        return r

    def vinv(self, a):
        return self.vpow(a, self.order - 2)

    def vtrace(self, a):
        x = np.asarray(a, dtype=np.uint64) & np.uint64(self.tr_mask)
        t = np.zeros_like(x)
        for i in range(self.n):
            t ^= (x >> np.uint64(i)) & np.uint64(1)
        return t

    def vhalf_trace(self, c):
        x = np.asarray(c, dtype=np.uint64)
        h = np.zeros_like(x)
        for _ in range((self.n + 1) // 2):
            h ^= x
            x = self.vsqr(self.vsqr(x))
        return h

    # F_2-linear maps on GF(2^n) as lists of images of the basis 1, z, ..., z^(n-1)
    def linear_map_kernel(self, images: list[int]) -> list[int]:
        rows = [sum(((images[i] >> r) & 1) << i for i in range(self.n)) for r in range(self.n)]
        return _f2_nullspace(rows, self.n)

    def linearized(self, g: int, x: int) -> int:
        """L_g(x) = sum g_i x^(2^i)."""
        r, y, i = 0, x, 0
        while g >> i:
            if (g >> i) & 1:
                r ^= y
            y = self.sqr(y)
            i += 1
        return r


# ------------------------------------------------------------------ subspaces
def echelon(vectors: list[int]) -> dict[int, int]:
    piv: dict[int, int] = {}
    for v in vectors:
        for c in sorted(piv, reverse=True):
            if (v >> c) & 1:
                v ^= piv[c]
        if v:
            piv[v.bit_length() - 1] = v
    return piv


def f2_rank(vectors: list[int]) -> int:
    return len(echelon(vectors))


def span_elements(basis: list[int]) -> np.ndarray:
    """All 2^k elements of span(basis), index i = sum bit_j(i) * basis_j."""
    out = np.zeros(1 << len(basis), dtype=np.uint64)
    for j, b in enumerate(basis):
        h = 1 << j
        out[h:2 * h] = out[:h] ^ np.uint64(b)
    return out


def sigma_stable_subspaces(F: GF2n) -> list[dict]:
    """Every sigma-stable subspace: kernels of L_g for each divisor g of x^n - 1."""
    n = F.n
    facs = factor_squarefree((1 << n) | 1)
    out = []
    for r in range(1, len(facs) + 1):
        for sub in itertools.combinations(range(len(facs)), r):
            g = 1
            for i in sub:
                g = pmul(g, facs[i])
            images = [F.linearized(g, 1 << i) for i in range(n)]
            basis = F.linear_map_kernel(images)
            assert len(basis) == pdeg(g)
            out.append({"g": g, "factors": [facs[i] for i in sub], "dim": len(basis), "basis": basis})
    return out


def is_sigma_stable(F: GF2n, basis: list[int]) -> bool:
    k = f2_rank(basis)
    return f2_rank(basis + [F.sqr(b) for b in basis]) == k


def random_subspace(F: GF2n, k: int, rng: random.Random, avoid_stable: bool = True) -> list[int]:
    while True:
        basis = []
        while len(basis) < k:
            v = rng.randrange(1, F.order)
            if f2_rank(basis + [v]) > len(basis):
                basis.append(v)
        if not avoid_stable or not is_sigma_stable(F, basis):
            return basis


def fppr_subspace(F: GF2n, k: int) -> list[int]:
    """Faugere-Perret-Petit-Renault's choice: span{1, z, ..., z^(k-1)} in the polynomial basis."""
    return [1 << i for i in range(k)]
