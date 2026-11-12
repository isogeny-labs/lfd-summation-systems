"""Factor base against degree of regularity: the point-decomposition systems of
index calculus on prime, binary and Koblitz curves, and their measured degrees.

Relation search writes a target R as R = P_1 + ... + P_m with every P_i in a
factor base F = {P : x(P) in V}. Algebraically that is S_{m+1}(x_1, ..., x_m,
x(R)) = 0 with every x_i in V (Semaev's summation polynomial). This module
builds that system for the three curve families the literature attacks this way
and measures how hard it is for Groebner / XL-style solvers, as a function of the
factor base.

* **binary** (Diem; Petit-Quisquater, ePrint 2012/146; Faugere-Perret-Petit-
  Renault): E: y^2 + xy = x^3 + a x^2 + b over F_{2^n}, V an F_2-subspace of
  dimension n' (the factor base dimension). x_i = sum_j c_ij v_j with Boolean
  c_ij; the coefficients of S_{m+1} in an F_2-basis of F_{2^n} give n Boolean
  equations in m n' unknowns, the field equations c^2 = c built in.
* **koblitz** (Galbraith-Granger-Merz-Petit, SAC 2020, ePrint 2020/1315): the
  binary system on y^2 + xy = x^3 + a x^2 + 1, a in {0, 1}. A factor base is
  Frobenius invariant iff V is: V is then a sum of the subspaces V_f attached to
  the irreducible factors f of x^n - 1 over F_2 (their Lemma 4.1), so its
  dimension is k ell + eps with ell = ord_n(2). `invariant_subspaces` builds
  them; `koblitz_invariant_dims` lists the dimensions a given n allows.
* **prime** (Petit-Kosters-Messeng, PKC 2016, section 3.3): E: y^2 = x^3 + ax + b
  over F_p with 2^{n1} | p - 1 and V the subgroup of order 2^{n1} of F_p^*, cut
  out by the chain x_{i,k+1} = x_{i,k}^2 (k < n1 - 1), x_{i,n1-1}^2 = 1. The
  chain is a Groebner basis by itself (pure-power leading terms in distinct
  variables), so every row of the Macaulay matrix reduces into the algebra
  F_p[x_1, ..., x_m]/(x_i^{2^{n1}} - 1), in which x^E has chain degree
  popcount(E): the same squarefree-monomial bookkeeping as the Boolean case,
  with carries. The Macaulay matrix at degree D of {S_3} + chain is the matrix
  of the reduced multiples NF(mu * S_3), mu squarefree of degree <= D - 4.

Measured degrees, per (curve, V, target):

* ``d_ff`` (binary / Koblitz): the first fall degree, the smallest D at which
  the top-degree forms of the system have a syzygy beyond the generic
  (semi-regular) count in the graded Boolean ring F_2[x]/(x_i^2): rank of the
  degree-D part of the top-form ideal below C(v, D) - [HS]_D, with
  HS = (1 + t)^v / prod(1 + t^{d_i}) (Bardet-Faugere-Salvy). Kosters-Yeo
  (arXiv 1503.08001, Cor. 4.11) prove it is 2 for m = 2 by a trace identity.
* ``d_reg`` (binary / Koblitz): the degree of regularity by its definition, the
  smallest D at which the top-form ideal contains every monomial of degree D.
* ``d_solve``: the solving degree of XL on the Macaulay matrix, the smallest D
  whose row space contains every affine-linear polynomial that vanishes on the
  solution set. The target dimension is computed from the solution set found
  by an independent enumeration (not from the system), so a run that never
  reaches it, or exceeds it, is detected. With at most one solution this is the
  degree at which the matrix contains the reduced Groebner basis.
* ``d_reg_sr`` (formula): the semi-regular degree of regularity of a system with
  the same number of equations, degrees and unknowns.

Also measured: |F| (exact count over V), the solution set and the genuine
decompositions of R (lifted to points and checked with the group law), and the
counted operations of the incremental XL run (kernel row operations, `elim.py`).
Nothing here is wall-clock.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import combinations

import numpy as np

from . import elim
from .gf2_ghs import GF2n, clmul, poly_mod
from .field import factor, is_prime, legendre, sqrt_mod

# ---------------------------------------------------------------------------
# F_2[x] helpers and moduli without PARI
# ---------------------------------------------------------------------------


def _f2_mulmod(a: int, b: int, f: int) -> int:
    return poly_mod(clmul(a, b), f)


def _f2_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


def f2_irreducible(f: int) -> bool:
    """Rabin's test over F_2: x^(2^N) = x mod f and gcd(x^(2^(N/q)) - x, f) = 1 for
    every prime q | N."""
    N = f.bit_length() - 1
    if N < 1:
        return False
    if N == 1:
        return True
    xp = [2]  # x^(2^k) mod f
    for _ in range(N):
        xp.append(_f2_mulmod(xp[-1], xp[-1], f))
    if xp[N] != 2 % f:
        return False
    for q in factor(N):
        if _f2_gcd(f, xp[N // q] ^ 2) != 1:
            return False
    return True


@lru_cache(maxsize=None)
def binary_modulus(n: int) -> int:
    """The lowest-weight irreducible: x^n + x^k + 1 with least k, else the least
    pentanomial x^n + x^k3 + x^k2 + x^k1 + 1."""
    for k in range(1, n):
        f = (1 << n) | (1 << k) | 1
        if f2_irreducible(f):
            return f
    for k3 in range(3, n):
        for k2 in range(2, k3):
            for k1 in range(1, k2):
                f = (1 << n) | (1 << k3) | (1 << k2) | (1 << k1) | 1
                if f2_irreducible(f):
                    return f
    raise ValueError(n)


def binary_field(n: int) -> GF2n:
    return GF2n(n, binary_modulus(n))


def _f2_solve(images: list[int], rhs: int, nbits: int) -> tuple[int, list[int]] | None:
    """All c in F_2^k with sum_j c_j images[j] = rhs (vectors of `nbits` bits):
    (a particular solution, a kernel basis), or None. c is returned as a k-bit int."""
    k = len(images)
    # columns = unknowns; eliminate on the augmented rows (one row per bit)
    rows = []
    for b in range(nbits):
        r = 0
        for j, im in enumerate(images):
            if im >> b & 1:
                r |= 1 << j
        rows.append((r, rhs >> b & 1))
    piv: list[tuple[int, int, int]] = []  # (pivot bit, row, rhs)
    for r, t in rows:
        for pb, pr, pt in piv:
            if r >> pb & 1:
                r ^= pr
                t ^= pt
        if r:
            pb = r.bit_length() - 1
            # keep the basis reduced on pb
            piv = [(qb, qr ^ r, qt ^ t) if qr >> pb & 1 else (qb, qr, qt) for qb, qr, qt in piv]
            piv.append((pb, r, t))
        elif t:
            return None
    pivbits = {pb for pb, _, _ in piv}
    part = 0
    for pb, pr, pt in piv:
        if pt:
            part |= 1 << pb
    kernel = []
    for f in range(k):
        if f in pivbits:
            continue
        v = 1 << f
        for pb, pr, pt in piv:
            if pr >> f & 1:
                v |= 1 << pb
        kernel.append(v)
    return part, kernel


# ---------------------------------------------------------------------------
# binary curves y^2 + xy = x^3 + a x^2 + b
# ---------------------------------------------------------------------------


class BinaryCurve:
    def __init__(self, F: GF2n, a: int, b: int):
        assert b != 0
        self.F, self.a, self.b = F, a, b
        n = F.N
        self._qimg = [F.sq(1 << i) ^ (1 << i) for i in range(n)]  # z -> z^2 + z

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        F = self.F
        x, y = P
        return F.mul(y, y) ^ F.mul(x, y) == F.mul(F.mul(x, x), x) ^ F.mul(self.a, F.mul(x, x)) ^ self.b

    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])

    def add(self, P, Q):
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if y2 == x1 ^ y1:  # Q = -P (also P = -P, the 2-torsion point x = 0)
                return None
            lam = x1 ^ F.mul(y1, F.inv(x1))  # doubling
            x3 = F.mul(lam, lam) ^ lam ^ self.a
            y3 = F.mul(x1, x1) ^ F.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.mul(lam, lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = F.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def mul(self, k: int, P):
        R, Q = None, P
        if k < 0:
            k, Q = -k, self.neg(P)
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def lift_x(self, x: int):
        """A point with abscissa x, or None if x is not the abscissa of an F_{2^n}-point."""
        F = self.F
        if x == 0:
            return (0, F.sqrt(self.b))
        x2 = F.mul(x, x)
        c = x ^ self.a ^ F.mul(self.b, F.inv(x2))  # y = x z, z^2 + z = c
        if F.trace(c):
            return None
        sol = _f2_solve(self._qimg, c, F.N)
        assert sol is not None
        return (x, F.mul(x, sol[0]))

    def random_point(self, rng: random.Random):
        while True:
            P = self.lift_x(rng.getrandbits(self.F.N))
            if P is not None:
                return P if rng.getrandbits(1) else self.neg(P)

    @property
    def q(self) -> int:
        return 1 << self.F.N


def koblitz_order(n: int, a: int) -> int:
    """#E(F_{2^n}) for y^2 + xy = x^3 + a x^2 + 1 by the Lucas sequence of
    T^2 - t T + 2, t = 2a - 1 (#E(F_2) = 4 for a = 0, 2 for a = 1)."""
    t = 2 * a - 1
    v0, v1 = 2, t
    for _ in range(n - 1):
        v0, v1 = v1, t * v1 - 2 * v0
    return (1 << n) + 1 - v1


# ---------------------------------------------------------------------------
# prime curves y^2 = x^3 + a x + b
# ---------------------------------------------------------------------------


class PrimeCurve:
    def __init__(self, p: int, a: int, b: int):
        assert (4 * a ** 3 + 27 * b * b) % p
        self.p, self.a, self.b = p, a % p, b % p

    def rhs(self, x: int) -> int:
        return (x * x * x + self.a * x + self.b) % self.p

    def on_curve(self, P) -> bool:
        return P is None or (P[1] * P[1] - self.rhs(P[0])) % self.p == 0

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, k: int, P):
        R, Q = None, P
        if k < 0:
            k, Q = -k, self.neg(P)
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def lift_x(self, x: int):
        r = self.rhs(x)
        if r == 0:
            return (x % self.p, 0)
        if legendre(r, self.p) != 1:
            return None
        return (x % self.p, sqrt_mod(r, self.p))

    def random_point(self, rng: random.Random):
        while True:
            P = self.lift_x(rng.randrange(self.p))
            if P is not None and P[1]:
                return P if rng.getrandbits(1) else self.neg(P)

    @property
    def q(self) -> int:
        return self.p


def group_order(E, rng: random.Random, tries: int = 24) -> int:
    """#E(F_q): the orders of random points (baby-step giant-step in the Hasse
    interval) until a single multiple of their lcm is left in the interval.
    Raises ValueError if the order stays ambiguous (exponent below 4 sqrt q)."""
    q = E.q
    w = math.isqrt(4 * q) + 1
    lo, hi = q + 1 - w, q + 1 + w
    s = math.isqrt(2 * w) + 1
    ex = 1
    cands: list[int] = []
    for _ in range(tries):
        P = E.random_point(rng)
        baby: dict = {}
        T = None
        for j in range(s):
            baby.setdefault(T, j)
            T = E.add(T, P)
        step, G, N = T, E.mul(lo, P), None   # T = s P
        for i in range(s + 1):
            j = baby.get(E.neg(G))           # (lo + i s) P + j P = O
            if j is not None:
                N = lo + i * s + j
                break
            G = E.add(G, step)
        assert N is not None and E.mul(N, P) is None
        o = _point_order(E, P, N)
        ex = ex * o // math.gcd(ex, o)
        cands = list(range(-(-lo // ex) * ex, hi + 1, ex))
        if len(cands) == 1:
            return cands[0]
    raise ValueError(f"order ambiguous: exponent {ex}, candidates {cands[:5]}")


def _point_order(E, P, N: int) -> int:
    """Order of P given N P = O."""
    o = N
    for ell in factor(N):
        while o % ell == 0 and E.mul(o // ell, P) is None:
            o //= ell
    return o


def largest_prime_factor(N: int) -> int:
    return max(factor(N))


# ---------------------------------------------------------------------------
# summation polynomials
# ---------------------------------------------------------------------------


def s3_binary_value(F: GF2n, b: int, x1: int, x2: int, x3: int) -> int:
    """S_3 for y^2 + xy = x^3 + a x^2 + b: (x1 x2 + x1 x3 + x2 x3)^2 + x1 x2 x3 + b."""
    t = F.mul(x1, x2) ^ F.mul(x1, x3) ^ F.mul(x2, x3)
    return F.mul(t, t) ^ F.mul(F.mul(x1, x2), x3) ^ b


def s4_binary_value(F: GF2n, b: int, x1: int, x2: int, x3: int, x4: int) -> int:
    """S_4 = Res_X(S_3(x1, x2, X), S_3(x3, x4, X)) with S_3(u, v, X) =
    A X^2 + B X + C, A = (u + v)^2, B = u v, C = (u v)^2 + b:
    Res = (A1 C2 + A2 C1)^2 + (A1 B2 + A2 B1)(B1 C2 + B2 C1) in characteristic 2."""
    m = F.mul
    A1, B1 = m(x1 ^ x2, x1 ^ x2), m(x1, x2)
    C1 = m(B1, B1) ^ b
    A2, B2 = m(x3 ^ x4, x3 ^ x4), m(x3, x4)
    C2 = m(B2, B2) ^ b
    u = m(A1, C2) ^ m(A2, C1)
    return m(u, u) ^ m(m(A1, B2) ^ m(A2, B1), m(B1, C2) ^ m(B2, C1))


def s3_prime_coeffs(E: PrimeCurve, X: int) -> dict[tuple[int, int], int]:
    """S_3(x1, x2, X) for y^2 = x^3 + a x + b as {(e1, e2): coefficient mod p}:
    (x1 - x2)^2 X^2 - 2((x1 + x2)(x1 x2 + a) + 2b) X + (x1 x2 - a)^2 - 4b(x1 + x2)."""
    p, a, b = E.p, E.a, E.b
    c = {
        (2, 0): X * X, (0, 2): X * X, (1, 1): -2 * X * X - 2 * a,
        (2, 1): -2 * X, (1, 2): -2 * X,
        (1, 0): -2 * a * X - 4 * b, (0, 1): -2 * a * X - 4 * b,
        (0, 0): -4 * b * X + a * a, (2, 2): 1,
    }
    return {e: v % p for e, v in c.items() if v % p}


def s3_prime_value(E: PrimeCurve, x1: int, x2: int, X: int) -> int:
    return sum(c * pow(x1, e1, E.p) * pow(x2, e2, E.p) for (e1, e2), c in s3_prime_coeffs(E, X).items()) % E.p


# ---------------------------------------------------------------------------
# subspaces of F_{2^n}
# ---------------------------------------------------------------------------


def standard_subspace(F: GF2n, d: int) -> list[int]:
    """<1, w, ..., w^{d-1}> (the 'standard' choice of Semaev / GGMP Table 2)."""
    return [1 << i for i in range(d)]


def random_subspace(F: GF2n, d: int, rng: random.Random) -> list[int]:
    while True:
        B = [rng.getrandbits(F.N) for _ in range(d)]
        if _f2_rank(B) == d:
            return B


def _f2_rank(vs: list[int]) -> int:
    basis: dict[int, int] = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return len(basis)


def _f2_poly_factor_xn1(n: int) -> list[int]:
    """Irreducible factors of x^n - 1 over F_2 for odd n (one per 2-cyclotomic coset)."""
    assert n % 2 == 1
    # x^n - 1 is squarefree for odd n; split it with gcd(g, T_D(r)), T_D the
    # F_{2^D}-trace map r + r^2 + ... + r^(2^(D-1)) (equal-degree splitting)
    f = (1 << n) | 1
    rng = random.Random(0)
    out = []

    def split(g: int):
        d = g.bit_length() - 1
        if d == 0:
            return
        if f2_irreducible(g):
            out.append(g)
            return
        # find the degree of the factors: distinct-degree on g
        while True:
            r = rng.getrandbits(d) | 2
            # trace map T(r) = r + r^2 + r^4 + ... + r^(2^(D-1)) for D = ell (unknown): try all D
            for D in range(1, d + 1):
                t, acc = r, r
                for _ in range(D - 1):
                    t = _f2_mulmod(t, t, g)
                    acc ^= t
                h = _f2_gcd(g, acc)
                if 0 < h.bit_length() - 1 < d:
                    split(h)
                    split(_f2_div(g, h))
                    return

    split(f)
    return sorted(out, key=lambda g: (g.bit_length(), g))


def _f2_div(a: int, b: int) -> int:
    q = 0
    db = b.bit_length() - 1
    while a.bit_length() - 1 >= db:
        s = a.bit_length() - 1 - db
        q |= 1 << s
        a ^= b << s
    assert a == 0
    return q


def multiplicative_order_mod(a: int, n: int) -> int:
    k, x = 1, a % n
    while x != 1:
        x = x * a % n
        k += 1
    return k


def koblitz_invariant_dims(n: int) -> list[int]:
    """Dimensions of the Frobenius-invariant F_2-subspaces of F_{2^n}, n odd:
    sums of the degrees of distinct irreducible factors of x^n - 1 (GGMP Lemma 4.1)."""
    degs = [g.bit_length() - 1 for g in _f2_poly_factor_xn1(n)]
    sums = {0}
    for d in degs:
        sums |= {s + d for s in sums}
    return sorted(sums)


def invariant_subspaces(F: GF2n, d: int) -> list[list[int]]:
    """Bases of every Frobenius-invariant subspace of dimension d. V_f = ker f(sigma)
    for f | x^n - 1, sigma the squaring map; a sum over a set of factors."""
    n = F.N
    facs = _f2_poly_factor_xn1(n)
    out = []
    for k in range(len(facs) + 1):
        for S in combinations(range(len(facs)), k):
            if sum(facs[i].bit_length() - 1 for i in S) != d:
                continue
            g = (1 << n) | 1
            for i in S:
                g = _f2_div(g, facs[i])
            # V = image of g(sigma) = kernel of prod_{i in S} f_i(sigma)
            imgs = [_apply_sigma_poly(F, g, 1 << j) for j in range(n)]
            basis = _f2_basis(imgs)
            assert len(basis) == d, (len(basis), d)
            out.append(basis)
    return out


def _apply_sigma_poly(F: GF2n, g: int, x: int) -> int:
    acc, cur, k = 0, x, 0
    while g >> k:
        if g >> k & 1:
            acc ^= cur
        cur = F.sq(cur)
        k += 1
    return acc


def _f2_basis(vs: list[int]) -> list[int]:
    basis: dict[int, int] = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return sorted(basis.values())


def is_frobenius_invariant(F: GF2n, basis: list[int]) -> bool:
    r = _f2_rank(basis)
    return all(_f2_rank(basis + [F.sq(v)]) == r for v in basis)


def span(basis: list[int]) -> list[int]:
    out = [0]
    for v in basis:
        out += [u ^ v for u in out]
    return out


# ---------------------------------------------------------------------------
# Boolean polynomials with F_{2^n} coefficients, and Weil descent
# ---------------------------------------------------------------------------

BPoly = dict  # {variable mask: F_{2^n} element}


def _bp_add(f: BPoly, g: BPoly) -> BPoly:
    out = dict(f)
    for k, c in g.items():
        v = out.get(k, 0) ^ c
        if v:
            out[k] = v
        else:
            out.pop(k, None)
    return out


def _bp_mul(F: GF2n, f: BPoly, g: BPoly) -> BPoly:
    out: BPoly = {}
    for k1, c1 in f.items():
        for k2, c2 in g.items():
            k = k1 | k2
            v = out.get(k, 0) ^ F.mul(c1, c2)
            if v:
                out[k] = v
            else:
                out.pop(k, None)
    return out


def _bp_sq(F: GF2n, f: BPoly) -> BPoly:
    # (sum c_M M)^2 = sum c_M^2 M in characteristic 2 with M^2 = M
    return {k: F.sq(c) for k, c in f.items()}


def _bp_const(c: int) -> BPoly:
    return {0: c} if c else {}


def _bp_block(basis: list[int], offset: int) -> BPoly:
    return {1 << (offset + j): v for j, v in enumerate(basis)}


def binary_descent(E: BinaryCurve, bases: list[list[int]], xR: int) -> list[list[int]]:
    """Weil descent of S_{m+1}(x_1, ..., x_m, xR), x_i in span(bases[i]), m = 2 or 3.
    Returns the n Boolean equations as sorted lists of monomial masks (0 = constant)."""
    F, b = E.F, E.b
    m = len(bases)
    offs = np.cumsum([0] + [len(B) for B in bases]).tolist()
    xs = [_bp_block(B, offs[i]) for i, B in enumerate(bases)]
    R = _bp_const(xR)
    if m == 2:
        x1, x2 = xs
        t = _bp_add(_bp_add(_bp_mul(F, x1, x2), _bp_mul(F, x1, R)), _bp_mul(F, x2, R))
        S = _bp_add(_bp_add(_bp_sq(F, t), _bp_mul(F, _bp_mul(F, x1, x2), R)), _bp_const(b))
    elif m == 3:
        x1, x2, x3 = xs
        s12 = _bp_add(x1, x2)
        A1, B1 = _bp_sq(F, s12), _bp_mul(F, x1, x2)
        C1 = _bp_add(_bp_sq(F, B1), _bp_const(b))
        A2, B2 = _bp_sq(F, _bp_add(x3, R)), _bp_mul(F, x3, R)
        C2 = _bp_add(_bp_sq(F, B2), _bp_const(b))
        u = _bp_add(_bp_mul(F, A1, C2), _bp_mul(F, A2, C1))
        S = _bp_add(_bp_sq(F, u), _bp_mul(F, _bp_add(_bp_mul(F, A1, B2), _bp_mul(F, A2, B1)),
                                          _bp_add(_bp_mul(F, B1, C2), _bp_mul(F, B2, C1))))
    else:
        raise ValueError(m)
    eqs = []
    for k in range(F.N):
        eq = sorted(mask for mask, c in S.items() if c >> k & 1)
        if eq:
            eqs.append(eq)
    return eqs


# ---------------------------------------------------------------------------
# independent enumeration of the solution set
# ---------------------------------------------------------------------------


def binary_solutions(E: BinaryCurve, bases: list[list[int]], xR: int) -> list[tuple[int, ...]]:
    """Every (c_1, ..., c_m) (coordinate vectors, ints of n' bits) with
    S_{m+1}(x_1, ..., x_m, xR) = 0, x_i = sum_j c_ij v_ij, found in F_{2^n}
    without the descended system. m = 2: for each x1, S_3 is
    (x1 + xR)^2 x2^2 + x1 xR x2 + (x1 xR)^2 + b, F_2-affine in x2, solved on V.
    m = 3: for each (x1, x2), brute force over x3."""
    F, b = E.F, E.b
    m = len(bases)
    if m == 2:
        B1, B2 = bases
        out = []
        V1 = span(B1)
        # coordinates of span elements: span() lists them in binary-counter order
        for c1, x1 in enumerate(V1):
            al = F.mul(x1 ^ xR, x1 ^ xR)
            be = F.mul(x1, xR)
            ga = F.mul(be, be) ^ b
            imgs = [F.mul(al, F.mul(v, v)) ^ F.mul(be, v) for v in B2]
            sol = _f2_solve(imgs, ga, F.N)
            if sol is None:
                continue
            part, ker = sol
            for c2 in _affine_span(part, ker):
                out.append((c1, c2))
        return out
    if m == 3:
        B1, B2, B3 = bases
        V1, V2, V3 = span(B1), span(B2), span(B3)
        out = []
        for c1, x1 in enumerate(V1):
            for c2, x2 in enumerate(V2):
                for c3, x3 in enumerate(V3):
                    if s4_binary_value(F, b, x1, x2, x3, xR) == 0:
                        out.append((c1, c2, c3))
        return out
    raise ValueError(m)


def _affine_span(part: int, ker: list[int]) -> list[int]:
    out = [part]
    for k in ker:
        out += [u ^ k for u in out]
    return out


def binary_decompositions(E: BinaryCurve, bases: list[list[int]], R, sols) -> int:
    """Number of distinct multisets {x_1, ..., x_m} among the solutions for which
    some choice of points P_i with x(P_i) = x_i sums to R exactly (group law)."""
    Vs = [span(B) for B in bases]
    seen = set()
    for c in sols:
        xs = tuple(sorted(Vs[i][ci] for i, ci in enumerate(c)))
        if xs in seen:
            continue
        pts = [E.lift_x(x) for x in xs]
        if any(P is None for P in pts):
            continue
        if _signs_hit(E, pts, R):
            seen.add(xs)
    return len(seen)


def _signs_hit(E, pts, R) -> bool:
    for s in range(1 << len(pts)):
        acc = None
        for i, P in enumerate(pts):
            acc = E.add(acc, E.neg(P) if s >> i & 1 else P)
        if acc == R:
            return True
    return False


def binary_factor_base_size(E: BinaryCurve, basis: list[int]) -> int:
    F = E.F
    n = 0
    for x in span(basis):
        if x == 0:
            n += 1
        elif F.trace(x ^ E.a ^ F.mul(E.b, F.inv(F.mul(x, x)))) == 0:
            n += 2
    return n


def subgroup_2power(p: int, n1: int) -> list[int]:
    """The subgroup of order 2^{n1} of F_p^* (2^{n1} | p - 1), listed as g^E, E = 0..2^{n1}-1,
    so that the chain coordinates of g^E are the bits of E."""
    assert (p - 1) % (1 << n1) == 0
    N = 1 << n1
    for h in range(2, p):
        g = pow(h, (p - 1) // N, p)
        if pow(g, N // 2, p) != 1:
            return [pow(g, E, p) for E in range(N)]
    raise ValueError(p)


def prime_solutions(E: PrimeCurve, V: list[int], xR: int) -> list[tuple[int, int]]:
    """Every (E1, E2) with S_3(V[E1], V[E2], xR) = 0, by evaluation over V x V."""
    p = E.p
    out = []
    coeffs = s3_prime_coeffs(E, xR)
    pw = {e: [pow(x, e, p) for x in V] for e in (0, 1, 2)}
    for i, x1 in enumerate(V):
        for j, x2 in enumerate(V):
            s = 0
            for (e1, e2), c in coeffs.items():
                s += c * pw[e1][i] * pw[e2][j]
            if s % p == 0:
                out.append((i, j))
    return out


def prime_decompositions(E: PrimeCurve, V: list[int], R, sols) -> int:
    seen = set()
    for i, j in sols:
        xs = tuple(sorted((V[i], V[j])))
        if xs in seen:
            continue
        pts = [E.lift_x(x) for x in xs]
        if any(P is None for P in pts):
            continue
        if _signs_hit(E, pts, R):
            seen.add(xs)
    return len(seen)


def prime_factor_base_size(E: PrimeCurve, V: list[int]) -> int:
    n = 0
    for x in V:
        r = E.rhs(x)
        n += 1 if r == 0 else (2 if legendre(r, E.p) == 1 else 0)
    return n


# ---------------------------------------------------------------------------
# Hilbert series formulas
# ---------------------------------------------------------------------------


def boolean_hilbert(v: int, degrees: list[int], dmax: int) -> list[int]:
    """Coefficients 0..dmax of (1 + t)^v / prod (1 + t^{d_i}) (not truncated)."""
    s = [math.comb(v, k) for k in range(dmax + 1)]
    for d in degrees:
        # divide by (1 + t^d): c_k <- c_k - c_{k-d}
        for k in range(d, dmax + 1):
            s[k] -= s[k - d]
    return s


def boolean_dreg(v: int, degrees: list[int]) -> int:
    """Semi-regular degree of regularity over F_2: index of the first non-positive
    coefficient of (1 + t)^v / prod (1 + t^{d_i})."""
    dmax = v + 1
    s = boolean_hilbert(v, degrees, dmax)
    for k, c in enumerate(s):
        if c <= 0:
            return k
    return dmax


def affine_dreg(v: int, degrees: list[int]) -> int:
    """Semi-regular degree of regularity of an affine system over a large field:
    first non-positive coefficient of prod (1 - t^{d_i}) / (1 - t)^v."""
    dmax = sum(d - 1 for d in degrees) + 2
    s = [math.comb(v - 1 + k, k) for k in range(dmax + 1)]
    for d in degrees:
        for k in range(dmax, d - 1, -1):
            s[k] -= s[k - d]
    for k, c in enumerate(s):
        if c <= 0:
            return k
    return dmax


# ---------------------------------------------------------------------------
# Macaulay matrices in squarefree-monomial algebras
# ---------------------------------------------------------------------------


@lru_cache(maxsize=64)
def _masks_upto(v: int, D: int) -> np.ndarray:
    """All masks of popcount <= D in v bits, sorted by (popcount desc, mask asc)."""
    out = []
    for k in range(D, -1, -1):
        out.extend(sum(1 << i for i in c) for c in combinations(range(v), k))
    return np.asarray(out, dtype=np.uint64)


@lru_cache(maxsize=64)
def _col_index(v: int, D: int):
    cols = _masks_upto(v, D)
    order = np.argsort(cols, kind="stable")
    return cols, cols[order], order


def _cols_of(v: int, D: int, masks: np.ndarray) -> np.ndarray:
    cols, sorted_masks, order = _col_index(v, D)
    pos = np.searchsorted(sorted_masks, masks)
    assert np.all(sorted_masks[np.minimum(pos, len(sorted_masks) - 1)] == masks), "monomial out of range"
    return order[pos]


def n_cols(v: int, D: int) -> int:
    return sum(math.comb(v, k) for k in range(D + 1))


def _popcount(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.uint64)
    c = np.zeros(a.shape, dtype=np.int64)
    while np.any(a):
        c += (a & np.uint64(1)).astype(np.int64)
        a = a >> np.uint64(1)
    return c


def boolean_macaulay(eqs: list[list[int]], v: int, D: int, chunk: int = 1 << 22) -> tuple[np.ndarray, int]:
    """Bit-packed Macaulay matrix at degree D of Boolean equations (rows mu * f,
    mu squarefree with deg mu + deg f <= D, field equations built in). Rows are
    written equation by equation in chunks of about `chunk` entries."""
    cols = n_cols(v, D)
    plan = []
    for f in eqs:
        fa = np.asarray(f, dtype=np.uint64)
        df = int(_popcount(fa).max()) if len(f) else 0
        if df <= D:
            plan.append((fa, _masks_upto(v, D - df)))
    nrows = sum(len(mus) for _, mus in plan)
    words = max(1, (cols + 63) // 64)
    M = np.zeros((nrows, words), dtype=np.uint64)
    row = 0
    one = np.uint64(1)
    for fa, mus in plan:
        step = max(1, chunk // max(1, len(fa)))
        for s0 in range(0, len(mus), step):
            mu = mus[s0:s0 + step]
            prod = (mu[:, None] | fa[None, :]).ravel()
            ri = np.repeat(np.arange(row + s0, row + s0 + len(mu)), len(fa))
            ci = _cols_of(v, D, prod).astype(np.int64)
            np.bitwise_xor.at(M, (ri, ci >> 6), np.left_shift(one, (ci & 63).astype(np.uint64)))
        row += len(mus)
    return M, cols


def boolean_top_rank(eqs: list[list[int]], v: int, D: int) -> int:
    """Rank of the degree-D part of the ideal of the top-degree forms in the graded
    Boolean ring F_2[x]/(x_i^2) (products with a repeated variable vanish)."""
    combs = [sum(1 << i for i in c) for c in combinations(range(v), D)]
    if not combs:
        return 0
    allD = np.asarray(sorted(combs), dtype=np.uint64)
    rows = []
    for f in eqs:
        fa = np.asarray(f, dtype=np.uint64)
        pc = _popcount(fa)
        d = int(pc.max())
        top = fa[pc == d]
        if d > D:
            continue
        mus = np.asarray([sum(1 << i for i in c) for c in combinations(range(v), D - d)], dtype=np.uint64)
        for mu in mus:
            t = top[(top & mu) == 0] | mu
            if t.size:
                rows.append(np.searchsorted(allD, t).tolist())  # distinct: mu is disjoint from t
    if not rows:
        return 0
    M = elim.pack_rows(rows, len(allD))
    r, _, _ = elim.gf2_echelon(M, len(allD))
    return r


def _lin_target_f2(points: list[int], v: int) -> int:
    """Dimension of the affine-linear Boolean polynomials vanishing on the points (v-bit ints)."""
    if not points:
        return v + 1
    rows = [(pt << 1) | 1 for pt in points]
    return v + 1 - _f2_rank(rows)


@dataclass
class DegreeProfile:
    v: int
    n_eq: int
    degrees: list[int]
    d_reg_sr: int
    d_ff: int | None = None
    d_reg: int | None = None
    d_solve: int | None = None
    lin_target: int = 0
    ranks: dict = field(default_factory=dict)       # D -> (rows, cols, rank, lin_dim, top_rank, top_generic)
    ops: dict = field(default_factory=dict)         # D -> counted kernel ops at D (XL)
    d_mxl: int | None = None
    mxl: dict = field(default_factory=dict)         # D -> (rank after closure, lin_dim, rounds, rows added)
    ops_mxl: dict = field(default_factory=dict)     # D -> counted kernel ops at D (MutantXL)
    capped: bool = False

    def xl_ops(self) -> float | None:
        """Counted ops of an XL run stepping D up to d_solve (sum over its degrees)."""
        if self.d_solve is None:
            return None
        return float(sum(o for D, o in self.ops.items() if D <= self.d_solve))

    def mxl_ops(self) -> float | None:
        """Counted ops of a MutantXL run stepping D up to d_mxl."""
        if self.d_mxl is None:
            return None
        return float(sum(o for D, o in self.ops_mxl.items() if D <= self.d_mxl))


# -- MutantXL ------------------------------------------------------------------
#
# MutantXL (Mohamed, Cabarcas, Ding, Buchmann, Bulygin, PQCrypto 2009 / ICISC 2008):
# at degree D, echelonise the Macaulay matrix; every row-space element of degree
# below D (a mutant: its pivot lies in a column of degree < D) is multiplied by
# every variable and the products are added; repeat until the rank stops growing.
# Unlike XL it multiplies the polynomials produced by degree falls, as F4 does.
# Here the scan reports the smallest D at which the closure contains every
# affine-linear polynomial vanishing on the solutions.


def _bool_mulvar(masks: np.ndarray, j: int) -> np.ndarray:
    return masks | np.uint64(1 << j)


def _chain_mulvar(n1: int, m: int):
    N = 1 << n1

    def mul(masks: np.ndarray, j: int) -> np.ndarray:
        blk, k = divmod(j, n1)
        sh = blk * n1
        a = masks.astype(np.int64)
        e = (a >> sh) & (N - 1)
        e = (e + (1 << k)) % N
        return ((a & ~((N - 1) << sh)) | (e << sh)).astype(np.uint64)

    return mul


def _support_gf2(row: np.ndarray, cols: int) -> np.ndarray:
    bits = np.unpackbits(row.view(np.uint8), bitorder="little")[:cols]
    return np.flatnonzero(bits)


def _mxl_closure(M: np.ndarray, cols: int, v: int, D: int, p: int, mulvar, max_rounds: int = 64):
    """Mutant closure at degree D. M is bit-packed (p = 2) or dense mod p. Returns
    (rank, pivots, ops, rounds, rows_added)."""
    masks = _masks_upto(v, D)
    colpop = _popcount(masks)
    first_low = int(np.searchsorted(-colpop, -(D - 1)))  # first column of degree <= D - 1
    ops = 0.0
    rounds = added = 0
    out: list = []
    if p == 2:
        rank, piv, o = elim.gf2_echelon(M, cols, out)
    else:
        rank, piv, o = elim.fp_echelon(M, p, out)
    ops += o
    E = out[-1]
    done: set[bytes] = set()
    while rounds < max_rounds:
        low = [i for i, c in enumerate(piv) if c >= first_low and colpop[c] >= 1]
        todo = []
        for i in low:
            key = E[i].tobytes()
            if key in done:
                continue
            done.add(key)
            todo.append(i)
        nrow = len(todo) * v
        if nrow == 0:
            break
        # products x_j * h written row by row (h in todo order, j inner): the same matrix
        # as building all index lists at once, without their memory
        if p == 2:
            N_ = np.zeros((nrow, E.shape[1]), dtype=np.uint64)
        else:
            N_ = np.zeros((nrow, cols), dtype=np.int64)
        one = np.uint64(1)
        for t, i in enumerate(todo):
            if p == 2:
                sup = _support_gf2(E[i], cols)
            else:
                sup = np.flatnonzero(E[i])
                coef = E[i][sup]
            sm = masks[sup]
            for j in range(v):
                ci = _cols_of(v, D, mulvar(sm, j)).astype(np.int64)
                rr = t * v + j
                if p == 2:
                    np.bitwise_xor.at(N_[rr], ci >> 6, np.left_shift(one, (ci & 63).astype(np.uint64)))
                else:
                    np.add.at(N_[rr], ci, coef)
        if p == 2:
            stack = np.vstack([E, N_])
            del N_
            out = []
            r2, piv2, o = elim.gf2_echelon(stack, cols, out)
        else:
            stack = np.vstack([E, N_ % p])
            del N_
            out = []
            r2, piv2, o = elim.fp_echelon(stack, p, out)
        ops += o
        rounds += 1
        added += nrow
        if r2 == rank:
            break
        rank, piv, E = r2, piv2, out[-1]
    return rank, piv, ops, rounds, added


def _mxl_scan(prof: DegreeProfile, build, v: int, p: int, mulvar, dmin: int, dmax: int, max_cols: int):
    for D in range(dmin, dmax + 1):
        if n_cols(v, D) > max_cols:
            prof.capped = True
            return
        M, cols = build(D)
        rank, piv, ops, rounds, added = _mxl_closure(M, cols, v, D, p, mulvar)
        lin = sum(1 for c in piv if c >= cols - (v + 1))
        assert lin <= prof.lin_target, "mutant closure has a linear form that does not vanish on the solutions"
        prof.ops_mxl[D] = ops
        prof.mxl[D] = (rank, lin, rounds, added)
        if lin == prof.lin_target:
            prof.d_mxl = D
            return


def boolean_profile(eqs: list[list[int]], v: int, points: list[int], dmax: int,
                    max_cols: int = 60000, mutant: bool = True) -> DegreeProfile:
    """Measure d_ff, d_reg and d_solve of a Boolean system; `points` is the solution
    set from the independent enumeration, as v-bit ints."""
    degs = [int(_popcount(np.asarray(f, dtype=np.uint64)).max()) for f in eqs]
    prof = DegreeProfile(v=v, n_eq=len(eqs), degrees=degs, d_reg_sr=boolean_dreg(v, degs))
    prof.lin_target = _lin_target_f2(points, v)
    hs = boolean_hilbert(v, degs, dmax)
    for D in range(min(degs), dmax + 1):
        if n_cols(v, D) > max_cols:
            prof.capped = True
            break
        need_solve = prof.d_solve is None
        if not need_solve:
            break
        need_top = True
        rows = cols = rank = lin = None
        if need_solve:
            M, cols = boolean_macaulay(eqs, v, D)
            rank, piv, ops = elim.gf2_echelon(M, cols)
            rows = M.shape[0]
            lin = sum(1 for c in piv if c >= cols - (v + 1))
            prof.ops[D] = ops
            assert lin <= prof.lin_target, "row space has a linear form that does not vanish on the solutions"
            if lin == prof.lin_target:
                prof.d_solve = D
        top = gen = None
        if need_top:
            top = boolean_top_rank(eqs, v, D)
            full = math.comb(v, D)
            gen = full - hs[D] if hs[D] > 0 else full
            gen = min(gen, full)
            if prof.d_ff is None and top < gen:
                prof.d_ff = D
            if prof.d_reg is None and top == full:
                prof.d_reg = D
        prof.ranks[D] = (rows, cols, rank, lin, top, gen)
    if mutant:
        _mxl_scan(prof, lambda D: boolean_macaulay(eqs, v, D), v, 2, _bool_mulvar, min(degs), dmax, max_cols)
    return prof


# -- PKM chain algebra over F_p ------------------------------------------------


def pkm_generator(E: PrimeCurve, xR: int, n1: int) -> dict[int, int]:
    """g = NF(S_3(x1, x2, xR)) in F_p[x1, x2]/(x_i^{2^{n1}} - 1), as {chain mask: coefficient}:
    x_i^e is the mask of the bits of e mod 2^{n1} in block i. S_3 and g differ by a
    combination of chain equations, so {g} + chain generates the same ideal; g has
    chain degree <= 2 where S_3 has degree 4 (x_i^2 = x_{i,1})."""
    N = 1 << n1
    g: dict[int, int] = {}
    for (e1, e2), c in s3_prime_coeffs(E, xR).items():
        k = (e1 % N) | ((e2 % N) << n1)
        g[k] = (g.get(k, 0) + c) % E.p
    return {k: c for k, c in g.items() if c}


def prime_macaulay(g: dict[int, int], n1: int, D: int, p: int) -> tuple[np.ndarray, int]:
    """Macaulay matrix at degree D of {g} + chain (m = 2), reduced by the chain: rows
    NF(mu * g) for mu squarefree in the 2 n1 chain variables, deg mu + deg g <= D.
    Monomials multiply by adding the exponents E_i blockwise mod 2^{n1}."""
    v = 2 * n1
    tm = np.asarray(list(g), dtype=np.int64)
    tc = np.asarray(list(g.values()), dtype=np.int64)
    dg = int(_popcount(tm.astype(np.uint64)).max())
    colsD = n_cols(v, D)
    if D < dg:
        return np.zeros((0, colsD), dtype=np.int64), colsD
    mus = _masks_upto(v, D - dg).astype(np.int64)
    N = 1 << n1
    lo = N - 1
    a1 = ((mus[:, None] & lo) + (tm[None, :] & lo)) % N
    a2 = (((mus[:, None] >> n1) & lo) + ((tm[None, :] >> n1) & lo)) % N
    ci = _cols_of(v, D, (a1 | (a2 << n1)).ravel().astype(np.uint64)).reshape(a1.shape)
    M = np.zeros((len(mus), colsD), dtype=np.int64)
    rows = np.repeat(np.arange(len(mus)), len(tm))
    np.add.at(M, (rows, ci.ravel()), np.tile(tc, len(mus)))
    return M % p, colsD


def _rank_mod_p(rows: list[list[int]], p: int) -> int:
    if not rows:
        return 0
    r, _, _ = elim.fp_echelon(np.asarray(rows, dtype=np.int64), p)
    return r


def prime_profile(E: PrimeCurve, V: list[int], xR: int, n1: int, sols: list[tuple[int, int]],
                  dmax: int, max_cols: int = 40000, mutant: bool = True,
                  max_cells: int = 2 * 10 ** 7) -> DegreeProfile:
    """d_solve of the PKM system {g} + chain (g = NF(S_3), `pkm_generator`). d_reg_sr is
    the affine semi-regular value for that system: one equation of degree deg g and
    2 n1 chain equations of degree 2 in 2 n1 unknowns."""
    v = 2 * n1
    p = E.p
    g = pkm_generator(E, xR, n1)
    dg = int(_popcount(np.asarray(list(g), dtype=np.uint64)).max())
    degs = [dg] + [2] * v
    prof = DegreeProfile(v=v, n_eq=1 + v, degrees=degs, d_reg_sr=affine_dreg(v, degs))
    # target: affine-linear functions in the chain variables x_{i,k} = x_i^{2^k} vanishing on W
    pts = []
    for i, j in sols:
        pts.append([1] + [pow(V[i], 1 << k, p) for k in range(n1)] + [pow(V[j], 1 << k, p) for k in range(n1)])
    prof.lin_target = v + 1 - _rank_mod_p(pts, p)
    # column order of _masks_upto: popcount descending, so the last v + 1 columns are
    # the degree-1 masks then the constant; the chain variable x_{i,k} is bit i*n1 + k
    for D in range(dg, dmax + 1):
        if n_cols(v, D) > max_cols or n_cols(v, D) * n_cols(v, D - dg) > max_cells:
            prof.capped = True
            break
        M, cols = prime_macaulay(g, n1, D, p)
        rank, piv, ops = elim.fp_echelon(M, p)
        lin = sum(1 for c in piv if c >= cols - (v + 1))
        prof.ops[D] = ops
        prof.ranks[D] = (M.shape[0], cols, rank, lin, None, None)
        assert lin <= prof.lin_target, "row space has a linear form that does not vanish on the solutions"
        if lin == prof.lin_target:
            prof.d_solve = D
            break
    if mutant:
        _mxl_scan(prof, lambda D: prime_macaulay(g, n1, D, p), v, p, _chain_mulvar(n1, 2), dg, dmax, max_cols)
    return prof
