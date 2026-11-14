"""Generators of the m-summand direct systems (m = 2, 3, 4, 5, ...) by EVALUATION, solution sets, instance builders.

Why evaluation.  S_{m+1}(x_1, ..., x_m, t) is too large to expand symbolically for m >= 4.  It is evaluated at a
point through Semaev's resultant recursion

    S_n(y_1..y_n) = Res_X( S_{j+1}(y_1..y_j, X), S_{n-j+1}(y_{j+1}..y_n, X) ),     S_3 = A X^2 + B X + C,

with numbers substituted for every variable except X.  The normalization is FIXED: each factor is a polynomial in X of
FORMAL degree 2^(j-1) (the degree of S_{j+1} in its last variable), its coefficient list is padded with zeros, and the
resultant is the determinant of the Sylvester matrix of that formal size (dg rows of f, df rows of g, coefficients
listed from the highest formal degree down).  The result is therefore a polynomial in the remaining variables (the
resultant of formal-degree polynomials), not a quantity that changes with the vanishing of a leading coefficient.
The coefficient list of S_{j+1}(y_1..y_j, X) for j >= 3 is obtained by interpolation of values at 2^(j-1) + 1 distinct
points (the polynomial has formal degree 2^(j-1) in X).  Verified in ``scripts/ic_bounds_msum_run.py validate``.

Prime family.  H = {g0^e}, 2^k | p - 1, chains u_{i,j}^2 = u_{i,j+1}.  The normal form of a polynomial modulo the chains
is determined by its values on the grid H^m (the chain quotient is the algebra of functions on the grid): the
coefficient of the monomial with exponent vector a = (a_1..a_m), 0 <= a_i < 2^k, is
N^(-m) sum_e f(g0^e) g0^(-a.e)  (an m-dimensional DFT over F_p, N = 2^k).  The chain mask has the binary digits of a_i in
block i (bit j of block i = variable u_{i,j}); chain degree = number of one bits.  This is exactly the convention of
``csd.ic_bounds_m3.chain_nf``.

Binary family.  The coordinate functions of S_{m+1}(x_1..x_m, x_R) over F_(2^n) on the 2^(mk) assignments of the Boolean
variables (x_i = sum_j c_{i,j} v_j, variable index i*k + j) have algebraic normal form given by the Moebius transform
of the truth table (XOR transform over all mk axes; all n coordinates at once on uint64 words).  This is the Boolean
normal form (modulo x^2 = x) that ``csd.ic_dreg.binary_descent`` computes for m = 2, 3.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from itertools import product

import numpy as np

from . import ic_dreg as I
from .ic_curve_sweep import subgroup_generator_value
from .ic_bounds_m3 import subgroup

# ---------------------------------------------------------------------------
# field operations
# ---------------------------------------------------------------------------


class PrimeOps:
    char = 0

    def __init__(self, p):
        self.p = p
        self.zero, self.one = 0, 1

    def add(self, a, b):
        return (a + b) % self.p

    def sub(self, a, b):
        return (a - b) % self.p

    def mul(self, a, b):
        return a * b % self.p

    def inv(self, a):
        return pow(a, -1, self.p)

    def neg(self, a):
        return (-a) % self.p

    def const(self, c):
        return c % self.p


class BinOps:
    """F_(2^n) with log / antilog tables (n <= 22).  The tables are built with the same modulus as ic_dreg's field."""
    char = 2

    def __init__(self, F):
        self.F, self.n = F, F.N
        self.zero, self.one = 0, 1
        order = (1 << self.n) - 1
        self.order = order
        for cand in range(2, 1 << self.n):
            exp = [0] * (2 * order)
            x = 1
            ok = True
            for i in range(order):
                exp[i] = x
                x = F.mul(x, cand)
                if x == 1 and i + 1 < order:
                    ok = False
                    break
            if ok and x == 1:
                break
        else:  # pragma: no cover
            raise ValueError('no primitive element found')
        for i in range(order):
            exp[order + i] = exp[i]
        log = [0] * (1 << self.n)
        for i in range(order):
            log[exp[i]] = i
        self.exp, self.log = exp, log

    def add(self, a, b):
        return a ^ b

    sub = add

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        return self.exp[self.order - self.log[a]]

    def neg(self, a):
        return a

    def const(self, c):
        return c


# ---------------------------------------------------------------------------
# curve-specific S_3 and the resultant recursion
# ---------------------------------------------------------------------------


class Curve:
    """Bundles the field operations and the S_3 coefficient formulas of one curve."""

    def __init__(self, ops, kind, a=0, b=0):
        self.ops, self.kind, self.a, self.b = ops, kind, a, b

    def s3_abc(self, u, v):
        """(A, B, C) with S_3(u, v, X) = A X^2 + B X + C."""
        o = self.ops
        uv = o.mul(u, v)
        if self.kind == 'prime':
            d = o.sub(u, v)
            s = o.add(u, v)
            A = o.mul(d, d)
            B = o.neg(o.mul(o.const(2), o.add(o.mul(s, o.add(uv, self.a)), o.const(2 * self.b))))
            w = o.sub(uv, self.a)
            C = o.sub(o.mul(w, w), o.mul(o.const(4 * self.b), s))
            return A, B, C
        s = o.add(u, v)
        return o.mul(s, s), uv, o.add(o.mul(uv, uv), self.b)

    def s3_value(self, u, v, w):
        A, B, C = self.s3_abc(u, v)
        o = self.ops
        return o.add(o.add(o.mul(A, o.mul(w, w)), o.mul(B, w)), C)


def det(o, M):
    """Determinant by Gaussian elimination over the field (list of lists, modified copy)."""
    n = len(M)
    M = [row[:] for row in M]
    d = o.one
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            return o.zero
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = o.neg(d)
        pv = M[c][c]
        d = o.mul(d, pv)
        ip = o.inv(pv)
        rowc = M[c]
        for r in range(c + 1, n):
            f = M[r][c]
            if f != 0:
                f = o.mul(f, ip)
                rr = M[r]
                for j in range(c, n):
                    if rowc[j] != 0:
                        rr[j] = o.sub(rr[j], o.mul(f, rowc[j]))
    return d


def sylvester_det(o, f, g):
    """Resultant of formal-degree polynomials: f, g coefficient lists LOW -> HIGH of formal lengths df+1, dg+1."""
    df, dg = len(f) - 1, len(g) - 1
    n = df + dg
    fh, gh = f[::-1], g[::-1]
    M = []
    for i in range(dg):
        row = [o.zero] * n
        for j, c in enumerate(fh):
            row[i + j] = c
        M.append(row)
    for i in range(df):
        row = [o.zero] * n
        for j, c in enumerate(gh):
            row[i + j] = c
        M.append(row)
    return det(o, M)


def default_split(n):
    """number j of variables in the first factor for S_n (n variables): 2 up to n = 5, then balanced."""
    return 2 if n <= 5 else n // 2


class Summation:
    """Evaluator of S_n at points of a curve; ``fast`` uses the closed form for two quadratics (n = 4)."""

    def __init__(self, curve, fast=True):
        self.c, self.o, self.fast = curve, curve.ops, fast
        self._pts = {}

    def points(self, count):
        if count not in self._pts:
            self._pts[count] = [self.o.const(i) for i in range(count)]
        return self._pts[count]

    def poly_X(self, ys, split=None):
        """Coefficients LOW -> HIGH (formal length 2^(j-1) + 1, j = len(ys) >= 2) of S_{j+1}(ys, X)."""
        j = len(ys)
        o = self.o
        if j == 2:
            A, B, C = self.c.s3_abc(*ys)
            return [C, B, A]
        d = 1 << (j - 1)
        xs = self.points(d + 1)
        vals = [self.value(list(ys) + [x], split) for x in xs]
        return self.interpolate(xs, vals)

    def interpolate(self, xs, vals):
        o = self.o
        d = len(xs)
        out = [o.zero] * d
        for i in range(d):
            num = [o.one]
            den = o.one
            for j in range(d):
                if j == i:
                    continue
                new = [o.zero] * (len(num) + 1)
                for t, c in enumerate(num):
                    new[t + 1] = o.add(new[t + 1], c)
                    new[t] = o.sub(new[t], o.mul(xs[j], c))
                num = new
                den = o.mul(den, o.sub(xs[i], xs[j]))
            f = o.mul(vals[i], o.inv(den))
            for t in range(d):
                out[t] = o.add(out[t], o.mul(f, num[t]))
        return out

    def value(self, ys, split=None):
        """S_n(ys), n = len(ys) >= 3, as the fixed-formal-degree resultant recursion."""
        o = self.o
        n = len(ys)
        if n == 3:
            return self.c.s3_value(*ys)
        j = default_split(n) if split is None else split
        if self.fast and n == 4 and j == 2:
            A, B, C = self.c.s3_abc(ys[0], ys[1])
            A2, B2, C2 = self.c.s3_abc(ys[2], ys[3])
            m = o.mul
            u = o.sub(m(A, C2), m(A2, C))
            v = o.sub(m(A, B2), m(A2, B))
            w = o.sub(m(B, C2), m(B2, C))
            return o.sub(m(u, u), m(v, w))
        f = self.poly_X(ys[:j])
        g = self.poly_X(ys[j:])
        return sylvester_det(o, f, g)


def prime_curve_ops(p, a, b):
    return Curve(PrimeOps(p), 'prime', a % p, b % p)


def binary_curve_ops(F, a, b):
    return Curve(BinOps(F), 'binary', a, b)


# ---------------------------------------------------------------------------
# prime family
# ---------------------------------------------------------------------------


def prime_grid(E, t, m, k, split=None, fast=True):
    """Values of S_{m+1}(g0^e_1, ..., g0^e_m, t) on the grid, array of shape (N,)*m (axis i = block i)."""
    p = E.p
    g0, H = subgroup(p, k)
    N = 1 << k
    S = Summation(prime_curve_ops(p, E.a, E.b), fast=fast)
    out = np.zeros((N,) * m, dtype=np.int64)
    flat = out.reshape(-1)
    for idx, es in enumerate(product(range(N), repeat=m)):
        flat[idx] = S.value([H[e] for e in es] + [t % p], split)
    return out, g0


def inverse_dft(grid, p, g0, k):
    """Coefficient array c[a_1..a_m] with grid[e] = sum_a c[a] g0^(a.e)."""
    N = 1 << k
    m = grid.ndim
    assert p < (1 << 26)
    W = np.array([[pow(g0, (-a * e) % N, p) for e in range(N)] for a in range(N)], dtype=np.int64)
    c = grid % p
    for ax in range(m):
        c = np.moveaxis(np.tensordot(W, c, axes=([1], [ax])), 0, ax) % p
    return c * pow(N, -m, p) % p


def forward_dft(coef, p, g0, k):
    N = 1 << k
    W = np.array([[pow(g0, (a * e) % N, p) for a in range(N)] for e in range(N)], dtype=np.int64)
    c = coef % p
    for ax in range(coef.ndim):
        c = np.moveaxis(np.tensordot(W, c, axes=([1], [ax])), 0, ax) % p
    return c


def coef_to_masks(coef, k):
    """{mask: c}: block i of the mask holds the binary digits of the exponent a_i (0 coefficients dropped)."""
    m = coef.ndim
    out = {}
    for idx, a in enumerate(product(range(1 << k), repeat=m)):
        c = int(coef.reshape(-1)[idx])
        if c:
            out[sum(ai << (i * k) for i, ai in enumerate(a))] = c
    return out


def prime_generator(E, t, m, k, split=None, fast=True):
    """(g, grid, g0): g = {mask: c} chain normal form of S_{m+1}(x_1..x_m, t) (raw normalization of the recursion)."""
    grid, g0 = prime_grid(E, t, m, k, split, fast)
    coef = inverse_dft(grid, E.p, g0, k)
    g = coef_to_masks(coef, k)
    return g, grid, g0


def prime_zero_set(grid, m, k):
    """Z as chain masks (E_1 | E_2 << k | ...), E_i the exponent of g0 in x_i."""
    N = 1 << k
    out = []
    flat = grid.reshape(-1)
    for idx, es in enumerate(product(range(N), repeat=m)):
        if flat[idx] == 0:
            out.append(sum(e << (i * k) for i, e in enumerate(es)))
    return out


def nf_values_at(g, m, k, p, g0, masks_points):
    """Value of the normal form {mask: c} at chain points (given as masks of the exponent vectors)."""
    N = 1 << k
    vals = []
    for pt in masks_points:
        E = [(pt >> (i * k)) & (N - 1) for i in range(m)]
        tot = 0
        for mask, c in g.items():
            s = sum(((mask >> (i * k)) & (N - 1)) * E[i] for i in range(m)) % N
            tot += c * pow(g0, s, p)
        vals.append(tot % p)
    return vals


def generator_degree(g):
    return max((ms.bit_count() for ms in g), default=None)


def generator_hash(p, k, m, g):
    return hashlib.sha256(json.dumps(dict(p=p, k=k, m=m, g=sorted(g.items())), separators=(',', ':')).encode()).hexdigest()


def carry_hypotheses(g, m, k, p):
    """Hypotheses of Proposition 1 (experiments/ic-bounds-m3/README.md), tested per variable u_{i,j}:
    (a) deg NF(u g) <= d = deg g;  (b') NF(u g) is not a scalar multiple of g (equivalent to (b): x_i^(2^j) g is not
    proportional to g as a function on H^m).  Returns the sorted list of (i, j) satisfying both (empty: the
    proposition does not apply); proposition 1 then gives LFD >= d + 1."""
    N = 1 << k
    d = generator_degree(g)
    out = []
    for i in range(m):
        for j in range(k):
            h = {}
            for mask, c in g.items():
                a = (mask >> (i * k)) & (N - 1)
                nm = (mask & ~((N - 1) << (i * k))) | (((a + (1 << j)) % N) << (i * k))
                h[nm] = (h.get(nm, 0) + c) % p
            h = {ms: c for ms, c in h.items() if c}
            if max((ms.bit_count() for ms in h), default=0) > d:
                continue
            proportional = False
            if set(h) == set(g):
                ms0 = next(iter(g))
                lam = h[ms0] * pow(g[ms0], -1, p) % p
                proportional = all(h[ms] == lam * g[ms] % p for ms in g)
            if not proportional:
                out.append((i, j))
    return out


# ---------------------------------------------------------------------------
# prime instances
# ---------------------------------------------------------------------------


def _rng(*parts):
    return random.Random(int.from_bytes(hashlib.sha256(':'.join(map(str, parts)).encode()).digest()[:8], 'big'))


def rational_decompositions_prime(E, H, R, m, limit=None):
    """Ordered m-tuples in H^m with signed lifts summing to R (group law).  Brute force, small m, k only."""
    lifts = {}
    for x in H:
        P = E.lift_x(x)
        lifts[x] = [] if P is None else list({P, E.neg(P)})
    out = []
    # meet in the middle over partial sums: dictionary sum -> list of prefixes
    def sums(prefix_len):
        acc = {None: [()]}
        for pos in range(prefix_len):
            nxt = {}
            for S, tuples in acc.items():
                for x in H:
                    for P in lifts[x]:
                        T = E.add(S, P)
                        nxt.setdefault(T, []).extend(tp + (x,) for tp in tuples)
            acc = nxt
        return acc
    h1 = m // 2
    left = sums(h1)
    right = sums(m - h1)
    for S, tl in left.items():
        need = E.add(R, E.neg(S))
        if need in right:
            for a in tl:
                for b in right[need]:
                    out.append(a + b)
    return sorted(set(out))


def make_prime_instance(p, m, k, seed, planted, a=None, b=None, label=''):
    """Seeded random nonsingular curve and target.  planted: target = P_1 + ... + P_m with x(P_i) in the order-4 subgroup
    (so it has a rational decomposition over H for every k >= 2); else a random rational target."""
    rng = _rng('msum-prime', p, m, k, seed, planted)
    N4 = min(k, 2)
    while True:
        aa = rng.randrange(1, p) if a is None else a
        bb = rng.randrange(1, p) if b is None else b
        if (4 * aa ** 3 + 27 * bb ** 2) % p == 0:
            continue
        E = I.PrimeCurve(p, aa, bb)
        if planted:
            H2 = subgroup(p, N4)[1]
            pts = [E.lift_x(x) for x in H2]
            pts = [P for P in pts if P is not None]
            if not pts:
                continue
            Ps = [rng.choice(pts) for _ in range(m)]
            Ps = [P if rng.random() < .5 else E.neg(P) for P in Ps]
            R = None
            for P in Ps:
                R = E.add(R, P)
            if R is None:
                continue
            t, y = R
            plant = [list(P) for P in Ps]
        else:
            t = rng.randrange(p)
            P = E.lift_x(t)
            if P is None:
                continue
            y = P[1]
            plant = None
        return dict(family='prime', p=p, a=aa, b=bb, t=t, y=y, m=m, k=k, planted=plant, seed=seed,
                    target_class='planted-H%d' % (1 << N4) if planted else 'random-rational', label=label)


# ---------------------------------------------------------------------------
# binary family
# ---------------------------------------------------------------------------


def binary_generator(E, basis, m, xR, split=None, fast=True, ops=None):
    """(equations, truth) for S_{m+1}(x_1..x_m, xR), x_i in span(basis), Boolean variable index i*k + j.
    equations: sorted list of sorted lists of monomial masks (one per nonzero coordinate), the canonical form of
    ic_bounds_binary.m3_direct_system; truth: uint64 array of the 2^(mk) values."""
    F = E.F
    k = len(basis)
    r = m * k
    ops = BinOps(F) if ops is None else ops
    S = Summation(Curve(ops, 'binary', E.a, E.b), fast=fast)
    V = I.span(basis)
    T = np.zeros(1 << r, dtype=np.uint64)
    msk = (1 << k) - 1
    for c in range(1 << r):
        xs = [V[(c >> (i * k)) & msk] for i in range(m)]
        T[c] = S.value(xs + [xR], split)
    truth = T.copy()
    for i in range(r):
        T = T.reshape(-1, 2, 1 << i)
        T[:, 1, :] ^= T[:, 0, :]
    T = T.reshape(-1)
    n = F.N
    eqs = []
    for bit in range(n):
        sel = np.flatnonzero((T >> np.uint64(bit)) & np.uint64(1))
        if len(sel):
            eqs.append(sorted(int(x) for x in sel))
    return sorted(set(tuple(e) for e in eqs)), truth


def binary_zero_set(truth):
    return sorted(int(c) for c in np.flatnonzero(truth == 0))


def binary_system_hash(kind, variables, equations):
    system = dict(kind=kind, variables=variables, equations=sorted(tuple(sorted(e)) for e in equations if e))
    return hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def boolean_degrees(equations):
    return [max(mk.bit_count() for mk in e) for e in equations]


def rational_decompositions_binary(E, V, R, m):
    """Number of distinct multisets {x_1..x_m} in V^m with signed lifts summing to R (group law), meet in the middle."""
    lifts = {}
    for x in V:
        P = E.lift_x(x)
        lifts[x] = [] if P is None else list({P, E.neg(P)})
    def sums(cnt):
        acc = {None: [()]}
        for _ in range(cnt):
            nxt = {}
            for S, tuples in acc.items():
                for x in V:
                    for P in lifts[x]:
                        nxt.setdefault(E.add(S, P), []).extend(tp + (x,) for tp in tuples)
            acc = nxt
        return acc
    h1 = m // 2
    left, right = sums(h1), sums(m - h1)
    seen = set()
    for S, tl in left.items():
        need = E.add(R, E.neg(S))
        if need in right:
            for a in tl:
                for b in right[need]:
                    seen.add(tuple(sorted(a + b)))
    return seen


def planted_target_m(E, V, m, rng):
    pts = [P for x in V if x and (P := E.lift_x(x)) is not None]
    for _ in range(400):
        if not pts:
            return None
        Ps = [rng.choice(pts) for _ in range(m)]
        Ps = [P if rng.randrange(2) else E.neg(P) for P in Ps]
        R = None
        for P in Ps:
            R = E.add(R, P)
        if R is not None and R[0]:
            return R, [list(P) for P in Ps]
    return None


# ---------------------------------------------------------------------------
# closures: prime ladder (CompactClosure of scripts/ic_bounds_m3_depth.py, generic in m through MixedSpace(m, k, D))
# ---------------------------------------------------------------------------

_DEPTH = None


def load_depth_module():
    """scripts/ic_bounds_m3_depth.py (not a package member) loaded by path; its CompactClosure takes any MixedSpace."""
    global _DEPTH
    if _DEPTH is None:
        import importlib.util
        import os
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts', 'ic_bounds_m3_depth.py')
        spec = importlib.util.spec_from_file_location('ic_bounds_m3_depth_loaded', path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _DEPTH = mod
    return _DEPTH


def compact_closure_m(gens, m, k, p, D, *, target_rank=None, batch=1024, chunk=4096, withhold_row=False):
    from .ic_bounds_m3 import MixedSpace
    depth = load_depth_module()
    sp = MixedSpace(m, k, D, aux=False)
    cl = depth.CompactClosure(sp, p, batch=batch, chunk=chunk)
    rows = [sp.row(g, p) for g in gens if g and generator_degree(g) <= D]
    cl.add_rows(rows)
    how, sweeps = cl.close(target_rank)
    if not cl.contains_one:
        for r in rows:
            assert cl.contains(r), 'input not in the closure'
    return cl, how, sweeps


def prime_ladder(p, m, k, g, pts, g0, *, cap=None, engine='compact', sandwich='auto', dense_max_columns=1500,
                 log=None, true_pts=None):
    """Exact LFD of {chains, g} on m blocks by the rule of experiments/ic-bounds-m3/README.md.
    unit ideal: least D with 1 in V_D.  proper ideal: first D with rank ev(M_<D on Z) = |Z| and dim V_D = |M_<=D| - |Z|;
    LFD = last strict fall at or below D.  pts = Z (chain masks) computed independently of the closure.
    true_pts (negative controls only): if given, V_D is also required to lie inside the true ideal (profile <= ideal).
    Returns dict(status, exact_lfd, certificate, strict_falls, profile rows)."""
    import time
    from . import ic_fall_certificates as C
    from .ic_bounds_m3 import MixedSpace, evaluation_rank
    r = m * k
    cap = r + 1 if cap is None else cap
    d = generator_degree(g) if g else None
    rows, falls = [], []
    prev_rank = None
    status, cert, unit = None, None, None
    reason = 'configured degree cap'
    inclusion_failures = []
    for D in range(cap + 1):
        t0 = time.time()
        sp = MixedSpace(m, k, D, aux=False)
        masks = [ms for ms, _ in sp.mons]
        n = len(masks)
        Mle = [sum(1 for ms in masks if ms.bit_count() <= j) for j in range(D + 1)]
        ev_all = evaluation_rank(pts, masks, k, g0, p, blocks=m) if pts else 0
        ideal_D = n - ev_all
        no_input = (not g) or d > D
        use_sand = (sandwich is True) or (sandwich == 'auto' and n > 1200)
        how = 'no-input'
        if no_input:
            rank, low, one = 0, [0] * (D + 1), False
        elif engine == 'dense' and n <= dense_max_columns:
            cl = C.chain_degree_closure([g], m, k, p, D, max_columns=dense_max_columns, max_cells=200_000_000)
            rank, one = cl.quotient_rank, cl.contains_one
            if one:
                low = [Mle[j] for j in range(D + 1)]
                how = 'unit'
            else:
                assert C.check_closed_exclusion_certificate(cl, [g])
                low = [sum(cl.masks[j].bit_count() <= dd for j in cl.pivots) for dd in range(D + 1)]
                how = 'check_closed_exclusion_certificate'
        else:
            cl, how, _ = compact_closure_m([g], m, k, p, D, target_rank=(ideal_D if use_sand else None))
            rank, one = cl.rank, cl.contains_one
            low = cl.low_profile()
        if true_pts is not None and pts is not None and not one:
            ev_true = [evaluation_rank(true_pts, [ms for ms in masks if ms.bit_count() <= j], k, g0, p, blocks=m)
                       if true_pts else 0 for j in range(D + 1)]
            if any(low[j] > Mle[j] - ev_true[j] for j in range(D + 1)):
                inclusion_failures.append(D)
        fall = prev_rank is not None and D >= 1 and low[D - 1] > prev_rank
        if fall:
            falls.append(D)
        row = dict(D=D, monomials_le_D=n, quotient_rank=rank, profile=low, contains_one=bool(one), strict_fall=bool(fall),
                   closure_check=how, seconds=round(time.time() - t0, 2))
        if not one and pts and D >= 1:
            low_masks = [ms for ms in masks if ms.bit_count() < D]
            er = evaluation_rank(pts, low_masks, k, g0, p, blocks=m)
            row.update(evaluation_rank_lt_D=er, evaluation_rank_le_D=ev_all, ideal_dim_le_D=ideal_D,
                       codim=n - rank)
            assert all(low[j] <= Mle[j] - evaluation_rank(pts, [ms for ms in masks if ms.bit_count() <= j], k, g0, p,
                                                           blocks=m) for j in range(D + 1)) or true_pts is not None, \
                'V_D not inside I cap M_<=j: generator and Z inconsistent'
            if er == len(pts) == n - rank:
                status = 'exact-proper-local'
                cert = dict(kind='local-to-global', degree=D, solutions=len(pts))
        rows.append(row)
        if log:
            log(f'  m={m} k={k} D={D} n={n} rank={rank} fall={fall} one={one} {how} {row["seconds"]}s')
        prev_rank = rank
        if one:
            assert not pts, 'unit ideal with a nonempty zero set'
            status, unit = 'exact-unit', D
            cert = dict(kind='unit', degree=D)
            break
        if status == 'exact-proper-local':
            break
    last = max(falls, default=0)
    out = dict(status=status or 'capped', exact_lfd=(unit if unit is not None else (last if status else None)),
               lfd_lower_bound=(unit if unit is not None else last), kernel_upper_bound=r + 1, strict_falls=falls,
               certificate=cert, profile=rows, solutions=len(pts), generator_degree=d,
               generator_terms=len(g or {}), cap_reason=None if status else reason)
    if true_pts is not None:
        out['inclusion_failures_vs_true_Z'] = inclusion_failures
    return out


# ---------------------------------------------------------------------------
# binary ladder
# ---------------------------------------------------------------------------


def binary_ladder(equations, r, points, *, cap=None, engine='py', max_columns=70000, max_bytes=900_000_000, work_gb=1.5,
                  batch=512, log=None):
    """Exact LFD of the Boolean system by csd.ic_local_certificate.local_profile (engine 'py') or the bit-packed
    csd.ic_bounds_bin_big.local_profile_bin on the cpu (engine 'big'); same definition, same certificate."""
    cap = r + 1 if cap is None else cap
    if engine == 'py':
        from .ic_local_certificate import local_profile
        prof = local_profile([list(e) for e in equations], r, points, max_degree=cap, max_columns=max_columns,
                             max_bytes=max_bytes)
        return dict(status=prof['status'], exact_lfd=prof['exact_lfd'], lfd_lower_bound=prof.get('lfd_lower_bound'),
                    strict_falls=prof['strict_falls'], quotient_ranks=prof['quotient_ranks'],
                    local_checks=prof['local_checks'], certificate_degree=prof.get('local_upper_bound',
                                                                                 prof['exact_lfd']),
                    cap_reason=prof.get('cap_reason'), engine='csd.ic_local_certificate.local_profile')
    from . import ic_bounds_bin_big as EB
    kw = dict(device='cpu', batch=batch, tile_rows=None, tile_mb=256.0, work_gb=work_gb, col_block=None,
              compact_frac=0.25, tiles_on='device', withhold=None, repair=True, status=None, quiet=True)
    prof = EB.local_profile_bin([list(e) for e in equations], r, list(points), cap, kw, log=log)
    return dict(status=prof['status'], exact_lfd=prof['exact_lfd'], strict_falls=prof['strict_falls'],
                quotient_ranks=prof['quotient_ranks'], local_checks=prof['local_checks'],
                certificate_degree=prof['certificate_degree'], cap_reason=None,
                rref_sha256=[x['rref_sha256'] for x in prof['per_degree']],
                engine='csd.ic_bounds_bin_big.local_profile_bin(cpu)')


def binary_degree_and_zeros(E, basis, m, xR, split=None, fast=True):
    """(input degree d_in, |Z|, number of monomials over all coordinates) of the direct binary system, without building the equation lists
    (large r): same evaluation and Moebius transform as binary_generator.  Added after the predictions freeze (new function, no change of the others)."""
    F = E.F
    k = len(basis)
    r = m * k
    S_ = Summation(Curve(BinOps(F), 'binary', E.a, E.b), fast=fast)
    V = I.span(basis)
    T = np.zeros(1 << r, dtype=np.uint64)
    msk = (1 << k) - 1
    for c in range(1 << r):
        T[c] = S_.value([V[(c >> (i * k)) & msk] for i in range(m)] + [xR], split)
    zeros = int(np.count_nonzero(T == 0))
    for i in range(r):
        T = T.reshape(-1, 2, 1 << i)
        T[:, 1, :] ^= T[:, 0, :]
    T = T.reshape(-1)
    nz = np.flatnonzero(T)
    deg = int(max((bin(int(x)).count('1') for x in nz), default=0)) if len(nz) < 3_000_000 else int(
        np.max(np.array([bin(int(x)).count('1') for x in nz[:: max(1, len(nz) // 2_000_000)]])))
    return deg, zeros, int(sum(bin(int(v)).count('1') for v in T[nz])) if len(nz) < 3_000_000 else None
