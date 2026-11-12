"""Exact ordinary iterated LFD of three-summand (m = 3) prime subgroup-chain systems.

Setting.  E: y^2 = x^3 + a x + b over F_p, H the subgroup of order 2^k of F_p^* (2^k | p-1),
H = {g0^e : 0 <= e < 2^k} with g0 = ic_curve_sweep.subgroup_generator_value(p, k).  Each summand
x_i (i = 1, 2, 3) is written through the chain u_{i,0} = x_i, u_{i,j}^2 = u_{i,j+1},
u_{i,k-1}^2 = 1; modulo the chains a monomial is a squarefree mask whose block i holds the bits of
the exponent e_i of x_i (mod 2^k), and multiplication by u_{i,j} adds 2^j to e_i modulo 2^k.

Encodings.

* ``direct``: inputs = the 3k chain equations and g = NF(S_4(x1, x2, x3, t)), where t is the
  x-coordinate of the target R and S_4(x1,x2,x3,x4) = Res_X(S_3(x1,x2,X), S_3(x3,x4,X)) (Semaev).
  S_4 has degree 4 in each variable; NF(S_4) has chain degree <= 6 (exponents 0..4 reduce to
  masks of weight <= 2 per block; <= 3 for k = 1).  The closure is csd.ic_fall_certificates
  .chain_degree_closure with m = 3 (dense, modulo chains); the exact value uses the
  local-to-global certificate with Z enumerated over H^3.

* ``split``: inputs = the 3k chain equations and
      g1 = NF(S_3(x1, x2, z)),   g2 = NF(S_3(z, x3, t)),
  with ONE auxiliary coordinate z (the x-coordinate of P1 + P2) that is an ordinary polynomial
  variable, NOT a chain variable: no subgroup restriction and no chain is imposed on z.  Two
  ideals share the same closures V_D for every D < p:
      I_free  = (chains, g1, g2)               (z ranges over the algebraic closure),
      I_field = (chains, g1, g2, z^p - z)      (z ranges over the full field F_p).
  The field equation has degree p and never enters V_D for D < p, so the closure is computed once
  and the local-to-global certificate is checked against both quotients.  This differs from
  csd.ic_research_campaign.prime_split, which writes the auxiliary as a chain variable in H.

Local-to-global certificate (general form, no radicality needed).  Let A = R/I be finite
dimensional.  If the images of the monomials of degree < D span A and dim V_D = |M_<=D| - dim A
(in the chain-quotient representation), then V_D = I cap M_<=D and V_D' = I cap M_<=D' for every
D' >= D; hence the exact LFD is the last strict fall at or below D (falls are computed from the
exact closures).  Proof (written in full in experiments/ic-bounds-m3/README.md): every monomial
of degree D'+1 is x * mu with deg mu = D'; mu = rho(mu) + v with deg rho < D and v in V_D'; x*v is
an admissible product in V_D'+1 and x*rho(mu) has degree <= D, so every element of M_<=D'+1 is
congruent mod V_D'+1 to an element of M_<D; an ideal element congruent to rho in M_<D has
rho in I cap M_<D subset V_D.
For I_free, A = (+)_{x in H^3} F_p[z]/(h_x) with h_x = gcd(g1(x,z), g2(x,z)) (A_chain = F_p^{H^3});
for I_field, h_x is replaced by gcd(h_x, z^p - z).  For the direct system A = F_p^Z.
"""
from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations, product

import numpy as np

from . import elim
from . import ic_dreg as I
from . import ic_fall_certificates as C
from .ic_curve_sweep import subgroup_generator_value

# ---------------------------------------------------------------------------
# multivariate polynomials over F_p as {exponent tuple: coefficient}
# ---------------------------------------------------------------------------


def padd(p, *terms):
    out = {}
    for scale, f in terms:
        for e, c in f.items():
            out[e] = (out.get(e, 0) + scale * c) % p
    return {e: c for e, c in out.items() if c}


def pmul(p, f, g):
    out = {}
    for e1, c1 in f.items():
        for e2, c2 in g.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = (out.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in out.items() if c}


def s3_quadratic(E, u, v, nvars):
    """Coefficients (A, B, C) of S_3(u, v, X) = A X^2 + B X + C, u and v polynomials in nvars
    variables: A = (u - v)^2, B = -2((u + v)(u v + a) + 2b), C = (u v - a)^2 - 4b(u + v)."""
    p, a, b = E.p, E.a, E.b
    one = {(0,) * nvars: 1}
    uv = pmul(p, u, v)
    d = padd(p, (1, u), (-1, v))
    A = pmul(p, d, d)
    s = padd(p, (1, u), (1, v))
    B = padd(p, (-2, pmul(p, s, padd(p, (1, uv), (a, one)))), (-4 * b, one))
    w = padd(p, (1, uv), (-a, one))
    Cc = padd(p, (1, pmul(p, w, w)), (-4 * b, s))
    return A, B, Cc


def s4_polynomial(E, t):
    """S_4(x1, x2, x3, t) = Res_X(S_3(x1, x2, X), S_3(x3, t, X)) as {(e1, e2, e3): c mod p}.
    Res(A X^2 + B X + C, A' X^2 + B' X + C') = (A C' - A' C)^2 - (A B' - A' B)(B C' - B' C)."""
    p = E.p
    x1, x2, x3 = {(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}
    tc = {(0, 0, 0): t % p}
    A, B, Cc = s3_quadratic(E, x1, x2, 3)
    A2, B2, C2 = s3_quadratic(E, x3, tc, 3)
    m = lambda f, g: pmul(p, f, g)
    u = padd(p, (1, m(A, C2)), (-1, m(A2, Cc)))
    v = padd(p, (1, m(A, B2)), (-1, m(A2, B)))
    w = padd(p, (1, m(B, C2)), (-1, m(B2, Cc)))
    return padd(p, (1, m(u, u)), (-1, m(v, w)))


def s4_value(E, x1, x2, x3, x4):
    """Numerical S_4 through the resultant formula (independent of s4_polynomial's expansion)."""
    p, a, b = E.p, E.a, E.b

    def abc(u, v):
        return ((u - v) ** 2, -2 * ((u + v) * (u * v + a) + 2 * b), (u * v - a) ** 2 - 4 * b * (u + v))
    A, B, Cc = abc(x1, x2)
    A2, B2, C2 = abc(x3, x4)
    return ((A * C2 - A2 * Cc) ** 2 - (A * B2 - A2 * B) * (B * C2 - B2 * Cc)) % p


def poly_value(p, f, xs):
    return sum(c * math.prod(pow(x, e, p) for x, e in zip(xs, ex)) for ex, c in f.items()) % p


def chain_nf(p, f, k, chain_vars):
    """Reduce the exponents of the first chain_vars variables modulo 2^k into block masks.
    Returns {(mask, rest exponents...): c} when extra variables exist, else {mask: c}."""
    N = 1 << k
    out = {}
    for e, c in f.items():
        mask = sum((e[i] % N) << (i * k) for i in range(chain_vars))
        key = mask if len(e) == chain_vars else (mask,) + tuple(e[chain_vars:])
        out[key] = (out.get(key, 0) + c) % p
    return {key: c for key, c in out.items() if c}


def direct_generator(E, t, k):
    return chain_nf(E.p, s4_polynomial(E, t), k, 3)


def split_generators(E, t, k):
    """(g1, g2) as {(mask, z exponent): c}; variables (x1, x2, x3, z)."""
    p = E.p
    x1, x2, x3 = {(1, 0, 0, 0): 1}, {(0, 1, 0, 0): 1}, {(0, 0, 1, 0): 1}
    z = {(0, 0, 0, 1): 1}
    tc = {(0, 0, 0, 0): t % p}
    A, B, Cc = s3_quadratic(E, x1, x2, 4)
    g1 = padd(p, (1, pmul(p, A, pmul(p, z, z))), (1, pmul(p, B, z)), (1, Cc))
    A, B, Cc = s3_quadratic(E, z, x3, 4)
    tt = {(0, 0, 0, 0): t * t % p}
    g2 = padd(p, (1, pmul(p, A, tt)), (1, pmul(p, B, {(0, 0, 0, 0): t % p})), (1, Cc))
    return chain_nf(p, g1, k, 3), chain_nf(p, g2, k, 3), (g1, g2)


# ---------------------------------------------------------------------------
# factor base, scalar solution sets and rational decompositions
# ---------------------------------------------------------------------------


def subgroup(p, k):
    g0 = subgroup_generator_value(p, k)
    H = [pow(g0, e, p) for e in range(1 << k)]
    assert len(set(H)) == 1 << k and pow(g0, 1 << k, p) == 1
    return g0, H


def nf_value(p, k, nf, E3, g0):
    """Value of a chain normal form {mask: c} at the chain point (E1, E2, E3)."""
    N = 1 << k
    total = 0
    for mask, c in nf.items():
        s = sum(((mask >> (i * k)) & (N - 1)) * E3[i] for i in range(3)) % N
        total += c * pow(g0, s, p)
    return total % p


def direct_points(E, t, k, nf):
    """Z = {(E1, E2, E3) in (Z/2^k)^3 : S_4(g0^E1, g0^E2, g0^E3, t) = 0}, as chain masks.
    Computed from the resultant formula and asserted equal to the zero set of NF(S_4)."""
    p = E.p
    g0, H = subgroup(p, k)
    N = 1 << k
    pts = []
    for E3 in product(range(N), repeat=3):
        xs = [H[e] for e in E3]
        zero = s4_value(E, *xs, t) == 0
        assert zero == (nf_value(p, k, nf, E3, g0) == 0)
        if zero:
            pts.append(E3[0] | (E3[1] << k) | (E3[2] << (2 * k)))
    return pts


def rational_decompositions(E, H, R):
    """Ordered triples (x1, x2, x3) in H^3 with signed lifts P_i and P1 + P2 + P3 = R."""
    lifts = {x: [P for P in {E.lift_x(x), E.neg(E.lift_x(x))}] if E.lift_x(x) is not None else []
             for x in H}
    out = []
    for xs in product(H, repeat=3):
        hit = False
        for P1 in lifts[xs[0]]:
            for P2 in lifts[xs[1]]:
                S = E.add(P1, P2)
                for P3 in lifts[xs[2]]:
                    if E.add(S, P3) == R:
                        hit = True
        if hit:
            out.append(xs)
    return out


def evaluation_rank(points, masks, k, g0, p, blocks=3):
    """Rank of the evaluation of squarefree chain monomials at chain points (rows = points)."""
    if not points or not masks:
        return 0
    N = 1 << k
    P = np.array([[(pt >> (i * k)) & (N - 1) for i in range(blocks)] for pt in points], dtype=np.int64)
    M = np.array([[(ms >> (i * k)) & (N - 1) for i in range(blocks)] for ms in masks], dtype=np.int64)
    expo = (P @ M.T) % N
    table = np.array([pow(g0, e, p) for e in range(N)], dtype=np.int64)
    return elim.fp_echelon(table[expo], p)[0]


# ---------------------------------------------------------------------------
# direct profile (chain_degree_closure, m = 3)
# ---------------------------------------------------------------------------


def masks_upto(r, D):
    return [sum(1 << j for j in c) for d in range(min(r, D) + 1) for c in combinations(range(r), d)]


def direct_profile(E, t, k, *, max_degree=None, max_columns=5000, max_cells=60_000_000,
                   check_exclusion=True, generator=None, points=None):
    """Exact LFD of {chains, NF(S_4(x1,x2,x3,t))}: unit ideal -> least D with 1 in V_D;
    proper ideal -> last strict fall once the local-to-global certificate holds at D."""
    p = E.p
    r = 3 * k
    g0, H = subgroup(p, k)
    g = direct_generator(E, t, k) if generator is None else generator
    pts = direct_points(E, t, k, g) if points is None else points
    maxdeg = r + 1 if max_degree is None else max_degree
    rows = []
    prev = None
    falls = []
    status = None
    cert = None
    reason = 'configured degree cap'
    for D in range(maxdeg + 1):
        try:
            cl = C.chain_degree_closure([g] if g else [], 3, k, p, D, max_columns=max_columns,
                                        max_cells=max_cells)
        except ValueError as error:
            reason = str(error)
            break
        fall = prev is not None and cl.low_degree_dimension(D - 1) > prev.rank
        if fall:
            falls.append(D)
        Mle = len(cl.masks)
        row = dict(D=D, monomials_le_D=Mle, quotient_rank=cl.quotient_rank, dim_V_D=cl.rank,
                   dim_V_D_low=cl.low_degree_dimension(D - 1), contains_one=cl.contains_one,
                   strict_fall=fall, rounds=len(cl.rank_history))
        if not cl.contains_one and check_exclusion:
            ok = C.check_closed_exclusion_certificate(cl, [g] if g else [])
            assert ok, (p, E.a, E.b, t, k, D)
            row['closed_exclusion_certificate'] = ok
        if pts and D >= 1:
            low = [ms for ms in cl.masks if ms.bit_count() < D]
            er = evaluation_rank(pts, low, k, g0, p)
            full = evaluation_rank(pts, list(cl.masks), k, g0, p)
            row.update(evaluation_rank_lt_D=er, evaluation_rank_le_D=full,
                       ideal_dim_le_D=Mle - full, deficiency=(Mle - full) - cl.quotient_rank)
        rows.append(row)
        prev = cl
        if cl.contains_one:
            assert not pts
            status = 'exact-unit'
            cert = dict(kind='unit', degree=D)
            break
        if pts and D >= 1 and row['evaluation_rank_lt_D'] == len(pts) == Mle - cl.quotient_rank:
            status = 'exact-proper-local'
            cert = dict(kind='local-to-global', degree=D, solutions=len(pts))
            break
    if status is None:
        out = dict(status='capped', exact_lfd=None, lfd_lower_bound=max(falls, default=0),
                   lfd_upper_bound=r + 1, cap_reason=reason)
    elif status == 'exact-unit':
        out = dict(status=status, exact_lfd=rows[-1]['D'], lfd_lower_bound=rows[-1]['D'],
                   lfd_upper_bound=rows[-1]['D'])
    else:
        last = max(falls, default=0)
        out = dict(status=status, exact_lfd=last, lfd_lower_bound=last, lfd_upper_bound=last)
    out.update(strict_falls=falls, certificate=cert, profile=rows, solutions=len(pts),
               generator_terms=len(g), generator_degree=max((ms.bit_count() for ms in g), default=None))
    return out


# ---------------------------------------------------------------------------
# mixed closure: chain blocks plus one ordinary auxiliary variable z
# ---------------------------------------------------------------------------


class MixedSpace:
    """Monomials (mask, e) = u^mask z^e with weight(mask) + e <= D, sorted by degree descending."""

    def __init__(self, m, k, D, aux=True):
        self.m, self.k, self.D, self.aux = m, k, D, aux
        r = m * k
        mons = []
        for d in range(D + 1):
            for e in range(d + 1 if aux else 1):
                w = d - e
                if w > r:
                    continue
                for c in combinations(range(r), w):
                    mons.append((sum(1 << j for j in c), e))
        mons.sort(key=lambda me: (me[0].bit_count() + me[1], me[1], me[0]), reverse=True)
        self.mons = mons
        self.index = {me: i for i, me in enumerate(mons)}
        self.deg = np.array([ms.bit_count() + e for ms, e in mons], dtype=np.int64)
        N = 1 << k
        nvar = r + (1 if aux else 0)
        cols = np.full((nvar, len(mons)), -1, dtype=np.int64)
        for s, (ms, e) in enumerate(mons):
            if ms.bit_count() + e >= D:
                continue
            for j in range(r):
                blk, bit = divmod(j, k)
                sh = blk * k
                ex = ((ms >> sh) + (1 << bit)) % N
                cols[j, s] = self.index[((ms & ~((N - 1) << sh)) | (ex << sh), e)]
            if aux:
                cols[r, s] = self.index[(ms, e + 1)]
        self.cols = cols

    def __len__(self):
        return len(self.mons)

    def row(self, poly, p):
        v = np.zeros(len(self.mons), dtype=np.int64)
        for key, c in poly.items():
            me = key if isinstance(key, tuple) else (key, 0)
            v[self.index[me]] = (v[self.index[me]] + c) % p
        return v

    def ordinary_count(self):
        nvar = self.m * self.k + (1 if self.aux else 0)
        return math.comb(nvar + self.D, self.D)


def _mixed_products(space, rows, pivots, D):
    low = [row for row, pv in zip(rows, pivots) if space.deg[pv] < D]
    if not low:
        return np.zeros((0, len(space)), dtype=np.int64)
    out = np.zeros((len(low) * len(space.cols), len(space)), dtype=np.int64)
    for i, row in enumerate(low):
        sup = np.flatnonzero(row)
        for j in range(len(space.cols)):
            dest = space.cols[j, sup]
            assert np.all(dest >= 0)
            out[i * len(space.cols) + j, dest] = row[sup]
    return out


def mixed_closure(generators, m, k, p, D, aux=True, max_cells=60_000_000):
    """Exact V_D of {chains} + generators in F_p[u, z], represented modulo the chains.
    generators: dicts {(mask, e): c} (or {mask: c} when aux is False), input degree = degree of
    the normal form.  Saturates by variable products of rows of degree < D, then VERIFIES closure:
    all products of the final low rows lie in the span.  Returns dict with basis and pivots."""
    sp = MixedSpace(m, k, D, aux)
    n = len(sp)
    rows = []
    for g in generators:
        if not g:
            continue
        deg = max((key[0].bit_count() + key[1]) if isinstance(key, tuple) else key.bit_count()
                  for key in g)
        if deg <= D:
            rows.append(sp.row(g, p))
    mat = np.asarray(rows, dtype=np.int64).reshape((-1, n))
    rank_hist = []
    while True:
        out = []
        rank, piv, _ = elim.fp_echelon(mat, p, out)
        basis = out[0] if rank else np.zeros((0, n), dtype=np.int64)
        rank_hist.append(rank)
        if n - 1 in piv:
            basis = np.eye(n, dtype=np.int64)
            piv = list(range(n))
            rank = n
            break
        if len(rank_hist) > 1 and rank_hist[-1] == rank_hist[-2]:
            break
        prods = _mixed_products(sp, basis, piv, D)
        if not len(prods):
            break
        if (len(basis) + len(prods)) * n > max_cells:
            raise ValueError('mixed closure cell budget exceeded')
        mat = np.vstack([basis, prods])
    contains_one = n - 1 in piv
    verified = True
    if not contains_one:
        # independent closure verification: inputs in span and all low products in span
        prods = _mixed_products(sp, basis, piv, D)
        gens = np.asarray(rows, dtype=np.int64).reshape((-1, n))
        verified = elim.fp_echelon(np.vstack([basis, gens, prods]), p)[0] == rank
        assert verified
    return dict(space=sp, basis=basis, pivots=list(piv), rank=rank, contains_one=contains_one,
                verified_closed=verified, rounds=len(rank_hist))


def low_dimension(cl, d):
    sp = cl['space']
    return int(sum(1 for pv in cl['pivots'] if sp.deg[pv] <= d))


def ordinary_dimension(cl):
    sp = cl['space']
    standard = len(sp)
    return sp.ordinary_count() - standard + cl['rank']


# ---------------------------------------------------------------------------
# split encoding: quotient algebra and spanning check
# ---------------------------------------------------------------------------


def upoly_trim(f, p):
    f = [c % p for c in f]
    while f and f[-1] == 0:
        f.pop()
    return f


def upoly_mod(f, g, p):
    f = upoly_trim(f, p)
    g = upoly_trim(g, p)
    inv = pow(g[-1], -1, p)
    while len(f) >= len(g):
        c = f[-1] * inv % p
        s = len(f) - len(g)
        for i, gc in enumerate(g):
            f[s + i] = (f[s + i] - c * gc) % p
        f = upoly_trim(f, p)
    return f


def upoly_gcd(f, g, p):
    f, g = upoly_trim(f, p), upoly_trim(g, p)
    while g:
        f, g = g, upoly_mod(f, g, p)
    if not f:
        return None  # gcd(0, 0): not zero-dimensional
    inv = pow(f[-1], -1, p)
    return [c * inv % p for c in f]


def upoly_mulmod(f, g, h, p):
    out = [0] * (len(f) + len(g))
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i + j] = (out[i + j] + a * b) % p
    return upoly_mod(out, h, p)


def upoly_powmod_z(e, h, p):
    result, base = [1], [0, 1]
    result = upoly_mod(result, h, p) if len(h) > 1 else []
    base = upoly_mod(base, h, p)
    while e:
        if e & 1:
            result = upoly_mulmod(result, base, h, p)
        base = upoly_mulmod(base, base, h, p)
        e >>= 1
    return result


def split_local_algebra(g1, g2, k, p, g0, field_equation):
    """For every chain point x in H^3: h_x = gcd(g1(x, z), g2(x, z)) (and z^p - z if requested).
    Returns {point mask: monic h_x} with deg h_x >= 1, or None if some gcd(0, 0) occurs."""
    N = 1 << k
    out = {}
    for E3 in product(range(N), repeat=3):
        pt = E3[0] | (E3[1] << k) | (E3[2] << (2 * k))
        polys = []
        for g in (g1, g2):
            f = [0, 0, 0, 0, 0]
            for (ms, e), c in g.items():
                s = sum(((ms >> (i * k)) & (N - 1)) * E3[i] for i in range(3)) % N
                f[e] = (f[e] + c * pow(g0, s, p)) % p
            polys.append(f)
        h = upoly_gcd(polys[0], polys[1], p)
        if h is None:
            return None
        if field_equation and len(h) > 1:
            zp = upoly_powmod_z(p, h, p)
            diff = upoly_trim([(zp[i] if i < len(zp) else 0) - (1 if i == 1 else 0)
                               for i in range(max(len(zp), 2))], p)
            h = upoly_gcd(h, diff, p) if diff else h
        if len(h) > 1:
            out[pt] = h
    return out


def split_evaluation_rank(local, mons, k, g0, p):
    """Rank of the images of the monomials u^mask z^e in (+)_x F_p[z]/(h_x)."""
    N = 1 << k
    pts = sorted(local)
    offsets = []
    width = 0
    for pt in pts:
        offsets.append(width)
        width += len(local[pt]) - 1
    if width == 0 or not mons:
        return 0, width
    maxe = max(e for _, e in mons)
    # z^e mod h_x for e <= maxe
    zpow = {}
    for pt in pts:
        h = local[pt]
        cur = upoly_mod([1], h, p) if len(h) > 1 else []
        lst = [cur]
        for _ in range(maxe):
            cur = upoly_mulmod(cur, [0, 1], h, p)
            lst.append(cur)
        zpow[pt] = lst
    mat = np.zeros((len(mons), width), dtype=np.int64)
    for i, (ms, e) in enumerate(mons):
        for pt, off in zip(pts, offsets):
            s = sum(((ms >> (j * k)) & (N - 1)) * ((pt >> (j * k)) & (N - 1)) for j in range(3)) % N
            val = pow(g0, s, p)
            for d, c in enumerate(zpow[pt][e]):
                mat[i, off + d] = val * c % p
    return elim.fp_echelon(mat, p)[0], width


def split_profile(E, t, k, *, max_degree=None, max_cells=60_000_000, max_columns=12000):
    p = E.p
    g0, H = subgroup(p, k)
    g1, g2, raw = split_generators(E, t, k)
    loc = {fe: split_local_algebra(g1, g2, k, p, g0, fe) for fe in (False, True)}
    dims = {fe: (None if loc[fe] is None else sum(len(h) - 1 for h in loc[fe].values()))
            for fe in (False, True)}
    r = 3 * k
    maxdeg = r + 4 if max_degree is None else max_degree
    rows = []
    prev = None
    falls = []
    certs = {False: None, True: None}
    unit = None
    reason = 'configured degree cap'
    for D in range(maxdeg + 1):
        if D >= p:
            reason = 'D >= p: the field equation would enter V_D'
            break
        sp_size = len(MixedSpace(3, k, D, True).mons)
        if sp_size > max_columns:
            reason = f'column budget {sp_size} > {max_columns}'
            break
        try:
            cl = mixed_closure([g1, g2], 3, k, p, D, True, max_cells=max_cells)
        except ValueError as error:
            reason = str(error)
            break
        sp = cl['space']
        ordinary = ordinary_dimension(cl)
        # strict fall: dim(V_D cap M_<=D-1) > dim V_{D-1}, in the ordinary ring
        if prev is not None:
            # ordinary dims: chain-ideal part of degree <= D-1 is the same on both sides
            low_q = low_dimension(cl, D - 1)
            fall = low_q > prev['rank']
        else:
            fall = False
        if fall:
            falls.append(D)
        row = dict(D=D, monomials_le_D=len(sp), quotient_rank=cl['rank'], dim_V_D=ordinary,
                   quotient_low_rank=low_dimension(cl, D - 1), contains_one=cl['contains_one'],
                   strict_fall=fall, verified_closed=cl['verified_closed'], rounds=cl['rounds'])
        for fe in (False, True):
            key = 'field' if fe else 'free'
            if loc[fe] is None:
                row[f'{key}_note'] = 'not zero-dimensional'
                continue
            low = [me for me in sp.mons if me[0].bit_count() + me[1] < D]
            rk, width = split_evaluation_rank(loc[fe], low, k, g0, p) if D >= 1 else (0, dims[fe])
            row[f'{key}_eval_rank_lt_D'] = rk
            row[f'{key}_ideal_dim_le_D_if_spanning'] = len(sp) - dims[fe]
            if (certs[fe] is None and not cl['contains_one'] and D >= 1 and rk == dims[fe]
                    and len(sp) - cl['rank'] == dims[fe]):
                certs[fe] = dict(kind='local-to-global', degree=D, quotient_dimension=dims[fe],
                                 exact_lfd=max(falls, default=0))
        rows.append(row)
        prev = cl
        if cl['contains_one']:
            unit = D
            break
        if all(certs[fe] is not None or loc[fe] is None for fe in (False, True)):
            break
    res = dict(quotient_dimension_free=dims[False], quotient_dimension_field=dims[True],
               strict_falls=falls, profile=rows,
               g1_degree=max(ms.bit_count() + e for ms, e in g1), g2_degree=max(ms.bit_count() + e for ms, e in g2))
    if unit is not None:
        assert dims[False] == 0 and dims[True] == 0
        for key in ('free', 'field'):
            res[key] = dict(status='exact-unit', exact_lfd=unit, certificate=dict(kind='unit', degree=unit))
    else:
        for fe, key in ((False, 'free'), (True, 'field')):
            if certs[fe] is not None:
                res[key] = dict(status='exact-proper-local', exact_lfd=certs[fe]['exact_lfd'],
                                certificate=certs[fe])
            else:
                res[key] = dict(status='capped', exact_lfd=None, lfd_lower_bound=max(falls, default=0),
                                cap_reason=reason if loc[fe] is not None else 'not zero-dimensional')
    return res


def system_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


# ---------------------------------------------------------------------------
# unified exact profile (dense engines for small spaces, batched engine otherwise)
# ---------------------------------------------------------------------------


class _Snap:
    """Engine-independent view of one closure V_D in the chain-quotient representation."""

    def __init__(self, D, n, q_rank, q_low, one, closed_check, engine):
        self.D, self.n, self.q_rank, self.q_low, self.one = D, n, q_rank, q_low, one
        self.closed_check, self.engine = closed_check, engine


def _closure_snapshot(encoding, gens, k, p, D, engine, target_rank=None):
    from . import ic_bounds_m3_batched as BC
    if encoding == 'direct' and engine == 'dense':
        cl = C.chain_degree_closure(gens, 3, k, p, D, max_columns=5000, max_cells=60_000_000)
        check = None
        if not cl.contains_one:
            check = C.check_closed_exclusion_certificate(cl, gens)
            assert check
        lows = [sum(cl.masks[j].bit_count() <= d for j in cl.pivots) for d in range(D + 1)]
        return _Snap(D, len(cl.masks), cl.quotient_rank, lows, cl.contains_one,
                     'check_closed_exclusion_certificate' if check else 'unit', 'chain_degree_closure')
    if encoding == 'split' and engine == 'dense':
        try:
            cl = mixed_closure(gens, 3, k, p, D, True)
        except ValueError:
            # dense cell budget exceeded: same V_D from the memory-bounded engine
            return _closure_snapshot(encoding, gens, k, p, D, 'batched', target_rank)
        lows = [low_dimension(cl, d) for d in range(D + 1)]
        return _Snap(D, len(cl['space']), cl['rank'], lows, cl['contains_one'],
                     'all-low-shifts-reduce' if not cl['contains_one'] else 'unit', 'mixed_closure')
    cl = BC.batched_closure(gens, 3, k, p, D, encoding == 'split', batch=1024, target_rank=target_rank)
    lows = [cl.low_dimension(d) for d in range(D + 1)]
    check = ('unit' if cl.contains_one else
             'rank-equals-dim(I cap M_<=D)' if cl.stopped_at_ideal_dimension else 'all-low-shifts-reduce')
    return _Snap(D, len(cl.sp), cl.rank, lows, cl.contains_one, check, 'batched_closure')


def _ordinary(nvars, D, n_quot, q_rank):
    return math.comb(nvars + D, D) - n_quot + q_rank


def exact_profile(encoding, E, t, k, *, engine=None, max_degree=None, max_columns=None, log=None):
    """Exact ordinary iterated LFD of the direct or split m = 3 system with certificates.

    Returns status 'exact-unit' (1 in V_D, no 1 below), 'exact-proper-local' (local-to-global
    certificate at D: monomials of degree < D span R/I and dim V_D = |M_<=D| - dim R/I; exact LFD =
    last strict fall), or 'capped' (lower bound = last strict fall of exact closures, never a
    stall).  For the split encoding both ideals I_free and I_field are certified separately."""
    p = E.p
    g0, H = subgroup(p, k)
    if engine is None:
        engine = 'dense' if k <= 3 else 'batched'
    if encoding == 'direct':
        g = direct_generator(E, t, k)
        gens = [g] if g else []
        pts = direct_points(E, t, k, g)
        algebras = {'direct': ('points', pts, len(pts))}
        nvars = 3 * k
        gen_degrees = [max((ms.bit_count() for ms in g), default=None)]
        gen_hash = system_hash(dict(p=p, k=k, m=3, encoding='direct', g=sorted(g.items())))
        cap = 3 * k + 1
    else:
        g1, g2, _ = split_generators(E, t, k)
        gens = [g1, g2]
        algebras = {}
        for fe, key in ((False, 'free'), (True, 'field')):
            loc = split_local_algebra(g1, g2, k, p, g0, fe)
            algebras[key] = ('local', loc, None if loc is None else sum(len(h) - 1 for h in loc.values()))
        nvars = 3 * k + 1
        gen_degrees = [max(ms.bit_count() + e for ms, e in g) for g in gens]
        gen_hash = system_hash(dict(p=p, k=k, m=3, encoding='split-free-aux',
                                    g=[sorted([list(kk), c] for kk, c in g.items()) for g in gens]))
        cap = 3 * k + 6
    if max_degree is not None:
        cap = max_degree
    if max_columns is None:
        max_columns = 6000 if encoding == 'split' else 5000
    rows, falls = [], []
    certs = {key: None for key in algebras}
    prev = None
    unit = None
    reason = 'configured degree cap'
    for D in range(cap + 1):
        if encoding == 'split' and D >= p:
            reason = 'D >= p'
            break
        n_est = len(MixedSpace(3, k, D, encoding == 'split').mons)
        if n_est > max_columns:
            reason = f'column budget: |M_<=D| = {n_est} > {max_columns}'
            break
        target_rank = None
        if engine == 'batched' and D >= 1:
            # exact dim(I cap M_<=D) from the quotient algebra (I_free for the split); used only to
            # stop the saturation early, which is sound by S subset V_D subset I cap M_<=D.
            spD = MixedSpace(3, k, D, encoding == 'split')
            kind, data, dim = algebras['direct' if encoding == 'direct' else 'free']
            if dim is not None:
                if kind == 'points':
                    rk_all = evaluation_rank(data, [ms for ms, _ in spD.mons], k, g0, p)
                else:
                    rk_all = split_evaluation_rank(data, spD.mons, k, g0, p)[0]
                target_rank = len(spD) - rk_all
        snap = _closure_snapshot(encoding, gens, k, p, D, engine, target_rank)
        fall = prev is not None and snap.q_low[D - 1] > prev.q_rank
        if fall:
            falls.append(D)
        row = dict(D=D, monomials_le_D=snap.n, quotient_rank=snap.q_rank,
                   dim_V_D=_ordinary(nvars, D, snap.n, snap.q_rank),
                   quotient_low_rank=snap.q_low[D - 1] if D else 0,
                   contains_one=snap.one, strict_fall=fall, closure_check=snap.closed_check,
                   engine=snap.engine)
        if not snap.one and D >= 1:
            sp = MixedSpace(3, k, D, encoding == 'split')
            low = [me for me in sp.mons if me[0].bit_count() + me[1] < D]
            for key, (kind, data, dim) in algebras.items():
                if dim is None:
                    row[f'{key}_note'] = 'not zero-dimensional'
                    continue
                if kind == 'points':
                    rk = evaluation_rank(data, [ms for ms, _ in low], k, g0, p)
                else:
                    rk = split_evaluation_rank(data, low, k, g0, p)[0]
                row[f'{key}_eval_rank_lt_D'] = rk
                row[f'{key}_quotient_dim'] = dim
                if certs[key] is None and rk == dim and snap.n - snap.q_rank == dim:
                    certs[key] = dict(kind='local-to-global', degree=D, quotient_dimension=dim,
                                      exact_lfd=max(falls, default=0))
        rows.append(row)
        if log:
            log(f'    {encoding} k={k} D={D} n={snap.n} q_rank={snap.q_rank} fall={fall} one={snap.one}')
        prev = snap
        if snap.one:
            unit = D
            break
        if encoding == 'split' and certs.get('free') is not None and certs['field'] is None \
                and algebras['field'][2] is not None and algebras['field'][2] < algebras['free'][2]:
            # I_free strictly inside I_field: the closures for D < p are those of I_free, so
            # LFD(F_field) >= p (README, Remark 3); no further degree can certify I_field below p.
            certs['field'] = dict(kind='strict-containment-lower-bound', lower_bound=p,
                                  free_dimension=algebras['free'][2], field_dimension=algebras['field'][2])
        if all(certs[key] is not None or algebras[key][2] is None for key in algebras):
            break
    results = {}
    for key, (kind, data, dim) in algebras.items():
        if unit is not None:
            assert dim == 0, (key, dim)
            results[key] = dict(status='exact-unit', exact_lfd=unit, lfd_lower_bound=unit,
                                lfd_upper_bound=unit, certificate=dict(kind='unit', degree=unit))
        elif certs[key] is not None and certs[key]['kind'] == 'strict-containment-lower-bound':
            results[key] = dict(status='lower-bound-p', exact_lfd=None, lfd_lower_bound=p,
                                lfd_upper_bound=None, certificate=certs[key])
        elif certs[key] is not None:
            v = certs[key]['exact_lfd']
            results[key] = dict(status='exact-proper-local', exact_lfd=v, lfd_lower_bound=v,
                                lfd_upper_bound=v, certificate=certs[key])
        else:
            results[key] = dict(status='capped', exact_lfd=None,
                                lfd_lower_bound=max(falls, default=0),
                                lfd_upper_bound=(3 * k + 1 if encoding == 'direct' else None),
                                cap_reason=reason if dim is not None else 'not zero-dimensional')
        results[key]['quotient_dimension'] = dim
    return dict(encoding=encoding, generator_hash=gen_hash, generator_degrees=gen_degrees,
                strict_falls=falls, profile=rows, results=results, engine=engine)
