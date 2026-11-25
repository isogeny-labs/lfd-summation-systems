"""Small-field S4 / affine-split correctness oracle and matched Boolean inputs.

No speed assumption. Raw resultant zeros and actual rational decompositions are
different objects; the PDP split includes its missing point-at-infinity branch.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import product

from . import ic_dreg as I


def rational_points(E):
    """Scan both coordinates; do not use lift_x or Semaev polynomials."""
    q = 1 << E.F.N
    return [(x, y) for x in range(q) for y in range(q) if E.on_curve((x, y))]


def decomposition_map(E, V):
    """Actual group-law oracle, keyed by affine target (x,y)."""
    points = [P for P in rational_points(E) if P[0] in set(V)]
    out = defaultdict(set)
    for P, Q, T in product(points, repeat=3):
        R = E.add(E.add(P, Q), T)
        if R is not None:
            out[R].add((P[0], Q[0], T[0]))
    return dict(out)


def root_table(F, b):
    q = 1 << F.N
    return {(x, y): sum(1 << u for u in range(q)
                        if I.s3_binary_value(F, b, x, y, u) == 0)
            for x in range(q) for y in range(q)}


def coordinate_sets(E, V, r, roots=None):
    """Raw S4, finite split, and corrected PDP sets, as ordered x triples."""
    F = E.F
    roots = roots if roots is not None else root_table(F, E.b)
    liftable = {P[0] for P in rational_points(E)}
    direct, finite, infinity = set(), set(), set()
    for x, y, z in product(V, repeat=3):
        t = (x, y, z)
        if I.s4_binary_value(F, E.b, x, y, z, r) == 0:
            direct.add(t)
        if roots[x, y] & roots[z, r]:
            finite.add(t)
        if x == y and z == r:
            infinity.add(t)
    admissible = lambda ts: {t for t in ts if all(x in liftable for x in t)}
    return dict(raw_direct=direct, raw_finite=finite,
                pdp_direct=admissible(direct), pdp_finite=admissible(finite),
                infinity=admissible(infinity),
                pdp_split=admissible(finite | infinity))


def quadratic_coefficients(F, b, x, y):
    return [F.sq(F.mul(x, y)) ^ b, F.mul(x, y), F.sq(x ^ y)]


def gcd_degree(F, a, b):
    def trim(p):
        p = list(p)
        while p and p[-1] == 0:
            p.pop()
        return p
    a, b = trim(a), trim(b)
    while b:
        while len(a) >= len(b):
            c = F.mul(a[-1], F.inv(b[-1]))
            off = len(a) - len(b)
            for j, v in enumerate(b):
                a[off + j] ^= F.mul(c, v)
            a = trim(a)
        a, b = b, a
    return len(a) - 1


def polynomial_stats(eqs, nvars):
    return dict(variables=nvars, coordinate_equations=len(eqs),
                field_equations=nvars,
                coordinate_terms=sum(map(len, eqs)),
                input_degree=max([2] + [m.bit_count() for p in eqs for m in p]))


def build_system(E, basis, r, formulation):
    """Squarefree inputs in an ordinary GF(2) ring; field equations added by solver."""
    F, ell = E.F, len(basis)
    V = I.span(basis)
    add, mul, sq = I._bp_add, lambda a, b: I._bp_mul(F, a, b), lambda a: I._bp_sq(F, a)
    const = I._bp_const
    xs = [I._bp_block(basis, i * ell) for i in range(3)]
    R = const(r)

    def s3(x, y, z):
        return add(add(sq(add(add(mul(x, y), mul(x, z)), mul(y, z))),
                       mul(mul(x, y), z)), const(E.b))

    if formulation == 'original':
        eqs = I.binary_descent(E, [basis] * 3, r)
        nvars = 3 * ell
    elif formulation == 'split':
        u = I._bp_block([1 << i for i in range(F.N)], 3 * ell)
        polys = [s3(xs[0], xs[1], u), s3(u, xs[2], R)]
        eqs = [sorted(m for m, c in p.items() if c >> j & 1)
               for p in polys for j in range(F.N)]
        eqs = [p for p in eqs if p]
        nvars = 3 * ell + F.N
    else:
        raise ValueError(formulation)

    # Truth-table -> algebraic normal form for non-liftability on V.
    liftable = {P[0] for P in rational_points(E)}
    anf = [int(x not in liftable) for x in V]
    for j in range(ell):
        for mask in range(1 << ell):
            if mask >> j & 1:
                anf[mask] ^= anf[mask ^ (1 << j)]
    for block in range(3):
        member = [m << (block * ell) for m, c in enumerate(anf) if c]
        if member:
            eqs.append(member)
    return eqs, nvars


def eval_boolean(poly, assignment):
    return sum((assignment & m) == m for m in poly) & 1


def decode(mask, basis, nblocks=3):
    V, ell = I.span(basis), len(basis)
    return tuple(V[(mask >> (ell * i)) & ((1 << ell) - 1)] for i in range(nblocks))


def infinity_solutions(E, basis, r):
    V = I.span(basis)
    if r not in V:
        return set()
    liftable = {P[0] for P in rational_points(E)}
    return {(x, x, r) for x in V if x in liftable}
