#!/usr/bin/env python3
"""Verifier (independent of csd/ic_bounds_m3*.py): symbolic facts about S_4.

S_3(u,v,w) = (u-v)^2 w^2 - 2((u+v)(uv+a)+2b) w + (uv-a)^2 - 4b(u+v)   (Semaev)
S_4(x1,x2,x3,x4) = Res_X(S_3(x1,x2,X), S_3(x3,x4,X)).
Checks: symmetry, degree per variable, total degree, coefficient of x1^3x2^3x3^3 at x4=t,
coefficient of x3^4 at x4=t, and that S_4 vanishes at x-coordinates of random decompositions
P1+P2+P3 = R computed with an independent affine group law written here.
Writes experiments/ic-bounds-m3/verify/s4_coeffs.json (exact integer coefficients in a,b,x1..x4)
and verify/symbolic_verify.json.
"""
import json, random, hashlib, itertools
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parents[1] / 'experiments' / 'ic-bounds-m3' / 'verify'
a, b, x1, x2, x3, x4, X, t = sp.symbols('a b x1 x2 x3 x4 X t')


def S3(u, v, w):
    return (u - v)**2 * w**2 - 2 * ((u + v) * (u * v + a) + 2 * b) * w + (u * v - a)**2 - 4 * b * (u + v)


def main():
    rep = {}
    S4 = sp.Poly(sp.resultant(S3(x1, x2, X), S3(x3, x4, X), X), a, b, x1, x2, x3, x4)
    gens = (x1, x2, x3, x4)
    P = sp.Poly(S4.as_expr(), *gens)
    rep['deg_per_var'] = [P.degree(g) for g in gens]
    rep['total_degree_in_x'] = P.total_degree()
    rep['terms_in_x_generic'] = len(P.terms())
    e = S4.as_expr()
    sym = all(sp.expand(e - e.subs({gens[i]: gens[j], gens[j]: gens[i]}, simultaneous=True)) == 0
              for i, j in [(0, 1), (0, 2), (0, 3), (2, 3)])
    rep['symmetric'] = bool(sym)
    Et = sp.Poly(e.subs(x4, t), x1, x2, x3)
    c333 = sp.factor(Et.coeff_monomial(x1**3 * x2**3 * x3**3))
    rep['c333'] = str(c333)
    rep['c333_identity'] = bool(sp.expand(c333 + 8 * (5 * t**3 + 3 * a * t + 6 * b)) == 0)
    # coefficient of x3^4 (as polynomial in x3 with coefficients in Z[a,b,x1,x2,t])
    c34 = sp.Poly(e.subs(x4, t), x3).coeff_monomial(x3**4)
    rep['x3^4_identity'] = bool(sp.expand(c34 - S3(x1, x2, t)**2) == 0)
    # exponent range check per monomial of Et (needed for the 'unique degree-6 mask' argument)
    rep['max_exp_per_var_at_x4=t'] = [Et.degree(g) for g in (x1, x2, x3)]
    # dump coefficients
    coeffs = [[list(m), int(c)] for m, c in S4.terms()]
    blob = json.dumps(coeffs, sort_keys=True).encode()
    (OUT / 's4_coeffs.json').write_text(json.dumps(dict(order=['a', 'b', 'x1', 'x2', 'x3', 'x4'], terms=coeffs)))
    rep['s4_coeffs_sha256'] = hashlib.sha256(blob).hexdigest()
    rep['s4_terms_in_abx'] = len(coeffs)
    # numeric decompositions with an independent group law
    rng = random.Random(7351)
    f = sp.lambdify((a, b, x1, x2, x3, x4), e, 'math')
    ok = tot = 0
    for p in [10007, 65537, 1000003, 999983]:
        for _ in range(25):
            while True:
                A, B = rng.randrange(p), rng.randrange(p)
                if (4 * A**3 + 27 * B**2) % p: break
            pts = []
            while len(pts) < 3:
                xx = rng.randrange(p); r = (xx**3 + A * xx + B) % p
                if r == 0 or pow(r, (p - 1) // 2, p) != 1: continue
                yy = sqrt_mod(r, p); pts.append((xx, yy if rng.random() < .5 else (-yy) % p))
            R = add(add(pts[0], pts[1], A, p), pts[2], A, p)
            if R is None: continue
            tot += 1
            val = int(eval_s4(coeffs, A, B, [pts[0][0], pts[1][0], pts[2][0], R[0]], p))
            ok += (val == 0)
            # negative control: perturb target
            val2 = eval_s4(coeffs, A, B, [pts[0][0], pts[1][0], pts[2][0], (R[0] + 1) % p], p)
            rep.setdefault('neg_control_nonzero', 0)
            rep['neg_control_nonzero'] += (val2 != 0)
    rep['decompositions_checked'] = tot
    rep['decompositions_vanish'] = ok
    (OUT / 'symbolic_verify.json').write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


def eval_s4(coeffs, A, B, xs, p):
    s = 0
    for (ea, eb, e1, e2, e3, e4), c in coeffs:
        s += c * pow(A, ea, p) * pow(B, eb, p) * pow(xs[0], e1, p) * pow(xs[1], e2, p) * pow(xs[2], e3, p) * pow(xs[3], e4, p)
    return s % p


def sqrt_mod(r, p):
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, tt, R = s, pow(z, q, p), pow(r, q, p), pow(r, (q + 1) // 2, p)
    while tt != 1:
        i, t2 = 0, tt
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        bb = pow(c, 1 << (m - i - 1), p); m, c = i, bb * bb % p; tt = tt * c % p; R = R * bb % p
    return R


def add(P, Q, A, p):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0: return None
    if P == Q: lam = (3 * P[0]**2 + A) * pow(2 * P[1], -1, p) % p
    else: lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x = (lam * lam - P[0] - Q[0]) % p
    return (x, (lam * (P[0] - x) - P[1]) % p)


if __name__ == '__main__':
    main()
