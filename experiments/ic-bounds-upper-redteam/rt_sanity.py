"""Sanity of the red-team builder: (1) S_3 is the summation polynomial (brute force over curve
points); (2) Lemma 2.2 by brute force: no fibre polynomial vanishes identically, all alpha in F_p;
(3) the variable-wise engine closure equals the literal Lean rule (closure_naive) on edge cases."""
import sys, itertools, time
from rt_lib import *

def curve_points(p, a, b):
    return [(x, y) for x in range(p) for y in range(p) if (y * y - x ** 3 - a * x - b) % p == 0]

def add(p, a, P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        l = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)

# (1) summation polynomial property
n1 = 0
for p in (5, 7, 11, 13, 17):
    for a in range(p):
        for b in range(p):
            if not nonsingular(p, a, b):
                continue
            pts = curve_points(p, a, b)
            sums = set()
            for P in pts:
                for Q in pts:
                    R = add(p, a, P, Q)
                    if R is not None:
                        sums.add((P[0], Q[0], R[0]))
            xs = sorted({P[0] for P in pts})
            for x1 in xs:
                for x2 in xs:
                    for x3 in xs:
                        z = s3_eval(p, a, b, x3, x1, x2) == 0
                        assert z == ((x1, x2, x3) in sums), (p, a, b, x1, x2, x3)
                        n1 += 1
print('S3 is the summation polynomial on rational abscissae: triples checked', n1)

# (2) fibre condition, all alpha, t in F_p, all nonsingular curves, incl. characteristic 3
n2 = 0; lowdeg = {0: 0, 1: 0, 2: 0}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
    for a in range(p):
        for b in range(p):
            if not nonsingular(p, a, b):
                continue
            for t in range(p):
                co = s3_coeffs(p, a, b, t)
                for al in range(p):
                    q = [0, 0, 0]
                    for (e1, e2), c in co.items():
                        q[e2] = (q[e2] + c * pow(al, e1, p)) % p
                    assert any(q), ('ZERO FIBRE', p, a, b, t, al)
                    lowdeg[2 if q[2] else (1 if q[1] else 0)] += 1
                    # discriminant identity disc = 16 f(al) f(t)
                    f = lambda x: (x ** 3 + a * x + b) % p
                    assert (q[1] * q[1] - 4 * q[0] * q[2] - 16 * f(al) * f(t)) % p == 0
                    n2 += 1
print('Lemma 2.2 brute force: fibres checked', n2, 'by degree', lowdeg)

# singular curves DO have a vanishing fibre (control: the check can fail)
found = None
for p in (5, 7, 13):
    for a in range(p):
        for b in range(p):
            if nonsingular(p, a, b):
                continue
            for t in range(p):
                co = s3_coeffs(p, a, b, t)
                q = [0, 0, 0]
                for (e1, e2), c in co.items():
                    q[e2] = (q[e2] + c * pow(t, e1, p)) % p
                if not any(q):
                    found = (p, a, b, t)
print('control: singular curve with identically vanishing fibre at alpha = t:', found)

# (3) engine closure vs literal Lean rule on edge cases
n3 = 0
for (p, a, b, t, k1, k2, s1, s2, raw) in [
        (13, 1, 2, 3, 2, 2, 1, 1, False), (13, 1, 2, 3, 2, 2, 2, 2, False), (5, 1, 1, 2, 2, 2, 1, 1, False),
        (17, 0, 1, 4, 2, 2, 1, 3, False), (17, 1, 0, 0, 2, 2, 1, 1, False), (17, 2, 3, 5, 1, 3, 1, 1, False),
        (17, 2, 3, 5, 3, 1, 1, 1, False), (11, 2, 1, 3, 1, 1, 1, 1, False), (13, 1, 2, 3, 2, 2, 1, 1, True),
        (17, 3, 5, 16, 1, 2, 1, 1, False), (13, 4, 4, 1, 1, 2, 1, 1, True)]:
    inputs, pts, r = build(p, a, b, t, k1, k2, s1, s2, raw)
    for D in range(0, 6):
        c1 = E.closure(inputs, r, D, p)
        c2 = E.closure_naive(inputs, r, D, p)
        assert (c1['dim'], c1['dim_low'], c1['contains_one']) == (c2['dim'], c2['dim_low'], c2['contains_one']), (p, a, b, t, k1, k2, D)
        n3 += 1
print('engine closure == literal Lean-rule closure: (instance, D) pairs', n3)
