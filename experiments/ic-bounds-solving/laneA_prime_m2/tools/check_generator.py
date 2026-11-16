"""Checks: (1) README g == NF(S_3) numerically, (2) matches lfd_engine zero set, (3) S_3 vanishes on
x-coords of P1+P2+P3=O on a small curve. [pcc:unverified]"""
import random, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
rng = random.Random(1)
# (1) polynomial identity at random points for several p
for p in (97, 12289, 2013265921 % (2**31) if False else 40961):
    for _ in range(2000):
        a, b = random_curve(p, rng); t = rng.randrange(p); x1, x2 = rng.randrange(p), rng.randrange(p)
        assert S3(x1, x2, t, a, b, p) == g_eval(a, b, t, p, x1, x2)
print("g == S_3 as polynomial: ok")
# (2) lfd_engine.prime_points zero set vs direct
for (p, k) in ((97, 3), (113, 4), (193, 3)):
    for _ in range(5):
        a, b = random_curve(p, rng); t = rng.randrange(p)
        items = []
        for (e1, e2), c in gen_coeffs(a, b, t, p).items():
            items.append((e1 | (e2 << k), c))
        inputs, G = LE.prime_inputs(items, k, p)
        pts = LE.prime_points(G, k, p)
        z = zeros_count(a, b, t, p, k)
        assert len(pts) == z, (len(pts), z)
print("prime_inputs/prime_points agree with direct zero count: ok")
# (3) Semaev property on y^2 = x^3+ax+b mod 97
p = 97; a, b = 2, 3
assert (4*a**3 + 27*b*b) % p
pts = [(x, y) for x in range(p) for y in range(p) if (y*y - x**3 - a*x - b) % p == 0]
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    lam = (3*x1*x1 + a) * pow(2*y1, -1, p) % p if P == Q else (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam*lam - x1 - x2) % p
    return (x3, (lam*(x1 - x3) - y1) % p)
cnt = 0
for P in pts:
    for Q in pts:
        R = add(P, Q)
        if R is None: continue
        assert S3(P[0], Q[0], R[0], a, b, p) == 0; cnt += 1
print("Semaev S_3(x(P),x(Q),x(P+Q))=0 on", cnt, "pairs: ok")
