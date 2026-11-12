#!/usr/bin/env python3
"""Verifier step 3: machine checks of the identities used in the producer's written proofs
(order-two square, t = 0 b = 0 factorisation, weighted homogeneity, the weight-2 matrix of g
and its determinant) and an explicit instance check of the tangent-line construction
h1 = (u11 + 2t u10) g, h2 = (u21 + 2t u20) g at several (p, k, a, t) with 3t^2 + a = 0,
plus a direct check of the H-scaling map on generators (x -> x/mu, mu in H)."""
import json
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ic_bounds_locus_verify import generator, subgroup  # own code

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "experiments", "ic-bounds-locus", "verify", "proof-checks.json")


def mulvar(poly, k, i, j, p):
    N = 1 << k
    lo = N - 1
    out = {}
    for m, c in poly.items():
        e = [m & lo, (m >> k) & lo]
        e[i] = (e[i] + (1 << j)) % N
        nm = e[0] | (e[1] << k)
        out[nm] = (out.get(nm, 0) + c) % p
    return {m: c for m, c in out.items() if c}


def add(f, g, p, s=1):
    out = dict(f)
    for m, c in g.items():
        out[m] = (out.get(m, 0) + s * c) % p
    return {m: c for m, c in out.items() if c}


def scal(f, s, p):
    return {m: c * s % p for m, c in f.items() if c * s % p}


def wt(f):
    return max((bin(m).count("1") for m in f), default=-1)


def rank3(polys, p):
    import numpy as np
    sys.path.insert(0, os.path.join(ROOT, "experiments", "ic-curve-validation", "code"))
    import lfd_engine
    keys = sorted(set().union(*polys))
    M = np.array([[f.get(m, 0) for m in keys] for f in polys], dtype=np.int64)
    return lfd_engine.rank_mod(M, p)


def main():
    res = {}
    x1, x2, t, a, b, mu = sp.symbols("x1 x2 t a b mu")
    S3 = (x1 - x2) ** 2 * t ** 2 - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)
    q = x1 * x2 - t * (x1 + x2) - a - 2 * t ** 2
    res["order_two_square"] = sp.expand(S3.subs(b, -t ** 3 - a * t) - q ** 2) == 0
    res["t0_b0_factor"] = str(sp.factor(S3.subs({t: 0, b: 0})))
    res["weighted_homogeneous"] = sp.expand(S3.subs({x1: mu * x1, x2: mu * x2, t: mu * t, a: mu ** 2 * a, b: mu ** 3 * b},
                                                    simultaneous=True) - mu ** 4 * S3) == 0
    d = dict(sp.Poly(sp.expand(S3), x1, x2).terms())
    M = sp.Matrix([[d[(2, 2)], d[(2, 1)]], [d[(1, 2)], d[(1, 1)]]])
    res["weight2_matrix"] = str(M)
    res["weight2_det"] = str(sp.factor(M.det()))
    # tangent instances
    inst = []
    for p, k in ((97, 3), (97, 5), (193, 5), (193, 6), (257, 7), (257, 8)):
        rng = random.Random(p * 100 + k)
        for _ in range(3):
            tt = rng.randrange(1, p)
            aa = (-3 * tt * tt) % p
            bb = rng.randrange(p)
            g = generator(p, k, aa, bb, tt)
            u = lambda i, j, f: mulvar(f, k, i, j, p)  # noqa: E731
            h1 = add(u(0, 1, g), scal(u(0, 0, g), 2 * tt, p), p)
            h2 = add(u(1, 1, g), scal(u(1, 0, g), 2 * tt, p), p)
            m1221 = (1 << 2) | (1 << (k + 1))   # u12 u21
            m2211 = (1 << 1) | (1 << (k + 2))   # u11 u22
            inst.append(dict(p=p, k=k, a=aa, b=bb, t=tt, wt_g=wt(g), wt_h1=wt(h1), wt_h2=wt(h2),
                             h1_u12u21=h1.get(m1221, 0), g_u12u21=g.get(m1221, 0), h2_u12u21=h2.get(m1221, 0),
                             h2_u11u22=h2.get(m2211, 0), g_u11u22=g.get(m2211, 0), h1_u11u22=h1.get(m2211, 0),
                             rank=rank3([g, h1, h2], p)))
    res["tangent_instances"] = inst
    res["tangent_all_ok"] = all(r["wt_h1"] <= 2 and r["wt_h2"] <= 2 and r["rank"] == 3 for r in inst)
    # control: off the tangent line h1 has weight 3
    ctrl = []
    for p, k in ((193, 6), (257, 7)):
        tt, aa, bb = 5, 11, 7
        g = generator(p, k, aa, bb, tt)
        h1 = add(mulvar(g, k, 0, 1, p), scal(mulvar(g, k, 0, 0, p), 2 * tt, p), p)
        ctrl.append(dict(p=p, k=k, a=aa, b=bb, t=tt, tangent=(3 * tt * tt + aa) % p == 0, wt_h1=wt(h1)))
    res["tangent_control_off_line"] = ctrl
    # H-scaling: generator of (mu^2 a, mu^3 b, mu t) equals mu^4 * phi(g_{a,b,t}), phi(x^e) = mu^-e x^e
    hs = []
    for p, k in ((193, 6), (257, 7)):
        H = subgroup(p, k)
        rng = random.Random(p)
        N = 1 << k
        lo = N - 1
        for _ in range(5):
            m_ = rng.choice(H)
            aa, bb, tt = rng.randrange(p), rng.randrange(p), rng.randrange(1, p)
            g = generator(p, k, aa, bb, tt)
            inv = pow(m_, -1, p)
            phi = {m: c * pow(inv, (m & lo) + ((m >> k) & lo), p) * pow(m_, 4, p) % p for m, c in g.items()}
            g2 = generator(p, k, aa * m_ ** 2 % p, bb * m_ ** 3 % p, tt * m_ % p)
            hs.append(dict(p=p, k=k, mu=m_, a=aa, b=bb, t=tt, equal=phi == g2))
    res["H_scaling_generator"] = hs
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("tangent_instances",)}, indent=1))
    print("tangent instances:", [(r["p"], r["k"], r["wt_h1"], r["wt_h2"], r["rank"]) for r in inst])


if __name__ == "__main__":
    main()
