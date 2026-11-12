#!/usr/bin/env python3
"""Independent verification of ITEM 2 (exceptional locus) of the bounds program.

Written by the verifier; imports NONE of the producer's modules (csd/ic_bounds_locus*.py,
scripts/ic_bounds_locus_*.py) and none of the repository closure engines
(ic_fall_certificates, ic_chain_shift_closure, ic_chain_block_closure, ic_curve_sweep).
Linear algebra: the `Ech` reduced-echelon class of the independent validator engine
experiments/ic-curve-validation/code/lfd_engine.py (and its ordinary-ring `closure` for
a second, ordinary-ring cross-check).  Everything else (generator, chain-quotient closure,
zero set, evaluation ranks, certificates) is own code here.

System: g = NF(S_3(x1, x2, t)) in F_p[x1,x2]/(x1^N - 1, x2^N - 1), N = 2^k, written in the
chain variables u_{i,j} = x_i^(2^j): a chain monomial is a 2k-bit mask (bits 0..k-1 = binary
digits of the exponent e1 mod N, bits k..2k-1 those of e2); its weight is its popcount.
Multiplication by u_{i,j} adds 2^j to e_i mod N.

V_D (modulo chains) = least subspace of span(M_<=D) containing g (weight 2) and closed under
h -> u*h for every chain variable u whenever weight(h) <= D-1.  The closure loop below stops
only when the shifts of a BASIS of V_D cap M_<=D-1 all reduce to zero against V_D, so every
value it returns is a verified closed fixed point (not a stall).

I_D = kernel of evaluation of M_<=D on Z (the chain algebra is a product of copies of F_p, so
every ideal is the vanishing ideal of its zero set Z inside H x H).
Certificate of LFD <= D0 (local-to-global, as in lfd_engine.exact_lfd): dim V_D0 = dim I_D0
and evaluation of M_<D0 onto F_p^Z; unit ideal: 1 in V_D0.  Exact LFD = largest D <= D0 with
V_{D-1} != I_{D-1} ... equivalently the last strict fall.
"""
import argparse
import hashlib
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "experiments", "ic-curve-validation", "code"))
import lfd_engine  # noqa: E402  (independent engine: Ech, closure, monomials)


# ---------------------------------------------------------------------------------------
# field / curve helpers (own)
# ---------------------------------------------------------------------------------------

def subgroup(p, k):
    N = 1 << k
    assert (p - 1) % N == 0
    H = sorted(x for x in range(1, p) if pow(x, N, p) == 1)
    assert len(H) == N
    return H


def is_square(v, p):
    v %= p
    return v == 0 or pow(v, (p - 1) // 2, p) == 1


def s3_coeffs(p, a, b, t):
    """S_3(x1, x2, t) = (x1-x2)^2 t^2 - 2((x1+x2)(x1x2+a) + 2b) t + (x1x2-a)^2 - 4b(x1+x2),
    as {(e1, e2): c mod p}, expanded by hand (own code)."""
    c = {}

    def add(e1, e2, v):
        c[(e1, e2)] = (c.get((e1, e2), 0) + v) % p
    # (x1-x2)^2 t^2
    add(2, 0, t * t); add(1, 1, -2 * t * t); add(0, 2, t * t)
    # -2t[(x1+x2)(x1x2+a) + 2b]
    add(2, 1, -2 * t); add(1, 2, -2 * t); add(1, 0, -2 * a * t); add(0, 1, -2 * a * t); add(0, 0, -4 * b * t)
    # (x1x2-a)^2
    add(2, 2, 1); add(1, 1, -2 * a); add(0, 0, a * a)
    # -4b(x1+x2)
    add(1, 0, -4 * b); add(0, 1, -4 * b)
    return {e: v for e, v in c.items() if v}


def s3_value(p, a, b, x1, x2, x3):
    return ((x1 - x2) ** 2 * x3 ** 2 - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * x3
            + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p


def generator(p, k, a, b, t):
    N = 1 << k
    g = {}
    for (e1, e2), v in s3_coeffs(p, a, b, t).items():
        m = (e1 % N) | ((e2 % N) << k)
        g[m] = (g.get(m, 0) + v) % p
    return {m: v for m, v in g.items() if v}


def gen_hash(g):
    s = json.dumps(sorted([m, v] for m, v in g.items()))
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def zero_set(p, k, a, b, t):
    H = subgroup(p, k)
    return [(x1, x2) for x1 in H for x2 in H if s3_value(p, a, b, x1, x2, t) == 0]


# ---------------------------------------------------------------------------------------
# chain-quotient monomial space (own)
# ---------------------------------------------------------------------------------------

class ChainSpace:
    def __init__(self, k, D):
        self.k, self.D, self.v = k, D, 2 * k
        N = 1 << k
        masks = [m for m in range(1 << self.v) if bin(m).count("1") <= D]
        masks.sort(key=lambda m: (-bin(m).count("1"), m))   # weight descending
        self.cols = masks
        self.n = len(masks)
        self.idx = {m: i for i, m in enumerate(masks)}
        self.w = np.array([bin(m).count("1") for m in masks])
        low = np.flatnonzero(self.w <= D - 1)
        self.low_start = int(low[0]) if low.size else self.n
        lo = N - 1
        self.mul = []
        for i in range(2):
            for j in range(k):
                arr = np.empty(self.n - self.low_start, dtype=np.int64)
                for c in range(self.low_start, self.n):
                    m = masks[c]
                    e = [m & lo, (m >> k) & lo]
                    e[i] = (e[i] + (1 << j)) % N
                    nm = e[0] | (e[1] << k)
                    arr[c - self.low_start] = self.idx[nm]   # weight(nm) <= weight(m)+1 <= D
                self.mul.append(arr)

    def shift(self, rows, var):
        out = np.zeros((rows.shape[0], self.n), dtype=rows.dtype)
        # distinct source columns map to distinct targets (multiplication by a unit x^(2^j))
        out[:, self.mul[var]] = rows[:, self.low_start:]
        return out

    def row(self, g, p):
        r = np.zeros(self.n)
        for m, c in g.items():
            r[self.idx[m]] = c % p
        return r


def chain_closure(g, k, D, p, sp=None):
    sp = sp or ChainSpace(k, D)
    V = lfd_engine.Ech(sp.n, p)
    done = lfd_engine.Ech(sp.n, p, cap=max(1, sp.n - sp.low_start))
    if max(bin(m).count("1") for m in g) <= D:
        V.add(sp.row(g, p).reshape(1, -1))
    rounds = 0
    while True:
        low = V.E[V.piv >= sp.low_start]
        if low.shape[0] == 0:
            break
        resid = done.reduce(low)
        resid = resid[np.any(resid != 0, axis=1)]
        if resid.shape[0] == 0:
            break
        new = done.add(resid)
        rounds += 1
        for blk in new:
            for var in range(2 * k):
                V.add(sp.shift(blk, var))
    # explicit closedness re-check: every shift of every basis row of V cap M_<=D-1 reduces to 0
    low = V.E[V.piv >= sp.low_start]
    closed = True
    for var in range(2 * k):
        if low.shape[0] and np.any(V.reduce(sp.shift(low, var))):
            closed = False
    lowdim = int(np.sum(V.piv >= sp.low_start))
    one = (sp.n - 1) in set(V.piv.tolist())   # last column = mask 0 = the constant
    assert sp.cols[-1] == 0
    return dict(dim=V.rk, low=lowdim, one=one, closed=closed, rounds=rounds)


def eval_ranks(p, k, Z, Dmax):
    """rank of evaluation of M_<=d on Z, d = 0..Dmax."""
    if not Z:
        return [0] * (Dmax + 1)
    N = 1 << k
    lo = N - 1
    sp = ChainSpace(k, Dmax)
    M = np.empty((len(Z), sp.n), dtype=np.int64)
    for r, (x1, x2) in enumerate(Z):
        p1 = [pow(x1, e, p) for e in range(N)]
        p2 = [pow(x2, e, p) for e in range(N)]
        M[r] = [p1[m & lo] * p2[(m >> k) & lo] % p for m in sp.cols]
    out = []
    for d in range(Dmax + 1):
        out.append(lfd_engine.rank_mod(M[:, sp.w <= d], p))
    return out


def nmono(k, d):
    from math import comb
    return sum(comb(2 * k, i) for i in range(d + 1))


def analyse(p, k, a, b, t, Dmax, spaces=None):
    t0 = time.time()
    g = generator(p, k, a, b, t)
    Z = zero_set(p, k, a, b, t)
    er = eval_ranks(p, k, Z, Dmax)
    per = {}
    cert = None
    for D in range(2, Dmax + 1):
        sp = spaces.get(D) if spaces else None
        r = chain_closure(g, k, D, p, sp)
        dimI = nmono(k, D) - er[D]
        assert r["dim"] <= dimI, "soundness violated"
        r["dimI"] = dimI
        per[D] = r
        if r["one"]:
            cert = ("unit", D)
            break
        if r["dim"] == dimI and Z and er[D - 1] == len(Z):
            cert = ("local", D)
            break
    falls = [D for D in per if D - 1 in per and per[D]["low"] > per[D - 1]["dim"]]
    if 3 in per and per[3]["low"] > 1:
        falls = sorted(set(falls) | {3})   # V_2 = span(g), dim 1
    lfd = None
    if cert:
        if cert[0] == "unit":
            lfd = cert[1]
        else:
            neq = [D for D in per if per[D]["dim"] != per[D]["dimI"]]
            lfd = max(neq) + 1 if neq else 2
    return dict(prime=p, k=k, a=a, b=b, t=t, Z=len(Z), gen_hash=gen_hash(g), generator=sorted([m, c] for m, c in g.items()),
                per={str(D): v for D, v in per.items()}, eval_ranks=er, cert=cert, falls=falls, lfd=lfd,
                seconds=round(time.time() - t0, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--dmax", type=int, default=5)
    ap.add_argument("--targets", required=True, help="a,b,t;a,b,t;...")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    spaces = {D: ChainSpace(args.k, D) for D in range(2, args.dmax + 1)}
    with open(args.out, "a") as f:
        for tok in args.targets.split(";"):
            a, b, t = map(int, tok.split(","))
            res = analyse(args.prime, args.k, a, b, t, args.dmax, spaces)
            f.write(json.dumps(res) + "\n")
            f.flush()
            print(a, b, t, res["Z"], {D: (v["dim"], v["low"], v["one"], v["closed"]) for D, v in res["per"].items()},
                  res["cert"], res["falls"], res["lfd"], res["seconds"], flush=True)


if __name__ == "__main__":
    main()
