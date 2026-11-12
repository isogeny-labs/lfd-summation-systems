#!/usr/bin/env python3
"""Verifier step 2: re-derive, with own code, the exact point-set claims made from the
producer's p = 193, k = 6 plane sweep (coverage, class counts, the C1 cubic's zero set,
the degree-18 rank of C2, the per-line counts) and the F_257 / F_449 cubic fits; then
write the plane-point sample that ic_bounds_locus_verify.py recomputes (`--sample`)."""
import argparse
import collections
import glob
import json
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "experiments", "ic-curve-validation", "code"))
import lfd_engine  # noqa: E402

EXP = os.path.join(ROOT, "experiments", "ic-bounds-locus")
CUBIC193 = {(3, 0): 1, (2, 1): 30, (1, 2): 72, (0, 3): 50, (2, 0): 96, (1, 1): 131, (0, 2): 8,
            (1, 0): 166, (0, 1): 162, (0, 0): 15}


def sig(r):
    if r.get("saturated4"):
        return ("sat",)
    return (r["V3"], r["low3"], r["V4"], r["low4"], r["one4"])


def ev(poly, A, B, p):
    return sum(c * pow(A, i, p) * pow(B, j, p) for (i, j), c in poly.items()) % p


def monos(d):
    return [(i, j) for s in range(d + 1) for i in range(s + 1) for j in [s - i]]


def nullspace_cubic(pts, p):
    """all cubics through pts: returns (rank, a basis vector if nullity 1)."""
    ms = monos(3)
    M = np.array([[pow(A, i, p) * pow(B, j, p) % p for (i, j) in ms] for A, B in pts], dtype=np.int64)
    rk = lfd_engine.rank_mod(M, p)
    sol = None
    if rk == len(ms) - 1:
        # solve with the A^3 coefficient = 1 by brute elimination (own small Gauss)
        Mat = [[int(x) for x in r] for r in M]
        n = len(ms)
        # unknowns c_m; normalize index of (3,0)
        i30 = ms.index((3, 0))
        cols = [c for c in range(n) if c != i30]
        rows = [[r[c] for c in cols] + [(-r[i30]) % p] for r in Mat]
        # gaussian elimination
        piv = []
        R = rows
        ri = 0
        for c in range(len(cols)):
            pr = next((i for i in range(ri, len(R)) if R[i][c] % p), None)
            if pr is None:
                continue
            R[ri], R[pr] = R[pr], R[ri]
            inv = pow(R[ri][c], -1, p)
            R[ri] = [x * inv % p for x in R[ri]]
            for i in range(len(R)):
                if i != ri and R[i][c]:
                    f = R[i][c]
                    R[i] = [(x - f * y) % p for x, y in zip(R[i], R[ri])]
            piv.append(c)
            ri += 1
        sol = {(3, 0): 1}
        for i, c in enumerate(piv):
            sol[ms[cols[c]]] = R[i][-1]
    return rk, sol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(EXP, "verify", "plane-checks.json"))
    ap.add_argument("--sample-out", default=os.path.join(EXP, "verify", "plane-sample.json"))
    args = ap.parse_args()
    p = 193
    rows = []
    for fn in sorted(glob.glob(os.path.join(EXP, "plane", "p193-k6-plane-c*.jsonl"))):
        rows += [json.loads(l) for l in open(fn)]
    pts = {}
    dup = 0
    for r in rows:
        key = (r["A"], r["B"])
        dup += key in pts
        pts[key] = r
    out = dict(rows=len(rows), distinct=len(pts), duplicates=dup, complete=len(pts) == p * p,
               all_t1=all(r["t"] == 1 for r in rows))
    cls = collections.Counter(sig(r) for r in pts.values())
    out["classes"] = [[list(k), v] for k, v in cls.most_common()]
    gen = (13, 1, 605, 111, False)
    # cubic zero set
    zeros = [(A, B) for A in range(p) for B in range(p) if ev(CUBIC193, A, B, p) == 0]
    zc = collections.Counter(sig(pts[z]) for z in zeros)
    c602 = [k for k, r in pts.items() if sig(r) == (13, 1, 602, 111, False)]
    out["cubic193"] = dict(zeros=len(zeros), zero_classes=[[list(k), v] for k, v in zc.items()],
                           class602=len(c602), class602_all_on_cubic=all(ev(CUBIC193, A, B, p) == 0 for A, B in c602),
                           special_on_cubic={str(z): [ev(CUBIC193, z[0], z[1], p) == 0, list(sig(pts[z]))]
                                             for z in [(190, 2), (190, 170), (162, 58), (32, 103)]},
                           generic_on_cubic=sum(1 for z in zeros if sig(pts[z]) == gen))
    # own refit of the cubic on the 602 class
    rk, sol = nullspace_cubic(c602, p)
    out["cubic193"]["refit_rank_of_10_monomials"] = rk
    out["cubic193"]["refit_equals_claimed"] = sol == CUBIC193 if sol else None
    # C2
    c2 = sorted(k for k, r in pts.items() if sig(r) == ("sat",))
    out["C2"] = dict(points=len(c2))
    for d in (17, 18, 19):
        ms = monos(d)
        M = np.array([[pow(A, i, p) * pow(B, j, p) % p for (i, j) in ms] for A, B in c2], dtype=np.int64)
        out["C2"][f"rank_deg{d}"] = [lfd_engine.rank_mod(M, p), len(ms)]
    per_line = collections.Counter(A for A, B in c2)
    dist = collections.Counter(per_line.get(A, 0) for A in range(p))
    out["C2"]["points_per_line_A"] = dict(sorted(dist.items()))
    out["C2"]["V4_plus_Z"] = sorted({r["V4"] + r["Z"] for k, r in pts.items() if sig(r) == ("sat",)})
    # lines
    out["line_1+A+B"] = collections.Counter(str(list(sig(pts[(A, (-1 - A) % p)]))) for A in range(p))
    out["line_A=-3"] = collections.Counter(str(list(sig(pts[(p - 3, B)]))) for B in range(p))
    # other primes: refit cubics on the k = 6 V4 = 602 class
    other = {}
    for q, files in ((257, ["p257-k6-lines-A1-20.jsonl", "p257-k6-B0.jsonl"]), (449, ["p449-k6-lines-A1-6.jsonl", "p449-k6-B0.jsonl"]),
                     (577, ["p577-k6-lines-A1-6.jsonl"])):
        P = {}
        for fn in files:
            for l in open(os.path.join(EXP, "plane", fn)):
                r = json.loads(l)
                P[(r["A"], r["B"])] = r
        if q == 257:
            for r in json.load(open(os.path.join(EXP, "plane", "p257-crossk.json"))):
                P.setdefault((r["A"], r["B"]), dict(r["k6"], A=r["A"], B=r["B"]))
        c1 = sorted(k for k, r in P.items() if (r["V3"], r["V4"]) == (13, 602))
        rk, sol = nullspace_cubic(c1, q) if len(c1) >= 9 else (None, None)
        onc = None
        if sol:
            onc = collections.Counter(str(list(sig(r)) if r.get("saturated4") is not None else (r["V3"], r["V4"]))
                                      for k, r in P.items() if ev(sol, k[0], k[1], q) == 0)
        other[q] = dict(points_computed=len(P), class602=len(c1), cubic_rank=rk,
                        cubic={f"{i},{j}": c for (i, j), c in sol.items()} if sol else None, points_on_cubic_by_class=onc)
    out["other_primes"] = other
    json.dump(out, open(args.out, "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:6000])
    # sample for recomputation: all specials, all C2, all C1, 25 order-two, 25 tangent, 300 generic
    rng = random.Random(4242)
    by = collections.defaultdict(list)
    for k, r in pts.items():
        by[sig(r)].append(k)
    samp = []
    for s, L in by.items():
        L = sorted(L)
        if s == gen:
            samp += [("generic", x) for x in rng.sample(L, 300)]
        elif len(L) <= 190 and s in (("sat",), (13, 1, 602, 111, False)):
            samp += [(str(s), x) for x in L]
        elif len(L) > 20:
            samp += [(str(s), x) for x in rng.sample(L, 25)]
        else:
            samp += [(str(s), x) for x in L]
    json.dump(samp, open(args.sample_out, "w"))
    print("sample", len(samp))


if __name__ == "__main__":
    main()
