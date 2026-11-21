#!/usr/bin/env python3
"""Baseline: relation yield per linear system for LRF, m=3, monomial V, prime n.

Funnel (C-kernel Stats): targets -> lin_solvable -> lin_candidates -> found.
  yield          = found / targets           (relations per linear system)
  solvable_rate  = lin_solvable / targets
  cand_per_solv  = lin_candidates / lin_solvable
  found_per_cand = found / lin_candidates
Cross-check: lrf found-set == brute found-set on every target.
dim V2 and delta = k + dimV2 - n are recorded (the determinacy knob).
"""
import argparse, json, random, sys
sys.path.insert(0, ".")
from csd import ic_rootfind as rf
from csd import ic_dreg as I


def monomial_basis(k):
    return [1 << i for i in range(k)]  # <1, w, ..., w^(k-1)>


def run(n, k, m, ntargets, seed, a=0):
    F = I.GF2n(n, (1 << n) | _irr_tail(n))
    rng = random.Random(seed)
    # random b giving an ordinary curve; a=0 Koblitz-like unless overridden
    b = rng.randrange(1, 1 << n)
    E = I.BinaryCurve(F, a, b)
    basis = monomial_basis(k)
    S = rf.Solver(E, basis)
    dimV2 = len(rf.product_space(F, basis))
    agg = dict(targets=0, lin_solvable=0, lin_candidates=0, found=0)
    rel_total = 0
    mismatches = 0
    for _ in range(ntargets):
        R = E.random_point(rng)
        lrf_rows, stats, _ = S.lrf(m, R)
        brute_rows, _ = S.brute(m, R)
        if sorted(lrf_rows) != sorted(brute_rows):
            mismatches += 1
        for key in agg:
            agg[key] += stats[key]
        rel_total += len(lrf_rows)
    t = agg["targets"]
    sv = agg["lin_solvable"]
    cd = agg["lin_candidates"]
    fnd = agg["found"]
    return dict(
        n=n, k=k, m=m, a=a, dimV2=dimV2, delta=k + dimV2 - n,
        fb_log2=k, ntargets=ntargets, seed=seed,
        targets=t, lin_solvable=sv, lin_candidates=cd, found=fnd,
        relations=rel_total, mismatches=mismatches,
        yield_per_system=fnd / t if t else 0.0,
        solvable_rate=sv / t if t else 0.0,
        cand_per_solvable=cd / sv if sv else 0.0,
        found_per_cand=fnd / cd if cd else 0.0,
        relations_per_target_point=rel_total / ntargets,
    )


_IRR = {11: 0b101, 13: 0b11011, 17: 0b1001, 19: 0b1100011, 23: 0b100001,
        29: 0b100000101, 31: 0b1001, 37: 0b1000001100001, 41: 0b1001}
def _irr_tail(n):
    # low bits of a known irreducible trinomial/pentanomial x^n + ... + 1
    tails = {11: (1 << 2) | 1, 13: (1 << 4) | (1 << 3) | (1 << 1) | 1,
             17: (1 << 3) | 1, 19: (1 << 5) | (1 << 2) | (1 << 1) | 1,
             23: (1 << 5) | 1, 29: (1 << 2) | 1, 31: (1 << 3) | 1,
             37: (1 << 6) | (1 << 4) | (1 << 1) | 1, 41: (1 << 3) | 1}
    return tails[n]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    # sweep: prime n, k around n/3, m=3; a few m=4 rows
    plan = []
    for n in (11, 13, 17, 19, 23):
        for k in range(max(2, n // 3 - 1), n // 3 + 3):
            if 2 * k - 1 <= n:  # monomial product space fits
                plan.append((n, k, 3))
    for n, k in ((13, 4), (17, 5), (19, 5)):
        plan.append((n, k, 4))
    rows = []
    for (n, k, m) in plan:
        r = run(n, k, m, args.targets, args.seed)
        rows.append(r)
        print(json.dumps(r))
    if args.out:
        with open(args.out, "w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
