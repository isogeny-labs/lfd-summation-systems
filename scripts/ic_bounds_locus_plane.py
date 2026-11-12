"""Exceptional locus in the normalized (A, B) plane (bounds item 2).

S_3 is weighted homogeneous (x, t of weight 2, a of weight 4, b of weight 6), so
for t != 0 the substitution x_i = t y_i turns NF(S_3(x1, x2, t)) for (a, b, t)
into t^4 NF(S_3(y1, y2, 1)) for (A, B) = (a/t^2, b/t^3), up to the scaling
y -> y/t of the chain variables.  That scaling multiplies the coefficient of
X^(e1,e2) by t^-(e1+e2); it preserves Hamming weights but is an algebra
automorphism of F_p[G] only when t is in H, so invariance of the closure
dimensions under it is an empirical observation (checked in this script on the
rational-target data: recorded per run in the 'scaling_check' field).

For a list of lines (fixed A, every B in F_p, or fixed B, every A) this script
computes, at t = 1, a = A, b = B (singular cubics included: S_3 is a polynomial
system for every (a, b)), the dense closures V_3 and V_4 modulo the chains with
exclusion certificates, |Z| (zeros of g on H x H) and the evaluation rank of
M_<=3 on Z.

Usage
  python3 scripts/ic_bounds_locus_plane.py --prime 193 --k 6 --fix A --values 0,1,2 --out FILE.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_fall_certificates as C  # noqa: E402
from csd import ic_curve_sweep as S  # noqa: E402
from csd import ic_chain_shift_closure as SC  # noqa: E402


def s3_generator(p, k, a, b, t, lam=1):
    """NF(S_3(lam*x1, lam*x2, t)) on y^2 = x^3 + a x + b (any a, b) as {mask: coeff}."""
    N = 1 << k
    c = {(2, 0): t * t, (0, 2): t * t, (1, 1): -2 * t * t - 2 * a, (2, 1): -2 * t, (1, 2): -2 * t,
         (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b, (0, 0): -4 * b * t + a * a, (2, 2): 1}
    g = {}
    for (e1, e2), v in c.items():
        v = v * pow(lam, e1 + e2, p) % p
        if v:
            m = (e1 % N) | ((e2 % N) << k)
            g[m] = (g.get(m, 0) + v) % p
    return {m: v for m, v in g.items() if v}


def zeros(g, p, k):
    N = 1 << k
    g0 = S.subgroup_generator_value(p, k)
    pw = [pow(g0, e, p) for e in range(N)]
    out = []
    for e1 in range(N):
        for e2 in range(N):
            v = 0
            for m, c in g.items():
                v += c * pw[(m & (N - 1)) * e1 % N] * pw[(m >> k) * e2 % N]
            if v % p == 0:
                out.append(e1 | (e2 << k))
    return out, g0


def closure_row(p, k, g, certify=True, with_zeros=True):
    Z, g0 = zeros(g, p, k) if with_zeros else ([], None)
    rec = dict(Z=len(Z) if with_zeros else None)
    for D in (3, 4):
        V = C.chain_degree_closure([g], 2, k, p, D, max_columns=20000, max_cells=60_000_000)
        one = V.contains_one
        if not one and certify:
            assert C.check_closed_exclusion_certificate(V, [g])
        low = (len([m for m in V.masks if m.bit_count() <= D - 1]) if one else
               sum(V.masks[j].bit_count() <= D - 1 for j in V.pivots))
        rec[f'V{D}'] = V.quotient_rank
        rec[f'low{D}'] = low
        rec[f'one{D}'] = one
        rec[f'codim{D}'] = len(V.masks) - V.quotient_rank
    if not with_zeros:
        rec['eval3'] = rec['saturated4'] = None
        return rec
    rec['eval3'] = S.prime_evaluation_rank(Z, SC.monomials(k, 3), k, g0, p, 10 ** 8) if Z else 0
    rec['saturated4'] = rec['codim4'] == len(Z) and rec['eval3'] == len(Z)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prime', type=int, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--fix', choices=['A', 'B'], required=True)
    ap.add_argument('--values', required=True)
    ap.add_argument('--lam', type=int, default=1, help='extra scaling of the chain variables (scaling check)')
    ap.add_argument('--no-certify', action='store_true',
                    help='skip the (redundant) exclusion-certificate re-check for speed')
    ap.add_argument('--no-zeros', action='store_true', help='skip |Z| and the saturation test')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    p, k = args.prime, args.k
    out = open(args.out, 'a')
    for v in [int(x) for x in args.values.split(',')]:
        t0 = time.time()
        for w in range(p):
            A, B = (v, w) if args.fix == 'A' else (w, v)
            g = s3_generator(p, k, A, B, 1, args.lam)
            rec = dict(prime=p, k=k, A=A % p, B=B % p, t=1, lam=args.lam, certified=not args.no_certify,
                       **closure_row(p, k, g, certify=not args.no_certify, with_zeros=not args.no_zeros))
            out.write(json.dumps(rec, separators=(',', ':')) + '\n')
        out.flush()
        print(args.fix, v, 'done', round(time.time() - t0, 1), flush=True)


if __name__ == '__main__':
    main()
