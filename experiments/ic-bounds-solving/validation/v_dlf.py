#!/usr/bin/env python3
"""Validator check 6: own implementation of the iterated degree-D closure V_{F,D} in the ORDINARY ring and the
last fall degree d_lf (unit ideal: least D with 1 in V_{F,D}); d_ff (first strict fall) recorded as well.

V_{F,D} (HKYY Def. 2.1 / CG Def. 1.5): the smallest subspace of R_{<=D} containing F cap R_{<=D} and closed under
multiplication by a variable whenever the product stays in degree <= D.  Implementation: incremental echelon
basis with pivots on the largest monomial in a degree-compatible order (degree, then lex).  In such an order
V cap R_{<=D-1} is spanned by the basis rows whose lead has degree <= D-1, so multiplying every such basis row by
every variable, inserting the product, and repeating until no insertion succeeds gives the closure.
Pivot rows are never modified after insertion, so each is multiplied exactly once.

Instances: rows of lane A's all_instances.json that carry a recorded d_lf (5 per k in {3,4,5}, chosen with my
own RNG), and optionally my own instances (--own JSONL from v_static.py).
"""
import argparse, json, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import MonoTable, SparsePivots, chain_polys, lift_generator, mul_monomial, generator_and_inverse

_TABLES = {}


def table_for(n, Dmax):
    if (n, Dmax) not in _TABLES:
        _TABLES[(n, Dmax)] = MonoTable(n, Dmax)
    return _TABLES[(n, Dmax)]


def insert_get_lead(piv, row):
    row = piv.reduce_lead(dict(row))
    if not row:
        return None
    lead = max(row)
    inv = pow(row[lead], -1, piv.p)
    piv.piv[lead] = {a: c * inv % piv.p for a, c in row.items()}
    return lead


def closure(m, k, p, g, delta, D):
    """returns (one_in_V, dim V_D, dim of V_D cap R_{<=D-1})."""
    n = m * k
    T = table_for(n, max(D, 2, delta))
    chains = chain_polys(m, k, T)
    ghat, dg = lift_generator(g, m, k, T, p)
    assert dg == delta
    F = [{a: c % p for a, c in r.items()} for r in chains] + [ghat]
    piv = SparsePivots(p)
    queue = []
    for f in F:
        if max(T.deg[i] for i in f) <= D:
            L = insert_get_lead(piv, f)
            if L is not None:
                queue.append(L)
    while queue:
        L = queue.pop()
        if T.deg[L] > D - 1:
            continue
        row = piv.piv[L]
        for v in range(n):
            L2 = insert_get_lead(piv, mul_monomial(row, T.var(v), T, p))
            if L2 is not None:
                queue.append(L2)
    dimV = len(piv.piv)
    dim_low = sum(1 for L in piv.piv if T.deg[L] <= D - 1)
    return (0 in piv.piv), dimV, dim_low


def dlf_dff(m, k, p, a, b, t, Dcap):
    info, g, ginv = generator_and_inverse(m, k, p, a, b, t)
    assert info['unit'], 'not a unit instance'
    delta = info['delta']
    prof = []
    prev_dim = 0
    d_lf = None
    d_ff = None
    for D in range(1, Dcap + 1):
        t0 = time.time()
        one, dimV, dim_low = closure(m, k, p, g, delta, D)
        prof.append(dict(D=D, dimV=dimV, dimV_low=dim_low, one_in_V=one, sec=round(time.time() - t0, 2)))
        if d_ff is None and dim_low > prev_dim:
            d_ff = D
        prev_dim = dimV
        if one:
            d_lf = D
            break
    return dict(m=m, k=k, p=p, a=a, b=b, t=t, delta=delta, maxwt=info['maxwt'], sd_predicted=info['sd_predicted'],
                d_lf=d_lf, d_ff=d_ff, profile=prof)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--laneA', default='/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneA_prime_m2/data/all_instances.json')
    ap.add_argument('--per-k', type=int, default=5)
    ap.add_argument('--ks', type=str, default='3,4,5')
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--own', type=str, default=None, help='JSONL of my own instances (m,k,p,a,b,t) to run too')
    ap.add_argument('--own-limit', type=int, default=5)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = []
    with open(a.out, 'w') as fo:
        if a.laneA:
            rows = json.load(open(a.laneA))
            for k in [int(x) for x in a.ks.split(',')]:
                cand = [r for r in rows if r['k'] == k and r.get('d_lf') is not None]
                pick = rng.sample(cand, min(a.per_k, len(cand)))
                for r in pick:
                    res = dlf_dff(2, k, r['p'], r['a'], r['b'], r['t'], Dcap=k + 3)
                    res['laneA_d_lf'] = r['d_lf']
                    res['laneA_d_ff'] = r.get('d_ff')
                    res['laneA_d_lf_source'] = r.get('d_lf_source')
                    res['laneA_sd_inverse'] = r['sd_inverse']
                    res['agree_d_lf'] = (res['d_lf'] == r['d_lf'])
                    res['agree_d_ff'] = (r.get('d_ff') is None) or (res['d_ff'] == r.get('d_ff'))
                    fo.write(json.dumps(res) + '\n')
                    fo.flush()
                    print(f"laneA k={k} p={r['p']} a={r['a']} b={r['b']} t={r['t']}: my d_lf={res['d_lf']} d_ff={res['d_ff']} | laneA d_lf={r['d_lf']} d_ff={r.get('d_ff')} "
                          f"sd={res['sd_predicted']} gap={res['sd_predicted'] - res['d_lf']} agree={res['agree_d_lf']} "
                          f"profile={[(x['D'], x['dimV'], x['dimV_low'], x['sec']) for x in res['profile']]}", flush=True)
        if a.own:
            own = [json.loads(l) for l in open(a.own)][:a.own_limit]
            for r in own:
                res = dlf_dff(r['m'], r['k'], r['p'], r['a'], r['b'], r['t'], Dcap=r['m'] * r['k'] + 3)
                res['own_sd_static'] = r.get('sd_static')
                fo.write(json.dumps(res) + '\n')
                fo.flush()
                print(f"own m={r['m']} k={r['k']} p={r['p']}: d_lf={res['d_lf']} d_ff={res['d_ff']} sd_static={r.get('sd_static')} "
                      f"profile={[(x['D'], x['dimV'], x['dimV_low'], x['sec']) for x in res['profile']]}", flush=True)


if __name__ == '__main__':
    main()
