#!/usr/bin/env python3
"""Verification pass: static (non-iterated) Boolean Macaulay degree d_full for sampled rows.

Rows of the degree-D matrix: mu * f (Boolean normal form) for every input f and every
squarefree monomial mu with deg mu <= D - deg f.  d_full = least D at which the row
space has codimension |Z| in M_<=D and M_<D evaluates onto Z.  Uses the system builder
and GF(2) rank kernel of scripts/ic_bounds_binary_verify.py (no producer code).

  python3 scripts/ic_bounds_binary_verify_static.py --family m3-direct --case-ids ... --max-columns 20000 --out F
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ic_bounds_binary_verify as Vf  # noqa: E402


def static_rank(r, eqs, D, with_one=False):
    S = Vf.Space(r, min(D, r))
    polys = []
    for f in eqs:
        df = max(bin(m).count('1') for m in f)
        if df > D:
            continue
        for mu in S.cols[S.deg <= min(D - df, r)].tolist():
            acc = {}
            for m in f:
                t = m | mu
                acc[t] = acc.get(t, 0) ^ 1
            polys.append([t for t, c in acc.items() if c])
    if with_one:
        polys.append([0])
    if not polys:
        return S.C, 0, 0
    P = S.pack(polys)
    work = np.zeros((S.C + 1, S.W), dtype=np.uint64)
    rk = Vf.lib().rank(S.C, len(polys), np.ascontiguousarray(P), work, S.C + 1)
    return S.C, len(polys), rk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--case-ids', nargs='+', required=True)
    ap.add_argument('--max-columns', type=int, default=20000)
    ap.add_argument('--max-degree', type=int, default=12)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rows = {r['case_id']: r for r in Vf.load(a.family)}
    rel = {}
    p = os.path.join(Vf.EXP, f'related-{a.family}.jsonl')
    for l in open(p):
        x = json.loads(l)
        rel[x['case_id']] = x
    for cid in a.case_ids:
        row = rows[cid]
        t0 = time.time()
        r, eqs, Z, h = Vf.build(row)
        per = []
        dfull = None
        dmax_in = max(bin(m).count('1') for e in eqs for m in e)
        for D in range(1, min(r + dmax_in, a.max_degree) + 1):
            C = sum(math.comb(r, i) for i in range(min(D, r) + 1))
            if C > a.max_columns:
                break
            C, nrows, rk = static_rank(r, eqs, D)
            per.append(dict(D=D, columns=C, rows=nrows, rank=rk))
            if not Z:
                # unit ideal: the certificate is 1 in the row space
                if static_rank(r, eqs, D, with_one=True)[2] == rk and nrows:
                    dfull = D
                    break
                continue
            er = Vf.eval_rank(r, Z, min(D, r + 1))
            if rk == C - len(Z) and er == len(Z):
                dfull = D
                break
        stored = rel.get(cid, {})
        sp = {x['D']: x for x in stored.get('per_degree', [])}
        rank_match = all(sp[x['D']]['rank'] == x['rank'] for x in per if x['D'] in sp)
        res = dict(case_id=cid, family=a.family, n=row['n'], k=row['k'], r=r, hash_match=h == row['system_hash'],
                   solutions=len(Z), dstar=row['certificate_degree'], d_full=dfull,
                   stored_d_full=stored.get('d_full'), stored_capped_at=stored.get('capped_at'),
                   rank_match_where_stored=rank_match, per_degree=per, seconds=round(time.time() - t0, 2))
        with open(a.out, 'a') as f:
            f.write(json.dumps(res, sort_keys=True) + '\n')
        print(cid, a.family, 'k', row['k'], 'D*', row['certificate_degree'], 'd_full', dfull, '/', stored.get('d_full'),
              'ranks', rank_match, res['seconds'], 's', flush=True)


if __name__ == '__main__':
    main()
