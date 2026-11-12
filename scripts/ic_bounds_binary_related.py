#!/usr/bin/env python3
"""Static Macaulay (XL) degrees for every item-6 row: how far plain linear algebra on the
Macaulay matrix is from the iterated closure certificate.

For each row (system rebuilt from its specification and hash-checked) and each degree D
from the input degree upward, the bit-packed Boolean Macaulay matrix Mac_D (every
equation times every squarefree monomial of degree <= D - deg, field equations built in;
csd.ic_dreg.boolean_macaulay) is echelonised over GF(2) and two degrees are recorded:

* d_lin  : least D at which the row space of Mac_D contains every affine-linear
           polynomial vanishing on Z (csd.ic_dreg's d_solve convention);
* d_full : least D at which rank Mac_D = |M_<=D| - |Z| and the monomials of degree < D
           evaluate onto Z, i.e. the static matrix alone already meets the local-to-global
           certificate (row space = I cap M_<=D); for a unit ideal, least D with 1 in the
           row space.

Both are exact (GF(2) ranks). The scan stops at the column cap; a capped degree is
reported as a lower bound "> D" only. The closure certificate degree D* <= d_full always
(row space of Mac_D is contained in V_D); the difference is the cost of NOT iterating.

  python3 scripts/ic_bounds_binary_related.py --results experiments/ic-bounds-binary/results-m2.jsonl \
      --out experiments/ic-bounds-binary/related-m2.jsonl
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from csd import elim  # noqa: E402
from csd import ic_dreg as I  # noqa: E402
from csd.ic_local_certificate import evaluation_rank  # noqa: E402
from ic_bounds_binary_validate import rebuild  # noqa: E402


def static_degrees(system, max_columns, max_degree):
    eqs, r, pts = system['equations'], system['variables'], system['points']
    Z = len(pts)
    lin_target = I._lin_target_f2(pts, r)
    d_in = max(m.bit_count() for e in eqs for m in e)
    out = dict(d_lin=None, d_full=None, per_degree=[], capped_at=None)
    for D in range(d_in, max_degree + 1):
        cols = I.n_cols(r, D)
        if cols > max_columns:
            out['capped_at'] = D
            break
        M, cols = I.boolean_macaulay(eqs, r, D)
        rows = M.shape[0]
        rank, piv, _ = elim.gf2_echelon(M, cols)
        del M
        lin = sum(1 for c in piv if c >= cols - (r + 1))
        assert lin <= lin_target
        one = (cols - 1) in set(piv)
        if Z:
            masks = I._masks_upto(r, D)
            ev = evaluation_rank(pts, [int(m) for m in masks if int(m).bit_count() < D]) if D else 0
            full = rank == cols - Z and ev == Z
        else:
            full = one
        out['per_degree'].append(dict(D=D, rows=rows, columns=cols, rank=rank, linear_forms=lin, one=one))
        if out['d_lin'] is None and lin == lin_target:
            out['d_lin'] = D
        if full:
            out['d_full'] = D
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', nargs='+', required=True)
    ap.add_argument('--max-columns', type=int, default=32000)
    ap.add_argument('--max-degree', type=int, default=12)
    ap.add_argument('--per-group', type=int, default=0, help='0 = every row')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)['case_id'] for l in open(args.out) if l.strip()}
    rows = [json.loads(l) for f in args.results for l in open(f) if l.strip()]
    if args.per_group:
        # deterministic subsample: the first per_group rows (by case_id) of every
        # (family, n, k, curve_class, base_kind, target_class) group
        groups = {}
        for r in sorted(rows, key=lambda r: r['case_id']):
            groups.setdefault((r['family'], r['n'], r['k'], r['curve_class'], r['base_kind'], r['target_class']), []).append(r)
        rows = [r for g in groups.values() for r in g[:args.per_group]]
    rows = [r for r in rows if r['case_id'] not in done]
    print(f'{len(rows)} rows to process', flush=True)
    t0 = time.time()
    with open(args.out, 'a') as fh:
        for i, row in enumerate(rows):
            system = rebuild(row)
            assert system['system_hash'] == row['system_hash']
            t1 = time.time()
            res = static_degrees(system, args.max_columns, args.max_degree)
            rec = dict(case_id=row['case_id'], family=row['family'], n=row['n'], k=row['k'],
                       system_hash=row['system_hash'], solutions=row['solutions'], exact_lfd=row['exact_lfd'],
                       certificate_degree=row.get('certificate_degree'), seconds=round(time.time() - t1, 2), **res)
            if rec['d_full'] is not None and rec['certificate_degree'] is not None:
                assert rec['certificate_degree'] <= rec['d_full'], rec
            fh.write(json.dumps(rec, sort_keys=True) + '\n')
            fh.flush()
            print(f"[{i + 1}/{len(rows)} {time.time() - t0:.0f}s] {row['family']} n={row['n']} k={row['k']} "
                  f"lfd={row['exact_lfd']} D*={row.get('certificate_degree')} d_lin={res['d_lin']} "
                  f"d_full={res['d_full']} capped_at={res['capped_at']}", flush=True)


if __name__ == '__main__':
    main()
