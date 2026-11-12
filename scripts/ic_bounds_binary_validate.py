#!/usr/bin/env python3
"""Independent re-derivation of a sample of item-6 rows.

For each sampled row the system is rebuilt from its specification (field modulus,
curve, basis, target), its hash is compared with the stored one, the zero set is
re-enumerated by brute force over all 2^r Boolean assignments (r <= --max-brute),
and the exact LFD is recomputed with the independent validator engine
experiments/ic-curve-validation/code/lfd_engine.py, which works in the ORDINARY
ring GF(2)[z_1..z_r] with the field equations z_i^2 + z_i as inputs (no Boolean
normal form, no shared closure code). The engine's certificate is the same
local-to-global rank criterion (dim V_D = dim I_D, evaluation of degree < D onto Z).

  python3 scripts/ic_bounds_binary_validate.py --out experiments/ic-bounds-binary/validation.json
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'ic-curve-validation', 'code'))

import lfd_engine as V  # noqa: E402
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_binary as B  # noqa: E402
from csd.ic_research_campaign import field  # noqa: E402


def rebuild(row):
    F = field(row['n'])
    assert F.modulus == row['modulus']
    E = I.BinaryCurve(F, row['a'], row['b'])
    target = tuple(row['target'])
    planted = [tuple(P) for P in row['planted']] if row['planted'] else None
    if row['family'] == 'm2':
        return B.m2_system(E, row['basis'], row['base_kind'], target, row['target_class'], planted)
    if row['family'] == 'm3-direct':
        return B.m3_direct_system(E, row['basis'], target)
    return B.m3_split_system(E, row['basis'], target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', nargs='+', default=sorted(glob.glob(os.path.join(ROOT, 'experiments/ic-bounds-binary/results-*.jsonl'))))
    ap.add_argument('--per-family', type=int, default=8)
    ap.add_argument('--max-vars', type=int, default=12)
    ap.add_argument('--max-brute', type=int, default=16)
    ap.add_argument('--max-degree', type=int, default=8)
    ap.add_argument('--max-ordinary-columns', type=int, default=6000,
                    help='the validator engine is dense float32: keep C(r+D, r) small')
    ap.add_argument('--seed', type=int, default=20261108)
    ap.add_argument('--hash-all', action='store_true', help='also rebuild and hash-check every row (no LFD)')
    ap.add_argument('--mxl-per-k', type=int, default=0,
                    help='engine-B check: rows per (family, k) whose closure ranks are recomputed with the '
                         'MutantXL closure of csd.ic_dreg (_mxl_closure on the full Macaulay matrix)')
    ap.add_argument('--mxl-max-columns', type=int, default=30000)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    rows = [json.loads(l) for f in args.results for l in open(f) if l.strip()]
    rng = random.Random(args.seed)
    sample = []
    for fam in sorted({r['family'] for r in rows}):
        pool = [r for r in rows if r['family'] == fam and r['exact_lfd'] is not None
                and r['variables'] <= args.max_vars and r['exact_lfd'] <= args.max_degree
                and math.comb(r['variables'] + max(r['exact_lfd'], r['certificate_degree']) + 1, r['variables'])
                <= args.max_ordinary_columns]
        # prefer distinct LFD values, then random
        pool.sort(key=lambda r: r['case_id'])
        rng.shuffle(pool)
        by = {}
        for r in pool:
            by.setdefault(r['exact_lfd'], []).append(r)
        chosen = []
        while len(chosen) < args.per_family and any(by.values()):
            for v in sorted(by):
                if by[v] and len(chosen) < args.per_family:
                    chosen.append(by[v].pop())
        sample += chosen
    checks = []
    for row in sample:
        t0 = time.time()
        system = rebuild(row)
        assert system['system_hash'] == row['system_hash'], row['case_id']
        r = system['variables']
        rec = dict(case_id=row['case_id'], family=row['family'], n=row['n'], k=row['k'], variables=r,
                   system_hash=row['system_hash'], stored_lfd=row['exact_lfd'], stored_solutions=row['solutions'])
        if r <= args.max_brute:
            brute = V.binary_points(system['equations'], r)
            as_masks = sorted(sum(b << i for i, b in enumerate(p)) for p in brute)
            rec['brute_force_solutions'] = len(brute)
            assert as_masks == sorted(system['points']), row['case_id']
        else:
            brute = [tuple((p >> i) & 1 for i in range(r)) for p in system['points']]
        inputs = V.binary_inputs(system['equations'], r)
        res = V.exact_lfd(inputs, r, 2, brute, max(row['exact_lfd'], row.get('certificate_degree') or 0) + 1)
        rec.update(validator_lfd=res['lfd'], validator_largest_fall=res.get('lfd_largest_fall'),
                   validator_certified_at=res['certified_at'], seconds=round(time.time() - t0, 2))
        # ordinary-ring and Boolean presentations agree for LFD >= 2 (field equations of degree 2)
        rec['agree'] = res['lfd'] == row['exact_lfd'] or (row['exact_lfd'] < 2 and res['lfd'] is not None and res['lfd'] <= 2)
        print(rec, flush=True)
        checks.append(rec)
    mxl = []
    if args.mxl_per_k:
        from csd.ic_dreg import _mxl_closure, _bool_mulvar, boolean_macaulay, n_cols
        groups = {}
        for row in rows:
            if row['exact_lfd'] is not None:
                groups.setdefault((row['family'], row['k']), []).append(row)
        for key in sorted(groups):
            pool = sorted(groups[key], key=lambda r: r['case_id'])
            rng.shuffle(pool)
            # prefer proper ideals (|Z| > 0) with the largest LFD in the group
            pool.sort(key=lambda r: (-(r['solutions'] > 0), -r['exact_lfd']))
            for row in pool[:args.mxl_per_k]:
                system = rebuild(row)
                assert system['system_hash'] == row['system_hash']
                r = system['variables']
                t0 = time.time()
                per = []
                for D in range(0, row['certificate_degree'] + 1):
                    if n_cols(r, D) > args.mxl_max_columns:
                        per.append(dict(D=D, skipped='column cap'))
                        continue
                    if not any(max(m.bit_count() for m in e) <= D for e in system['equations']):
                        per.append(dict(D=D, mxl_rank=0, stored_rank=row['quotient_ranks'][D],
                                        agree=row['quotient_ranks'][D] == 0))
                        continue
                    M, cols = boolean_macaulay(system['equations'], r, D)
                    rank, piv, _, rounds, _ = _mxl_closure(M, cols, r, D, 2, _bool_mulvar, max_rounds=10 ** 6)
                    stored = row['quotient_ranks'][D]
                    if row['status'] == 'exact-unit' and D == row['certificate_degree']:
                        # MutantXL does not multiply a constant row: compare "1 in V_D" instead of the rank
                        per.append(dict(D=D, mxl_rank=rank, stored_rank=stored, rounds=rounds,
                                        compared='contains-one', agree=(cols - 1) in set(piv)))
                        continue
                    per.append(dict(D=D, mxl_rank=rank, stored_rank=stored, rounds=rounds, agree=rank == stored))
                rec = dict(case_id=row['case_id'], family=row['family'], n=row['n'], k=row['k'], variables=r,
                           solutions=row['solutions'], stored_lfd=row['exact_lfd'], per_degree=per,
                           agree=all(p.get('agree', True) for p in per),
                           checked_degrees=sum('agree' in p for p in per), seconds=round(time.time() - t0, 2))
                print(rec, flush=True)
                mxl.append(rec)
    hashed = 0
    if args.hash_all:
        for row in rows:
            assert rebuild(row)['system_hash'] == row['system_hash'], row['case_id']
            hashed += 1
    out = dict(engine='experiments/ic-curve-validation/code/lfd_engine.py (ordinary ring, field equations as inputs)',
               sampled=len(checks), agree=sum(c['agree'] for c in checks), rows_rehashed=hashed, checks=checks,
               engine_b='csd.ic_dreg._mxl_closure (MutantXL closure of the full Boolean Macaulay matrix, run to '
                        'a fixed point); compared rank by rank with the stored local_profile quotient ranks',
               mxl_sampled=len(mxl), mxl_agree=sum(c['agree'] for c in mxl), mxl_checks=mxl)
    json.dump(out, open(args.out, 'w'), indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in out.items() if k not in ('checks', 'mxl_checks')}))


if __name__ == '__main__':
    main()
