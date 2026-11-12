"""Validation of the split closure engine (csd/ic_bounds_d5.py) at small depth.

Checks, each recorded in experiments/ic-bounds-d5/validation.json:
  1. instance helpers: prime_points_fast == ic_curve_sweep.prime_points and evaluation_rank_fast ==
     prime_evaluation_rank (k <= 8);
  2. k = 6, 7 (and --dense-k8): the split engine run to a verified closure (no early stop) against the
     dense exact engine ic_fall_certificates.chain_degree_closure at D = 3, 4, 5 -- dimension,
     dim(V_D cap M_<=w) for every w, 1 in V_D -- and against the depth-10 probe's shift engine
     (ic_chain_shift_closure.shift_closure, verified) at D = 4, 5;
  3. k = 8: every depth-8 system of experiments/ic-chain-depth8/systems.jsonl: the split engine at
     D = 4 (verified closure) reproduces the saved dim V_4 ('previous_dimension') and at D = 5
     reaches the saved rank target, with the saved low dimension.
Three instance families at k = 6, 7: the depth-10 probe instance (p = 12289, a = 1, b = 0, seed 10,
unit ideal), the second depth-10 instance (p = 13313, seed 11, proper ideal) and the third
(p = 15361, a = 4021, b = 9917, seed 12, proper ideal).
"""
import argparse, hashlib, json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from csd import ic_bounds_d5 as B
from csd import ic_fall_certificates as C
from csd import ic_chain_shift_closure as SC
from csd import ic_curve_sweep as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--depths', default='6,7')
ap.add_argument('--dense-k8', action='store_true', help='also run the dense engine at k = 8 (D <= 5)')
ap.add_argument('--depth8-limit', type=int, default=None, help='only the first N depth-8 systems')
ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-d5/validation.json')
args = ap.parse_args()
INSTANCES = [dict(name='probe-unit', p=12289, a=1, b=0, seed=10, want=False),
             dict(name='proper-13313', p=13313, a=1, b=0, seed=11, want=True),
             dict(name='random-15361', p=15361, a=4021, b=9917, seed=12, want=True)]
report = dict(helpers=[], small_depth=[], depth8=[])
t0 = time.time()


def dense_profile(V, D):
    w = [V.masks[j].bit_count() for j in V.pivots]
    return [sum(1 for x in w if x <= d) for d in range(D + 1)]


# 1. helpers
for spec in INSTANCES:
    for k in (6, 7, 8):
        inst = B.make_instance(spec['p'], spec['a'], spec['b'], spec['seed'], k, spec['want'])
        E = inst['E']
        slow = S.prime_points(E, k, inst['g0'], 1, 1, inst['target'][0])
        masks = [m for m in B.monomials(k, 5) if m.bit_count() < 5]
        er_fast = B.evaluation_rank_fast(inst['points'], masks, k, inst['g0'], E.p)
        er_slow = S.prime_evaluation_rank(slow, masks, k, inst['g0'], E.p, 10 ** 9) if slow else 0
        row = dict(instance=spec['name'], k=k, points=len(slow), points_equal=slow == inst['points'],
                   eval_rank_equal=er_fast == er_slow, eval_rank=er_fast)
        assert row['points_equal'] and row['eval_rank_equal'], row
        report['helpers'].append(row)
print(json.dumps(dict(stage='helpers', ok=True, rows=len(report['helpers']))), flush=True)

# 2. small depth: dense vs split vs shift
depths = [int(x) for x in args.depths.split(',')] + ([8] if args.dense_k8 else [])
for k in depths:
    for spec in INSTANCES:
        inst = B.make_instance(spec['p'], spec['a'], spec['b'], spec['seed'], k, spec['want'])
        g, Z, p = inst['generator'], inst['points'], inst['E'].p
        for D in (3, 4, 5):
            t = time.time()
            V = C.chain_degree_closure([g], 2, k, p, D, max_columns=20000, max_cells=1_200_000_000)
            td = time.time() - t
            dense = dict(dim=V.quotient_rank, contains_one=V.contains_one, profile=dense_profile(V, D))
            if not V.contains_one:
                dense['exclusion_certificate'] = C.check_closed_exclusion_certificate(V, [g])
            t = time.time()
            r = B.split_closure(g, k, p, D, B.column_count(k, D) - len(Z), early_stop=False)
            ts = time.time() - t
            split = dict(dim=r['rank'], contains_one=r['contains_one'], profile=r['profile'], status=r['status'],
                         sweeps=r['sweeps'], echelon_ok=r['echelon_ok'], ranks=r['ranks'])
            row = dict(instance=spec['name'], k=k, D=D, points=len(Z), dense=dense, split=split,
                       seconds=dict(dense=round(td, 2), split=round(ts, 2)))
            if D >= 4:
                t = time.time()
                sh = SC.shift_closure(g, k, p, D, 10 ** 9, product_cap=6000, max_cells=1_200_000_000, verify=True)
                row['shift'] = dict(dim=sh['rank'], closed=sh['closed'], low_dimension=sh['low_dimension'])
                row['seconds']['shift'] = round(time.time() - t, 2)
            ok = (split['status'] == 'closed' and split['echelon_ok'] and split['dim'] == dense['dim']
                  and split['profile'] == dense['profile'] and split['contains_one'] == dense['contains_one']
                  and dense.get('exclusion_certificate', True)
                  and ('shift' not in row or (row['shift']['dim'] == dense['dim'] and row['shift']['closed'])))
            row['agree'] = ok
            report['small_depth'].append(row)
            print(json.dumps(dict(stage='small', k=k, D=D, instance=spec['name'], dim=dense['dim'], agree=ok,
                                  profile=split['profile'], seconds=row['seconds'])), flush=True)
            assert ok, row

# 3. depth 8 against the saved results
systems = [json.loads(l) for l in (ROOT / 'experiments/ic-chain-depth8/systems.jsonl').read_text().splitlines()]
results = {r['system_hash']: r for r in map(json.loads, (ROOT / 'experiments/ic-chain-depth8/results.jsonl').read_text().splitlines())}
systems = [s for s in systems if s['k'] == 8][:args.depth8_limit]
for sysm in systems:
    g = {int(m): int(c) for m, c in sysm['generators'][0]}
    k, p, Z = sysm['k'], sysm['prime'], sysm['coordinate_solutions']
    saved = results[sysm['system_hash']]
    chk = saved['certificate_checks'][0]
    t = time.time()
    r4 = B.split_closure(g, k, p, 4, B.column_count(k, 4) - len(Z), early_stop=False)
    r5 = B.split_closure(g, k, p, 5, B.column_count(k, 5) - len(Z), early_stop=True)
    row = dict(system_hash=sysm['system_hash'], prime=p, points=len(Z), saved_status=saved['status'],
               saved_exact_lfd=saved['exact_lfd'], saved_dim_V4=chk['previous_dimension'], split_dim_V4=r4['rank'],
               split_V4_status=r4['status'], saved_target=chk['target_rank'], saved_rank=chk['rank'],
               split_rank_5=r5['rank'], split_status_5=r5['status'], saved_low=chk['low_dimension'],
               split_low_5=r5['low_dimension'], seconds=round(time.time() - t, 1))
    row['agree'] = (r4['status'] == 'closed' and r4['rank'] == chk['previous_dimension']
                    and r5['status'] == 'reached' and r5['rank'] == chk['rank'] == chk['target_rank']
                    and r5['low_dimension'] == chk['low_dimension'])
    report['depth8'].append(row)
    print(json.dumps(dict(stage='depth8', **row)), flush=True)
    assert row['agree'], row

report['summary'] = dict(helpers=len(report['helpers']), small_depth=len(report['small_depth']),
                         depth8=len(report['depth8']),
                         all_agree=all(r['agree'] for r in report['small_depth'] + report['depth8']),
                         seconds=round(time.time() - t0, 1))
report['source_sha256'] = {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest()
                           for n in ['scripts/ic_bounds_d5_validate.py', 'csd/ic_bounds_d5.py',
                                     'csd/ic_chain_shift_closure.py', 'csd/ic_fall_certificates.py',
                                     'csd/ic_curve_sweep.py', 'csd/ic_dreg.py', 'csd/elim.py', 'csd/_elim.c']}
report['data_sha256'] = {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest()
                         for n in ['experiments/ic-chain-depth8/systems.jsonl', 'experiments/ic-chain-depth8/results.jsonl']}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(report, indent=1) + '\n')
print(json.dumps(report['summary']))
