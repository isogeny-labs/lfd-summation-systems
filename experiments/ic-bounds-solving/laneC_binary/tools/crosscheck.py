#!/usr/bin/env python3
"""Cross-check tools/static_scan.py against the repository's own static_degrees (scripts/ic_bounds_binary_related.py)
on a few small systems: d_full must agree exactly (same definition, different elimination kernel and independently
brute-forced Z).  Also checks the new rank kernel against csd.elim.gf2_echelon on the actual Macaulay matrices."""
import sys, json, os
sys.path.insert(0, '/home/user/ml-cryptanalysis'); sys.path.insert(0, '/home/user/ml-cryptanalysis/scripts')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import static_scan as S
from csd import elim
from csd import ic_dreg as I
from ic_bounds_binary_related import static_degrees
LANE = os.path.dirname(HERE)
inst = {x['system_hash']: x for x in S.load_instances({5, 6, 7})}
done = [json.loads(l) for l in open(os.path.join(LANE, 'data/static_scan_k2_7.jsonl'))]
res = []
import random
random.Random(1).shuffle(done)
cnt = 0
for d in done:
    if d['k'] not in (5, 6, 7) or d['Z'] == 0 and d['k'] == 7 or cnt >= int(sys.argv[1]):
        continue
    x = inst[d['system_hash']]
    eqs, r = x['eqs'], 2 * x['k']
    pts = S.brute_points(eqs, r)
    ref = static_degrees(dict(equations=eqs, variables=r, points=pts), 70000, 12)
    # kernel check at D = 4
    M, cols = I.boolean_macaulay(eqs, r, 4)
    r1, p1, _ = elim.gf2_echelon(M, cols)
    r2, p2 = S.echelon(M.copy(), cols, 1)
    ok = (ref['d_full'] == d['d_full'] and ref['d_lin'] == None or True)
    rec = dict(k=d['k'], n=d['n'], hash=d['system_hash'][:12], ref_d_full=ref['d_full'], scan_d_full=d['d_full'], ref_d_lin=ref['d_lin'],
               kernel_rank_equal=(r1 == r2 and list(p1) == list(p2)), match=(ref['d_full'] == d['d_full']))
    res.append(rec); cnt += 1; print(rec, flush=True)
json.dump(res, open(os.path.join(LANE, 'data/crosscheck.json'), 'w'), indent=1)
print('all match:', all(r['match'] and r['kernel_rank_equal'] for r in res), len(res))
