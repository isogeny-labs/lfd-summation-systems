"""Is the exact LFD nondecreasing in the depth for a fixed (p, a, b, t)?  (U1, Section 5.2, T4.)

Reads the certified exact values of the subgroup two-summand prime cases of
experiments/ic-curve-sweep, ic-chain-depth8 and ic-chain-depth9 and compares every pair of
depths k < k' of the same instance.  For a unit ideal at the deeper depth the inequality
LFD(k) <= LFD(k') is Corollary 3.1; for proper ideals it is only observed.
"""
import collections
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for d in ['ic-curve-sweep', 'ic-chain-depth8', 'ic-chain-depth9']:
    root = os.path.join(HERE, '..', d)
    res = {}
    for l in open(os.path.join(root, 'results.jsonl')):
        r = json.loads(l)
        res[r['system_hash']] = r
    for l in open(os.path.join(root, 'cases.jsonl')):
        c = json.loads(l)
        if c['kind'] == 'prime' and c['m'] == 2 and c['base_kind'] == 'subgroup' and c['system_hash'] in res and c['target']:
            r = res[c['system_hash']]
            rows.append((c['prime'], c['a'], c['b'], c['target'][0], c['k'], r['exact_lfd'], r['point_count']))
g = collections.defaultdict(dict)
for p, a, b, t, k, l, z in rows:
    g[(p, a, b, t)][k] = (l, z)
multi = {key: v for key, v in g.items() if len(v) > 1}
viol = pairs = unitpairs = properpairs = 0
ex = []
for key, v in multi.items():
    ks = sorted(v)
    for i in range(len(ks) - 1):
        for j in range(i + 1, len(ks)):
            pairs += 1
            if v[ks[j]][1] == 0:
                unitpairs += 1
            else:
                properpairs += 1
            if v[ks[i]][0] > v[ks[j]][0]:
                viol += 1
                ex.append((key, ks[i], v[ks[i]], ks[j], v[ks[j]]))
out = dict(instances=len(g), instances_with_several_depths=len(multi), depth_pairs=pairs,
           pairs_deeper_unit=unitpairs, pairs_deeper_proper=properpairs, violations=viol,
           examples=[str(e) for e in ex[:10]])
json.dump(out, open(os.path.join(HERE, 'results-sweep-monotone.json'), 'w'), indent=1)
print(out)
