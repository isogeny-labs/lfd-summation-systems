"""Cross-check of the red-team builder + validator engine against the certified sweep values:
every two-summand prime case of experiments/ic-curve-sweep with k <= 4 (subgroup and random coset,
x_i = c*u_{i,0}), exact LFD and |Z| recomputed in the ordinary ring."""
import json, os, collections
from rt_lib import *
root = os.path.join(os.path.dirname(os.path.abspath(E.__file__)), '..', '..', 'ic-curve-sweep')
cases = [json.loads(l) for l in open(os.path.join(root, 'cases.jsonl'))]
res = {json.loads(l)['system_hash']: json.loads(l) for l in open(os.path.join(root, 'results.jsonl'))}
agree = collections.Counter(); bad = []
seen = set()
for c in cases:
    if c['kind'] != 'prime' or c['m'] != 2 or c['k'] > 4 or c['system_hash'] in seen:
        continue
    seen.add(c['system_hash'])
    r = res[c['system_hash']]
    s = c.get('base_scale', 1)
    out, bound, z = exact(c['prime'], c['a'], c['b'], c['target'][0], c['k'], c['k'], s, s)
    ok = (out['lfd'] == r['exact_lfd'] and z == r['point_count'])
    agree[(c['base_kind'], c['k'], ok)] += 1
    if not ok:
        bad.append((c['case_id'], c['prime'], c['a'], c['b'], c['target'], c['k'], s, out['lfd'], r['exact_lfd'], z, r['point_count']))
print('distinct systems compared', len(seen))
for k in sorted(agree): print('  ', k, agree[k])
print('mismatches', len(bad), bad[:10])
