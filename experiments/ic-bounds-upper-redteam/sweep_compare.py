import json, sys, collections
import os
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..') + os.sep
def load(path):
    return [json.loads(l) for l in open(path) if l.strip()]
for name in ['ic-curve-sweep', 'ic-chain-depth8', 'ic-chain-depth9']:
    cases = load(root + name + '/cases.jsonl')
    results = {r['system_hash']: r for r in load(root + name + '/results.jsonl')}
    print('==', name, 'cases', len(cases), 'results', len(results))
    kinds = collections.Counter((c.get('kind'), c.get('m'), c.get('base_kind')) for c in cases)
    print(' kinds', dict(kinds))
    statuses = collections.Counter(r.get('status') for r in results.values())
    print(' statuses', dict(statuses))
    hist = collections.defaultdict(collections.Counter)
    viol = []
    nores = 0
    seen = set()
    for c in cases:
        if c.get('kind') != 'prime' or c.get('m') != 2:
            continue
        r = results.get(c['system_hash'])
        if r is None:
            nores += 1
            continue
        k = c['k']
        key = (c['base_kind'], k)
        ex = r.get('exact_lfd')
        lo, up = r.get('lfd_lower_bound'), r.get('lfd_upper_bound')
        hist[key][(r.get('status'), ex, lo, up)] += 1
        # violation tests: exact value above k+2, or a certified lower bound above k+2
        if (ex is not None and ex > k + 2) or (lo is not None and lo > k + 2):
            viol.append((c['case_id'], c['prime'], c['a'], c['b'], c['target'], k, c['base_kind'], c.get('base_scale'), ex, lo, up))
    print(' prime m=2 cases without result:', nores)
    for key in sorted(hist):
        print(' ', key, 'k+2 =', key[1] + 2, dict(hist[key]))
    print(' VIOLATIONS of LFD <= k+2:', viol)
