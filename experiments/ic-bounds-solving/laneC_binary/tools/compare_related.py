#!/usr/bin/env python3
"""Agreement of tools/static_scan.py with the repository's related-m2.jsonl on every system present in both:
exact vs exact must be equal; repository cap 'capped_at = c' (d_full >= c) must be consistent with ours (exact >= c or lower bound >= c)."""
import json, glob, os, collections
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rel = {}
for l in open('/home/user/ml-cryptanalysis/experiments/ic-bounds-binary/related-m2.jsonl'):
    x = json.loads(l); rel.setdefault(x['system_hash'], x)
res = collections.Counter(); bad = []; rows = []
for f in sorted(glob.glob(os.path.join(LANE, 'data/static_scan_*.jsonl'))):
    for l in open(f):
        x = json.loads(l)
        q = rel.get(x['system_hash'])
        if not q: continue
        mine_lo = x['d_full'] if x['d_full'] is not None else x['d_full_lb']
        if q['d_full'] is not None:
            ok = (x['d_full'] == q['d_full']); kind = 'exact-vs-exact'
        else:
            ok = (mine_lo >= q['capped_at']); kind = 'cap-vs-ours'
        res[(kind, ok)] += 1
        rows.append(dict(k=x['k'], n=x['n'], hash=x['system_hash'][:12], repo_d_full=q['d_full'], repo_capped_at=q['capped_at'],
                         ours_d_full=x['d_full'], ours_lb=x['d_full_lb'], ok=ok))
        if not ok: bad.append(rows[-1])
print(dict(res)); print('disagreements:', bad)
json.dump(dict(summary={str(k): v for k, v in res.items()}, rows=rows), open(os.path.join(LANE, 'data/related_agreement.json'), 'w'), indent=1)
