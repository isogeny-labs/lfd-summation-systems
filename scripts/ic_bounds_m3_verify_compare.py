#!/usr/bin/env python3
"""Compare verify/recheck.jsonl (independent engine) with the producer's rows.jsonl -> verify/compare.json."""
import json
from pathlib import Path
from collections import Counter, defaultdict
EXP = Path(__file__).resolve().parents[1] / 'experiments' / 'ic-bounds-m3'
rows = {json.loads(l)['case_id']: json.loads(l) for l in open(EXP / 'rows.jsonl')}
rec = [json.loads(l) for l in open(EXP / 'verify' / 'recheck.jsonl')]
stat = defaultdict(Counter); mism = []
for r in rec:
    P = rows[r['case_id']]; e = r['enc']; key = f"k{r['k']}-{e}"
    if r['status'] == 'capped':
        stat[key]['capped'] += 1; continue
    if e == 'direct':
        theirs = P['direct_lfd']; dimA = P['solutions']
        ranks = {d[0]: d[2] for d in P['direct_dims']}
        status_ok = (r['status'] == 'unit') == (P['direct_status'] == 'exact-unit')
        extra = dict(deg=(r['input_degrees'][0] == P['generator_degree']))
    elif e == 'm2':
        theirs = P['m2_lfd']; dimA = P['m2_solutions']; ranks = {}
        status_ok = (r['status'] == 'unit') == (P['m2_status'] == 'exact-unit')
        extra = {}
    else:
        theirs = P['split_free_lfd']; dimA = P['split_free_dim']
        ranks = {d[0]: d[2] for d in P['split_dims']}
        status_ok = (r['status'] == 'unit') == ('unit' in P['split_free_status'])
        extra = dict(field_dim=(r['dimA_field'] == P['split_field_dim']),
                     field_status=((r['dimA_field'] < r['dimA']) == (P['split_field_status'].startswith('lower-bound-p'))))
    ok_l = r['lfd'] == theirs; ok_d = r['dimA'] == dimA
    ok_r = all(ranks.get(q['D'], q['rank']) == q['rank'] for q in r['profile'])
    checks = dict(lfd=ok_l, dimA=ok_d, ranks=ok_r, status=status_ok, **extra)
    for c, v in checks.items(): stat[key][c + ('_ok' if v else '_BAD')] += 1
    stat[key]['n'] += 1
    stat[key][f"lfd={r['lfd']}"] += 1
    if not all(checks.values()):
        mism.append(dict(case_id=r['case_id'], enc=e, k=r['k'], ours=r['lfd'], theirs=theirs, checks=checks,
                         our_profile=r['profile'], their_ranks=ranks))
out = dict(n_rechecked=len(rec), per_stratum={k: dict(v) for k, v in sorted(stat.items())}, mismatches=mism)
(EXP / 'verify' / 'compare.json').write_text(json.dumps(out, indent=1))
for k, v in sorted(stat.items()): print(k, dict(v))
print('mismatches', len(mism))
for m in mism[:10]: print(json.dumps(m)[:600])
