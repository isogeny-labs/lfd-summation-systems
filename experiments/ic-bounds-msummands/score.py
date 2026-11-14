#!/usr/bin/env python3
"""Scores predictions-2026-11-12.json against records/*.json -> scores.json, results.tsv.  Read-only on both inputs."""
import json, os, glob, hashlib, collections
HERE = os.path.dirname(os.path.abspath(__file__))
pred = json.load(open(os.path.join(HERE, 'predictions-2026-11-12.json')))
recs = {}
for f in glob.glob(os.path.join(HERE, 'records', '*.json')):
    r = json.load(open(f)); recs[r['case']] = r
rows = []; sc = collections.Counter()
for name, p in sorted(pred['cases'].items()):
    r = recs.get(name)
    row = dict(case=name, family=p['family'], m=p['m'], k=p['k'])
    if r is None:
        row.update(outcome='no-record'); sc['no-record'] += 1; rows.append(row); continue
    ex = r.get('exact_lfd')
    row.update(status=r['status'], lfd=ex, solutions=r.get('solutions'), ideal_class=r.get('ideal_class'),
               certificate_degree=r.get('certificate_degree'), strict_falls=r.get('strict_falls'), seconds=r.get('seconds'))
    if p['family'] == 'prime':
        row.update(generator_degree=p['generator_degree'], predicted=p['predicted_lfd'], kernel_bound=p['kernel_upper_bound'])
        if ex is None:
            row['hit'] = False; sc['prime-nonexact'] += 1
        else:
            row['hit'] = ex == p['predicted_lfd']
            row['falls_hit'] = r['strict_falls'] == p['predicted_strict_falls']
            row['within_kernel'] = ex <= p['kernel_upper_bound']
            sc['prime-exact'] += 1; sc['prime-hit'] += row['hit']; sc['prime-falls-hit'] += row['falls_hit']; sc['prime-within-kernel'] += row['within_kernel']
            if p.get('lfd_lower_bound_if_prop1') is not None:
                sc['prime-prop1-lower-ok'] += ex >= p['lfd_lower_bound_if_prop1']; sc['prime-prop1-cases'] += 1
    else:
        row.update(variables=p['variables'], input_degree=p['input_degree'], pred_upper=p['pred_c_upper'], pred_point=p['pred_c_point'], scored_point=p['point_scored'])
        if ex is None:
            row['hit_upper'] = None; row['hit_point'] = None
            sc['binary-nonexact-or-not-run'] += 1
        else:
            row['hit_upper'] = ex <= p['pred_c_upper']; row['hit_point'] = ex == p['pred_c_point']
            row['below_kernel'] = ex < p['kernel_upper_bound']
            sc['binary-exact'] += 1; sc['binary-upper-hit'] += row['hit_upper']; sc['binary-below-r+1'] += row['below_kernel']
            if p['point_scored']:
                sc['binary-point-cases'] += 1; sc['binary-point-hit'] += row['hit_point']
    rows.append(row)
json.dump(dict(counts=dict(sc), rows=rows), open(os.path.join(HERE, 'scores.json'), 'w'), indent=1)
print(json.dumps(dict(sc), indent=1))
