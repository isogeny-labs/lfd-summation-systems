#!/usr/bin/env python3
"""P5 first half and d_full vs LFD counts over ALL rows of data/binary_table.jsonl (all families, rows with exact d_full),
plus lower-bound statements; writes data/p5_check.json and prints a table."""
import json, os, collections
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = [json.loads(l) for l in open(os.path.join(LANE, 'data/binary_table.jsonl'))]
out = {}
tot = collections.Counter()
for fam in ['m2', 'm3-direct', 'm3-split']:
    rows = [r for r in R if r['family'] == fam]
    ex = [r for r in rows if r['d_full'] is not None]
    lb = [r for r in rows if r['d_full'] is None and r['d_full_lb'] is not None]
    c = collections.defaultdict(collections.Counter)
    for r in ex:
        c[r['k']]['gt' if r['d_full'] > r['lfd'] else ('eq' if r['d_full'] == r['lfd'] else 'LT')] += 1
    for r in lb:
        c[r['k']]['lb>LFD' if r['d_full_lb'] > r['lfd'] else 'lb<=LFD(open)'] += 1
    viol_dstar = [r['case_id'] for r in ex if r['dstar'] is not None and r['dstar'] > r['d_full']]
    viol_lt = [r['case_id'] for r in ex if r['d_full'] < r['lfd']]
    viol_lb = [r['case_id'] for r in lb if r['dstar'] is not None and r['d_full_lb'] > 0 and False]
    out[fam] = dict(rows=len(rows), exact_dfull_rows=len(ex), lower_bound_rows=len(lb), by_k={k: dict(v) for k, v in sorted(c.items())},
                    dstar_gt_dfull=viol_dstar, dfull_lt_lfd=viol_lt)
    print(fam, out[fam]['rows'], 'exact', len(ex), 'lb', len(lb), 'D*>d_full:', len(viol_dstar), 'd_full<LFD:', len(viol_lt))
    for k, v in sorted(c.items()):
        print('   k=%d' % k, dict(v))
json.dump(out, open(os.path.join(LANE, 'data/p5_check.json'), 'w'), indent=1, sort_keys=True)
