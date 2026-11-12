#!/usr/bin/env python3
"""Score the preregistered predictions (predictions.json) against results-prereg.jsonl.

Per distinct affine system (system_hash): hit iff the certified exact LFD equals the
law's prediction.  Also checks that each result row's system_hash and |Z| equal the
preregistered ones (the run rebuilt the same systems).  Output: score.json.
"""
import json, os, collections

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    pre = json.load(open(os.path.join(HERE, 'predictions.json')))
    res = {}
    for l in open(os.path.join(HERE, 'results-prereg.jsonl')):
        x = json.loads(l); res[x['case_id']] = x
    mismatch, missing = [], []
    seen = {}
    for c in pre['cases']:
        r = res.get(c['case_id'])
        if r is None:
            missing.append(c['case_id']); continue
        if r['system_hash'] != c['system_hash'] or r['solutions'] != c['Z']:
            mismatch.append(c['case_id'])
        if r['exact_lfd'] is None:
            missing.append(c['case_id']); continue
        seen.setdefault(c['system_hash'], (c, r))
    cells = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    detail = []
    for h, (c, r) in seen.items():
        L = r['exact_lfd']
        row = dict(cell=c['cell'], cls=c['cls'], n=c['n'], k=c['k'], e=c['e'], Z=c['Z'], lfd=L,
                   status=r['status'], falls=r['strict_falls'], cert=r.get('certificate_degree'),
                   seconds=r['closure_seconds'], rss=r['peak_rss_mb'])
        if c['cls'] == 'affine':
            for law in ('A', 'B', 'A_Z'):
                p = c.get('pred_' + law)
                if p is None:
                    continue
                row['pred_' + law] = p
                cells[c['cell']][law]['hit' if L == p else 'miss'] += 1
        elif c['cls'] == 'const':
            cells[c['cell']]['const<=2']['hit' if L <= 2 else 'miss'] += 1
        else:
            cells[c['cell']]['order-two (unscored)'][L] += 1
        detail.append(row)
    out = dict(missing_or_capped=missing, mismatched=mismatch, distinct_systems=len(seen),
               by_cell={c: {law: dict(v) for law, v in d.items()} for c, d in sorted(cells.items())},
               totals={law: dict(sum((cells[c][law] for c in cells if c in ('P1', 'P2', 'P3', 'P4')),
                                     collections.Counter())) for law in ('A', 'B', 'A_Z')},
               control_P5={law: dict(cells['P5'][law]) for law in ('A', 'B', 'A_Z')},
               systems=sorted(detail, key=lambda d: (d['cell'], d['cls'], d['Z'])))
    json.dump(out, open(os.path.join(HERE, 'score.json'), 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('missing_or_capped', 'mismatched', 'distinct_systems', 'by_cell',
                                          'totals', 'control_P5')}, indent=1))


if __name__ == '__main__':
    main()
