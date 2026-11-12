#!/usr/bin/env python3
"""Compare every producer validation row (k=5,6,7; both batches) with the independent closure
values in experiments/ic-bounds-exact10/verify/grid-*.jsonl.  Reads only; writes verify/compare.json."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
E = ROOT / 'experiments/ic-bounds-exact10'
mine = {}
for f in sorted((E / 'verify').glob('grid-*.jsonl')):
    for line in open(f):
        r = json.loads(line)
        for rec in r['recs']:
            mine[(r['instance'], r['k'], rec['D'])] = dict(rec, points=r['points'], x=r['x'])
out = dict(batches={}, mismatches=[])
tot = dict(rows=0, dim=0, low=0, one=0, prev_dim=0, prev_one=0, block_one=0, block_eq=0, neg_expected=0, neg_detected=0, neg_rows=0, undecided=0, unit_not_one=0, proper=0)
for f in sorted(E.glob('validation-*.jsonl')):
    b = dict(rows=0, agree=0)
    for line in open(f):
        r = json.loads(line)
        m = mine.get((r['instance'], r['k'], r['degree']))
        mp = mine.get((r['instance'], r['k'], r['degree'] - 1))
        if m is None:
            continue
        b['rows'] += 1
        ss, sp, be, nc = r['split_scratch'], r['split_from_prev_sparse'], r['block_engine'], r.get('negative_control')
        ok = {
            'points': m['points'] == r['points'] and m['x'] == r['target_x'],
            'dim': m['dimV'] == r['ref_dim'] == ss['exact_dimension'],
            'low': m['lowV'] == r['ref_low_dim'] == ss['low_dimension'],
            'one': m['one'] == r['ref_contains_one'] == ss['contains_one'],
            'prev_dim': m['dimV'] == sp['exact_dimension'] and m['lowV'] == sp['low_dimension'],
            'prev_one': m['one'] == sp['contains_one'],
            'block_one': m['one'] == be['contains_one'],
            'block_eq': m['one'] or (be['sym_rank'] + be['anti_rank'] == m['dimV']),
            'target': m['target'] == ss['target'],
        }
        if nc is not None and mp is not None and not nc.get('one_before'):
            # producer's control space contains V_{D-1}; it is a proper subspace of V_D at least when
            # dim V_{D-1} < dim V_D (independent value); rank_before is the producer's own count
            ok['neg_expected'] = nc['expected_missing'] == (mp['dimV'] < m['dimV'])
            ok['neg_detected'] = nc['detected'] == nc['expected_missing']
            if 'v2' not in f.name:
                tot['neg_rows'] += 1
        if all(ok.values()):
            b['agree'] += 1
        else:
            out['mismatches'].append(dict(file=f.name, instance=r['instance'], k=r['k'], D=r['degree'],
                                          bad=[k for k, v in ok.items() if not v]))
        if 'v2' not in f.name:
            tot['rows'] += 1
            for key in ('dim', 'low', 'one', 'prev_dim', 'prev_one', 'block_one', 'block_eq',
                        'neg_expected', 'neg_detected'):
                tot[key] += int(ok.get(key, False))
            tot['undecided'] += int(not m['one'])
            tot['unit_not_one'] += int(not m['one'] and m['points'] == 0)
            tot['proper'] += int(m['points'] > 0)
    out['batches'][f.name] = b
out['totals_excluding_v2'] = tot
(E / 'verify/compare.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
