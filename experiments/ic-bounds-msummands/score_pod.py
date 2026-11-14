#!/usr/bin/env python3
"""Score the pod runs (gpu/RUN.md) against the frozen predictions (G4, G5, G7, G8 of predictions-2026-11-12.json).
Reads pod-records/<chain>/out/*.json (unpacked from fetched-p2/msum-*.tgz) and gpu/cost.json (frozen hashes, d_in, r).
Checks every record's generator/system hash against cost.json.  Writes pod-scores.json and prints a Markdown table."""
import glob, hashlib, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, 'predictions-2026-11-12.json')))
assert hashlib.sha256(open(os.path.join(HERE, 'predictions-2026-11-12.json'), 'rb').read()).hexdigest().startswith('deb4b304')
cost = json.load(open(os.path.join(HERE, 'gpu', 'cost.json')))
meta = {r['name']: dict(r, family=f) for f in ('binary', 'prime') for r in cost[f]}
rows = []
for chain_status in sorted(glob.glob(os.path.join(HERE, 'pod-records', '*', 'out', 'results.jsonl'))):
    chain = chain_status.split(os.sep)[-3]
    for line in open(chain_status):
        res = json.loads(line); name = res['name']; m = meta[name]
        rec_path = os.path.join(os.path.dirname(chain_status), name + '.json')
        rec = json.load(open(rec_path)) if os.path.exists(rec_path) else {}
        # binary: the ladder record's system_hash is the GPU engine's own hash, a different definition from the frozen
        # msum system_hash; chain_msum.py compared the rebuilt instance's msum_system_hash with the frozen one before the
        # run and would have recorded 'instance mismatch' instead of running.  So: frozen hash present -> checked on the pod.
        if m['family'] == 'binary':
            want = m.get('system_hash'); got = want if (want and res['status'] != 'instance mismatch or build failure') else None
        else:
            want = m.get('generator_hash'); got = rec.get('generator_hash')
        row = dict(name=name, chain=chain, family=m['family'], m=m['m'], k=m['k'], run_status=res['status'],
                   engine_status=res.get('engine_status'), lfd=res.get('exact_lfd'), strict_falls=res.get('strict_falls'),
                   seconds=res.get('seconds'), hash_equal=('no frozen hash' if want is None else (got == want)) if rec else None,
                   lower_bound=rec.get('lfd_lower_bound'))
        if m['family'] == 'prime':
            d = m.get('generator_degree') or rec.get('generator_degree'); row.update(generator_degree=d, mk1=m['m'] * m['k'] + 1, solutions=m.get('solutions'))
            if m['m'] == 4 and m['k'] == 4:
                row['G7_range_ok'] = row['lfd'] is not None and 13 <= row['lfd'] <= 17
                row['G7_point_b'] = row['lfd'] == 13
            if m['m'] == 5 and m['k'] == 3:
                row['G8_hit'] = row['lfd'] == 16
        else:
            r, din = m['r'], m['input_degree']
            row.update(r=r, d_in=din, d_in_source=m.get('input_degree_source'), solutions=m.get('solutions'))
            if row['lfd'] is not None:
                row['G4_upper_ok'] = row['lfd'] <= min(r, din + 2)
                row['point_rule'] = min(r, din + 1); row['point_hit'] = row['lfd'] == row['point_rule']
        rows.append(row)
json.dump(rows, open(os.path.join(HERE, 'pod-scores.json'), 'w'), indent=1)
print('| case | family | m | k | deg g or d_in | r | \\|Z\\| | status | LFD | strict falls | s | hash equal | predictions |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for x in rows:
    pred = []
    if 'G7_point_b' in x: pred.append('G7 range %s, (b)=13 %s' % (x['G7_range_ok'], x['G7_point_b']))
    if 'G8_hit' in x: pred.append('G8 (16) %s' % x['G8_hit'])
    if 'G4_upper_ok' in x: pred.append('upper %s, point %s %s' % (x['G4_upper_ok'], x['point_rule'], x['point_hit']))
    print('| %s | %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        x['name'], x['family'], x['m'], x['k'], x.get('generator_degree', x.get('d_in')), x.get('r', ''), x.get('solutions'),
        x['engine_status'] or x['run_status'], x['lfd'] if x['lfd'] is not None else ('>= %s' % x['lower_bound'] if x['lower_bound'] else '-'),
        x['strict_falls'], x['seconds'], x['hash_equal'], '; '.join(pred)))
