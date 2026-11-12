#!/usr/bin/env python3
import json, glob, csv, collections
rows = []
for f in sorted(glob.glob('final_*.jsonl')) + sorted(glob.glob('bonus_k*.jsonl')):
    grp = 'bonus' if f.startswith('bonus') else 'required'
    for l in open(f):
        r = json.loads(l)
        r['group'] = grp
        rows.append(r)
with open('all_compare.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['group', 'system_hash', 'family', 'vars', 'repo_lfd', 'mine_lfd', 'repo_points', 'mine_points', 'repo_status', 'mine_status',
                'repo_strict_falls', 'mine_falls', 'certified_at_D', 'match_all', 'seconds', 'maxrss_mb'])
    for r in rows:
        w.writerow([r['group'], r['system_hash'], r['label'], r['vars'], r['repo']['exact_lfd'], r['mine']['lfd'], r['repo']['point_count'],
                    r['mine']['points'], r['repo']['status'], r['mine']['status'], r['repo']['strict_falls'], r['mine']['falls'],
                    r['mine']['certified_at'], all(r['match'].values()), r['seconds'], r['maxrss_mb']])
cnt = collections.Counter()
for r in rows:
    fam = r['label'].split(' p=')[0].replace('binary(binary)', 'binary').replace('binary(koblitz)', 'koblitz')
    cnt[(r['group'], fam, r['vars'])] += 1
    if not all(r['match'].values()):
        cnt[('MISMATCH', r['system_hash'])] += 1
for k, v in sorted(cnt.items(), key=str):
    print(k, v)
print('total', len(rows), 'mismatching', sum(1 for r in rows if not all(r['match'].values())))
print('max seconds', max(r['seconds'] for r in rows), 'max rss MB', max(r['maxrss_mb'] for r in rows))
