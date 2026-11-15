"""Collect data/scan_*.json, dlf_crosscheck_*.json, dff_*.json into data/all_instances.json and data/summary.json.
Law/recorded d_lf values are labelled as such.  [pcc:unverified]"""
import json, glob, os, collections
D = os.path.join(os.path.dirname(__file__), '..', 'data')
rows, groups = [], []
for f in sorted(glob.glob(D + '/scan_*.json')):
    d = json.load(open(f))
    groups.append(dict(k=d['k'], p=d['p'], seed=d['seed'], tries=d['tries'], n_unit=d['n_unit'], unit_rate=d['unit_rate'],
                       wall_sec=d['wall_sec'], sec_inverse=d['sec_inverse'], sec_bruteforce=d['sec_bruteforce'], sec_dlf=d['sec_dlf'], file=os.path.basename(f)))
    rows += d['rows']
json.dump(rows, open(D + '/all_instances.json', 'w'))
# recorded values: ml-cryptanalysis/experiments/ic-bounds-exact10/summary.json profile_lfd_by_k (5:4, 6:5, 7:5); law floor((k+14)/4)
def law(k): return (k + 14) // 4
by_k = collections.defaultdict(list)
for r in rows: by_k[r['k']].append(r)
summ = {}
for k in sorted(by_k):
    rs = by_k[k]
    meas = [r for r in rs if r['d_lf'] is not None]
    ffs = [r for r in rs if r['d_ff'] is not None]
    # P2 per prime
    per_p = collections.defaultdict(lambda: [0, 0, collections.Counter()])
    for r in rs:
        x = per_p[r['p']]; x[0] += 1; x[1] += (r['sd_inverse'] == 2 * k + 2); x[2][r['sd_inverse']] += 1
    summ[k] = dict(n=len(rs), p2={p: dict(n=v[0], n_generic=v[1], frac=v[1] / v[0], others=dict(v[2]), one_minus_3_over_p=1 - 3 / p,
                                          one_minus_10_over_p=1 - 10 / p) for p, v in per_p.items()},
                   sd_values=dict(collections.Counter(r['sd_inverse'] for r in rs)),
                   d_lf_measured=dict(collections.Counter(r['d_lf'] for r in meas)), d_ff_measured=dict(collections.Counter(r['d_ff'] for r in ffs)),
                   d_lf_law=law(k), P1_checked=sum(1 for r in rs if r['sd_bruteforce'] is not None),
                   P1_mismatch=sum(1 for r in rs if r.get('P1_match') is False),
                   inverse_unverified=sum(1 for r in rs if not r['inverse_verified']))
json.dump(dict(groups=groups, per_k=summ), open(D + '/summary.json', 'w'), indent=1, default=str)
print('rows', len(rows), 'groups', len(groups))
for k, s in summ.items():
    print(k, s['n'], 'sd', s['sd_values'], 'dlf', s['d_lf_measured'], 'dff', s['d_ff_measured'], 'law', s['d_lf_law'], 'P1', s['P1_checked'], s['P1_mismatch'], 'unverif', s['inverse_unverified'])
