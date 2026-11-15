#!/usr/bin/env python3
import glob, json, os, collections
D = os.path.join(os.path.dirname(__file__), '..', 'data')
rows = []
for f in sorted(glob.glob(D + '/scan_k*.jsonl')):
    for l in open(f):
        rows.append(json.loads(l))
_seen = set(); _u = []
for r in rows:   # small primes: sampling with replacement repeats (p,a,b,t); keep distinct instances only
    kk = (r['k'], r['p'], r['a'], r['b'], r['t'])
    if kk in _seen: continue
    _seen.add(kk); _u.append(r)
rows = _u
summ = {}
exc = []
for k in sorted({r['k'] for r in rows}):
    for key in ('generic', 'locus'):
        for p in sorted({r['p'] for r in rows if r['k'] == k}) + ['all']:
            sel = [r for r in rows if r['k'] == k and (p == 'all' or r['p'] == p) and (r['locus5'] if key == 'locus' else not r['locus5'])]
            if not sel: continue
            units = [r for r in sel if r['unit']]
            n = len(units)
            full = [r for r in units if r['maxwt_ginv'] == 3 * k]
            e = [r for r in units if r['maxwt_ginv'] != 3 * k]
            summ[f'k{k}_{key}_p{p}'] = dict(k=k, kind=key, p=p, instances=len(sel), units=n, unit_fraction=round(n / len(sel), 3),
                delta_counts=dict(collections.Counter(r['delta'] for r in units)),
                n_sd_eq_delta_plus_3k=len(full), frac=(len(full) / n if n else None),
                n_exceptions=len(e), frac_vs_bounds=(None if not n else dict(frac=round(len(full)/n,4), ge_1_minus_3_over_p=(len(full)/n >= 1-3/p) if p!='all' else None, ge_1_minus_10_over_p=(len(full)/n >= 1-10/p) if p!='all' else None)), pred_sd_set=sorted({r['predicted_sd'] for r in units}),
                n_top_coeff_zero=sum(r['top_coeff_zero'] for r in units),
                max_seconds=max(r['seconds'] for r in sel))
            if p != 'all':
                for r in e:
                    exc.append(dict(k=k, p=r['p'], a=r['a'], b=r['b'], t=r['t'], kind=key, delta=r['delta'], maxwt=r['maxwt_ginv'],
                                    predicted_sd=r['predicted_sd'], generic_sd=r['generic_sd'], top_coeff_zero=r['top_coeff_zero']))
bf = []
for f in sorted(glob.glob(D + '/bruteforce_k*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        bf.append({x: r[x] for x in ('k', 'p', 'a', 'b', 't', 'delta', 'locus5', 'maxwt_ginv', 'predicted_sd', 'sd_bruteforce', 'd_ff', 'd_lf', 'sd_matches_prediction')})
json.dump(dict(summary=summ, exceptions=exc, bruteforce=bf), open(D + '/summary.json', 'w'), indent=1)
print('| k | prime | kind | instances | units | delta set | sd=delta+3k | exceptions | pred sd values | >=1-3/p / >=1-10/p |')
print('|---|---|---|---|---|---|---|---|---|---|')
for key, s in summ.items():
    if s['p'] == 'all': continue
    print(f"| {s['k']} | {s['p']} | {s['kind']} | {s['instances']} | {s['units']} | {s['delta_counts']} | {s['n_sd_eq_delta_plus_3k']}/{s['units']} | {s['n_exceptions']} | {s['pred_sd_set']} | {(s['frac_vs_bounds'] or {}).get('ge_1_minus_3_over_p')} / {(s['frac_vs_bounds'] or {}).get('ge_1_minus_10_over_p')} |")
print()
print('n exceptions (distinct):', len(exc)); print('exception maxwt vs 3k-1:', collections.Counter((x['k'], x['maxwt'] - (3*x['k']-1)) for x in exc))
print('brute force:', len(bf), 'matches', sum(b['sd_matches_prediction'] for b in bf))
