#!/usr/bin/env python3
"""Inventory: data/binary_table.jsonl (lane C).  Read-only on the repository records.

Rows: every case of results-{m2,m3-direct,m3-split}.jsonl (ic-bounds-binary, 1871 cases) plus the
two-summand rows k = 2..7 of ic-bounds-binary-law/table.jsonl (source ic-curve-sweep; no related-*.jsonl exists
for them).  Columns: family, n, k, r, e (m2 only), Z, class, LFD, D*, d_lin, d_full / d_full_lb (lower bound),
d_ff (= min of the recorded strict falls), d_reg (semi-regular, m2 only), kernel bound r+1, provenance.
d_full / d_lin come from related-*.jsonl (cap) and are replaced by tools/static_scan.py results when available
(data/static_scan_*.jsonl, keyed by system hash).
"""
import json, glob, os, math, collections, sys
REPO = '/home/user/ml-cryptanalysis'
HERE = os.path.dirname(os.path.abspath(__file__)); LANE = os.path.dirname(HERE)
B = os.path.join(REPO, 'experiments/ic-bounds-binary'); LAW = os.path.join(REPO, 'experiments/ic-bounds-binary-law')


def jl(p):
    return [json.loads(l) for l in open(p) if l.strip()]


def coeffs(r, m, N):
    return [sum(math.comb(r, i - 2 * j) * (-1) ** j * math.comb(m + j - 1, j) for j in range(i // 2 + 1)) for i in range(N + 1)]


def d_reg(r, m):   # same as analysis.semi_regular_dreg (Boolean series (1+t)^r/(1+t^2)^m, first non-positive coefficient)
    for i, c in enumerate(coeffs(r, m, r + 2)):
        if c <= 0:
            return i
    return r + 1


def main():
    struct = {x['system_hash']: x for x in jl(os.path.join(LAW, 'structure.jsonl'))}
    scans = {}
    for f in sorted(glob.glob(os.path.join(LANE, 'data/static_scan_*.jsonl'))):
        for x in jl(f):
            old = scans.get(x['system_hash'])
            if old is None or (x['d_full'] is not None) or ((x['d_full_lb'] or 0) > (old['d_full_lb'] or 0) and old['d_full'] is None):
                scans[x['system_hash']] = x
    out, seen_case = [], set()
    for fam in ['m2', 'm3-direct', 'm3-split']:
        rel = {x['case_id']: x for x in jl(os.path.join(B, f'related-{fam}.jsonl'))}
        for r in jl(os.path.join(B, f'results-{fam}.jsonl')):
            s = struct.get(r['system_hash']) if fam == 'm2' else None
            row = dict(family=fam, case_id=r['case_id'], system_hash=r['system_hash'], n=r['n'], k=r['k'], r=r['variables'],
                       e=(s['r_eff'] - s['m_eff']) if s else None, r_eff=s['r_eff'] if s else None, m_eff=s['m_eff'] if s else None,
                       Z=r['solutions'], status=r['status'], unit=(r['solutions'] == 0),
                       cls=r['target_class'], curve_class=r['curve_class'], base_kind=r['base_kind'],
                       lfd=r['exact_lfd'], dstar=r['certificate_degree'], d_ff=(min(r['strict_falls']) if r['strict_falls'] else None),
                       strict_falls=r['strict_falls'], d_lin=None, d_full=None, d_full_lb=None, d_full_src=None,
                       d_reg=d_reg(s['r_eff'], s['m_eff']) if s else None, kernel_r_plus_1=r['variables'] + 1, source='ic-bounds-binary')
            if r['case_id'] in rel:
                q = rel[r['case_id']]
                row['d_lin'] = q['d_lin']; row['d_full'] = q['d_full']
                row['d_full_lb'] = q['capped_at'] if q['d_full'] is None else None
                row['d_full_src'] = f'related-{fam}.jsonl (cap {q["capped_at"]})' if q['d_full'] is None else f'related-{fam}.jsonl'
            out.append(row)
    # k = 2..7 two-summand rows
    for t in jl(os.path.join(LAW, 'table.jsonl')):
        if t['source'] != 'ic-curve-sweep':
            continue
        s = struct[t['system_hash']]
        out.append(dict(family='m2', case_id=t['case_id'], system_hash=t['system_hash'], n=t['n'], k=t['k'], r=t['r'],
                        e=s['r_eff'] - s['m_eff'], r_eff=s['r_eff'], m_eff=s['m_eff'], Z=t['Z'], status=t['status'], unit=t['unit'],
                        cls=t['target_class'], curve_class=t['curve_class'], base_kind=t['base_kind'], lfd=t['lfd'], dstar=t['cert_degree'],
                        d_ff=(min(t['falls']) if t['falls'] else None), strict_falls=t['falls'], d_lin=None, d_full=None, d_full_lb=None,
                        d_full_src=None, d_reg=d_reg(s['r_eff'], s['m_eff']), kernel_r_plus_1=t['r'] + 1, source='ic-curve-sweep'))
    for row in out:
        sc = scans.get(row['system_hash'])
        if sc and row['family'] == 'm2':
            lin = [p['D'] for p in sc['per_degree'] if False]
            if sc['d_full'] is not None:
                row['d_full'], row['d_full_lb'] = sc['d_full'], None
                row['d_full_src'] = 'tools/static_scan.py'
            elif row['d_full'] is None and (sc['d_full_lb'] or 0) > (row['d_full_lb'] or 0):
                row['d_full_lb'] = sc['d_full_lb']; row['d_full_src'] = 'tools/static_scan.py (lower bound)'
            row['scan_deficiency_last'] = sc['per_degree'][-1]['deficiency'] if sc['per_degree'] else None
            row['scan_last_D'] = sc['per_degree'][-1]['D'] if sc['per_degree'] else None
    with open(os.path.join(LANE, 'data/binary_table.jsonl'), 'w') as fh:
        for row in out:
            fh.write(json.dumps(row, sort_keys=True) + '\n')
    c = collections.Counter(r['family'] for r in out)
    print(len(out), dict(c))


if __name__ == '__main__':
    main()
