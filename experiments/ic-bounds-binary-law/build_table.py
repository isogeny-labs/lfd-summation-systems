#!/usr/bin/env python3
"""Build one table of the certified two-summand binary/Koblitz LFD data.

Sources (read only, never modified):
  experiments/ic-curve-sweep/{cases,results,systems}.jsonl   k = 2..7 (kind binary/koblitz, m = 2)
  experiments/ic-bounds-binary/results-m2.jsonl              k = 8..10
  experiments/ic-binary-inverse/profiles.json               k = n, order-two target, full field

Output: table.jsonl (one row per *case*; `system_hash` identifies the system; the
equations depend on (n, V, b, x_R) only, so Koblitz (0,1) and (1,1) share systems).

Per-degree data: dims[D] = dim V_D (closure at budget D, Boolean ring), for
D = 0..D* (certificate degree); codims[D] = |M_<=D| - dims[D].  For the inverse
family the recorded ranks are of a different (Moore-coordinate) model and are
not copied (dims = None).
"""
import json, math, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'table.jsonl')


def M(r, D):
    return sum(math.comb(r, i) for i in range(0, min(r, D) + 1))


def rows_sweep():
    p = lambda f: os.path.join(ROOT, 'experiments/ic-curve-sweep', f)
    res = {}
    for l in open(p('results.jsonl')):
        x = json.loads(l); res[x['system_hash']] = x
    sysz = {}
    for l in open(p('systems.jsonl')):
        x = json.loads(l)
        if x['kind'] != 'prime':
            sysz[x['system_hash']] = (len(x['coordinate_solutions']), x['variables'])
    for l in open(p('cases.jsonl')):
        c = json.loads(l)
        if c['kind'] == 'prime' or c['m'] != 2:
            continue
        r = res[c['system_hash']]
        z, nv = sysz[c['system_hash']]
        assert z == r['point_count'] and nv == 2 * c['k']
        yield dict(source='ic-curve-sweep', case_id=c['case_id'], system_hash=c['system_hash'],
                   n=c['n'], k=c['k'], a=c['a'], b=c['b'],
                   curve_class='koblitz' if c['kind'] == 'koblitz' else 'ordinary',
                   base_kind=c['base_kind'], target_class=c['target_class'],
                   target_x=c['target'][0], Z=r['point_count'],
                   rational_unordered=c['rational_unordered_multisets'],
                   status=r['status'], lfd=r['exact_lfd'], falls=r['strict_falls'],
                   cert_degree=r['scanned_through'], dims=r['quotient_ranks'])


def rows_m2():
    for l in open(os.path.join(ROOT, 'experiments/ic-bounds-binary/results-m2.jsonl')):
        x = json.loads(l)
        yield dict(source='ic-bounds-binary', case_id=x['case_id'], system_hash=x['system_hash'],
                   n=x['n'], k=x['k'], a=x['a'], b=x['b'],
                   curve_class='koblitz' if x['curve_class'] == 'koblitz' else 'ordinary',
                   base_kind=x['base_kind'], target_class=x['target_class'],
                   target_x=x['target'][0], Z=x['solutions'],
                   rational_unordered=x['rational_unordered_multisets'],
                   status=x['status'], lfd=x['exact_lfd'], falls=x['strict_falls'],
                   cert_degree=x['scanned_through'], dims=x['quotient_ranks'])


def rows_inverse():
    for x in json.load(open(os.path.join(ROOT, 'experiments/ic-binary-inverse/profiles.json'))):
        yield dict(source='ic-binary-inverse', case_id='inverse-n%d' % x['n'],
                   system_hash='inverse-n%d' % x['n'], n=x['n'], k=x['n'], a=0, b=1,
                   curve_class='koblitz', base_kind='full-field', target_class='order-two',
                   target_x=0, Z=x['solutions'], rational_unordered=None,
                   status='exact-proper-model', lfd=x['exact_lfd'], falls=x['strict_falls'],
                   cert_degree=x['upper_criterion_at'], dims=None)


def main():
    out = []
    for row in list(rows_sweep()) + list(rows_m2()) + list(rows_inverse()):
        r = 2 * row['k']
        row['r'] = r
        row['unit'] = row['Z'] == 0
        row['proper_subspace'] = row['k'] < row['n']
        row['codim_n_minus_k'] = row['n'] - row['k']
        row['two_k_minus_n'] = 2 * row['k'] - row['n']
        if row['dims'] is not None:
            row['codims'] = [M(r, D) - d for D, d in enumerate(row['dims'])]
            if not row['unit']:
                assert row['codims'][-1] == row['Z'], row['case_id']
        else:
            row['codims'] = None
        out.append(row)
    with open(OUT, 'w') as f:
        for row in out:
            f.write(json.dumps(row, sort_keys=True) + '\n')
    print('rows', len(out), 'systems', len({r['system_hash'] for r in out}))


if __name__ == '__main__':
    main()
