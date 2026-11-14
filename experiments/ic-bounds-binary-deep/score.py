#!/usr/bin/env python3
"""Score the binary deep runs against the frozen predictions.

Reads predictions-2026-11-12.json (checks its sha256 against the value frozen in PLAN.md) and the per-system
records copied from the pod (fetched/<name>.json).  For every case it checks that the record is about the predicted
system (system_hash, points_sha256, n, k), that the engine file hash in the record equals csd/ic_bounds_bin_big.py,
and that the record is certified (status exact-*, closed_all, a certificate degree).  Then it scores each law on the
affine class (hit iff the certified LFD equals the prediction), units (|Z| = 0) reported separately, and the
inverse-family value on the scaled-subfield cases.  Writes outcomes-2026-11-12.json and prints a Markdown table.
usage: score.py [--fetched DIR] [--out FILE]
"""
import argparse
import collections
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
FROZEN = '7feee5dbafbde72b4b8ddf94b310ee4b8b56a45446b886bf61b8f98c038e33e1'
LAWS = ['B', 'posthoc', 'A', 'A_Z']


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fetched', default=os.path.join(HERE, 'fetched'))
    ap.add_argument('--out', default=os.path.join(HERE, 'outcomes-2026-11-12.json'))
    a = ap.parse_args()
    ppath = os.path.join(HERE, 'predictions-2026-11-12.json')
    psha = sha(ppath)
    assert psha == FROZEN, 'predictions file changed: ' + psha
    pred = json.load(open(ppath))
    engine_sha = sha(os.path.join(ROOT, 'csd', 'ic_bounds_bin_big.py'))
    rows, problems = [], []
    for c in pred['cases']:
        path = os.path.join(a.fetched, c['name'] + '.json')
        if not os.path.exists(path):
            rows.append(dict(name=c['name'], cell=c['cell'], status='no record'))
            problems.append(c['name'] + ': no record')
            continue
        r = json.load(open(path))
        checks = dict(system_hash=r.get('system_hash') == c['system_hash'],
                      points_sha256=r.get('points_sha256') == c['points_sha256'],
                      n_k=(r.get('n'), r.get('k')) == (c['n'], c['k']), n_points=r.get('n_points') == c['Z'],
                      engine_sha256=r.get('engine_sha256') == engine_sha,
                      certified=str(r.get('status', '')).startswith('exact') and bool(r.get('closed_all'))
                      and r.get('certificate_degree') is not None and r.get('exact_lfd') is not None)
        for key, ok in checks.items():
            if not ok:
                problems.append('%s: check %s failed' % (c['name'], key))
        row = dict(name=c['name'], cell=c['cell'], n=c['n'], k=c['k'], e=c['e'], k_plus_e=c['k_plus_e'], Z=c['Z'],
                   cls=c['cls'], base=c['base'], curve=c['curve'], subfield=bool(c.get('subfield')),
                   trace_zero=bool(c.get('trace_zero')), status=r.get('status'), lfd=r.get('exact_lfd'),
                   certificate_degree=r.get('certificate_degree'), strict_falls=r.get('strict_falls'),
                   repairs=r.get('repairs'), seconds=r.get('seconds'), columns=r.get('columns'),
                   rref_sha256=r.get('rref_sha256'), record_sha256=sha(path), checks=checks)
        if c['cls'] == 'affine':
            row['pred'] = {law: c['pred_' + law] for law in LAWS}
            row['hit'] = {law: (r.get('exact_lfd') == c['pred_' + law]) for law in LAWS}
        if c.get('subfield'):
            inv = max(3, c['k'] // 2 + 1)
            row['pred_inverse_family'] = inv
            row['hit_inverse_family'] = r.get('exact_lfd') == inv
        rows.append(row)
    done = [r for r in rows if r.get('lfd') is not None]
    aff = [r for r in done if r['cls'] == 'affine']
    score = {}
    for label, sub in (('affine, all', aff), ('affine, |Z| > 0', [r for r in aff if r['Z'] > 0]),
                       ('affine, |Z| = 0 (units)', [r for r in aff if r['Z'] == 0])):
        score[label] = dict(systems=len(sub), **{law: sum(r['hit'][law] for r in sub) for law in LAWS})
    sub = [r for r in done if r.get('subfield')]
    score['scaled subfield, inverse-family value'] = dict(systems=len(sub), hits=sum(r['hit_inverse_family'] for r in sub))
    cells = collections.OrderedDict()
    for r in rows:
        cells.setdefault(r['cell'], []).append(r)
    out = dict(predictions_sha256=psha, engine_sha256=engine_sha, cases=len(pred['cases']), records=len(done),
               problems=problems, score=score, rows=rows)
    json.dump(out, open(a.out, 'w'), indent=1)
    print('predictions sha256 %s (frozen value: equal)' % psha)
    print('records %d of %d, problems %d' % (len(done), len(pred['cases']), len(problems)))
    for p in problems:
        print('PROBLEM', p)
    print()
    print('| cell | (n, k) | e | k+e | base / curve | class | systems | \\|Z\\| | certified LFD | cert. degree | B | post-hoc | A | A_Z | inverse |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for cell, rs in cells.items():
        groups = collections.OrderedDict()
        for r in rs:
            groups.setdefault((r.get('cls'), r.get('subfield')), []).append(r)
        for (cls, subf), g in groups.items():
            g = [r for r in g if r.get('lfd') is not None]
            if not g:
                continue
            def span(vals):
                vals = sorted(set(vals))
                return str(vals[0]) if len(vals) == 1 else '%s-%s' % (vals[0], vals[-1])
            lf = collections.Counter(r['lfd'] for r in g)
            lfs = ', '.join('%d (x%d)' % (v, n) for v, n in sorted(lf.items()))
            pr = lambda law: ('%s: %d/%d' % (span([r['pred'][law] for r in g]), sum(r['hit'][law] for r in g), len(g))) if cls == 'affine' else '-'
            inv = ('%s: %d/%d' % (span([r['pred_inverse_family'] for r in g]), sum(r['hit_inverse_family'] for r in g), len(g))) if subf else '-'
            print('| %s | (%d, %d) | %s | %s | %s / %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                cell, g[0]['n'], g[0]['k'], span([r['e'] for r in g]), span([r['k_plus_e'] for r in g]),
                '+'.join(sorted(set(r['base'] for r in g))), '+'.join(sorted(set(r['curve'] for r in g))), cls, len(g),
                span([r['Z'] for r in g]), lfs, span([r['certificate_degree'] for r in g]),
                pr('B'), pr('posthoc'), pr('A'), pr('A_Z'), inv))
    print()
    for label, sc in score.items():
        print(label, sc)


if __name__ == '__main__':
    main()
