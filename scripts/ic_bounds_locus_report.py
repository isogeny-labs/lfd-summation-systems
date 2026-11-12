"""Summary tables, figure and hashes for the exceptional-locus map (bounds item 2).

Reads experiments/ic-bounds-locus/{data,heldout,plane,predictions,validation}
and writes experiments/ic-bounds-locus/summary.json, plane-map.svg and
plane-map.png, and verification.json (sha256 of every source and data file).

Every target is assigned a rule label from the frozen rules:
  t0-j0, t0-j1728, t0-other, tangent (3t^2+a = 0, t != 0), order-two (rhs = 0),
  C1 (cubic in the normalized plane, V_4 three below generic),
  C2 (normalized point in the C2 class of the p = 193 plane, k = 6 only),
  none.
"""
from __future__ import annotations

import glob
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from ic_bounds_locus_classify import summarize, deviations, generic_values  # noqa: E402

EXP = ROOT / 'experiments/ic-bounds-locus'


def load_dir(sub, prime, k):
    recs = []
    for f in sorted(glob.glob(str(EXP / sub / f'p{prime}-k{k}-*.jsonl'))):
        recs += [json.loads(line) for line in open(f)]
    return recs


def plane_classes(prime, k):
    out = {}
    for f in sorted(glob.glob(str(EXP / 'plane' / f'p{prime}-k{k}-plane-c*.jsonl'))):
        for line in open(f):
            r = json.loads(line)
            out[(r['A'], r['B'])] = r
    return out


def cubic(prime):
    lift = json.loads((EXP / 'c1-lift.json').read_text())
    c = lift['cubics'].get(str(prime))
    return None if c is None else {tuple(int(x) for x in key.split(',')): v for key, v in c['coefficients'].items()}


def label(p, a, b, t, C1=None, plane=None):
    if t == 0:
        return 't0-j0' if a == 0 else 't0-j1728' if b == 0 else 't0-other'
    if (3 * t * t + a) % p == 0:
        return 'tangent'
    if (t ** 3 + a * t + b) % p == 0:
        return 'order-two'
    it = pow(t, -1, p)
    A, B = a * it * it % p, b * it ** 3 % p
    if C1 is not None and sum(c * pow(A, i, p) * pow(B, j, p) for (i, j), c in C1.items()) % p == 0:
        return 'C1'
    if plane is not None and (A, B) in plane:
        r = plane[(A, B)]
        if r['V3'] == 13 and r['V4'] > 700:
            return 'C2'
    return 'none'


def table(recs, k, C1, plane):
    rows = []
    for r in recs:
        s = summarize(r)
        dev = deviations(s, k)
        lab = label(r['prime'], r['a'], r['b'], r['t'], C1, plane)
        rows.append(dict(a=r['a'], b=r['b'], t=r['t'], label=lab, deviation=dev, V3=s['V3'], V4=s['V4'], low4=s['low4'],
                         one4=s['one4'], Z=s['Z'], V5=s.get('V5'), exact_lfd=s['exact_lfd'], status=s['status']))
    by = defaultdict(list)
    for x in rows:
        by[x['label']].append(x)
    classes = {}
    for lab, xs in sorted(by.items()):
        classes[lab] = dict(count=len(xs), deviating=sum(bool(x['deviation']) for x in xs),
                            signatures=Counter(json.dumps(dict(V3=x['V3'], V4=x['V4'], low4=x['low4'],
                                                               one4=x['one4'], lfd=x['exact_lfd']))
                                               for x in xs).most_common(8))
    unexplained = [x for x in rows if x['deviation'] and x['label'] == 'none']
    false_alarms = [x for x in rows if not x['deviation'] and x['label'] != 'none']
    return rows, classes, unexplained, false_alarms


def figure(plane, prime, k, out_svg, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups = defaultdict(list)
    for (A, B), r in plane.items():
        g = generic_values(k)
        if (r['V3'], r['V4'], r['low4']) == (g['V3'], g['V4'], g['low4']):
            continue
        if r['V3'] != g['V3']:
            name = 'A = -3 (3t^2+a = 0): V3 = 6k-2'
        elif r['V4'] == g['V4'] + 2:
            name = '1+A+B = 0 (order two): V4 +2'
        elif r['V4'] == g['V4'] - 3:
            name = 'C1 cubic: V4 -3'
        elif r['V4'] > 700:
            name = 'C2: V4 = |M<=4|-|Z| (saturated)'
        else:
            name = f'other V4={r["V4"]}'
        groups[name].append((A, B))
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    colors = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#8c564b', '#ff7f0e']
    for (name, pts), c in zip(sorted(groups.items()), colors):
        ax.scatter([a for a, _ in pts], [b for _, b in pts], s=7, color=c, label=f'{name} ({len(pts)})')
    lines = sorted({A for A, _ in plane})
    ax.set_xlim(-2, prime + 1)
    ax.set_ylim(-2, prime + 1)
    ax.set_xlabel('A = a / t^2  (mod p)')
    ax.set_ylabel('B = b / t^3  (mod p)')
    ax.set_title(f'Exceptional points of the closure, normalized plane, p = {prime}, k = {k}\n'
                 f'({len(plane)} points computed on {len(lines)} lines A = const; generic points not drawn)', fontsize=9)
    ax.legend(fontsize=7, loc='upper right', framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=110)


def main():
    summary = {}
    plane193 = plane_classes(193, 6)
    C1_193, C1_257 = cubic(193), cubic(257)
    for prime, k, C1, plane in ((193, 6, C1_193, plane193), (257, 7, C1_257, None)):
        for sub in ('data', 'heldout'):
            recs = load_dir(sub, prime, k)
            if not recs:
                continue
            rows, classes, unexplained, false_alarms = table(recs, k, C1, plane)
            key = f'p{prime}-k{k}-{sub}'
            summary[key] = dict(targets=len(rows), curves=sorted({(r['a'], r['b']) for r in rows}),
                                deviating=sum(bool(r['deviation']) for r in rows),
                                lfd_histogram=dict(Counter(str(r['exact_lfd']) for r in rows)),
                                status_histogram=dict(Counter(r['status'] for r in rows)),
                                classes=classes,
                                unexplained=[dict(a=x['a'], b=x['b'], t=x['t'], V3=x['V3'], V4=x['V4'], low4=x['low4'],
                                                  lfd=x['exact_lfd']) for x in unexplained],
                                false_alarms=[dict(a=x['a'], b=x['b'], t=x['t'], label=x['label']) for x in false_alarms])
            print(key, len(rows), 'deviating', summary[key]['deviating'], 'unexplained', len(unexplained),
                  'false alarms', len(false_alarms))
    if plane193:
        # (A, B)-invariance test: every k = 6 target with t != 0 against the t = 1 plane point
        inv = Counter()
        mism = []
        for sub in ('data', 'heldout'):
            for r in load_dir(sub, 193, 6):
                if r['t'] == 0:
                    continue
                s = summarize(r)
                it = pow(r['t'], -1, 193)
                q = plane193[(r['a'] * it * it % 193, r['b'] * it ** 3 % 193)]
                sat_r = s['V4'] + s['Z'] == 794 and s['V4'] > 700
                sat_q = q['V4'] + q['Z'] == 794 and q['V4'] > 700
                same = (s['V3'], s['low3']) == (q['V3'], q['low3']) and sat_r == sat_q and (
                    sat_r or (s['V4'], s['low4']) == (q['V4'], q['low4']))
                t_in_H = pow(r['t'], 64, 193) == 1
                inv[(sub, 't in H' if t_in_H else 't not in H', 'agree' if same else 'DISAGREE')] += 1
                if not same:
                    mism.append(dict(sub=sub, a=r['a'], b=r['b'], t=r['t']))
        summary['AB_invariance_k6'] = dict(counts={' / '.join(k): v for k, v in sorted(inv.items())}, mismatches=mism)
        print('AB invariance', dict(inv), mism[:5])
        figure(plane193, 193, 6, EXP / 'plane-map.svg', EXP / 'plane-map.png')
        summary['plane_p193_k6'] = dict(points=len(plane193), lines=len({A for A, _ in plane193}),
                                        classes=Counter(json.dumps(dict(V3=r['V3'], V4=r['V4'] if r['V4'] < 700 else '>700'))
                                                        for r in plane193.values()).most_common())
    (EXP / 'summary.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    files = sorted(glob.glob(str(ROOT / 'scripts/ic_bounds_locus*')) + [str(ROOT / 'csd/ic_bounds_locus.py')])
    data = sorted(f for f in glob.glob(str(EXP / '**/*'), recursive=True)
                  if Path(f).is_file() and not f.endswith('verification.json'))
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
    engines = ['csd/ic_fall_certificates.py', 'csd/ic_chain_shift_closure.py', 'csd/ic_curve_sweep.py', 'csd/ic_dreg.py',
               'csd/elim.py', 'csd/_elim.c', 'experiments/ic-curve-validation/code/lfd_engine.py']
    (EXP / 'verification.json').write_text(json.dumps(dict(
        sources={str(Path(f).relative_to(ROOT)): sha(f) for f in files},
        engines={f: sha(ROOT / f) for f in engines},
        data={str(Path(f).relative_to(ROOT)): sha(f) for f in data}), indent=1) + '\n')


if __name__ == '__main__':
    main()
