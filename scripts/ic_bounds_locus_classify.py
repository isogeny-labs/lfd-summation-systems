"""Classify the deviations of the two-summand prime chain closure (bounds item 2).

Reads experiments/ic-bounds-locus/data/p<P>-k<K>-*.jsonl (scripts/ic_bounds_locus.py),
marks every target whose closure data differ from the generic values

  dim V_3 = 2k+1, dim V_4 = 4k^2+172k-571, dim(V_4 ∩ M_<=3) = 4k+87, 1 not in V_4,
  and at k = 6: dim V_5 = |M_<=5| - |Z| (V_5 = I ∩ M_<=5) and exact LFD 5,

and searches for polynomial/class conditions in (a, b, t) (csd.ic_bounds_locus.features)
that describe the deviating set exactly: unions of conjunctions of at most two
atoms (atoms and their negations), with zero false positives and zero misses
on every rational target of every curve in the run.

Usage
  python3 scripts/ic_bounds_locus_classify.py --prime 193 --k 6 [--holdout a,b ...] --out FILE.json
"""
from __future__ import annotations

import argparse
import glob
import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd.ic_bounds_locus import features  # noqa: E402


def generic_values(k):
    return dict(V3=2 * k + 1, V4=4 * k * k + 172 * k - 571, low4=4 * k + 87)


def load(prime, k, data_dir):
    recs = []
    for f in sorted(glob.glob(str(Path(data_dir) / f'p{prime}-k{k}-*.jsonl'))):
        for line in open(f):
            recs.append(json.loads(line))
    return recs


def summarize(r):
    d = {row['D']: row for row in r['rows']}
    s = dict(a=r['a'], b=r['b'], t=r['t'], j=r['j'], y_zero=r['y_zero'], Z=r['Z'],
             V3=d[3]['dim'], low3=d[3]['low'], V4=d[4]['dim'], low4=d[4]['low'], one4=d[4]['one'],
             status=r['status'], exact_lfd=r['exact_lfd'], falls=r['falls'])
    if 5 in d:
        s.update(V5=d[5]['dim'], low5=d[5]['low'], codim5=d[5]['codim'], one5=d[5]['one'])
    if 6 in d:
        s.update(V6=d[6]['dim'], low6=d[6]['low'], codim6=d[6]['codim'])
    return s


def deviations(s, k):
    g = generic_values(k)
    dev = {}
    if s['V3'] != g['V3']:
        dev['V3'] = s['V3']
    if s['V4'] != g['V4']:
        dev['V4'] = s['V4']
    if s['low4'] != g['low4']:
        dev['low4'] = s['low4']
    if s['one4']:
        dev['one4'] = True
    if 'V5' in s:
        if s['codim5'] != s['Z']:
            dev['V5_deficiency'] = s['codim5'] - s['Z']
        if s['exact_lfd'] != 5:
            dev['LFD'] = s['exact_lfd']
    return dev


def search(rows, positive_key, max_terms=4):
    """Exact descriptions of {rows with positive_key} as unions of <= max_terms
    conjunctions of <= 2 literals; returns all minimal-size solutions (capped)."""
    atoms = sorted({name for r in rows for name in r['features']})
    literals = [(n, True) for n in atoms] + [(n, False) for n in atoms]
    pos = frozenset(i for i, r in enumerate(rows) if r[positive_key])
    def holds(lits, r):
        return all(r['features'].get(n, False) == v for n, v in lits)
    clauses = []
    for size in (1, 2):
        for lits in itertools.combinations(literals, size):
            if size == 2 and lits[0][0] == lits[1][0]:
                continue
            cover = frozenset(i for i, r in enumerate(rows) if holds(lits, r))
            if cover and cover <= pos:
                clauses.append((lits, cover))
    # keep only maximal covers per cover set (shortest literal list)
    best = {}
    for lits, cover in clauses:
        if cover not in best or len(lits) < len(best[cover]):
            best[cover] = lits
    items = sorted(best.items(), key=lambda kv: -len(kv[0]))
    # drop covers strictly contained in another
    maximal = [(c, l) for c, l in items if not any(c < c2 for c2, _ in items)]
    solutions = []
    for n in range(1, max_terms + 1):
        for combo in itertools.combinations(maximal, n):
            if frozenset().union(*(c for c, _ in combo)) == pos:
                solutions.append([[f'{name}' if v else f'not({name})' for name, v in l] for _, l in combo])
        if solutions:
            break
    return dict(positives=len(pos), total=len(rows), minimal_terms=(len(solutions[0]) if solutions else None),
                solutions=solutions[:40], solution_count=len(solutions))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prime', type=int, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--data', default=str(ROOT / 'experiments/ic-bounds-locus/data'))
    ap.add_argument('--holdout', nargs='*', default=[])
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    hold = {tuple(int(v) for v in c.split(',')) for c in args.holdout}
    recs = load(args.prime, args.k, args.data)
    rows = []
    for r in recs:
        s = summarize(r)
        s['deviation'] = deviations(s, args.k)
        s['features'] = features(r['prime'], r['a'], r['b'], r['t'], r['k'])
        s['holdout'] = (r['a'], r['b']) in hold
        for key in ('V3', 'V4', 'low4', 'one4', 'V5_deficiency', 'LFD'):
            s['dev_' + key] = key in s['deviation']
        s['dev_any'] = bool(s['deviation'])
        rows.append(s)
    train = [r for r in rows if not r['holdout']]
    curves = sorted({(r['a'], r['b']) for r in rows})
    per_curve = []
    for a, b in curves:
        sub = [r for r in rows if (r['a'], r['b']) == (a, b)]
        per_curve.append(dict(a=a, b=b, j=sub[0]['j'], holdout=sub[0]['holdout'], targets=len(sub),
                              order_two_targets=sum(r['y_zero'] for r in sub),
                              deviating=sum(r['dev_any'] for r in sub),
                              deviating_t=[r['t'] for r in sub if r['dev_any']]))
    signatures = defaultdict(list)
    for r in rows:
        if r['dev_any']:
            key = json.dumps({x: r.get(x) for x in ('V3', 'low3', 'V4', 'low4', 'one4', 'exact_lfd')}, sort_keys=True)
            signatures[key].append(r)
    sig_out = []
    for key, members in sorted(signatures.items(), key=lambda kv: -len(kv[1])):
        common = set.intersection(*({n for n, v in m['features'].items() if v} for m in members))
        generic_rate = {n: sum(r['features'].get(n, False) for r in rows if not r['dev_any']) for n in sorted(common)}
        sig_out.append(dict(signature=json.loads(key), count=len(members),
                            members=[dict(a=m['a'], b=m['b'], t=m['t'], Z=m['Z'], deviation=m['deviation'],
                                          true_features=sorted(n for n, v in m['features'].items() if v))
                                     for m in members],
                            common_features=sorted(common), common_feature_counts_among_generic=generic_rate))
    searches = {}
    for key in ('dev_any', 'dev_V3', 'dev_V4', 'dev_low4', 'dev_one4', 'dev_V5_deficiency', 'dev_LFD'):
        if any(r[key] for r in train):
            searches[key] = search(train, key)
    out = dict(prime=args.prime, k=args.k, generic=generic_values(args.k), records=len(rows),
               training_records=len(train), per_curve=per_curve, signatures=sig_out, searches=searches,
               signature_histogram=[[list(key), n] for key, n in
                                     Counter((r['V3'], r['V4'], r['low4']) for r in rows).most_common()])
    Path(args.out).write_text(json.dumps(out, indent=1, default=str) + '\n')
    for pc in per_curve:
        print(pc)
    for s in sig_out:
        print(s['count'], s['signature'], [(m['a'], m['b'], m['t']) for m in s['members']][:12])
        print('   common:', s['common_features'])
    for key, res in searches.items():
        print(key, res['positives'], '/', res['total'], 'terms', res['minimal_terms'], res['solutions'][:5])


if __name__ == '__main__':
    main()
