"""Separate the subfield and solution-count confounders in the binary/Koblitz base-class comparison.

Reads experiments/ic-curve-sweep/{cases,results}.jsonl and writes
experiments/ic-curve-sweep/base-class-analysis.json: unit-ideal rates and mean
solution counts per class, paired-cell comparisons for all systems and for
proper ideals only, and a comparison stratified by depth k and solution-count
bucket (unit, 1-2, 3-8, 9-32, 33+).
"""
from collections import defaultdict
from pathlib import Path
import json


def cls(c):
    if c['base_kind'] == 'random-subspace':
        return 'random'
    return 'subfield' if c['subfield'] else 'invariant-nonsubfield'


def bucket(z):
    return 'unit' if z == 0 else '1-2' if z <= 2 else '3-8' if z <= 8 else '9-32' if z <= 32 else '33+'


def paired(sub, res, a_cl, b_cl, proper_only=False):
    A = defaultdict(list)
    B = defaultdict(list)
    for c in sub:
        r = res[c['system_hash']]
        if proper_only and r['status'] == 'exact-unit':
            continue
        key = (c['n'], c['a'], c['b'], c['k'], c['target_class'])
        if cls(c) == a_cl:
            A[key].append(r['exact_lfd'])
        if cls(c) == b_cl:
            B[key].append(r['exact_lfd'])
    common = sorted(set(A) & set(B))
    d = [max(A[k]) - max(B[k]) for k in common]
    return dict(cells=len(common), a_higher=sum(x > 0 for x in d), a_lower=sum(x < 0 for x in d), equal=sum(x == 0 for x in d))


def main():
    root = Path(__file__).resolve().parents[1]
    sweep = root / 'experiments/ic-curve-sweep'
    cases = [json.loads(l) for l in (sweep / 'cases.jsonl').read_text().splitlines()]
    res = {r['system_hash']: r for r in map(json.loads, (sweep / 'results.jsonl').read_text().splitlines())}
    report = {}
    for kind in ('binary', 'koblitz'):
        sub = [c for c in cases if c['kind'] == kind]
        classes = {}
        for cl in ('random', 'invariant-nonsubfield', 'subfield'):
            seen = set()
            units = tot = 0
            zs = []
            for c in sub:
                h = c['system_hash']
                if cls(c) != cl or h in seen:
                    continue
                seen.add(h)
                tot += 1
                units += res[h]['status'] == 'exact-unit'
                zs.append(c['polynomial_tuples'])
            classes[cl] = dict(systems=tot, unit_rate=round(units / tot, 3) if tot else None,
                               mean_solutions=round(sum(zs) / len(zs), 2) if zs else None)
        strata = defaultdict(lambda: defaultdict(list))
        seen = set()
        for c in sub:
            h = c['system_hash']
            if h in seen:
                continue
            seen.add(h)
            strata[(c['k'], bucket(c['polynomial_tuples']))][cls(c)].append(res[h]['exact_lfd'])
        stratified = []
        tot = wdiff = 0
        higher = lower = equal = 0
        for key in sorted(strata):
            d = strata[key]
            if d['random'] and d['invariant-nonsubfield']:
                mr = sum(d['random']) / len(d['random'])
                mi = sum(d['invariant-nonsubfield']) / len(d['invariant-nonsubfield'])
                w = min(len(d['random']), len(d['invariant-nonsubfield']))
                tot += w
                wdiff += w * (mi - mr)
                higher += mi > mr + 1e-9
                lower += mi < mr - 1e-9
                equal += abs(mi - mr) <= 1e-9
                stratified.append(dict(k=key[0], solutions=key[1], random_n=len(d['random']), random_mean=round(mr, 3),
                                       invariant_n=len(d['invariant-nonsubfield']), invariant_mean=round(mi, 3)))
        report[kind] = dict(classes=classes,
                            paired_all={f'{a} vs {b}': paired(sub, res, a, b) for a, b in
                                        [('random', 'invariant-nonsubfield'), ('random', 'subfield'), ('invariant-nonsubfield', 'subfield')]},
                            paired_proper_only={f'{a} vs {b}': paired(sub, res, a, b, True) for a, b in
                                                [('random', 'invariant-nonsubfield'), ('random', 'subfield')]},
                            stratified_by_k_and_solutions=stratified,
                            stratified_summary=dict(strata=higher + lower + equal, invariant_higher=higher, invariant_lower=lower,
                                                    equal=equal, weighted_mean_invariant_minus_random=round(wdiff / tot, 3) if tot else None))
    (sweep / 'base-class-analysis.json').write_text(json.dumps(report, indent=1) + '\n')
    for kind in report:
        print(kind, json.dumps(report[kind]['classes']), json.dumps(report[kind]['paired_proper_only']),
              json.dumps(report[kind]['stratified_summary']))


if __name__ == '__main__':
    main()
