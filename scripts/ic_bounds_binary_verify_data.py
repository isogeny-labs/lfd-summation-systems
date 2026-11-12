#!/usr/bin/env python3
"""Recompute the item-6 headline counts from the stored rows (verification pass).

Reads experiments/ic-bounds-binary/results-*.jsonl and related-*.jsonl only; does not
import the producer's summarizer.  Prints a JSON object of recomputed quantities.

  python3 scripts/ic_bounds_binary_verify_data.py > experiments/ic-bounds-binary/verify/data-checks.json
"""
import collections
import hashlib
import json
import math
import os
import statistics

EXP = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'ic-bounds-binary')


def load(name):
    return [json.loads(l) for l in open(os.path.join(EXP, name))]


def sha(name):
    return hashlib.sha256(open(os.path.join(EXP, name), 'rb').read()).hexdigest()


def M(r, D):
    return sum(math.comb(r, i) for i in range(0, min(r, D) + 1))


def main():
    fams = {f: load(f'results-{f}.jsonl') for f in ('m2', 'm3-direct', 'm3-split')}
    allrows = [r for f in fams.values() for r in f]
    out = {}
    out['file_sha256'] = {n: sha(n) for n in ['results-m2.jsonl', 'results-m3-direct.jsonl', 'results-m3-split.jsonl',
                                               'related-m2.jsonl', 'related-m3-direct.jsonl', 'related-m3-split.jsonl']}
    out['cases'] = len(allrows)
    out['cases_by_family'] = {f: len(v) for f, v in fams.items()}
    out['distinct_system_hashes'] = len({r['system_hash'] for r in allrows})
    out['status'] = dict(collections.Counter(r['status'] for r in allrows))
    out['exact_lfd_none'] = sum(r['exact_lfd'] is None for r in allrows)
    # duplicates: are they exactly the (0,1)/(1,1) pairs?
    byh = collections.defaultdict(list)
    for r in allrows:
        byh[r['system_hash']].append(r)
    dup_groups = [g for g in byh.values() if len(g) > 1]
    out['duplicate_groups'] = len(dup_groups)
    out['duplicate_groups_all_koblitz_pairs'] = all(
        len(g) == 2 and {(x['a'], x['b']) for x in g} == {(0, 1), (1, 1)} for g in dup_groups)
    out['duplicate_groups_lfd_consistent'] = all(len({x['exact_lfd'] for x in g}) == 1 for g in dup_groups)
    out['duplicate_target_classes'] = dict(collections.Counter(g[0]['target_class'] for g in dup_groups))
    # per family/k
    per = {}
    for f, rows in fams.items():
        for k in sorted({r['k'] for r in rows}):
            rr = [r for r in rows if r['k'] == k]
            hist = collections.Counter(r['exact_lfd'] for r in rr)
            per[f'{f}:k={k}'] = dict(cases=len(rr), systems=len({r['system_hash'] for r in rr}),
                                    max_lfd=max(hist), hist=dict(sorted(hist.items())),
                                    unit=sum(r['status'] == 'exact-unit' for r in rr),
                                    r=sorted({r['variables'] for r in rr}),
                                    n=sorted({r['n'] for r in rr}))
    out['per_family_k'] = per
    # proper / unit counts over cases and systems
    out['proper_cases'] = sum(r['status'] == 'exact-proper-local' for r in allrows)
    out['unit_cases'] = sum(r['status'] == 'exact-unit' for r in allrows)
    # LFD-5 binary m2
    m2 = fams['m2']
    l5 = [r for r in m2 if r['exact_lfd'] == 5]
    out['m2_lfd5_cases'] = len(l5)
    out['m2_lfd5_systems'] = len({r['system_hash'] for r in l5})
    out['m2_lfd5_systems_by_n'] = dict(collections.Counter(
        n for n, _ in {(r['n'], r['system_hash']) for r in l5}))
    k10n11 = [r for r in m2 if r['k'] == 10 and r['n'] == 11]
    out['m2_k10_n11_proper_systems'] = len({r['system_hash'] for r in k10n11 if r['status'] == 'exact-proper-local'})
    out['m2_k10_n11_cases'] = len(k10n11)
    # separation at k=10 (systems)
    sep = {}
    for k in (8, 9, 10):
        prop = {r['system_hash']: r for r in m2 if r['k'] == k and r['status'] == 'exact-proper-local'}
        by = collections.defaultdict(list)
        for r in prop.values():
            by[r['exact_lfd']].append(r['solutions'])
        sep[k] = {l: dict(systems=len(v), zmin=min(v), zmax=max(v)) for l, v in sorted(by.items())}
    out['m2_lfd_vs_Z_systems'] = sep
    # |Z| vs 2^(2k-n), random bases, non-order-two
    tr = {}
    ratio = []
    for n in (11, 12, 13, 14):
        for k in (8, 9, 10):
            rr = [r for r in m2 if r['n'] == n and r['k'] == k and r['base_kind'] == 'random-subspace'
                  and r['target_class'] != 'order-two']
            if rr:
                tr[f'n={n},k={k}'] = dict(cases=len(rr), mean_Z=round(statistics.mean(r['solutions'] for r in rr), 1),
                                         pred=2 ** (2 * k - n))
    for k in (8, 9, 10):
        rr = [r for r in m2 if r['k'] == k and r['solutions'] > 0]
        tr[f'unordered_over_Z_k={k}'] = round(sum(r['rational_unordered_multisets'] for r in rr) /
                                              sum(r['solutions'] for r in rr), 3)
    out['m2_Z_vs_2^(2k-n)'] = tr
    # D* - LFD, counting bound, HKYY
    dd = collections.Counter(r['certificate_degree'] - r['exact_lfd'] for r in allrows)
    out['dstar_minus_lfd'] = dict(dd)
    out['max_dstar'] = {f: max(r['certificate_degree'] for r in rows) for f, rows in fams.items()}
    cnt_eq = collections.Counter()
    cnt_viol = 0
    hk = collections.Counter()
    hk_tot = collections.Counter()
    macmax = {}
    for f, rows in fams.items():
        for r in rows:
            R = r['variables']
            if r['status'] == 'exact-proper-local':
                Z = r['solutions']
                lb = min(D for D in range(0, R + 2) if M(R, D - 1) >= Z)
                if lb > r['certificate_degree']:
                    cnt_viol += 1
                if lb == r['certificate_degree']:
                    cnt_eq[f] += 1
                hk_tot[f] += 1
                if max(r['exact_lfd'], Z) >= R:
                    hk[f] += 1
            mc = r['macaulay_at_certificate']
            assert mc['columns'] == M(R, mc['degree'])
            macmax[f] = max(macmax.get(f, (0, 0)), (mc['columns'], mc['rows']))
    out['counting_bound_equals_dstar_proper'] = dict(cnt_eq, total=sum(cnt_eq.values()))
    out['counting_bound_violations'] = cnt_viol
    out['hkyy_ge_r'] = dict(hk, total=sum(hk.values()))
    out['hkyy_proper_totals'] = dict(hk_tot, total=sum(hk_tot.values()))
    n11 = [r for r in m2 if r['n'] == 11 and r['status'] == 'exact-proper-local']
    out['hkyy_every_m2_n11_proper'] = all(max(r['exact_lfd'], r['solutions']) >= r['variables'] for r in n11)
    out['largest_certificate_macaulay'] = {f: dict(columns=v[0], rows=v[1]) for f, v in macmax.items()}
    out['max_rows_m2_k10_D5'] = max(r['macaulay_at_certificate']['rows'] for r in m2
                                    if r['k'] == 10 and r['certificate_degree'] == 5)
    # Macaulay row formula check on m2 rows (all equations of degree 2 in m2)
    # paired order-two direct vs split
    def key(r):
        return (r['n'], r['k'], r['a'], r['b'], tuple(r['basis']))
    d2 = {key(r): r for r in fams['m3-direct'] if r['target_class'] == 'order-two'}
    s2 = {key(r): r for r in fams['m3-split'] if r['target_class'] == 'order-two'}
    common = sorted(set(d2) & set(s2))
    cmp = collections.Counter(
        'direct>split' if d2[c]['exact_lfd'] > s2[c]['exact_lfd'] else
        'equal' if d2[c]['exact_lfd'] == s2[c]['exact_lfd'] else 'direct<split' for c in common)
    out['paired_order_two'] = dict(pairs=len(common), **cmp)
    # related: d_full
    rel = {}
    for f in ('m2', 'm3-direct', 'm3-split'):
        rr = load(f'related-{f}.jsonl')
        keys = sorted(rr[0].keys())
        res = fams[f]
        byid = {r['case_id']: r for r in res}
        met = [x for x in rr if x.get('d_full') is not None]
        viol = [x for x in met if byid[x['case_id']]['certificate_degree'] > x['d_full']]
        rel[f] = dict(rows=len(rr), keys=keys, d_full_met=len(met), not_met=len(rr) - len(met), dstar_gt_dfull=len(viol))
        if f == 'm3-direct':
            for k in (2, 3, 4, 5):
                ks = [x['d_full'] for x in met if byid[x['case_id']]['k'] == k]
                if ks:
                    rel[f][f'd_full_range_k={k}'] = [min(ks), max(ks)]
    out['related'] = rel
    print(json.dumps(out, indent=1, sort_keys=True, default=str))


if __name__ == '__main__':
    main()
