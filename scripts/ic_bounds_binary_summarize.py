#!/usr/bin/env python3
"""Summaries, tables and figures for item 6 (experiments/ic-bounds-binary/).

  python3 scripts/ic_bounds_binary_summarize.py

Reads results-*.jsonl, related-*.jsonl and validation.json in the experiment directory,
and the existing k <= 7 sweep (experiments/ic-curve-sweep) for the figure only.
Writes summary.json, tables.md, lfd_vs_k.{png,svg}, cost_proxy.{png,svg}.
Every number is a count or an exact certified value; the cost proxies are column/row
counts, not running-time or ECDLP claims.
"""
from __future__ import annotations

import collections
import glob
import hashlib
import json
import math
import os
import statistics

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
EXP = os.path.join(ROOT, 'experiments', 'ic-bounds-binary')


def load(pattern):
    rows = []
    for f in sorted(glob.glob(os.path.join(EXP, pattern))):
        rows += [json.loads(l) for l in open(f) if l.strip()]
    return rows


def squarefree_count(r, D):
    return sum(math.comb(r, i) for i in range(0, min(r, D) + 1)) if D >= 0 else 0


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def hist(values):
    items = sorted(collections.Counter(values).items(), key=lambda kv: (kv[0] is None, kv[0] if kv[0] is not None else 0))
    return {('none' if k is None else str(k)): v for k, v in items}


def group_table(rows, keys):
    groups = collections.defaultdict(list)
    for r in rows:
        groups[tuple(r[k] for k in keys)].append(r)
    out = []
    for key in sorted(groups, key=lambda t: tuple(str(x) if not isinstance(x, int) else f'{x:04d}' for x in t)):
        g = groups[key]
        exact = [r for r in g if r['exact_lfd'] is not None]
        Z = [r['solutions'] for r in g]
        rec = dict(zip(keys, key))
        rec.update(systems=len({r['system_hash'] for r in g}), cases=len(g), exact=len(exact),
                   capped=len(g) - len(exact), lfd_histogram=hist([r['exact_lfd'] for r in g]),
                   max_lfd=max((r['exact_lfd'] for r in exact), default=None),
                   certificate_degree_histogram=hist([r.get('certificate_degree') for r in exact]),
                   unit_ideals=sum(r['status'] == 'exact-unit' for r in g),
                   solutions_min=min(Z), solutions_median=statistics.median(Z), solutions_max=max(Z),
                   variables=sorted({r['variables'] for r in g}),
                   max_columns_at_certificate=max((r['macaulay_at_certificate']['columns'] for r in exact), default=None),
                   max_rows_at_certificate=max((r['macaulay_at_certificate']['rows'] for r in exact), default=None),
                   kernel_r_plus_1=max(r['bounds']['kernel_r_plus_1'] for r in g),
                   written_k_plus_2=(g[0]['k'] + 2) if g[0]['family'] == 'm2' else None,
                   max_closure_seconds=max(r['closure_seconds'] for r in g),
                   max_peak_rss_mb=max(r['peak_rss_mb'] for r in g))
        out.append(rec)
    return out


def fmt_hist(h):
    return ', '.join(f'{k}:{v}' for k, v in h.items())


def md_table(recs, cols, headers):
    lines = ['| ' + ' | '.join(headers) + ' |', '|' + '---|' * len(headers)]
    for r in recs:
        cells = []
        for c in cols:
            v = r.get(c)
            cells.append(fmt_hist(v) if isinstance(v, dict) else ('' if v is None else str(v)))
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def main():
    rows = load('results-*.jsonl')
    related = {r['case_id']: r for r in load('related-*.jsonl')}
    for r in rows:
        rel = related.get(r['case_id'])
        r['d_lin'] = rel['d_lin'] if rel else None
        r['d_full'] = rel['d_full'] if rel else None
        r['static_capped_at'] = rel['capped_at'] if rel else None
    for r in rows:
        rv = r['variables']
        if r['exact_lfd'] is not None:
            # HKYY Prop. 2.8 (as summarised in docs/reviews-curve-lfd-novelty-20261105.md): solving needs
            # V_{max(d_F, e)}, e = number of solutions; in the Boolean ring M_<=D is everything once D >= r.
            r['hkyy_degree'] = max(r['exact_lfd'], r['solutions'])
            r['hkyy_columns'] = squarefree_count(rv, min(rv, r['hkyy_degree']))
            # counting lower bound for the certificate degree: |M_<D| >= |Z| is needed for evaluation onto Z
            r['counting_certificate_lower'] = next(D for D in range(rv + 2) if squarefree_count(rv, D - 1) >= r['solutions'])
            assert r['certificate_degree'] >= r['counting_certificate_lower']
    validation = {}
    for f in sorted(glob.glob(os.path.join(EXP, 'validation-*.json'))):
        v = json.load(open(f))
        validation[os.path.basename(f)] = dict(sha256=sha(f), **{k: x for k, x in v.items() if k not in ('checks', 'mxl_checks')})

    fams = sorted({r['family'] for r in rows})
    summary = dict(rows=len(rows), unique_systems=len({r['system_hash'] for r in rows}),
                   status=dict(collections.Counter(r['status'] for r in rows)),
                   by_family={f: dict(cases=sum(r['family'] == f for r in rows),
                                      exact=sum(r['family'] == f and r['exact_lfd'] is not None for r in rows),
                                      lfd_histogram=hist([r['exact_lfd'] for r in rows if r['family'] == f]))
                              for f in fams},
                   input_files={os.path.basename(f): sha(f) for f in sorted(glob.glob(os.path.join(EXP, 'results-*.jsonl')))
                                + sorted(glob.glob(os.path.join(EXP, 'related-*.jsonl')))})
    summary['by_family_k'] = group_table(rows, ['family', 'k'])
    summary['by_family_n_k'] = group_table(rows, ['family', 'n', 'k'])
    summary['by_family_k_curve_base'] = group_table(rows, ['family', 'k', 'curve_class', 'base_kind'])
    summary['by_family_k_target'] = group_table(rows, ['family', 'k', 'target_class'])
    # bound comparison
    comp = []
    for f in fams:
        for k in sorted({r['k'] for r in rows if r['family'] == f}):
            g = [r for r in rows if r['family'] == f and r['k'] == k and r['exact_lfd'] is not None]
            if not g:
                continue
            mx = max(r['exact_lfd'] for r in g)
            rmax = max(r['variables'] for r in g)
            rec = dict(family=f, k=k, exact_systems=len({r['system_hash'] for r in g}), max_lfd=mx,
                       kernel_r_plus_1=rmax + 1, slack_to_kernel=rmax + 1 - mx)
            if f == 'm2':
                rec.update(kernel_k_l_1=2 * k + 1, written_min_k_l_2=k + 2, slack_to_written=k + 2 - mx)
            comp.append(rec)
    summary['bound_comparison'] = comp
    # cost proxy (exact counts) per family, n, k
    cost = []
    groups = collections.defaultdict(list)
    for r in rows:
        if r['exact_lfd'] is not None:
            groups[(r['family'], r['n'], r['k'])].append(r)
    for key in sorted(groups):
        g = groups[key]
        proper = [r for r in g if r['solutions'] > 0]
        rec = dict(family=key[0], n=key[1], k=key[2], variables=g[0]['variables'], input_degree=max(r['input_degree'] for r in g),
                   systems=len(g), max_lfd=max(r['exact_lfd'] for r in g),
                   max_cost_degree=max(r['cost_degree'] for r in g),
                   max_certificate_degree=max(r['certificate_degree'] for r in g),
                   log2_columns_at_cost_degree=round(math.log2(max(r['macaulay_at_cost_degree']['columns'] for r in g)), 2),
                   log2_columns_at_certificate=round(math.log2(max(r['macaulay_at_certificate']['columns'] for r in g)), 2),
                   log2_rows_at_certificate=round(math.log2(max(r['macaulay_at_certificate']['rows'] for r in g)), 2),
                   solutions_median=statistics.median([r['solutions'] for r in g]),
                   solutions_max=max(r['solutions'] for r in g),
                   proper=len(proper), unit=len(g) - len(proper),
                   d_full_histogram=hist([r['d_full'] for r in g if r['d_full'] is not None or r['static_capped_at']]),
                   d_lin_histogram=hist([r['d_lin'] for r in g if r['d_full'] is not None or r['static_capped_at']]),
                   max_hkyy_degree=max(r['hkyy_degree'] for r in g),
                   hkyy_degree_at_least_r=sum(r['hkyy_degree'] >= r['variables'] for r in g),
                   log2_columns_hkyy_max=round(math.log2(max(r['hkyy_columns'] for r in g)), 2),
                   counting_lower_histogram=hist([r['counting_certificate_lower'] for r in g]),
                   certificate_minus_lfd_histogram=hist([r['certificate_degree'] - r['exact_lfd'] for r in g]),
                   static_scanned=sum(r['d_full'] is not None or r['static_capped_at'] is not None for r in g),
                   static_capped=sum(r['static_capped_at'] is not None and r['d_full'] is None for r in g))
        cost.append(rec)
    summary['cost_proxy'] = cost
    if validation:
        summary['validation'] = validation
    # empirical: LFD against |Z| for m = 2 (proper ideals), per k
    lz = []
    for k in sorted({r['k'] for r in rows if r['family'] == 'm2'}):
        g = [r for r in rows if r['family'] == 'm2' and r['k'] == k and r['exact_lfd'] is not None and r['solutions'] > 0]
        by = collections.defaultdict(list)
        for r in g:
            by[r['exact_lfd']].append(r['solutions'])
        lz.append(dict(k=k, by_lfd={str(v): dict(count=len(z), min_solutions=min(z), max_solutions=max(z))
                                     for v, z in sorted(by.items())}))
    summary['m2_lfd_vs_solutions'] = lz
    # empirical: paired order-two targets (the only target shared across families), same (n, k, curve, basis)
    key = lambda r: (r['n'], r['k'], r['a'], r['b'], tuple(r['basis']))
    two = collections.defaultdict(dict)
    for r in rows:
        if r['target_class'] == 'order-two' and r['exact_lfd'] is not None:
            two[key(r)][r['family']] = r
    pairs = collections.Counter()
    for v in two.values():
        if 'm3-direct' in v and 'm3-split' in v:
            pairs[('direct', v['m3-direct']['exact_lfd'], 'split', v['m3-split']['exact_lfd'])] += 1
    summary['paired_order_two_direct_vs_split'] = {str(k): c for k, c in sorted(pairs.items())}
    json.dump(summary, open(os.path.join(EXP, 'summary.json'), 'w'), indent=1, sort_keys=True)

    lines = ['# Item 6 tables (generated by scripts/ic_bounds_binary_summarize.py)', '']
    lines += ['## Exact LFD by family and k', '',
              md_table(summary['by_family_k'], ['family', 'k', 'variables', 'cases', 'exact', 'capped', 'unit_ideals',
                                                'lfd_histogram', 'max_lfd', 'certificate_degree_histogram',
                                                'written_k_plus_2', 'kernel_r_plus_1'],
                       ['family', 'k', 'r', 'cases', 'exact', 'capped', 'unit', 'LFD histogram', 'max LFD',
                        'D* histogram', 'written k+2', 'kernel r+1']), '']
    lines += ['## Exact LFD by family, n, k', '',
              md_table(summary['by_family_n_k'], ['family', 'n', 'k', 'cases', 'exact', 'lfd_histogram',
                                                  'certificate_degree_histogram', 'solutions_min', 'solutions_median',
                                                  'solutions_max', 'max_columns_at_certificate', 'max_rows_at_certificate',
                                                  'max_closure_seconds', 'max_peak_rss_mb'],
                       ['family', 'n', 'k', 'cases', 'exact', 'LFD histogram', 'D* histogram', '|Z| min', '|Z| median',
                        '|Z| max', 'max cols at D*', 'max rows at D*', 'max s', 'max RSS MB']), '']
    lines += ['## By curve class and base kind', '',
              md_table(summary['by_family_k_curve_base'], ['family', 'k', 'curve_class', 'base_kind', 'cases', 'lfd_histogram', 'max_lfd'],
                       ['family', 'k', 'curve', 'base', 'cases', 'LFD histogram', 'max']), '']
    lines += ['## By target class', '',
              md_table(summary['by_family_k_target'], ['family', 'k', 'target_class', 'cases', 'lfd_histogram', 'solutions_median'],
                       ['family', 'k', 'target', 'cases', 'LFD histogram', '|Z| median']), '']
    lines += ['## Bound comparison (max certified LFD against proved bounds)', '',
              md_table(comp, ['family', 'k', 'exact_systems', 'max_lfd', 'written_min_k_l_2', 'kernel_k_l_1', 'kernel_r_plus_1'],
                       ['family', 'k', 'exact systems', 'max LFD', 'written min(k,l)+2', 'kernel k+l+1', 'kernel r+1']), '']
    lines += ['## Cost-related quantities (exact counts; no running-time or ECDLP claim)', '',
              'cost degree = max(LFD, input degree); columns = squarefree monomials of degree <= D; rows = standard '
              'Boolean Macaulay rows sum_e |M_<=D-deg e|; d_full / d_lin = static Macaulay (no iteration) degrees from '
              'scripts/ic_bounds_binary_related.py.', '',
              md_table(cost, ['family', 'n', 'k', 'variables', 'input_degree', 'systems', 'unit', 'max_lfd', 'max_cost_degree',
                              'max_certificate_degree', 'log2_columns_at_cost_degree', 'log2_columns_at_certificate',
                              'log2_rows_at_certificate', 'solutions_median', 'solutions_max', 'static_scanned',
                              'd_full_histogram', 'd_lin_histogram', 'static_capped', 'max_hkyy_degree',
                              'hkyy_degree_at_least_r', 'log2_columns_hkyy_max', 'certificate_minus_lfd_histogram',
                              'counting_lower_histogram'],
                       ['family', 'n', 'k', 'r', 'd_in', 'systems', 'unit', 'max LFD', 'max cost deg', 'max D*',
                        'log2 cols(cost)', 'log2 cols(D*)', 'log2 rows(D*)', '|Z| med', '|Z| max', 'static scanned',
                        'd_full hist', 'd_lin hist', 'static capped', 'max HKYY deg max(LFD,|Z|)',
                        '#HKYY deg >= r', 'log2 cols(HKYY)', 'D*-LFD hist', 'counting D* lower hist'])]
    lines += ['', '## m = 2: LFD against the number of solutions (proper ideals)', '']
    for rec in lz:
        lines.append(f"- k={rec['k']}: " + '; '.join(f"LFD {v}: {d['count']} systems, |Z| in [{d['min_solutions']}, {d['max_solutions']}]"
                                                    for v, d in rec['by_lfd'].items()))
    lines += ['', '## m = 3: paired order-two targets, direct S_4 against split S_3 (same n, k, curve, basis)', '']
    lines += [f'- {k}: {c}' for k, c in summary['paired_order_two_direct_vs_split'].items()]
    open(os.path.join(EXP, 'tables.md'), 'w').write('\n'.join(lines) + '\n')

    # figures
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    old = []
    sweep = os.path.join(ROOT, 'experiments', 'ic-curve-sweep')
    if os.path.exists(os.path.join(sweep, 'cases.jsonl')):
        res = {}
        for l in open(os.path.join(sweep, 'results.jsonl')):
            x = json.loads(l)
            res[x['system_hash']] = x['exact_lfd']
        for l in open(os.path.join(sweep, 'cases.jsonl')):
            c = json.loads(l)
            if c['kind'] in ('binary', 'koblitz') and res.get(c['system_hash']) is not None:
                old.append((c['k'], res[c['system_hash']]))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ks = list(range(2, 11))
    ax.plot(ks, [k + 2 for k in ks], '--', color='#888', label='written min(k,l)+2 = k+2 (m=2 only)')
    ax.plot(ks, [2 * k + 1 for k in ks], ':', color='#444', label='kernel k+l+1 = 2k+1 (m=2)')
    ax.plot(ks, [3 * k + 1 for k in ks], '-.', color='#c0392b', alpha=0.5, label='kernel r+1 = 3k+1 (m=3 direct)')
    if old:
        oc = collections.Counter(old)
        ax.scatter([k - 0.12 for (k, v) in oc], [v for (k, v) in oc], s=[8 + 2 * c for c in oc.values()],
                   color='#9aa', label='m=2 sweep k<=7 (ic-curve-sweep)')
    m2 = collections.Counter((r['k'], r['exact_lfd']) for r in rows if r['family'] == 'm2' and r['exact_lfd'] is not None)
    ax.scatter([k + 0.12 for (k, v) in m2], [v for (k, v) in m2], s=[8 + 2 * c for c in m2.values()],
               color='#1f5fa8', label='m=2, this item (k=8..10)')
    for f, col, off in (('m3-direct', '#c0392b', 0.0), ('m3-split', '#27864a', 0.24)):
        c3 = collections.Counter((r['k'], r['exact_lfd']) for r in rows if r['family'] == f and r['exact_lfd'] is not None)
        if c3:
            ax.scatter([k + off for (k, v) in c3], [v for (k, v) in c3], s=[8 + 2 * c for c in c3.values()],
                       marker='s', color=col, label=f'm=3 {f.split("-")[1]}, this item')
    ax.set_xlabel('coordinate-space dimension k')
    ax.set_ylabel('exact LFD (certified)')
    ax.set_title('Binary / Koblitz summation systems: exact LFD vs k (marker area = count)')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_ylim(-0.3, 13)
    fig.tight_layout()
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(os.path.join(EXP, f'lfd_vs_k.{ext}'), dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for f, col in (('m2', '#1f5fa8'), ('m3-direct', '#c0392b'), ('m3-split', '#27864a')):
        g = [r for r in rows if r['family'] == f and r['exact_lfd'] is not None and r['solutions'] > 0]
        if g:
            ax.scatter([math.log2(r['solutions']) for r in g],
                       [math.log2(r['macaulay_at_certificate']['columns']) for r in g], s=10, alpha=0.5, color=col, label=f)
    ax.set_xlabel('log2 |Z| (coordinate solutions)')
    ax.set_ylabel('log2 columns |M_<=D*| at the certificate degree')
    ax.set_title('Certificate Macaulay size against the number of solutions (proper ideals)')
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(os.path.join(EXP, f'cost_proxy.{ext}'), dpi=140)
    plt.close(fig)
    print(json.dumps({k: summary[k] for k in ('rows', 'unique_systems', 'status', 'by_family')}, indent=1))


if __name__ == '__main__':
    main()
