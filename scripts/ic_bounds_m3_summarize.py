#!/usr/bin/env python3
"""Summarize experiments/ic-bounds-m3/results*.jsonl into summary.json, table.md and figures.

Also re-checks, per direct case, the hypotheses of the carry lower bound (README, Proposition 1):
deg NF(S_4) = 6 iff 5t^3 + 3at + 6b != 0 (mod p), and NF(u_{1,0} * g) is not a scalar multiple of
g (instance check of the support condition), and that every exact value satisfies the
kernel-checked upper bound LFD <= mk + 1 = 3k + 1.
"""
import hashlib
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_m3 as M  # noqa: E402

OUT = ROOT / 'experiments' / 'ic-bounds-m3'


def shift_block1(g, k):
    N = 1 << k
    out = {}
    for ms, c in g.items():
        e1 = ((ms & (N - 1)) + 1) % N
        out[(ms & ~(N - 1)) | e1] = c
    return out


def carry_check(case):
    p, k = case['prime'], case['k']
    E = I.PrimeCurve(p, case['a'], case['b'])
    t = case['target'][0]
    g = M.direct_generator(E, t, k)
    deg = max((ms.bit_count() for ms in g), default=None)
    c333 = (-8 * (5 * t ** 3 + 3 * E.a * t + 6 * E.b)) % p
    top = (3 | (3 << k) | (3 << (2 * k))) if k >= 2 else None
    if k >= 2:
        assert g.get(top, 0) == c333, 'coefficient of the (3,3,3) mask differs from -8(5t^3+3at+6b)'
        assert (deg == 6) == (c333 != 0)
    h = shift_block1(g, k) if g else {}
    hdeg = max((ms.bit_count() for ms in h), default=None)
    proportional = False
    if g:
        ms0, c0 = next(iter(g.items()))
        lam = h.get(ms0, 0) * pow(c0, -1, p) % p
        proportional = set(h) == set(g) and all(h[ms] == lam * g[ms] % p for ms in g)
    hyp = k >= 2 and c333 != 0 and hdeg is not None and hdeg <= 6 and not proportional
    return dict(generator_degree=deg, c333=c333, shifted_degree=hdeg, shift_proportional=proportional,
                carry_lower_bound_7_applies=hyp)


def main():
    cases = {json.loads(l)['case_id']: json.loads(l) for l in open(OUT / 'cases.jsonl')}
    res = defaultdict(dict)
    for path in sorted(OUT.glob('results*.jsonl')):
        for l in open(path):
            r = json.loads(l)
            res[r['case_id']][r['encoding']] = r
    rows = []
    for cid, c in cases.items():
        if cid not in res:
            continue
        row = dict(case_id=cid, prime=c['prime'], a=c['a'], b=c['b'], curve_class=c['curve_class'], k=c['k'],
                   target=c['target'], target_class=c['target_class'])
        d = res[cid].get('direct')
        if d:
            v = d['results']['direct']
            row.update(direct_status=v['status'], direct_lfd=v['exact_lfd'], direct_lower=v['lfd_lower_bound'],
                       solutions=v['quotient_dimension'], generator_degree=d['generator_degrees'][0],
                       direct_falls=d['strict_falls'],
                       direct_dims=[(x['D'], x['monomials_le_D'], x['quotient_rank'], x['dim_V_D']) for x in d['profile']],
                       rational_ordered=d['geometry']['rational_ordered_triples'],
                       rational_unordered=d['geometry']['rational_unordered'],
                       direct_hash=d['generator_hash'], direct_seconds=d['seconds'])
            if v['exact_lfd'] is not None:
                assert v['exact_lfd'] <= 3 * c['k'] + 1, 'kernel upper bound violated'
            row.update(carry_check(c))
            if row['carry_lower_bound_7_applies'] and v['exact_lfd'] is not None:
                assert v['exact_lfd'] >= 7, 'carry lower bound contradicted'
        s = res[cid].get('split')
        if s:
            f, fe = s['results']['free'], s['results']['field']
            row.update(split_free_status=f['status'], split_free_lfd=f['exact_lfd'],
                       split_free_lower=f.get('lfd_lower_bound'), split_free_dim=f['quotient_dimension'],
                       split_field_dim=fe['quotient_dimension'], split_falls=s['strict_falls'],
                       split_degrees=s['generator_degrees'], split_hash=s['generator_hash'],
                       split_dims=[(x['D'], x['monomials_le_D'], x['quotient_rank'], x['dim_V_D']) for x in s['profile']],
                       split_seconds=s['seconds'])
            if fe['exact_lfd'] is not None:
                row['split_field_status'], row['split_field_lfd'] = fe['status'], fe['exact_lfd']
            elif f['quotient_dimension'] is not None and fe['quotient_dimension'] is not None \
                    and fe['quotient_dimension'] < f['quotient_dimension']:
                # I_free strictly inside I_field: some h in I_field of degree < p is missing from every
                # V_D with D < p (those equal the I_free closures), so LFD(F_field) >= p.
                row['split_field_status'], row['split_field_lfd'] = 'lower-bound-p (proof)', None
                row['split_field_lower'] = c['prime']
            else:
                row['split_field_status'], row['split_field_lfd'] = fe['status'], None
        m2 = res[cid].get('m2')
        if m2:
            v = m2['results']['m2']
            row.update(m2_status=v['status'], m2_lfd=v['exact_lfd'], m2_solutions=v['quotient_dimension'],
                       m2_falls=m2['strict_falls'])
        rows.append(row)
    rows.sort(key=lambda r: (r['k'], r['prime'], r['curve_class'], r['a'], r['b'], r['target_class'], r['target']))

    def hist(vals):
        return dict(sorted(Counter(str(v) for v in vals).items()))

    groups = []
    for k in sorted({r['k'] for r in rows}):
        sub = [r for r in rows if r['k'] == k]
        nondeg = [r for r in sub if r.get('generator_degree') == 6]
        g = dict(k=k, cases=len(sub),
                 direct=dict(n=sum('direct_status' in r for r in sub),
                             exact=sum(r.get('direct_lfd') is not None for r in sub),
                             hist=hist(r['direct_lfd'] for r in sub if 'direct_status' in r),
                             hist_deg6=hist(r['direct_lfd'] for r in nondeg),
                             hist_lfd_minus_deg=hist((r['direct_lfd'] - r['generator_degree'])
                                                     for r in sub if r.get('direct_lfd') is not None
                                                     and r.get('generator_degree')),
                             unit=sum(r.get('direct_status') == 'exact-unit' for r in sub),
                             kernel_upper=3 * k + 1),
                 split_free=dict(n=sum('split_free_status' in r for r in sub),
                                 exact=sum(r.get('split_free_lfd') is not None for r in sub),
                                 hist=hist(r.get('split_free_lfd') for r in sub if 'split_free_status' in r),
                                 capped_lower=hist(r.get('split_free_lower') for r in sub
                                                   if r.get('split_free_status') == 'capped')),
                 split_field=dict(hist=hist(r.get('split_field_status') + ':' + str(r.get('split_field_lfd'))
                                            for r in sub if 'split_field_status' in r)),
                 m2=dict(n=sum('m2_status' in r for r in sub),
                         hist=hist(r.get('m2_lfd') for r in sub if 'm2_status' in r)),
                 paired_direct_minus_m2=hist(r['direct_lfd'] - r['m2_lfd'] for r in sub
                                             if r.get('direct_lfd') is not None and r.get('m2_lfd') is not None),
                 paired_split_minus_m2=hist(r['split_free_lfd'] - r['m2_lfd'] for r in sub
                                            if r.get('split_free_lfd') is not None and r.get('m2_lfd') is not None),
                 carry_lower_bound_applies=sum(bool(r.get('carry_lower_bound_7_applies')) for r in sub))
        groups.append(g)
    by_class = defaultdict(Counter)
    for r in rows:
        if 'direct_status' in r:
            by_class[(r['k'], r['curve_class'], r['target_class'], 'direct')][str(r['direct_lfd'])] += 1
        if 'split_free_status' in r:
            by_class[(r['k'], r['curve_class'], r['target_class'], 'split-free')][str(r['split_free_lfd'])] += 1
        if 'm2_status' in r:
            by_class[(r['k'], r['curve_class'], r['target_class'], 'm2')][str(r['m2_lfd'])] += 1
    # dimensions of V_D (direct, generator degree 6, proper ideals): per k and D, min/median/max
    dims = defaultdict(list)
    for r in rows:
        for D, n, q, o in r.get('direct_dims', []):
            dims[('direct', r['k'], D, n)].append(q)
        for D, n, q, o in r.get('split_dims', []):
            dims[('split', r['k'], D, n)].append(q)
    dim_table = [dict(encoding=e, k=k, D=D, monomials_le_D=n, cases=len(v), quotient_rank_min=min(v),
                      quotient_rank_median=statistics.median(v), quotient_rank_max=max(v))
                 for (e, k, D, n), v in sorted(dims.items())]
    summary = dict(seed=20261107, cases_with_results=len(rows), groups=groups,
                   by_class={' | '.join(map(str, k)): dict(v) for k, v in sorted(by_class.items())},
                   dimension_table=dim_table,
                   files={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(OUT.glob('*.jsonl'))})
    summary['files'] = {name: h for name, h in summary['files'].items() if name != 'rows.jsonl'}
    summary['rows_sha256'] = hashlib.sha256(''.join(json.dumps(r, sort_keys=True) + '\n'
                                                    for r in rows).encode()).hexdigest()
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1))
    with open(OUT / 'rows.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(r, sort_keys=True) + '\n')
    # markdown table
    lines = ['| k | direct m=3 LFD (all) | direct, deg NF = 6 | LFD - deg NF | split I_free LFD | split I_field | m=2 LFD | direct - m2 | split - m2 |',
             '|---|---|---|---|---|---|---|---|---|']
    for g in groups:
        lines.append(f"| {g['k']} | {g['direct']['hist']} | {g['direct']['hist_deg6']} | {g['direct']['hist_lfd_minus_deg']} "
                     f"| {g['split_free']['hist']} | {g['split_field']['hist']} | {g['m2']['hist']} "
                     f"| {g['paired_direct_minus_m2']} | {g['paired_split_minus_m2']} |")
    lines += ['', '| encoding | k | D | |M_<=D| (quotient) | cases | dim V_D/chains min | median | max |', '|---|---|---|---|---|---|---|---|']
    for d in dim_table:
        lines.append(f"| {d['encoding']} | {d['k']} | {d['D']} | {d['monomials_le_D']} | {d['cases']} | "
                     f"{d['quotient_rank_min']} | {d['quotient_rank_median']} | {d['quotient_rank_max']} |")
    (OUT / 'table.md').write_text('\n'.join(lines) + '\n')
    # figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        series = [('m2_lfd', 'm=2 (S_3)', 'o', -0.15), ('split_free_lfd', 'm=3 split, z free', 's', 0.0),
                  ('direct_lfd', 'm=3 direct NF(S_4)', '^', 0.15)]
        for key, label, mk, dx in series:
            pts = Counter((r['k'], r[key]) for r in rows if r.get(key) is not None)
            if not pts:
                continue
            xs = [kk + dx for kk, _ in pts]
            ys = [v for _, v in pts]
            ss = [12 + 4 * n for n in pts.values()]
            ax.scatter(xs, ys, s=ss, marker=mk, label=label, alpha=0.7)
        ks = sorted({r['k'] for r in rows})
        ax.plot(ks, [3 * k + 1 for k in ks], 'k--', lw=1, label='kernel upper bound 3k+1 (direct)')
        ax.plot(ks, [2 * k + 1 for k in ks], 'k:', lw=1, label='kernel upper bound 2k+1 (m=2)')
        ax.set_xlabel('chain depth k (|H| = 2^k)')
        ax.set_ylabel('exact certified LFD (marker area ~ cases)')
        ax.set_xticks(ks)
        ax.legend(fontsize=7)
        ax.set_title('Exact LFD, prime subgroup-chain systems, m=3 vs m=2')
        fig.tight_layout()
        fig.savefig(OUT / 'lfd_vs_k.png', dpi=130)
        fig.savefig(OUT / 'lfd_vs_k.svg')
    except ImportError:
        pass
    print(json.dumps(groups, indent=1))


if __name__ == '__main__':
    main()
