#!/usr/bin/env python3
"""Analyze experiments/ic-bounds-m3-depth/closures.jsonl: exact LFD per (instance, k), falling-behind
depths per degree, and the three-block free-extension recursion p_j(k+1) = sum_w C(3,w) p_{j-w}(k).

Exact LFD rule (README of experiments/ic-bounds-m3, local-to-global criterion):
  * V_D = 0 modulo the chains for D < deg NF(S_4) = 6 (no input of degree <= D); asserted
    generator_degree == 6 for every row used;
  * strict fall at D: dim(V_D cap M_<=D-1) > dim V_{D-1};
  * unit: least D with 1 in V_D;
  * proper: at the first D with eval rank(M_<D on Z) = |Z| and dim V_D = |M_<=D| - |Z|
    (equivalently V_D = I cap M_<=D), exact LFD = last strict fall <= D.
"""
from __future__ import annotations

import json
import os
from math import comb

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'ic-bounds-m3-depth')


def load():
    rows = {}
    for fn in sorted(os.listdir(HERE)):
        if fn.startswith('closures') and fn.endswith('.jsonl'):
            with open(os.path.join(HERE, fn)) as f:
                for line in f:
                    if line.strip():
                        r = json.loads(line)
                        rows[(r['instance'], r['k'], r['D'])] = r
    return rows


def rank_of(rows, inst, k, D):
    """dim V_D modulo chains; 0 for D < generator degree."""
    r = rows.get((inst, k, D))
    if r is not None:
        return r['rank']
    any_row = next((v for (i, kk, _), v in rows.items() if i == inst and kk == k), None)
    if any_row is not None and D < any_row['generator_degree']:
        return 0
    return None


def exact_lfd(rows, inst, k):
    Ds = sorted(D for (i, kk, D) in rows if i == inst and kk == k)
    falls = []
    for D in range(1, max(Ds) + 1):
        r = rows.get((inst, k, D))
        prev = rank_of(rows, inst, k, D - 1)
        if r is None:
            if D < 6:
                continue
            return dict(status='capped', lower_bound=max(falls, default=0), falls=falls, missing_D=D)
        assert r['generator_degree'] == 6
        if prev is None:
            return dict(status='capped', lower_bound=max(falls, default=0), falls=falls, missing_D=D - 1)
        if r['profile'][D - 1] > prev:
            falls.append(D)
        if r['contains_one']:
            assert r['solutions'] == 0
            return dict(status='exact-unit', lfd=D, falls=falls, certificate=dict(kind='unit', degree=D))
        if (r['solutions'] > 0 and r['eval_rank_le_j'][D - 1] == r['solutions']
                and r['rank'] == r['monomials_le_j'][D] - r['solutions']):
            return dict(status='exact-proper-local', lfd=max(falls, default=0), falls=falls,
                        certificate=dict(kind='local-to-global', degree=D, solutions=r['solutions'],
                                         closure_check=r['closure_check']))
    lb = max(falls, default=0)
    top = rows.get((inst, k, max(Ds)))
    if top is not None and top['closure_check'] == 'closed-sweep' and not top['equals_ideal']:
        # V_D computed exactly and strictly inside I cap M_<=D: some ideal element of degree <= D
        # is missing from V_D, so LFD >= D + 1 (exact lower bound).
        return dict(status='lower-bound-deficient-closure', lower_bound=max(lb, top['D'] + 1), falls=falls,
                    certificate=dict(kind='closed V_D strictly inside I cap M_<=D', degree=top['D'],
                                     dim_V_D=top['rank'], dim_I_le_D=top['ideal_dim_le_j'][top['D']]))
    return dict(status='capped', lower_bound=lb, falls=falls)


def recursion_pred(prof, D):
    return [sum(comb(3, w) * prof[j - w] for w in range(4) if j - w >= 0) for j in range(D + 1)]


def main():
    rows = load()
    insts = sorted({i for (i, _, _) in rows})
    ks = sorted({k for (_, k, _) in rows})
    Dset = sorted({D for (_, _, D) in rows})
    out = dict(lfd={}, behind={}, recursion=[], ideal_recursion=[])
    for inst in insts:
        for k in ks:
            if any(i == inst and kk == k for (i, kk, _) in rows):
                out['lfd'][f'{inst}|{k}'] = exact_lfd(rows, inst, k)
        for D in Dset:
            first = None
            seen = []
            for k in ks:
                r = rows.get((inst, k, D))
                if r is None:
                    continue
                eq = r['profile'][D] == r['ideal_dim_le_j'][D]
                seen.append([k, eq])
                if not eq and first is None:
                    first = k
            out['behind'][f'{inst}|{D}'] = dict(first_depth_behind=first, depths_equal=seen)
            for k in ks:
                a, b = rows.get((inst, k, D)), rows.get((inst, k + 1, D))
                if a is None or b is None:
                    continue
                pred = recursion_pred(a['profile'], D)
                bad = [j for j in range(D + 1) if pred[j] != b['profile'][j]]
                full = lambda r: r['contains_one']
                gonly = lambda r: r['profile'] == [0] * D + [1]
                trivial = (full(a) and full(b)) or (gonly(a) and gonly(b))
                out['recursion'].append(dict(instance=inst, D=D, k=k, k1=k + 1, holds=not bad, trivial=trivial,
                                             first_bad_j=bad[0] if bad else None,
                                             predicted=pred, actual=b['profile']))
                ip = recursion_pred(a['ideal_dim_le_j'], D)
                out['ideal_recursion'].append(dict(instance=inst, D=D, k=k, holds=ip == b['ideal_dim_le_j'],
                                                   predicted=ip, actual=b['ideal_dim_le_j']))
    with open(os.path.join(HERE, 'summary.json'), 'w') as f:
        json.dump(out, f, indent=1)

    L = ['# Closure profiles (generated by scripts/ic_bounds_m3_depth_analyze.py)', '',
         'p_j = dim(V_D cap M_<=j) modulo the chains; ideal_j = dim(I cap M_<=j) = |M_<=j| - rank ev(M_<=j -> F_p^Z).', '',
         '| instance | k | D | |M_<=D| | |Z| | profile p_0..p_D | ideal dim (I cap M_<=D) | V_D = I cap M_<=D | closure check | s |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for (inst, k, D) in sorted(rows):
        r = rows[(inst, k, D)]
        L.append(f"| {inst} | {k} | {D} | {r['monomials_le_j'][D]} | {r['solutions']} | {' '.join(map(str, r['profile']))} | "
                 f"{r['ideal_dim_le_j'][D]} | {'yes' if r['equals_ideal'] else 'no'} | {r['closure_check']} | {r['seconds']} |")
    L += ['', '## Exact LFD', '', '| instance | k | status | LFD | strict falls | certificate |', '|---|---|---|---|---|---|']
    for key, v in out['lfd'].items():
        inst, k = key.split('|')
        L.append(f"| {inst} | {k} | {v['status']} | {v.get('lfd', '>= %s' % v.get('lower_bound'))} | {v['falls']} | "
                 f"{json.dumps(v.get('certificate'))} |")
    L += ['', '## Recursion p_j(k+1) = sum_w C(3,w) p_{j-w}(k)', '',
          'trivial = both closures are the whole space M_<=D (unit ideal; |M_<=j| obeys the recursion by Vandermonde) or both are span(g) (D = 6).', '',
          '| instance | D | k -> k+1 | holds | trivial | first failing j | predicted | actual |', '|---|---|---|---|---|---|---|---|']
    for r in out['recursion']:
        L.append(f"| {r['instance']} | {r['D']} | {r['k']} -> {r['k1']} | {r['holds']} | {r['trivial']} | {r['first_bad_j']} | "
                 f"{' '.join(map(str, r['predicted']))} | {' '.join(map(str, r['actual']))} |")
    L += ['', '## First depth where V_D != I cap M_<=D', '', '| instance | D | first depth behind | (k, equal) |', '|---|---|---|---|']
    for key, v in out['behind'].items():
        inst, D = key.split('|')
        L.append(f"| {inst} | {D} | {v['first_depth_behind']} | {v['depths_equal']} |")
    nt = os.path.join(HERE, 'nosandwich-k5.jsonl')
    if os.path.exists(nt):
        L += ['', '## Full closedness sweep (no dimension sandwich) vs sandwich runs', '',
              '| instance | k | D | profile (full sweep) | same as sandwich run | check |', '|---|---|---|---|---|---|']
        for line in open(nt):
            r = json.loads(line)
            q = rows.get((r['instance'], r['k'], r['D']))
            L.append(f"| {r['instance']} | {r['k']} | {r['D']} | {' '.join(map(str, r['profile']))} | "
                     f"{q is not None and q['profile'] == r['profile']} | {r['closure_check']} |")
    with open(os.path.join(HERE, 'table.md'), 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L[-60:]))


if __name__ == '__main__':
    main()
