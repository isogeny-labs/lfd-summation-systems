"""Exact LFD through the filtered-generation criterion (docs/proofs/U1-upper-bounds.md, Prop. 4):

    LFD <= D   iff   V_D = I_D   and   I_{D'+1} = I_{D'} + sum_v v I_{D'}  for all D' >= D,

where I_D = I ∩ M_<=D is the kernel of evaluation on Z.  By Theorem 1 the generation condition
holds for every D' >= k+2, so only D' <= k+1 is computed.  For each instance the script records

  dimV[D], dimI[D]       for D = 2 .. first D with V_D = I_D (closures verified closed),
  gen[D']                whether I_{D'+1} = I_{D'} + sum_v v I_{D'}, D' = 2 .. k+1,
  gen_degree             least D0 with gen[D'] for all D' in [D0, k+1],
  lfd                    least D >= gen_degree with V_D = I_D.

Modes:
  small    exhaustive planes at k = 1, 2, 3 (small p) and random instances at k = 4, 5
  sweep    the subgroup two-summand prime cases of experiments/ic-curve-sweep (k <= 5 all,
           k = 6 a deterministic sample), compared with the recorded exact LFD
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import random
import sys
import time

import numpy as np

from chain_algebra import (Window, Echelon, closure, generator_row, kernel_basis, profile,
                           solutions)

HERE = os.path.dirname(os.path.abspath(__file__))


def embed(rows, src, dst):
    out = np.zeros((len(rows), dst.n), dtype=np.int64)
    cols = np.array([dst.index[m] for m in src.mons], dtype=np.int64)
    out[:, cols] = rows
    return out


def analyse(p, a, b, t, k, with_generation=True, dmax=None):
    Z = solutions(p, a, b, t, k)
    dmax = (k + 2) if dmax is None else dmax
    wins = {}

    def win(D):
        if D not in wins:
            wins[D] = Window(k, k, D)
        return wins[D]

    out = dict(p=p, a=a, b=b, t=t, k=k, points=len(Z), dimV={}, dimI={}, profile={}, gen={})
    # evaluation ranks and ideal parts
    Ibasis = {}
    for D in range(2, dmax + 1):
        w = win(D)
        if Z:
            Kb = kernel_basis(w.eval_matrix(Z, p), p)
        else:
            Kb = None
        Ibasis[D] = Kb
        out['dimI'][D] = w.n if Kb is None else len(Kb)
    # filtered generation of the vanishing ideal
    gen_degree = 2
    if Z and with_generation:
        for D in range(2, k + 2):
            if D + 1 not in Ibasis:
                break
            src, dst = win(D), win(D + 1)
            target = out['dimI'][D + 1]
            E = Echelon(dst.n, p)
            base = embed(Ibasis[D], src, dst)
            E.add(base)
            for vi in range(len(dst.vars)):
                if E.rank == target:
                    break
                E.add(dst.shift(base, vi))
            assert E.rank <= target
            out['gen'][D] = (E.rank == target)
            if not out['gen'][D]:
                gen_degree = D + 1
    out['gen_degree'] = gen_degree
    # closures until V_D = I_D at a degree D >= gen_degree
    first_eq = lfd = None
    for D in range(2, dmax + 1):
        w = win(D)
        grow, g = generator_row(w, p, a, b, t)
        gdeg = max((bin(m[0]).count('1') + bin(m[1]).count('1') for m in g), default=0)
        E = closure(w, grow.reshape(1, -1) if gdeg <= D else np.zeros((0, w.n), dtype=np.int64), p)
        out['dimV'][D] = E.rank
        out['profile'][D] = profile(E, w)
        if Ibasis[D] is not None:
            assert E.rank <= out['dimI'][D]
        if E.rank == out['dimI'][D]:
            if first_eq is None:
                first_eq = D
            if D >= gen_degree:
                lfd = D
                break
    out['first_equal'] = first_eq
    out['lfd'] = lfd
    if first_eq is not None and first_eq < gen_degree:
        out['note'] = 'first_equal < gen_degree'
    return out


def mode_small(args):
    rng = random.Random(20261110)
    todo = []
    for p, k in [(5, 1), (7, 1), (13, 1), (5, 2), (13, 2), (17, 2), (17, 3)]:
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                for t in range(p):
                    todo.append((p, a, b, t, k))
    for p, k, count in [(41, 3, 200), (97, 3, 100), (17, 4, 200), (97, 4, 100), (113, 4, 60), (97, 5, 40), (193, 5, 40)]:
        for _ in range(count):
            while True:
                a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
                if (4 * a ** 3 + 27 * b * b) % p:
                    break
            todo.append((p, a, b, t, k))
    hist = collections.defaultdict(collections.Counter)
    pair = collections.defaultdict(collections.Counter)
    notes = 0
    t0 = time.time()
    with open(os.path.join(HERE, 'results-structure-small.jsonl'), 'w') as fh:
        for inst in todo:
            r = analyse(*inst)
            assert r['lfd'] is not None and r['lfd'] <= r['k'] + 2, r
            notes += 'note' in r
            hist[f"k={r['k']}"][r['lfd']] += 1
            if r['points']:
                pair[f"k={r['k']}"][f"first_equal={r['first_equal']},gen_degree={r['gen_degree']}"] += 1
            fh.write(json.dumps(r, sort_keys=True) + '\n')
    summary = dict(instances=len(todo), seconds=round(time.time() - t0, 1),
                   lfd_histogram={k: dict(sorted(v.items())) for k, v in hist.items()},
                   proper_first_equal_vs_generation={k: dict(v) for k, v in pair.items()},
                   first_equal_below_generation=notes, bound_k_plus_2_violations=0)
    json.dump(summary, open(os.path.join(HERE, 'results-structure-small-summary.json'), 'w'), indent=1, sort_keys=True)
    print(json.dumps(summary, indent=1, sort_keys=True))


def mode_sweep(args):
    root = os.path.join(HERE, '..', 'ic-curve-sweep')
    res = {}
    for line in open(os.path.join(root, 'results.jsonl')):
        r = json.loads(line)
        res[r['system_hash']] = r
    cases = [json.loads(l) for l in open(os.path.join(root, 'cases.jsonl'))]
    cases = [c for c in cases if c['kind'] == 'prime' and c['m'] == 2 and c['base_kind'] == 'subgroup'
             and c.get('base_scale', 1) == 1 and c['target'] is not None]
    todo = [c for c in cases if c['k'] <= 5]
    k6 = [c for c in cases if c['k'] == 6]
    todo += k6[::max(1, len(k6) // args.k6)][:args.k6]
    agree = mismatch = 0
    rows = []
    hist = collections.defaultdict(collections.Counter)
    t0 = time.time()
    for c in todo:
        r = analyse(c['prime'], c['a'], c['b'], c['target'][0], c['k'])
        ref = res[c['system_hash']]
        ok = (r['lfd'] == ref['exact_lfd'] and r['points'] == ref['point_count'])
        agree += ok
        mismatch += not ok
        hist[f"k={c['k']}"][r['lfd']] += 1
        rows.append(dict(case_id=c['case_id'], k=c['k'], prime=c['prime'], lfd=r['lfd'], sweep_lfd=ref['exact_lfd'],
                         points=r['points'], sweep_points=ref['point_count'], first_equal=r['first_equal'],
                         gen_degree=r['gen_degree'], agree=ok))
    summary = dict(cases=len(todo), agree=agree, mismatch=mismatch, seconds=round(time.time() - t0, 1),
                   lfd_histogram={k: dict(sorted(v.items())) for k, v in hist.items()})
    with open(os.path.join(HERE, 'results-structure-sweep.jsonl'), 'w') as fh:
        for r in rows:
            fh.write(json.dumps(r, sort_keys=True) + '\n')
    json.dump(summary, open(os.path.join(HERE, 'results-structure-sweep-summary.json'), 'w'), indent=1, sort_keys=True)
    print(json.dumps(summary, indent=1, sort_keys=True))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['small', 'sweep'])
    ap.add_argument('--k6', type=int, default=8)
    args = ap.parse_args()
    sys.exit(dict(small=mode_small, sweep=mode_sweep)[args.mode](args))
