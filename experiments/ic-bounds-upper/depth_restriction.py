"""Depth restriction pi: A_{k+1} -> A_k (u_{i,k} -> 1) and the closures (Prop. 3 of U1).

For an instance with 2^(k+1) | p-1 and a degree D the script computes V_D(k), V_D(k+1) and

  contained   pi(V_D(k+1)) ⊆ V_D(k)                (proved; must be True)
  dim_pi      dim pi(V_D(k+1))
  wrap_gain   dim V_D(k) - dim pi(V_D(k+1))        (rows that need the wrap at level k)
  low_part    dim(V_D(k+1) ∩ M_{levels < k})       (elements not touching the new level)
  one         1 in V_D(k), 1 in V_D(k+1)
  pi_filtered dim pi(V_D(k+1) ∩ M_<=j), j = 0..D, to compare with profile_k = dim(V_D(k) ∩ M_<=j)

Usage: python3 depth_restriction.py
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from chain_algebra import Window, Echelon, closure, generator_row, profile, solutions

HERE = os.path.dirname(os.path.abspath(__file__))


def run(p, a, b, t, D, k, cache):
    def get(kk):
        key = (D, kk)
        if key not in cache:
            w = Window(kk, kk, D)
            grow, _ = generator_row(w, p, a, b, t)
            cache[key] = (w, closure(w, grow.reshape(1, -1), p))
        return cache[key]
    w0, E0 = get(k)
    w1, E1 = get(k + 1)
    N = 1 << k
    cols = np.array([w0.index[(e1 % N, e2 % N)] for (e1, e2) in w1.mons], dtype=np.int64)
    proj = np.zeros((E1.rank, w0.n), dtype=np.int64)
    np.add.at(proj.T, cols, E1.B.T)
    proj %= p
    contained = E0.contains(proj)
    Ep = Echelon(w0.n, p)
    Ep.add(proj)
    # filtered version: dim pi(V_D(k+1) ∩ M_<=j) against dim(V_D(k) ∩ M_<=j)
    filt = []
    pw = w1.weights[E1.piv]
    for j in range(D + 1):
        Ej = Echelon(w0.n, p)
        Ej.add(proj[pw <= j])
        filt.append(Ej.rank)
    high = np.array([i for i, (e1, e2) in enumerate(w1.mons) if e1 >= N or e2 >= N], dtype=np.int64)
    Eh = Echelon(len(high), p)
    Eh.add(E1.B[:, high])
    one0 = (w0.n - 1) in E0.piv
    one1 = (w1.n - 1) in E1.piv
    return dict(D=D, k=k, dimV_k=E0.rank, dimV_k1=E1.rank, window_k=w0.n, window_k1=w1.n,
                points_k=len(solutions(p, a, b, t, k)), points_k1=len(solutions(p, a, b, t, k + 1)),
                contained=bool(contained), dim_pi=Ep.rank, wrap_gain=E0.rank - Ep.rank,
                low_part=E1.rank - Eh.rank, one_k=bool(one0), one_k1=bool(one1),
                profile_k=profile(E0, w0), profile_k1=profile(E1, w1), pi_filtered=filt)


def main():
    insts = [dict(name='unit-12289', p=12289, a=3, b=4, t=5, kmax={3: 8, 4: 8, 5: 7}),
             dict(name='unit-40961', p=40961, a=21223, b=9887, t=4748, kmax={3: 8, 4: 8, 5: 7}),
             dict(name='proper-257', p=257, a=1, b=3, t=190, kmax={3: 8, 4: 8, 5: 7}),
             dict(name='proper-193', p=193, a=5, b=7, t=11, kmax={3: 6, 4: 6, 5: 6})]
    out = []
    t0 = time.time()
    for inst in insts:
        cache = {}
        for D in (3, 4, 5):
            for k in range(2, inst['kmax'][D]):
                r = run(inst['p'], inst['a'], inst['b'], inst['t'], D, k, cache)
                r.update(instance=inst['name'], p=inst['p'], a=inst['a'], b=inst['b'], t=inst['t'])
                assert r['contained'], r
                out.append(r)
                print(json.dumps({x: r[x] for x in ('instance', 'D', 'k', 'dimV_k', 'window_k', 'dimV_k1', 'window_k1',
                                                    'dim_pi', 'wrap_gain', 'low_part', 'one_k', 'one_k1',
                                                    'points_k', 'points_k1', 'profile_k', 'pi_filtered')}), flush=True)
            for key in [key for key in cache if key[0] == D]:
                del cache[key]
    # unit-case monotonicity: 1 in V_D(k+1) implies 1 in V_D(k)
    viol = [r for r in out if r['one_k1'] and not r['one_k']]
    json.dump(dict(rows=out, all_contained=all(r['contained'] for r in out), monotonicity_violations=len(viol),
                   seconds=round(time.time() - t0, 1)),
              open(os.path.join(HERE, 'results-depth-restriction.json'), 'w'), indent=1, sort_keys=True)
    print('violations', len(viol), 'seconds', round(time.time() - t0, 1))


if __name__ == '__main__':
    main()
