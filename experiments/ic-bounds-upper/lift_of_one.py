"""Least degree j such that V_D(k+1) ∩ M_<=j contains an element restricting to 1 on H_k^2
(i.e. 1 in pi(V_D(k+1) ∩ M_<=j)), at the transition D = 4, k = 5 -> 6 (U1, Section 5.2, T2)."""
import json
import os

import numpy as np

from chain_algebra import Window, Echelon, closure, generator_row

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
for (p, a, b, t) in [(12289, 3, 4, 5), (40961, 21223, 9887, 4748), (12289, 1, 0, 9362)]:
    D, k = 4, 5
    w0, w1 = Window(k, k, D), Window(k + 1, k + 1, D)
    grow, _ = generator_row(w1, p, a, b, t)
    E1 = closure(w1, grow.reshape(1, -1), p)
    N = 1 << k
    cols = np.array([w0.index[(e1 % N, e2 % N)] for (e1, e2) in w1.mons], dtype=np.int64)
    proj = np.zeros((E1.rank, w0.n), dtype=np.int64)
    np.add.at(proj.T, cols, E1.B.T)
    proj %= p
    pw = w1.weights[E1.piv]
    one = np.zeros((1, w0.n), dtype=np.int64)
    one[0, w0.index[(0, 0)]] = 1
    least = None
    for j in range(D + 1):
        Ej = Echelon(w0.n, p)
        Ej.add(proj[pw <= j])
        if Ej.contains(one):
            least = j
            break
    rec = dict(p=p, a=a, b=b, t=t, D=D, k=k, dimV_k1=E1.rank, one_in_V_k1=bool((w1.n - 1) in E1.piv),
               least_degree_of_a_lift_of_one=least)
    print(rec)
    out.append(rec)
json.dump(out, open(os.path.join(HERE, 'results-lift-of-one.json'), 'w'), indent=1)
