"""Spot check of P3 / E1 in the ORDINARY ring: sigma = (u_{i,k} -> 1) maps the degree-D closure at
depth k+1 into the degree-D closure at depth k (P3), and onto it (E1)."""
import numpy as np
from rt_lib import *
for (p, a, b, t) in ((12289, 3, 4, 5), (40961, 21223, 9887, 4748), (257, 1, 3, 190), (193, 5, 7, 11)):
    for D in (3, 4):
        for k in (2, 3, 4):
            if (p - 1) % (1 << (k + 1)):
                continue
            inp1, pts1, r1 = build(p, a, b, t, k + 1, k + 1)
            inp0, pts0, r0 = build(p, a, b, t, k, k)
            c1 = E.closure(inp1, r1, D, p)
            c0 = E.closure(inp0, r0, D, p)
            sp1, sp0 = c1['space'], c0['space']
            keep = [i for i in range(r1) if i not in (k, 2 * k + 1)]
            cols = np.array([sp0.idx[tuple(m[i] for i in keep)] for m in sp1.cols], dtype=np.int64)
            img = np.zeros((c1['V'].rk, sp0.n), dtype=np.float64)
            np.add.at(img.T, cols, c1['V'].E.T)
            img = np.mod(img, p)
            inside = not np.any(c0['V'].reduce(img))
            rk = E.rank_mod(img.astype(np.int64), p)
            print(f'p={p} D={D} k={k}->{k+1}: |Z_k|={len(pts0)} |Z_k+1|={len(pts1)} dim V(k+1)={c1["dim"]} dim V(k)={c0["dim"]} '
                  f'image inside V(k): {inside}  rank image={rk}  onto: {rk == c0["dim"]}  1 in V(k+1)/V(k): {c1["contains_one"]}/{c0["contains_one"]}', flush=True)
