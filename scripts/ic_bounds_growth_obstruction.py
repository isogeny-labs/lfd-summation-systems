"""N8-growth, Section 6.2: a basis-independent measure of the failure of FEI_k(D-1).

For each coordinate i and level j the top-run map  tau_{i,j} : S_{D-1}(k) -> F,  f -> (top-run part of f) on
the curve, is linear; FEI_k(D-1) says that all tau_{i,j} vanish.  The script reports rank(tau_{i,j}) (the
dimension of the image of the top-run map modulo the ideal), computed as the rank of the evaluations of the
top-run parts of a basis of S_{D-1}(k) on more than 4(2^k - 1) points of C.  Unlike the number of failing
basis rows printed by ic_bounds_growth_ideal.py, this does not depend on the basis.
usage: ic_bounds_growth_obstruction.py p a b t D kmin kmax [out.json]
"""
import json, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import *

p, a, b, t, D, kmin, kmax = map(int, sys.argv[1:8])
pts = curve_points(p, a, b, t); out = []
for k in range(kmin, kmax + 1):
    L = Layout(k, D - 1); N = L.N
    assert len(pts) >= 4 * N + 40
    E = eval_matrix(L, pts[:4 * N + 40], p)
    S = nullspace(E, p)
    ranks = {}
    for i in range(2):
        ai = L.a1 if i == 0 else L.a2
        for j in range(k):
            mask = (ai >= N - (1 << j)).astype(np.int64)
            v = ((S * mask) @ E.T) % p
            r = len(rref(v, p)[1]) if v.any() else 0
            if r: ranks[f'{i},{j}'] = r
    rec = dict(p=p, a=a, b=b, t=t, D=D, k=k, dim_S_Dminus1=int(S.shape[0]), FEI=not ranks, top_run_map_ranks=ranks)
    out.append(rec); print(rec, flush=True)
if len(sys.argv) > 8: json.dump(out, open(sys.argv[8], 'w'), indent=1)
