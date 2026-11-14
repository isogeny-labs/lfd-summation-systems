"""N8-growth, Section 7: the low part V_D(k) ∩ M_<D against the no-wrap ideal space S_{D-1}(k).

For each k: dim(V_D(k) ∩ M_<=D-1), dim S_{D-1}(k), and the two inclusions (exact: S from more than
4(2^k - 1) points of C); and the criterion of Proposition 7.1: the space S# = S_{D-1} + sum_u u S_{D-1}
(products in A_k, with wrap) contains nothing of weight < D outside S_{D-1}.
usage: ic_bounds_growth_lowpart.py p a b t D kmin kmax [out.json]
"""
import json, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import *

p, a, b, t, D, kmin, kmax = map(int, sys.argv[1:8])
g = s3_coeffs(p, a, b, t); pts = curve_points(p, a, b, t); out = []
for k in range(kmin, kmax + 1):
    L = Layout(k, D); N = L.N
    assert len(pts) >= 4 * N + 40
    B, piv = closure(L, g, p)
    lowcols = np.nonzero(L.wt < D)[0]
    Vlow = B[[r for r, c in enumerate(piv) if L.wt[c] < D]][:, lowcols]      # basis of V ∩ M_<D
    E = eval_matrix(L, pts[:4 * N + 40], p)[:, lowcols]
    S = nullspace(E, p); S, ps = rref(S, p)                                    # S_{D-1}(k)
    V_in_S = not ((Vlow @ E.T) % p).any()
    Rv, pv = rref(Vlow, p) if Vlow.shape[0] else (Vlow, [])
    S_in_V = not reduce_against(Rv, pv, S, p).any() if Vlow.shape[0] else S.shape[0] == 0
    # criterion of Proposition 7.1: S# = S_{D-1} + sum_u u S_{D-1} (products in A_k) has low part S_{D-1}
    Sfull = np.zeros((S.shape[0], L.n), dtype=np.int64); Sfull[:, lowcols] = S
    rows = [Sfull]
    for i in range(2):
        for j in range(k):
            pm = L.shift_perm(i, j)
            Y = np.zeros_like(Sfull); src = np.nonzero(pm >= 0)[0]
            Y[:, pm[src]] = Sfull[:, src]
            assert not Sfull[:, pm < 0].any()
            rows.append(Y)
    Sh, psh = rref(np.vstack(rows), p)
    sh_low = Sh[[r_ for r_, c in enumerate(psh) if L.wt[c] < D]][:, lowcols]
    crit = sh_low.shape[0] == S.shape[0]          # S_{D-1} is inside the low part, so equality of dimensions
    r = dict(p=p, a=a, b=b, t=t, D=D, k=k, dim_V=len(piv), one_in_V=bool(L.idx[(0, 0)] in piv),
             dim_V_low=int(Vlow.shape[0]), dim_S_Dminus1=int(S.shape[0]), Vlow_in_S=bool(V_in_S), S_in_Vlow=bool(S_in_V),
             dim_Ssharp=len(psh), dim_Ssharp_low=int(sh_low.shape[0]), criterion=bool(crit),
             one_in_Ssharp=bool(L.idx[(0, 0)] in psh))
    out.append(r); print(r, flush=True)
if len(sys.argv) > 8: json.dump(out, open(sys.argv[8], 'w'), indent=1)
