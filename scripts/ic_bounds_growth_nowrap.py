"""N8-growth, Section 7: is V_D(k) the low part of a wrap-free closure?

In A_L a shift x_i^(2^j) with j <= L - D never wraps on an element of weight < D (a wrapping monomial has the
L - j >= D bits j..L-1 set).  So the closure T^-(L) of g in A_L under the shifts with j <= L - D only is a
derivation in the polynomial ring K[x1, x2] (integer exponents, no reduction mod x^N - 1).  The script
compares, for k <= L - D + 1,

    dim( T^-(L) restricted to exponents < 2^k )   with   dim V_D(k)   (closure in A_k, with wrap).

The inclusion "low part of a wrap-free closure, reduced mod 2^k, lies in V_D(k)" always holds; equality of
the dimensions means that every element of V_D(k) is obtained without any wrap-around.
usage: ic_bounds_growth_nowrap.py p a b t D k L [out.json]
"""
import json, sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import *

def main():
    p, a, b, t, D, k, Ltop = map(int, sys.argv[1:8])
    g = s3_coeffs(p, a, b, t)
    Lk = Layout(k, D); Bk, pk = closure(Lk, g, p)
    out = dict(p=p, a=a, b=b, t=t, D=D, k=k, dim_V_k=len(pk), profile_V_k=profile(Lk, pk), rows=[])
    print(out, flush=True)
    for L in range(max(k, D + 1), Ltop + 1):
        t0 = time.time()
        LL = Layout(L, D)
        shifts = [(i, j) for i in range(2) for j in range(L - D + 1)]
        B, piv = closure(LL, g, p, shifts=shifts)
        # low part: elements supported on exponents < 2^k.  Eliminate the other columns first.
        N = 1 << k
        lowcol = (LL.a1 < N) & (LL.a2 < N)
        order = np.concatenate([np.nonzero(~lowcol)[0], np.nonzero(lowcol)[0]])
        R, pv = rref(B[:, order], p)
        nhigh = int((~lowcol).sum())
        lowrows = R[[r for r, c in enumerate(pv) if c >= nhigh]][:, nhigh:]
        # express in the layout of depth k and compare with V_D(k)
        cols_k = [Lk.idx[LL.mons[c]] for c in order[nhigh:]]
        M = np.zeros((lowrows.shape[0], Lk.n), dtype=np.int64); M[:, cols_k] = lowrows
        res = reduce_against(Bk, pk, M, p)
        inside = not res.any()
        _, pm = rref(M, p) if M.shape[0] else (None, [])
        row = dict(L=L, shifts_levels=L - D + 1, dim_Tminus=len(piv), profile_Tminus=profile(LL, piv),
                   dim_low_part=len(pm), low_part_inside_V=inside, equal=(len(pm) == len(pk)) and inside,
                   profile_low_part=profile(Lk, pm), seconds=round(time.time() - t0, 1))
        out['rows'].append(row); print(row, flush=True)
    if len(sys.argv) > 8: json.dump(out, open(sys.argv[8], 'w'), indent=1)

main()
