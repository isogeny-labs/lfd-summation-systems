"""N8-growth, exact tests of the no-wrap ideal space S_D(k) = {f in M_<=D(k) : lift(f) vanishes on the curve C}.

For each (D, k): dim V_D(k), dim S_D(k), whether V_D(k) is inside S_D(k), whether S_D(k) is closed under the
degree-D rule of A_k, and the splitting property FEI_k(w), w = D-1 (every top-run part of every element of
S of weight <= w is again in the ideal).  A polynomial of bidegree < 2^k vanishing at more than 4(2^k - 1)
points of the irreducible (2,2)-curve C lies in its ideal, so with enough points the tests are exact.

usage: ic_bounds_growth_ideal.py p a b t D kmin kmax [out.json]
"""
import json, sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import *

def run(p, a, b, t, D, k, with_closure=True):
    L = Layout(k, D); N = L.N
    g = s3_coeffs(p, a, b, t)
    pts = curve_points(p, a, b, t)
    need = 4 * (N - 1) + 1
    rec = dict(p=p, a=a, b=b, t=t, D=D, k=k, columns=L.n, points=len(pts), points_needed=need,
               exact=len(pts) >= need)
    E = eval_matrix(L, pts[:max(need + 40, 0) if len(pts) >= need + 40 else len(pts)], p)
    R, pe = rref(E, p)
    rec['rank_eval'] = len(pe); rec['dim_S'] = L.n - len(pe)
    S = nullspace(E, p)                      # RREF basis of S, columns in weight-descending order
    # nullspace basis: each row has a 1 at a free column and entries at pivot columns; bring to RREF
    S, ps = rref(S, p)
    rec['profile_S'] = profile(L, ps)
    rec['one_in_S'] = bool(L.idx[(0, 0)] in ps)
    # closedness of S in A_k
    low = S[[r for r, c in enumerate(ps) if L.wt[c] < D]]
    bad = {}
    for i in range(2):
        for j in range(k):
            pm = L.shift_perm(i, j)
            Y = np.zeros_like(low); src = np.nonzero(pm >= 0)[0]
            Y[:, pm[src]] = low[:, src]
            res = reduce_against(S, ps, Y, p)
            nb = int(np.any(res != 0, axis=1).sum())
            if nb: bad[f'{i},{j}'] = nb
    rec['S_closed'] = not bad; rec['S_not_closed_shifts'] = bad
    # FEI: top-run parts of low elements of S vanish on C
    fei_bad = {}
    for i in range(2):
        ai = L.a1 if i == 0 else L.a2
        for j in range(k):
            mask = (ai >= N - (1 << j)).astype(np.int64)
            hi = low * mask
            v = (hi @ E.T) % p
            nb = int(np.any(v != 0, axis=1).sum())
            if nb: fei_bad[f'{i},{j}'] = nb
    rec['FEI'] = not fei_bad; rec['FEI_fail'] = fei_bad
    if with_closure:
        t0 = time.time()
        B, piv = closure(L, g, p)
        rec['dim_V'] = len(piv); rec['profile_V'] = profile(L, piv)
        rec['one_in_V'] = bool(L.idx[(0, 0)] in piv)
        rec['V_sweep_failures'] = verify_closed(L, B, piv, p)
        v = (B @ E.T) % p
        rec['V_in_S'] = not v.any()
        rec['V_rows_outside_S'] = int(np.any(v != 0, axis=1).sum())
        rec['closure_seconds'] = round(time.time() - t0, 1)
    return rec

if __name__ == '__main__':
    p, a, b, t, D, kmin, kmax = map(int, sys.argv[1:8])
    out = []
    for k in range(kmin, kmax + 1):
        r = run(p, a, b, t, D, k)
        out.append(r)
        print({x: r[x] for x in ('D', 'k', 'columns', 'exact', 'rank_eval', 'dim_S', 'dim_V', 'one_in_V', 'V_sweep_failures', 'V_in_S',
                                 'S_closed', 'FEI', 'profile_V', 'profile_S')}, flush=True)
        if r['S_not_closed_shifts']: print('   S fails at shifts', r['S_not_closed_shifts'], ' FEI fails', r['FEI_fail'], flush=True)
    if len(sys.argv) > 8:
        json.dump(out, open(sys.argv[8], 'w'), indent=1)
