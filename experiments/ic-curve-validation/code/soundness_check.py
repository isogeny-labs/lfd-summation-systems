#!/usr/bin/env python3
"""Direct soundness / certificate check for every proper (Z non-empty) system of the audit:
at the certified degree D0 recompute V_{D0} and verify (a) every basis row of V_{D0} vanishes on every point of Z
(so V_{D0} is inside the vanishing ideal, independent of the closure code), (b) dim V_{D0} = |M_{<=D0}| - |Z|,
(c) the monomials of degree < D0 evaluate onto K^Z (rank = |Z|).  Together: V_{D0} = I_{D0} and LFD <= D0."""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lfd_engine as L
import run_checks as RC
S, R, C = RC.load()
n = bad = 0
for f in sorted(glob.glob('final_*.jsonl')):
    for line in open(f):
        rec = json.loads(line)
        if rec['mine']['points'] == 0:
            continue
        h = rec['system_hash']; s = S[h]; D0 = rec['mine']['certified_at']
        if s['kind'] == 'prime':
            k, p = s['k'], s['prime']; r = 2 * k
            inputs, G = L.prime_inputs([(int(m), int(c)) for m, c in s['generators'][0]], k, p)
            pts = L.prime_points(G, k, p)
        else:
            r = s['variables']; p = 2
            inputs = L.binary_inputs(s['equations'], r); pts = L.binary_points(s['equations'], r)
        res = L.closure(inputs, r, D0, p)
        sp = res['space']
        M = L.eval_matrix(pts, sp, p)                       # |Z| x n
        vanish = not np.any(np.mod(M @ res['V'].E.T.astype(np.int64), p))
        dim_ok = res['dim'] == sp.n - len(pts)
        cols_lt = sp.deg <= D0 - 1
        onto = L.rank_mod(M[:, cols_lt], p) == len(pts)
        n += 1
        if not (vanish and dim_ok and onto):
            bad += 1
            print('FAIL', h, vanish, dim_ok, onto)
print('proper systems checked', n, 'failures', bad)
