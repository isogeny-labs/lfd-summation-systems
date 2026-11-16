"""Task 5: proper instances (|Z|>0), group-algebra static analogue
   d_full = least D with  dim(g * M_{<=D-2}) = dim(I cap M_{<=D})  and  eval on M_{<D} onto K^Z,
   dim I_D = |M_<=D| - rank(eval of M_<=D on Z),  compared with the exact d_lf (quotient closure, certificate as
   in lfd_engine.exact_lfd:  dimV_D = dimI_D and eval of M_{<D} onto K^Z).
Usage: python3 task5_dfull.py k p n [seed]   [pcc:unverified]"""
import sys, os, json, random, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import lf_quotient as Q


def zero_points(a, b, t, p, k):
    N = 1 << k
    w = root_of_unity(p, N)
    H = [pow(w, i, p) for i in range(N)]
    return [(x, y) for x in H for y in H if S3(x, y, t, a, b, p) == 0]


def analyse(a, b, t, p, k, Z, Dmax):
    N = 1 << k
    pc = popcount_arr(N)
    W = (pc[:, None] + pc[None, :]).reshape(-1)
    e1 = np.arange(N * N) // N; e2 = np.arange(N * N) % N
    ev = np.array([[pow(x, int(u), p) * pow(y, int(v), p) % p for u, v in zip(e1, e2)] for (x, y) in Z], dtype=np.int64)
    gc = gen_coeffs(a, b, t, p)
    out = dict(per=[])
    rank_eval = {}
    DM = 2 * k + 2   # static scan range; d_lf scan stops at its certificate
    for D in range(0, DM + 1):
        cols = np.flatnonzero(W <= D)
        rank_eval[D] = LE.rank_mod(ev[:, cols], p)
    # Mac quotient: rows tau_mu g, wt(mu) <= D-2
    d_full = None; V_prof = None
    for D in range(2, DM + 1):
        sp = Q.QSpace(k, D)
        mus = np.flatnonzero(W <= D - 2)
        rows = np.zeros((len(mus), sp.n))
        for i, mu in enumerate(mus):
            m1, m2 = mu // N, mu % N
            for (u, v), c in gc.items():
                rows[i, sp.colof[((u + m1) % N) * N + (v + m2) % N]] = c % p
        ech = LE.Ech(sp.n, p, cap=min(len(mus), sp.n)); ech.add(rows)
        dimMac = ech.rk
        dimI = sp.n - rank_eval[D]
        onto = (rank_eval[D - 1] == len(Z))
        rec = dict(D=D, dimMac=dimMac, dimI=dimI, onto_below=onto, mac_eq_I=(dimMac == dimI))
        out['per'].append(rec)
        if dimMac == dimI and onto and d_full is None: d_full = D
    out['d_full'] = d_full
    out['d_full_stable'] = (d_full is not None and all(r['mac_eq_I'] and r['onto_below'] for r in out['per'] if r['D'] >= d_full))
    # exact d_lf via quotient closure
    neq = []; cert = None; lf = []
    for D in range(2, Dmax + 1):
        sp = Q.QSpace(k, D)
        r = Q.closure_q(sp, a, b, t, p)
        dimI = sp.n - rank_eval[D]
        assert r['dim'] <= dimI
        lf.append(dict(D=D, dimV=r['dim'], dimI=dimI))
        if r['dim'] != dimI: neq.append(D)
        if r['dim'] == dimI and rank_eval[D - 1] == len(Z):
            cert = D; break
    out['d_lf'] = (max(neq) + 1) if (cert is not None and neq) else (2 if cert is not None else None)
    out['d_lf_certified_at'] = cert; out['lf_per'] = lf
    return out


if __name__ == '__main__':
    k, p, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    rng = random.Random(seed * 99991 + p + k)
    rows, tries = [], 0
    while len(rows) < n and tries < 5000:
        tries += 1
        a, b = random_curve(p, rng); t = rng.randrange(p)
        Z = zero_points(a, b, t, p, k)
        if not (0 < len(Z) <= 40): continue
        t0 = time.time()
        r = analyse(a, b, t, p, k, Z, Dmax=k + 4)
        r.update(p=p, a=a, b=b, t=t, k=k, nZ=len(Z), sec=round(time.time() - t0, 1))
        rows.append(r)
        print(dict(k=k, p=p, a=a, b=b, t=t, nZ=len(Z), d_full=r['d_full'], stable=r['d_full_stable'], d_lf=r['d_lf'], sec=r['sec']), flush=True)
    json.dump(rows, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'task5_k%d_p%d.json' % (k, p)), 'w'), indent=1)
