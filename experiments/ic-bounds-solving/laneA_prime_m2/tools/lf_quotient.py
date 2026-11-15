"""Last/first fall degree of the two-summand prime chain system, computed in the quotient algebra
A = F_p[(Z/N)^2] (U1-upper-bounds.md section 1: V_D = least subspace of M_{<=D} containing g and closed
under u_{i,j} (translation by 2^j in block i) on elements of weight <= D-1).  Uses lfd_engine.Ech for the
exact mod-p linear algebra and the same round structure as lfd_engine.closure.  Unit instances:
d_lf = least D with 1 in V_D; d_ff = least D >= 2 with dim(V_D cap M_{<=D-1}) > dim V_{D-1}.
By the Bridge (U1 sec. 1, kernel-checked) this equals the ordinary-ring value for D >= 2; cross-checked
against lfd_engine.exact_lfd at k=3,4 in data/dlf_ordinary_vs_quotient.json.   [pcc:unverified]
"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import *


class QSpace:
    def __init__(self, k, D):
        N = 1 << k
        self.k, self.N, self.D = k, N, D
        pc = popcount_arr(N)
        W = (pc[:, None] + pc[None, :]).reshape(-1)           # weight of exponent pair index e1*N+e2
        idx = np.flatnonzero(W <= D)
        order = np.lexsort((idx, -W[idx]))                    # weight descending
        self.cols = idx[order]                                # exponent-pair index of each column
        self.n = len(self.cols)
        self.W = W[self.cols]
        self.colof = -np.ones(N * N, dtype=np.int64); self.colof[self.cols] = np.arange(self.n)
        below = np.flatnonzero(self.W <= D - 1)
        self.low_start = int(below[0]) if below.size else self.n
        e1 = self.cols // N; e2 = self.cols % N
        self.shift = []
        for i in range(2):
            for j in range(k):
                if i == 0: tgt = ((e1 + (1 << j)) % N) * N + e2
                else: tgt = e1 * N + (e2 + (1 << j)) % N
                self.shift.append(self.colof[tgt[self.low_start:]])
        assert (np.concatenate(self.shift) >= 0).all()

    def row_of_g(self, a, b, t, p):
        row = np.zeros(self.n)
        for (e1, e2), c in gen_coeffs(a, b, t, p).items():
            col = self.colof[e1 * self.N + e2]
            assert col >= 0
            row[col] = c % p
        return row


def closure_q(sp, a, b, t, p, deadline=None):
    V = LE.Ech(sp.n, p)
    done = LE.Ech(sp.n, p, cap=max(1, sp.n - sp.low_start))
    V.add(sp.row_of_g(a, b, t, p).reshape(1, -1))
    rounds = 0
    while True:
        lowmask = V.piv >= sp.low_start
        low = V.E[lowmask]
        if low.shape[0] == 0: break
        resid = done.reduce(low); resid = resid[np.any(resid != 0, axis=1)]
        if resid.shape[0] == 0: break
        new = done.add(resid); rounds += 1
        for blk in new:
            for sh in sp.shift:
                out = np.zeros((blk.shape[0], sp.n)); out[:, sh] = blk[:, sp.low_start:]
                V.add(out)
        if deadline and time.time() > deadline: return None
    lowcount = int(np.sum(V.piv >= sp.low_start))
    one = bool(sp.n - 1 in set(V.piv.tolist()))
    return dict(dim=V.rk, dim_low=lowcount, one=one, rounds=rounds, nmono=sp.n)


def profile(a, b, t, p, k, Dmax=None, max_seconds=1e9):
    """per-D data until 1 in V_D (unit instance).  returns d_lf, d_ff, per-D list."""
    Dmax = Dmax or k + 3
    t_end = time.time() + max_seconds
    per, dprev = [], 0
    d_lf = d_ff = None
    for D in range(2, Dmax + 1):
        t0 = time.time()
        r = closure_q(QSpace(k, D), a, b, t, p, deadline=t_end)
        if r is None: break
        r['D'] = D; r['sec'] = round(time.time() - t0, 2)
        per.append(r)
        if d_ff is None and r['dim_low'] > dprev: d_ff = D
        dprev = r['dim']
        if r['one']:
            d_lf = D; break
    return dict(d_lf=d_lf, d_ff=d_ff, per_degree=per)


if __name__ == '__main__':
    k, p, a, b, t = map(int, sys.argv[1:6])
    print(json.dumps(profile(a, b, t, p, k)))
