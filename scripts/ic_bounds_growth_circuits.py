"""N8-growth: Theorem 3.3 on actual low-weight elements of the ideal of C.

For (w, k): the space S_w(k) of polynomials of weight <= w, exponents < 2^k, vanishing on C is computed as the
kernel of the evaluation matrix; each row of a reduced echelon basis is a minimal vanishing sum over the
constants (a circuit of the linear matroid of the monomials).  Several random column orders give different
circuits.  For each circuit with m terms the script records spread = max over pairs of |a-a'| + |b-b'| and
checks spread <= 3(m-1)(m-2)  (Theorem 3.3 for p = 0; here p is a prime far above the exponents, and a circuit
over the constants need not be a circuit over F^p, so this is a consistency test, not a consequence).
It also reports the least number of terms of a circuit and the largest ratio spread / (number of terms).
usage: ic_bounds_growth_circuits.py p a b t w kmin kmax [out.json]
"""
import json, random, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import *

def run(p, a, b, t, w, k, orders=4, seed=3):
    L = Layout(k, w); N = L.N
    pts = curve_points(p, a, b, t)
    need = 4 * (N - 1) + 41
    assert len(pts) >= need
    E = eval_matrix(L, pts[:need], p)
    rng = random.Random(seed)
    rec = dict(p=p, a=a, b=b, t=t, w=w, k=k, columns=L.n, circuits=0, violations=0, min_terms=None,
               min_slack=None, max_spread_per_term=0.0, sparsest=None)
    for o in range(orders):
        perm = np.arange(L.n)
        if o: rng.shuffle(perm)
        Nb = nullspace(E[:, perm], p)
        if o == 0: rec['dim_S'] = int(Nb.shape[0])
        for row in Nb:
            sup = perm[np.nonzero(row)[0]]
            m = len(sup)
            A = L.a1[sup]; B = L.a2[sup]
            # max over pairs of |da| + |db| = max over the four sign patterns of the range of (+-a +- b)
            spread = int(max(np.ptp(A + B), np.ptp(A - B)))
            slack = 3 * (m - 1) * (m - 2) - spread
            rec['circuits'] += 1
            if slack < 0: rec['violations'] += 1
            if rec['min_slack'] is None or slack < rec['min_slack']: rec['min_slack'] = slack
            if rec['min_terms'] is None or m < rec['min_terms']:
                rec['min_terms'] = m; rec['sparsest'] = sorted((int(x), int(y)) for x, y in zip(A, B))
            rec['max_spread_per_term'] = max(rec['max_spread_per_term'], round(spread / m, 2))
    return rec

if __name__ == '__main__':
    p, a, b, t, w, kmin, kmax = map(int, sys.argv[1:8])
    out = []
    for k in range(kmin, kmax + 1):
        r = run(p, a, b, t, w, k); out.append(r); print(r, flush=True)
    if len(sys.argv) > 8: json.dump(out, open(sys.argv[8], 'w'), indent=1)
