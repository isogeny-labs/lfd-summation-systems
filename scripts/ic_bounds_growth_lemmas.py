"""N8-growth: brute-force tests of the combinatorial lemmas of docs/proofs/N8-growth.md (Section 4) and
sanity tests of the cluster bound (Section 3) on the curve S_3(x1, x2, t) = 0.

  lambda   : table of Lambda(w) and k_D = D + Lambda(D-1); check Lambda(w) <= 13 w + 23 for w <= 120 and the
             monotonicity used in Lemma 4.1
  window   : Lemma 4.2 (window lemma) by brute force over small (N, p, sigma)
  straddle : Lemma 4.4 (popcount of a straddling pair) by brute force
  sparse   : Theorem 3.3 (a vanishing sum with m terms that has no proper F_0-dependent subset has exponent
             spread <= 3(m-1)(m-2)); tested in the contrapositive form "few monomials with large spread are
             linearly independent on C" on random and structured monomial sets, exponents far below p
  semaev   : S_3(x(P), x(R-P), x(R)) = 0 on random points
usage: ic_bounds_growth_lemmas.py [lambda|window|straddle|sparse|semaev|all] [out.json]
"""
import json, random, sys
from math import comb
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import curve_points, rref, popcount

CONST = 7776

def Cle(n, w):
    return sum(comb(n, i) for i in range(w + 1))

def Lam(w, top=None):
    top = top or 40 * w + 200
    best = 0
    for l in range(1, top):
        if 2 ** (l - 1) < CONST * Cle(2 * l, w) ** 2:
            best = l
    return best

def t_lambda():
    rows = []; ok = True
    for w in range(1, 121):
        L = Lam(w)
        # initial-segment property on [4w, top): once false, stays false
        seen_false = False
        for l in range(4 * w, 40 * w + 200):
            holds = 2 ** (l - 1) < CONST * Cle(2 * l, w) ** 2
            if seen_false and holds: ok = False
            if not holds: seen_false = True
        if L > 13 * w + 23 or L < 4 * w: ok = False
        if w <= 20 or w % 20 == 0: rows.append(dict(w=w, D=w + 1, Lambda=L, k_D=w + 1 + L))
    return dict(test='lambda', ok=ok, rows=rows)

def dist(n, p):
    n %= p
    return min(n, p - n)

def t_window():
    """{a in [0,N): ||2(a-a0)||_p <= 2 sigma} lies in <= 6 aligned blocks of size 2^l, l = ceil(log2(2 sigma+1)),
    for p >= 2N+1."""
    checked = 0; worst = 0
    for k in range(2, 9):
        N = 1 << k
        for p in [q for q in range(2 * N + 1, 12 * N, 2) if all(q % d for d in range(3, int(q ** .5) + 1, 2))][:12]:
            for sigma in [0, 1, 2, 3, 5, 6, 9, 18, 36, 60]:
                l = max(0, (2 * sigma).bit_length())   # 2^l >= 2 sigma + 1
                assert (1 << l) >= 2 * sigma + 1
                for a0 in range(N):
                    blocks = {a >> l for a in range(N) if dist(2 * (a - a0), p) <= 2 * sigma}
                    worst = max(worst, len(blocks)); checked += 1
    return dict(test='window', ok=worst <= 6, cases=checked, max_blocks=worst)

def t_straddle():
    """For k, Lam <= k, j >= Lam-1, sigma < 2^(Lam-1), p >= 2N+1: if a >= N-2^j > a' and
    (a-a' <= sigma or a-a' >= (p-2 sigma)/2) then popcount(a') >= k-Lam or popcount(a) >= k-Lam+1."""
    checked = 0; ok = True
    for k in range(2, 9):
        N = 1 << k
        for p in (2 * N + 1, 3 * N + 1, 5 * N + 1):      # primality is not used by the lemma
            for Lm in range(1, k + 1):
                sigma = (1 << (Lm - 1)) - 1
                for j in range(max(Lm - 1, 0), k):
                    for a in range(N - (1 << j), N):
                        for a2 in range(0, N - (1 << j)):
                            d = a - a2
                            if d <= sigma:
                                checked += 1
                                if popcount(a2) < k - Lm: ok = False
                            elif 2 * d >= p - 2 * sigma:
                                checked += 1
                                if popcount(a) < k - Lm + 1: ok = False
    return dict(test='straddle', ok=ok, cases=checked)

def t_sparse(seed=1):
    """Sets of m <= 7 monomials with exponent spread > 3(m-1)(m-2) in some pair (|da|+|db| measured between the
    two extreme monomials) contain no vanishing sum with all m coefficients nonzero whose proper subsets are
    independent; in particular: if EVERY pair of the set has spread > 3(m-1)(m-2) the set is linearly
    independent on C.  Tested by rank over F_p on all affine F_p-points of C (exponents <= 400, p = 12289,
    40961; a nonzero polynomial of bidegree <= 400 has at most 1600 zeros on C, the curves have > 12000 points)."""
    rng = random.Random(seed); res = []
    for (p, a, b, t) in [(12289, 4021, 9917, 4408), (12289, 1, 0, 9362), (40961, 21223, 9887, 4748)]:
        pts = curve_points(p, a, b, t)
        x1 = np.array([q[0] for q in pts], dtype=object); x2 = np.array([q[1] for q in pts], dtype=object)
        def col(e):
            return np.array([pow(int(u), e[0], p) * pow(int(v), e[1], p) % p for u, v in zip(x1[:2000], x2[:2000])], dtype=np.int64)
        tested = dep = 0
        for trial in range(400):
            m = rng.choice([3, 4, 5, 6, 7]); s = 3 * (m - 1) * (m - 2)
            kind = trial % 3
            for attempt in range(200):
                if kind == 0:      # random
                    M = [(rng.randrange(401), rng.randrange(401)) for _ in range(m)]
                elif kind == 1:    # lattice (Graeffe-like): multiples of a step
                    st = rng.choice([s + 1, s + 2, 2 * s + 1])
                    M = [(st * rng.randrange(400 // st + 1), st * rng.randrange(400 // st + 1)) for _ in range(m)]
                else:              # powers of two only
                    M = [(rng.choice([0] + [1 << i for i in range(9)]), rng.choice([0] + [1 << i for i in range(9)])) for _ in range(m)]
                if len(set(M)) < m: continue
                if all(abs(u[0] - v[0]) + abs(u[1] - v[1]) > s for i, u in enumerate(M) for v in M[:i]): break
            else:
                continue
            E = np.stack([col(e) for e in M], axis=1)
            _, piv = rref(E, p)
            tested += 1
            if len(piv) < m: dep += 1
        res.append(dict(p=p, a=a, b=b, t=t, points=len(pts), sets_tested=tested, dependent=dep))
    return dict(test='sparse', ok=all(r['dependent'] == 0 for r in res), rows=res)

def t_semaev(seed=5):
    """S_3(x(P), x(R - P), x(R)) = 0 for the formula of the repository (chord-tangent arithmetic mod p)."""
    from ic_bounds_growth_lib import s3_coeffs
    rng = random.Random(seed); tested = 0; ok = True
    for (p, a, b) in [(12289, 4021, 9917), (12289, 1, 0), (40961, 21223, 9887), (10007, 3, 5)]:
        sq = {}
        for y in range(p): sq.setdefault(y * y % p, y)
        def rand_pt():
            while True:
                x = rng.randrange(p); r = (x * x * x + a * x + b) % p
                if r in sq: return (x, sq[r])
        def add(P, Q):
            if P is None: return Q
            if Q is None: return P
            if P[0] == Q[0]:
                if (P[1] + Q[1]) % p == 0: return None
                lam = (3 * P[0] * P[0] + a) * pow(2 * P[1], p - 2, p) % p
            else:
                lam = (Q[1] - P[1]) * pow(Q[0] - P[0], p - 2, p) % p
            x = (lam * lam - P[0] - Q[0]) % p
            return (x, (lam * (P[0] - x) - P[1]) % p)
        for _ in range(500):
            R = rand_pt(); P = rand_pt(); Q = add(R, (P[0], (-P[1]) % p))
            if Q is None: continue
            g = s3_coeffs(p, a, b, R[0])
            val = sum(c * pow(P[0], e1, p) * pow(Q[0], e2, p) for (e1, e2), c in g.items()) % p
            tested += 1; ok &= (val == 0)
    return dict(test='semaev', ok=bool(ok), cases=tested)

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    tests = dict(**{'lambda': t_lambda}, window=t_window, straddle=t_straddle, sparse=t_sparse, semaev=t_semaev)
    out = []
    for name, fn in tests.items():
        if which in ('all', name):
            r = fn(); out.append(r); print(json.dumps(r)[:1500], flush=True)
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], 'w'), indent=1)
