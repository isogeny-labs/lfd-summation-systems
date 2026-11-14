"""N8-growth: attack on the characteristic-p part of Theorem 3.3 (docs/proofs/N8-growth.md).

Claim tested (consequence of Theorem 3.3 with F_0 = F^p): let M be a set of m monomials x1^a x2^b such that
every two of them satisfy ||2(a-a')||_p + ||2(b-b')||_p > 6(m-1)(m-2)  (|| || = distance to pZ).  Then M is
linearly independent over the constants on the curve C: S_3(x1, x2, t) = 0, t != 0.
Exponents range up to 4p, where Frobenius relations (g^p = g(x1^p, x2^p), exponent differences divisible by p)
exist; such pairs are close modulo p and are excluded by the hypothesis.  Control: sets containing the
support of g(x1^p, x2^p) ARE dependent (the script checks that it detects this).

Points of C over F_{p^2} (so that x^p != x on the sample); independence = full column rank of the evaluation
matrix on 60 random points (rigorous for independence); a rank defect is re-examined on all points.
Limit of this test (review finding F7): it decides dependence over the CONSTANTS only.  Dependence over F^p with
non-constant coefficients is tested in ic_bounds_growth_fpdep.py.
usage: ic_bounds_growth_charp.py [out.json]
"""
import json, random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import s3_coeffs

class Fq2:
    def __init__(self, p):
        self.p = p
        self.nr = next(n for n in range(2, p) if pow(n, (p - 1) // 2, p) == p - 1)
    def mul(self, x, y):
        p = self.p
        return ((x[0] * y[0] + self.nr * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)
    def add(self, x, y):
        return ((x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p)
    def neg(self, x):
        return ((-x[0]) % self.p, (-x[1]) % self.p)
    def inv(self, x):
        p = self.p
        n = (x[0] * x[0] - self.nr * x[1] * x[1]) % p
        ni = pow(n, p - 2, p)
        return (x[0] * ni % p, (-x[1]) * ni % p)
    def pw(self, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = self.mul(r, x)
            x = self.mul(x, x); e >>= 1
        return r

def curve_points_q2(K, a, b, t):
    p = K.p; g = s3_coeffs(p, a, b, t)
    els = [(u, v) for u in range(p) for v in range(p)]
    sq = {}
    for y in els: sq.setdefault(K.mul(y, y), y)
    pts = []
    for x1 in els:
        co = [(0, 0)] * 3
        pw = [(1, 0), x1, K.mul(x1, x1)]
        for (e1, e2), c in g.items():
            co[e2] = K.add(co[e2], K.mul((c, 0), pw[e1]))
        c0, c1, c2 = co
        if c2 == (0, 0): continue
        disc = K.add(K.mul(c1, c1), K.neg(K.mul((4, 0), K.mul(c2, c0))))
        if disc not in sq: continue
        s = sq[disc]; inv = K.inv(K.mul((2, 0), c2))
        for r in {s, K.neg(s)}:
            pts.append((x1, K.mul(K.add(K.neg(c1), r), inv)))
    return pts

def rank_q2(K, rows):
    rows = [list(r) for r in rows]; rk = 0; ncol = len(rows[0])
    for c in range(ncol):
        s = next((i for i in range(rk, len(rows)) if rows[i][c] != (0, 0)), None)
        if s is None: continue
        rows[rk], rows[s] = rows[s], rows[rk]
        iv = K.inv(rows[rk][c]); rows[rk] = [K.mul(iv, v) for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != (0, 0):
                f = rows[i][c]
                rows[i] = [K.add(v, K.neg(K.mul(f, w))) for v, w in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def dist(n, p):
    n %= p
    return min(n, p - n)

def main():
    rng = random.Random(7); out = []
    for (p, a, b, t) in [(101, 3, 7, 5), (103, 1, 0, 11), (127, 5, 9, 2)]:
        K = Fq2(p); pts = curve_points_q2(K, a, b, t)
        def rank_of(M, sample):
            return rank_q2(K, [[K.mul(K.pw(P[0], e[0]), K.pw(P[1], e[1])) for e in M] for P in sample])
        tested = dep = 0
        for trial in range(300):
            m = rng.choice([3, 3, 4]); thr = 6 * (m - 1) * (m - 2)
            kind = trial % 3
            for attempt in range(500):
                if kind == 0:
                    M = [(rng.randrange(4 * p), rng.randrange(4 * p)) for _ in range(m)]
                elif kind == 1:     # near multiples of p (Frobenius-like), offsets just beyond the threshold
                    M = [(p * rng.randrange(4) + rng.randrange(-thr - 3, thr + 4), p * rng.randrange(4) + rng.randrange(-thr - 3, thr + 4)) for _ in range(m)]
                    if any(u < 0 or v < 0 for u, v in M): continue
                else:               # one coordinate fixed
                    c0 = rng.randrange(4 * p)
                    M = [(rng.randrange(4 * p), c0) for _ in range(m)]
                if len(set(M)) < m: continue
                if all(dist(2 * (u[0] - v[0]), p) + dist(2 * (u[1] - v[1]), p) > thr for i, u in enumerate(M) for v in M[:i]): break
            else:
                continue
            tested += 1
            if rank_of(M, rng.sample(pts, 60)) < m:
                if rank_of(M, pts[:4000]) < m: dep += 1
        # control: the Frobenius relation is detected
        g = s3_coeffs(p, a, b, t)
        Mf = [(p * e1, p * e2) for (e1, e2) in g]
        ctrl = rank_of(Mf, pts[:4000])
        out.append(dict(p=p, a=a, b=b, t=t, points_Fp2=len(pts), sets_tested=tested, dependent=dep,
                        control_frobenius_terms=len(Mf), control_frobenius_rank=ctrl))
        print(out[-1], flush=True)
    ok = all(r['dependent'] == 0 and r['control_frobenius_rank'] < r['control_frobenius_terms'] for r in out)
    print('ok', ok)
    if len(sys.argv) > 1:
        json.dump(dict(ok=ok, rows=out), open(sys.argv[1], 'w'), indent=1)

if __name__ == '__main__':
    main()
