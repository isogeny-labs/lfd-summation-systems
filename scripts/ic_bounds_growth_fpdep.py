"""N8-growth, Theorem 3.3 in characteristic p: tests of linear dependence over F_0 = F^p itself.

The older tests (ic_bounds_growth_charp.py, ..._circuits.py, lemmas.py sparse) decide dependence over the
CONSTANTS (ranks of evaluation matrices).  Theorem 3.3 is about circuits over F^p, whose coefficients may be
non-constant p-th powers.  This script decides F^p-dependence by the Wronskian criterion (Theorem 3.1 (a)):
for m <= p functions f_1..f_m, they are dependent over F^p = ker(d/dx1) iff det(D^(i) f_j)_{i<m} = 0 as a
function on the curve, where D^(i) is the i-th Hasse derivative w.r.t. x1 (D^i = i! D^(i), i < p).

Evaluation at a point Q = (alpha, beta) of C over F_{p^2} with dS_3/dx2 != 0: x1 - alpha = s is a uniformizer,
x2 = y(s) is the power-series root of S_3(alpha + s, y) = 0 (computed mod s^m and checked), and D^(i) f (Q) is
the coefficient of s^i of f(alpha + s, y(s)).  A nonzero determinant at ONE point proves F^p-independence.
A determinant vanishing at every sampled point is reported as "dependent" (evidence, not proof; used for the
controls, where dependence is known).

F^p-dependence of monomials depends only on the exponents modulo p (xi^(e + p u) = (xi^u)^p xi^e), so a scan
over the p^2 residues covers all monomials.

Tests
  controls : (1) the 9 monomials of S_3: dependent, every 8 of them independent;
             (2) twisted S_3: exponents e + p u_e with random u_e: dependent over F^p, although the rank over the
                 constants is 9 of 9 (this is the dependence the older tests cannot see);
             (3) a pair e, e + p u: dependent.
  three    : all three-term F^p-circuits {0, f1, f2} (exponent 0 without loss, by translation), exhaustively
             over all residues: Theorem 3.3 allows them only with pairwise ||2 da||_p + ||2 db||_p <= 12.
  far      : random sets B of m - 1 residues, pairwise ||2 da||_p + ||2 db||_p > 6(m-1)(m-2).  Theorem 3.3
             predicts that B is independent and that B + {f} is independent for EVERY residue f not in B
             (a circuit in B + {f} would contain two elements of B, or be a pair with equal residues).
             Every f of the p^2 residues is scanned.
usage: ic_bounds_growth_fpdep.py [out.json]
"""
import json, random, sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import s3_coeffs
from ic_bounds_growth_charp import Fq2, curve_points_q2, rank_q2

Z, ONE = (0, 0), (1, 0)

class Curve:
    def __init__(self, p, a, b, t, seed=1):
        self.p, self.K = p, Fq2(p)
        self.g = s3_coeffs(p, a, b, t)
        self.pts = curve_points_q2(self.K, a, b, t)
        self.rng = random.Random(seed)

    def smul(self, f, h, m):
        K = self.K; out = [Z] * m
        for i, u in enumerate(f):
            if u == Z: continue
            for j in range(m - i):
                if h[j] != Z: out[i + j] = K.add(out[i + j], K.mul(u, h[j]))
        return out

    def S_series(self, xs, ys, m):
        """S_3(xs, ys) for truncated series xs, ys."""
        K = self.K
        px = [[ONE] + [Z] * (m - 1), xs, self.smul(xs, xs, m)]
        py = [[ONE] + [Z] * (m - 1), ys, self.smul(ys, ys, m)]
        out = [Z] * m
        for (e1, e2), c in self.g.items():
            term = self.smul(px[e1], py[e2], m)
            out = [K.add(o, K.mul((c, 0), v)) for o, v in zip(out, term)]
        return out

    def local(self, m):
        """A random point with a usable local expansion: returns (alpha, beta, xs, ys)."""
        K = self.K
        while True:
            al, be = self.rng.choice(self.pts)
            if al == Z or be == Z: continue
            # dS/dx2 at the point
            d = Z
            for (e1, e2), c in self.g.items():
                if e2: d = K.add(d, K.mul((c * e2 % self.p, 0), K.mul(K.pw(al, e1), K.pw(be, e2 - 1))))
            if d == Z: continue
            dinv = K.inv(d)
            xs = [al, ONE] + [Z] * (m - 2) if m > 1 else [al]
            ys = [be] + [Z] * (m - 1)
            for n in range(1, m):
                r = self.S_series(xs, ys, m)[n]
                ys[n] = K.neg(K.mul(r, dinv))
            assert all(v == Z for v in self.S_series(xs, ys, m))
            return xs, ys

    def tables(self, xs, ys, m):
        """series of xs^a and ys^b for 0 <= a, b < p."""
        A = [[ONE] + [Z] * (m - 1)]; B = [[ONE] + [Z] * (m - 1)]
        for _ in range(self.p - 1):
            A.append(self.smul(A[-1], xs, m)); B.append(self.smul(B[-1], ys, m))
        return A, B

    def spow(self, f, e, m):
        r = [ONE] + [Z] * (m - 1)
        while e:
            if e & 1: r = self.smul(r, f, m)
            f = self.smul(f, f, m); e >>= 1
        return r

    def det(self, M):
        K = self.K; M = [list(r) for r in M]; n = len(M); d = ONE
        for c in range(n):
            s = next((i for i in range(c, n) if M[i][c] != Z), None)
            if s is None: return Z
            if s != c: M[c], M[s] = M[s], M[c]; d = K.neg(d)
            d = K.mul(d, M[c][c]); iv = K.inv(M[c][c])
            for i in range(c + 1, n):
                if M[i][c] != Z:
                    f = K.mul(M[i][c], iv)
                    M[i] = [K.add(v, K.neg(K.mul(f, w))) for v, w in zip(M[i], M[c])]
        return d

    def wronskian_nonzero(self, exps, points=6):
        """True if the Hasse Wronskian of the monomials (arbitrary exponents) is nonzero at some sampled point."""
        m = len(exps)
        for _ in range(points):
            xs, ys = self.local(m)
            cols = [self.smul(self.spow(xs, e[0], m), self.spow(ys, e[1], m), m) for e in exps]
            if self.det([[cols[j][i] for j in range(m)] for i in range(m)]) != Z: return True
        return False

def dist(n, p):
    n %= p
    return min(n, p - n)

def controls(C, rng):
    p = C.p; base = list(C.g)
    r = dict(S3_terms=len(base), S3_dependent=not C.wronskian_nonzero(base, 8))
    r['every_8_subset_independent'] = all(C.wronskian_nonzero(base[:i] + base[i + 1:], 8) for i in range(len(base)))
    tw = []; ranks = []
    for _ in range(4):
        M = [(e[0] + p * rng.randrange(4), e[1] + p * rng.randrange(4)) for e in base]
        tw.append(not C.wronskian_nonzero(M, 6))
        sample = rng.sample(C.pts, 80)
        ranks.append(rank_q2(C.K, [[C.K.mul(C.K.pw(P[0], e[0]), C.K.pw(P[1], e[1])) for e in M] for P in sample]))
    r['twisted_dependent_over_Fp'] = tw; r['twisted_rank_over_constants'] = ranks
    e = (rng.randrange(p), rng.randrange(p)); u = (rng.randrange(1, 4), rng.randrange(4))
    r['pair_e_e_plus_pu_dependent'] = not C.wronskian_nonzero([e, (e[0] + p * u[0], e[1] + p * u[1])], 8)
    r['pair_generic_independent'] = C.wronskian_nonzero([e, ((e[0] + 1) % p, e[1])], 8)
    return r

def three(C, points=5):
    """Residues f1 != f2, both != 0, such that {1, xi^f1, xi^f2} is F^p-dependent at every sampled point.
    With first column (1, 0, 0) the Wronskian is the 2x2 determinant of the vectors v(f) = (c_1(f), c_2(f))."""
    p, K = C.p, C.K
    res = [(a, b) for a in range(p) for b in range(p) if (a, b) != (0, 0)]
    groups = [res]
    for _ in range(points):
        xs, ys = C.local(3); A, B = C.tables(xs, ys, 3)
        v = {}
        for f in res:
            s = C.smul(A[f[0]], B[f[1]], 3)
            v[f] = (s[1], s[2])
        new = []
        for G in groups:
            cls = {}; wild = []
            for f in G:
                c1, c2 = v[f]
                if c1 == Z and c2 == Z: wild.append(f); continue
                key = ('inf',) if c1 == Z else K.mul(c2, K.inv(c1))
                cls.setdefault(key, []).append(f)
            for members in cls.values():
                if len(members) + len(wild) >= 2: new.append(members + wild)
            if len(wild) >= 2 and not cls: new.append(wild)
        groups = new
        if not groups: break
    survivors = sorted({f for G in groups for f in G})
    pairs = [(f1, f2) for G in groups for i, f1 in enumerate(G) for f2 in G[:i]]
    viol = [q for q in pairs if dist(2 * (q[0][0] - q[1][0]), p) + dist(2 * (q[0][1] - q[1][1]), p) > 12
            or dist(2 * q[0][0], p) + dist(2 * q[0][1], p) > 12 or dist(2 * q[1][0], p) + dist(2 * q[1][1], p) > 12]
    return dict(residues=len(res), surviving_groups=len(groups), surviving_residues=survivors[:20],
                dependent_triples=len(pairs), violations_of_theorem_3_3=len(viol))

def far(C, m, nsets, rng, points=4):
    p, K = C.p, C.K; thr = 6 * (m - 1) * (m - 2)
    done = bad_B = bad_ext = scanned = 0
    while done < nsets:
        B = []
        for _ in range(4000):
            f = (rng.randrange(p), rng.randrange(p))
            if all(dist(2 * (f[0] - h[0]), p) + dist(2 * (f[1] - h[1]), p) > thr for h in B): B.append(f)
            if len(B) == m - 1: break
        if len(B) < m - 1: continue
        done += 1
        cand = [(a, b) for a in range(p) for b in range(p) if (a, b) not in B]
        B_indep = False
        for _ in range(points):
            xs, ys = C.local(m); A, Bt = C.tables(xs, ys, m)
            cols = [C.smul(A[f[0]], Bt[f[1]], m) for f in B]
            # cofactors of the last column of the m x m matrix [cols | c(f)]
            cof = []
            for i in range(m):
                minor = [[cols[j][r] for j in range(m - 1)] for r in range(m) if r != i]
                d = C.det(minor)
                cof.append(d if (i + m - 1) % 2 == 0 else K.neg(d))
            if all(c == Z for c in cof): continue      # B itself degenerate at this point
            B_indep = True
            keep = []
            for f in cand:
                s = C.smul(A[f[0]], Bt[f[1]], m)
                w = Z
                for i in range(m):
                    if cof[i] != Z and s[i] != Z: w = K.add(w, K.mul(cof[i], s[i]))
                if w == Z: keep.append(f)
            scanned += len(cand); cand = keep
            if not cand: break
        if not B_indep: bad_B += 1
        elif cand: bad_ext += 1
    return dict(m=m, threshold=thr, sets=done, residues_scanned=scanned, sets_B_not_shown_independent=bad_B,
                sets_with_a_dependent_extension=bad_ext)

def main():
    rng = random.Random(5); out = []; t0 = time.time()
    C = Curve(131, 5, 9, 2)
    out.append(dict(test='controls', p=131, a=5, b=9, t=2, points_Fp2=len(C.pts), **controls(C, rng))); print(out[-1], flush=True)
    for (p, a, b, t) in [(61, 3, 7, 5), (131, 5, 9, 2), (103, 1, 0, 11)]:
        Cc = C if p == 131 else Curve(p, a, b, t)
        out.append(dict(test='three', p=p, a=a, b=b, t=t, **three(Cc))); print(out[-1], flush=True)
    out.append(dict(test='far', p=131, a=5, b=9, t=2, **far(C, 5, 40, rng))); print(out[-1], flush=True)
    C2 = Curve(311, 3, 7, 5)
    out.append(dict(test='far', p=311, a=3, b=7, t=5, **far(C2, 6, 6, rng))); print(out[-1], flush=True)
    c = out[0]
    ok = (c['S3_dependent'] and c['every_8_subset_independent'] and all(c['twisted_dependent_over_Fp'])
          and all(r == c['S3_terms'] for r in c['twisted_rank_over_constants']) and c['pair_e_e_plus_pu_dependent']
          and c['pair_generic_independent']
          and all(r['violations_of_theorem_3_3'] == 0 for r in out if r['test'] == 'three')
          and all(r['sets_B_not_shown_independent'] == 0 and r['sets_with_a_dependent_extension'] == 0 for r in out if r['test'] == 'far'))
    print('ok', ok, 'seconds', round(time.time() - t0))
    if len(sys.argv) > 1: json.dump(dict(ok=ok, rows=out), open(sys.argv[1], 'w'), indent=1)

if __name__ == '__main__':
    main()
