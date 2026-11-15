#!/usr/bin/env python3
"""Lane D small exact check (independent of repository engines).

For tiny (m, k, p) build the ORDINARY-ring static Macaulay space
    Mac_D = span{ mon * f : f in chains + [ghat], deg(mon*f) <= D }
in F_p[u_{i,j}] (i < m, j < k), by sparse Gaussian elimination mod p with
pivots on the highest monomial (degree, then lex), and find
    sd = least D with 1 in Mac_D.
Compare with  deg(ghat) + maxwt(supp g^{-1}),  g = NF(ghat) in A = F_p[(Z/N)^m],
where g^{-1} is computed twice: (i) by solving g*h = 1 with the multiplication
matrix of A, (ii) by Fourier inversion over H^m, H = mu_N in F_p.  Also checks
the top-coefficient character sum  c_top = N^{-m} sum_alpha prod(alpha_i)/g(alpha).

Usage: python3 check_small.py  (prints one line per instance and a summary)
Single thread, pure Python; k <= 3.
"""
import itertools, random, sys

def egcd_inv(a, p):
    return pow(a % p, p - 2, p)

def wt(e):
    return sum(bin(x).count("1") for x in e)

# ---------- polynomials as dict {exponent tuple: coeff} ----------
def padd(P, Q, p, s=1):
    R = dict(P)
    for e, c in Q.items():
        R[e] = (R.get(e, 0) + s * c) % p
        if R[e] == 0:
            del R[e]
    return R

def pmul(P, Q, p):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            R[e] = (R.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in R.items() if c}

def deg(P):
    return max((sum(e) for e in P), default=-1)

# ---------- setting ----------
class Inst:
    def __init__(s, m, k, p):
        assert (p - 1) % (2 ** k) == 0
        s.m, s.k, s.p, s.N = m, k, p, 2 ** k
        s.r = m * k
        # a generator of mu_N in F_p
        for z in range(2, p):
            if pow(z, s.N, p) == 1 and all(pow(z, s.N // q, p) != 1 for q in [2]):
                s.zeta = z
                break
        s.H = [pow(s.zeta, j, p) for j in range(s.N)]

    def var(s, i, j):
        e = [0] * s.r
        e[i * s.k + j] = 1
        return tuple(e)

    def chains(s):
        p, out = s.p, []
        one = tuple([0] * s.r)
        for i in range(s.m):
            for j in range(s.k):
                sq = list(s.var(i, j)); sq[i * s.k + j] = 2; sq = tuple(sq)
                nxt = s.var(i, j + 1) if j < s.k - 1 else one
                out.append({sq: 1, nxt: p - 1})
        return out

    def x_power(s, i, n):
        """ordinary reduced monomial u-representative of x_i^n, 0<=n<N"""
        e = [0] * s.r
        for j in range(s.k):
            if (n >> j) & 1:
                e[i * s.k + j] = 1
        return tuple(e)

    # NF: ordinary monomial -> group element e in (Z/N)^m
    def nf_mon(s, e):
        return tuple(sum(e[i * s.k + j] * 2 ** j for j in range(s.k)) % s.N for i in range(s.m))

    def NF(s, P):
        R = {}
        for e, c in P.items():
            g = s.nf_mon(e)
            R[g] = (R.get(g, 0) + c) % s.p
        return {g: c for g, c in R.items() if c}

    def lift(s, h):  # element of A -> reduced ordinary polynomial
        out = {}
        for g, c in h.items():
            e = [0] * s.r
            for i in range(s.m):
                for j in range(s.k):
                    e[i * s.k + j] = (g[i] >> j) & 1
            out[tuple(e)] = c
        return out

    # ---------- group algebra ----------
    def elems(s):
        return list(itertools.product(range(s.N), repeat=s.m))

    def amul(s, a, b):
        R = {}
        for g1, c1 in a.items():
            for g2, c2 in b.items():
                g = tuple((x + y) % s.N for x, y in zip(g1, g2))
                R[g] = (R.get(g, 0) + c1 * c2) % s.p
        return {g: c for g, c in R.items() if c}

    def ev(s, a, alpha):
        p, t = s.p, 0
        for g, c in a.items():
            v = c
            for ai, gi in zip(alpha, g):
                v = v * pow(ai, gi, p) % p
            t = (t + v) % p
        return t

    def inv_linear(s, g):
        """solve g*h = 1 by Gaussian elimination on the multiplication matrix"""
        p, E = s.p, s.elems()
        idx = {e: n for n, e in enumerate(E)}
        n = len(E)
        # column c = g * x^{E[c]}
        M = [[0] * (n + 1) for _ in range(n)]
        for c, e in enumerate(E):
            prod = s.amul(g, {e: 1})
            for gg, v in prod.items():
                M[idx[gg]][c] = v
        M[idx[tuple([0] * s.m)]][n] = 1
        row = 0
        piv = []
        for c in range(n):
            pr = next((r for r in range(row, n) if M[r][c]), None)
            if pr is None:
                return None
            M[row], M[pr] = M[pr], M[row]
            iv = egcd_inv(M[row][c], p)
            M[row] = [v * iv % p for v in M[row]]
            for r in range(n):
                if r != row and M[r][c]:
                    f = M[r][c]
                    M[r] = [(v - f * w) % p for v, w in zip(M[r], M[row])]
            row += 1
        return {E[c]: M[c][n] for c in range(n) if M[c][n]}

    def inv_fourier(s, g):
        p, N, m = s.p, s.N, s.m
        Hm = list(itertools.product(s.H, repeat=m))
        vals = {}
        for al in Hm:
            v = s.ev(g, al)
            if v == 0:
                return None
            vals[al] = egcd_inv(v, p)
        Nm_inv = egcd_inv(pow(N, m, p), p)
        out = {}
        for e in s.elems():
            t = 0
            for al in Hm:
                w = vals[al]
                for ai, ei in zip(al, e):
                    w = w * pow(ai, (N - ei) % N, p) % p
                t = (t + w) % p
            t = t * Nm_inv % p
            if t:
                out[e] = t
        return out

    def ctop_charsum(s, g):
        p, N, m = s.p, s.N, s.m
        t = 0
        for al in itertools.product(s.H, repeat=m):
            v = s.ev(g, al)
            pr = 1
            for ai in al:
                pr = pr * ai % p
            t = (t + pr * egcd_inv(v, p)) % p
        return t * egcd_inv(pow(N, m, p), p) % p

    # ---------- ordinary-ring Macaulay ----------
    def monomials_upto(s, d):
        out = []
        for tot in range(d + 1):
            for c in itertools.combinations_with_replacement(range(s.r), tot):
                e = [0] * s.r
                for v in c:
                    e[v] += 1
                out.append(tuple(e))
        return out

    def sd_bruteforce(s, ghat, Dmax):
        """least D <= Dmax with 1 in Mac_D (None if none). Incremental in D."""
        p = s.p
        gens = s.chains() + [ghat]
        key = lambda e: (-sum(e), tuple(-x for x in e))  # highest monomial first
        basis = {}  # pivot monomial -> row dict (leading coeff 1)
        one = tuple([0] * s.r)
        added = set()
        for D in range(0, Dmax + 1):
            for fi, f in enumerate(gens):
                df = deg(f)
                if df > D:
                    continue
                for mon in s.monomials_upto(D - df):
                    if (fi, mon) in added:
                        continue
                    added.add((fi, mon))
                    row = pmul({mon: 1}, f, p)
                    while row:
                        lead = min(row, key=key)
                        if lead in basis:
                            c = row[lead]
                            row = padd(row, basis[lead], p, -c)
                        else:
                            iv = egcd_inv(row[lead], p)
                            basis[lead] = {e: v * iv % p for e, v in row.items()}
                            break
            if one in basis:
                return D
        return None

# ---------- generators ----------
def s3_normalized(I, a, b, t):
    """chain normal-form representative (reduced, degree 2) of S_3(x1,x2,t)."""
    p = I.p
    # S3 as dict over (n1, n2) exponents of x1, x2
    S = {}
    def add(n1, n2, c):
        S[(n1, n2)] = (S.get((n1, n2), 0) + c) % p
    # (x1-x2)^2 t^2
    add(2, 0, t * t); add(1, 1, -2 * t * t); add(0, 2, t * t)
    # -2((x1+x2)(x1x2+a)+2b) t
    add(2, 1, -2 * t); add(1, 2, -2 * t); add(1, 0, -2 * a * t); add(0, 1, -2 * a * t); add(0, 0, -4 * b * t)
    # (x1x2-a)^2 - 4b(x1+x2)
    add(2, 2, 1); add(1, 1, -2 * a); add(0, 0, a * a); add(1, 0, -4 * b); add(0, 1, -4 * b)
    P = {}
    for (n1, n2), c in S.items():
        if c % p == 0:
            continue
        e = tuple(x + y for x, y in zip(I.x_power(0, n1 % I.N), I.x_power(1, n2 % I.N)))
        P[e] = (P.get(e, 0) + c) % p
    return {e: c for e, c in P.items() if c}

def s3_raw(I, a, b, t):
    """unreduced S_3(x1,x2,t) in x_i = u_{i,0}: degree 4 ordinary polynomial."""
    p = I.p
    x1, x2 = {I.var(0, 0): 1}, {I.var(1, 0): 1}
    one = {tuple([0] * I.r): 1}
    def sc(P, c): return {e: v * c % p for e, v in P.items() if v * c % p}
    d = padd(x1, x2, p, -1)
    T1 = sc(pmul(d, d, p), t * t)
    s12 = padd(x1, x2, p); pr = pmul(x1, x2, p)
    T2 = sc(padd(pmul(s12, padd(pr, sc(one, a), p), p), sc(one, 2 * b), p), -2 * t)
    q = padd(pr, sc(one, a), p, -1)
    T3 = padd(pmul(q, q, p), sc(s12, 4 * b), p, -1)
    return padd(padd(T1, T2, p), T3, p)

def random_reduced(I, delta, rng, density=0.6):
    """random reduced ordinary polynomial of degree exactly delta."""
    p = I.p
    mons = [e for e in I.monomials_upto(delta) if max(e, default=0) <= 1]
    P = {}
    top = [e for e in mons if sum(e) == delta]
    P[rng.choice(top)] = rng.randrange(1, p)
    for e in mons:
        if rng.random() < density:
            P[e] = rng.randrange(1, p)
    return {e: c for e, c in P.items() if c}

def run_case(I, ghat, label):
    g = I.NF(ghat)
    dhat = deg(ghat)
    hl = I.inv_linear(g)
    hf = I.inv_fourier(g)
    Dmax = dhat + I.m * I.k
    sd = I.sd_bruteforce(ghat, Dmax + 1)
    if hl is None:
        ok = (sd is None) and hf is None
        return dict(label=label, unit=False, sd=sd, pred=None, ok=ok)
    assert hl == hf, "linear and Fourier inverses differ"
    mw = max(wt(e) for e in hl)
    top = tuple([I.N - 1] * I.m)
    ct = hl.get(top, 0)
    cs = I.ctop_charsum(g)
    pred = dhat + mw
    ok = (sd == pred) and (ct == cs) and ((pred == dhat + I.m * I.k) == (ct != 0))
    return dict(label=label, unit=True, sd=sd, pred=pred, maxwt=mw, ctop=ct, charsum=cs, ok=ok)

def main():
    rng = random.Random(20261115)
    results = []
    plans = [  # (m, k, p, n_s3, n_raw, n_random)
        (2, 1, 3, 9, 4, 6),
        (2, 1, 7, 12, 4, 6),
        (2, 2, 5, 25, 6, 10),
        (2, 2, 13, 20, 4, 10),
        (2, 2, 17, 12, 0, 6),
        (2, 3, 17, 4, 0, 2),
        (2, 3, 41, 10, 2, 3),
        (2, 3, 73, 8, 0, 2),
        (3, 1, 5, 0, 0, 6),
    ]
    for (m, k, p, ns3, nraw, nrand) in plans:
        I = Inst(m, k, p)
        cases = []
        if m == 2:
            for _ in range(ns3):
                a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                cases.append((s3_normalized(I, a, b, t), f"S3norm a={a} b={b} t={t}"))
            for _ in range(nraw):
                a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                cases.append((s3_raw(I, a, b, t), f"S3raw a={a} b={b} t={t}"))
        for _ in range(nrand):
            delta = rng.choice([1, 2, 3]) if m * k >= 3 else rng.choice([1, 2])
            cases.append((random_reduced(I, delta, rng), f"random delta={delta}"))
        # deliberate low-weight-inverse case: ghat = u_{0,k-1} (g^2 = 1)
        cases.append(({I.var(0, k - 1): 1}, "ghat=u_{1,k-1}"))
        for ghat, lab in cases:
            r = run_case(I, ghat, lab)
            r.update(m=m, k=k, p=p)
            results.append(r)
            print(r, flush=True)
    units = [r for r in results if r["unit"]]
    bad = [r for r in results if not r["ok"]]
    print(f"SUMMARY cases={len(results)} unit={len(units)} nonunit={len(results)-len(units)} "
          f"mismatches={len(bad)}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
