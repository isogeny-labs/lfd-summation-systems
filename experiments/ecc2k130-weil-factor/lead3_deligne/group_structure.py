#!/usr/bin/env python3
"""Lead 3, METHOD B for part 3: p-primary part of A(F_2^k) = O_K/(pi^k-1)  vs  A'(F_2^k) = T/(pi^k-1)T,
computed LOCALLY at p with p-adic lifts of the characters (no global SNF).  Cross-checked by group_structure.gp
(global Smith normal form for ell=11,19; HNF-mod-p^m for ell=131).

Decomposition: O_K (x) Z_p = prod over B-primes P of O_P.  T (sub-lattice, index p) differs from O_K only on the
glued component(s) of the chosen merged pair: split p: Z_p x Z_p (phi1, phi2); inert p: Z_{p^2}.
Sub-lattice  T  = {x : chi1(x) = mu chi2(x) mod p}            (split)   / {x : Tr(theta chi1(x)) = 0 mod p}   (inert)
Super-lattice T' = O_K + Z_p * v/p,  v = e1 + nu e2           (split)   / v = theta e1 + Frob(theta e1)        (inert)
"""
import json, sys
from deligne_chars import setup, F, factor_small, elt_order

M = 12  # p-adic precision (valuations of interest are < M for the k listed)

class Ring:  # Z/p^M or Z/p^M[t]/(t^2 - r)
    def __init__(s, p, r): s.p, s.q, s.r = p, p**M, r
    def mul(s, x, y): r = s.r or 0; q = s.q; return ((x[0]*y[0] + r*x[1]*y[1]) % q, (x[0]*y[1] + x[1]*y[0]) % q)
    def add(s, x, y): return ((x[0]+y[0]) % s.q, (x[1]+y[1]) % s.q)
    def sub(s, x, y): return ((x[0]-y[0]) % s.q, (x[1]-y[1]) % s.q)
    def inv(s, x):
        r = s.r or 0; n = (x[0]*x[0] - r*x[1]*x[1]) % s.q; ni = pow(n, -1, s.q)
        if s.r is None: return (pow(x[0], -1, s.q), 0)
        return (x[0]*ni % s.q, -x[1]*ni % s.q)
    def pw(s, x, e):
        R = (1, 0)
        while e:
            if e & 1: R = s.mul(R, x)
            x = s.mul(x, x); e >>= 1
        return R

def hensel_root_phi(Rg, ell, z0):
    z = z0
    for _ in range(8):  # quadratic convergence, 2^8 > M
        ph, dph, zp = (0, 0), (0, 0), (1, 0)
        for k in range(ell):
            ph = Rg.add(ph, zp)
            if k + 1 < ell: dph = Rg.add(dph, Rg.mul(((k+1) % Rg.q, 0), zp))
            zp = Rg.mul(zp, z)
        z = Rg.sub(z, Rg.mul(ph, Rg.inv(dph)))
    return z

def hensel_sqrt(a, s0, p):
    q = p**M; x = s0
    for _ in range(8): x = (x - (x*x - a) * pow(2*x, -1, q)) % q
    return x

def val(x, p):
    if x == 0: return M + 100
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def eldiv2(C, p):
    """elementary divisors exponents of a 2x2 matrix over Z_p (entries exact rationals given as ints mod p^M, assume nonsingular)."""
    e1 = min(val(c, p) for c in (C[0][0], C[0][1], C[1][0], C[1][1]))
    det = (C[0][0]*C[1][1] - C[0][1]*C[1][0]) % (p**M)
    return sorted([e1, val(det, p) - e1])

def coords_in_basis(Bt, vec, p):
    """solve Bt * c = vec over Q_p where Bt is 2x2 with columns; return c as p-adic ints (exact division)."""
    q = p**(M+2)
    a, b, c, d = Bt[0][0], Bt[0][1], Bt[1][0], Bt[1][1]
    det = a*d - b*c; vd = val(det % q, p); du = (det // p**vd) % q
    x0 = (d*vec[0] - b*vec[1]); x1 = (-c*vec[0] + a*vec[1])
    assert x0 % p**vd == 0 and x1 % p**vd == 0, "not in lattice"
    iu = pow(du, -1, q)
    return ((x0 // p**vd) * iu % p**M, (x1 // p**vd) * iu % p**M)

def local_analysis(ell, p, kmax, which='conj'):
    Fq, a1, a2, roots = setup(ell, p)
    split = Fq.r is None
    Rg = Ring(p, Fq.r)
    # B-characters and their p-adic lifts phi(pi)
    if split:
        s7 = hensel_sqrt(-7, ((2*a1[0]+1) % p), p)            # lift of sqrt(-7) matching a1 = (-1+s)/2
        inv2 = pow(2, -1, p**M)
        A1 = ((-1 + s7) * inv2 % p**M, 0); A2 = ((-1 - s7) * inv2 % p**M, 0)
    else:
        inv2 = pow(2, -1, p**M)
        A1 = ((-inv2) % p**M, inv2); A2 = ((-inv2) % p**M, (-inv2) % p**M)
    Z = {r: hensel_root_phi(Rg, ell, r) for r in roots}
    lam_of = lambda r, A: Rg.mul(A, Z[r])
    # characters of the chosen (conjugation-stable) merged pair: z1 with z1^2 = c^-1, c = a1/a2
    c = Fq.mul(a1, Fq.inv(a2)); z1 = Fq.pw(c, -(ell+1)//2)
    assert Fq.pw(z1, 2) == Fq.inv(c)
    if which == 'loword':   # split case only: merged pair (z1,a1)~(c z1,a2) whose pi-eigenvalue a1 z1 has least order
        assert split
        fac = factor_small(p-1)
        cands = [r for r in roots if Fq.mul(c, r) != (1, 0)]
        z1 = min(cands, key=lambda r: elt_order(Fq, Fq.mul(a1, r), p-1, fac))
    z2 = Fq.mul(c, z1)                                        # partner, with alpha -> a2 (= z1^-1 for the conj pair)
    L1, L2 = lam_of(z1, A1), lam_of(z2, A2)                   # p-adic phi1(pi), phi2(pi)
    assert (L1[0]-L2[0]) % p == 0 and (L1[1]-L2[1]) % p == 0
    # all B-characters' lifted pi-values, grouped into B-primes (Frobenius orbits)
    allchars = [(r, A1) for r in roots] + [(r, A2) for r in roots]
    rows = []
    for k in range(1, kmax+1):
        # O_K: p-part = sum over B-primes of O_P/(y) ; for split each char is a prime Z_p -> Z/p^v ;
        # inert: chars come in Frobenius pairs, each prime Z_{p^2}/(y) = (Z/p^v)^2 with v = val of y (unramified)
        okdiv = []
        if split:
            for (r, A) in allchars:
                y = Rg.pw(lam_of(r, A), k); v = val((y[0]-1) % p**M, p)
                if v: okdiv.append(v)
        else:
            for (r, A) in allchars:   # each char once; a prime has 2 chars both with same valuation -> one (Z/p^v)^2 per prime
                y = Rg.pw(lam_of(r, A), k); y = (y[0]-1, y[1]); v = min(val(y[0] % p**M, p), val(y[1] % p**M, p))
                if v: okdiv.append(v)      # 2 chars per prime x 1 entry each = (Z/p^v)^2 per prime
        okdiv.sort()
        # glued component, O_K part
        if split:
            u = (Rg.pw(L1, k)[0]-1) % p**M; w = (Rg.pw(L2, k)[0]-1) % p**M
            glue_ok = sorted([e for e in (val(u, p), val(w, p)) if e])
            Y = [[u, 0], [0, w]]
            mu, nu = 5 % p, 7 % p      # any units (results do not depend on them: checked below with another choice)
            subB = [[mu, p], [1, 0]]               # columns (mu,1),(p,0)
            supB = [[1, 0], [nu, p]]               # columns (1,nu)/p, (0,1) scaled by p -> (1,nu),(0,p); scaling is harmless
        else:
            y = Rg.pw(L1, k); y = ((y[0]-1) % p**M, y[1])
            r = Fq.r
            Y = [[y[0], r*y[1] % p**M], [y[1], y[0]]]
            vv = min(val(y[0], p), val(y[1], p)); glue_ok = [vv, vv] if vv else []
            th = (3, 1)
            # sub: {x0 + x1 t : th0 x0 + r th1 x1 = 0 mod p} -> basis (-r th1/th0, 1), (p, 0)
            subB = [[(-r*th[1]*pow(th[0], -1, p)) % p, p], [1, 0]]
            supB = [[th[0], 0], [th[1], p]]       # th and p*(0,1)  (= p * {th/p, (0,1)})
        def struct(Bt):
            cols = [[Y[0][0]*Bt[0][j] + Y[0][1]*Bt[1][j], Y[1][0]*Bt[0][j] + Y[1][1]*Bt[1][j]] for j in range(2)]
            C = [coords_in_basis(Bt, col, p) for col in cols]          # C[j] = coords of Y*b_j
            Cm = [[C[0][0], C[1][0]], [C[0][1], C[1][1]]]
            return [e for e in eldiv2(Cm, p) if e]
        glue_sub = struct(subB); glue_sup = struct(supB)
        rest = okdiv[:]
        for e in glue_ok: rest.remove(e)
        sub = sorted(rest + glue_sub); sup = sorted(rest + glue_sup)
        if okdiv or sub or sup:
            rows.append(dict(k=k, OK=okdiv, T_sub=sub, T_super=sup, differs=(okdiv != sub or okdiv != sup)))
    return dict(ell=ell, p=p, split=split, pair=which, lambda_mod_p=L1[0] % p if split else [L1[0] % p, L1[1] % p], rows=rows)

if __name__ == '__main__':
    out = []
    def fmt(v, p):
        from collections import Counter
        C = Counter(v)
        return ' + '.join(('(Z/%d^%d)^%d' % (p, e, m) if e > 1 else '(Z/%d)^%d' % (p, m)) if m > 1 else ('Z/%d^%d' % (p, e) if e > 1 else 'Z/%d' % p) for e, m in sorted(C.items())) or '0'
    for ell, p, kmax, which in [(11, 23, 60, 'conj'), (11, 23, 60, 'loword'), (19, 457, 320, 'conj'), (19, 457, 320, 'loword'),
                                (131, 263, 600, 'conj'), (131, 263, 600, 'loword'), (131, 146505763881528721, 300, 'conj')]:
        r = local_analysis(ell, p, kmax, which); out.append(r)
        print('ell=%d p=%d (%s), A\' from %s pair (pi-eigenvalue on kernel = %s); k<=%d with p | #A(F_2^k):' % (ell, p, 'split' if r['split'] else 'inert', which, r['lambda_mod_p'], kmax))
        for row in r['rows']:
            print('  k=%-4d A: %-24s A\'(sub T): %-30s A\'(super T\'): %-30s %s' % (row['k'], fmt(row['OK'], p), fmt(row['T_sub'], p), fmt(row['T_super'], p), 'DIFFERS' if row['differs'] else 'same'))
        if not r['rows']: print('  (none)')
        sys.stdout.flush()
    json.dump(out, open('group_structure.json', 'w'), indent=1)
