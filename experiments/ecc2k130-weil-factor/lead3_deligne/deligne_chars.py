#!/usr/bin/env python3
"""Lead 3, METHOD B (independent of PARI): explicit Deligne-module construction via characters.

O_K = Z[z]/Phi_ell (x) Z[w]/(w^2-w+2); vector coords x[i + d*j] <-> z^i w^j (d = ell-1).
alpha = w-1, beta = -w, pi = alpha z, pibar = beta z^-1  (pi*pibar = 2).
For a prime p | U_ell, characters of B = O_K/p into Fbar_p are chi_(z0,a): z -> z0 (primitive ell-th root),
alpha -> a (a root of x^2+x+2).  chi(pi) = a z0, chi(pibar) = a' z0^-1 (a' = 2/a).
Everything claimed by character theory is then RE-CHECKED by plain integer arithmetic in O_K:
the lattices are written down by explicit Z-bases and pi/pibar/zeta are applied by polynomial multiplication.
Usage: python3 deligne_chars.py  (runs ell=11/p=23, ell=19/p=457, ell=131/p=263, ell=131/p=146505763881528721,
and controls). Output -> deligne_chars.json
"""
import json, random, sys
random.seed(20261127)

# ---------- F_p or F_p^2 = F_p[t]/(t^2 - r) ----------
class F:
    def __init__(s, p, r=None): s.p, s.r = p, r
    def e(s, a, b=0): return (a % s.p, b % s.p)
    def add(s, x, y): return ((x[0]+y[0]) % s.p, (x[1]+y[1]) % s.p)
    def sub(s, x, y): return ((x[0]-y[0]) % s.p, (x[1]-y[1]) % s.p)
    def mul(s, x, y):
        p, r = s.p, (s.r or 0)
        return ((x[0]*y[0] + r*x[1]*y[1]) % p, (x[0]*y[1] + x[1]*y[0]) % p)
    def pw(s, x, n):
        if n < 0: x, n = s.inv(x), -n
        R = (1, 0)
        while n:
            if n & 1: R = s.mul(R, x)
            x = s.mul(x, x); n >>= 1
        return R
    def conj(s, x): return (x[0], (-x[1]) % s.p)
    def norm(s, x): return (x[0]*x[0] - (s.r or 0)*x[1]*x[1]) % s.p
    def inv(s, x):
        if s.r is None: return (pow(x[0], s.p-2, s.p), 0)
        n = pow(s.norm(x), s.p-2, s.p); c = s.conj(x); return (c[0]*n % s.p, c[1]*n % s.p)
    def order_size(s): return s.p - 1 if s.r is None else s.p*s.p - 1

def sqrt_mod(a, p):  # Tonelli-Shanks
    a %= p
    if a == 0: return 0
    assert pow(a, (p-1)//2, p) == 1
    q, s_ = p-1, 0
    while q % 2 == 0: q //= 2; s_ += 1
    zz = 2
    while pow(zz, (p-1)//2, p) != p-1: zz += 1
    m, c, t, r = s_, pow(zz, q, p), pow(a, q, p), pow(a, (q+1)//2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2*t2 % p; i += 1
        b = pow(c, 1 << (m-i-1), p); m, c, t, r = i, b*b % p, t*b*b % p, r*b % p
    return r

def is_prime(n):
    if n < 2: return False
    for q in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n % q == 0: return n == q
    d, s_ = n-1, 0
    while d % 2 == 0: d //= 2; s_ += 1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s_-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def factor_small(n):
    f, q = {}, 2
    while q*q <= n and q < 10**7:
        while n % q == 0: f[q] = f.get(q, 0)+1; n //= q
        q += 1
    if n > 1:
        assert is_prime(n), n
        f[n] = f.get(n, 0)+1
    return f

def elt_order(Fq, x, N, fac):
    o = N
    for q in fac:
        while o % q == 0 and Fq.pw(x, o//q) == (1, 0): o //= q
    return o

# ---------- integer arithmetic in O_K ----------
def mul_OK(x, y, ell):
    d = ell-1
    def pm(a, b):
        r = [0]*(2*d)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj: r[i+j] += ai*bj
        # reduce mod Phi_ell: z^ell = 1, then z^(ell-1) = -(1+...+z^(ell-2))
        rr = [0]*ell
        for k, c in enumerate(r): rr[k % ell] += c
        top = rr[ell-1]
        return [rr[i]-top for i in range(d)]
    x0, x1, y0, y1 = x[:d], x[d:], y[:d], y[d:]
    t0, t2 = pm(x0, y0), pm(x1, y1)
    t1 = [u+v for u, v in zip(pm(x0, y1), pm(x1, y0))]
    return [u - 2*v for u, v in zip(t0, t2)] + [u + v for u, v in zip(t1, t2)]   # w^2 = w - 2

def unit(ell, k):
    v = [0]*(2*(ell-1)); v[k] = 1; return v

def gens(ell):
    d = ell-1
    z = unit(ell, 1); w = unit(ell, d); one = unit(ell, 0)
    alpha = [a-b for a, b in zip(w, one)]; beta = [-a for a in w]
    zinv = [0]*(2*d)                       # z^-1 = z^(ell-1) = -(1+z+...+z^(ell-2))
    for i in range(d): zinv[i] = -1
    return dict(z=z, w=w, one=one, pi=mul_OK(alpha, z, ell), pib=mul_OK(beta, zinv, ell), zinv=zinv)

# ---------- characters ----------
def setup(ell, p):
    d = ell-1
    r = None if pow(-7 % p, (p-1)//2, p) == 1 else -7
    Fq = F(p, r)
    s7 = (sqrt_mod(-7, p), 0) if r is None else (0, 1)
    half = pow(2, p-2, p)
    a1 = Fq.mul(Fq.sub(s7, (1, 0)), (half, 0))         # alpha -> (-1+sqrt-7)/2
    a2 = Fq.sub(Fq.e(-1), a1)
    N = Fq.order_size(); assert N % ell == 0, "mu_ell not in residue field"
    g = (2, 0) if r is None else (1, 1)
    while True:
        zeta0 = Fq.pw(g, N//ell)
        if zeta0 != (1, 0): break
        g = (g[0]+1, g[1])
    roots = [Fq.pw(zeta0, k) for k in range(1, ell)]
    return Fq, a1, a2, roots

def chi(Fq, x, ell, z0, a):
    d = ell-1; wv = Fq.add(a, (1, 0))
    acc0 = (0, 0); zp = (1, 0)
    acc1 = (0, 0)
    for i in range(d):
        acc0 = Fq.add(acc0, Fq.mul((x[i] % Fq.p, 0), zp))
        acc1 = Fq.add(acc1, Fq.mul((x[d+i] % Fq.p, 0), zp))
        zp = Fq.mul(zp, z0)
    return Fq.add(acc0, Fq.mul(acc1, wv))

def idempotent(Fq, ell, z0, a, aother):
    """coords in Fq of the primitive idempotent e_(z0,a) of B (x) Fq: Lagrange in z times Lagrange in w."""
    d = ell-1
    # E(z) = Phi(z)/((z-z0) Phi'(z0)); Phi(z)/(z-z0) = sum_{k} q_k z^k by synthetic division of 1+z+..+z^(ell-1)
    q = [(0, 0)]*d; acc = (0, 0)
    coeffs = [(1, 0)]*ell   # Phi coefficients high->low all 1
    for k in range(ell-1):  # quotient degree ell-2 ... 0
        acc = Fq.add(Fq.mul(acc, z0), coeffs[k]); q[ell-2-k] = acc
    dphi = (0, 0); zp = (1, 0)
    for k in range(1, ell):
        dphi = Fq.add(dphi, Fq.mul((k, 0), zp)); zp = Fq.mul(zp, z0)
    inv = Fq.inv(dphi); Ez = [Fq.mul(c, inv) for c in q]
    s, s2 = Fq.add(a, (1, 0)), Fq.add(aother, (1, 0))     # images of w
    den = Fq.inv(Fq.sub(s, s2))
    c0 = Fq.mul(Fq.sub((0, 0), s2), den); c1 = den       # (w - s2)/(s - s2)
    return [Fq.mul(e, c0) for e in Ez] + [Fq.mul(e, c1) for e in Ez]

def run(ell, p, full_checks=True):
    d, n = ell-1, 2*(ell-1)
    G = gens(ell)
    Fq, a1, a2, roots = setup(ell, p)
    out = dict(ell=ell, p=p, split_in_Q_sqrt_m7=(Fq.r is None), p_mod_ell=p % ell)
    c = Fq.mul(a1, Fq.inv(a2))
    out['c=a1/a2 has order ell'] = (Fq.pw(c, ell) == (1, 0) and c != (1, 0))
    chars = [(z0, a1) for z0 in roots] + [(z0, a2) for z0 in roots]
    # S-characters: (chi(pi), chi(pibar))
    def sval(ch):
        z0, a = ch; return (Fq.mul(a, z0), Fq.mul(Fq.mul((2, 0), Fq.inv(a)), Fq.inv(z0)))
    groups = {}
    for ch in chars: groups.setdefault(sval(ch), []).append(ch)
    merged = [g for g in groups.values() if len(g) == 2]
    out['#B-characters'] = len(chars); out['#S-characters'] = len(groups)
    out['#merged pairs (S-char with 2 B-chars)'] = len(merged)
    out['max fibre'] = max(len(g) for g in groups.values())
    # Frobenius x->x^p on characters: (z0,a) -> (z0^p, a^p)
    frob = lambda ch: (Fq.pw(ch[0], p), Fq.pw(ch[1], p))
    rat_merged = [g for g in merged if frob(g[0]) in g]   # S-character with values in F_p
    out['#merged pairs whose S-character is F_p-rational'] = len(rat_merged)
    out['#merged pairs F_p-rational with the two B-chars Frobenius-conjugate (inert type)'] = sum(1 for g in rat_merged if frob(g[0]) == g[1] and g[0] != g[1])
    # conj-stable pair: chi2 = chi1 o conj, conj: z->z^-1, alpha->beta
    cconj = lambda ch: (Fq.inv(ch[0]), Fq.sub(Fq.e(-1), ch[1]))
    conjstab = [g for g in merged if cconj(g[0]) in g]
    out['#merged pairs stable under complex conjugation'] = len(conjstab)
    g = conjstab[0]; ch1, ch2 = (g[0], g[1]) if g[0][1] == a1 else (g[1], g[0])
    lam = sval(ch1)[0]
    out['chosen pair: chi1(pi) = chi2(pi) = lambda'] = lam
    out['lambda^2 == 2'] = Fq.mul(lam, lam) == (2, 0)
    out['lambda in F_p'] = lam[1] == 0
    out['chosen pair is the F_p-rational one'] = g in rat_merged
    fac = factor_small(p-1)
    out['ord(lambda) in F_p^* (field of definition of kernel points: F_2^ord)'] = elt_order(Fq, lam, p-1, fac) if lam[1] == 0 else None
    # orders of pi-eigenvalues over all merged pairs which are F_p-rational
    ords = {}
    for gg in merged:
        lv = sval(gg[0])[0]
        if lv[1] == 0:
            o = elt_order(Fq, lv, p-1, fac); ords[o] = ords.get(o, 0)+1
    out['histogram ord(chi(pi)) over F_p-rational merged pairs'] = {str(k): v for k, v in sorted(ords.items())}

    # ---------- explicit SUB-lattice T = {x in O_K : f(x) = 0 mod p} ----------
    if Fq.r is None:
        mu = random.randrange(1, p)
        fcoef = [Fq.sub(chi(Fq, unit(ell, k), ell, *ch1), Fq.mul((mu, 0), chi(Fq, unit(ell, k), ell, *ch2)))[0] for k in range(n)]
        out['functional'] = 'f = chi1 - mu*chi2, mu=%d' % mu
    else:
        theta = (random.randrange(p), random.randrange(1, p))
        def tr(x): y = Fq.mul(theta, x); return (2*y[0]) % p
        fcoef = [tr(chi(Fq, unit(ell, k), ell, *ch1)) for k in range(n)]
        out['functional'] = 'f = Tr_{F_p^2/F_p}(theta*chi1), theta=%s' % (theta,)
    fval = lambda x: sum(a*b for a, b in zip(fcoef, x)) % p
    j0 = next(k for k in range(n) if fcoef[k]); ij = pow(fcoef[j0], p-2, p)
    Tbasis = []
    for k in range(n):
        if k == j0: v = [0]*n; v[j0] = p
        else: v = unit(ell, k); v[j0] = -(fcoef[k]*ij % p)
        Tbasis.append(v)
    # index = |det| = p (basis is triangular-up-to-permutation with one entry p)
    assert all(fval(v) == 0 for v in Tbasis)
    out['T: all basis vectors satisfy f=0'] = True
    out['T: index [O_K:T]'] = p
    pi_ok = all(fval(mul_OK(G['pi'], v, ell)) == 0 for v in Tbasis)
    pib_ok = all(fval(mul_OK(G['pib'], v, ell)) == 0 for v in Tbasis)
    zeta_bad = [k for k, v in enumerate(Tbasis) if fval(mul_OK(G['z'], v, ell)) != 0]
    out['T: pi*T in T (integer check on every basis vector)'] = pi_ok
    out['T: pibar*T in T'] = pib_ok
    out['T: #basis vectors t with zeta*t notin T'] = len(zeta_bad)
    # End(T) claim: E = {x : chi1(x) = chi2(x)} (split) or {x: chi1(x) in F_p} (inert); check zeta^k and random samples
    def in_E(x):
        v1 = chi(Fq, x, ell, *ch1)
        return v1 == chi(Fq, x, ell, *ch2) if Fq.r is None else v1[1] == 0
    zk = G['one']; zin = []
    for k in range(1, 2*ell+1):
        zk = mul_OK(zk, G['z'], ell)
        if in_E(zk): zin.append(k)
    out['k in 1..2ell with zeta^k in claimed End(T)'] = zin
    out['-1 in End(T)'] = True
    # sampled verification x in E => x T in T, x notin E => x T not in T (integer multiplication)
    okE, okNotE = 0, 0
    for trial in range(6 if ell > 100 else 20):
        x = [random.randrange(-5, 6) for _ in range(n)]
        t = random.choice(Tbasis)
        if trial % 2 == 0:   # project x into the claimed E (codim-1 condition mod p) by fixing one coordinate
            def gfun(y):
                v1 = chi(Fq, y, ell, *ch1)
                return Fq.sub(v1, chi(Fq, y, ell, *ch2))[0] if Fq.r is None else v1[1]
            jj = next(k for k in range(n) if gfun(unit(ell, k)) % p)
            x[jj] -= gfun(x) * pow(gfun(unit(ell, jj)), p-2, p) % p
        inE = in_E(x)
        stable = all(fval(mul_OK(x, tb, ell)) == 0 for tb in (Tbasis if ell < 30 else random.sample(Tbasis, 8)))
        if inE: okE += stable
        else: okNotE += (not stable)
        if inE != stable and ell < 30: raise SystemExit('End(T) mismatch')
    out['sampled x: (in E and xT<=T) count, (notin E and xT notin T) count'] = [okE, okNotE]
    # elements of E constructed deliberately: R generators and pO_K
    out['pi, pibar, 2, p*zeta in claimed End(T)'] = [in_E(G['pi']), in_E(G['pib']), in_E([2]+[0]*(n-1)), in_E([p*a for a in G['z']])]

    # ---------- explicit SUPER-lattice T' = O_K + Z*(v/p) : kernel of A -> A' inside A[p] = p^-1 O_K/O_K ----------
    if Fq.r is None:
        nu = random.randrange(1, p)
        e1 = idempotent(Fq, ell, ch1[0], ch1[1], ch2[1] if False else Fq.sub(Fq.e(-1), ch1[1]))
        e2 = idempotent(Fq, ell, ch2[0], ch2[1], Fq.sub(Fq.e(-1), ch2[1]))
        vF = [Fq.add(x, Fq.mul((nu, 0), y)) for x, y in zip(e1, e2)]
    else:
        th = (random.randrange(p), random.randrange(1, p))
        e1 = idempotent(Fq, ell, ch1[0], ch1[1], Fq.sub(Fq.e(-1), ch1[1]))
        vF = [Fq.add(Fq.mul(th, x), Fq.conj(Fq.mul(th, x))) for x in e1]   # th*e1 + Frob(th*e1)
    assert all(x[1] == 0 for x in vF), "kernel generator not F_p-rational"
    v = [x[0] for x in vF]
    def prop(u, lv):  # u == lv * v mod p ?
        return all((a - lv*b) % p == 0 for a, b in zip(u, v))
    lamp = lam[0]
    out["T': pi*v == lambda*v mod p (kernel point is a pi-eigenvector, lambda=chi1(pi))"] = prop(mul_OK(G['pi'], v, ell), lamp)
    out["T': pibar*v == (2/lambda)*v mod p"] = prop(mul_OK(G['pib'], v, ell), 2*pow(lamp, p-2, p) % p)
    zv = mul_OK(G['z'], v, ell)
    # zeta*v in F_p v + pO_K ?  find a coordinate where v nonzero, solve
    j = next(k for k in range(n) if v[k] % p); sc = zv[j]*pow(v[j], p-2, p) % p
    out["T': zeta*v in F_p*v + pO_K (would mean O_K-stable)"] = prop(zv, sc)
    out["T': v nonzero mod p"] = any(a % p for a in v)
    return out

if __name__ == '__main__':
    res = []
    cases = [(11, 23), (19, 457), (131, 263), (131, 146505763881528721)]
    for ell, p in cases:
        r = run(ell, p); res.append(r)
        print(json.dumps(r, default=str, indent=1)); sys.stdout.flush()
    json.dump(res, open('deligne_chars.json', 'w'), default=str, indent=1)
