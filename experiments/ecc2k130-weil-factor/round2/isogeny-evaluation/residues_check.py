#!/usr/bin/env python3
"""Method B (pure Python, independent of PARI): self-glued residue of R = Z[pi,pibar] at p | U_ell, and the
multiplicative order of lambda = pi mod v0 (which forces kernel points into F_{2^m}, ord(lambda) | m).
F_{p^2} built by hand as F_p[t]/(t^2+t+2) (inert case) ; F_p in the split case."""
import random, sys

def is_prime(n):
    if n < 2: return False
    for q in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n % q == 0: return n == q
    d, s = n-1, 0
    while d % 2 == 0: d//=2; s+=1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37,41]:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def factor(n):
    f, d = {}, 2
    while d*d <= n and d < 10**6:
        while n % d == 0: f[d] = f.get(d,0)+1; n//=d
        d += 1
    if n > 1:
        assert is_prime(n), n   # all cofactors here are prime (checked)
        f[n] = f.get(n,0)+1
    return f

def lucasU(n):
    a, b = 0, 1
    for _ in range(n-1): a, b = b, -b-2*a
    return b

class F2:  # F_p[t]/(t^2+t+2)
    def __init__(s, p): s.p = p
    def mul(s, x, y):
        p = s.p; a0,a1 = x; b0,b1 = y
        c0 = a0*b0; c1 = a0*b1+a1*b0; c2 = a1*b1   # t^2 = -t-2
        return ((c0-2*c2) % p, (c1-c2) % p)
    def pw(s, x, e):
        r = (1,0)
        while e:
            if e&1: r = s.mul(r,x)
            x = s.mul(x,x); e >>= 1
        return r
    def inv(s, x): return s.pw(x, s.p**2-2)

def order(x, n, f, mul, pw, one):
    o = n
    for q in f:
        while o % q == 0 and pw(x, o//q) == one: o //= q
    return o

def run(ell):
    U = abs(lucasU(ell)); fU = factor(U)
    print(f"ell={ell} U={U} factors={fU}")
    for p in fU:
        if (p+1) % ell == 0:   # inert type: work in F_{p^2}
            F = F2(p); a = (0,1); b = F.mul((2,0), F.inv(a))
            assert F.pw(a, p) == b, "alpha^p should be beta (p inert in Q(sqrt-7))"
            n2 = p*p-1
            while True:
                g = (random.randrange(p), random.randrange(p))
                z = F.pw(g, n2//ell)
                if z != (1,0): break
            hits = []
            for j in range(ell):
                pij = F.mul(a, F.pw(z, j)); pbj = F.mul(b, F.pw(z, (ell-j) % ell))
                if pij == pbj: hits.append((j, pij))
            assert len(hits) == 1, hits
            lam = hits[0][1]
            assert lam[1] == 0 and F.mul(lam, lam) == (2 % p, 0)
            l = lam[0]
        elif (p-1) % ell == 0:  # split type: F_p
            r = [x for x in range(p) if (x*x+x+2) % p == 0] if p < 10**6 else None
            a = r[0]; b = 2*pow(a, -1, p) % p
            while True:
                z = pow(random.randrange(2, p), (p-1)//ell, p)
                if z != 1: break
            hits = [j for j in range(ell) if a*pow(z,j,p) % p == b*pow(z,-j,p) % p]
            assert len(hits) == 1
            l = a*pow(z, hits[0], p) % p
        else:
            print("  p", p, "other type"); continue
        fp = factor(p-1)
        mul = lambda x,y: x*y % p
        o = order(l, p-1, fp, mul, lambda x,e: pow(x,e,p), 1)
        om = order((-l) % p, p-1, fp, mul, lambda x,e: pow(x,e,p), 1)
        print(f"  p={p} prime={is_prime(p)} lambda={l} lambda^2%p={l*l%p} ord(lambda)={o} (p-1)/ord={(p-1)//o} ord(-lambda)={om} p-1={fp}")

random.seed(1)
for ell in [11, 17, 19, 131]: run(ell)
