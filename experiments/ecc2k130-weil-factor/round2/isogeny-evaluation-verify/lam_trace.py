#!/usr/bin/env python3
"""Independent check (no F_{p^2}, no zeta): at a self-glued hom pi==pibar=lam, lam^2=2 and
lam^ell = alpha^ell = beta^ell, so t_ell = alpha^ell+beta^ell == 2*lam^ell (mod p). This picks the sign
of sqrt(2) mod p. Also recompute U_ell, its factorisation, primality, and the multiplicative orders."""
import random, math
def trace(n):  # t_k = alpha^k+beta^k, alpha+beta=-1, alpha*beta=2
    a,b=2,-1
    for _ in range(n-1): a,b=b,-b-2*a
    return b
def lucasU(n):  # (alpha^n-beta^n)/(alpha-beta)
    a,b=0,1
    for _ in range(n-1): a,b=b,-b-2*a
    return b
def mr(n,k=64):
    if n<2: return False
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n%p==0: return n==p
    d,s=n-1,0
    while d%2==0: d//=2;s+=1
    rng=random.Random(12345)
    for _ in range(k):
        a=rng.randrange(2,n-1);x=pow(a,d,n)
        if x in(1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True
def factor(n):
    f={};p=2
    while p*p<=n and p<10**6:
        while n%p==0: f[p]=f.get(p,0)+1;n//=p
        p+=1
    if n>1:
        assert mr(n), n
        f[n]=f.get(n,0)+1
    return f
def sqrtmod(a,p):  # Tonelli-Shanks
    assert pow(a,(p-1)//2,p)==1
    q,s=p-1,0
    while q%2==0:q//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1: z+=1
    m,c,t,r=s,pow(z,q,p),pow(a,q,p),pow(a,(q+1)//2,p)
    while t!=1:
        i,tt=0,t
        while tt!=1: tt=tt*tt%p;i+=1
        b=pow(c,1<<(m-i-1),p);m,c,t,r=i,b*b%p,t*b*b%p,r*b%p
    return r
def order(x,p):
    o=p-1
    for q in factor(p-1):
        while o%q==0 and pow(x,o//q,p)==1: o//=q
    return o
expected={(11,23):18,(17,271):96,(19,457):377,(131,263):110,(131,146505763881528721):14975904886999252}
for ell in (11,17,19,131):
    U=abs(lucasU(ell)); t=trace(ell); F=factor(U)
    print("ell=%d U=%d factors=%s t_ell=%d"%(ell,U,F,t))
    for p in F:
        s=sqrtmod(2,p); half=t*pow(2,-1,p)%p
        cands=[x for x in (s,p-s) if pow(x,ell,p)==half]
        lam=cands[0] if len(cands)==1 else None
        print("  p=%d  #sign-consistent sqrt2=%d lam=%s expected=%s match=%s ord(lam)=%d (p-1)/ord=%d ord(-lam)=%d ord(2)=%d"%(
            p,len(cands),lam,expected.get((ell,p)),lam==expected.get((ell,p)),order(lam,p),(p-1)//order(lam,p),order(p-lam,p),order(2,p)))
    # negative control: a prime not dividing U gives no consistent lam with lam^2=2 & t=2lam^ell? (just report)
p2=146505763881528721
print("p2-1 factorisation:",factor(p2-1)," (p2-1)/12 =",(p2-1)//12)
print("p2 mod 131 =",p2%131," kronecker(-7,p2) via Euler:",pow(-7%p2,(p2-1)//2,p2)==1)
# negative control: random primes q !| U_131 with sqrt2 existing: check t_131 == 2 lam^131 fails for both signs
rng=random.Random(7);bad=0;tested=0
while tested<200:
    q=rng.randrange(10**6,10**9)
    if not mr(q) or q in (263,p2) or pow(2,(q-1)//2,q)!=1: continue
    tested+=1; s=sqrtmod(2,q); h=trace(131)*pow(2,-1,q)%q
    if any(pow(x,131,q)==h for x in (s,q-s)): bad+=1
print("negative control: random primes with sqrt2 where t_131==2*lam^131 holds: %d of %d"%(bad,tested))
