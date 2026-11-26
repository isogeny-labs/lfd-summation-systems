#!/usr/bin/env python3
"""Measured ANF structure of the projected root-membership indicator.

At high delta the 2-summand linear system has a delta-dim affine solution space (2^delta
candidates); only a few give a quadratic with BOTH roots in V. The
indicator's ANF degree and support density are diagnostics, not a lower
bound on the cost of solving equations with auxiliary root coordinates.
Even a unique solution to linear equations has a full-degree indicator.

For many solvable targets we build the Boolean indicator g over the delta kernel coordinates and
measure: ANF degree of g, whether g^{-1}(1) is affine, and the density. Baseline control: a random
Boolean function of the same delta and density has ANF degree ~delta (full) with overwhelming prob.
"""
import argparse, json, random, sys
sys.path.insert(0, ".")
from csd import ic_rootfind as rf
from csd import ic_dreg as I

TAILS={17:(1<<3)|1,19:(1<<5)|(1<<2)|(1<<1)|1,23:(1<<5)|1,29:(1<<2)|1,31:(1<<3)|1}

def anf_degree(bits, d):
    """bits: list of 2^d values g(x) for x=0..2^d-1 (x bit i = coord i). Return (anf_degree, n_ones)."""
    f=list(bits)
    N=1<<d
    # Moebius transform (ANF) in place
    i=1
    while i<N:
        for j in range(N):
            if j & i:
                f[j]^=f[j^i]
        i<<=1
    deg=0
    for m in range(N):
        if f[m]:
            w=bin(m).count("1")
            if w>deg: deg=w
    return deg, sum(bits)

def is_affine_set(points, d):
    """points: set of ints in [0,2^d). affine iff closed: for any a,b,c in set, a^b^c in set."""
    pts=list(points)
    if len(pts)<=1: return True
    base=pts[0]
    diffs={p^base for p in pts}
    # diffs must be a linear subspace: closed under xor
    basis=[]
    for v in diffs:
        x=v
        for b in basis: x=min(x,x^b)
        if x: basis.append(x); basis.sort(reverse=True)
    return (1<<len(basis))==len(diffs)

def run(n,k,targets,seed,want_delta):
    F=I.GF2n(n,(1<<n)|TAILS[n]); rng=random.Random(seed)
    b=rng.randrange(1,1<<n); E=I.BinaryCurve(F,0,b)
    basis=[1<<i for i in range(k)]
    S=rf.Solver(E,basis)
    Qf=lambda x:rf._apply(S.Q,x); pif=lambda x:rf._apply(S.pi,x); top=S.w-1
    degs=[]; densities=[]; affine_cnt=0; systems=0; deltas=[]
    tried=0
    while systems<targets and tried<targets*200:
        tried+=1
        R=E.random_point(rng)
        # pick a random factor-base point P3 to form T=R-P3
        x3=rng.choice(list(range(1,1<<k)))
        e=0
        for i in range(k):
            if x3>>i&1: e^=basis[i]
        P3=E.lift_x(e)
        if P3 is None: continue
        T=E.add(R,E.neg(P3))
        xT,yT=T
        if xT==0: continue
        ym=yT^xT; c1=E.a^xT; c2=F.mul(F.sq(ym)^E.b,F.inv(xT))
        cols=[pif(F.mul(xT,Qf(v)^v))|(F.trace(v)<<top) for v in basis]+[pif(xT)]
        rhs=pif(F.mul(xT,Qf(c1)^c1)^c2)|(F.trace(c1)<<top)
        sol=I._f2_solve(cols,rhs,S.w)
        if sol is None: continue
        part,ker=sol; d=len(ker)
        if d!=want_delta: continue
        if d>16: continue
        systems+=1; deltas.append(d)
        # enumerate kernel coords 0..2^d-1, g=1 iff both roots in V
        good=[0]*(1<<d)
        goodset=set()
        for mask in range(1<<d):
            s=part
            mm=mask;
            while mm:
                s^=ker[(mm & -mm).bit_length()-1]; mm&=mm-1
            e1=0
            for i in range(k):
                if s>>i&1: e1^=basis[i]
            eps=s>>k&1
            lam=Qf(e1^c1)^eps
            e2=F.mul(xT,F.sq(lam))^c2
            if e1==0:
                r0=F.sqrt(e2); ok=(rf._apply(S.chk,r0)==0)
            else:
                dd=F.mul(e2,F.inv(F.sq(e1)))
                if F.trace(dd): ok=False
                else:
                    x1=F.mul(e1,Qf(dd)); ok=(rf._apply(S.chk,x1)==0 and rf._apply(S.chk,x1^e1)==0)
            if ok: good[mask]=1; goodset.add(mask)
        deg,ones=anf_degree(good,d)
        if ones==0:
            systems-=1; deltas.pop(); continue   # skip empty (no info on structure)
        degs.append(deg); densities.append(ones/(1<<d))
        if is_affine_set(goodset,d): affine_cnt+=1
    import statistics as st
    return dict(n=n,k=k,delta=want_delta,systems_with_good=len(degs),
                anf_degree_mean=st.mean(degs) if degs else None,
                anf_degree_max=max(degs) if degs else None,
                anf_degree_min=min(degs) if degs else None,
                affine_fraction=affine_cnt/len(degs) if degs else None,
                density_mean=st.mean(densities) if densities else None)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--seed",type=int,default=1)
    ap.add_argument("--targets",type=int,default=200); a=ap.parse_args()
    # choose (n,k,delta): delta=3k-1-n, keep dimV2=2k-1<n
    plan=[(23,10,6),(23,11,9),(29,12,6),(29,13,9),(29,14,12)]
    for (n,k,dlt) in plan:
        print(json.dumps(run(n,k,a.targets,a.seed,dlt)))
