#!/usr/bin/env python3
"""Measured Frobenius orbit reduction for LRF index calculus on a Koblitz curve.

Curve y^2+xy=x^3+1 (a=0,b=1) over F_(2^n), n prime; phi:(x,y)->(x^2,y^2) = [tau], tau^2-t tau+2=0.

Honest apples-to-apples: LRF searches ONE subspace V^(0)=<1,..,w^(k-1)> and finds decompositions
with all abscissas in V. (Searching the conjugates V^(j) only yields Frobenius images of these =
the SAME orbit rows, no new rank, so it is wasted.) From the SAME collected relations we build:
  (A) per-point columns  = the V-points that appear            -> saturated rank R_A
  (B) orbit columns      = Frobenius orbits those points touch -> saturated rank R_B (tau-weighted)
Real linear-algebra reduction = R_A / R_B. Also report the 'ideal' n_eff = |F|/|orbits| for the
fully Frobenius-closed factor base (what the textbook n-fold assumes but LRF-from-V^(0) does not reach).
"""
import argparse, json, random, sys
sys.path.insert(0, ".")
from csd import ic_dreg as I
from csd import ic_rootfind as rf

TAILS = {11:(1<<2)|1, 13:(1<<4)|(1<<3)|(1<<1)|1, 17:(1<<3)|1,
         19:(1<<5)|(1<<2)|(1<<1)|1, 23:(1<<5)|1}

def order_t(a,b,n):
    pts=1
    for x in (0,1):
        for y in (0,1):
            if (y*y^x*y)==(x*x*x^(a*x*x)^b): pts+=1
    t=3-pts; s=[2,t]
    for i in range(2,n+1): s.append(t*s[-1]-2*s[-2])
    return 2**n+1-s[n], t

def lpf(m):
    f,big=2,1
    while f*f<=m:
        while m%f==0: big=f; m//=f
        f+=1
    return big if m==1 else m

def span(basis):
    S={0}
    for v in basis: S|={s^v for s in S}
    return S

class ModRank:
    def __init__(self,ell): self.ell=ell; self.piv={}; self.rank=0
    def add(self,row):
        ell=self.ell; rr={c:w%ell for c,w in row.items() if w%ell}
        while rr:
            p=max(rr)
            if p in self.piv:
                f=rr[p]; pr=self.piv[p]
                for c,w in pr.items(): rr[c]=(rr.get(c,0)-f*w)%ell
                rr={c:w for c,w in rr.items() if w}
            else:
                inv=pow(rr[p],ell-2,ell)
                self.piv[p]={c:(w*inv)%ell for c,w in rr.items()}; self.rank+=1; return True
        return False

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,required=True); ap.add_argument("--k",type=int,required=True)
    ap.add_argument("--max-searches",type=int,default=40000); ap.add_argument("--seed",type=int,default=1)
    a=ap.parse_args(); n,k=a.n,a.k
    F=I.GF2n(n,(1<<n)|TAILS[n]); E=I.BinaryCurve(F,0,1)
    order,t=order_t(0,1,n); ell=lpf(order); sq=lambda x:F.mul(x,x)
    basis=[1<<i for i in range(k)]; V=span(basis)
    # full Frobenius closure (ideal factor base F) for n_eff
    absF=set()
    for x in V:
        y=x
        for _ in range(n): absF.add(y); y=sq(y)
    fbF=[]
    for x in sorted(absF):
        P=E.lift_x(x)
        if P is None: continue
        fbF.append(P); Q=E.neg(P)
        if Q!=P: fbF.append(Q)
    phi=lambda P:(sq(P[0]),sq(P[1]))
    seen=set(); orbits=[]
    for P in fbF:
        if P in seen: continue
        o=[]; Q=P
        while Q not in seen: seen.add(Q); o.append(Q); Q=phi(Q)
        orbits.append(o)
    norbF=len(orbits); rep_of={}
    for oid,o in enumerate(orbits):
        for s,P in enumerate(o): rep_of[P]=(oid,s)
    tau=next(c for c in range(1,ell) if (c*c-t*c+2)%ell==0)
    taupow=[1]
    for _ in range(n): taupow.append((taupow[-1]*tau)%ell)

    # collect V^(0) relations
    S=rf.Solver(E,basis); rng=random.Random(a.seed)
    Vpoints=set()
    for x in V:
        P=E.lift_x(x)
        if P is None: continue
        Vpoints.add(P); Vpoints.add(E.neg(P))
    ptcol={}; rk_pt=ModRank(ell); rk_orb=ModRank(ell)
    searches=0; systems=0; decomps=0
    sat_pt=0; sat_orb=0; stall=0
    hist=[]
    while searches < a.max_searches and stall < 3000:
        R=E.random_point(rng); rows,stats,_=S.lrf(3,R)
        searches+=1; systems+=stats["targets"]; grew=False
        for row in rows:
            pts=[(int(x),int(y)) for (x,y) in row]
            if any(P not in Vpoints for P in pts): continue
            decomps+=1
            vec={}
            for P in pts:
                if P not in ptcol: ptcol[P]=len(ptcol)
                vec[ptcol[P]]=vec.get(ptcol[P],0)+1
            if rk_pt.add(vec): grew=True
            orow={}
            for P in pts:
                oid,sft=rep_of[P]; orow[oid]=orow.get(oid,0)+taupow[sft]
            if rk_orb.add(orow): grew=True
        stall = 0 if grew else stall+1
    out=dict(n=n,k=k,ell=ell,tau=tau,t=t,order=order,
             fbF_points=len(fbF), n_orbits_F=norbF, n_eff_ideal=len(fbF)/norbF,
             Vpoints=len(Vpoints),
             searches=searches, systems=systems, decomps=decomps,
             cols_perpoint=rk_pt.rank, cols_orbit=rk_orb.rank,
             real_LA_reduction=rk_pt.rank/rk_orb.rank if rk_orb.rank else None)
    print(json.dumps(out))

if __name__=="__main__": main()
