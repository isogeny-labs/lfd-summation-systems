#!/usr/bin/env python3
"""Head-to-head: monomial V vs sliding-window Frobenius-span V, for LRF index calculus.

Frobenius-span V = <w, w^2, w^4, ..., w^(2^(k-1))> for a generic w (full-length Frobenius orbit).
phi(V) overlaps V in dim k-1, so V-points share Frobenius orbits -> a real column reduction,
unlike the monomial basis. Measured per basis: dim V2 (LRF feasibility), yield/system,
saturated per-point columns vs orbit columns (true linear-algebra reduction) over E[ell].
"""
import argparse, json, random, sys
sys.path.insert(0, ".")
from csd import ic_dreg as I
from csd import ic_rootfind as rf

TAILS={11:(1<<2)|1,13:(1<<4)|(1<<3)|(1<<1)|1,17:(1<<3)|1,19:(1<<5)|(1<<2)|(1<<1)|1,23:(1<<5)|1}

def order_t(n):
    pts=1
    for x in (0,1):
        for y in (0,1):
            if (y*y^x*y)==(x*x*x^1): pts+=1
    t=3-pts; s=[2,t]
    for i in range(2,n+1): s.append(t*s[-1]-2*s[-2])
    return 2**n+1-s[n], t
def lpf(m):
    f,big=2,1
    while f*f<=m:
        while m%f==0: big=f; m//=f
        f+=1
    return big if m==1 else m
def span(bs):
    S={0}
    for v in bs: S|={s^v for s in S}
    return S
class ModRank:
    def __init__(s,ell): s.ell=ell; s.piv={}; s.rank=0
    def add(s,row):
        ell=s.ell; rr={c:w%ell for c,w in row.items() if w%ell}
        while rr:
            p=max(rr)
            if p in s.piv:
                f=rr[p]
                for c,w in s.piv[p].items(): rr[c]=(rr.get(c,0)-f*w)%ell
                rr={c:w for c,w in rr.items() if w}
            else:
                inv=pow(rr[p],ell-2,ell); s.piv[p]={c:(w*inv)%ell for c,w in rr.items()}; s.rank+=1; return True
        return False

def frob_basis(F,n,k,rng):
    sq=lambda x:F.mul(x,x)
    while True:
        w=rng.randrange(2,1<<n)
        orb=[]; e=w; ok=True
        for i in range(k):
            orb.append(e); e=sq(e)
        # require the k window elements independent and orbit length >= k
        if len(span(orb))==1<<k: return orb

def run(n,k,basis_kind,targets,seed):
    F=I.GF2n(n,(1<<n)|TAILS[n]); E=I.BinaryCurve(F,0,1)
    order,t=order_t(n); ell=lpf(order); sq=lambda x:F.mul(x,x)
    rng=random.Random(seed)
    if basis_kind=="mono": basis=[1<<i for i in range(k)]
    else: basis=frob_basis(F,n,k,rng)
    dimV2=len(rf.product_space(F,basis))
    if dimV2>=n:  # LRF infeasible
        return dict(n=n,k=k,basis=basis_kind,dimV2=dimV2,feasible=False)
    V=span(basis)
    Vpoints=set()
    for x in V:
        P=E.lift_x(x)
        if P is None: continue
        Vpoints.add(P); Vpoints.add(E.neg(P))
    # orbit map over the Frobenius orbits of V-points
    phi=lambda P:(sq(P[0]),sq(P[1]))
    seen=set(); orbits=[]
    for P in list(Vpoints):
        if P in seen: continue
        o=[]; Q=P
        while Q not in seen and Q in Vpoints:
            seen.add(Q); o.append(Q); Q=phi(Q)
        # orbit may leave Vpoints; include the full cycle rep anyway
        orbits.append(o if o else [P])
    rep_of={}
    for oid,o in enumerate(orbits):
        for s,P in enumerate(o): rep_of[P]=(oid,s)
    tau=next(c for c in range(1,ell) if (c*c-t*c+2)%ell==0)
    taupow=[1]
    for _ in range(n): taupow.append((taupow[-1]*tau)%ell)
    S=rf.Solver(E,basis)
    ptcol={}; rk_pt=ModRank(ell); rk_orb=ModRank(ell)
    agg=dict(targets=0,found=0); decomps=0; searches=0; stall=0
    while searches<targets and stall<4000:
        R=E.random_point(rng); rows,st,_=S.lrf(3,R)
        searches+=1; agg["targets"]+=st["targets"]; agg["found"]+=st["found"]; grew=False
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
                if P not in rep_of: rep_of[P]=(len(orbits),0); orbits.append([P])
                oid,sft=rep_of[P]; orow[oid]=orow.get(oid,0)+taupow[sft%len(taupow)]
            if rk_orb.add(orow): grew=True
        stall=0 if grew else stall+1
    return dict(n=n,k=k,basis=basis_kind,dimV2=dimV2,feasible=True,delta=k+dimV2-n,
                Vpoints=len(Vpoints),searches=searches,systems=agg["targets"],decomps=decomps,
                yield_per_system=agg["found"]/agg["targets"] if agg["targets"] else 0,
                cols_perpoint=rk_pt.rank,cols_orbit=rk_orb.rank,
                real_LA_reduction=rk_pt.rank/rk_orb.rank if rk_orb.rank else None)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--targets",type=int,default=8000); ap.add_argument("--seed",type=int,default=1)
    a=ap.parse_args()
    for n in (11,13,17,19,23):
        for k in (n//3+1, n//3+2):
            for kind in ("mono","frob"):
                try:
                    r=run(n,k,kind,a.targets,a.seed); print(json.dumps(r))
                except Exception as e:
                    print(json.dumps(dict(n=n,k=k,basis=kind,error=str(e))))
