import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ct_lib import build, s_of
from collections import Counter
ell,n,N=19,18,int(sys.argv[1])
q,g,dlog,mul,pw=build(n)
assert (dlog>=0).sum()==q-1 and len(set(dlog[dlog>=0].tolist()))==q-1, "modulus not irreducible"
print("ord_19(2) =", next(e for e in range(1,19) if pow(2,e,19)==1), "; F_2^18 dlog table complete")
s=s_of(ell); xs=np.arange(q,dtype=np.int64)
L0,L1,Lg=dlog[xs],dlog[xs^1],dlog[xs^g]
rng=np.random.default_rng(7); coh=Counter(); smax=0; done=0
while done<N:
    lam=int(rng.integers(2,q))
    if lam==g: continue
    a2,a3=int(rng.integers(1,ell)),int(rng.integers(1,ell))
    a4=(-(1+a2+a3))%ell
    if a4==0: continue
    Ll=dlog[xs^lam]; good=(L0>=0)&(L1>=0)&(Lg>=0)&(Ll>=0)
    d=(L0+a2*L1+a3*Lg+a4*Ll)%ell
    h=np.bincount(d[good],minlength=ell).astype(np.int64); h[0]+=1
    srt=np.sort(h); smax=max(smax,int(srt[-1]-srt[0])); done+=1
    if srt[1]==srt[-1]: delta=int(srt[0]-srt[1])
    elif srt[0]==srt[-2]: delta=int(srt[-1]-srt[0])
    else: continue
    assert abs(delta)<=2*int(q**0.5)+1
    coh[-delta]+=1
print(f"ell=19 q={q}: A needs trace {s}; sampled 4-rational covers={done}; max spread={smax} (Weil 2sqrt(q)={2*int(q**0.5)}); coherent={dict(coh)}; MATCH={coh.get(s,0)}")
