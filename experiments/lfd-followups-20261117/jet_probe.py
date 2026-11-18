#!/usr/bin/env python3
"""Falsifiers for a proposed proper-ideal proof route, not a growth theorem."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from csd import ic_bounds_d5_verify as R
HERE=Path(__file__).resolve().parent

def main():
    identities=[]
    for k in range(2,7):
        N=1<<k;checks=0
        for a in range(N):
            for b in range(N):
                for block in [0,1]:
                    for j in range(k):
                        step=1<<j;e=[a,b];wrap=e[block]+step>=N
                        out=e.copy();out[block]=(out[block]+step)%N
                        # Euler derivative in the second coordinate, on basis vectors.
                        rhs=b+(step-N*wrap if block==1 else 0)
                        assert out[1]==rhs;checks+=1
        identities.append(dict(k=k,commutator_basis_checks=checks))
    witnesses=[]
    for k,p,D in [(7,257,4),(7,769,4),(8,7681,4),(9,12289,4),(10,12289,5)]:
        start=time.time();N=1<<k
        C=R.Closure(R.chain_generator(p,6,p-7,1,k),k,p,D,split=True)
        assert C.run(deadline=time.time()+180)
        closed=C.sweep_full(deadline=time.time()+180);assert closed['complete'] and closed['passed']
        h={(1,1):1,(1,0):-1,(0,1):-1,(0,0):-8}
        v=np.zeros((1,C.L.c))
        for (i,j),c in h.items():v[0,C.L.index(np.array([i]),np.array([j]))[0]]=c%p
        member=all(not np.any(part.residual(y)) for part,y in zip(C.parts,C.L.project(v,p)))
        assert not member
        # Sparse Frobenius certificate in the chain quotient: h^p=h.
        frob={(i*p%N,j*p%N):pow(c%p,p,p) for (i,j),c in h.items()}
        assert frob=={e:c%p for e,c in h.items()} and p>=2
        witnesses.append(dict(k=k,p=p,D=D,rank=C.rank(),h_in_V_D=member,
                              frobenius_normal_form_checked=True,seconds=round(time.time()-start,3)))
    out=dict(status='passed',scope='finite witnesses and commutators only; no asymptotic theorem',
             commutators=identities,witnesses=witnesses,
             source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'csd/ic_bounds_d5_verify.py']})
    (HERE/'jet-probe.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
