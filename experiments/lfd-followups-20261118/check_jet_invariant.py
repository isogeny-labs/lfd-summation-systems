#!/usr/bin/env python3
"""Falsify unqualified jet invariance; finite rational evaluations only.

A nonzero evaluation certifies nonvanishing. Zero at the sampled points
is not claimed to prove a rational-function identity.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from csd import ic_bounds_d5_verify as R
HERE=Path(__file__).resolve().parent


def compute(k,p,D):
    g=R.chain_generator(p,6,p-7,1,k)
    C=R.Closure(g,k,p,D);assert C.run()
    closed=C.sweep_full();assert closed['complete'] and closed['passed']
    samples=[x for x in range(2,min(p,18))]
    vals=[];derivs=[]
    for x in samples:
        y=(x+8)*pow(x-1,-1,p)%p
        v=np.array([pow(x,int(a),p)*pow(y,int(b),p)%p for a,b in zip(C.L.a1,C.L.a2)],dtype=np.int64)
        vals.append(v);derivs.append(v*C.L.a2%p)
    E=np.asarray(vals,dtype=np.int64).T;J=np.asarray(derivs,dtype=np.int64).T
    bad=None;checked=0
    for q,P in enumerate(C.parts):
        for a in range(0,P.rank,32):
            rows=np.mod(C.L.lift(q,P.rows(np.arange(a,min(P.rank,a+32)))),p).astype(np.int64)
            plain=rows@E%p;jet=rows@J%p
            nz=np.argwhere((plain!=0)|(jet!=0))
            checked+=len(rows)
            if len(nz):
                i,j=map(int,nz[0]);r=rows[i];nzrow=np.flatnonzero(r)
                bad=dict(part=q,row=a+i,x=samples[j],phi=int(plain[i,j]),phi_euler=int(jet[i,j]),
                         row_terms=[[int(C.L.a1[c]),int(C.L.a2[c]),int(r[c])] for c in nzrow])
                break
        if bad:break
    return dict(k=k,p=p,D=D,rank=C.rank(),rows_sampled=checked,sample_count=len(samples),
                status='counterexample-to-unqualified-invariance' if bad else 'no-counterexample-in-samples',witness=bad)


def main():
    cases=[compute(k,17 if k==2 else 65537,D) for k in range(2,7) for D in [3,4]]
    assert any(r['witness'] for r in cases),'falsifier failed to exercise wrap regime'
    result=dict(scope='finite falsifiers; no all-depth invariant or sampled-zero identity claim',cases=cases,
      source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'csd/ic_bounds_d5_verify.py']})
    (HERE/'jet-invariant-falsifiers.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['k','p','D','rank','status']} for r in cases]))
if __name__=='__main__':main()
