#!/usr/bin/env python3
"""A proof-carrying counterexample to dropping the wrapping hypothesis."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'experiments/lfd-followups-20261117')]
from csd.ic_bounds_lift_independent import s3_generator
from derivations import produce
from checker import verify
HERE=Path(__file__).resolve().parent


def main():
    p,k,D=17,2,3
    g={(a|(b<<k)):c for (a,b),c in s3_generator(p,6,-7,1,k).items()}
    cert=produce([g],2,k,p,D)
    target=cert['row_sha256']
    rows=verify(cert,generators=[g],parameters=(2,k,p,D),expected_output_sha256=target)
    x=2;y=(x+8)*pow(x-1,-1,p)%p
    for i,row in enumerate(rows):
        phi=sum(c*pow(x,e%4,p)*pow(y,e//4,p) for e,c in row.items())%p
        jet=sum(c*(e//4)*pow(x,e%4,p)*pow(y,e//4,p) for e,c in row.items())%p
        if phi or jet:break
    else:raise AssertionError('expected falsifier not found')
    assert phi!=0 or jet!=0
    (HERE/'counterexample-certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    files=[Path(__file__),ROOT/'experiments/lfd-followups-20261117/derivations.py',ROOT/'experiments/lfd-followups-20261117/checker.py']
    result=dict(k=k,p=p,D=D,curve=[6,-7],target=[1,0],rational_sample=[x,y],
        checked_output_sha256=target,row_index=i,row_terms=sorted(row.items()),phi=phi,phi_euler=jet,
        status='derived counterexample to unqualified jet invariance',
        scope='ordinary normalized chain closure modulo degree-compatible chain relations; no growth theorem',
        source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    (HERE/'certified-falsifier.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
