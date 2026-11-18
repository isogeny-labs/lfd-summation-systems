#!/usr/bin/env python3
"""Bounded certificate and arithmetic regressions with explicit negative controls."""
import copy
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from csd import ic_bounds_lift_independent as lift
from csd.ic_fall_certificates import chain_degree_closure
from derivations import produce
from checker import verify
HERE=Path(__file__).resolve().parent


def main():
    cases=[]; controls=[]
    for k,D,X in [(2,3,4271),(2,4,4271),(3,4,4271),(4,3,4271),(4,4,4271),(4,4,2)]:
        p=65537; a,b=(1,0) if X==4271 else (4021,9917)
        gen={(i|(j<<k)):c for (i,j),c in lift.s3_generator(p,a,b,X,k).items()}
        start=time.monotonic();cert=produce([gen],2,k,p,D)
        rows=verify(cert,generators=[gen],parameters=(2,k,p,D))
        ref=chain_degree_closure([gen],2,k,p,D)
        assert len(rows)==ref.quotient_rank
        assert all(ref.contains_normal_form(row) for row in rows)
        # Producer rows are echelon with distinct pivots, checked independently.
        piv=[max(row,key=lambda e:(e.bit_count(),e)) for row in rows]
        assert len(set(piv))==len(rows)
        for row,pivot in zip(rows,piv): assert row[pivot]==1
        case=dict(k=k,D=D,X=X,rank=len(rows),nodes=len(cert['nodes']),
                  reference_contains_one=ref.contains_one,seconds=round(time.monotonic()-start,3),
                  row_sha256=cert['row_sha256'])
        cases.append(case)
        if k==2 and D==4:
            (HERE/'certificate-k2-D4.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
            for name,mutate in [
                ('forged-output',lambda c:c['rows'][0].append([0,1])),
                ('forward-reference',lambda c:c['nodes'].append({'op':'lin','terms':[[len(c['nodes']),1]]})),
                ('wrong-input',lambda c:c['generators'][0].append([255,1])),
                ('wrong-field',lambda c:c.update(p=17)),
                ('unknown-operation',lambda c:c['nodes'].append({'op':'trust-me'}))]:
                bad=copy.deepcopy(cert);mutate(bad)
                try: verify(bad,generators=[gen],parameters=(2,k,p,D))
                except (ValueError,KeyError,IndexError): controls.append(name)
                else: raise AssertionError(f'accepted {name}')
    # A carry drops the output degree, but the pre-reduction degree is illegal.
    bad=dict(m=1,k=3,p=17,D=2,generators=[[[3,1]]],nodes=[{'op':'input','index':0},
        {'op':'shift','source':0,'variable':0}],outputs=[],rows=[],row_sha256='')
    try: verify(bad,generators=[{3:1}],parameters=(1,3,17,2))
    except ValueError as e:
        assert 'unreduced' in str(e);controls.append('carry-does-not-authorize-overdegree-shift')
    else: raise AssertionError('accepted inadmissible carry')
    arithmetic=[]
    for p in [65537,163841,16777213]:
        M=np.zeros((64,65));M[:63,:63]=np.eye(63);M[:63,64]=p-1
        M[63,:63]=p-1;M[63,63]=1;M[63,64]=p-1
        R,piv=lift.rref(M,p,np)
        expected=(p-1-63*(p-1)**2)%p
        assert len(piv)==64 and int(R[-1,64])==expected
        arithmetic.append(dict(p=p,last_free=int(R[-1,64]),expected=expected))
    for unsafe in [0,-1]:
        try: lift.mm(np.ones((1,1)),np.ones((1,1)),16777213,np,unsafe)
        except lift.Refuse: controls.append(f'unsafe-inner-{unsafe}')
        else: raise AssertionError('unsafe inner accepted')
    files=['csd/ic_bounds_lift_independent.py']+[str(p.relative_to(ROOT)) for p in HERE.glob('*.py')]
    out=dict(status='passed',scope='bounded derivation prototype and high-prime arithmetic; not a depth-14 certificate',
             cases=cases,negative_controls=controls,arithmetic=arithmetic,
             source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)})
    (HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
