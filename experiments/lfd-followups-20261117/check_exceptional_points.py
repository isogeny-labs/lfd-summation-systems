#!/usr/bin/env python3
"""Independent O(2^k) enumeration using the order-two target's involution."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from csd import ic_curve_sweep as S
from csd import ic_bounds_d5_verify as R
HERE=Path(__file__).resolve().parent

def main():
    records=[]
    for k,p in [(7,257),(7,769),(8,7681),(9,12289),(10,12289)]:
        E=S.I.PrimeCurve(p,6,-7);h=S.subgroup_generator_value(p,k)
        H={pow(h,j,p) for j in range(1<<k)};assert len(H)==1<<k
        linear={(1,1):1,(1,0):-1,(0,1):-1,(0,0):-8};square={}
        for (i,j),c in linear.items():
            for (u,v),d in linear.items():square[i+u,j+v]=(square.get((i+u,j+v),0)+c*d)%p
        square={e:c for e,c in square.items() if c}
        assert square==S.I.s3_prime_coeffs(E,1)
        actual=set()
        for x in H:
            if x==1:
                assert (-9)%p;continue
            y=(x+8)*pow(x-1,-1,p)%p
            if y in H:actual.add((x,y))
        hr,hrvalues=R.subgroup(p,k)
        reference={(hrvalues[i],hrvalues[j]) for i,j in R.solutions(p,6,p-7,1,k)}
        assert actual==reference
        records.append(dict(k=k,p=p,points=len(actual),points_sha256=hashlib.sha256(json.dumps(sorted(actual),separators=(',',':')).encode()).hexdigest(),
                            square_identity_checked=True,involution_equals_full_grid=True))
    result=dict(status='passed',records=records,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'exceptional-points-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
