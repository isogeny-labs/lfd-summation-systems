#!/usr/bin/env python3
"""Independent finite regressions for the proof's algebra, not a growth replay."""
import hashlib,json
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    x,y=s.symbols('x y');h=x*y-x-y-8
    S=(x-y)**2-2*((x+y)*(x*y+6)-14)+(x*y-6)**2+28*(x+y)
    assert s.expand(S-h*h)==0
    delta=lambda f:s.expand(y*s.diff(f,y))
    assert s.expand(delta(h)-y*(x-1))==0
    assert s.cancel(delta(h).subs(y,(x+8)/(x-1))-(x+8))==0
    assert s.expand(delta(h*h)-2*h*delta(h))==0
    cells=[]
    for k,p in [(2,17),(3,257),(4,65537),(5,641),(6,274177)]:
        N=1<<k;assert s.isprime(p) and (p-1)%(2*N)==0
        # Freshman's dream reduces each exponent separately in the chain quotient.
        poly=s.Poly(h,x,y,domain=s.ZZ)
        normal={(a*p%N,b*p%N):pow(int(c)%p,p,p) for (a,b),c in poly.terms()}
        expected={e:int(c)%p for e,c in poly.terms()}
        assert normal==expected
        # Explicit proper-ideal point on the Fermat-factor diagonal for these cases.
        assert pow(2,N,p)==p-1 and (4*4-4-4-8)%p==0 and pow(4,N,p)==1
        cells.append(dict(k=k,p=p,frobenius_factor_fixed=True,proper_grid_point=[4,4]))
    # Deliberate algebra changes must be detected; this is not an all-depth test.
    assert s.expand(S-(h+1)**2)!=0
    assert s.cancel(delta(h).subs(y,(x+8)/(x-1)))!=0
    out=dict(status='passed',scope='exact symbolic identities and finite Frobenius checks, not an asymptotic closure computation',
             cells=cells,negative_controls=['altered-factor-square','witness-derivative-does-not-vanish'],
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'algebra-regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
