"""Independent finite algebra checks for the E6 first-derivative proof."""
import hashlib,json
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'experiments/lfd-jet-polynomial-20261118'


def test_e6_factor_and_nonzero_derivative():
    x,y=sp.symbols('x y');h=x*y-x-y-8
    semaev=(x-y)**2-2*((x+y)*(x*y+6)-14)+(x*y-6)**2+28*(x+y)
    assert sp.expand(semaev-h**2)==0
    derivative=y*sp.diff(h,y)
    assert sp.cancel(derivative.subs(y,(x+8)/(x-1)))==x+8
    assert sp.expand(y*sp.diff(h**2,y)-2*h*derivative)==0
    assert sp.expand(semaev-(h+1)**2)!=0


def test_frobenius_proper_ideal_cells_are_preserved():
    result=json.loads((HERE/'algebra-regressions.json').read_text())
    for case in result['cells']:
        k,p=case['k'],case['p'];N=2**k
        assert sp.isprime(p) and (p-1)%(2*N)==0
        terms={(1,1):1,(1,0):-1,(0,1):-1,(0,0):-8}
        assert {(a*p%N,b*p%N):pow(c%p,p,p) for (a,b),c in terms.items()}=={e:c%p for e,c in terms.items()}
        assert pow(4,N,p)==1 and (4*4-4-4-8)%p==0


def test_new_algebra_receipt_source_matches():
    result=json.loads((HERE/'algebra-regressions.json').read_text())
    assert hashlib.sha256((HERE/'check_algebra.py').read_bytes()).hexdigest()==result['source_sha256']
