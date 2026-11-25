"""Matched finite-field screen of a two-torsion symmetrized S3 presentation.

The closure calculation is exact for unit ideals.  It charges all input
polynomials and iterates multiplication/reduction to a fixed degree cap.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.ic_dreg import PrimeCurve, prime_factor_base_size, prime_decompositions, prime_solutions  # noqa: E402


def add(a, b, p):
    out = dict(a)
    for m, c in b.items():
        out[m] = (out.get(m, 0) + c) % p
        if not out[m]:
            del out[m]
    return out


def mul(a, b, p):
    out = {}
    for (i, j), c in a.items():
        for (k, l), d in b.items():
            m = (i+k, j+l)
            out[m] = (out.get(m, 0) + c*d) % p
    return {m:c for m,c in out.items() if c}


def scale(a, c, p):
    return {m:(v*c)%p for m,v in a.items() if (v*c)%p}


def power(a, n, p):
    out = {(0, 0):1}
    for _ in range(n):
        out = mul(out, a, p)
    return out


def degree(a):
    return max((sum(m) for m in a), default=-1)


def value(a, x, y, p):
    return sum(c*pow(x,i,p)*pow(y,j,p) for (i,j),c in a.items())%p


def shift(a, idx, D):
    return {(i+(idx == 0), j+(idx == 1)):c for (i,j),c in a.items()
            if i+j+1 <= D}


def closure_unit_degree(gens, p, cap):
    """Exact iterated V_D, returning least D containing 1 and rank profile."""
    profile = []
    for D in range(cap+1):
        basis = {}
        queue = [g for g in gens if degree(g) <= D]
        at = 0
        while at < len(queue):
            row = dict(queue[at]); at += 1
            while row:
                lead = max(row, key=lambda m: (sum(m), m[0]))
                if lead not in basis:
                    inv = pow(row[lead], -1, p)
                    row = scale(row, inv, p)
                    basis[lead] = row
                    if degree(row) < D:
                        queue.append(shift(row, 0, D))
                        queue.append(shift(row, 1, D))
                    break
                row = add(row, scale(basis[lead], -row[lead], p), p)
        profile.append(len(basis))
        if (0, 0) in basis:
            return D, profile
    return None, profile


def pows(e, r, max_n, p):
    # u^j+v^j for roots of Z^2-eZ+r.
    vals = [{(0,0):2}, {(1,0):1}]
    for _ in range(2, max_n+1):
        vals.append(add(mul(e, vals[-1], p), scale(mul(r, vals[-2], p), -1, p), p))
    return vals


def symm_generators(p, N, t):
    # u=(x+1)/(x-1); x in mu_N minus {1}.  The numerator of x^N-1 is
    # F(u)=(u+1)^N-(u-1)^N, whose roots are exactly the chosen u-values.
    from math import comb
    c = {j:(comb(N,j)*(1-(-1)**(N-j)))%p for j in range(N+1)}
    c = {j:v for j,v in c.items() if v}
    e,r = {(1,0):1},{(0,1):1}
    ps = pows(e,r,N,p)
    fsum = {}
    fprod = {}
    for i,ci in c.items():
        fsum = add(fsum, scale(ps[i],ci,p),p)
        fprod = add(fprod,{(0,i):ci*ci%p},p)
        for j,cj in c.items():
            if j <= i: continue
            fprod = add(fprod,scale(mul(power(r,i,p),ps[j-i],p),ci*cj,p),p)
    relation = {(2,0):(1+t*t)%p,(0,1):(-4*t*t)%p,
                (0,2):(-2*t)%p,(0,0):(2*t)%p}
    return [relation, fsum, fprod]


def torsion_quotient_generators(p, N, t):
    # e -> -e is the simultaneous translation by T=(0,0).  Set q=e^2.
    # fsum=e*H(q,r), so fsum=0 is equivalent to q*H(q,r)=0.
    relation, fsum, fprod = symm_generators(p,N,t)
    assert all(i % 2 == 1 for (i,j) in fsum)
    assert all(i % 2 == 0 for (i,j) in fprod)
    qsum = {((i+1)//2,j):c for (i,j),c in fsum.items()}
    qprod = {(i//2,j):c for (i,j),c in fprod.items()}
    qrel = {(i//2,j):c for (i,j),c in relation.items()}
    return [qrel,qsum,qprod]


def direct_generators(p,N,t):
    x,y={(1,0):1},{(0,1):1}
    one={(0,0):1}
    xy=mul(x,y,p)
    a=add(x,scale(y,-1,p),p)
    b=add(x,y,p)
    s3=add(scale(power(a,2,p),t*t,p),
           scale(mul(b,add(xy,one,p),p),-2*t,p),p)
    s3=add(s3,power(add(xy,scale(one,-1,p),p),2,p),p)
    fx={(j,0):1 for j in range(N)}
    fy={(0,j):1 for j in range(N)}
    return [s3,fx,fy]


def main():
    cells=[]
    for p,k in ((97,3),(193,3),(257,3),(193,4),(257,4)):
        N=1<<k
        g=next(z for z in range(2,p) if pow(z,N,p)==1 and pow(z,N//2,p)!=1)
        H=sorted({pow(g,j,p) for j in range(N)})
        V=[x for x in H if x!=1]
        E=PrimeCurve(p,1,0)
        assert all(pow(x,-1,p) in V for x in V)
        T=(0,0)
        assert E.on_curve(T) and E.add(T,T) is None
        for x in V:
            P=E.lift_x(x)
            if P is not None:
                assert E.add(P,T)[0]==pow(x,-1,p)
        U=[(x+1)*pow(x-1,-1,p)%p for x in V]
        image={( (u+v)**2%p, u*v%p) for u in U for v in U}
        _, qsum, qprod = torsion_quotient_generators(p,N,0)
        roots={(q,r) for q in range(p) for r in range(p)
               if value(qsum,q,r,p)==value(qprod,q,r,p)==0}
        assert roots==image, (p,k,len(roots),len(image))
        targets=[]
        for x in range(p):
            R=E.lift_x(x)
            if R and R[1]:
                targets.append(R)
            if len(targets)==12: break
        rows=[]
        for R in targets:
            t=R[0]
            sols=prime_solutions(E,V,t)
            decomps=prime_decompositions(E,V,R,sols)
            # Every scalar solution maps to the invariant quadratic.
            image=set()
            for i,j in sols:
                u=(V[i]+1)*pow(V[i]-1,-1,p)%p
                v=(V[j]+1)*pow(V[j]-1,-1,p)%p
                e=(u+v)%p;r=u*v%p
                assert ((1+t*t)*e*e-4*t*t*r-2*t*(r*r-1))%p==0
                image.add((e*e%p,r))
            qrel=torsion_quotient_generators(p,N,t)[0]
            assert image == {(q,r) for q,r in roots if value(qrel,q,r,p)==0}
            direct=direct_generators(p,N,t)
            symm=symm_generators(p,N,t)
            torsion=torsion_quotient_generators(p,N,t)
            d_direct=degree(direct[0]);d_symm=degree(symm[0])
            row={'target':list(R),'scalar_ordered_pairs':len(sols),
                 'unordered_x_pairs_with_point_decomposition':decomps,
                 'invariant_solution_images':len(image),
                 'input_max_degree_direct':max(map(degree,direct)),
                 'input_max_degree_symmetrized':max(map(degree,symm)),
                 'input_max_degree_torsion_quotient':max(map(degree,torsion)),
                 'relation_degree_direct':d_direct,
                 'relation_degree_symmetrized':d_symm,
                 'relation_degree_torsion_quotient':degree(torsion[0])}
            if not sols:
                ld,pd=closure_unit_degree(direct,p,22)
                ls,ps=closure_unit_degree(symm,p,22)
                lt,pt=closure_unit_degree(torsion,p,22)
                row.update(unit_lfd_direct=ld,unit_lfd_symmetrized=ls,
                           unit_lfd_torsion_quotient=lt,
                           closure_rank_direct=pd,closure_rank_symmetrized=ps,
                           closure_rank_torsion_quotient=pt)
            rows.append(row)
        cells.append({'p':p,'k':k,'curve':'y^2=x^3+x',
                      'factor_base_rule':'x in mu_8 excluding x=1',
                      'x_values':V,'factor_base_points':prime_factor_base_size(E,V),
                      'torsion_quotient_factor_base_image_count':len(image),
                      'torsion_quotient_constraint_roots':len(roots),
                      'target_rule':'first 12 non-order-two rational points by ascending x; first square root',
                      'rows':rows})
    out={'status':'passed','scope':'matched small prime-field two-summand screen',
         'construction':'u=(x+1)/(x-1), e=u+v, r=uv, q=e^2; P=(1+t^2)q-4t^2r-2t(r^2-1)',
         'degree_note':'Unit-ideal LFD is the exact first degree containing 1 in iterated closure; relation input degree alone is not solver degree.',
         'source_sha256':{'screen.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                          'csd/ic_dreg.py':hashlib.sha256((ROOT/'csd/ic_dreg.py').read_bytes()).hexdigest()},
         'cells':cells}
    path=Path(__file__).with_name('results.json')
    path.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    for cell in cells:
        rows=cell['rows']
        unit=[r for r in rows if r['scalar_ordered_pairs']==0]
        print(cell['p'],len(unit),Counter((r['unit_lfd_direct'],r['unit_lfd_symmetrized'],r['unit_lfd_torsion_quotient']) for r in unit))

if __name__=='__main__': main()
