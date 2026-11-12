"""Exact binary-S3 experiments on spaces outside the scaled-subfield family."""
from __future__ import annotations
import hashlib,json
from . import ic_dreg as I


def canonical_basis(vectors):
    rows={}
    for value in vectors:
        for pivot in sorted(rows,reverse=True):
            if value>>pivot&1:value^=rows[pivot]
        if not value:continue
        pivot=value.bit_length()-1
        rows={p:(r^value if r>>pivot&1 else r) for p,r in rows.items()}
        rows[pivot]=value
    return [rows[p] for p in sorted(rows,reverse=True)]


def scaled_subfield(F,basis):
    k=len(basis)
    if not k or F.N%k:return False
    inverse=F.inv(basis[0])
    normalized=[F.mul(inverse,b) for b in basis]
    return all(F.frob(b,k)==b for b in normalized)


def boolean_points(equations,k,ell=None):
    ell=k if ell is None else ell;out=[];mask=(1<<k)-1
    if any((m&mask).bit_count()>1 or (m>>k).bit_count()>1 for eq in equations for m in eq):
        raise ValueError('Require a bilinear Boolean coordinate presentation.')
    for left in range(1<<k):
        images=[0]*ell;rhs=0
        for j,eq in enumerate(equations):
            for m in eq:
                a,b=m&mask,m>>k
                if left&a!=a:continue
                if b:images[b.bit_length()-1]^=1<<j
                else:rhs^=1<<j
        solution=I._f2_solve(images,rhs,len(equations))
        if solution is None:continue
        particular,kernel=solution
        if len(kernel)>1:raise ValueError('A scalar-S3 fiber has more than two roots.')
        out.extend(left|((particular^free)<<k) for free in I.span(kernel))
    return sorted(out)


def prepare(E,basis,target,target_class,planted=None):
    F=E.F;k=len(basis)
    if not k or any(not isinstance(v,int) or not 0<v<1<<F.N for v in basis):
        raise ValueError('Require nonzero field elements in the supplied basis.')
    basis=canonical_basis(basis)
    if len(basis)!=k:raise ValueError('The coordinate basis must be linearly independent.')
    assert E.on_curve(target) and not scaled_subfield(F,basis)
    equations=sorted(set(tuple(sorted(row)) for row in I.binary_descent(E,[basis,basis],target[0]) if row))
    system=dict(variables=2*k,equations=equations)
    key=hashlib.sha256(json.dumps(system,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    points=boolean_points(equations,k);values=I.span(basis);rational=0;unordered=set()
    for bits in points:
        x,y=values[bits&((1<<k)-1)],values[bits>>k]
        assert I.s3_binary_value(F,E.b,x,y,target[0])==0
        P=E.lift_x(x)
        if P is None:continue
        if any((Q:=E.add(target,E.neg(P0))) is not None and Q[0]==y for P0 in {P,E.neg(P)}):
            rational+=1;unordered.add(tuple(sorted((x,y))))
    if planted:
        P,Q=map(tuple,planted)
        assert E.add(P,Q)==target and P[0] in values and Q[0] in values
        assert rational>0
    case=dict(n=F.N,modulus=F.modulus,k=k,a=E.a,b=E.b,basis=basis,target=list(target),
        target_class=target_class,planted=planted,scaled_subfield=False,
        frobenius_invariant=I.is_frobenius_invariant(F,basis),
        polynomial_tuples=len(points),rational_ordered_tuples=rational,
        rational_unordered_multisets=len(unordered),independent_relation_rank=None,system_hash=key)
    case['case_id']=hashlib.sha256(json.dumps(case,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:20]
    return case,dict(system_hash=key,**system,coordinate_solutions=points)
