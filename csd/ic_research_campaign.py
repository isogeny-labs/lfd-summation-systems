"""Reproducible curve/target/base variation and explicitly named encodings."""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import product
import hashlib
import json
import math
import random

from . import ic_dreg as I
from . import ic_fall_certificates as C
from .ic_boolean_closure import closure_profile, check_containing_space
from .ic_binary_ffd import truncated_ffd


def field(n):
    return I.GF2n(1,3) if n==1 else I.binary_field(n)


def s3_boolean(F,b,x,y,z):
    add,mul,sq=I._bp_add,I._bp_mul,I._bp_sq
    pairs=add(add(mul(F,x,y),mul(F,x,z)),mul(F,y,z))
    return add(add(sq(F,pairs),mul(F,mul(F,x,y),z)),I._bp_const(b))


def binary_split(E,basis,m,target_x):
    if m<2 or not basis or I._f2_rank(basis)!=len(basis):
        raise ValueError('Require m>=2 and an independent nonempty factor-base basis.')
    k=len(basis);n=E.F.N
    xs=[I._bp_block(basis,i*k) for i in range(m)]
    auxiliary=[I._bp_block([1<<j for j in range(n)],m*k+i*n) for i in range(m-2)]
    triples=[(xs[0],xs[1],I._bp_const(target_x))] if m==2 else [
        (xs[0],xs[1],auxiliary[0])]+[
        (auxiliary[j-1],xs[j+1],auxiliary[j]) for j in range(1,m-2)]+[
        (auxiliary[-1],xs[-1],I._bp_const(target_x))]
    equations=[]
    for triple in triples:
        f=s3_boolean(E.F,E.b,*triple)
        equations += [sorted(mask for mask,c in f.items() if c>>i&1) for i in range(n)]
    return equations,m*k+(m-2)*n


def _prime_mul(f,g,blocks,k,p):
    output={};N=1<<k;mask=N-1
    for a,c in f.items():
        for b,d in g.items():
            monomial=sum(((((a>>(j*k))&mask)+((b>>(j*k))&mask))%N)<<(j*k) for j in range(blocks))
            output[monomial]=(output.get(monomial,0)+c*d)%p
    return {m:c for m,c in output.items() if c}


def prime_split(E,k,m,target_x):
    if m<2 or k<1 or (E.p-1)%(1<<k):
        raise ValueError('Require subgroup-admissible parameters.')
    blocks=2*m-2;p=E.p
    mul=lambda f,g:_prime_mul(f,g,blocks,k,p)
    def add(*terms):
        out={}
        for scale,term in terms:
            for mask,c in term.items():out[mask]=(out.get(mask,0)+scale*c)%p
        return {mask:c for mask,c in out.items() if c}
    xs=[{1<<(i*k):1} for i in range(m)]
    aux=[{1<<((m+i)*k):1} for i in range(m-2)]
    triples=[(xs[0],xs[1],{0:target_x})] if m==2 else [
        (xs[0],xs[1],aux[0])]+[(aux[j-1],xs[j+1],aux[j]) for j in range(1,m-2)]+[(aux[-1],xs[-1],{0:target_x})]
    generators=[]
    for x,y,z in triples:
        xy,xz,yz=mul(x,y),mul(x,z),mul(y,z)
        xyz=mul(xy,z)
        generators.append(add((1,mul(xy,xy)),(1,mul(xz,xz)),(1,mul(yz,yz)),
            (-2,mul(xyz,add((1,x),(1,y),(1,z)))),
            (-2*E.a,add((1,xy),(1,xz),(1,yz))),
            (-4*E.b,add((1,x),(1,y),(1,z))),(1,{0:E.a*E.a})))
    return generators,blocks


def _point_options(E,x):
    point=E.lift_x(x)
    return set() if point is None else {point,E.neg(point)}


def geometry(E,values,m,target,encoding,aux_values=None):
    """Independent scalar S3 evaluation and rational group-law validation.

    Split encodings require finite partial sums; prime split encodings also
    restrict intermediate x-coordinates to the chosen auxiliary subgroup.
    Counts of coordinate tuples are not independent relation ranks.
    """
    split=encoding=='split-s3'
    scalar=(lambda x,y,z:I.s3_binary_value(E.F,E.b,x,y,z)) if hasattr(E,'F') else (lambda x,y,z:I.s3_prime_value(E,x,y,z))
    if split:
        aux=list(aux_values)
        @lru_cache(None)
        def transition(x,y):
            return frozenset(z for z in aux if scalar(x,y,z)==0)
        counts=Counter(z for x in values for y in values for z in transition(x,y))
        for _ in range(2,m-1):
            new=Counter()
            for z,count in counts.items():
                for x in values:
                    for nxt in transition(z,x):new[nxt]+=count
            counts=new
        assignment_count=sum(count for z,count in counts.items() for x in values if scalar(z,x,target[0])==0)
    polynomial_endpoints=0;rational=0;unrestricted=0;unordered=set()
    for xs in product(values,repeat=m):
        if split:
            states=set(transition(xs[0],xs[1]))
            for x in xs[2:-1]:states={nxt for z in states for nxt in transition(z,x)}
            zero=any(scalar(z,xs[-1],target[0])==0 for z in states)
        elif m==2:zero=scalar(*xs,target[0])==0
        elif m==3:zero=I.s4_binary_value(E.F,E.b,*xs,target[0])==0
        else:raise ValueError('Direct expanded generation supports m=2,3 only.')
        polynomial_endpoints+=int(zero)
        points=[_point_options(E,x) for x in xs]
        def reachable(restrict):
            states={None}
            for i,options in enumerate(points):
                states={E.add(p,q) for p in states for q in options}
                if restrict and 1<=i<m-1:
                    states={p for p in states if p is not None and p[0] in aux_values}
            return target in states
        real=reachable(False)
        restricted=reachable(True) if split else real
        if restricted:
            assert zero,'A finite rational partial-sum path failed its scalar S3 equations.'
            rational+=1;unordered.add(tuple(sorted(xs)))
        unrestricted+=int(real)
    return dict(polynomial_assignment_count=assignment_count if split else polynomial_endpoints,
                polynomial_endpoint_tuples=polynomial_endpoints,
                rational_endpoint_tuples=rational,unrestricted_rational_endpoint_tuples=unrestricted,
                rational_unordered_coordinate_multisets=len(unordered),
                split_loses_rational_endpoint_tuples=unrestricted-rational,
                independent_relation_rank=None)


def binary_cases(seed=20261029):
    rng=random.Random(seed)
    cases=[]
    for n in range(2,7):
        F=field(n)
        curves=[(0,1),(1,1)]+[(rng.randrange(1<<n),rng.randrange(2,1<<n)) for _ in range(2)]
        bases=[]
        for k in range(1,min(3,n)+1):
            candidates=[I.standard_subspace(F,k),I.random_subspace(F,k,rng)]
            invariant=[]
            for element in range(1,1<<n):
                orbit=[];value=element
                for _ in range(n):orbit.append(value);value=F.sq(value)
                candidate=I._f2_basis(orbit)
                if len(candidate)==k:
                    invariant.append(candidate);break
            candidates+=invariant
            for basis in candidates:
                if basis not in bases:bases.append(basis)
        for a,b in curves:
            E=I.BinaryCurve(F,a,b)
            targets=[E.lift_x(x) for x in range(1<<n) if E.lift_x(x) is not None]
            targets=targets if n<=3 else targets[:1]+rng.sample(targets[1:],min(7,len(targets)-1))
            for basis in bases:
                for target in targets:
                    cases.append(dict(kind='binary',n=n,a=a,b=b,basis=basis,target=list(target),m=2,encoding='direct-s3'))
        # A broad proper-ideal inversion graph; one case at each full-field dimension.
        cases.append(dict(kind='binary',n=n,a=0,b=1,basis=[1<<j for j in range(n)],target=[0,1],m=2,encoding='direct-s3',stratum='full-field-inversion'))
    for n in (2,3,4):
        F=field(n)
        for a,b in [(0,1),(1,1),(2,3)]:
            E=I.BinaryCurve(F,a,b)
            targets=[E.lift_x(x) for x in range(1<<n) if E.lift_x(x) is not None][:5]
            for k in (1,2):
                basis=I.standard_subspace(F,k)
                for target in targets:
                    cases.append(dict(kind='binary',n=n,a=a,b=b,basis=basis,target=list(target),m=3,encoding='direct-s4'))
    for n in (1,2,3):
        F=field(n)
        for a,b in [(0,1),(1,1)]:
            E=I.BinaryCurve(F,a,b)
            targets=[E.lift_x(x) for x in range(1<<n) if E.lift_x(x) is not None][:4]
            for m in range(3,10):
                for target in targets:
                    cases.append(dict(kind='binary',n=n,a=a,b=b,basis=[1],target=list(target),m=m,encoding='split-s3'))
    # Deduplicate exact parameter cases while retaining the designated stratum.
    unique={}
    for case in cases:
        key=json.dumps({k:v for k,v in case.items() if k!='stratum'},sort_keys=True)
        unique[key]=case
    return list(unique.values())


def prime_cases(seed=20261029):
    rng=random.Random(seed+1);cases=[]
    for p in (17,41,97,257):
        curves=[(2,2),(5,3),(0,1)]
        for a,b in curves:
            if not (4*a**3+27*b*b)%p:continue
            E=I.PrimeCurve(p,a,b)
            targets=[E.lift_x(x) for x in range(p) if E.lift_x(x) is not None]
            if p>41:targets=targets[:1]+rng.sample(targets[1:],min(15,len(targets)-1))
            for k in range(1,5):
                if (p-1)%(1<<k):continue
                for target in targets:
                    cases.append(dict(kind='prime',prime=p,a=a,b=b,k=k,target=list(target),m=2,encoding='direct-s3'))
    for p in (5,17):
        for a,b in [(2,2),(1,1)]:
            if not (4*a**3+27*b*b)%p:continue
            E=I.PrimeCurve(p,a,b)
            targets=[E.lift_x(x) for x in range(p) if E.lift_x(x) is not None][:4]
            for m in range(3,10):
                for target in targets:
                    cases.append(dict(kind='prime',prime=p,a=a,b=b,k=1,target=list(target),m=m,encoding='split-s3',auxiliary_restriction='same coordinate subgroup; finite partial sums only'))
    return cases


def prepare(case):
    case=dict(case);m=case['m'];target=tuple(case['target']);encoding=case['encoding']
    if case['kind']=='binary':
        E=I.BinaryCurve(field(case['n']),case['a'],case['b']);basis=case['basis'];k=len(basis)
        if encoding=='split-s3':equations,r=binary_split(E,basis,m,target[0])
        else:equations=I.binary_descent(E,[basis]*m,target[0]);r=m*k
        clean=sorted(set(tuple(sorted(row)) for row in equations if row))
        system=dict(kind='binary',variables=r,equations=clean)
        case.update(modulus=E.F.modulus,k=k,variables=r,
                    frobenius_invariant=I.is_frobenius_invariant(E.F,basis),koblitz=E.a in (0,1) and E.b==1)
        values=I.span(basis)
        aux=set(range(1<<case['n'])) if encoding=='split-s3' else None
    else:
        E=I.PrimeCurve(case['prime'],case['a'],case['b']);k=case['k']
        if encoding=='split-s3':generators,blocks=prime_split(E,k,m,target[0])
        else:generators=[I.pkm_generator(E,target[0],k)];blocks=m
        system=dict(kind='prime',prime=E.p,m=blocks,k=k,variables=blocks*k,
                    generators=[sorted(g.items()) for g in generators])
        case.update(variables=blocks*k)
        values=I.subgroup_2power(E.p,k);aux=set(values) if encoding=='split-s3' else None
    assert E.on_curve(target)
    case['geometry']=geometry(E,values,m,target,'split-s3' if encoding=='split-s3' else 'direct',aux)
    case['system_hash']=hashlib.sha256(json.dumps(system,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    case['case_id']=hashlib.sha256(json.dumps(case,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:20]
    return case,system


def prime_profile(system,max_columns=4096,max_cells=8_000_000):
    generators=[dict(g) for g in system['generators']]
    m,k,p=system['m'],system['k'],system['prime'];r=m*k
    maximum=max((mask.bit_count() for g in generators for mask in g),default=0)
    upper=max(r+1,maximum,2)
    spaces=[];falls=[];cap=None
    for D in range(upper+1):
        try:space=C.chain_degree_closure(generators,m,k,p,D,max_columns=max_columns,max_cells=max_cells)
        except ValueError as error:
            cap=str(error);break
        if spaces and space.low_degree_dimension(D-1)>spaces[-1].rank:falls.append(D)
        spaces.append(space)
        if space.contains_one:
            label=dict(status='exact-unit',exact_lfd=D,lfd_lower_bound=D,lfd_upper_bound=D)
            break
    else:label=dict(status='exact-proper',exact_lfd=max(falls,default=0),lfd_lower_bound=max(falls,default=0),lfd_upper_bound=max(falls,default=0))
    if cap:label=dict(status='capped',exact_lfd=None,lfd_lower_bound=max(falls,default=0),lfd_upper_bound=upper)
    label.update(strict_falls=falls,first_affine_strict_fall=min(falls,default=None),
                 quotient_ranks=[s.quotient_rank for s in spaces],ordinary_ranks=[s.rank for s in spaces],
                 scanned_through=len(spaces)-1,cap_reason=cap)
    return label,spaces


def evaluate_system(system):
    if system['kind']=='binary':
        result,spaces=closure_profile(system['equations'],system['variables'],max_columns=4096,max_bytes=96_000_000)
        for space in spaces:
            assert check_containing_space(space,system['equations'],max_bytes=192_000_000)
        # FFD is a separate truncated-syzygy calculation, never inferred from LFD.
        ffd=truncated_ffd(system['equations'],system['variables'],max_columns=4096,max_bytes=64_000_000)
        result['truncated_ffd']=ffd
    else:
        result,spaces=prime_profile(system)
        generators=[dict(g) for g in system['generators']]
        for space in spaces:
            if not space.contains_one:assert C.check_closed_exclusion_certificate(space,generators)
    result['containing_spaces_checked']=True
    return result
