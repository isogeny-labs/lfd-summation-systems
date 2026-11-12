"""Finite evaluation-slice certificates for full ordinary ideal LFD bounds.

The general implication is proved in formal/ICFallDegrees/LocalReductionBound.
These concrete GF(2) matrices are checked in Python, not replayed by Lean.
"""
from __future__ import annotations
from . import elim
from .ic_boolean_closure import (boolean_degree_closure, check_containing_space,
                                ClosureBudgetExceeded)


def evaluation_rank(points, masks, max_bytes=256_000_000):
    if len(points)*max(1,(len(masks)+63)//64)*8 > max_bytes:
        raise ClosureBudgetExceeded('Evaluation matrix exceeds byte cap')
    matrix=elim.pack_rows([[j for j,m in enumerate(masks) if bits&m==m]
                          for bits in points],len(masks))
    return elim.gf2_echelon(matrix,len(masks))[0]


def local_profile(equations, variables, points, *, max_degree=6,
                  max_columns=16384, max_bytes=256_000_000):
    """A supplied set of zeros suffices: completeness is certified by ranks.

    If evaluation on degree <D is onto, and V_D has codimension #points,
    every degree <=D polynomial reduces at budget D to a section of degree
    <D. The local-to-global theorem gives full-ideal LFD <=D.
    """
    points=list(points)
    if len(points)!=len(set(points)) or any(not isinstance(p,int) or not 0<=p<1<<variables for p in points):
        raise ValueError('Evaluation points must be distinct Boolean assignments.')
    if any(sum((p&m)==m for m in row)%2 for p in points for row in equations):
        raise ValueError('A supplied point does not vanish on the inputs.')
    spaces=[];falls=[];checks=[];reason='configured degree cap'
    for D in range(max_degree+1):
        try:
            current=boolean_degree_closure(equations,variables,D,
                max_columns=max_columns,max_bytes=max_bytes)
            assert check_containing_space(current,equations,max_bytes=max_bytes)
            if spaces and current.low_degree_dimension(D-1)>spaces[-1].rank:falls.append(D)
            spaces.append(current)
            if current.contains_one:
                assert not points
                return dict(status='exact-unit',exact_lfd=D,lfd_lower_bound=D,
                    lfd_upper_bound=D,strict_falls=falls,scanned_through=D,
                    quotient_ranks=[s.quotient_rank for s in spaces],
                    point_count=0,local_checks=checks,containing_spaces_checked=True)
            if D>0 and points:
                rank=evaluation_rank(points,[m for m in current.masks if m.bit_count()<D],max_bytes)
                codim=len(current.masks)-current.quotient_rank
                checks.append(dict(degree=D,evaluation_rank=rank,codimension=codim))
                if rank==len(points)==codim:
                    last=max(falls,default=0)
                    return dict(status='exact-proper-local',exact_lfd=last,
                        lfd_lower_bound=last,lfd_upper_bound=last,strict_falls=falls,
                        scanned_through=D,local_upper_bound=D,point_count=len(points),
                        quotient_ranks=[s.quotient_rank for s in spaces],local_checks=checks,
                        containing_spaces_checked=True)
        except ClosureBudgetExceeded as error:
            reason=str(error);break
    return dict(status='capped',exact_lfd=None,lfd_lower_bound=max(falls,default=0),
        lfd_upper_bound=variables+1,strict_falls=falls,scanned_through=len(spaces)-1,
        quotient_ranks=[s.quotient_rank for s in spaces],point_count=len(points),
        local_checks=checks,cap_reason=reason,containing_spaces_checked=True)


def split_points(case, max_paths=100000):
    """Enumerate scalar S3 paths; auxiliary coordinates use the full field."""
    from itertools import product
    from . import ic_dreg as I
    from .ic_research_campaign import field
    F=field(case['n']);basis=case['basis'];k=len(basis);m=case['m'];n=F.N
    if case['encoding']!='split-s3' or m<3:
        raise ValueError('Require a binary split system with m>=3.')
    values=[0]*(1<<k)
    for bits in range(1<<k):
        for j,b in enumerate(basis):
            if bits>>j&1:values[bits]^=b
    transition={(x,y):[z for z in range(1<<n) if I.s3_binary_value(F,case['b'],x,y,z)==0]
                for x in range(1<<n) for y in values}
    paths=[(z,i|(j<<k)|(z<<(m*k))) for i,j in product(range(1<<k),repeat=2)
           for z in transition[values[i],values[j]]]
    for index in range(2,m-1):
        paths=[(z,bits|(i<<(index*k))|(z<<(m*k+(index-1)*n)))
               for previous,bits in paths for i in range(1<<k)
               for z in transition[previous,values[i]]]
        if len(paths)>max_paths:raise ClosureBudgetExceeded('Scalar path enumeration cap')
    points=[bits|(i<<((m-1)*k)) for previous,bits in paths for i in range(1<<k)
            if case['target'][0] in transition[previous,values[i]]]
    assert len(points)==len(set(points))
    assert len(points)==case['geometry']['polynomial_assignment_count']
    return sorted(points)


def zero_coordinate_profile(equations, variables, points, *, budget=4,
                            max_columns=32768, max_bytes=1_000_000_000):
    """Transport an upper bound through coordinate zeros derived at budget B.

    Telescoping gives p-p|_(zeros=0) at max(B,deg p). Each projected input
    is available at max(B,its original degree). Lifting every reduced
    closure step therefore costs max(B,D), and full-ideal LFD is at most
    max(B,LFD(reduced)). This is a written transport proof, not a Lean instance.
    """
    if any(m.bit_count()>budget for row in equations for m in row):
        raise ValueError('Coordinate transport requires every input degree <= budget.')
    space=boolean_degree_closure(equations,variables,budget,
                                max_columns=max_columns,max_bytes=max_bytes)
    assert check_containing_space(space,equations,max_bytes=max_bytes)
    zero=[j for j in range(variables) if space.contains([1<<j])]
    keep=[j for j in range(variables) if j not in zero]
    if not zero or not keep:return None
    killed=sum(1<<j for j in zero)
    compress=lambda mask:sum(((mask>>j)&1)<<i for i,j in enumerate(keep))
    projected=[[compress(m) for m in row if not m&killed] for row in equations]
    assert all(not p&killed for p in points)
    reduced=local_profile(projected,len(keep),[compress(p) for p in points],
                          max_degree=budget+1,max_columns=max_columns,max_bytes=max_bytes)
    if reduced['exact_lfd'] is None:return None
    return dict(lfd_upper_bound=max(budget,reduced['exact_lfd']),zero_coordinates=zero,
                retained_coordinates=keep,preprocessing_budget=budget,
                reduced_variables=len(keep),reduced_result=reduced)
