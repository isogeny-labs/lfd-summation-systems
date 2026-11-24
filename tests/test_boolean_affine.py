"""Assignment-level controls for the generic reducer and domain-specific traces."""
from itertools import product
import random

from csd.boolean_affine import reduce_affine, trace_constraints
from csd import ic_dreg as I
from csd import semaev_matched as M


def solutions(equations,n):
    return {v for v in range(1 << n) if all(M.eval_boolean(p,v)==0 for p in equations)}


def test_generic_reduction_exhaustively():
    rng=random.Random(20261124)
    for n in range(1,8):
        for _ in range(60):
            eqs=[rng.sample(range(1 << n),rng.randrange(min(10,1 << n)+1))
                 for _ in range(rng.randrange(1,12))]
            reduced=reduce_affine(eqs,n)
            got=set() if reduced.inconsistent else {reduced.recover(v) for v in
                                                    solutions(reduced.equations,reduced.variables)}
            assert got==solutions(eqs,n)
            if not reduced.inconsistent:
                assert len(got)==len(solutions(reduced.equations,reduced.variables))


def test_hidden_affine_row_and_contradiction():
    # xy+x+z=0 and xy+z+1=0 imply x=1; this is not an input linear row.
    eqs=[[3,1,4],[3,4,0]]
    r=reduce_affine(eqs,3)
    assert r.statistics['eliminated']>=1
    assert {r.recover(v) for v in solutions(r.equations,r.variables)}==solutions(eqs,3)
    assert reduce_affine([[1],[1,0]],1).inconsistent
    assert reduce_affine([],4).variables==4


def test_satisfiable_hidden_affine_maps_after_row_mixing():
    """Force nontrivial bijections, rather than relying on random unsatisfiable systems."""
    rng=random.Random(731)
    for n in range(2,9):
        for _ in range(20):
            free=n//2
            eqs=[]
            assignment=rng.randrange(1<<free)
            for j in range(free,n):
                rhs=[1<<i for i in range(free) if rng.randrange(2)]
                if rng.randrange(2):
                    rhs.append(0)
                if M.eval_boolean(rhs,assignment):
                    assignment |= 1<<j
                eqs.append(set([1<<j]+rhs))
            for _ in range(3):
                p=set(rng.sample(range(1<<n),min(7,1<<n)))
                if M.eval_boolean(p,assignment):
                    p.symmetric_difference_update({0})
                eqs.append(p)
            for _ in range(40):
                a,b=rng.sample(range(len(eqs)),2)
                eqs[a]=eqs[a]^eqs[b]
            reduced=reduce_affine(eqs,n)
            assert not reduced.inconsistent
            assert reduced.statistics['eliminated']>=n-free
            got={reduced.recover(v) for v in solutions(reduced.equations,reduced.variables)}
            assert got==solutions(eqs,n)


def test_trace_constraints_on_rational_point_sums():
    for n in (2,3,4):
        E=I.BinaryCurve(I.binary_field(n),1,3)
        points=M.rational_points(E)
        for P,Q in product(points,repeat=2):
            R=E.add(P,Q)
            if R is not None:
                assert E.F.trace(P[0]^Q[0]^R[0]^E.a)==0


def test_trace_augmentation_and_recovery_against_point_oracle():
    for n in (2,3,4):
        E=I.BinaryCurve(I.binary_field(n),1,1)
        basis=[1,2]
        oracle=M.decomposition_map(E,I.span(basis))
        for R in M.rational_points(E):
            for formulation in ('original','split'):
                eqs,nv=M.build_system(E,basis,R[0],formulation)
                augmented=eqs+trace_constraints(E,basis,R[0],formulation)
                assert solutions(eqs,nv)==solutions(augmented,nv)
                reduction=reduce_affine(augmented,nv)
                masks=set() if reduction.inconsistent else {
                    reduction.recover(v) for v in solutions(reduction.equations,reduction.variables)}
                got={M.decode(v,basis) for v in masks}
                if formulation=='split':
                    got |= M.infinity_solutions(E,basis,R[0])
                assert got==oracle.get(R,set())
