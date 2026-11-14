from dataclasses import replace

import numpy as np
import pytest
import sympy as s

from csd import ic_dreg as I
from csd import ic_fall_proofs as P
from csd import ic_fall_certificates as C
from csd.field import is_prime


def _ordinary_input(curve,target_x,k):
    variables = s.symbols(f"z0:{2*k}")
    chains = []
    for start in (0,k):
        chains.extend(variables[start+j]**2-variables[start+j+1] for j in range(k-1))
        chains.append(variables[start+k-1]**2-1)
    normal = I.pkm_generator(curve,target_x,k)
    generator = sum(c*s.prod(variables[j] for j in range(2*k) if mask >> j & 1)
                    for mask,c in normal.items())
    return variables,chains,generator,normal


@pytest.mark.parametrize("k,degree", [(1,d) for d in range(5)]+[(2,d) for d in range(6)])
def test_exact_quotient_matches_independent_ordinary_closure(k,degree):
    curve = I.PrimeCurve(17,2,2)
    variables,chains,generator,normal = _ordinary_input(curve,3,k)
    ordinary = P.degree_closure(chains+[generator],variables,17,degree)
    quotient = C.chain_degree_closure([normal],2,k,17,degree)
    assert ordinary.rank == quotient.rank
    assert ordinary.contains(1) == quotient.contains_one
    for d in range(degree+1):
        assert ordinary.low_degree_dimension(d) == quotient.low_degree_dimension(d)


def test_original_input_degree_is_honored_after_normalization():
    x,y = s.symbols("x y")
    for degree in range(5):
        ordinary = P.degree_closure([x*x-1,y*y-1,x**4+2],(x,y),17,degree)
        quotient = C.chain_degree_closure([{0:3}],2,1,17,degree,input_degrees=[4])
        assert ordinary.rank == quotient.rank
        assert ordinary.contains(1) == quotient.contains_one == (degree==4)


def test_generic_counterexample_exact_ffd4_lfd5_with_prime_order_target():
    curve = I.PrimeCurve(2689,1739,1574)
    target,k = (797,1184),6
    assert curve.on_curve(target)
    assert (4*curve.a**3+27*curve.b**2) % curve.p
    assert is_prime(1307) and curve.mul(1307,target) is None
    points = [curve.lift_x(x) for x in range(curve.p)]
    assert 1+sum(1 if q[1]==0 else 2 for q in points if q is not None) == 2614
    subgroup = I.subgroup_2power(curve.p,k)
    assert target[0] not in subgroup
    assert I.prime_solutions(curve,subgroup,target[0]) == []
    variables,chains,generator,normal = _ordinary_input(curve,target[0],k)
    assert normal == {2:605,128:605,65:690,66:1095,129:1095,
                      1:2164,64:2164,0:1447,130:1}
    assert P.m2_normalized_homogeneous_ffd(curve.p,curve.a,target[0],k) == 4
    spaces = [C.chain_degree_closure([normal],2,k,curve.p,d) for d in range(6)]
    assert [v.quotient_rank for v in spaces] == [0,0,1,13,605,1586]
    assert [v.rank for v in spaces] == [0,0,13,169,1631,6188]
    assert not spaces[4].contains_one and spaces[5].contains_one
    falls = [d for d in range(1,6)
             if spaces[d].low_degree_dimension(d-1) > spaces[d-1].rank]
    assert falls == [4,5]  # Affine first fall also 4; unit ideal gives exact LFD 5.
    assert C.check_closed_exclusion_certificate(spaces[4],[normal])
    dual = C.exclusion_dual(spaces[4])
    assert dual[-1] == 1
    assert np.all((spaces[4].basis.astype(object) @ dual.astype(object)) % curve.p == 0)
    # A separate implementation in the ordinary ring confirms the lower bound.
    ordinary = P.degree_closure(chains+[generator],variables,curve.p,4)
    assert ordinary.rank == 1631 and not ordinary.contains(1)


def test_degenerate_counterexample_exact_ffd3_lfd4():
    curve = I.PrimeCurve(257,105,196)
    variables,chains,generator,normal = _ordinary_input(curve,42,3)
    assert curve.on_curve((42,56))
    assert (3*42**2+curve.a) % curve.p == 0
    assert P.m2_normalized_homogeneous_ffd(curve.p,curve.a,42,3) == 3
    for degree,rank,unit in [(2,7,False),(3,58,False),(4,210,True)]:
        ordinary = P.degree_closure(chains+[generator],variables,curve.p,degree)
        quotient = C.chain_degree_closure([normal],2,3,curve.p,degree)
        assert ordinary.rank == quotient.rank == rank
        assert ordinary.contains(1) == quotient.contains_one == unit


def test_lower_certificate_rejects_missing_input_and_unstable_or_unit_spaces():
    curve = I.PrimeCurve(17,2,2)
    normal = I.pkm_generator(curve,3,2)
    low = C.chain_degree_closure([normal],2,2,17,2)
    assert C.check_closed_exclusion_certificate(low,[normal])
    assert not C.check_closed_exclusion_certificate(low,[{0:1}])
    assert not C.check_closed_exclusion_certificate(low,[normal],input_degrees=[])
    wrong_columns = replace(low,masks=low.masks[:-1],basis=low.basis[:,:-1])
    assert not C.check_closed_exclusion_certificate(wrong_columns,[normal])
    # The degree-2 row alone fails closure at degree 3.
    high = C.chain_degree_closure([normal],2,2,17,3)
    index = {mask:j for j,mask in enumerate(high.masks)}
    only_generator = np.zeros((1,len(high.masks)),dtype=np.int64)
    for mask,c in normal.items():
        only_generator[0,index[mask]] = c
    unstable = replace(high,basis=only_generator,pivots=())
    assert not C.check_closed_exclusion_certificate(unstable,[normal])
    unit = C.chain_degree_closure([{0:1}],2,2,17,3)
    assert not C.check_closed_exclusion_certificate(unit,[{0:1}])
    with pytest.raises(ValueError,match="contains it"):
        C.exclusion_dual(unit)


def test_empty_zero_and_constant_generators_and_resource_guards():
    for generators in ([],[{}],[{0:17}]):
        space = C.chain_degree_closure(generators,2,2,17,3)
        assert space.quotient_rank == 0 and space.rank == 20
        assert not space.contains_one
    unit = C.chain_degree_closure([{0:3}],2,2,17,1)
    assert unit.contains_one and unit.rank == 5
    assert unit.contains_normal_form({1:2,2:5})
    assert not unit.contains_normal_form({3:1})
    for args in [(2,0,17,3),(2,2,15,3),(2,5,17,3),(0,2,17,3)]:
        with pytest.raises(ValueError):
            C.chain_degree_closure([], *args)
    with pytest.raises(ValueError,match="masks"):
        C.chain_degree_closure([{16:1}],2,2,17,3)
    with pytest.raises(ValueError,match="Input degrees"):
        C.chain_degree_closure([{3:1}],2,2,17,3,input_degrees=[1])
    with pytest.raises(ValueError,match="column budget"):
        C.chain_degree_closure([],2,6,2689,5,max_columns=100)
    with pytest.raises(ValueError,match="budget"):
        C.chain_degree_closure([{1:1}],2,2,17,3,max_cells=1)
