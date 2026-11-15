"""Construct or independently check exact FFD/LFD counterexample certificates."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
import sympy as s

from csd import ic_dreg as I
from csd import ic_fall_certificates as C
from csd import ic_fall_proofs as P
from csd import elim
from csd.field import is_prime,factor


CASES = [
    dict(name="generic-prime-order",p=2689,a=1739,b=1574,target=[797,1184],
         m=2,k=6,exact_lfd=5,exact_ffd=4),
    dict(name="derivative-zero",p=257,a=105,b=196,target=[42,56],
         m=2,k=3,exact_lfd=4,exact_ffd=3),
]


def ordinary_exclusion(curve,target_x,k,degree):
    variables = s.symbols(f"z0:{2*k}")
    chains = []
    for start in (0,k):
        chains.extend(variables[start+j]**2-variables[start+j+1] for j in range(k-1))
        chains.append(variables[start+k-1]**2-1)
    normal = I.pkm_generator(curve,target_x,k)
    generator = sum(c*s.prod(variables[j] for j in range(2*k) if mask >> j & 1)
                    for mask,c in normal.items())
    space = P.degree_closure(chains+[generator],variables,curve.p,degree)
    assert not space.contains(1)
    return dict(degree=degree,ordinary_rank=space.rank,contains_one=False)


def validate_case(case):
    curve = I.PrimeCurve(case['p'],case['a'],case['b'])
    target = tuple(case['target'])
    assert curve.on_curve(target)
    assert (4*curve.a**3+27*curve.b**2) % curve.p
    subgroup = I.subgroup_2power(curve.p,case['k'])
    assert I.prime_solutions(curve,subgroup,target[0]) == []
    points = [curve.lift_x(x) for x in range(curve.p)]
    group_order = 1+sum(1 if q[1]==0 else 2 for q in points if q is not None)
    target_order = I._point_order(curve,target,group_order)
    assert curve.mul(target_order,target) is None
    if case['name']=='generic-prime-order':
        assert group_order==2614 and target_order==1307 and is_prime(target_order)
        assert target[0] not in subgroup
    generator = I.pkm_generator(curve,target[0],case['k'])
    assert P.m2_normalized_homogeneous_ffd(curve.p,curve.a,target[0],case['k']) == case['exact_ffd']
    return curve,generator,dict(group_order=group_order,target_order=target_order,
        target_order_factorization=factor(target_order),target_has_prime_order=is_prime(target_order),
        target_outside_factor_base=target[0] not in subgroup,
        subgroup_grid_size=len(subgroup)**2,subgroup_grid_zeros=0,
        derivative=(3*target[0]**2+curve.a) % curve.p)


def profile(case,generator):
    spaces = [C.chain_degree_closure([generator],2,case['k'],case['p'],d)
              for d in range(case['exact_lfd']+1)]
    first_unit = next(d for d,space in enumerate(spaces) if space.contains_one)
    assert first_unit==case['exact_lfd']
    falls = [d for d in range(1,len(spaces))
             if spaces[d].low_degree_dimension(d-1)>spaces[d-1].rank]
    assert falls[0]==case['exact_ffd'] and falls[-1]==case['exact_lfd']
    rows = [dict(degree=d,standard_columns=len(space.masks),quotient_rank=space.quotient_rank,
                 ordinary_rank=space.rank,lower_degree_dimension=space.low_degree_dimension(d-1),
                 contains_one=space.contains_one,rank_history=list(space.rank_history))
            for d,space in enumerate(spaces)]
    return spaces,rows,falls


def save_certificate(path,space,generator):
    assert C.check_closed_exclusion_certificate(space,[generator])
    dual = C.exclusion_dual(space)
    np.savez_compressed(path,basis=space.basis.astype(np.uint16),
                        masks=np.asarray(space.masks,dtype=np.uint64),
                        dual=dual.astype(np.uint16),prime=space.prime,m=space.m,k=space.k,
                        degree=space.degree_bound)


def check_saved_certificate(path,case,generator):
    with np.load(path,allow_pickle=False) as arrays:
        assert int(arrays['prime'])==case['p'] and int(arrays['m'])==2
        assert int(arrays['k'])==case['k'] and int(arrays['degree'])==case['exact_lfd']-1
        basis = arrays['basis'].astype(np.int64)
        masks = tuple(int(mask) for mask in arrays['masks'])
        dual = arrays['dual'].astype(np.int64)
    rank,pivots,_ = elim.fp_echelon(basis,case['p'])
    assert rank==len(basis) and dual.shape==(len(masks),) and dual[-1]==1
    assert np.all((basis.astype(object) @ dual.astype(object)) % case['p']==0)
    space = C.ChainDegreeClosure(2,case['k'],case['p'],case['exact_lfd']-1,
                                masks,basis,tuple(pivots),())
    assert C.check_closed_exclusion_certificate(space,[generator])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path('experiments/ic-fall-counterexamples'))
    ap.add_argument('--check-only',action='store_true',help='Check saved exclusion matrices, duals, and replay upper derivations.')
    ap.add_argument('--ordinary-check',action='store_true',help='Also verify exclusion with the independent ordinary-ring implementation.')
    args = ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    report_path = args.out/'counterexamples.json'
    previous = json.loads(report_path.read_text()) if args.check_only else None
    results = []
    for case in CASES:
        curve,generator,context = validate_case(case)
        spaces,rows,falls = profile(case,generator)
        filename = case['name']+'-lower-bound.npz'
        certificate = args.out/filename
        if not args.check_only:
            save_certificate(certificate,spaces[-2],generator)
        check_saved_certificate(certificate,case,generator)
        item = dict(**case,**context,generator={str(mask):c for mask,c in sorted(generator.items())},
                    exact_affine_first_strict_fall=falls[0],affine_strict_falls=falls,
                    degree_closures=rows,exclusion_certificate=filename,
                    exclusion_certificate_sha256=hashlib.sha256(certificate.read_bytes()).hexdigest())
        if args.ordinary_check:
            check = ordinary_exclusion(curve,case['target'][0],case['k'],case['exact_lfd']-1)
            assert check['ordinary_rank']==spaces[-2].rank
            item['independent_ordinary_ring_exclusion']=check
        if previous is not None:
            old = next(c for c in previous['counterexamples'] if c['name']==case['name'])
            for key in ('generator','degree_closures','exclusion_certificate_sha256',
                        'exact_ffd','exact_lfd','exact_affine_first_strict_fall'):
                assert old[key]==item[key],key
        results.append(item)
    root = Path(__file__).resolve().parents[1]
    sources = ['docs/IC_FALL_COUNTEREXAMPLES.md','docs/IC_FALL_PROOFS.md','csd/ic_fall_certificates.py',
               'scripts/verify_ic_fall_counterexamples.py','tests/test_ic_fall_certificates.py',
               'csd/ic_fall_proofs.py','csd/ic_dreg.py','csd/elim.py','csd/_elim.c','csd/field.py']
    hashes = {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in sources}
    if previous is not None:
        assert previous['source_sha256']==hashes
    else:
        report = dict(disproved_universal_claims=[
            'Normalized prime-field PKM FFD equals LFD for every instance.',
            'Normalized prime-field PKM LFD is always at most FFD.',
            'Normalized two-summand PKM LFD is always at most 4.'],
            counterexamples=results,source_sha256=hashes,
            versions=dict(numpy=np.__version__,sympy=s.__version__),
            limitations=['Finite exact counterexamples disprove universal statements; they do not prove an unbounded LFD family.',
                         'No claim about a high-probability asymptotic heuristic, other encodings, binary curves, or cryptographic speedup.',
                         'This report checks exact Python arithmetic. A separate Lean audit verifies F257 ordinary Koszul FFD3, affine first strict fall3, and LFD4; the F2689 gap remains outside that audit.'])
        report_path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(mode='check-only' if args.check_only else 'construct-and-check',
                         counterexamples=[{k:c[k] for k in ['name','p','k','exact_ffd','exact_affine_first_strict_fall','exact_lfd','target_order']} for c in results]),indent=2))


if __name__=='__main__':
    main()
