"""Regression checks for the high-prime arithmetic fix and proof-carrying rows."""
import copy
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from csd.ic_bounds_lift_independent import Refuse, max_inner, mm, rref
from csd.ic_bounds_split_m import SplitSpace
from csd.ic_bounds_m3_batched import BatchedClosure

HERE=Path(__file__).resolve().parents[1]/'experiments/lfd-followups-20261117'

def module(name):
    spec=importlib.util.spec_from_file_location('lfd_'+name,HERE/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

@pytest.mark.parametrize('p',[65537,163841,16777213])
def test_rref_exact_when_inner_block_exceeds_float_limit(p):
    M=np.zeros((64,65));M[:63,:63]=np.eye(63);M[:63,64]=p-1
    M[63,:63]=p-1;M[63,63]=1;M[63,64]=p-1
    R,piv=rref(M,p,np)
    assert len(piv)==64
    assert int(R[-1,64])==(p-1-63*(p-1)**2)%p

@pytest.mark.parametrize('inner',[0,-1])
def test_unsafe_accumulation_override_is_refused(inner):
    with pytest.raises(Refuse):mm(np.ones((1,1)),np.ones((1,1)),16777213,np,inner)

def test_derivation_checker_binds_inputs_and_outputs():
    producer,checker=module('derivations'),module('checker')
    gens=[{3:1,0:2}];params=(1,3,17,3)
    cert=producer.produce(gens,*params[:3],params[3])
    assert checker.verify(cert,generators=gens,parameters=params,expected_output_sha256=cert["row_sha256"])
    bad=copy.deepcopy(cert);bad['nodes'].append({'op':'lin','terms':[[len(bad['nodes']),1]]})
    with pytest.raises(ValueError):checker.verify(bad,generators=gens,parameters=params,expected_output_sha256=cert["row_sha256"])
    with pytest.raises(ValueError):checker.verify(cert,generators=[{0:1}],parameters=params,expected_output_sha256=cert["row_sha256"])
    bad=copy.deepcopy(cert);bad['rows'][0]=[[0,1]]
    with pytest.raises(ValueError):checker.verify(bad,generators=gens,parameters=params,expected_output_sha256=cert["row_sha256"])

def test_carry_cannot_hide_illegal_unreduced_degree():
    cert=dict(m=1,k=3,p=17,D=2,generators=[[[3,1]]],nodes=[{'op':'input','index':0},
         {'op':'shift','source':0,'variable':0}],outputs=[],rows=[],row_sha256='')
    with pytest.raises(ValueError,match='unreduced'):
        module('checker').verify(cert,generators=[{3:1}],parameters=(1,3,17,2),expected_output_sha256="")

@pytest.mark.parametrize('D',[3,4])
def test_tiled_queue_matches_dense_closure(D):
    sp=SplitSpace(3,2,D,1);p=17
    # Includes auxiliary powers and a chain carry in the shifted polynomials.
    g={(3,(0,)):1,(0,(1,)):2,(0,(0,)):3}
    rows=[sp.row(g,p)]
    dense=BatchedClosure(sp,p,batch=16);dense.add_inputs(rows);dense.saturate()
    tiled=module('tiled_split').TiledClosure(sp,p,batch=16,tile_rows=3,tile_cols=7)
    tiled.add_inputs(rows);assert tiled.saturate()
    assert tiled.rank==dense.rank and tiled.contains_one==dense.contains_one
    assert [tiled.low_dimension(d) for d in range(D+1)]==[dense.low_dimension(d) for d in range(D+1)]
    assert tiled.verify_closed() and tiled.contains(rows[0])


def test_requested_product_chunk_is_capped_to_exact_bound():
    p=16777213
    A=np.full((1,65),p-1,dtype=np.float64)
    B=np.full((65,1),p-1,dtype=np.float64)
    assert int(mm(A,B,p,np,1000)[0,0])==65


def test_saved_certificate_replays_against_reconstructed_curve_input():
    import json
    from csd.ic_bounds_lift_independent import s3_generator
    cert=json.loads((HERE/'certificate-k2-D4.json').read_text())
    gen={(i|(j<<2)):c for (i,j),c in s3_generator(65537,1,0,4271,2).items()}
    rows=module('checker').verify(cert,generators=[gen],parameters=(2,2,65537,4),expected_output_sha256="96dcd59adb15c96fe3fb77202647e3b93640b71885b7846e8f7f3eba0e7e322a")
    assert len(rows)==16


def test_followup_records_bind_current_sources():
    import hashlib
    import json
    root=HERE.parents[1]
    for name in ['validation.json','split-tiled-controls.json']:
        record=json.loads((HERE/name).read_text())
        assert record['status']=='passed'
        for path,digest in record['source_sha256'].items():
            assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest


def test_exceptional_zero_sets_match_independent_involution():
    # Asserts the polynomial square identity and equality with full-grid enumeration.
    module('check_exceptional_points').main()


def test_certificate_cannot_substitute_an_empty_output_set():
    import hashlib
    producer,checker=module("derivations"),module("checker")
    gens=[{3:1,0:2}];params=(1,3,17,3)
    cert=producer.produce(gens,*params[:3],params[3])
    expected=cert["row_sha256"]
    cert.update(outputs=[],rows=[],row_sha256=hashlib.sha256(b"[]").hexdigest())
    with pytest.raises(ValueError,match="external target"):
        checker.verify(cert,generators=gens,parameters=params,expected_output_sha256=expected)
