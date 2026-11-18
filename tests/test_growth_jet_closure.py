"""Protect the explicit hypotheses of the proposed first-derivative route."""
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'experiments/lfd-followups-20261118'


def load_checker():
    path=ROOT/'experiments/lfd-followups-20261117/checker.py'
    spec=importlib.util.spec_from_file_location('jet_independent_checker',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def test_proof_carrying_counterexample_replays():
    from csd.ic_bounds_lift_independent import s3_generator
    record=json.loads((HERE/'certified-falsifier.json').read_text())
    cert=json.loads((HERE/'counterexample-certificate.json').read_text())
    g={(a|(b<<2)):c for (a,b),c in s3_generator(17,6,-7,1,2).items()}
    rows=load_checker().verify(cert,generators=[g],parameters=(2,2,17,3),
        expected_output_sha256=record['checked_output_sha256'])
    row=rows[record['row_index']]
    assert sorted(row.items())==[tuple(t) for t in record['row_terms']]
    x,y=record['rational_sample']
    assert (x*y-x-y-8)%17==0
    plain=sum(c*pow(x,e%4,17)*pow(y,e//4,17) for e,c in row.items())%17
    derivative=sum(c*(e//4)*pow(x,e%4,17)*pow(y,e//4,17) for e,c in row.items())%17
    assert (plain,derivative)==(8,14)


def test_saved_falsifier_source_manifest():
    for path,digest in json.loads((HERE/'certified-falsifier.json').read_text())['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest


def test_cyclic_commutator_requires_wrapping_term():
    # Exact integer check of the concrete coefficient convention used by Lean.
    # Includes characteristic dividing N: no inverse of N is assumed.
    detected=False
    for N in [4,8,16]:
        for s in [1,2]:
            for p in [2,3,17]:
                f=[(i*i+3*i+1)%p for i in range(N)]
                for b in range(N):
                    q=(b+N-s)%N
                    lhs=b*f[q]%p
                    uncorrected=(q*f[q]+s*f[q])%p
                    correction=N*(f[q] if N-s<=q else 0)
                    assert lhs==(uncorrected-correction)%p
                    detected=detected or lhs!=uncorrected
    assert detected
