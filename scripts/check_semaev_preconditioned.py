#!/usr/bin/env python3
"""Check frozen-source records and independently enumerate all PDP oracle tuples.

No timing replay or performance threshold is asserted by this integrity check.
"""
from collections import Counter
import contextlib
import hashlib
from itertools import product
import io
import json
import math
from pathlib import Path
import sys
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_semaev_preconditioned as RUN
from analyze_semaev_timing_bounds import analyze
from csd import ic_dreg as I


def main():
    out=RUN.OUT
    env=json.loads((out/'environment.json').read_text())
    replaced=[]
    for path,digest in env['source_sha256'].items():
        current=RUN.ROOT/path
        if hashlib.sha256(current.read_bytes()).hexdigest()!=digest:
            # The native rootfinder was subsequently extended. Preserve the
            # measured source bytes rather than rewriting the frozen record.
            assert path in {'csd/_rootfind.c','csd/ic_rootfind.py'},path
            snapshot=out/'frozen-source'/path
            assert hashlib.sha256(snapshot.read_bytes()).hexdigest()==digest,path
            replaced.append(path)
    for path,key in [('PROTOCOL.md','protocol_sha256'),('instances.json','instances_sha256')]:
        assert hashlib.sha256((out/path).read_bytes()).hexdigest()==env[key],path
    inputs=json.loads((out/'instances.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        target=Path(tmp)
        with contextlib.redirect_stdout(io.StringIO()):
            regenerated=RUN.generate(target)
        assert json.loads(json.dumps(regenerated))==inputs
        # Enumerate y for each x in V, independently of lift_x and LRF.
        for c in inputs['cases']:
            E=I.BinaryCurve(I.binary_field(c['n']),c['a'],c['b'])
            fb=[(x,y) for x in I.span(c['basis']) for y in range(1<<c['n'])
                if E.on_curve((x,y))]
            R=tuple(c['target'])
            assert E.on_curve(R)
            want={(A[0],B[0],C[0]) for A,B,C in product(fb,repeat=3)
                  if E.add(E.add(A,B),C)==R}
            assert want=={tuple(t) for t in c['oracle']},c['id']
        runs=[json.loads(s) for s in (out/'runs.jsonl').read_text().splitlines()]
        keys=[(r['case_id'],r['method'],r['repetition']) for r in runs]
        expected={(c['id'],m,rep) for c in inputs['cases'] for m in RUN.METHODS for rep in range(3)}
        assert len(keys)==len(set(keys)) and set(keys)==expected
        oracle_counts={c['id']:len(c['oracle']) for c in inputs['cases']}
        for r in runs:
            if r['status']!='ok':
                continue
            assert r['projected_solutions']==oracle_counts[r['case_id']]
            phase=sum(r[k] for k in ('setup_seconds','preprocess_seconds','gb_seconds','recovery_seconds'))
            phase+=r.get('native_seconds',0)
            assert math.isclose(phase,r['total_seconds'],rel_tol=1e-10,abs_tol=1e-10)
            if 'residual_variables' in r:
                assert r['residual_variables']<=r['input']['variables']
                if r['method'].endswith('_preconditioned'):
                    info=r['preprocessing']
                    assert len(info['row_ranks'])==len(info['row_columns'])==info['row_rounds']
                    assert all(rank<=cols for rank,cols in zip(info['row_ranks'],info['row_columns']))
        (target/'instances.json').write_bytes((out/'instances.json').read_bytes())
        (target/'runs.jsonl').write_bytes((out/'runs.jsonl').read_bytes())
        with contextlib.redirect_stdout(io.StringIO()):
            RUN.summarize(target)
        assert json.loads((target/'summary.json').read_text())==json.loads((out/'summary.json').read_text())
        assert analyze(out)==json.loads((out/'timing-bounds.json').read_text())
    print(json.dumps(dict(status='passed',cases=len(inputs['cases']),records=len(runs),
                          statuses=dict(Counter(r['status'] for r in runs)),
                          frozen_source_snapshots=sorted(replaced),
                          scope='independent point-oracle scan, source/record/summary integrity; not timing replay')))


if __name__=='__main__':
    main()
