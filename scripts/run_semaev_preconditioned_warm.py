#!/usr/bin/env python3
"""Same held-out candidate and cases, confirmed in one warm SymPy process."""
from collections import Counter,defaultdict
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from csd import ic_dreg as I
from csd import semaev_preconditioned as P
from analyze_semaev_timing_bounds import median_bounds,gain_bounds,geometric,interval

PARENT=ROOT/'experiments/semaev-preconditioned-20261124'
OUT=PARENT/'warm'
METHODS=['original','split','original_preconditioned']


def write(path,value):
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def run():
    import sympy as sp
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/'runs.jsonl'
    if path.exists():
        raise FileExistsError('Refusing to overwrite warm measurements')
    cases=[c for c in json.loads((PARENT/'instances.json').read_text())['cases'] if c['cohort']=='heldout']
    write(OUT/'instances.json',dict(cases=cases))
    x=sp.symbols('warmup_x')
    sp.groebner([x*x+x],x,modulus=2,method='buchberger',order='grevlex')
    sources=['csd/boolean_affine.py','csd/semaev_preconditioned.py',
             'scripts/run_semaev_preconditioned_warm.py','scripts/analyze_semaev_timing_bounds.py']
    write(OUT/'environment.json',dict(python=sys.version,sympy=sp.__version__,platform=platform.platform(),
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        protocol_sha256=hashlib.sha256((PARENT/'WARM-PROTOCOL.md').read_bytes()).hexdigest(),
        instances_sha256=hashlib.sha256((OUT/'instances.json').read_bytes()).hexdigest(),
        mode='one process; trivial SymPy Boolean GB warmup before timed calls'))
    begin=time.monotonic()
    with path.open('x') as stream:
        for ci,case in enumerate(cases):
            E=I.BinaryCurve(I.binary_field(case['n']),case['a'],case['b'])
            want={tuple(t) for t in case['oracle']}
            for rep in range(3):
                offset=(ci+rep)%3
                order=METHODS[offset:]+METHODS[:offset]
                for method in order:
                    row=dict(case_id=case['id'],method=method,repetition=rep)
                    if time.monotonic()-begin>=900:
                        row['status']='unrun_budget'
                    else:
                        try:
                            tuples,metrics=P.solve(E,case['basis'],tuple(case['target']),
                                method.split('_')[0],method.endswith('_preconditioned'),gb_timeout=10)
                            if tuples!=want:
                                raise AssertionError(dict(missing=sorted(want-tuples),extra=sorted(tuples-want)))
                            row.update(status='ok',**metrics)
                        except P.GBTimeout:
                            row['status']='gb_timeout'
                        except Exception as error:
                            row.update(status='error',error=repr(error))
                    stream.write(json.dumps(row,sort_keys=True)+'\n')
                    stream.flush()
                    print(case['id'],method,rep,row['status'],round(row.get('total_seconds',0),5),flush=True)
    summarize()


def summary_data():
    cases=json.loads((OUT/'instances.json').read_text())['cases']
    runs=[json.loads(s) for s in (OUT/'runs.jsonl').read_text().splitlines()]
    grouped=defaultdict(lambda:defaultdict(list))
    for r in runs:
        grouped[r['case_id']][r['method']].append(r)
    result=[]
    for c in cases:
        methods=grouped[c['id']]
        bounds={m:median_bounds(methods.get(m,[])) for m in METHODS}
        best=(min(bounds[m][0] for m in ('original','split')),
              min(bounds[m][1] for m in ('original','split')))
        gain=gain_bounds(best,bounds['original_preconditioned'])
        complete=len(methods['original_preconditioned'])==3 and all(
            r['status']=='ok' for r in methods['original_preconditioned'])
        result.append(dict(id=c['id'],baseline=interval(*best),candidate=interval(*bounds['original_preconditioned']),
                           gain=interval(*gain),primary_complete=complete))
    lower=geometric([r['gain']['lower'] for r in result])
    return dict(cases=result,statuses=dict(Counter(r['status'] for r in runs)),
                geometric_mean_gain_lower_bound=lower,
                primary_complete_cases=sum(r['primary_complete'] for r in result),
                target_supported=len(result)==12 and all(r['primary_complete'] for r in result) and lower>=2)


def summarize():
    data=summary_data()
    write(OUT/'summary.json',data)
    print(json.dumps({k:v for k,v in data.items() if k!='cases'},indent=2))


def check():
    env=json.loads((OUT/'environment.json').read_text())
    for p,digest in env['source_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest,p
    assert hashlib.sha256((PARENT/'WARM-PROTOCOL.md').read_bytes()).hexdigest()==env['protocol_sha256']
    assert hashlib.sha256((OUT/'instances.json').read_bytes()).hexdigest()==env['instances_sha256']
    cs=json.loads((OUT/'instances.json').read_text())['cases']
    expected_cs=[c for c in json.loads((PARENT/'instances.json').read_text())['cases'] if c['cohort']=='heldout']
    assert cs==expected_cs
    runs=[json.loads(s) for s in (OUT/'runs.jsonl').read_text().splitlines()]
    keys=[(r['case_id'],r['method'],r['repetition']) for r in runs]
    expected={(c['id'],m,rep) for c in cs for m in METHODS for rep in range(3)}
    assert len(keys)==len(set(keys)) and set(keys)==expected
    counts={c['id']:len(c['oracle']) for c in cs}
    for r in runs:
        if r['status']=='ok':
            assert r['projected_solutions']==counts[r['case_id']]
    assert summary_data()==json.loads((OUT/'summary.json').read_text())
    print(json.dumps(dict(status='passed',records=len(runs),scope='source/record/summary integrity; not timing replay')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['run','summary','check'])
    a=p.parse_args()
    {'run':run,'summary':summarize,'check':check}[a.command]()
