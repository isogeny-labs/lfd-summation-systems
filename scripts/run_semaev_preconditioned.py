#!/usr/bin/env python3
"""Fresh controlled comparison of exact affine preprocessing and existing LRF."""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import hashlib
from itertools import permutations,product
import json
import math
import os
from pathlib import Path
import platform
import random
import resource
import statistics
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from csd import ic_dreg as I
from csd import semaev_preconditioned as P

OUT=ROOT/'experiments/semaev-preconditioned-20261124'
METHODS=['original','split','original_preconditioned','split_preconditioned','lrf']


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')


def generate(out):
    prior=json.loads((ROOT/'experiments/semaev-matched-20261124/instances.json').read_text())
    cases=[dict(c,cohort='prior') for c in prior['cases']]
    missing=[]
    for n,ell in [(7,3),(9,3),(11,3)]:
        for seed in (3,4):
            rng=random.Random((n<<16)+(ell<<8)+seed)
            E=I.BinaryCurve(I.binary_field(n),rng.randrange(2),rng.randrange(1,1<<n))
            basis=[1<<i for i in range(ell)]
            points=P.lifted_points(E,I.span(basis))
            oracle=defaultdict(set)
            for A,B,C in product(points,repeat=3):
                R=E.add(E.add(A,B),C)
                if R is not None:
                    oracle[R].add((A[0],B[0],C[0]))
            targets=P.lifted_points(E,range(1<<n))
            rng.shuffle(targets)
            for cat in ('yes','no'):
                R=next((R for R in targets if bool(oracle.get(R))==(cat=='yes')),None)
                if R is None:
                    missing.append(dict(n=n,ell=ell,seed=seed,category=cat))
                    continue
                cases.append(dict(id=f'n{n}-l{ell}-s{seed}-{cat}',n=n,ell=ell,seed=seed,
                                  a=E.a,b=E.b,modulus=E.F.modulus,target=R,basis=basis,
                                  category=cat,oracle=sorted(oracle.get(R,set())),cohort='heldout'))
    data=dict(cases=cases,missing_heldout_categories=missing,
              prior_instances_sha256=hashlib.sha256((ROOT/'experiments/semaev-matched-20261124/instances.json').read_bytes()).hexdigest())
    write(out/'instances.json',data)
    print(len(cases),'cases;',len(missing),'missing held-out categories')
    return data


def worker(payload):
    case,method=payload['case'],payload['method']
    E=I.BinaryCurve(I.binary_field(case['n']),case['a'],case['b'])
    R=tuple(case['target'])
    if method=='lrf':
        from csd.ic_rootfind import Solver
        start=time.perf_counter()
        solver=Solver(E,case['basis'])
        setup=time.perf_counter()-start
        start=time.perf_counter()
        decompositions,stats,_=solver.lrf(3,R,cap=65536)
        native=time.perf_counter()-start
        start=time.perf_counter()
        raw={t for dec in decompositions for t in permutations([p[0] for p in dec])}
        tuples=raw|P.infinity_solutions(E,case['basis'],R[0])
        recovery=time.perf_counter()-start
        metrics=dict(setup_seconds=setup,preprocess_seconds=0,gb_seconds=0,
                     native_seconds=native,recovery_seconds=recovery,
                     total_seconds=setup+native+recovery,lrf_statistics=stats,
                     projected_solutions=len(tuples),raw_projected_solutions=len(raw),
                     preprocessing_solved=False)
    else:
        pre=method.endswith('_preconditioned')
        formulation=method.split('_')[0]
        started=time.perf_counter()
        try:
            tuples,metrics=P.solve(E,case['basis'],R,formulation,pre,gb_timeout=10)
        except P.GBTimeout:
            print(json.dumps(dict(status='gb_timeout',partial_wall_seconds=time.perf_counter()-started,
                                  peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
            return
    start=time.perf_counter()
    want={tuple(t) for t in case['oracle']}
    if tuples!=want:
        raise AssertionError(dict(case=case['id'],method=method,
                                  missing=sorted(want-tuples),extra=sorted(tuples-want)))
    metrics['verification_seconds']=time.perf_counter()-start
    metrics['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps(dict(status='ok',**metrics)))


def invoke(case,method):
    start=time.perf_counter()
    try:
        r=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],
                         input=json.dumps(dict(case=case,method=method)),text=True,
                         capture_output=True,timeout=40,
                         env={**os.environ,'PYTHONHASHSEED':'0','OPENBLAS_NUM_THREADS':'1'})
        result=json.loads(r.stdout) if r.returncode==0 else dict(status='error',stderr=r.stderr[-6000:])
    except subprocess.TimeoutExpired:
        result=dict(status='process_timeout')
    result['process_wall_seconds']=time.perf_counter()-start
    return result


def benchmark(out):
    import sympy
    from csd.ic_rootfind import _lib
    start=time.perf_counter()
    _lib()
    kernel_ready=time.perf_counter()-start
    cases=json.loads((out/'instances.json').read_text())['cases']
    sources=['csd/boolean_affine.py','csd/semaev_preconditioned.py',
             'scripts/run_semaev_preconditioned.py','csd/semaev_matched.py',
             'csd/ic_dreg.py','csd/gf2_ghs.py','csd/ic_rootfind.py','csd/_rootfind.c']
    write(out/'environment.json',dict(python=sys.version,sympy=sympy.__version__,
          platform=platform.platform(),kernel_cache_ready_seconds=kernel_ready,
          kernel_preparation='readiness check; may reuse compiled cache; excluded once from repeated-target timing',
          source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
          protocol_sha256=hashlib.sha256((out/'PROTOCOL.md').read_bytes()).hexdigest(),
          instances_sha256=hashlib.sha256((out/'instances.json').read_bytes()).hexdigest()))
    jobs=[]
    for ci,case in enumerate(cases):
        for rep in range(3):
            offset=(ci+rep)%len(METHODS)
            order=METHODS[offset:]+METHODS[:offset]
            jobs.extend((case,method,rep) for method in order)
    path=out/'runs.jsonl'
    if path.exists():
        raise FileExistsError('Existing record; use a fresh --out directory')
    start=time.monotonic()
    with path.open('x') as stream:
        for case,method,rep in jobs:
            result=invoke(case,method) if time.monotonic()-start<1800 else dict(status='unrun_budget')
            row=dict(case_id=case['id'],cohort=case['cohort'],method=method,repetition=rep,**result)
            stream.write(json.dumps(row,sort_keys=True)+'\n')
            stream.flush()
            print(case['id'],method,rep,row['status'],round(row.get('total_seconds',0),5),flush=True)
    summarize(out)


def geometric(values):
    return math.exp(statistics.mean(map(math.log,values))) if values else None


def summarize(out):
    cases=json.loads((out/'instances.json').read_text())['cases']
    runs=[json.loads(s) for s in (out/'runs.jsonl').read_text().splitlines()]
    grouped=defaultdict(lambda:defaultdict(list))
    for r in runs:
        grouped[r['case_id']][r['method']].append(r)
    rows=[]
    for case in cases:
        row={k:case[k] for k in ('id','cohort','n','ell','seed','category')}
        for method in METHODS:
            samples=grouped[case['id']].get(method,[])
            info=dict(statuses=[s['status'] for s in samples])
            if len(samples)==3 and all(s['status']=='ok' for s in samples):
                info.update({key:statistics.median(s[key] for s in samples) for key in
                             ('setup_seconds','preprocess_seconds','gb_seconds','recovery_seconds',
                              'total_seconds','peak_rss_kib')})
                info['preprocessing_solved']=samples[0]['preprocessing_solved']
                if 'residual_variables' in samples[0]:
                    info['residual_variables']=samples[0]['residual_variables']
                    info['input_variables']=samples[0]['input']['variables']
            row[method]=info
        if all('total_seconds' in row[m] for m in ('original','split')):
            best=min(row[m]['total_seconds'] for m in ('original','split'))
            for method in ('original_preconditioned','split_preconditioned'):
                if 'total_seconds' not in row[method]:
                    continue
                row[method]['gain_vs_best_gb']=best/row[method]['total_seconds']
                if 'total_seconds' in row['lrf']:
                    row[method]['gain_vs_lrf']=row['lrf']['total_seconds']/row[method]['total_seconds']
        rows.append(row)
    aggregate={}
    for cohort in ('prior','heldout','all'):
        selected=[r for r in rows if cohort=='all' or r['cohort']==cohort]
        aggregate[cohort]={}
        for method in ('original_preconditioned','split_preconditioned'):
            cs=[r for r in selected if 'gain_vs_best_gb' in r[method]]
            gains=[r[method]['gain_vs_best_gb'] for r in cs]
            lrf=[r[method]['gain_vs_lrf'] for r in cs if 'gain_vs_lrf' in r[method]]
            aggregate[cohort][method]=dict(complete_pairs=len(cs),
                geometric_mean_gain=geometric(gains),median_gain=statistics.median(gains) if gains else None,
                wins=sum(g>1 for g in gains),losses=sum(g<1 for g in gains),
                minimum_gain=min(gains,default=None),maximum_gain=max(gains,default=None),
                geometric_mean_gain_vs_lrf=geometric(lrf),
                preprocessing_solved=sum(r[method]['preprocessing_solved'] for r in cs))
    write(out/'summary.json',dict(cases=rows,aggregate=aggregate,
                                 statuses=dict(Counter(r['status'] for r in runs))))
    print(json.dumps(aggregate,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['instances','worker','benchmark','summary'])
    p.add_argument('--out',type=Path,default=OUT)
    a=p.parse_args()
    if a.command=='instances':
        generate(a.out)
    elif a.command=='worker':
        worker(json.load(sys.stdin))
    elif a.command=='benchmark':
        benchmark(a.out)
    else:
        summarize(a.out)


if __name__=='__main__':
    main()
