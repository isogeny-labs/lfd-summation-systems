#!/usr/bin/env python3
"""Conservative charged-time intervals; never treat a timeout as exact timing."""
from collections import defaultdict
import json
import math
from pathlib import Path
import statistics

OUT=Path(__file__).resolve().parents[1]/'experiments/semaev-preconditioned-20261124'


def median_bounds(rows):
    intervals=[]
    for r in rows:
        if r['status']=='ok':
            t=r['total_seconds']
            intervals.append((t,t))
        elif r['status']=='gb_timeout':
            intervals.append((10.0,math.inf))
        else:
            intervals.append((0.0,math.inf))
    intervals += [(0.0,math.inf)]*(3-len(intervals))
    if len(intervals)!=3:
        raise ValueError('expected three unique repetitions')
    return statistics.median(a for a,b in intervals),statistics.median(b for a,b in intervals)


def gain_bounds(base,candidate):
    bl,bu=base
    cl,cu=candidate
    lower=0.0 if math.isinf(cu) else bl/cu
    upper=math.inf if cl==0 else bu/cl
    return lower,upper


def geometric(values):
    if any(v==0 for v in values):
        return 0.0
    return math.exp(statistics.mean(map(math.log,values))) if values else None


def interval(a,b):
    return dict(lower=a,upper=None if math.isinf(b) else b)


def analyze(out=OUT):
    inputs=json.loads((out/'instances.json').read_text())['cases']
    runs=[json.loads(s) for s in (out/'runs.jsonl').read_text().splitlines()]
    groups=defaultdict(lambda:defaultdict(list))
    for r in runs:
        groups[r['case_id']][r['method']].append(r)
    result=[]
    for c in inputs:
        methods=groups[c['id']]
        bounds={m:median_bounds(methods.get(m,[])) for m in
                ('original','split','original_preconditioned')}
        best=(min(bounds[m][0] for m in ('original','split')),
              min(bounds[m][1] for m in ('original','split')))
        ratio=gain_bounds(best,bounds['original_preconditioned'])
        result.append(dict(id=c['id'],cohort=c['cohort'],
                      baseline=interval(*best),candidate=interval(*bounds['original_preconditioned']),
                      gain=interval(*ratio),
                      primary_complete=len(methods.get('original_preconditioned',[]))==3 and
                      all(r['status']=='ok' for r in methods.get('original_preconditioned',[])),
                      baseline_censored=any(r['status']!='ok' for m in ('original','split')
                                            for r in methods.get(m,[]))))
    aggregate={}
    for cohort in ('prior','heldout','all'):
        rows=[r for r in result if cohort=='all' or r['cohort']==cohort]
        aggregate[cohort]=dict(cases=len(rows),
            geometric_mean_gain_lower_bound=geometric([r['gain']['lower'] for r in rows]),
            cases_with_censored_baseline=sum(r['baseline_censored'] for r in rows),
            primary_complete=sum(r['primary_complete'] for r in rows),
            exactly_resolved_best_baselines=sum(r['baseline']['lower']==r['baseline']['upper'] for r in rows))
    target=all(aggregate[c]['primary_complete']==aggregate[c]['cases'] and
               aggregate[c]['geometric_mean_gain_lower_bound']>=2 for c in ('prior','heldout'))
    return dict(analysis='secondary conservative intervals; frozen complete-pair gate reported separately',
                target_supported=target,cases=result,aggregate=aggregate)


def main():
    result=analyze()
    (OUT/'timing-bounds.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(target_supported=result['target_supported'],aggregate=result['aggregate']),indent=2))


if __name__=='__main__':
    main()
