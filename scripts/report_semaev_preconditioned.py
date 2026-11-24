#!/usr/bin/env python3
"""Render the prespecified candidate against the optimistic GB baseline and LRF."""
import json
import math
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parents[1]/'experiments/semaev-preconditioned-20261124'


def gm(values):
    return math.exp(statistics.mean(map(math.log,values)))


def main():
    data=json.loads((OUT/'summary.json').read_text())
    bounds=json.loads((OUT/'timing-bounds.json').read_text())
    warm=json.loads((OUT/'warm/summary.json').read_text())
    by_id={c['id']:c for c in bounds['cases']}
    method='original_preconditioned'
    cases=[c for c in data['cases'] if 'total_seconds' in c[method] and
           by_id[c['id']]['baseline']['upper']==by_id[c['id']]['baseline']['lower']]
    fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for cohort,color,marker in [('prior','#1976a3','o'),('heldout','#ca6b24','s')]:
        cs=[c for c in cases if c['cohort']==cohort]
        x=[by_id[c['id']]['baseline']['lower'] for c in cs]
        y=[c[method]['total_seconds'] for c in cs]
        ax.scatter(x,y,c=color,marker=marker,label=f'{cohort} ({len(cs)})',s=32,alpha=.8)
    values=[v for c in cases for v in (by_id[c['id']]['baseline']['lower'],c[method]['total_seconds'])]
    lo,hi=min(values)/1.5,max(values)*1.5
    ax.plot([lo,hi],[lo,hi],color='#555555',linestyle='--',label='Equal charged time')
    ax.plot([lo,hi],[lo/2,hi/2],color='#999999',linestyle=':',label='2x improvement')
    ax.set(xscale='log',yscale='log',xlabel='Faster original/split baseline (seconds)',
           ylabel='Preconditioned original (seconds)',title='Complete setup + solve + recovery')
    ax.legend(frameon=False,fontsize=8)
    cohorts=[([by_id[c['id']]['gain']['lower'] for c in cases if c['cohort']==cohort],color)
             for cohort,color in [('prior','#1976a3'),('heldout','#ca6b24')]]
    cohorts.append(([c['gain']['lower'] for c in warm['cases']],'#39824c'))
    for i,(ratios,color) in enumerate(cohorts):
        bx.scatter([i]*len(ratios),ratios,c=color,alpha=.7,s=28)
        bx.scatter([i],[gm(ratios)],c='black',s=110,marker='_',zorder=3)
    bx.axhline(1,color='#555555',linestyle='--')
    bx.axhline(2,color='#999999',linestyle=':')
    bx.set(xticks=[0,1,2],xticklabels=['Prior fresh','Held-out fresh','Held-out warm'],yscale='log',
           ylabel='Charged gain over faster GB baseline',title='Fixed primary candidate; no per-case selection')
    bx.text(.02,.02,'Black marks: geometric means',transform=bx.transAxes,fontsize=8)
    fig.suptitle('Trace-aware affine preconditioning · matched charged timings · 3 repetitions',fontsize=12)
    fig.savefig(OUT/'comparison.svg',metadata={'Date':None})
    fig.savefig(OUT/'comparison.png',dpi=150)
    plt.close(fig)
    rows=[]
    for cohort in ('prior','heldout','all'):
        for candidate in ('original_preconditioned','split_preconditioned'):
            a=data['aggregate'][cohort][candidate]
            rows.append(f"| {cohort} | {candidate} | {a['complete_pairs']} | {a['geometric_mean_gain']:.3f} | {a['wins']}/{a['losses']} | {a['preprocessing_solved']} | {a['geometric_mean_gain_vs_lrf']:.4f} |")
    (OUT/'aggregate-table.md').write_text(
        '| Cohort | Candidate | Complete cases | Geometric mean gain vs faster GB | Wins/losses | GB skipped | Gain vs existing LRF |\n'
        '| --- | --- | ---: | ---: | ---: | ---: | ---: |\n'+'\n'.join(rows)+'\n')
    bound_rows=[]
    for cohort in ('prior','heldout','all'):
        a=bounds['aggregate'][cohort]
        bound_rows.append(f"| {cohort} | {a['cases']} | {a['primary_complete']} | {a['geometric_mean_gain_lower_bound']:.3f} | {a['cases_with_censored_baseline']} | {a['exactly_resolved_best_baselines']} |")
    (OUT/'bounds-table.md').write_text(
        '| Cohort | All cases | Primary complete | Conservative geometric mean gain | Cases with censored baseline | Best medians resolved |\n'
        '| --- | ---: | ---: | ---: | ---: | ---: |\n'+'\n'.join(bound_rows)+'\n')
    primary=data['aggregate']['heldout']['original_preconditioned']
    prior=data['aggregate']['prior']['original_preconditioned']
    write=dict(primary_candidate=method,
               target='at least 2x charged geometric-mean gain vs faster original/split GB on each cohort',
               met=primary['complete_pairs']==12 and prior['complete_pairs']==33 and
                   primary['geometric_mean_gain']>=2 and prior['geometric_mean_gain']>=2,
               prior_gain=prior['geometric_mean_gain'],heldout_gain=primary['geometric_mean_gain'],
               heldout_complete_pairs=primary['complete_pairs'],
               secondary_interval_target_supported=bounds['target_supported'],
               warm_confirmation_target_supported=warm['target_supported'],
               warm_gain_lower_bound=warm['geometric_mean_gain_lower_bound'],
               explanation='Strict complete-pair gate retains all baseline censoring; secondary interval analysis is separate.',
               caveat='observed finite-grid practical significance gate, not an asymptotic or statistical theorem')
    (OUT/'target.json').write_text(json.dumps(write,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
