#!/usr/bin/env python3
"""records/*.json + validation/*.json -> results.md (the tables of the README)."""
import json, glob, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(HERE, 'records', '*.json')))]
out = []
P = [r for r in recs if r['family'] == 'prime']
out.append('### Prime chains (all cases)\n')
out.append('| m | k | case | p | target | deg g | terms | \\|Z\\| | class | LFD | cert. degree | strict falls | kernel bound | status | s |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(P, key=lambda r: (r['m'], r['k'], r['case'])):
    out.append(f"| {r['m']} | {r['k']} | {r['case']} | {r['p']} | {r['target_class']} | {r['generator_degree']} | {r['generator_terms']} | {r['solutions']} | {r['ideal_class']} | {r['exact_lfd']} | {r['certificate_degree']} | {r['strict_falls']} | {r['kernel_upper_bound']} | {r['status']} | {r['seconds']} |")
B = [r for r in recs if r['family'] == 'binary']
out.append('\n### Binary (all cases)\n')
out.append('| m | n | k | r | case | curve | base | target | d_in | eqs | \\|Z\\| | class | LFD | cert. degree | strict falls | r+1 | status | s |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(B, key=lambda r: (r['m'], r['k'], r['n'], r['case'])):
    if 'curve' not in r:
        out.append(f"| {r['m']} | {r['n']} | {r['k']} | {r['m']*r['k']} | {r['case']} | | | | {r.get('input_degree')} | | {r.get('solutions')} | | | | | | {r['status']} | |")
        continue
    out.append(f"| {r['m']} | {r['n']} | {r['k']} | {r['variables']} | {r['case']} | {r['curve']} | {r['base_kind'][:12]} | {r['target_class'][:10]} | {r['input_degree']} | {r['equations']} | {r['solutions']} | {r['ideal_class']} | {r['exact_lfd']} | {r.get('certificate_degree')} | {r['strict_falls']} | {r['kernel_upper_bound']} | {r['status']} | {r['seconds']} |")
# aggregate
agg = collections.defaultdict(collections.Counter)
for r in P:
    agg[('prime', r['m'], r['k'])][(r['exact_lfd'])] += 1
for r in B:
    agg[('binary', r['m'], r['k'])][(r.get('exact_lfd') if r.get('exact_lfd') is not None else r['status'][:20])] += 1
out.append('\n### Histogram of exact LFD per (family, m, k)\n')
out.append('| family | m | k | LFD: count |')
out.append('|---|---|---|---|')
for key in sorted(agg):
    out.append(f"| {key[0]} | {key[1]} | {key[2]} | " + ', '.join(f'{v}: {c}' for v, c in sorted(agg[key].items(), key=lambda t: str(t[0]))) + ' |')
open(os.path.join(HERE, 'results.md'), 'w').write('\n'.join(out) + '\n')
print(len(recs))
