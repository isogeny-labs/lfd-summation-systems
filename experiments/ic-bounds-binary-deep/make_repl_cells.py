#!/usr/bin/env python3
"""Write cells_cpu.json: the systems of cells.json whose closure at the certificate degree of the first run has at
most --max-columns columns, with max_degree set to that certificate degree (the ladder stops there anyway).
The subset is chosen from the first run's outcomes (outcomes-2026-11-12.json); that is a choice of what is cheap
enough to repeat on the CPU, not a choice by result: every system is repeated on the GPU in any case."""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--max-columns', type=int, default=60000)
ap.add_argument('--time-limit', type=int, default=3600)
a = ap.parse_args()
cells = json.load(open(os.path.join(HERE, 'cells.json')))
rows = {r['name']: r for r in json.load(open(os.path.join(HERE, 'outcomes-2026-11-12.json')))['rows']}
out = dict(cells)
out['order'] = []
kept = dropped = 0
for cell in cells['order']:
    systems = []
    for s in cell['systems']:
        r = rows[s['name']]
        cols = r['columns'][r['certificate_degree']] if isinstance(r['columns'], list) else r['columns']
        if cols <= a.max_columns:
            systems.append(dict(s, max_degree=r['certificate_degree'], est_seconds=0, time_limit_s=a.time_limit))
            kept += 1
        else:
            dropped += 1
    if systems:
        out['order'].append(dict(cell, systems=systems))
json.dump(out, open(os.path.join(HERE, 'cells_cpu.json'), 'w'), indent=1)
print('cpu subset: %d systems kept, %d left to the GPU only' % (kept, dropped))
