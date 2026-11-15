#!/usr/bin/env python3
"""Scan: for each (k, p) draw random instances until n_units unit instances (or max_tries / time budget).
usage: scan.py k p n_units max_tries budget_s [locus]"""
import sys, os, json, time, random
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *
from sd_inverse import analyse
k, p, nu, mt, budget = map(int, sys.argv[1:6])
locus = len(sys.argv) > 6 and sys.argv[6] == 'locus'
rng = random.Random(hash((k, p, locus)) & 0xffffffff ^ 20261115)
rng = random.Random(f'{k}-{p}-{locus}-20261115')
tag = f'scan_k{k}_p{p}' + ('_locus' if locus else '')
path = os.path.join(os.path.dirname(__file__), '..', 'data', tag + '.jsonl')
t0 = time.time(); tries = 0; units = 0
with open(path, 'w') as fh:
    while units < nu and tries < mt and time.time() - t0 < budget:
        a, b, t = (locus_instance if locus else random_instance)(rng, p)
        rec = analyse(p, a, b, t, k, check_conv=(k <= 3 and units < 3))
        tries += 1
        if rec['unit']:
            units += 1
        fh.write(json.dumps(rec) + '\n'); fh.flush()
summ = dict(k=k, p=p, locus=locus, tries=tries, units=units, wall_s=round(time.time() - t0, 1))
print(json.dumps(summ))
with open(os.path.join(os.path.dirname(__file__), '..', 'data', 'scan_summaries.jsonl'), 'a') as fh:
    fh.write(json.dumps(summ) + '\n')
