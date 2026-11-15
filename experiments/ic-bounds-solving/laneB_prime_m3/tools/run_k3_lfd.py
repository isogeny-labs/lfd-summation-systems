#!/usr/bin/env python3
"""k=3 exact d_lf / d_ff with lfd_engine (9 variables, Dmax=7: 11440 columns) on one unit instance."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *
from static_bruteforce import system, lfd_ff
D = os.path.join(os.path.dirname(__file__), '..', 'data')
rec = [json.loads(l) for l in open(D + '/scan_k3_p65537.jsonl')]
rec = [r for r in rec if r['unit'] and not r['locus5']][0]
chains, G, delta = system(rec['p'], rec['a'], rec['b'], rec['t'], 3)
t0 = time.time()
L = lfd_ff(chains, G, rec['p'], 9, 7)
rec.update(L); rec['lfd_seconds'] = round(time.time() - t0, 1)
json.dump(rec, open(D + '/lfd_k3.json', 'w'), indent=1)
print({k: rec[k] for k in ('lfd', 'd_ff', 'falls', 'certified_at', 'lfd_seconds')})
