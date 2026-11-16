#!/usr/bin/env python3
"""Brute-force sd (ordinary ring, k=2) on exceptional instances (maxwt(g^-1) < 3k) found by the scan."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *
from sd_inverse import analyse
from static_bruteforce import system, sd_sparse_int
D = os.path.join(os.path.dirname(__file__), '..', 'data')
S = json.load(open(D + '/summary.json'))
seen = set(); out = []
for e in S['exceptions']:
    if e['k'] != 2: continue
    key = (e['p'], e['a'], e['b'], e['t'])
    if key in seen: continue
    seen.add(key)
    if len(out) >= 12: break
    rec = analyse(e['p'], e['a'], e['b'], e['t'], 2)
    chains, G, delta = system(e['p'], e['a'], e['b'], e['t'], 2)
    sd, prof = sd_sparse_int(chains, G, e['p'], 6, 2 * 3 + delta)
    rec['sd_bruteforce'] = sd; rec['sd_matches_prediction'] = (sd == rec['predicted_sd']); rec['k'] = 2
    out.append({x: rec[x] for x in ('k', 'p', 'a', 'b', 't', 'delta', 'locus5', 'maxwt_ginv', 'predicted_sd', 'sd_bruteforce', 'sd_matches_prediction')})
    print(out[-1], flush=True)
json.dump(out, open(D + '/bruteforce_k2_exceptions.json', 'w'), indent=1)
