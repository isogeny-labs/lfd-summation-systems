#!/usr/bin/env python3
"""k=3 ordinary-ring static brute force (9 variables, D up to 15: 1,307,504 monomials). One unit instance from the scan."""
import sys, os, json, time, resource
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *
from static_bruteforce import system, sd_sparse_int
D = os.path.join(os.path.dirname(__file__), '..', 'data')
idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
rec = [json.loads(l) for l in open(D + '/scan_k3_p65537.jsonl')]
rec = [r for r in rec if r['unit'] and not r['locus5']][idx]
k = 3; r = 9
chains, G, delta = system(rec['p'], rec['a'], rec['b'], rec['t'], k)
log = open(D + f'/bruteforce_k3_{idx}.log', 'w')
t0 = time.time()
out = dict(rec)
try:
    sd, prof = sd_sparse_int(chains, G, rec['p'], r, 3 * k + delta, log=log)
    out.update(sd_bruteforce=sd, sd_profile=prof, status='ok')
except MemoryError:
    out.update(status='MemoryError')
out['seconds'] = round(time.time() - t0, 1)
out['maxrss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
out['sd_matches_prediction'] = out.get('sd_bruteforce') == out['predicted_sd']
json.dump(out, open(D + f'/bruteforce_k3_{idx}.json', 'w'), indent=1)
print({k_: out[k_] for k_ in ('status', 'predicted_sd', 'sd_bruteforce', 'seconds', 'maxrss_MB')})
