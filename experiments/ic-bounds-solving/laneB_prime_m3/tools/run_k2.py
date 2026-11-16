#!/usr/bin/env python3
"""k=2: brute-force sd (ordinary ring) vs delta+maxwt(supp g^-1), plus exact d_ff, d_lf (lfd_engine)."""
import sys, os, json, time, random
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *
from sd_inverse import analyse
from static_bruteforce import system, sd_sparse, lfd_ff

k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
out = os.path.join(os.path.dirname(__file__), '..', 'data', f'bruteforce_k{k}.jsonl')
rng = random.Random(20261115 + k)
per_prime = {17: 7, 4013: 7, 65537: 7}
locus_per_prime = {4013: 3, 65537: 3}
done = []
def unit_search(p, locus):
    for _ in range(100000):
        a, b, t = (locus_instance if locus else random_instance)(rng, p)
        rec = analyse(p, a, b, t, k)
        if rec['unit']:
            return rec
    return None
fh = open(out, 'w')
for locus, table in ((False, per_prime), (True, locus_per_prime)):
    for p, n in table.items():
        for i in range(n):
            rec = unit_search(p, locus)
            if rec is None:
                print('no unit instance', p, locus); continue
            r = 3 * k
            chains, G, delta = system(p, rec['a'], rec['b'], rec['t'], k)
            assert delta == rec['delta']
            t0 = time.time()
            sd, prof = sd_sparse(chains, G, p, r, 3 * k + delta)
            rec['sd_bruteforce'] = sd
            rec['sd_bruteforce_seconds'] = round(time.time() - t0, 2)
            rec['sd_profile'] = prof
            t0 = time.time()
            L = lfd_ff(chains, G, p, r, 9)
            rec['d_lf'] = L['lfd']; rec['d_ff'] = L['d_ff']; rec['falls'] = L['falls']; rec['lfd_certified_at'] = L['certified_at']
            rec['lfd_seconds'] = round(time.time() - t0, 2)
            rec['lfd_per_degree'] = L['per_degree']
            rec['sd_matches_prediction'] = (sd == rec['predicted_sd'])
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print(p, 'locus' if locus else 'generic', 'delta', delta, 'maxwt', rec['maxwt_ginv'], 'pred', rec['predicted_sd'], 'sd', sd,
                  'd_ff', rec['d_ff'], 'd_lf', rec['d_lf'], rec['sd_bruteforce_seconds'], rec['lfd_seconds'], flush=True)
