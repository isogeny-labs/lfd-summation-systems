"""d_lf (lfd_engine, ordinary ring) and d_ff on the exceptional instances (sd = 2k+1) at k = 3, 4.  [pcc:unverified]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import static_bruteforce as SB
R = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'all_instances.json')))
out = []
for k, n in ((3, 25), (4, 20)):
    ex = [r for r in R if r['k'] == k and r['sd_inverse'] != 2 * k + 2][:n]
    for r in ex:
        o = SB.lfd_per_degree(r['a'], r['b'], r['t'], r['p'], k, Dmax=k + 3)
        out.append(dict(k=k, p=r['p'], a=r['a'], b=r['b'], t=r['t'], sd=r['sd_inverse'], sd_bruteforce=r['sd_bruteforce'], d_lf=o['d_lf'], d_ff=o['d_ff']))
        print(out[-1], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'exceptions_dlf.json'), 'w'), indent=1)
