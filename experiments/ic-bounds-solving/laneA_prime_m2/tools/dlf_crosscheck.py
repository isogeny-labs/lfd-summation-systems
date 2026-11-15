"""Ordinary-ring lfd_engine vs quotient engine: d_lf, d_ff and dimV per degree.  [pcc:unverified]
Usage: python3 dlf_crosscheck.py k p n"""
import sys, os, json, random, time
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import lf_quotient as Q, static_bruteforce as SB
k, p, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(77 + k)
inst, _ = unit_instances(p, k, n, rng)
rows = []
for (a, b, t) in inst:
    t0 = time.time()
    q = Q.profile(a, b, t, p, k, Dmax=k + 3); tq = time.time() - t0
    t0 = time.time()
    o = SB.lfd_per_degree(a, b, t, p, k, Dmax=q['d_lf'] or (k + 3)); to = time.time() - t0
    # ordinary dims include chain span; compare lfd + fall
    agree = (o['d_lf'] == q['d_lf']) and (o['d_ff'] == q['d_ff'])
    rows.append(dict(p=p, a=a, b=b, t=t, k=k, d_lf_ord=o['d_lf'], d_ff_ord=o['d_ff'], d_lf_quot=q['d_lf'], d_ff_quot=q['d_ff'],
                     falls_ord=o['falls'], agree=agree, sec_ord=round(to, 1), sec_quot=round(tq, 1)))
    print(rows[-1], flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'dlf_crosscheck_k%d_p%d.json' % (k, p)), 'w'), indent=1)
