"""d_ff (first strict fall) for k = 7..10 from the quotient closure at D = 2..4 only (no d_lf search).
By the Bridge, = ordinary value.  Usage: python3 dff_only.py k n   [pcc:unverified]"""
import sys, os, json, random, time
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import lf_quotient as Q
k, n = int(sys.argv[1]), int(sys.argv[2])
p = prime_above(3 * 4 ** k, k) if k <= 9 else prime_above(3 * 4 ** k, k)
rng = random.Random(4242 + k)
rows = []
while len(rows) < n:
    a, b = random_curve(p, rng); t = rng.randrange(p)
    if zeros_count(a, b, t, p, k): continue
    t0 = time.time(); per = []; prev = 0; dff = None
    for D in range(2, 5):
        r = Q.closure_q(Q.QSpace(k, D), a, b, t, p)
        per.append(dict(D=D, dim=r['dim'], dim_low=r['dim_low'], one=r['one']))
        if dff is None and r['dim_low'] > prev: dff = D
        prev = r['dim']
    rows.append(dict(p=p, a=a, b=b, t=t, k=k, d_ff=dff, per=per, sec=round(time.time() - t0, 1)))
    print(rows[-1], flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'dff_k%d.json' % k), 'w'), indent=1)
