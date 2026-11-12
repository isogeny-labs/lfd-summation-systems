#!/usr/bin/env python3
"""Second, fully independent implementation path: exact_lfd() driven by the pure-Python literal-rule closure
(closure_naive: all monomial multipliers a with deg a + deg f <= D, own integer echelon) instead of the numpy
variable-wise closure.  Run on a stratified sample; compare with the repository."""
import json, os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lfd_engine as L
import run_checks as RC
S, R, C = RC.load()
rng = random.Random(99)
def pick(pred, n):
    pool = sorted(h for h, s in S.items() if pred(s))
    rng.shuffle(pool)
    strata = {}
    for h in pool:
        strata.setdefault((C[h][0]['kind'], R[h]['status'], R[h]['exact_lfd']), []).append(h)
    out, keys = [], sorted(strata)
    while len(out) < n and any(strata.values()):
        for k in keys:
            if strata[k] and len(out) < n:
                out.append(strata[k].pop())
    return out
hs = pick(lambda s: s['kind'] == 'prime' and s['k'] == 2, 10) + pick(lambda s: s['kind'] == 'prime' and s['k'] == 3, 20) \
   + pick(lambda s: s['kind'] == 'binary' and s['variables'] == 4, 10) + pick(lambda s: s['kind'] == 'binary' and s['variables'] == 6, 10) \
   + pick(lambda s: s['kind'] == 'binary' and s['variables'] == 8, 10)
bad = 0
out = open('naive_lfd.jsonl', 'a')
for h in hs:
    s, rep = S[h], R[h]
    t = time.time()
    if s['kind'] == 'prime':
        k, p = s['k'], s['prime']; r = 2 * k
        inputs, G = L.prime_inputs([(int(m), int(c)) for m, c in s['generators'][0]], k, p)
        pts = L.prime_points(G, k, p); Dmax = 2 * k + 1
    else:
        r = s['variables']; p = 2
        inputs = L.binary_inputs(s['equations'], r); pts = L.binary_points(s['equations'], r); Dmax = r // 2 + 2
    res = L.exact_lfd(inputs, r, p, pts, Dmax, closure_fn=L.closure_naive)
    ok = res['lfd'] == rep['exact_lfd'] and len(pts) == rep['point_count'] and res['falls'] == rep['strict_falls']
    bad += not ok
    out.write(json.dumps(dict(system_hash=h, kind=s['kind'], vars=r, mine_lfd=res['lfd'], repo_lfd=rep['exact_lfd'], points=len(pts), falls=res['falls'], repo_falls=rep['strict_falls'], ok=ok)) + '\n'); out.flush()
    print(h[:12], s['kind'], 'vars', r, 'naive-LFD', res['lfd'], 'repo', rep['exact_lfd'], 'falls', res['falls'], rep['strict_falls'], 'OK' if ok else 'MISMATCH', '%.0fs' % (time.time() - t), flush=True)
print('done', len(hs), 'mismatches', bad)
