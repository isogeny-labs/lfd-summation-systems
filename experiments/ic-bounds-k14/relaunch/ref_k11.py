#!/usr/bin/env python3
"""Reference engine (csd/ic_bounds_d5_verify.py, pushed as refeng_d5_verify.py) on the unit instance over
p = 65537: closure with a full sweep, same record fields as scripts/ic_bounds_d5_verify.py `closure`.

usage: python3 ref_k11.py [K [D]]     (default 11 6)
The target abscissa 4271 is fixed (the instance of experiments/ic-bounds-k14/instances), so the repository's
target-drawing helpers are not needed on the pod.  WORKSPACE overrides /workspace (tests).
"""
import hashlib
import json
import os
import resource
import sys
import time
from math import comb

import refeng_d5_verify as V

P, A, B, X, TARGET = 65537, 1, 0, 4271, [4271, 11720]
K = int(sys.argv[1]) if len(sys.argv) > 1 else 11
D = int(sys.argv[2]) if len(sys.argv) > 2 else 6
WS = os.environ.get('WORKSPACE', '/workspace')
T0 = time.time()


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['rss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(json.dumps(kw), flush=True)


def main():
    Z = V.solutions(P, A, B, X, K)
    g = V.chain_generator(P, A, B, X, K)
    log(stage='instance', target=TARGET, points=len(Z), k=K, D=D)
    t = time.time()
    C = V.Closure(g, K, P, D, split=True, gens='all', log=log)
    C.run(deadline=None)
    prof, parts = C.profile()
    rec = dict(k=K, D=D, p=P, split=True, gens='all', columns=C.L.c, finished=True, rank=C.rank(), profile=prof,
               part_profiles=parts, part_ranks=[Q.rank for Q in C.parts], contains_one=C.contains_one(), rounds=C.rounds,
               closure_seconds=round(time.time() - t, 1), madds=sum(Q.madds for Q in C.parts),
               peak_entries=[Q.peak_entries for Q in C.parts], rref_ok=all(Q.check_rref() for Q in C.parts),
               basis_sha256=[Q.sha() for Q in C.parts])
    rec['sweep_random'] = C.sweep_random(t=6)
    rec['sweep_full'] = C.sweep_full(deadline=None)
    cols = sum(comb(2 * K, w) for w in range(D + 1))
    rec.update(points=len(Z), target=cols - len(Z))
    rec['deficient'] = rec['rank'] < rec['target']
    rec['seconds'] = round(time.time() - t, 1)
    rec['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    rec.update(a=A, b=B, seed=10, want=False, instance_target=TARGET, driver='experiments/ic-bounds-k14/relaunch/ref_k11.py',
               source_sha256={n: hashlib.sha256(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), n), 'rb').read()).hexdigest()
                              for n in ('refeng_d5_verify.py', 'ref_k11.py')})
    out = os.path.join(WS, 'out', 'reference-engine-unit-p65537-k%d-D%d.json' % (K, D))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + '.tmp', 'w') as f:
        json.dump(rec, f, indent=1)
        f.write('\n')
    os.replace(out + '.tmp', out)
    log(stage='saved', file=out, rank=rec['rank'], contains_one=rec['contains_one'], basis_sha256=rec['basis_sha256'])


if __name__ == '__main__':
    main()
