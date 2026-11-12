#!/usr/bin/env python3
"""Literal-rule (full monomial multipliers, pure-Python exact arithmetic) vs variable-wise numpy closure on larger instances."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lfd_engine as L
ROOT = '/home/user/ml-cryptanalysis'
S = {}
for l in open(ROOT + '/experiments/ic-curve-sweep/systems.jsonl'):
    s = json.loads(l); S[s['system_hash']] = s
jobs = [(h, int(D)) for h, D in (a.split(':') for a in sys.argv[1:])]
for pref, Dmax in jobs:
    h = [x for x in S if x.startswith(pref)][0]
    s = S[h]
    if s['kind'] == 'prime':
        inputs, G = L.prime_inputs([(int(m), int(c)) for m, c in s['generators'][0]], s['k'], s['prime'])
        r, p = 2 * s['k'], s['prime']
    else:
        r, p = s['variables'], 2
        inputs = L.binary_inputs(s['equations'], r)
    for D in range(Dmax + 1):
        t = time.time()
        a = L.closure(inputs, r, D, p)
        b = L.closure_naive(inputs, r, D, p)
        ok = (a['dim'], a['dim_low'], a['contains_one']) == (b['dim'], b['dim_low'], b['contains_one'])
        print(h[:12], 'r=%d p=%d D=%d' % (r, p, D), 'variable-wise', (a['dim'], a['dim_low']), 'literal', (b['dim'], b['dim_low']), 'EQUAL' if ok else 'DIFFER', '%.0fs' % (time.time() - t), flush=True)
