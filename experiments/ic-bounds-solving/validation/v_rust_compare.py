#!/usr/bin/env python3
"""Validator check 3b: run the producers' Rust binary (built from the frozen commit) as a black box on instances
whose sd_static I computed by brute force, and compare maxwt / delta / sd / top-coefficient flag."""
import json, subprocess, sys, os

BIN = sys.argv[1]
OUT = sys.argv[2]
FILES = sys.argv[3:]
HERE = os.path.dirname(os.path.abspath(__file__))

picks = []
for f in FILES:
    rows = [json.loads(l) for l in open(f)]
    # take the first row and, if present, the first exceptional row (c_top_zero) of each file
    picks.append(rows[0])
    ex = [r for r in rows if r.get('c_top_zero')]
    if ex:
        picks.append(ex[0])
res = []
for r in picks:
    cmd = [BIN, '--m', str(r['m']), '--k', str(r['k']), '--p', str(r['p']), '--a', str(r['a']), '--b', str(r['b']), '--t', str(r['t'])]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    j = json.loads(out)
    rec = dict(cmd=' '.join(cmd), rust=dict(unit=j['unit'], maxwt=j.get('maxwt'), delta=j.get('delta'), sd_static=j.get('sd_static'),
               top_coeff_zero=j.get('top_coeff_zero')),
               mine=dict(unit=r['unit'], maxwt=r['maxwt'], delta=r['delta'], sd_inverse=r['sd_predicted'], sd_bruteforce=r['sd_static'],
                         c_top_zero=r['c_top_zero']))
    rec['agree'] = (j['unit'] == r['unit'] and j.get('maxwt') == r['maxwt'] and j.get('delta') == r['delta']
                    and j.get('sd_static') == r['sd_static'] and j.get('top_coeff_zero') == r['c_top_zero'])
    res.append(rec)
    print(rec['cmd'].split('/')[-1], '| rust', rec['rust'], '| mine', rec['mine'], '| agree', rec['agree'])
json.dump(res, open(OUT, 'w'), indent=1)
print('agree', sum(r['agree'] for r in res), '/', len(res))
