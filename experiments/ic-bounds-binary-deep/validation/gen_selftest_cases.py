#!/usr/bin/env python3
"""Generate the SELFTEST_CASES literal of csd/ic_bounds_bin_big.py from the EXISTING engine (reference values: LFD,
certificate degree, dim V_D, falls from csd.ic_local_certificate.local_profile; RREF hashes from the existing basis
re-reduced by the new canonical hash).  Appends to the engine file after the marker (idempotent)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, HERE)
import validate as V
import ic_bounds_bin_big_run as R
from csd import ic_bounds_bin_big as EB
from csd.ic_local_certificate import local_profile
recs = {json.loads(l)['name']: json.loads(l) for l in open(os.path.join(HERE, 'records.jsonl'))}
sel = {c['name']: c for c in json.load(open(os.path.join(HERE, 'selection.json')))}
want = []
for lab, k in (('affine', 5), ('affine', 6), ('unit', 5), ('const', 5), ('order-two', 6), ('affine', 7)):
    for n, r in recs.items():
        if r['cls'] == lab and r['k'] == k and n in sel and n not in [w for w, _ in want]:
            want.append((n, 'sel')); break
cases = []
for name, kind in want + [('sf12_6_a0', 'sub')]:
    inst = R.build_instance(V.spec_of(sel[name], name) if kind == 'sel' else [s for s in V.SUBFIELD if s['name'] == name][0])
    r = inst['variables']
    old = local_profile(inst['equations'], r, inst['points'], max_degree=7, max_columns=200000, max_bytes=10**9)
    D = old['exact_lfd'] if old['status'] == 'exact-unit' else old['local_upper_bound']
    hashes = [V.old_hash(inst['equations'], r, d, 10**9)[0] for d in range(D + 1)]
    cases.append(dict(name=name, system_hash=inst['system_hash'], variables=r, equations=inst['equations'], points=inst['points'],
                      max_degree=D + 1, lfd=old['exact_lfd'], cert=D, dims=old['quotient_ranks'], falls=old['strict_falls'], hashes=hashes))
text = '\nSELFTEST_CASES = ' + json.dumps(cases, separators=(',', ':')) + '\n'
p = os.path.join(ROOT, 'csd', 'ic_bounds_bin_big.py')
s = open(p).read()
mark = '# ---- generated selftest cases (scripts/../validation/gen_selftest_cases.py) ----\n'
head = s[:s.index(mark) + len(mark)]
open(p, 'w').write(head + text)
print(len(cases), 'cases', len(text), 'bytes', [c['name'] for c in cases])
