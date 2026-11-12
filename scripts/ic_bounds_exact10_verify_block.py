#!/usr/bin/env python3
"""Re-execute the existing block engine csd/ic_chain_block_closure.block_closure_contains_one on a
sample of rows where 1 is NOT in V_D, fed with a V_{D-1} computed by the independent closure in
scripts/ic_bounds_exact10_verify.py, and compare with the independent dim V_D."""
import hashlib, json, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
import ic_bounds_exact10_verify as V
from csd import ic_chain_block_closure as BC
SAMPLE = [('j1728-p7681-x936', 6, 4, 936), ('rand-p257-x58', 6, 4, 58), ('proper-p12289-seed43', 6, 5, 6061),
          ('rand-p12289-x11294', 7, 4, 11294), ('rand-p257-x58', 7, 4, 58)]
out = []
inst = {g['name']: g for g in V.GRID}
for name, k, D, x in SAMPLE:
    g = inst[name]; p = g['p']
    g0 = V.my_subgroup_generator(p, k)
    gen = V.my_generator(p, g['a'], g['b'], x, k, g['scale'])
    prev = None
    for d in range(2, D):
        sp, cl, st = V.closure_degree(gen, k, d, p, prev=prev)
        prev = (sp, cl)
    sp6, cl6, st6 = V.closure_degree(gen, k, D, p, prev=prev)
    sp, cl = prev
    t = time.time()
    res = BC.block_closure_contains_one(k, p, D, cl.R[:cl.r].astype(np.int64), [int(m) for m in sp.masks], seed=7)
    rec = dict(instance=name, k=k, D=D, x=x, my_dimV=cl6.r, my_one=cl6.contains_one, my_status=st6,
               block_one=bool(res['contains_one']), block_rank=int(res.get('sym_rank', 0)) + int(res.get('anti_rank', 0)),
               block_keys=sorted(res), sec=round(time.time() - t, 1))
    rec['agree'] = rec['block_one'] == rec['my_one'] and (rec['my_one'] or rec['block_rank'] == rec['my_dimV'])
    out.append(rec); print(json.dumps(rec), flush=True)
h = hashlib.sha256((ROOT / 'csd/ic_chain_block_closure.py').read_bytes()).hexdigest()
(ROOT / 'experiments/ic-bounds-exact10/verify/block-rerun.json').write_text(json.dumps(dict(block_engine_sha256=h, rows=out), indent=1) + '\n')
