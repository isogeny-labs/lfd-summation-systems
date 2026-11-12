"""Exact dim V_4 (and its filtration profile) at depths k = 6..13 with the split engine.

Tests the degree-4 identity d_4(k) = 4k^2 + 172k - 571 (data so far: k = 6..10) at new depths, and
records dim(V_4 cap M_<=w).  Each closure is run to a verified closure (no early stop).  Instances:
the depth-10 probe family (a = 1, b = 0, seed 10) and a random curve (a = 4021, b = 9917, seed 12,
--want-solutions), over p = 12289 for k <= 12 and p = 40961 = 5 * 2^13 + 1 for k = 13.
Writes experiments/ic-bounds-d5/d4-check.json.
"""
import json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_bounds_d5 as B

rows = []
for k in range(6, 14):
    p = 12289 if k <= 12 else 40961
    for a, b, seed, want in ((1, 0, 10, False), (4021, 9917, 12, True)):
        t = time.time()
        inst = B.make_instance(p, a % p, b % p, seed, k, want)
        g, Z = inst['generator'], inst['points']
        r = B.split_closure(g, k, p, 4, B.column_count(k, 4) - len(Z), early_stop=False)
        assert r['status'] == 'closed' and r['echelon_ok']
        fit = 4 * k * k + 172 * k - 571
        row = dict(k=k, prime=p, a=a % p, b=b % p, seed=seed, target=inst['target'], points=len(Z),
                   dim_V4=r['rank'], fit=fit, matches_fit=r['rank'] == fit, profile=r['profile'],
                   contains_one=r['contains_one'], columns=r['columns'], seconds=round(time.time() - t, 1))
        rows.append(row)
        print(json.dumps(row), flush=True)
out = ROOT / 'experiments/ic-bounds-d5/d4-check.json'
out.write_text(json.dumps(dict(rows=rows, all_match=all(r['matches_fit'] for r in rows)), indent=1) + '\n')
print(json.dumps(dict(all_match=all(r['matches_fit'] for r in rows))))
