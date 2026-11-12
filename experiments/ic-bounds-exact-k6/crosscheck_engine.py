import sys, json
sys.path.insert(0, '.')
from csd import ic_bounds_fall as F
p, a, b, t = 1153, 139, 1094, 185
co = {(2, 2): 1, (2, 1): -2 * t, (1, 2): -2 * t, (2, 0): t * t, (0, 2): t * t, (1, 1): -2 * (t * t + a),
      (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b, (0, 0): a * a - 4 * b * t}
gm = {e1 | (e2 << 6): c % p for (e1, e2), c in co.items()}
out = {}
for D in (4, 5):
    grid = F.make_grid(6, D)
    r = F.closure(gm, grid, p, verify=True, reduced=True)
    w = [int(x) for x in grid.weights[r['pivots']]]
    out[D] = dict(rank=len(r['pivots']), closed=bool(r['closed']), contains_one=bool(r['contains_one']),
                  profile=[sum(1 for x in w if x <= e) for e in range(D + 1)])
    print(D, out[D], flush=True)
json.dump(out, open('experiments/ic-bounds-exact-k6/p1153-a139-b1094-t185/crosscheck-engine.json', 'w'), indent=1)
