#!/usr/bin/env python3
"""How much would the degree lift of experiments/ic-bounds-k14/LIFT.md save on the split presentation?

For a closed V_D without 1 (this engine, split), count with the definitions of LIFT.md section 1:
  regular weight-D monomials     the orbit is a pivot column of both parts (self-paired: of the + part)
  regular weight-(D+1) monomials every weight-D divisor is regular
  special weight-(D+1) monomials the others
and dim T = special top columns + free columns of V_D, the space the lift has to eliminate in,
against the number of columns of the direct closure at D + 1.
"""
import json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
from csd import ic_bounds_split_m as S, ic_bounds_split_gpu as G   # noqa: E402
import ic_split_gpu_gate as GATE                                    # noqa: E402


def main():
    cells = [tuple(map(int, c.split(','))) for c in sys.argv[1:]] or [(4, 2, 5), (5, 2, 5), (4, 3, 6), (6, 2, 6)]
    for m, k, D in cells:
        p, gens, _ = GATE.instance(m, k)
        rec, eng = G.closure(gens, m, k, p, D, device='cpu', split=True, verify=False)
        assert rec['status'] == 'closed', rec['status']
        sp, co = eng.sp, eng.co
        n = len(sp)
        pivot = []
        for q, pt in enumerate(eng.parts):
            f = np.zeros(n, dtype=bool)
            f[co.parts[q]['rep'][pt.piv]] = True
            pivot.append(f | f[co.sigma])
        selfp = co.sigma == np.arange(n)
        regular = pivot[0] & (pivot[1] | selfp)
        top = sp.deg == D
        sp1 = S.SplitSpace(m, k, D + 1, m - 2)
        to1 = np.array([sp1.index[mon] for mon in sp.mons])
        n1 = len(sp1)
        special = np.zeros(n1, dtype=bool)
        src = np.flatnonzero(top & ~regular)                 # free-type weight-D monomials
        for v in range(len(sp1.cols)):
            img = sp1.cols[v, to1[src]]
            img = img[sp1.deg[img] == D + 1]
            special[img] = True
        ntop1 = int(np.sum(sp1.deg == D + 1))
        nspecial = int(special.sum())
        T = sum(pt.nfree for pt in eng.parts)        # free columns of V_D, both parts
        rep1 = np.minimum(np.arange(n1), G.sigma_permutation(sp1))
        self1 = G.sigma_permutation(sp1) == np.arange(n1)
        sp_plus = int(np.sum(special & (rep1 == np.arange(n1))))
        sp_minus = int(np.sum(special & (rep1 == np.arange(n1)) & ~self1))
        T += sp_plus + sp_minus
        line = dict(m=m, k=k, D=D, p=p, columns_D=n, dim_V_D=rec['dim_V_D'], weight_D_monomials=int(top.sum()),
                    regular_weight_D=int((top & regular).sum()), columns_D_plus_1=n1, weight_D1_monomials=ntop1,
                    special_weight_D1=nspecial, regular_weight_D1=ntop1 - nspecial,
                    regular_fraction_D1=round((ntop1 - nspecial) / ntop1, 5),
                    dim_T_both_parts=int(T), T_over_columns_D_plus_1=round(T / n1, 4))
        print(json.dumps(line), flush=True)
        with open(ROOT / 'experiments/gpu-split-m/raw/lift_regular_count.jsonl', 'a') as f:
            f.write(json.dumps(line) + '\n')


if __name__ == '__main__':
    main()
