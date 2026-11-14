"""Test M of docs/proofs/N9-block-four.md: is 1 in V_D(k1, k2) for unequal depths k1 <= k2?

For each pair: closure (fast engine of ic_bounds_block4_lib, block RREF of csd/ic_bounds_d5_verify.Part),
deterministic full sweep (every shift of every final low row reduces to zero), whether 1 is in V_D, dimension,
profile, and the number of solutions on mu_{2^k1} x mu_{2^k2} (0 = unit ideal).

usage: python3 ic_bounds_block4_mixed.py OUT.jsonl D k1min k1max k2max [maxcols] [unit|unit2] [--slow-check]
       pairs k1min <= k1 <= k1max, k1 <= k2 <= k2max with |M_<=D| <= maxcols (default 30000).
       --slow-check also recomputes each pair with the slow engine of ic_bounds_growth_lib (small pairs only)
       and compares dimension and profile.
"""
import json
import sys
import time
from math import comb

import ic_bounds_block4_lib as B4


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    slow = '--slow-check' in sys.argv
    out, D, k1min, k1max, k2max = args[0], int(args[1]), int(args[2]), int(args[3]), int(args[4])
    maxcols = int(args[5]) if len(args) > 5 else 30000
    inst = B4.UNIT2 if (len(args) > 6 and args[6] == 'unit2') else B4.UNIT
    p = inst['p']
    g = B4.generator(inst)
    with open(out, 'a') as fh:
        for k1 in range(k1min, k1max + 1):
            for k2 in range(k1, k2max + 1):
                ncols = sum(comb(k1 + k2, j) for j in range(D + 1))
                if ncols > maxcols:
                    continue
                t0 = time.time()
                C = B4.FastClosure(k1, k2, D, p, [g])
                C.run()
                bad = C.sweep_full()
                rec = dict(test='mixed', inst=inst, D=D, k1=k1, k2=k2, columns=C.L.c, dim=C.dim(),
                           profile=C.profile(), one=C.contains_one(), sweep_failures=bad,
                           solutions=B4.unit_check(inst, k1, k2), seconds=round(time.time() - t0, 2),
                           basis_sha=C.P.sha()[:16])
                if slow:
                    L = B4.Layout2(k1, k2, D)
                    Bm, piv = B4.closure(L, g, p)
                    rec['slow_engine'] = dict(dim=len(piv), profile=B4.profile(L, piv),
                                              one=B4.has_one(L, piv), sweep_failures=B4.verify_closed(L, Bm, piv, p))
                    rec['engines_agree'] = (rec['slow_engine']['dim'] == rec['dim'] and
                                            rec['slow_engine']['profile'] == rec['profile'])
                print(json.dumps(rec), flush=True)
                fh.write(json.dumps(rec) + '\n'); fh.flush()


if __name__ == '__main__':
    main()
