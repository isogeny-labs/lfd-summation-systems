#!/usr/bin/env python3
"""Two-summand reference on the SAME (p, curve, target) as the three-summand depth runs.

System {2k chains, NF(S_3(x1, x2, t))}; exact LFD from csd.ic_curve_sweep.prime_local_profile
(unit or local-to-global certificate, the engine used for the m = 2 reference in
experiments/ic-bounds-m3).  Output: experiments/ic-bounds-m3-depth/m2.jsonl (resumable).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from csd import ic_curve_sweep as S  # noqa: E402
from csd import ic_dreg as I  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ic_bounds_m3_depth import INSTANCES, OUTDIR, make_instance  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--instances', nargs='+', default=list(INSTANCES))
    ap.add_argument('--k', nargs='+', type=int, default=[2, 3, 4, 5, 6, 7])
    ap.add_argument('--max-columns', type=int, default=8000)
    ap.add_argument('--max-cells', type=int, default=40_000_000)
    a = ap.parse_args()
    path = os.path.join(OUTDIR, 'm2.jsonl')
    done = set()
    if os.path.exists(path):
        done = {(r['instance'], r['k']) for r in map(json.loads, open(path))}
    for name in a.instances:
        inst = make_instance(*INSTANCES[name])
        p = inst['p']
        E = I.PrimeCurve(p, inst['a'], inst['b'])
        for k in a.k:
            if (p - 1) % (1 << k) or (name, k) in done:
                continue
            t0 = time.time()
            g0 = S.subgroup_generator_value(p, k)
            gen = S.prime_coset_generator(E, inst['t'], k, 1, 1)
            pts = S.prime_points(E, k, g0, 1, 1, inst['t'])
            res = S.prime_local_profile(gen, k, p, pts, g0, max_degree=2 * k + 1,
                                        max_columns=a.max_columns, max_cells=a.max_cells)
            row = dict(instance=name, p=p, k=k, status=res['status'], exact_lfd=res['exact_lfd'],
                       lfd_lower_bound=res['lfd_lower_bound'], lfd_upper_bound=res['lfd_upper_bound'],
                       solutions=len(pts), strict_falls=res['strict_falls'],
                       quotient_ranks=res['quotient_ranks'], seconds=round(time.time() - t0, 1),
                       engine='ic_curve_sweep.prime_local_profile')
            print(json.dumps({k_: row[k_] for k_ in ('instance', 'k', 'status', 'exact_lfd', 'lfd_lower_bound', 'solutions', 'seconds')}), flush=True)
            with open(path, 'a') as f:
                f.write(json.dumps(row) + '\n')


if __name__ == '__main__':
    main()
