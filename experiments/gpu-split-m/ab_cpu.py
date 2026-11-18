#!/usr/bin/env python3
"""Interleaved CPU A/B: csd.ic_bounds_split_m.closure (InPlaceClosure) against
csd.ic_bounds_split_gpu.closure with the numpy backend (split and plain), same instances, same
degrees, verification pass included in both.  Arm order rotates every round; the 1-minute load
average is recorded before every run.  One JSON line per run."""
import argparse, json, os, pathlib, sys, time
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
from csd import ic_bounds_split_m as S, ic_bounds_split_gpu as G   # noqa: E402
import ic_split_gpu_gate as GATE                                    # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', default='4,2:5;4,2:6;5,2:5;4,3:5;4,3:6')
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--arms', default='reference,split,plain')
    ap.add_argument('--batch', type=int, default=1024)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    arms = a.arms.split(',')
    for item in a.cells.split(';'):
        mk, D = item.split(':')
        m, k = map(int, mk.split(','))
        D = int(D)
        p, gens, _ = GATE.instance(m, k)
        for r in range(a.rounds):
            for i in range(len(arms)):
                arm = arms[(i + r) % len(arms)]
                load = os.getloadavg()[0]
                t0 = time.perf_counter()
                if arm == 'reference':
                    rec = S.closure(gens, m, k, p, D)
                    extra = {}
                else:
                    rec, _ = G.closure(gens, m, k, p, D, device='cpu', split=(arm == 'split'), batch=a.batch)
                    extra = dict(madds=rec['ops']['total'], stats=rec['stats'])
                sec = time.perf_counter() - t0
                line = dict(m=m, k=k, D=D, p=p, arm=arm, round=r, seconds=round(sec, 3), load_before=load,
                            dim_V_D=rec['dim_V_D'], contains_one=rec['contains_one'],
                            verified_closed=rec['verified_closed'], low_profile=rec['low_profile'], **extra)
                with open(a.out, 'a') as f:
                    f.write(json.dumps(line) + '\n')
                print(json.dumps({k_: v for k_, v in line.items() if k_ not in ('stats', 'low_profile')}), flush=True)


if __name__ == '__main__':
    main()
