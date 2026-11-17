#!/usr/bin/env python3
"""Profile of the reference dense engine (csd.ic_bounds_split_m.InPlaceClosure) on certified cells.

The engine is not modified: a subclass times its three kernels and reads densities off its state.
One JSON line per closure: rank at stop, density of the stored basis (all columns / free columns
only), density of raw products and of reduced products, products reduced per pivot, time shares.
"""
import argparse, json, os, pathlib, sys, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd import ic_dreg as I                      # noqa: E402
from csd import ic_bounds_split_m as S            # noqa: E402


class Probe(S.InPlaceClosure):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.t = dict(products=0.0, reduce=0.0, absorb=0.0)
        self.c = dict(products=0, zero_rows=0, pivots=0, nnz_raw=0, nnz_red=0, red_rows=0)
        self.trace = []

    def saturate(self, log=None, target_rank=None):
        rounds = 0
        while self.queue and not self.contains_one:
            items = self.queue[:self.batch]
            del self.queue[:self.batch]
            t0 = time.perf_counter()
            P = self._product_rows(items)
            t1 = time.perf_counter()
            self.c['nnz_raw'] += int(np.count_nonzero(P))
            R = self._reduce(P)
            t2 = time.perf_counter()
            nzrow = np.any(R != 0, axis=1)
            self.c['zero_rows'] += int(len(items) - nzrow.sum())
            self.c['nnz_red'] += int(np.count_nonzero(R))
            self.c['red_rows'] += int(nzrow.sum())
            before = len(self.piv)
            t3 = time.perf_counter()
            self._absorb(R)
            t4 = time.perf_counter()
            self.c['products'] += len(items)
            self.c['pivots'] += len(self.piv) - before
            self.t['products'] += t1 - t0
            self.t['reduce'] += t2 - t1
            self.t['absorb'] += t4 - t3
            rounds += 1
            if rounds % 4 == 0 or not self.queue:
                self.trace.append(self.snap())
        return rounds

    def snap(self):
        r = len(self.piv)
        if not r:
            return dict(rank=0)
        nz = int(np.count_nonzero(self.B))
        free = np.ones(self.n, bool)
        free[self.piv] = False
        nzf = nz - r
        return dict(rank=r, queue=len(self.queue), density_all=round(nz / (r * self.n), 5),
                    density_free=round(nzf / max(1, r * int(free.sum())), 5),
                    low_rank=int(np.sum(self.sp.deg[self.piv] < self.sp.D)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cell', required=True, help='m,k')
    ap.add_argument('--degrees', required=True, help='e.g. 5,6')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    m, k = map(int, a.cell.split(','))
    rec0 = json.loads((ROOT / f'experiments/ic-bounds-split-m/cell-m{m}-k{k}.jsonl').read_text().splitlines()[0]) \
        if (ROOT / f'experiments/ic-bounds-split-m/cell-m{m}-k{k}.jsonl').exists() else None
    if rec0 is None:
        raise SystemExit('no certified record for this cell; use profile_engine.py')
    p = rec0['p']
    E = I.PrimeCurve(p, rec0['a'], rec0['b'])
    gens, _ = S.split_generators_m(E, rec0['target'][0], k, m)
    for D in map(int, a.degrees.split(',')):
        sp = S.SplitSpace(m, k, D, m - 2)
        cl = Probe(sp, p, 768)
        t0 = time.time()
        cl.add_inputs([sp.row(g, p) for g in gens if S.generator_degree(g) <= D])
        cl.saturate()
        sec = time.time() - t0
        c = cl.c
        out = dict(m=m, k=k, p=p, D=D, n=len(sp), rank_at_stop=int(len(cl.piv)), contains_one=bool(cl.contains_one),
                   seconds=round(sec, 2), loadavg=os.getloadavg()[0], t=cl.t, **c,
                   products_per_pivot=round(c['products'] / max(1, c['pivots']), 3),
                   zero_fraction=round(c['zero_rows'] / max(1, c['products']), 4),
                   raw_product_density=round(c['nnz_raw'] / max(1, c['products'] * len(sp)), 6),
                   reduced_product_density=round(c['nnz_red'] / max(1, c['red_rows'] * len(sp)), 5),
                   final=cl.snap(), trace=cl.trace)
        with open(a.out, 'a') as f:
            f.write(json.dumps(out) + '\n')
        print(json.dumps({k_: v for k_, v in out.items() if k_ != 'trace'}), flush=True)


if __name__ == '__main__':
    main()
