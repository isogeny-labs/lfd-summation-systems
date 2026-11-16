#!/usr/bin/env python3
"""Exact last fall degree of the split presentation with m summands (csd/ic_bounds_split_m.py).

  control   positive control: the recorded unit cases of experiments/ic-bounds-m3 (m = 3, split,
            free auxiliary) are recomputed with the new space and generators; generator equality,
            the exact value and dim V_D at every recorded degree must agree.
            negative control: planted targets (proper ideals) must never give 1.
  cell      seeded unit instances of one cell (m, k): random curve, random rational target redrawn
            until the grid has no zero of S_{m+1}; one record per instance.

Records go to experiments/ic-bounds-split-m/ (see PROTOCOL.md there).
"""
import argparse
import hashlib
import json
import pathlib
import random
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I                      # noqa: E402
from csd import ic_bounds_m3 as M3                # noqa: E402
from csd import ic_bounds_msum as MS              # noqa: E402
from csd import ic_bounds_split_m as S            # noqa: E402

OUT = ROOT / 'experiments/ic-bounds-split-m'


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def source_hashes():
    files = ['csd/ic_bounds_split_m.py', 'csd/ic_bounds_m3_batched.py', 'csd/ic_bounds_m3.py', 'csd/elim.py',
             'scripts/run_ic_split_m.py']
    return {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in files}


def control(args):
    cases = {}
    for line in (ROOT / 'experiments/ic-bounds-m3/cases.jsonl').read_text().splitlines():
        c = json.loads(line)
        cases[c['case_id']] = c
    recorded = {}
    for path in sorted((ROOT / 'experiments/ic-bounds-m3').glob('results*.jsonl')):
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if row['encoding'] == 'split' and row['results']['free']['status'] == 'exact-unit':
                recorded.setdefault(row['case_id'], row)
    rng = random.Random(20261116)
    chosen = []
    for k in args.depths:
        ids = sorted(i for i in recorded if cases[i]['k'] == k)
        rng.shuffle(ids)
        chosen += ids[:args.per_depth]
    positive, failures = [], []
    for cid in chosen:
        c, row = cases[cid], recorded[cid]
        E = I.PrimeCurve(c['prime'], c['a'], c['b'])
        t, k = c['target'][0], c['k']
        gens, _ = S.split_generators_m(E, t, k, 3)
        g1, g2, _ = M3.split_generators(E, t, k)
        same_gens = [{(key[0], (key[1],)): v for key, v in g.items()} for g in (g1, g2)] == gens
        want = row['results']['free']['exact_lfd']
        res = S.exact_unit_lfd(gens, 3, k, E.p, d_min=2, d_max=want + 1)
        dims = {r['D']: r['dim_V_D'] for r in res['profile']}
        old = {r['D']: r['quotient_rank'] for r in row['profile'] if r['D'] in dims}
        ok = same_gens and res['exact_lfd'] == want and dims == old
        positive.append(dict(case_id=cid, p=E.p, k=k, recorded_lfd=want, recomputed_lfd=res['exact_lfd'],
                             generators_equal=same_gens, dims_compared=len(old), dims_equal=dims == old, ok=ok))
        log(f'positive {cid} p={E.p} k={k} recorded={want} new={res["exact_lfd"]} dims={len(old)} ok={ok}')
        if not ok:
            failures.append(cid)
    negative = []
    for m, k, p, seed in args_negative():
        inst = MS.make_prime_instance(p, m, k, seed, planted=True)
        E = I.PrimeCurve(p, inst['a'], inst['b'])
        gens, _ = S.split_generators_m(E, inst['t'], k, m)
        res = S.exact_unit_lfd(gens, m, k, p, d_min=3, d_max=args.negative_cap)
        ok = res['exact_lfd'] is None and all(r['verified_closed'] for r in res['profile'])
        negative.append(dict(m=m, k=k, p=p, seed=seed, target=inst['t'], closures_without_one=[r['D'] for r in res['profile']],
                             status=res['status'], ok=ok))
        log(f'negative m={m} k={k} p={p} seed={seed} closures={[r["D"] for r in res["profile"]]} ok={ok}')
        if not ok:
            failures.append(f'negative-{m}-{k}-{p}-{seed}')
    out = dict(status='failed' if failures else 'passed', date=time.strftime('%Y-%m-%d'), positive=positive,
               negative=negative, failures=failures, source_sha256=source_hashes())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'control.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], positive=len(positive), negative=len(negative), failures=failures)))
    return bool(failures)


def args_negative():
    return [(3, 2, 7681, 1), (3, 2, 12289, 2), (3, 3, 7681, 1), (4, 2, 7681, 1), (4, 2, 12289, 2)]


def unit_instance(p, m, k, seed):
    """Random curve and random rational target with no zero of S_{m+1} on the grid H^m."""
    for attempt in range(200):
        inst = MS.make_prime_instance(p, m, k, seed * 1000 + attempt, planted=False)
        E = I.PrimeCurve(p, inst['a'], inst['b'])
        grid, _ = MS.prime_grid(E, inst['t'], m, k)
        if not MS.prime_zero_set(grid, m, k):
            inst['attempt'] = attempt
            return E, inst
    raise RuntimeError('no unit instance found')


def cell(args):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f'cell-m{args.m}-k{args.k}.jsonl'
    done = set()
    if path.exists():
        done = {(r['p'], r['seed']) for r in map(json.loads, path.read_text().splitlines())}
    hashes = source_hashes()
    for p in args.primes:
        assert (p - 1) % (1 << args.k) == 0
        for seed in range(args.seeds):
            if (p, seed) in done:
                continue
            E, inst = unit_instance(p, args.m, args.k, seed)
            gens, _ = S.split_generators_m(E, inst['t'], args.k, args.m)
            start = time.time()
            res = S.exact_unit_lfd(gens, args.m, args.k, p, d_min=3, d_max=args.cap, max_columns=args.max_columns, log=log)
            rec = dict(m=args.m, k=args.k, p=p, seed=seed, a=E.a, b=E.b, target=[inst['t'], inst['y']],
                       attempt=inst['attempt'], grid_zeros=0, generator_degrees=[S.generator_degree(g) for g in gens],
                       system_hash=S.system_hash(p, args.m, args.k, gens), status=res['status'],
                       exact_lfd=res['exact_lfd'], lfd_lower_bound=res['lfd_lower_bound'], profile=res['profile'],
                       seconds=round(time.time() - start, 2), source_sha256=hashes)
            with path.open('a') as f:
                f.write(json.dumps(rec) + '\n')
            log(f'cell m={args.m} k={args.k} p={p} seed={seed}: {res["status"]} {res["exact_lfd"]} ({rec["seconds"]}s)')
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('control')
    c.add_argument('--depths', type=int, nargs='+', default=[2, 3])
    c.add_argument('--per-depth', type=int, default=6)
    c.add_argument('--negative-cap', type=int, default=7)
    c.set_defaults(fn=control)
    c = sub.add_parser('cell')
    c.add_argument('--m', type=int, required=True)
    c.add_argument('--k', type=int, required=True)
    c.add_argument('--primes', type=int, nargs='+', default=[7681, 12289, 40961])
    c.add_argument('--seeds', type=int, default=2)
    c.add_argument('--cap', type=int, default=9)
    c.add_argument('--max-columns', type=int, default=None)
    c.set_defaults(fn=cell)
    args = ap.parse_args()
    sys.exit(args.fn(args))


if __name__ == '__main__':
    main()
