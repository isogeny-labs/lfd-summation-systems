#!/usr/bin/env python3
"""Item 6 runner: exact LFD of binary ordinary / Koblitz systems (m = 2 large k, m = 3 small k).

Usage (resource rules of the workflow: ulimit -v 3000000, one BLAS thread):
  python3 scripts/ic_bounds_binary_run.py --family m2 --n 11 12 13 14 --k 8 9 10 \
      --out experiments/ic-bounds-binary/results-m2.jsonl
  python3 scripts/ic_bounds_binary_run.py --family m3-direct --n 5 6 7 8 9 --k 2 3 4 5 \
      --out experiments/ic-bounds-binary/results-m3-direct.jsonl
  python3 scripts/ic_bounds_binary_run.py --family m3-split --n 5 6 7 --k 2 3 \
      --out experiments/ic-bounds-binary/results-m3-split.jsonl

Each row carries the full case specification (field modulus, curve, basis, target,
target class, planted pair), the system hash, |Z|, the certificate and the
Macaulay sizes. Rows already in --out (same case_id) are skipped, so a run can be
continued with the same command.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_binary as B  # noqa: E402
from csd.ic_research_campaign import field  # noqa: E402


def specs(family, ns, ks, n_random, n_invariant, n_curves, seed, curve_indices=None):
    out = []
    for n in ns:
        F = field(n)
        curves = B.curves_for(n, seed)
        curves = [curves[i] for i in curve_indices] if curve_indices else curves[:n_curves]
        for k in ks:
            if not k < n:
                continue
            bases = B.bases_for(F, k, n_random, n_invariant, seed)
            for a, b, curve_class in curves:
                E = I.BinaryCurve(F, a, b)
                for base_kind, basis in bases:
                    values = I.span(basis)
                    for target, target_class, planted in B.targets_for(E, values, (seed, family, n, k, a, b, basis)):
                        spec = dict(family=family, n=n, modulus=F.modulus, k=k, a=a, b=b, curve_class=curve_class,
                                    base_kind=base_kind, basis=basis,
                                    frobenius_invariant=I.is_frobenius_invariant(F, basis),
                                    target=list(target), target_class=target_class,
                                    planted=[list(P) for P in planted] if planted else None, seed=seed)
                        spec['case_id'] = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:20]
                        out.append(spec)
    return out


def run(args):
    spec, max_degree, max_columns, max_bytes = args
    F = field(spec['n'])
    assert F.modulus == spec['modulus']
    E = I.BinaryCurve(F, spec['a'], spec['b'])
    target = tuple(spec['target'])
    planted = [tuple(P) for P in spec['planted']] if spec['planted'] else None
    t0 = time.time()
    if spec['family'] == 'm2':
        system = B.m2_system(E, spec['basis'], spec['base_kind'], target, spec['target_class'], planted)
    elif spec['family'] == 'm3-direct':
        system = B.m3_direct_system(E, spec['basis'], target)
    elif spec['family'] == 'm3-split':
        system = B.m3_split_system(E, spec['basis'], target)
    else:
        raise ValueError(spec['family'])
    prep = time.time() - t0
    cert = B.certify(system, max_degree, max_columns, max_bytes)
    row = dict(spec)
    row.update(variables=system['variables'], system_hash=system['system_hash'], solutions=len(system['points']),
               points_sha256=hashlib.sha256(json.dumps(system['points']).encode()).hexdigest()[:16],
               prep_seconds=round(prep, 3), **{k: v for k, v in system.items()
                                               if k.startswith('rational') or k in ('scaled_subfield', 'subfield')})
    row.update(cert)
    row['bounds'] = B.bounds(spec['family'], spec['k'], spec['n'], system['variables'])
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True, choices=['m2', 'm3-direct', 'm3-split'])
    ap.add_argument('--n', type=int, nargs='+', required=True)
    ap.add_argument('--k', type=int, nargs='+', required=True)
    ap.add_argument('--random-bases', type=int, default=2)
    ap.add_argument('--invariant-bases', type=int, default=2)
    ap.add_argument('--curves', type=int, default=4)
    ap.add_argument('--curve-indices', type=int, nargs='+', default=None,
                    help='subset of the four curves (0: (0,1), 1: (1,1), 2, 3: random ordinary)')
    ap.add_argument('--seed', type=int, default=B.SEED)
    ap.add_argument('--max-degree', type=int, default=7)
    ap.add_argument('--max-columns', type=int, default=70000)
    ap.add_argument('--max-bytes', type=int, default=900_000_000)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--limit', type=int, default=0, help='stop after this many new rows (0 = all)')
    ap.add_argument('--list', action='store_true', help='print the case count and exit')
    ap.add_argument('--submit-seconds', type=float, default=0,
                    help='stop starting new cases after this many seconds (0 = no limit)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    todo = specs(args.family, args.n, args.k, args.random_bases, args.invariant_bases, args.curves, args.seed,
                 args.curve_indices)
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)['case_id'] for l in open(args.out) if l.strip()}
    todo = [s for s in todo if s['case_id'] not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f'{len(todo)} cases to run ({len(done)} already in {args.out})', flush=True)
    if args.list or not todo:
        return
    work = [(s, args.max_degree, args.max_columns, args.max_bytes) for s in todo]
    t0 = time.time()

    def results(pool):
        # submit at most `workers` tasks at a time; stop submitting after --submit-seconds so that a
        # wall-clock-limited run ends with every started case finished (no half-done work is lost)
        pending, queue = [], list(work)
        while queue or pending:
            while queue and len(pending) < args.workers and (not args.submit_seconds
                                                             or time.time() - t0 < args.submit_seconds):
                pending.append(pool.apply_async(run, (queue.pop(0),)))
            if not pending:
                print(f'submit window closed; {len(queue)} cases left for the next run', flush=True)
                return
            for job in list(pending):
                if job.ready():
                    pending.remove(job)
                    yield job.get()
            time.sleep(0.2)

    with open(args.out, 'a') as fh, mp.get_context('fork').Pool(args.workers, maxtasksperchild=1) as pool:
        for i, row in enumerate(results(pool)):
            fh.write(json.dumps(row, sort_keys=True) + '\n')
            fh.flush()
            print(f"[{i + 1}/{len(work)} {time.time() - t0:.0f}s] {row['family']} n={row['n']} k={row['k']} "
                  f"{row['curve_class']} {row['base_kind']} {row['target_class']} |Z|={row['solutions']} "
                  f"lfd={row['exact_lfd']} {row['status']} D*={row.get('certificate_degree')} "
                  f"{row['closure_seconds']}s rss={row['peak_rss_mb']}MB", flush=True)


if __name__ == '__main__':
    main()
