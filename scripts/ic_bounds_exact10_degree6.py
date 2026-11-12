"""Degree-D split closure for one two-summand prime chain instance (default D = 6): unit or proper.

Memory-lean replacement for scripts/run_ic_chain_depth10_degree6.py.  That runner embeds the whole
V_5 basis as a full-width int64 matrix (18,569 x 60,460 x 8 B = 9.0 GB) and then forms int64
projections of it (4.5 GB each) and int64 V_5 (3.2 GB): more than the 15 GB machine, and four such
processes (anon-rss 12.7-13.3 GB) were OOM-killed at 22:34 on 2026-11-07.  Here the V_5 file is
memory-mapped (uint16), inserted 500 rows at a time, and every float64 tile is bounded (--tile).

Instance: y^2 = x^3 + a x + b over F_p, depth k (2^k | p-1), target abscissa --target-x (or drawn
as in scripts/run_ic_chain_depth10_probe.py from --seed, redrawn with --want-solutions), subgroup
or coset scale --scale (both summands, so the generator is swap-symmetric).
Input (optional): --v5 DIR holding basis.npy (rows of a verified V_{D-1} over the squarefree masks of
weight <= D-1, the probe's format) and state.json; --v5-sha256 is checked when given.  Without --v5
the closure starts from the generator rows alone.
Certificates: csd/ic_bounds_exact10.py ('reached': V_D = I cap M_<=D, LFD <= D; 'closed-deficient':
verified closed below target, LFD >= D + 1; otherwise 'undecided', no bound).  For proper ideals the
rows of the input basis are checked to vanish on the solution set.

Checkpoints the two parts after every round in --ckpt; rerunning the same command resumes.
--dry-run prints the instance and the memory model and exits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from csd import ic_dreg as I
from csd import ic_curve_sweep as S
from csd import ic_chain_shift_closure as SC
from csd import ic_bounds_exact10 as X


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--k', type=int, default=10)
    ap.add_argument('--degree', type=int, default=6)
    ap.add_argument('--prime', type=int, default=12289)
    ap.add_argument('--a', type=int, default=1)
    ap.add_argument('--b', type=int, default=0)
    ap.add_argument('--scale', type=int, default=1)
    ap.add_argument('--target-x', type=int, default=None)
    ap.add_argument('--seed', type=int, default=10, help='target seed (probe convention) when --target-x is absent')
    ap.add_argument('--want-solutions', action='store_true')
    ap.add_argument('--v5', type=Path, default=None)
    ap.add_argument('--v5-sha256', default=None)
    ap.add_argument('--ckpt', type=Path, default=None)
    ap.add_argument('--tag', default=None)
    ap.add_argument('--out-dir', type=Path, default=ROOT / 'experiments/ic-bounds-exact10')
    ap.add_argument('--combos', type=int, default=2000)
    ap.add_argument('--terms', type=int, default=4, help='sparse mixing terms per batch and shift (0 = dense)')
    ap.add_argument('--tile', type=int, default=1024)
    ap.add_argument('--round-seed', type=int, default=20261107)
    ap.add_argument('--max-rounds', type=int, default=400, help='stop (and checkpoint) after this many rounds in total')
    ap.add_argument('--no-verify', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    p, k, D = args.prime, args.k, args.degree
    started = time.time()

    def log(**kw):
        kw['elapsed'] = round(time.time() - started)
        print(json.dumps(kw), flush=True)

    assert (p - 1) % (1 << k) == 0 and (4 * args.a ** 3 + 27 * args.b ** 2) % p
    E = I.PrimeCurve(p, args.a, args.b)
    g0 = S.subgroup_generator_value(p, k)
    c = args.scale
    if args.target_x is not None:
        x = args.target_x
        assert E.lift_x(x) is not None
        Z = S.prime_points(E, k, g0, c, c, x)
        attempts = 0
    else:
        rng = random.Random(args.seed)
        for attempts in range(1, 201):
            x = S.random_rational_target(E, rng, p)[0]
            Z = S.prime_points(E, k, g0, c, c, x)
            if Z or not args.want_solutions:
                break
        else:
            raise SystemExit('no target with solutions')
    g = S.prime_coset_generator(E, x, k, c, c)
    tag = args.tag or f'k{k}-D{D}-p{p}-x{x}'
    gsha = hashlib.sha256(json.dumps(sorted(g.items())).encode()).hexdigest()
    mem = X.memory_model(k, D, prime=p, prev_rank=0, combos=args.combos, tile=args.tile)
    log(stage='instance', prime=p, a=args.a, b=args.b, k=k, degree=D, scale=c, target_x=x, attempts=attempts,
        points=len(Z), generator_sha256=gsha, memory_model_gb=round(mem['peak_bytes'] / 1e9, 2),
        cols=mem['cols'], sym_cols=mem['sym_cols'])
    if args.dry_run:
        return

    sc = X.SplitClosure(g, k, p, D, points=Z, generator_value=g0, tile=args.tile)
    log(stage='target', target=sc.target, eval_rank_le=sc.eval_rank_le, eval_rank_lt=sc.eval_rank_lt,
        evaluation_onto_lt=sc.onto_lt)
    ckpt = args.ckpt or Path(os.environ.get('TMPDIR', '/tmp')) / f'ic-bounds-exact10-{tag}'
    ckpt.mkdir(parents=True, exist_ok=True)
    v5_hash = None
    resumed = False
    prior_seconds = 0.0
    if (ckpt / 'state.json').exists():
        meta = json.loads((ckpt / 'state.json').read_text())
        assert meta['generator_sha256'] == gsha and meta['degree'] == D and meta['k'] == k
        for name in ('sym', 'anti'):
            part = getattr(sc, name)
            part.basis = np.load(ckpt / f'{name}_basis.npy')
            part.pivots = np.load(ckpt / f'{name}_pivots.npy')
            assert part.basis.shape == (len(part.pivots), part.cols)
        sc.rounds, sc.history, sc.passes = meta['round'], meta['history'], meta['passes']
        v5_hash = meta.get('v5_sha256')
        prior_seconds = meta.get('seconds_total', 0.0)
        resumed = True
        log(stage='resumed', round=sc.rounds, rank=sc.rank)
    else:
        sc.add_generator_rows()
        log(stage='generator-rows', rank=sc.rank)
        if args.v5 is not None:
            state = json.loads((args.v5 / 'state.json').read_text())
            assert state.get('phase') in ('verify', 'reached'), state.get('phase')
            v5_hash = sha256_file(args.v5 / 'basis.npy')
            if args.v5_sha256:
                assert v5_hash == args.v5_sha256, ('V_5 basis hash mismatch', v5_hash)
            v5 = np.load(args.v5 / 'basis.npy', mmap_mode='r')
            masks_prev = SC.monomials(k, D - 1)
            assert v5.shape[1] == len(masks_prev) and v5.shape[0] == state['rank'], (v5.shape, state.get('rank'))
            if Z:  # every input row must vanish on the solution set
                N = 1 << k
                ev = np.asarray([[pow(g0, ((m & (N - 1)) * (z & (N - 1)) + (m >> k) * (z >> k)) % N, p)
                                  for z in Z] for m in masks_prev], dtype=np.int64)
                for s0 in range(0, len(v5), 2000):
                    blk = np.asarray(v5[s0:s0 + 2000], dtype=np.int64)
                    assert not np.any((blk @ ev) % p), 'an input row does not vanish on Z'
            sc.add_previous(v5, masks_prev, chunk=500)
            del v5
            log(stage='v5-inserted', v5_rows=state['rank'], v5_sha256=v5_hash, rank=sc.rank)

    def checkpoint():
        for name in ('sym', 'anti'):
            part = getattr(sc, name)
            np.save(ckpt / f'{name}_basis.tmp.npy', part.basis)
            os.replace(ckpt / f'{name}_basis.tmp.npy', ckpt / f'{name}_basis.npy')
            np.save(ckpt / f'{name}_pivots.tmp.npy', part.pivots)
            os.replace(ckpt / f'{name}_pivots.tmp.npy', ckpt / f'{name}_pivots.npy')
        meta = dict(k=k, degree=D, prime=p, generator_sha256=gsha, round=sc.rounds, history=sc.history,
                    passes=sc.passes, rank=sc.rank, v5_sha256=v5_hash,
                    seconds_total=round(prior_seconds + time.time() - started, 1))
        (ckpt / 'state.tmp.json').write_text(json.dumps(meta))
        os.replace(ckpt / 'state.tmp.json', ckpt / 'state.json')

    if not resumed:
        checkpoint()

    def on_event(d):
        log(**d)
        checkpoint()

    res = sc.run(seed=args.round_seed + sc.rounds, combos=args.combos, max_rounds=args.max_rounds,
                 verify=not args.no_verify, terms=args.terms or None, log=on_event)
    checkpoint()
    record = dict(instance=dict(prime=p, a=args.a, b=args.b, k=k, scale=c, target_x=x, seed=args.seed,
                                points=len(Z), generator_sha256=gsha),
                  degree=D, v5_sha256=v5_hash, combos=args.combos, terms=args.terms, round_seed=args.round_seed,
                  result={key: v for key, v in res.items() if key != 'rank_history'}, rank_history=res['rank_history'],
                  seconds=round(time.time() - started), seconds_total=round(prior_seconds + time.time() - started))
    out = args.out_dir / f'degree{D}-{tag}.json'
    out.write_text(json.dumps(record, indent=1) + '\n')
    log(stage='done', status=res['status'], rank=res['rank'], target=res['target'], out=str(out))


if __name__ == '__main__':
    main()
