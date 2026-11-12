"""Item 4 (a)+(b): validate the symmetry-split block engine against the dense engine at small depth.

For every instance and depth k in --ks, and every degree D in --degrees:
  * reference: the dense exact closure csd.ic_fall_certificates.chain_degree_closure (exclusion-
    certificate checked when 1 is not in V_D).  Where the dense closure exceeds its budget (k = 7,
    D = 6) the reference is the local-to-global theorem applied to the dense certificate at D - 1
    (V_{D-1} = I cap M_<=D-1 with ev(M_<D-1) onto Z implies V_D = I cap M_<=D), labelled as such;
  * (a) csd.ic_chain_block_closure.block_closure_contains_one, the engine the depth-10 runner uses,
    started from the dense V_{D-1} basis: compare "1 in V_D" and sym_rank + anti_rank with the
    reference dimension;
  * (b) csd.ic_bounds_exact10.SplitClosure (unit and proper ideals, deterministic verification),
    from the generator rows alone (dense random mixing) and from the dense V_{D-1} (sparse mixing):
    compare status, exact dimension, dim(V_D cap M_<D) and 1 in V_D with the reference;
  * negative control: the generator rows plus the dense V_{D-1}, verified as if they were V_D, must
    fail the closure check whenever dim V_{D-1} < dim V_D;
  * profile: csd.ic_bounds_exact10.split_profile must reproduce the exact LFD and the strict falls of
    csd.ic_curve_sweep.prime_local_profile.
Writes experiments/ic-bounds-exact10/validation.jsonl, profiles.jsonl and validation-summary.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from csd import ic_dreg as I
from csd import ic_curve_sweep as S
from csd import ic_fall_certificates as C
from csd import ic_chain_block_closure as B
from csd import ic_bounds_exact10 as X

INSTANCES = [
    dict(name='j1728-p7681-x936', p=7681, a=1, b=0, x=936, scale=1),
    dict(name='rand-p257-x58', p=257, a=2, b=3, x=58, scale=1),
    dict(name='rand-p7681-x2244', p=7681, a=5, b=7, x=2244, scale=1),
    dict(name='rand-p12289-x11294', p=12289, a=3, b=11, x=11294, scale=1),
    dict(name='j0-p7681-seed41', p=7681, a=0, b=7, target_seed=41, scale=1),
    dict(name='coset5-p12289-seed42', p=12289, a=3, b=11, target_seed=42, scale=5),
    dict(name='proper-p12289-seed43', p=12289, a=1, b=0, target_seed=43, scale=1, want_points=True),
    dict(name='proper-coset7-p7681-seed44', p=7681, a=5, b=7, target_seed=44, scale=7, want_points=True),
]


def maxrss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


def build(inst, k):
    p = inst['p']
    E = I.PrimeCurve(p, inst['a'], inst['b'])
    assert (4 * inst['a'] ** 3 + 27 * inst['b'] ** 2) % p and (p - 1) % (1 << k) == 0
    if 'x' in inst:
        x = inst['x']
        assert E.lift_x(x) is not None
    c = inst['scale']
    if c != 1:
        assert pow(c, 1 << k, p) != 1, 'the coset scale lies in H'
    g0 = S.subgroup_generator_value(p, k)
    if 'x' not in inst:
        rng = random.Random(inst['target_seed'])
        for _ in range(500):  # want_points: redraw (same seeded stream) until the system has scalar solutions
            x = S.random_rational_target(E, rng, p)[0]
            if not inst.get('want_points') or S.prime_points(E, k, g0, c, c, x):
                break
        else:
            raise RuntimeError('no target with solutions')
    g = S.prime_coset_generator(E, x, k, c, c)
    Z = S.prime_points(E, k, g0, c, c, x)
    if c == 1:
        assert g == I.pkm_generator(E, x, k)
    return g, Z, g0, x


def ghash(g):
    return hashlib.sha256(json.dumps(sorted(g.items())).encode()).hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ks', type=int, nargs='+', default=[5, 6, 7])
    ap.add_argument('--degrees', type=int, nargs='+', default=[4, 5, 6])
    ap.add_argument('--instances', nargs='*', default=None)
    ap.add_argument('--seed', type=int, default=20261107)
    ap.add_argument('--dense-max-cells', type=int, default=80_000_000)
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-exact10')
    ap.add_argument('--tag', default='')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    vpath = args.out / f'validation{args.tag}.jsonl'
    ppath = args.out / f'profiles{args.tag}.jsonl'
    vf, pf = vpath.open('w'), ppath.open('w')
    started = time.time()
    insts = [i for i in INSTANCES if args.instances is None or i['name'] in args.instances]
    for inst in insts:
        for k in args.ks:
            p = inst['p']
            g, Z, g0, x = build(inst, k)
            base = dict(instance=inst['name'], prime=p, a=inst['a'], b=inst['b'], target_x=x, scale=inst['scale'],
                        k=k, points=len(Z), unit=not Z, generator_sha=ghash(g))
            dense = {}
            for D in range(2, max(args.degrees) + 1):
                try:
                    V = C.chain_degree_closure([g], 2, k, p, D, max_columns=20000, max_cells=args.dense_max_cells)
                except ValueError as error:
                    dense[D] = str(error)
                    continue
                if not V.contains_one:
                    assert C.check_closed_exclusion_certificate(V, [g])
                dense[D] = V
            ref = S.prime_local_profile(g, k, p, Z, g0, max_degree=max(args.degrees), max_columns=20000,
                                        max_cells=args.dense_max_cells)
            for D in args.degrees:
                row = dict(base, degree=D, cols=X.monomial_count(k, D))
                V = dense.get(D)
                if isinstance(V, C.ChainDegreeClosure):
                    row.update(ref_source='dense exact closure', ref_dim=V.quotient_rank, ref_contains_one=V.contains_one,
                               ref_low_dim=int(sum(V.masks[j].bit_count() < D for j in V.pivots)))
                else:
                    W = dense.get(D - 1)
                    sc0 = X.SplitClosure(g, k, p, D - 1, points=Z, generator_value=g0)
                    cert_prev = (isinstance(W, C.ChainDegreeClosure) and
                                 (W.contains_one or (W.quotient_rank == sc0.target and sc0.onto_lt)))
                    if cert_prev:
                        sc1 = X.SplitClosure(g, k, p, D, points=Z, generator_value=g0)
                        row.update(ref_source='local-to-global from the dense certificate at D-1 (proof)',
                                   ref_dim=sc1.target, ref_contains_one=not Z,
                                   ref_low_dim=X.monomial_count(k, D - 1) - sc1.eval_rank_lt)
                    else:
                        row.update(ref_source=f'none ({V})', ref_dim=None, ref_contains_one=None, ref_low_dim=None)
                # (a) the original block engine, from the dense V_{D-1}
                W = dense.get(D - 1)
                if isinstance(W, C.ChainDegreeClosure):
                    t = time.time()
                    res = B.block_closure_contains_one(k, p, D, W.basis, list(W.masks), seed=args.seed + D,
                                                       combos_per_round=4000, max_rounds=60)
                    row['block_engine'] = dict(contains_one=res['contains_one'],
                                               rank=res['sym_rank'] + res['anti_rank'], sym_rank=res['sym_rank'],
                                               anti_rank=res['anti_rank'], rounds=res['rounds'],
                                               seconds=round(time.time() - t, 2))
                    del res
                # (b) the new split engine: from scratch (dense mixing) and from V_{D-1} (sparse mixing)
                t = time.time()
                r1 = X.split_degree(g, k, p, D, points=Z, generator_value=g0, seed=args.seed + 10 * D)
                row['split_scratch'] = {key: r1[key] for key in ('status', 'rank', 'exact_dimension', 'low_dimension',
                                                                'contains_one', 'closed', 'reached', 'target',
                                                                'evaluation_onto_lt', 'rounds', 'verify_passes')}
                row['split_scratch']['seconds'] = round(time.time() - t, 2)
                if isinstance(W, C.ChainDegreeClosure):
                    t = time.time()
                    r2 = X.split_degree(g, k, p, D, points=Z, generator_value=g0, prev=(W.basis, list(W.masks)),
                                        seed=args.seed + 10 * D + 1, terms=4)
                    row['split_from_prev_sparse'] = {key: r2[key] for key in ('status', 'rank', 'exact_dimension',
                                                                             'low_dimension', 'contains_one', 'closed',
                                                                             'rounds', 'verify_passes')}
                    row['split_from_prev_sparse']['seconds'] = round(time.time() - t, 2)
                    # negative control: V_{D-1} + generator rows presented as V_D must fail verification
                    sc = X.SplitClosure(g, k, p, D, points=Z, generator_value=g0)
                    sc.add_generator_rows()
                    sc.add_previous(W.basis, list(W.masks))
                    rank0 = sc.rank
                    missing = sc.verify_pass()
                    row['negative_control'] = dict(rank_before=rank0, missing=missing, one_before=sc.contains_one,
                                                   expected_missing=(row['ref_dim'] is not None and rank0 < row['ref_dim']),
                                                   detected=missing > 0)
                # comparisons
                chk = {}
                if row['ref_dim'] is not None:
                    chk['scratch_dim'] = r1['exact_dimension'] == row['ref_dim']
                    chk['scratch_low_dim'] = r1['low_dimension'] == row['ref_low_dim']
                    chk['scratch_one'] = r1['contains_one'] == row['ref_contains_one']
                    chk['scratch_exact'] = r1['status'] in ('reached', 'closed-deficient', 'reached-not-onto')
                    if 'split_from_prev_sparse' in row:
                        q = row['split_from_prev_sparse']
                        chk['prev_dim'] = q['exact_dimension'] == row['ref_dim']
                        chk['prev_one'] = q['contains_one'] == row['ref_contains_one']
                    if 'block_engine' in row:
                        be = row['block_engine']
                        chk['block_one'] = be['contains_one'] == row['ref_contains_one']
                        chk['block_rank_le_ref'] = be['rank'] <= row['ref_dim']
                        chk['block_rank_eq_ref_when_not_one'] = (be['contains_one'] or be['rank'] == row['ref_dim'])
                    if 'negative_control' in row:
                        nc = row['negative_control']
                        if not nc['one_before']:
                            chk['negative_control'] = nc['detected'] == nc['expected_missing']
                row['checks'] = chk
                row['all_pass'] = all(chk.values()) if chk else None
                row['maxrss_mb'] = maxrss_mb()
                vf.write(json.dumps(row) + '\n')
                vf.flush()
                print(json.dumps(dict(instance=inst['name'], k=k, D=D, points=len(Z), ref=row['ref_dim'],
                                      scratch=r1['status'], dim=r1['exact_dimension'], pass_=row['all_pass'],
                                      sec=round(time.time() - started), rss=row['maxrss_mb'])), flush=True)
            # profile comparison
            t = time.time()
            prof = X.split_profile(g, k, p, Z, g0, max_degree=max(args.degrees), seed=args.seed)
            prow = dict(base, reference=dict(status=ref['status'], exact_lfd=ref['exact_lfd'],
                                             strict_falls=ref['strict_falls'], quotient_ranks=ref.get('quotient_ranks')),
                        split=dict(status=prof['status'], exact_lfd=prof['exact_lfd'], strict_falls=prof['strict_falls'],
                                   dims=prof['dims'],
                                   statuses=[q['status'] for q in prof['per_degree']]),
                        seconds=round(time.time() - t, 2))
            prow['lfd_agree'] = prof['exact_lfd'] == ref['exact_lfd']
            prow['falls_agree'] = [f for f in prof['strict_falls'] if f >= 3] == [f for f in ref['strict_falls'] if f >= 3]
            pf.write(json.dumps(prow) + '\n')
            pf.flush()
            print(json.dumps(dict(profile=inst['name'], k=k, ref_lfd=ref['exact_lfd'], split_lfd=prof['exact_lfd'],
                                  falls=prof['strict_falls'], agree=prow['lfd_agree'] and prow['falls_agree'])), flush=True)
    vf.close()
    pf.close()


if __name__ == '__main__':
    main()
