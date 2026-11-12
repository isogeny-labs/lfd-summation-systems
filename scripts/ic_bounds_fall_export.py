"""Export the degree-4 closures V_4 (and their low parts) of a two-summand prime chain system for a
range of depths k, and V_5 ∩ M_{<=4} for small k, on one generic random instance (random curve,
random target of large order).  Bounds program, item 5 (fall spaces).

Group-algebra engine: csd/ic_bounds_fall.py (once-per-pivot shift closure, verified closed, RREF
with weight-sorted columns, birth round of every pivot).  One instance is used for every depth:
the generator g = NF(S_3(x1, x2, t)) has exponents <= 2 and is depth-independent as exponent
pairs; only its chain-mask encoding depends on k.

Outputs, under experiments/ic-bounds-fall/<tag>/ (tag = p<prime>-s<seed>):
  instance.json          prime, curve, target, target order, solution set, generator, hashes
  k<k>-D4.npz            RREF basis of V_4 (uint16/uint32), pivots, pivot weights, birth rounds
  k<k>-D4.json           dims, low profile (dim V_4 ∩ M_{<=d}), rank history, closure flag,
                         comparison with 4k^2+172k-571 and 4k+87, basis sha256
  k<k>-D5-low.npz/.json  V_5 ∩ M_{<=4} (basis rows of weight <= 4) and the V_5 summary

Example:  ulimit -v 3000000; python3 scripts/ic_bounds_fall_export.py --prime 40961 --seed 7 --k 6 7 8 9 10 --k5 6 7 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_bounds_fall as F  # noqa: E402
from csd import ic_curve_sweep as S  # noqa: E402


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_basis(path: Path, result: dict, grid: F.Grid, p: int, rows=None) -> dict:
    basis = result['basis'] if rows is None else rows
    dtype = np.uint16 if p < 65536 else np.uint32
    np.savez(path, basis=basis.astype(dtype), pivots=np.asarray(result['pivots'], dtype=np.int64),
             pivot_weights=np.asarray(result['pivot_weights'], dtype=np.int64),
             birth=np.asarray(result['birth'], dtype=np.int64), masks=np.asarray(grid.masks, dtype=np.int64),
             prime=np.asarray([p]), k=np.asarray([grid.k]), D=np.asarray([grid.D]))
    return dict(file=path.name, sha256=sha256_file(path), basis_sha256=F.row_digest(basis), shape=list(basis.shape))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prime', type=int, default=40961)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--k', type=int, nargs='*', default=[6, 7, 8, 9, 10], help='depths for the degree-4 export')
    ap.add_argument('--k5', type=int, nargs='*', default=[], help='depths for the degree-5 low-part export')
    ap.add_argument('--no-basis', action='store_true', help='write only the JSON summaries (no .npz)')
    ap.add_argument('--product-cap', type=int, default=20000)
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-fall')
    args = ap.parse_args()
    p = args.prime
    tag = f'p{p}-s{args.seed}'
    out = args.out / tag
    out.mkdir(parents=True, exist_ok=True)
    inst_path = out / 'instance.json'
    if inst_path.exists():
        inst = json.loads(inst_path.read_text())
    else:
        inst = F.random_instance(p, 6, args.seed)
        inst['instance_hash'] = F.instance_hash(inst)
        inst['note'] = ('random curve y^2=x^3+ax+b and random rational target of large order drawn from '
                        'random.Random(seed) at depth 6; the solution set listed is the depth-6 one, '
                        'per-depth solution counts are in the k*.json files')
        inst_path.write_text(json.dumps(inst, indent=1))
    pairs = inst['generator_pairs']
    E = F.I.PrimeCurve(p, inst['a'], inst['b'])
    print(json.dumps(dict(tag=tag, a=inst['a'], b=inst['b'], target=inst['target'], target_order=inst['target_order'])))

    for k in args.k:
        t0 = time.time()
        g = F.encode_generator(pairs, k)
        grid = F.make_grid(k, 4)
        g0 = S.subgroup_generator_value(p, k)  # the depth-k subgroup generator (depth-specific)
        Z = S.prime_points(E, k, g0, 1, 1, inst['target'][0])
        res = F.closure(g, grid, p, product_cap=args.product_cap)
        profile = F.low_profile(res, grid)
        count4 = len(grid.masks)
        summary = dict(prime=p, k=k, D=4, seed=args.seed, solutions=len(Z), monomials_le4=count4,
                       ideal_dimension_le4_if_onto=count4 - len(Z), dim_V4=res['rank'], low_profile=profile,
                       dim_V4_le3=profile[3], dim_V4_le2=profile[2], formula_4k2=4 * k * k + 172 * k - 571,
                       formula_4k87=4 * k + 87, matches_formulas=(res['rank'] == 4 * k * k + 172 * k - 571 and profile[3] == 4 * k + 87),
                       contains_one=res['contains_one'], closed_verified=res['closed'], rounds=res['rounds'],
                       rank_history=res['history'], generator_rows=len(res['generator_mus']),
                       birth_histogram={str(r): int(sum(1 for b in res['birth'] if b == r)) for r in range(res['rounds'] + 1)},
                       low_birth_histogram={str(r): int(sum(1 for b, w in zip(res['birth'], res['pivot_weights']) if b == r and w <= 3))
                                            for r in range(res['rounds'] + 1)},
                       seconds=round(time.time() - t0, 1))
        if not args.no_basis:
            summary['basis_file'] = save_basis(out / f'k{k}-D4.npz', res, grid, p)
        else:
            summary['basis_sha256'] = F.row_digest(res['basis'])
        (out / f'k{k}-D4.json').write_text(json.dumps(summary, indent=1))
        print(json.dumps({a: summary[a] for a in ('k', 'dim_V4', 'low_profile', 'formula_4k2', 'matches_formulas', 'closed_verified', 'rounds', 'solutions', 'seconds')}))

    for k in args.k5:
        t0 = time.time()
        g = F.encode_generator(pairs, k)
        grid = F.make_grid(k, 5)
        g0 = S.subgroup_generator_value(p, k)
        Z = S.prime_points(E, k, g0, 1, 1, inst['target'][0])
        res = F.closure(g, grid, p, product_cap=min(args.product_cap, 6000))
        profile = F.low_profile(res, grid)
        count = lambda d: sum(math.comb(2 * k, j) for j in range(d + 1))
        onto = True
        if Z:
            onto = S.prime_evaluation_rank(Z, [m for m in grid.masks if m.bit_count() < 5], k, g0, p, 10 ** 8) == len(Z)
        target = count(5) - len(Z)
        low = F.low_part(res, grid, 4)
        summary = dict(prime=p, k=k, D=5, seed=args.seed, solutions=len(Z), monomials_le5=count(5), monomials_le4=count(4),
                       dim_V5=res['rank'], low_profile=profile, dim_V5_le4=profile[4], ideal_dimension_le5=target,
                       evaluation_onto_below_5=onto, reached_ideal=(res['rank'] == target and onto),
                       certificate=('V_5 = I ∩ M_<=5 (rank = |M_<=5| - |Z|, M_<5 onto Z): LFD <= 5, V_5 ∩ M_<=4 = I ∩ M_<=4'
                                    if res['rank'] == target and onto else 'closure below the ideal dimension or evaluation not onto'),
                       contains_one=res['contains_one'], closed_verified=res['closed'], rounds=res['rounds'],
                       rank_history=res['history'], seconds=round(time.time() - t0, 1))
        if not args.no_basis:
            sub = dict(res)
            keep = [i for i, w in enumerate(res['pivot_weights']) if w <= 4]
            sub['pivots'] = [res['pivots'][i] for i in keep]
            sub['pivot_weights'] = [res['pivot_weights'][i] for i in keep]
            sub['birth'] = [res['birth'][i] for i in keep]
            summary['low_basis_file'] = save_basis(out / f'k{k}-D5-low.npz', sub, grid, p, rows=low)
        (out / f'k{k}-D5.json').write_text(json.dumps(summary, indent=1))
        print(json.dumps({a: summary[a] for a in ('k', 'dim_V5', 'low_profile', 'ideal_dimension_le5', 'reached_ideal', 'contains_one', 'rounds', 'seconds')}))


if __name__ == '__main__':
    main()
