"""Exceptional locus of the two-summand prime chain closure (bounds program, item 2).

For a prime p, a depth k with 2^k | p-1 and a curve y^2 = x^3 + a x + b, every
rational target abscissa t (rhs(t) a square, including rhs(t) = 0) gives the
PKM system {g_t} + chain, g_t = NF(S_3(x1, x2, t)) (csd.ic_dreg.pkm_generator).
For each system this script computes, with the exact dense engine
(csd.ic_fall_certificates.chain_degree_closure, every closure that excludes 1
re-checked by check_closed_exclusion_certificate):

  dim V_D modulo the chains for D = 0..dmax,
  dim (V_D ∩ M_{<=D-1}) modulo the chains (the 'low' dimension),
  whether 1 is in V_D, |Z| (scalar chain solutions, csd.ic_curve_sweep.prime_points),
  the evaluation rank of M_{<D} on Z, the strict falls and, when the local-to-global
  certificate (or 1 in V_D) closes, the exact LFD (same rule as
  csd.ic_curve_sweep.prime_local_profile, which is re-run on a seeded subsample).

All dimensions are quotient (modulo chains) dimensions; the ordinary-ring
dimension is C(2k+D, D) - |M_{<=D}| + quotient (equal for D >= 2 by the
chain normal-form argument).

Usage
  python3 scripts/ic_bounds_locus.py run --prime 193 --k 6 --curve 0,1 --dmax 5 --out FILE.jsonl
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_curve_sweep as S  # noqa: E402
from csd import ic_fall_certificates as C  # noqa: E402
from csd import ic_chain_shift_closure as SC  # noqa: E402


def j_invariant(p, a, b):
    num = 1728 * 4 * a ** 3 % p
    den = (4 * a ** 3 + 27 * b * b) % p
    return num * pow(den, -1, p) % p


def rational_targets(E):
    """Every abscissa t with rhs(t) a square in F_p (rhs(t) = 0 included)."""
    out = []
    for t in range(E.p):
        if E.lift_x(t) is not None:
            out.append(t)
    return out


def generator_hash(g):
    return hashlib.sha256(json.dumps(sorted(g.items()), separators=(',', ':')).encode()).hexdigest()[:16]


DENSE_MAX_COLUMNS = 1500   # dense engine above this width needs ~0.6 GB; use the shift engine instead


def profile(E, t, k, dmax, extend_to=None, max_cells=60_000_000, check_certificates=True, product_cap=2400):
    """Closures V_0..V_dmax of {g_t}+chain, plus the exact LFD when certified.

    Degrees with at most DENSE_MAX_COLUMNS squarefree monomials use the dense
    engine and its exclusion certificate.  Wider degrees use
    csd.ic_chain_shift_closure.shift_closure (every row lies in V_D by
    construction): its result is exact when it either reaches the local target
    |M_<=D| - |Z| (then V_D = I ∩ M_<=D) or is verified closed under every
    variable shift of every row of weight < D; anything else is recorded as
    inexact and never used.  If extend_to > dmax and nothing is certified by
    dmax, the profile continues up to extend_to."""
    p = E.p
    g = I.pkm_generator(E, t, k)
    g0 = S.subgroup_generator_value(p, k)
    Z = S.prime_points(E, k, g0, 1, 1, t)
    r = 2 * k
    rows = []
    falls = []
    prev_dim = None
    status, exact, local_D = 'capped', None, None
    top = max(dmax, extend_to or dmax)
    for D in range(top + 1):
        if D > dmax and status != 'capped':
            break
        masks = SC.monomials(k, D)
        count = len(masks)
        low_count = sum(math.comb(r, d) for d in range(min(r, D - 1) + 1)) if D > 0 else 0
        ev = None
        if D > 0 and Z:
            ev = S.prime_evaluation_rank(Z, [m for m in masks if m.bit_count() < D], k, g0, p, 10 ** 8)
        if count <= DENSE_MAX_COLUMNS:
            V = C.chain_degree_closure([g], 2, k, p, D, max_columns=20000, max_cells=max_cells)
            one = V.contains_one
            dim = V.quotient_rank
            low = low_count if one else (sum(V.masks[j].bit_count() <= D - 1 for j in V.pivots) if D > 0 else 0)
            cert = None
            if not one and check_certificates:
                cert = C.check_closed_exclusion_certificate(V, [g])
                assert cert, ('exclusion certificate failed', p, E.a, E.b, t, k, D)
            engine, how = 'dense', ('one-derived' if one else 'closed-exclusion-certificate')
        else:
            target = count - len(Z)
            res = SC.shift_closure(g, k, p, D, target, product_cap=product_cap, max_cells=max_cells, verify=True)
            if res['reached']:
                how = 'reached-local-target'
                assert res['rank'] == target
            elif res['closed']:
                how = 'verified-closed'
            else:
                raise RuntimeError(('shift closure neither reached nor closed', p, E.a, E.b, t, k, D))
            dim, low = res['rank'], res['low_dimension']
            one = (not Z) and res['reached']
            if one:
                low = low_count
            cert, engine = None, 'shift'
        if prev_dim is not None and D >= 2 and low > prev_dim:
            falls.append(D)
        elif prev_dim is not None and D < 2:
            # degrees 0, 1: the generator has degree 2, so V_0 = V_1 = 0 (no falls)
            assert dim == 0 and low == 0
        rows.append(dict(D=D, monomials=count, dim=dim, low=low, one=one, codim=count - dim, eval_rank=ev,
                         engine=engine, exactness=how, exclusion_certificate=cert))
        prev_dim = dim
        if status == 'capped':
            if one:
                assert not Z
                status, exact = 'exact-unit', D
            elif D > 0 and Z and ev == len(Z) == count - dim:
                status, exact, local_D = 'exact-proper-local', max(falls, default=0), D
        if status != 'capped' and D >= dmax:
            break
    return dict(g=g, g_hash=generator_hash(g), Z=Z, rows=rows, falls=falls, status=status, exact_lfd=exact,
                local_upper_bound=local_D)


def run(args):
    p, k = args.prime, args.k
    a, b = [int(v) for v in args.curve.split(',')]
    E = I.PrimeCurve(p, a, b)
    assert (p - 1) % (1 << k) == 0
    targets = rational_targets(E)
    if args.targets:
        targets = [int(v) for v in args.targets.split(',')]
    rng = random.Random(args.crosscheck_seed * 1000003 + p * 1009 + a * 31 + b)
    cross = set(rng.sample(targets, min(args.crosscheck, len(targets))))
    out = open(args.out, 'a')
    j = j_invariant(p, a, b)
    g0 = S.subgroup_generator_value(p, k)
    for t in targets:
        t0 = time.time()
        res = profile(E, t, k, args.dmax, extend_to=args.extend_to)
        rec = dict(prime=p, k=k, a=a % p, b=b % p, j=j, t=t, rhs=E.rhs(t), y_zero=E.rhs(t) == 0,
                   subgroup_generator=g0, generator_hash=res['g_hash'], generator=sorted(res['g'].items()),
                   Z=len(res['Z']), Z_points=res['Z'], rows=res['rows'], falls=res['falls'],
                   status=res['status'], exact_lfd=res['exact_lfd'], local_upper_bound=res['local_upper_bound'])
        if t in cross:
            # reference profile through degree 4 only (dense degree 5 at k = 6 needs ~0.6 GB)
            ref = S.prime_local_profile(res['g'], k, p, res['Z'], g0, max_degree=min(4, args.dmax),
                                        max_columns=20000, max_cells=60_000_000)
            mine = [row['dim'] for row in res['rows']]
            agree = (ref['quotient_ranks'] == mine[:len(ref['quotient_ranks'])]
                     and (ref['exact_lfd'] is None or ref['exact_lfd'] == res['exact_lfd']))
            rec['crosscheck_prime_local_profile'] = dict(agree=agree, exact_lfd=ref['exact_lfd'],
                                                         status=ref['status'], quotient_ranks=ref['quotient_ranks'])
            assert agree, (rec, ref)
        rec['seconds'] = round(time.time() - t0, 2)
        out.write(json.dumps(rec, separators=(',', ':')) + '\n')
        out.flush()
        dims = {row['D']: (row['dim'], row['low']) for row in res['rows']}
        print(p, k, a, b, t, len(res['Z']), dims, res['status'], res['exact_lfd'], rec['seconds'], flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--prime', type=int, required=True)
    r.add_argument('--k', type=int, required=True)
    r.add_argument('--curve', required=True, help='a,b')
    r.add_argument('--dmax', type=int, default=4)
    r.add_argument('--extend-to', type=int, default=None)
    r.add_argument('--targets', default=None, help='comma list (default: every rational abscissa)')
    r.add_argument('--crosscheck', type=int, default=3, help='targets re-run through prime_local_profile')
    r.add_argument('--crosscheck-seed', type=int, default=20261107)
    r.add_argument('--out', required=True)
    args = ap.parse_args()
    if args.cmd == 'run':
        run(args)


if __name__ == '__main__':
    main()
