#!/usr/bin/env python3
"""Item 1 driver: exact certified LFD of three-summand prime subgroup-chain systems.

  python3 scripts/ic_bounds_m3_run.py generate            # writes cases.jsonl (seeded grid)
  python3 scripts/ic_bounds_m3_run.py run --k 1 2 3 [--encodings direct split m2] [--budget-min 22]
  python3 scripts/ic_bounds_m3_run.py summarize           # summary.json + table.md

Every case: prime p (2^4 | p-1), curve (j=0, j=1728, random), target R (order two, random rational,
planted P1+P2+P3 with x_i in H), depth k.  Encodings:
  direct : {chains, NF(S_4(x1,x2,x3,t))}, Z enumerated over H^3;
  split  : {chains, NF(S_3(x1,x2,z)), NF(S_3(z,x3,t))}, z an ordinary variable (I_free) and
           additionally z^p - z (I_field), both certified from the same closures (D < p);
  m2     : {chains, NF(S_3(x1,x2,t))} at the same (p, curve, t, k), through
           csd.ic_curve_sweep.prime_local_profile (the m = 2 reference).
Results are appended to results.jsonl (one line per case and encoding) and the run is resumable;
each run stops after --budget-min minutes.
"""
import argparse
import hashlib
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_m3 as M  # noqa: E402
from csd import ic_curve_sweep as S  # noqa: E402

OUT = ROOT / 'experiments' / 'ic-bounds-m3'
SEED = 20261107
PRIMES = (17, 97, 257, 641, 7681)
KS = (1, 2, 3, 4)


def curves_for(p, rng):
    out = [('j0', 0, 1), ('j1728', 1, 0)]
    while True:
        b = rng.randrange(2, p)
        if (27 * b * b) % p and ('j0', 0, b) not in out:
            out.append(('j0', 0, b))
            break
    while True:
        a = rng.randrange(2, p)
        if (4 * a ** 3) % p and ('j1728', a, 0) not in out:
            out.append(('j1728', a, 0))
            break
    while len(out) < 6:
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        if (4 * a ** 3 + 27 * b * b) % p and ('random', a, b) not in out:
            out.append(('random', a, b))
    return out


def planted(E, H, rng):
    pts = [P for x in H if (P := E.lift_x(x)) is not None]
    if not pts:
        return None
    for _ in range(500):
        P = [rng.choice(pts) for _ in range(3)]
        P = [E.neg(Q) if rng.randrange(2) else Q for Q in P]
        R = E.add(E.add(P[0], P[1]), P[2])
        if R is not None:
            return R, [list(Q) for Q in P]
    return None


def generate():
    rng = random.Random(SEED)
    cases, skips = [], []
    for p in PRIMES:
        assert (p - 1) % 16 == 0
        for cls, a, b in curves_for(p, rng):
            E = I.PrimeCurve(p, a, b)
            base_targets = []
            two = S.order_two_prime_target(E)
            if two is not None:
                base_targets.append((two, 'order-two', None))
            else:
                skips.append(dict(p=p, a=a, b=b, reason='no rational 2-torsion point'))
            used = {two[0]} if two else set()
            for _ in range(2):
                R = S.random_rational_target(E, rng, p, exclude_x=used)
                used.add(R[0])
                base_targets.append((R, 'random-rational', None))
            for k in KS:
                g0, H = M.subgroup(p, k)
                targets = list(base_targets)
                pl = planted(E, H, rng)
                if pl is None:
                    skips.append(dict(p=p, a=a, b=b, k=k, reason='no point with x in H: no planted target'))
                else:
                    targets.append((pl[0], 'planted-P1+P2+P3', pl[1]))
                for R, tclass, plant in targets:
                    assert E.on_curve(R)
                    case = dict(prime=p, a=a, b=b, curve_class=cls, k=k, target=list(R), target_class=tclass,
                                planted=plant, subgroup_generator=g0)
                    case['case_id'] = hashlib.sha256(json.dumps(case, sort_keys=True).encode()).hexdigest()[:16]
                    cases.append(case)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / 'cases.jsonl', 'w') as f:
        for c in cases:
            f.write(json.dumps(c, sort_keys=True) + '\n')
    (OUT / 'skips.json').write_text(json.dumps(dict(seed=SEED, skips=skips), indent=1))
    print(len(cases), 'cases;', len(skips), 'skips')


def geometry(case):
    """Scalar and rational counts; asserts every rational decomposition lies in Z."""
    p, k = case['prime'], case['k']
    E = I.PrimeCurve(p, case['a'], case['b'])
    R = tuple(case['target'])
    g0, H = M.subgroup(p, k)
    g = M.direct_generator(E, R[0], k)
    Z = M.direct_points(E, R[0], k, g)
    dec = M.rational_decompositions(E, H, R)
    Zset = set(Z)
    idx = {h: e for e, h in enumerate(H)}
    for xs in dec:
        pt = idx[xs[0]] | (idx[xs[1]] << k) | (idx[xs[2]] << (2 * k))
        assert pt in Zset, 'rational decomposition outside Z(S_4)'
    if case['planted']:
        P = [tuple(Q) for Q in case['planted']]
        assert E.add(E.add(P[0], P[1]), P[2]) == R and all(Q[0] in idx for Q in P)
        assert tuple(Q[0] for Q in P) in set(dec)
    return dict(scalar_solutions_H3=len(Z), rational_ordered_triples=len(dec),
                rational_unordered=len({tuple(sorted(x)) for x in dec}), H_size=len(H),
                factor_base_points=S.I.prime_factor_base_size(E, H))


def run_m2(case):
    p, k = case['prime'], case['k']
    E = I.PrimeCurve(p, case['a'], case['b'])
    t = case['target'][0]
    g0 = S.subgroup_generator_value(p, k)
    gen = S.prime_coset_generator(E, t, k, 1, 1)
    pts = S.prime_points(E, k, g0, 1, 1, t)
    res = S.prime_local_profile(gen, k, p, pts, g0, max_degree=2 * k + 1)
    return dict(encoding='m2', generator_hash=M.system_hash(dict(p=p, k=k, m=2, g=sorted(gen.items()))),
                results={'m2': dict(status=res['status'], exact_lfd=res['exact_lfd'],
                                    lfd_lower_bound=res['lfd_lower_bound'], lfd_upper_bound=res['lfd_upper_bound'],
                                    quotient_dimension=len(pts))},
                strict_falls=res['strict_falls'], quotient_ranks=res['quotient_ranks'],
                local_checks=res.get('local_checks'), engine='ic_curve_sweep.prime_local_profile')


def run(args):
    cases = [json.loads(l) for l in open(OUT / 'cases.jsonl')]
    done = set()
    rpath = OUT / args.results
    for path in sorted(OUT.glob('results*.jsonl')):
        for l in open(path):
            r = json.loads(l)
            done.add((r['case_id'], r['encoding']))
    t_start = time.time()
    n = 0
    order = [c for c in cases if c['k'] in args.k and (not args.primes or c['prime'] in args.primes)]
    if args.select:
        keep = set(json.loads((OUT / args.select).read_text())['case_ids'])
        order = [c for c in order if c['case_id'] in keep]
    with open(rpath, 'a') as f:
        for c in order:
            for enc in args.encodings:
                if (c['case_id'], enc) in done:
                    continue
                if time.time() - t_start > args.budget_min * 60:
                    print('budget reached; rerun to resume', flush=True)
                    return
                E = I.PrimeCurve(c['prime'], c['a'], c['b'])
                t0 = time.time()
                if enc == 'm2':
                    res = run_m2(c)
                else:
                    res = M.exact_profile(enc, E, c['target'][0], c['k'],
                                          max_degree=args.max_degree, max_columns=args.max_columns)
                res['geometry'] = geometry(c) if enc == 'direct' else None
                res.update(case_id=c['case_id'], seconds=round(time.time() - t0, 2))
                f.write(json.dumps(res, sort_keys=True) + '\n')
                f.flush()
                n += 1
                brief = {k: (v['status'], v['exact_lfd'], v['quotient_dimension']) for k, v in res['results'].items()}
                print(c['prime'], c['curve_class'], c['a'], c['b'], c['k'], c['target_class'], enc, brief,
                      res['seconds'], flush=True)
    print('done', n, 'new results', flush=True)


def select_k4_split():
    """Deterministic k = 4 split subset (the split closure at D = 6 has 5282 columns, ~2 min per
    case): per prime, the curves (0, 1) [j = 0], (1, 0) [j = 1728] and the first random curve; for each
    the first random-rational target and the planted target, plus the order-two target of (1, 0)."""
    cases = [json.loads(l) for l in open(OUT / 'cases.jsonl')]
    keep = []
    for p in PRIMES:
        sub = [c for c in cases if c['prime'] == p and c['k'] == 4]
        rnd = [c for c in sub if c['curve_class'] == 'random']
        first_random = (rnd[0]['a'], rnd[0]['b']) if rnd else None
        for ab in [(0, 1), (1, 0), first_random]:
            cs = [c for c in sub if (c['a'], c['b']) == ab]
            rr = [c for c in cs if c['target_class'] == 'random-rational'][:1]
            pl = [c for c in cs if c['target_class'] == 'planted-P1+P2+P3']
            keep += [c['case_id'] for c in rr + pl]
            if ab == (1, 0):
                keep += [c['case_id'] for c in cs if c['target_class'] == 'order-two']
    (OUT / 'k4-split-selection.json').write_text(json.dumps(dict(rule=select_k4_split.__doc__.strip(),
                                                                 case_ids=keep), indent=1))
    print(len(keep), 'cases selected')


def select_k4_split2():
    """Second k = 4 split subset (curve coverage): per prime, the three curves not in the first
    subset (j = 0 with random b, j = 1728 with random a, second random curve), first random-rational
    target each."""
    cases = [json.loads(l) for l in open(OUT / 'cases.jsonl')]
    first = set(json.loads((OUT / 'k4-split-selection.json').read_text())['case_ids'])
    keep = []
    for p in PRIMES:
        sub = [c for c in cases if c['prime'] == p and c['k'] == 4]
        curves = []
        for c in sub:
            if (c['a'], c['b']) not in curves:
                curves.append((c['a'], c['b']))
        for ab in curves:
            cs = [c for c in sub if (c['a'], c['b']) == ab]
            if any(c['case_id'] in first for c in cs):
                continue
            keep += [c['case_id'] for c in cs if c['target_class'] == 'random-rational'][:1]
    (OUT / 'k4-split-selection2.json').write_text(json.dumps(dict(rule=select_k4_split2.__doc__.strip(),
                                                                  case_ids=keep), indent=1))
    print(len(keep), 'cases selected')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['generate', 'run', 'select-k4-split', 'select-k4-split2'])
    ap.add_argument('--k', type=int, nargs='+', default=[1, 2, 3])
    ap.add_argument('--primes', type=int, nargs='*', default=None)
    ap.add_argument('--encodings', nargs='+', default=['direct', 'split', 'm2'])
    ap.add_argument('--budget-min', type=float, default=22)
    ap.add_argument('--max-degree', type=int, default=None)
    ap.add_argument('--max-columns', type=int, default=None)
    ap.add_argument('--results', default='results.jsonl', help='output file name (results*.jsonl)')
    ap.add_argument('--select', default=None, help='JSON file in the output directory with case_ids')
    args = ap.parse_args()
    if args.cmd == 'generate':
        generate()
    elif args.cmd == 'select-k4-split':
        select_k4_split()
    elif args.cmd == 'select-k4-split2':
        select_k4_split2()
    else:
        run(args)


if __name__ == '__main__':
    main()
