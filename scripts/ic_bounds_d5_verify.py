"""Independent verification driver for item 3 (degree-5 closure growth).

Uses only csd/ic_bounds_d5_verify.py (own closure engine) plus the instance helpers
csd/ic_dreg.PrimeCurve and csd/ic_curve_sweep.random_rational_target (target drawing, so the
recorded targets can be reproduced), and, for cross-checks only, the existing dense engine
csd/ic_fall_certificates.chain_degree_closure and the independent ordinary-ring engine
experiments/ic-curve-validation/code/lfd_engine.py.

Subcommands (each prints JSON lines and writes experiments/ic-bounds-d5/verify/<name>.json):
  instance  --k K --prime P --a A --b B --seed S [--want]     target, |Z|, eval onto
  closure   --k K --D D ... [--nosplit] [--gens g|all] [--full-sweep] [--minutes M]
  small                                                      18 closures k = 6, 7 vs validation.json
  depth8    [--limit N]                                      saved depth-8 systems
  ordinary                                                   chain closure vs ordinary-ring closure
  d4                                                         dim V_4 at k = 6..13, two instances
"""
import argparse, hashlib, json, os, resource, sys, time, random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from csd import ic_bounds_d5_verify as V
from csd import ic_dreg as I
from csd import ic_curve_sweep as S

OUT = ROOT / 'experiments/ic-bounds-d5/verify'
T0 = time.time()


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['rss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(json.dumps(kw), flush=True)


def instance(p, a, b, seed, k, want):
    """target drawn as in the depth-10 probe: random.Random(seed), ic_curve_sweep.random_rational_target,
    redrawn (want) until the scalar solution set is nonempty; solutions by own enumeration."""
    rng = random.Random(seed)
    E = I.PrimeCurve(p, a, b)
    for attempt in range(200):
        R = S.random_rational_target(E, rng, p)
        Z = V.solutions(p, a % p, b % p, R[0], k)
        if Z or not want:
            break
    g = V.chain_generator(p, a % p, b % p, R[0], k)
    # cross-check the coefficient expansion against ic_dreg.s3_prime_coeffs
    assert V.s3_coeffs(p, a % p, b % p, R[0]) == I.s3_prime_coeffs(E, R[0])
    return dict(target=list(R), Z=Z, g=g, attempts=attempt + 1)


def count(k, D):
    from math import comb
    return sum(comb(2 * k, w) for w in range(D + 1))


def run_closure(g, k, p, D, split=True, gens='all', full_sweep=False, minutes=None, Z=None, tsweep=6):
    t = time.time()
    deadline = t + 60 * minutes if minutes else None
    C = V.Closure(g, k, p, D, split=split, gens=gens, log=log)
    done = C.run(deadline=deadline)
    rec = dict(k=k, D=D, p=p, split=split, gens=gens, columns=C.L.c, finished=done)
    if not done:
        rec['status'] = 'timeout'
        return rec, C
    prof, parts = C.profile()
    rec.update(rank=C.rank(), profile=prof, part_profiles=parts, part_ranks=[P.rank for P in C.parts],
               contains_one=C.contains_one(), rounds=C.rounds, closure_seconds=round(time.time() - t, 1),
               madds=sum(P.madds for P in C.parts), peak_entries=[P.peak_entries for P in C.parts],
               rref_ok=all(P.check_rref() for P in C.parts), basis_sha256=[P.sha() for P in C.parts])
    rec['sweep_random'] = C.sweep_random(t=tsweep)
    if full_sweep:
        rem = (deadline - time.time()) / 60 if deadline else None
        rec['sweep_full'] = C.sweep_full(deadline=deadline)
    if Z is not None:
        rec['points'] = len(Z)
        rec['target'] = count(k, D) - len(Z)
        if Z:
            X = C.random_elements_full(4, seed=7)
            ev = V.evaluate_full(C.L, X, Z, p, k)
            rec['soundness_vanish_on_Z'] = bool(np.all(ev == 0))
            rec['eval_rank_below_D'] = V.eval_rank(C.L, Z, p, k, D - 1)
        rec['deficient'] = rec['rank'] < rec['target']
    rec['seconds'] = round(time.time() - t, 1)
    rec['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    return rec, C


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    src = {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest()
           for n in ('csd/ic_bounds_d5_verify.py', 'scripts/ic_bounds_d5_verify.py')}
    obj = dict(obj, source_sha256=src)
    (OUT / f'{name}.json').write_text(json.dumps(obj, indent=1) + '\n')
    log(stage='saved', file=str(OUT / f'{name}.json'))


INST = dict(unit=(12289, 1, 0, 10, False), p13313=(13313, 1, 0, 11, True), p15361=(15361, 4021, 9917, 12, True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('--k', type=int)
    ap.add_argument('--D', type=int, default=5)
    ap.add_argument('--prime', type=int, default=12289)
    ap.add_argument('--a', type=int, default=1)
    ap.add_argument('--b', type=int, default=0)
    ap.add_argument('--seed', type=int, default=10)
    ap.add_argument('--want', action='store_true')
    ap.add_argument('--nosplit', action='store_true')
    ap.add_argument('--gens', default='all')
    ap.add_argument('--full-sweep', action='store_true')
    ap.add_argument('--minutes', type=float, default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--name', default=None)
    args = ap.parse_args()

    if args.cmd == 'instance':
        inst = instance(args.prime, args.a, args.b, args.seed, args.k, args.want)
        L = V.Layout(args.k, args.D, False)
        er = V.eval_rank(L, inst['Z'], args.prime, args.k, args.D - 1)
        rec = dict(k=args.k, prime=args.prime, a=args.a, b=args.b, seed=args.seed, want=args.want,
                   target=inst['target'], points=len(inst['Z']), attempts=inst['attempts'],
                   eval_rank_below_D=er, onto=er == len(inst['Z']), columns=count(args.k, args.D),
                   Z_bound_2_2k=2 << args.k, generator={f'{e1},{e2}': v for (e1, e2), v in sorted(inst['g'].items())})
        log(**rec)
        save(args.name or f'instance-k{args.k}-p{args.prime}-s{args.seed}', rec)

    elif args.cmd == 'closure':
        inst = instance(args.prime, args.a, args.b, args.seed, args.k, args.want)
        log(stage='instance', target=inst['target'], points=len(inst['Z']))
        rec, C = run_closure(inst['g'], args.k, args.prime, args.D, split=not args.nosplit, gens=args.gens,
                             full_sweep=args.full_sweep, minutes=args.minutes, Z=inst['Z'])
        rec.update(a=args.a, b=args.b, seed=args.seed, want=args.want, instance_target=inst['target'])
        log(stage='result', **{k: v for k, v in rec.items() if k != 'basis_sha256'})
        name = args.name or f'k{args.k}-p{args.prime}-s{args.seed}-D{args.D}' + ('-nosplit' if args.nosplit else '')
        save(name, rec)

    elif args.cmd == 'small':
        val = json.load(open(ROOT / 'experiments/ic-bounds-d5/validation.json'))
        names = {'probe-unit': 'unit', 'proper-13313': 'p13313', 'random-15361': 'p15361'}
        rows = []
        for ref in val['small_depth']:
            p, a, b, seed, want = INST[names[ref['instance']]]
            inst = instance(p, a, b, seed, ref['k'], want)
            res = {}
            for split in (True, False):
                rec, C = run_closure(inst['g'], ref['k'], p, ref['D'], split=split, gens='g', full_sweep=True, Z=inst['Z'])
                res['split' if split else 'nosplit'] = rec
            mine = res['split']
            agree = (mine['rank'] == ref['dense']['dim'] == res['nosplit']['rank']
                     and mine['profile'] == ref['dense']['profile'] == res['nosplit']['profile']
                     and mine['contains_one'] == ref['dense']['contains_one'] == res['nosplit']['contains_one']
                     and len(inst['Z']) == ref['points'])
            row = dict(instance=ref['instance'], k=ref['k'], D=ref['D'], points=len(inst['Z']), ref_points=ref['points'],
                       mine_split=dict(rank=mine['rank'], profile=mine['profile'], one=mine['contains_one'],
                                       sweep=mine['sweep_full'].get('passed'), rsweep=mine['sweep_random']['passed']),
                       mine_nosplit=dict(rank=res['nosplit']['rank'], profile=res['nosplit']['profile'],
                                         sweep=res['nosplit']['sweep_full'].get('passed')),
                       ref_dense=ref['dense'], agree=agree)
            log(**row)
            rows.append(row)
        save('small', dict(rows=rows, all_agree=all(r['agree'] for r in rows)))

    elif args.cmd == 'depth8':
        systems = [json.loads(l) for l in open(ROOT / 'experiments/ic-chain-depth8/systems.jsonl')]
        results = {json.loads(l)['system_hash']: json.loads(l) for l in open(ROOT / 'experiments/ic-chain-depth8/results.jsonl')}
        rows = []
        for s in systems[:args.limit]:
            k, p = s['k'], s['prime']
            N = 1 << k
            g = {}
            for mask, c in s['generators'][0]:
                key = (mask & (N - 1), mask >> k)
                g[key] = (g.get(key, 0) + c) % p
            sym = all(g.get((e2, e1), 0) == v for (e1, e2), v in g.items())
            Z = [(m & (N - 1), m >> k) for m in s['coordinate_solutions']]
            cert = results[s['system_hash']]['certificate_checks'][0]
            r4, _ = run_closure(g, k, p, 4, split=sym, gens='all', full_sweep=True, Z=Z)
            r5, _ = run_closure(g, k, p, 5, split=sym, gens='all', full_sweep=True, Z=Z)
            row = dict(system_hash=s['system_hash'], k=k, prime=p, points=len(Z), symmetric=sym,
                       dim_V4=r4['rank'], saved_dim_V4=cert['previous_dimension'],
                       dim_V5=r5['rank'], saved_rank5=cert['rank'], target5=r5['target'], saved_target5=cert['target_rank'],
                       low5=r5['profile'][4], saved_low5=cert['low_dimension'],
                       sweeps_ok=bool(r4['sweep_full']['passed'] and r5['sweep_full']['passed']),
                       sound=r5.get('soundness_vanish_on_Z'), onto=r5.get('eval_rank_below_D') == len(Z) if Z else None)
            row['agree'] = (row['dim_V4'] == row['saved_dim_V4'] and row['dim_V5'] == row['saved_rank5']
                            and row['target5'] == row['saved_target5'] and row['low5'] == row['saved_low5'])
            log(**row)
            rows.append(row)
        save('depth8', dict(rows=rows, all_agree=all(r['agree'] for r in rows)))

    elif args.cmd == 'ordinary':
        sys.path.insert(0, str(ROOT / 'experiments/ic-curve-validation/code'))
        import lfd_engine as LE
        rows = []
        for name, k, Ds in (('unit', 4, (3, 4, 5)), ('p13313', 4, (3, 4, 5)), ('unit', 5, (3, 4, 5)), ('p15361', 5, (4, 5))):
            p, a, b, seed, want = INST[name]
            inst = instance(p, a, b, seed, k, want)
            N = 1 << k
            gen_items = [((e1 & (N - 1)) | (e2 << k), v) for (e1, e2), v in inst['g'].items()]
            # lfd_engine uses variable i*k + j for u_{i,j}: mask bit b -> variable b
            inputs, G = LE.prime_inputs(gen_items, k, p)
            for D in Ds:
                t = time.time()
                res = LE.closure(inputs, 2 * k, D, p)
                sp, Vb = res['space'], res['V']
                # image under NF: ordinary monomial -> (sum_j e_{1,j} 2^j mod N, sum_j e_{2,j} 2^j mod N)
                L = V.Layout(k, D, False)
                img = np.zeros((Vb.rk, L.c), dtype=np.float64)
                cols = []
                for m in sp.cols:
                    a1 = sum(m[j] << j for j in range(k)) % N
                    a2 = sum(m[k + j] << j for j in range(k)) % N
                    cols.append(int(L.index(np.array([a1]), np.array([a2]))[0]))
                cols = np.array(cols)
                E = Vb.E
                for j in range(sp.n):
                    img[:, cols[j]] += E[:, j]
                img = np.mod(img, p)
                P = V.Part(L.w, p)
                for a0 in range(0, img.shape[0], 512):
                    P.add(img[a0:a0 + 512])
                rec, C = run_closure(inst['g'], k, p, D, split=True, gens='g', full_sweep=True)
                row = dict(instance=name, k=k, D=D, ordinary_dim=res['dim'], ordinary_one=res['contains_one'],
                           image_rank=P.rank, chain_rank=rec['rank'], chain_one=rec['contains_one'],
                           agree=P.rank == rec['rank'] and res['contains_one'] == rec['contains_one'],
                           seconds=round(time.time() - t, 1))
                log(**row)
                rows.append(row)
        save('ordinary', dict(rows=rows, all_agree=all(r['agree'] for r in rows)))

    elif args.cmd == 'd4':
        rows = []
        for k in range(6, 14):
            p = 12289 if k <= 12 else 40961
            for a, b, seed, want in ((1, 0, 10, False), (4021, 9917, 12, True)):
                inst = instance(p, a, b, seed, k, want)
                rec, C = run_closure(inst['g'], k, p, 4, split=True, gens='all', full_sweep=True, Z=inst['Z'])
                row = dict(k=k, prime=p, a=a, b=b, seed=seed, target=inst['target'], points=len(inst['Z']),
                           dim_V4=rec['rank'], fit=4 * k * k + 172 * k - 571, low3=rec['profile'][3], fit_low3=4 * k + 87,
                           profile=rec['profile'], one=rec['contains_one'], sweep_full=rec['sweep_full'],
                           sweep_random=rec['sweep_random']['passed'], sound=rec.get('soundness_vanish_on_Z'),
                           seconds=rec['seconds'])
                row['match'] = row['dim_V4'] == row['fit'] and row['low3'] == row['fit_low3'] and rec['sweep_full']['passed']
                log(**row)
                rows.append(row)
        save('d4', dict(rows=rows, all_match=all(r['match'] for r in rows)))


if __name__ == '__main__':
    main()
