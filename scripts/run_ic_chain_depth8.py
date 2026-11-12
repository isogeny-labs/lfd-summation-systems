"""Prime chain arm at depth k >= 8 with the shift-structured closure engine (plus k = 7 controls)."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse, hashlib, json, random, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from csd import ic_dreg as I
from csd import ic_curve_sweep as S
from csd import ic_chain_shift_closure as SC

SEED = 20261104
SOURCES = ['csd/ic_chain_shift_closure.py', 'scripts/run_ic_chain_depth8.py', 'csd/ic_curve_sweep.py',
           'csd/ic_fall_certificates.py', 'csd/ic_dreg.py', 'csd/elim.py', 'csd/_elim.c']


def generate(seed=SEED, primes=(257, 769, 7681), depths=(8,), controls=2):
    rng = random.Random(seed)
    cases, systems = [], {}
    for p in primes:
        curves = [(1, 0)] if (4 * 1 + 0) % p else []
        while len(curves) < 2:
            a, b = rng.randrange(1, p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p and (a, b) not in curves:
                curves.append((a, b))
        for a, b in curves:
            E = I.PrimeCurve(p, a, b)
            for k in depths:
                if (p - 1) % (1 << k):
                    continue
                N = 1 << k
                g0 = S.subgroup_generator_value(p, k)
                H = [pow(g0, e, p) for e in range(N)]
                bases = [('subgroup', 1)]
                if N < p - 1:
                    while len(bases) < 2:
                        c = rng.randrange(2, p)
                        if c not in H:
                            bases.append(('random-coset', c))
                for base_kind, scale in bases:
                    values = [scale * h % p for h in H]
                    targets = [(S.random_rational_target(E, rng, p), 'random-rational', None)]
                    plant = S.planted_target(E, values, rng)
                    if plant:
                        targets.append((plant[0], 'planted-rational-decomposition', plant[1]))
                    for target, kind, planted in targets:
                        case, system = S.prepare_prime(E, k, base_kind, scale, target, kind, planted)
                        case['arm'] = f'depth-{k}' if k >= 8 else 'control'
                        cases.append(case)
                        systems.setdefault(system['system_hash'], system)
    # k = 7 controls taken from the main sweep's grid (same seed logic, first two systems of p = 7681)
    E = I.PrimeCurve(7681, 1, 0)
    for _ in range(controls):
        target = S.random_rational_target(E, rng, 7681)
        case, system = S.prepare_prime(E, 7, 'subgroup', 1, target, 'random-rational', None)
        case['arm'] = 'control'
        cases.append(case)
        systems.setdefault(system['system_hash'], system)
    return cases, systems


PRODUCT_CAP = 16000
SHIFT_MAX_CELLS = 200_000_000


def evaluate(system):
    generator = {int(m): int(c) for m, c in system['generators'][0]}
    started = time.time()
    result = SC.shift_profile(generator, system['k'], system['prime'], system['coordinate_solutions'],
                              system['generator_value'], dense_degree=4, certificate_degrees=(5, 6),
                              product_cap=PRODUCT_CAP if system['k'] >= 8 else 9000,
                              shift_max_cells=SHIFT_MAX_CELLS)
    if result['exact_lfd'] is not None:
        assert result['exact_lfd'] <= 2 * system['k'] + 1
    result['seconds'] = round(time.time() - started, 1)
    return dict(system_hash=system['system_hash'], **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--out', type=Path, default=Path('experiments/ic-chain-depth8'))
    parser.add_argument('--depths', default='8', help='comma-separated depths, e.g. 9')
    parser.add_argument('--primes', default='257,769,7681', help='comma-separated primes with 2^k | p-1')
    parser.add_argument('--product-cap', type=int, default=16000)
    parser.add_argument('--shift-max-cells', type=int, default=200_000_000)
    args = parser.parse_args()
    global PRODUCT_CAP, SHIFT_MAX_CELLS
    PRODUCT_CAP = args.product_cap
    SHIFT_MAX_CELLS = args.shift_max_cells
    root = Path(__file__).resolve().parents[1]
    out = root / args.out
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    depths = tuple(int(d) for d in args.depths.split(','))
    primes = tuple(int(q) for q in args.primes.split(','))
    cases, systems = generate(primes=primes, depths=depths)
    print(json.dumps(dict(cases=len(cases), systems=len(systems))), flush=True)
    if args.check_only:
        report = json.loads((out / 'verification.json').read_text())
        for section in ['source_sha256', 'data_sha256']:
            for name, expected in report[section].items():
                assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
        saved = {r['system_hash']: r for r in map(json.loads, (out / 'results.jsonl').read_text().splitlines())}
        assert set(saved) == set(systems)
    else:
        for name, rows in [('cases.jsonl', cases), ('systems.jsonl', list(systems.values()))]:
            (out / name).write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in rows))
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(evaluate, systems[key]): key for key in sorted(systems)}
        for i, job in enumerate(as_completed(jobs), 1):
            key = jobs[job]
            results[key] = job.result()
            r = results[key]
            print(json.dumps(dict(done=i, total=len(jobs), k=systems[key]['k'], prime=systems[key]['prime'],
                                  status=r['status'], lfd=r['exact_lfd'], falls=r['strict_falls'],
                                  seconds=r['seconds'], elapsed=round(time.time() - started))), flush=True)
            if args.check_only:
                a = {k: v for k, v in r.items() if k != 'seconds'}
                b = {k: v for k, v in saved[key].items() if k != 'seconds'}
                assert a == b, key
    if args.check_only:
        print(json.dumps(dict(status='passed: every result replayed')))
        return
    ordered = [results[key] for key in sorted(results)]
    (out / 'results.jsonl').write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in ordered))
    summary = dict(seed=SEED, cases=len(cases), unique_systems=len(systems),
                   elapsed_seconds=round(time.time() - started),
                   by_depth={str(k): dict(systems=len([s for s in systems.values() if s['k'] == k]),
                                          status=dict(Counter(results[h]['status'] for h, s in systems.items() if s['k'] == k)),
                                          lfd=dict(Counter(str(results[h]['exact_lfd']) for h, s in systems.items() if s['k'] == k)),
                                          falls=dict(Counter(str(results[h]['strict_falls']) for h, s in systems.items() if s['k'] == k)))
                             for k in sorted({s['k'] for s in systems.values()})},
                   by_base={bk: dict(Counter(str(results[c['system_hash']]['exact_lfd']) for c in cases if c['base_kind'] == bk and c['k'] >= 8))
                            for bk in ('subgroup', 'random-coset')},
                   engine='shift-structured closure: exact dense profile through degree 4, once-per-pivot shift closure at degrees 5 and 6',
                   limitations=['Small seeded arm; no population inference.',
                                'An exact value needs the certificate at degree 5 (then V_4 exact decides the fall at 5); a certificate only at degree 6 gives an interval.'])
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    report = dict(status='passed: exact dense profiles through degree 4 with exclusion certificates; shift-closure rank certificates at degree 5/6',
                  source_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in SOURCES},
                  data_sha256={str(q.relative_to(root)): hashlib.sha256(q.read_bytes()).hexdigest()
                               for q in sorted(out.iterdir()) if q.name not in ('verification.json', 'run.log')})
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(summary['by_depth'], indent=1))


if __name__ == '__main__':
    main()
