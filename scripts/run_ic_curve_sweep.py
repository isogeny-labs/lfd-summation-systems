"""Exact LFD sweep: prime, binary ordinary and Koblitz curves; random and symmetric factor bases."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse, csv, hashlib, json, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from csd import ic_curve_sweep as S

SEED = 20261104
SOURCES = ['csd/ic_curve_sweep.py', 'scripts/run_ic_curve_sweep.py', 'csd/ic_arbitrary_bases.py',
           'csd/ic_local_certificate.py', 'csd/ic_boolean_closure.py', 'csd/ic_fall_certificates.py',
           'csd/ic_dreg.py', 'csd/ic_research_campaign.py', 'csd/elim.py', 'csd/_elim.c', 'csd/gf2_ghs.py',
           'csd/field.py', 'csd/gf2n.py']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true',
                        help='Regenerate every case and replay every saved exact closure without changing files.')
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--limit', type=int, help='Evaluate only the first N systems (smoke test; no files written).')
    parser.add_argument('--out', type=Path, default=Path('experiments/ic-curve-sweep'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / args.out
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    cases, systems, skips = S.generate(SEED)
    print(json.dumps(dict(cases=len(cases), systems=len(systems), skips=len(skips),
                          generation_seconds=round(time.time() - started, 1))), flush=True)
    results = {}
    if args.check_only:
        report = json.loads((out / 'verification.json').read_text())
        for section in ['source_sha256', 'data_sha256']:
            for name, expected in report[section].items():
                assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
        saved_cases = [json.loads(line) for line in (out / 'cases.jsonl').read_text().splitlines()]
        saved_systems = [json.loads(line) for line in (out / 'systems.jsonl').read_text().splitlines()]
        assert json.loads(json.dumps(cases)) == saved_cases
        assert json.loads(json.dumps(list(systems.values()))) == saved_systems
        saved_results = {r['system_hash']: r for r in map(json.loads, (out / 'results.jsonl').read_text().splitlines())}
        assert set(saved_results) == set(systems)
    keys = sorted(systems)
    if args.limit:
        keys = keys[:args.limit]
    if not args.check_only and not args.limit:
        for name, rows in [('cases.jsonl', cases), ('systems.jsonl', list(systems.values()))]:
            (out / name).write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in rows))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(S.evaluate, systems[key]): key for key in keys}
        for i, job in enumerate(as_completed(jobs), 1):
            key = jobs[job]
            results[key] = job.result()
            if args.check_only:
                assert results[key] == saved_results[key], key
            if i % 25 == 0 or i == len(keys):
                print(json.dumps(dict(computed=i, systems=len(keys), elapsed_seconds=round(time.time() - started))),
                      flush=True)
    if args.limit:
        print(json.dumps([dict(hash=k[:12], kind=systems[k]['kind'], status=results[k]['status'],
                               lfd=results[k]['exact_lfd']) for k in keys], indent=1))
        return
    ordered = [results[key] for key in sorted(results)]
    groups = S.summarize(cases, results)
    if args.check_only:
        summary = json.loads((out / 'summary.json').read_text())
        assert len(cases) == summary['cases'] and len(systems) == summary['unique_systems']
        assert skips == summary['skips']
        assert dict(Counter(str(r['exact_lfd']) for r in ordered)) == summary['lfd_histogram']
        assert groups == summary['groups']
        print(json.dumps(dict(status='passed: cases regenerated and every exact result replayed',
                              cases=len(cases), systems=len(systems))))
        return
    (out / 'results.jsonl').write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in ordered))
    summary = dict(seed=SEED, cases=len(cases), unique_systems=len(systems), skips=skips,
                   exact_systems=sum(r['exact_lfd'] is not None for r in ordered),
                   capped_systems=sum(r['exact_lfd'] is None for r in ordered),
                   lfd_histogram=dict(Counter(str(r['exact_lfd']) for r in ordered)),
                   by_kind={kind: dict(cases=sum(c['kind'] == kind for c in cases),
                                       systems=len({c['system_hash'] for c in cases if c['kind'] == kind}))
                            for kind in ('prime', 'binary', 'koblitz')},
                   groups=groups,
                   elapsed_seconds=round(time.time() - started),
                   limitations=[
                       'Seeded, structured case grid: no population inference and no claim about other bases.',
                       'Caps retain intervals (lower bound from strict falls, upper bound from the Lean theorem) and never certify an exact value.',
                       'Exact values are ordinary iterated LFDs of the specified presentations; relation yield, independent relation rank and ECDLP cost are not measured.',
                       'A planted decomposition proves one rational decomposition exists; it is not a relation-rate statement.'])
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (out / 'cases.csv').open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['case_id', 'system_hash', 'kind', 'field', 'a', 'b', 'k', 'base_kind', 'symmetric',
                         'target_x', 'target_class', 'polynomial_tuples', 'rational_ordered', 'status', 'lfd',
                         'lower', 'upper'])
        for c in cases:
            r = results[c['system_hash']]
            fieldname = c['prime'] if c['kind'] == 'prime' else f"2^{c['n']}"
            writer.writerow([c['case_id'], c['system_hash'], c['kind'], fieldname, c['a'], c['b'], c['k'],
                             c['base_kind'], c['symmetric'], c['target'][0], c['target_class'],
                             c['polynomial_tuples'], c['rational_ordered_tuples'], r['status'], r['exact_lfd'],
                             r['lfd_lower_bound'], r['lfd_upper_bound']])
    report = dict(status='passed: exact arithmetic, containing spaces and local certificates checked, '
                         'scalar roots and rational decompositions verified',
                  source_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in SOURCES},
                  data_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(out.iterdir()) if p.name not in ('verification.json', 'run.log')})
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ['skips', 'limitations', 'groups']}, indent=2))


if __name__ == '__main__':
    main()
