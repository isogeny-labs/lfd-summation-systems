#!/usr/bin/env python3
"""Bounded matched GB benchmark and exhaustive small-field correctness audit.

Run --help. No numerical result is assumed; workers are isolated for peak RSS.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import platform
import random
import resource
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_dreg as I
from csd import semaev_matched as M

OUT = ROOT / 'experiments/semaev-matched-20261124'
GRID = [(3, 1), (3, 2), (4, 2), (5, 2), (5, 3), (7, 2)]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def curve(n, a, b):
    return I.BinaryCurve(I.binary_field(n), a, b)


def audit(out):
    counts = Counter()
    witnesses = {}
    for n in (2, 3, 4):
        q = 1 << n
        for a, b in product(range(q) if n < 4 else (0, 1),
                            range(1, q) if n < 4 else (1, 3)):
            E = curve(n, a, b)
            V = list(range(q))
            oracle = M.decomposition_map(E, V)
            roots = M.root_table(E.F, b)
            counts['curves'] += 1
            for R in M.rational_points(E):
                sets = M.coordinate_sets(E, V, R[0], roots)
                want = oracle.get(R, set())
                for key in ('pdp_direct', 'pdp_split'):
                    if sets[key] != want:
                        raise AssertionError((n, a, b, R, key, sets[key] ^ want))
                missing = sets['raw_direct'] - sets['raw_finite']
                extra = sets['raw_finite'] - sets['raw_direct']
                counts['targets'] += 1
                counts['coordinate_triples'] += q ** 3
                counts['raw_direct_zeros'] += len(sets['raw_direct'])
                counts['raw_split_zeros'] += len(sets['raw_finite'])
                counts['raw_split_missing'] += len(missing)
                counts['raw_split_extra'] += len(extra)
                counts['nonliftable_direct_zeros'] += len(sets['raw_direct'] - sets['pdp_direct'])
                counts['pdp_affine_missing'] += len(want - sets['pdp_finite'])
                counts['pdp_corrected_mismatches'] += int(sets['pdp_split'] != want)
                for t in sorted(missing):
                    x, y, z = t
                    coeffs = [M.quadratic_coefficients(E.F, b, x, y),
                              M.quadratic_coefficients(E.F, b, z, R[0])]
                    label = ('degree_drop' if not all(p[-1] for p in coeffs)
                             else 'quadratic_gcd_no_k_root')
                    counts[label] += 1
                    if label not in witnesses:
                        witnesses[label] = dict(n=n, modulus=E.F.modulus, a=a, b=b,
                                                target=R, triple=t, quadratics=coeffs,
                                                gcd_degree=M.gcd_degree(E.F, *coeffs),
                                                actual_decomposition=t in want)
                if want - sets['pdp_finite'] and 'pdp_infinity' not in witnesses:
                    t = min(want - sets['pdp_finite'])
                    witnesses['pdp_infinity'] = dict(n=n, a=a, b=b, target=R, triple=t)
    result = dict(status='passed', counts=dict(counts), witnesses=witnesses,
                  scope='Exhaustive for the specified finite grid, not a universal proof')
    write(out / 'correctness.json', result)
    print(json.dumps(result, indent=2))


def instances(out):
    cases, missing = [], []
    for n, ell in GRID:
        for seed in (0, 1, 2):
            rng = random.Random((n << 16) + (ell << 8) + seed)
            E = curve(n, rng.randrange(2), rng.randrange(1, 1 << n))
            basis = [1 << i for i in range(ell)]
            oracle = M.decomposition_map(E, I.span(basis))
            targets = M.rational_points(E)
            rng.shuffle(targets)
            for category in ('yes', 'no'):
                selected = next((R for R in targets
                                 if bool(oracle.get(R)) == (category == 'yes')), None)
                if selected is None:
                    missing.append(dict(n=n, ell=ell, seed=seed, category=category))
                    continue
                case = dict(id=f'n{n}-l{ell}-s{seed}-{category}', n=n, ell=ell,
                            seed=seed, a=E.a, b=E.b, modulus=E.F.modulus,
                            target=selected, basis=basis, category=category,
                            oracle=sorted(oracle.get(selected, set())))
                cases.append(case)
    result = dict(cases=cases, missing_categories=missing)
    write(out / 'instances.json', result)
    return result


class GBTimeout(Exception):
    pass


def worker(payload):
    import sympy as sp
    from sympy.polys import groebnertools as gt
    case, formulation = payload['case'], payload['formulation']
    E = curve(case['n'], case['a'], case['b'])
    begin = time.perf_counter()
    eqs, nv = M.build_system(E, case['basis'], case['target'][0], formulation)
    variables = sp.symbols(f'z0:{nv}')
    exprs = [sp.Add(*(sp.Mul(*(variables[j] for j in range(nv) if m >> j & 1))
                       for m in eq)) for eq in eqs]
    exprs += [x*x-x for x in variables]
    setup = time.perf_counter() - begin
    diag = dict(spoly_calls=0, max_spoly_degree=0, max_spoly_terms=0,
                matrix_rows=None, matrix_columns=None, matrix_rank=None,
                matrix_sparsity=None)
    if payload.get('diagnostic'):
        original = gt.spoly
        def counted(p, q, ring):
            result = original(p, q, ring)
            diag['spoly_calls'] += 1
            diag['max_spoly_terms'] = max(diag['max_spoly_terms'], len(result))
            diag['max_spoly_degree'] = max(diag['max_spoly_degree'],
                                          max((sum(m) for m in result), default=0))
            return result
        gt.spoly = counted
    def alarm(_sig, _frame):
        raise GBTimeout()
    oldhandler = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, payload.get('gb_timeout', 5))
    start = time.perf_counter()
    try:
        G = sp.groebner(exprs, *variables, modulus=2, order='grevlex', method='buchberger')
    except GBTimeout:
        print(json.dumps(dict(status='gb_timeout', setup_seconds=setup,
                              gb_seconds=time.perf_counter()-start,
                              input=M.polynomial_stats(eqs, nv), diagnostics=diag,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
        return
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, oldhandler)
    gb = time.perf_counter() - start
    gb_degree = max((p.total_degree() for p in G.polys), default=0)
    if payload.get('diagnostic'):
        print(json.dumps(dict(status='ok', diagnostics=diag, gb_degree=gb_degree)))
        return

    start = time.perf_counter()
    lex = G if list(G) == [1] else G.fglm(order='lex')
    # Reverse-variable search with exact pruning using the lex basis. No oracle
    # enters this phase. Every full Boolean assignment satisfying G is recovered.
    constraints = []
    for p in lex.polys:
        terms = [sum(1 << j for j, e in enumerate(m) if e) for m, c in p.terms() if int(c) & 1]
        support = 0
        for m in terms:
            support |= m
        constraints.append((support, terms))
    solutions = []
    def visit(j, mask, assigned):
        if any((support & assigned) == support and M.eval_boolean(poly, mask)
               for support, poly in constraints):
            return
        if j < 0:
            solutions.append(mask)
            return
        visit(j-1, mask, assigned | (1 << j))
        visit(j-1, mask | (1 << j), assigned | (1 << j))
    visit(nv-1, 0, 0)
    projected = {M.decode(m, case['basis']) for m in solutions}
    infinity_start = time.perf_counter()
    extra = M.infinity_solutions(E, case['basis'], case['target'][0]) if formulation == 'split' else set()
    projected |= extra
    infinity_seconds = time.perf_counter() - infinity_start
    recovery = time.perf_counter() - start
    start = time.perf_counter()
    want = {tuple(t) for t in case['oracle']}
    if projected != want:
        raise AssertionError(dict(case=case['id'], formulation=formulation,
                                  missing=sorted(want-projected), extra=sorted(projected-want)))
    verify = time.perf_counter() - start
    print(json.dumps(dict(status='ok', setup_seconds=setup, gb_seconds=gb,
                          recovery_seconds=recovery, infinity_seconds=infinity_seconds,
                          verification_seconds=verify, total_seconds=setup+gb+recovery,
                          peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          input=M.polynomial_stats(eqs, nv), gb_degree=gb_degree,
                          gb_polynomials=len(G.polys), gb_terms=sum(len(p.terms()) for p in G.polys),
                          full_boolean_solutions=len(solutions), projected_solutions=len(projected),
                          infinity_branch_solutions=len(extra))))


def invoke(payload):
    start = time.perf_counter()
    try:
        run = subprocess.run([sys.executable, str(Path(__file__).resolve()), 'worker'],
                             input=json.dumps(payload), text=True, capture_output=True, timeout=20,
                             env={**os.environ, 'PYTHONHASHSEED': '0', 'OPENBLAS_NUM_THREADS': '1'})
        result = (json.loads(run.stdout) if run.returncode == 0 else
                  dict(status='error', stderr=run.stderr[-6000:], returncode=run.returncode))
    except subprocess.TimeoutExpired:
        result = dict(status='process_timeout')
    result['process_wall_seconds'] = time.perf_counter() - start
    return result


def benchmark(out):
    frozen = json.loads((out / 'instances.json').read_text())
    inputs = ['csd/semaev_matched.py', 'scripts/run_semaev_matched.py',
              'csd/ic_dreg.py', 'csd/gf2_ghs.py']
    import sympy
    write(out / 'environment.json', dict(python=sys.version, sympy=sympy.__version__,
          platform=platform.platform(), cpu=Path('/proc/cpuinfo').read_text().split('model name')[1].split('\n')[0].strip(),
          source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
          protocol_sha256=hashlib.sha256((out/'PROTOCOL.md').read_bytes()).hexdigest(),
          instances_sha256=hashlib.sha256((out/'instances.json').read_bytes()).hexdigest(),
          gb_timeout_seconds=5, process_timeout_seconds=20,
          matrix_diagnostics='not applicable to Buchberger; no matrices constructed'))
    began = time.monotonic()
    jobs = []
    for ci, case in enumerate(frozen['cases']):
        for rep in range(3):
            order = ['original', 'split'] if (ci+rep) % 2 == 0 else ['split', 'original']
            jobs += [(case, f, rep, False) for f in order]
        jobs += [(case, f, None, True) for f in ('original', 'split')]
    result_path = out / 'runs.jsonl'
    if result_path.exists():
        raise FileExistsError('Refusing to overwrite an existing benchmark record; choose a fresh --out')
    with result_path.open('x') as stream:
        for case, form, rep, diagnostic in jobs:
            if time.monotonic() - began >= 1200:
                result = dict(status='unrun_budget')
            else:
                result = invoke(dict(case=case, formulation=form, diagnostic=diagnostic))
            row = dict(case_id=case['id'], formulation=form, repetition=rep,
                       diagnostic=diagnostic, **result)
            stream.write(json.dumps(row, sort_keys=True)+'\n')
            stream.flush()
            print(case['id'], form, 'diagnostic' if diagnostic else rep, row['status'],
                  round(row.get('gb_seconds', 0), 4), flush=True)
    summarize(out)


def summarize(out):
    cases = {c['id']:c for c in json.loads((out/'instances.json').read_text())['cases']}
    runs = [json.loads(line) for line in (out/'runs.jsonl').read_text().splitlines()]
    pairs = defaultdict(dict)
    for row in runs:
        if not row['diagnostic']:
            pairs[row['case_id']].setdefault(row['formulation'], []).append(row)
    rows = []
    for cid, methods in pairs.items():
        c = cases[cid]
        row = {k:c[k] for k in ('id', 'n', 'ell', 'seed', 'category')}
        for method, samples in methods.items():
            row[method] = dict(statuses=[s['status'] for s in samples])
            if len(samples) == 3 and all(s['status']=='ok' for s in samples):
                row[method].update({key:statistics.median(s[key] for s in samples)
                                    for key in ('setup_seconds','gb_seconds','recovery_seconds',
                                                'total_seconds','peak_rss_kib')})
                row[method]['input'] = samples[0]['input']
        if all('gb_seconds' in row.get(m,{}) for m in ('original','split')):
            row['gb_speedup'] = row['original']['gb_seconds']/row['split']['gb_seconds']
            row['total_speedup'] = row['original']['total_seconds']/row['split']['total_seconds']
        rows.append(row)
    write(out/'summary.json', dict(cases=rows, statuses=dict(Counter(r['status'] for r in runs)),
          complete_pairs=sum('gb_speedup' in r for r in rows),
          split_gb_wins=sum(r.get('gb_speedup',0)>1 for r in rows),
          split_gb_losses=sum(0<r.get('gb_speedup',0)<1 for r in rows)))
    print('summary:', out/'summary.json')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['audit','instances','benchmark','summary','worker'])
    p.add_argument('--out', type=Path, default=OUT)
    a = p.parse_args()
    if a.command == 'worker':
        worker(json.load(sys.stdin))
    elif a.command == 'audit':
        audit(a.out)
    elif a.command == 'instances':
        result=instances(a.out)
        print(len(result['cases']), 'cases;', len(result['missing_categories']), 'missing categories')
    elif a.command == 'benchmark':
        benchmark(a.out)
    else:
        summarize(a.out)


if __name__ == '__main__':
    main()
