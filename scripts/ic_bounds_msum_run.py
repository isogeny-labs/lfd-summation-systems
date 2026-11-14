#!/usr/bin/env python3
"""Driver for csd/ic_bounds_msum.py: last fall degree (LFD) of the DIRECT m-summand systems, m = 2..5.

Subcommands
  survey     build every planned instance (generators by evaluation, Z) and print degrees, |Z|, class (NO closure)
  instance   write instance files (prime: parameters + hashes [+ generator grid for GPU cells]; binary: ic-bin-instance/1)
  validate   the validation suite (generators against the existing code; engines against each other; negative controls)
  lfd        exact LFD with certificates for the named cells -> records/<case>.json (resumable)
  plan       GPU plan for the cells that do not fit (instances, chain script, cost) -> gpu/

Layout: experiments/ic-bounds-msummands/{README.md, predictions-2026-11-12.json, records/, validation/, gpu/, instances/}.
Memory: run every closure under ``ulimit -v 7000000`` and the heavy lock for anything above 1.2 GB resident.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import resource
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, ROOT)

from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_msum as S  # noqa: E402
from csd import ic_curve_sweep as CS  # noqa: E402
from csd.ic_arbitrary_bases import canonical_basis  # noqa: E402
from csd.ic_research_campaign import field  # noqa: E402

OUT = os.path.join(ROOT, 'experiments', 'ic-bounds-msummands')
SEED = 20261112


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def jhash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=str)
    os.replace(tmp, path)


# =============================================================================== cells

PRIME_PRIMES = [(97, True), (193, False), (257, True), (449, False), (641, True), (7681, False)]
PRIME_SEED_BASE = 10          # seeds < 10 are the pre-freeze smoke tests and the controls


def prime_cells():
    """(m, k) -> list of specs; 6 instances per cell (3 planted-H4/H2 targets, 3 random rational)."""
    out = []
    for m, k in [(3, 2), (3, 3), (3, 4), (4, 2), (4, 3), (5, 2)]:
        for j, (p, planted) in enumerate(PRIME_PRIMES):
            out.append(dict(family='prime', m=m, k=k, p=p, planted=planted, seed=PRIME_SEED_BASE + j,
                            name=f'P-m{m}k{k}-p{p}-{"pl" if planted else "r"}'))
    return out


BIN_CURVES = [('koblitz0', 0, 1), ('koblitz1', 1, 1), ('ordA', None, None), ('ordB', None, None)]
BIN_PLAN = [  # (m, n, k): local cells r = m k <= 16 whose estimated cost fits
    (3, 6, 2), (3, 8, 3), (3, 9, 4), (3, 8, 5),
    (4, 6, 2), (4, 7, 2), (4, 8, 3), (4, 9, 3),
    (5, 6, 2), (5, 7, 2), (5, 8, 3), (5, 9, 3),
    (4, 9, 4),
]


def binary_cells():
    out = []
    for m, n, k in BIN_PLAN:
        for j, (cname, a, b) in enumerate(BIN_CURVES):
            base = 'invariant' if (j % 2 == 1) else 'random'
            target = ['order-two', 'random', 'planted', 'random'][j]
            out.append(dict(family='binary', m=m, n=n, k=k, curve=cname, a=a, b=b, base=base, target=target,
                            seed=SEED, idx=j, name=f'B-m{m}n{n}k{k}-{cname}-{base[:3]}-{target[:3]}'))
    return out


def all_cells():
    return prime_cells() + binary_cells()


# =============================================================================== builders


def build_prime(spec):
    p, m, k = spec['p'], spec['m'], spec['k']
    inst = S.make_prime_instance(p, m, k, spec['seed'], spec['planted'])
    E = I.PrimeCurve(p, inst['a'], inst['b'])
    t0 = time.time()
    g, grid, g0 = S.prime_generator(E, inst['t'], m, k)
    pts = S.prime_zero_set(grid, m, k)
    # independence checks: Z from the grid equals the zero set of the normal form (forward DFT), planted decomposition in Z
    coef_back = S.forward_dft(S.inverse_dft(grid, p, g0, k), p, g0, k)
    assert np.array_equal(coef_back % p, grid % p)
    H = S.subgroup(p, k)[1]
    info = dict(spec=spec, instance=inst, g0=g0, generator_terms=len(g), generator_degree=S.generator_degree(g),
                generator_hash=S.generator_hash(p, k, m, g), solutions=len(pts),
                points_sha256=hashlib.sha256(json.dumps(pts).encode()).hexdigest(),
                build_seconds=round(time.time() - t0, 2))
    if inst['planted']:
        xs = [P[0] for P in inst['planted']]
        N = 1 << k
        if all(x in H for x in xs):
            pt = sum(H.index(x) << (i * k) for i, x in enumerate(xs))
            info['planted_in_Z'] = pt in set(pts)
            assert info['planted_in_Z'], 'planted decomposition not in Z'
        else:
            info['planted_in_Z'] = None   # planted points lie in H_2 which is inside H for k >= 2 (else not applicable)
    info['carry_prop1_variables'] = len(S.carry_hypotheses(g, m, k, p))
    info['instance_hash'] = jhash(dict(spec=spec, instance=inst, g0=g0))
    return dict(E=E, inst=inst, g=g, grid=grid, g0=g0, pts=pts, info=info)


def build_binary(spec):
    m, n, k = spec['m'], spec['n'], spec['k']
    F = field(n)
    rng = random.Random(int.from_bytes(hashlib.sha256(
        f"{spec['seed']}:msum-binary:{m}:{n}:{k}:{spec['curve']}:{spec['idx']}".encode()).digest()[:8], 'big'))
    a, b = spec['a'], spec['b']
    if a is None:
        a, b = rng.randrange(1 << n), rng.randrange(2, 1 << n)
    E = I.BinaryCurve(F, a, b)
    base_note = spec['base']
    if spec['base'] == 'invariant':
        cands = CS.invariant_candidates(F, k, 2, rng)
        if cands:
            basis = canonical_basis(cands[spec['idx'] % len(cands)])
        else:
            basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
            base_note = 'random(no invariant space of this dimension)'
    else:
        basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
    assert len(basis) == k
    V = I.span(basis)
    tk = spec['target']
    plant = None
    if tk == 'order-two':
        target = (0, F.sqrt(b))
    elif tk == 'planted':
        pl = S.planted_target_m(E, V, m, rng)
        if pl is None:
            target = CS.random_rational_target(E, rng, 1 << n, exclude_x={0})
            tk = 'random(planting failed)'
        else:
            target, plant = pl
    else:
        target = CS.random_rational_target(E, rng, 1 << n, exclude_x={0})
    assert E.on_curve(target)
    t0 = time.time()
    eqs, truth = S.binary_generator(E, basis, m, target[0])
    pts = S.binary_zero_set(truth)
    r = m * k
    degs = S.boolean_degrees([list(e) for e in eqs]) if eqs else [0]
    info = dict(spec=spec, n=n, modulus=F.modulus, a=a, b=b, basis=basis, target=list(target), target_class=tk,
                base_kind=base_note, frobenius_invariant=bool(I.is_frobenius_invariant(F, basis)),
                variables=r, equations=len(eqs), input_degree=max(degs), equation_degrees=degs,
                solutions=len(pts), points_sha256=hashlib.sha256(json.dumps(pts).encode()).hexdigest(),
                system_hash=S.binary_system_hash(f'binary-m{m}-direct-eval', r, eqs),
                const_equation=any(0 in e and len(e) == 1 for e in eqs), planted=plant,
                build_seconds=round(time.time() - t0, 2))
    if plant:
        xs = [P[0] for P in plant]
        pos = [V.index(x) for x in xs]
        pt = sum(c << (i * k) for i, c in enumerate(pos))
        info['planted_in_Z'] = pt in set(pts)
        assert info['planted_in_Z']
    info['instance_hash'] = jhash(dict(spec=spec, a=a, b=b, basis=basis, target=list(target), n=n))
    return dict(E=E, eqs=eqs, truth=truth, pts=pts, info=info, F=F, basis=basis, target=target)


def build_binary_light(spec):
    """as build_binary but only (d_in, |Z|) from the transformed truth table (no equation lists): for r = 19..20."""
    m, n, k = spec['m'], spec['n'], spec['k']
    F = field(n)
    rng = random.Random(int.from_bytes(hashlib.sha256(
        f"{spec['seed']}:msum-binary:{m}:{n}:{k}:{spec['curve']}:{spec['idx']}".encode()).digest()[:8], 'big'))
    a, b = spec['a'], spec['b']
    if a is None:
        a, b = rng.randrange(1 << n), rng.randrange(2, 1 << n)
    E = I.BinaryCurve(F, a, b)
    if spec['base'] == 'invariant':
        cands = CS.invariant_candidates(F, k, 2, rng)
        basis = canonical_basis(cands[spec['idx'] % len(cands)]) if cands else canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
    else:
        basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
    V = I.span(basis)
    if spec['target'] == 'order-two':
        target = (0, F.sqrt(b))
    elif spec['target'] == 'planted':
        pl = S.planted_target_m(E, V, m, rng)
        target = pl[0] if pl else CS.random_rational_target(E, rng, 1 << n, exclude_x={0})
    else:
        target = CS.random_rational_target(E, rng, 1 << n, exclude_x={0})
    d, z, terms = S.binary_degree_and_zeros(E, basis, m, target[0])
    return dict(d_in=d, Z=z, terms=terms, a=a, b=b, basis=basis, target=list(target))


def build(spec):
    return build_prime(spec) if spec['family'] == 'prime' else build_binary(spec)


# =============================================================================== survey / instance


def cmd_survey(a):
    rows = []
    for spec in all_cells():
        if a.family and spec['family'] != a.family:
            continue
        if a.m and spec['m'] not in a.m:
            continue
        c = build(spec)
        i = c['info']
        if spec['family'] == 'prime':
            row = dict(name=spec['name'], m=spec['m'], k=spec['k'], p=spec['p'], deg_g=i['generator_degree'],
                       terms=i['generator_terms'], Z=i['solutions'], prop1_vars=i['carry_prop1_variables'],
                       planted_in_Z=i.get('planted_in_Z'), cls='unit?' if i['solutions'] == 0 else 'proper',
                       gen_hash=i['generator_hash'][:16], secs=i['build_seconds'])
        else:
            row = dict(name=spec['name'], m=spec['m'], n=spec['n'], k=spec['k'], r=i['variables'],
                       d_in=i['input_degree'], eqs=i['equations'], Z=i['solutions'], const=i['const_equation'],
                       target=i['target_class'], base=i['base_kind'], inv=i['frobenius_invariant'],
                       cls='unit?' if i['solutions'] == 0 else 'proper', secs=i['build_seconds'])
        rows.append(row)
        log(json.dumps(row))
    if a.out:
        write_json(a.out, rows)


# =============================================================================== lfd


def record_prime(spec, c, engine, sandwich, cap, log_fn=None):
    i = c['info']
    p, m, k = spec['p'], spec['m'], spec['k']
    t0 = time.time()
    res = S.prime_ladder(p, m, k, c['g'], c['pts'], c['g0'], cap=cap, engine=engine, sandwich=sandwich, log=log_fn)
    rec = dict(case=spec['name'], family='prime', m=m, k=k, p=p, a=c['inst']['a'], b=c['inst']['b'], t=c['inst']['t'],
               target_class=c['inst']['target_class'], generator_degree=i['generator_degree'],
               generator_terms=i['generator_terms'], generator_hash=i['generator_hash'],
               instance_hash=i['instance_hash'], points_sha256=i['points_sha256'], solutions=i['solutions'],
               ideal_class='unit' if i['solutions'] == 0 else 'proper', engine=engine, sandwich=sandwich,
               status=res['status'], exact_lfd=res['exact_lfd'], lfd_lower_bound=res['lfd_lower_bound'],
               kernel_upper_bound=res['kernel_upper_bound'], strict_falls=res['strict_falls'],
               certificate=res['certificate'],
               certificate_degree=(res['certificate'] or {}).get('degree'),
               prop1_variables=i['carry_prop1_variables'], profile=res['profile'],
               cap_reason=res['cap_reason'], seconds=round(time.time() - t0, 2), max_rss_gb=round(rss_gb(), 2))
    return rec


def record_binary(spec, c, engine, cap, log_fn=None):
    i = c['info']
    t0 = time.time()
    r = i['variables']
    res = S.binary_ladder(c['eqs'], r, c['pts'], cap=cap, engine=engine, log=log_fn)
    rec = dict(case=spec['name'], family='binary', m=spec['m'], n=spec['n'], k=spec['k'], variables=r,
               curve=spec['curve'], a=i['a'], b=i['b'], basis=i['basis'], target=i['target'],
               target_class=i['target_class'], base_kind=i['base_kind'], frobenius_invariant=i['frobenius_invariant'],
               input_degree=i['input_degree'], equations=i['equations'], const_equation=i['const_equation'],
               system_hash=i['system_hash'], instance_hash=i['instance_hash'], points_sha256=i['points_sha256'],
               solutions=i['solutions'], ideal_class='unit' if i['solutions'] == 0 else 'proper',
               status=res['status'], exact_lfd=res['exact_lfd'], strict_falls=res['strict_falls'],
               certificate_degree=res.get('certificate_degree'), quotient_ranks=res['quotient_ranks'],
               local_checks=res['local_checks'], kernel_upper_bound=r + 1, cap_reason=res.get('cap_reason'),
               engine=res['engine'], seconds=round(time.time() - t0, 2), max_rss_gb=round(rss_gb(), 2))
    return rec


def est_binary_cost(r, d_in):
    """columns at the cap degree D = d_in + 2 and the 0.28 N^3 flop model (PLAN.md of ic-bounds-binary-deep)."""
    D = min(r, d_in + 2)
    N = sum(math.comb(r, j) for j in range(D + 1))
    return D, N, 0.28 * N ** 3


def cmd_lfd(a):
    cells = [s for s in all_cells() if (not a.cells or s['name'] in a.cells or any(s['name'].startswith(x) for x in a.prefix))]
    if a.family:
        cells = [s for s in cells if s['family'] == a.family]
    if a.m:
        cells = [s for s in cells if s['m'] in a.m]
    if a.k:
        cells = [s for s in cells if s['k'] in a.k]
    t_start = time.time()
    os.makedirs(os.path.join(OUT, 'records'), exist_ok=True)
    for spec in cells:
        path = os.path.join(OUT, 'records', spec['name'] + '.json')
        if os.path.exists(path) and not a.force:
            continue
        if time.time() - t_start > a.budget_min * 60:
            log('budget exhausted; rerun to continue')
            return
        c = build(spec)
        i = c['info']
        if spec['family'] == 'binary':
            D, N, fl = est_binary_cost(i['variables'], i['input_degree'])
            if N > a.max_columns:
                rec = dict(case=spec['name'], family='binary', m=spec['m'], n=spec['n'], k=spec['k'],
                           status='not-run: estimated columns %d at D=%d exceed --max-columns %d' % (N, D, a.max_columns),
                           input_degree=i['input_degree'], solutions=i['solutions'], system_hash=i['system_hash'],
                           estimated_columns=N, estimated_flops=fl)
                write_json(path, rec)
                log('skip', spec['name'], rec['status'])
                continue
            log('start', spec['name'], 'r', i['variables'], 'd_in', i['input_degree'], '|Z|', i['solutions'], 'N', N)
            rec = record_binary(spec, c, a.bin_engine or ('py' if i['variables'] <= 12 else 'big'), a.cap, None)
        else:
            log('start', spec['name'], 'deg g', i['generator_degree'], '|Z|', i['solutions'])
            rec = record_prime(spec, c, a.engine, a.sandwich, a.cap, None)
        write_json(path, rec)
        log('done', spec['name'], rec['status'], 'LFD', rec['exact_lfd'], 'falls', rec['strict_falls'], rec['seconds'], 's',
            'rss', rec['max_rss_gb'])


def cmd_validate(a):
    from csd import ic_bounds_msum_checks as C
    checks = dict(summation_prime=C.check_summation_prime, summation_binary=C.check_summation_binary,
                  prime_generators=C.check_prime_generators, binary_generators=C.check_binary_generators,
                  binary_m45=C.check_binary_m45_consistency, control_m2_prime=C.control_m2_prime,
                  control_m2_binary=C.control_m2_binary, control_m3_prime=C.control_m3_prime,
                  control_m3_binary=C.control_m3_binary, engines_prime=C.engines_prime, engines_binary=C.engines_binary,
                  negative_prime=C.negative_prime, negative_binary=C.negative_binary)
    names = a.which or list(checks)
    summary = {}
    for nm in names:
        t0 = time.time()
        log('check', nm)
        r = checks[nm]()
        r['seconds'] = round(time.time() - t0, 1)
        write_json(os.path.join(OUT, 'validation', r['name'] + '.json'), r)
        summary[nm] = dict(ok=r.get('ok'), seconds=r['seconds'])
        log(nm, 'ok=', r.get('ok'), r['seconds'], 's')
    print(json.dumps(summary, indent=1))


# =============================================================================== instance files


def driver_bin_hash(r, equations):
    """same definition as scripts/ic_bounds_bin_big_run.py::system_hash_bin (kind 'binary')."""
    return S.binary_system_hash('binary', r, equations)


def instance_binary(spec, c):
    i = c['info']
    r = i['variables']
    eqs = [list(e) for e in c['eqs']]
    F = c['F']
    inst = dict(format='ic-bin-instance/1', name=spec['name'], spec=spec, m=spec['m'], n=i['n'], modulus=i['modulus'], k=spec['k'],
                variables=r, basis=i['basis'], curve=dict(a=i['a'], b=i['b'], family='koblitz' if i['a'] in (0, 1) and i['b'] == 1 else 'ordinary'),
                target=i['target'], frobenius_invariant=i['frobenius_invariant'], equations=eqs,
                system_hash=driver_bin_hash(r, eqs), msum_system_hash=i['system_hash'], points=c['pts'], n_points=len(c['pts']),
                points_sha256=i['points_sha256'], input_degree=i['input_degree'],
                unit_class=bool(len(c['pts']) == 0), **{'class': 'unit' if not c['pts'] else 'proper'})
    return inst


def instance_prime(spec, c):
    i = c['info']
    return dict(format='ic-prime-msum-instance/1', name=spec['name'], spec=spec, m=spec['m'], k=spec['k'], p=spec['p'], g0=c['g0'],
                curve=dict(a=c['inst']['a'], b=c['inst']['b']), target=[c['inst']['t'], c['inst']['y']], planted=c['inst']['planted'],
                generator=sorted([int(ms), int(cf)] for ms, cf in c['g'].items()), generator_degree=i['generator_degree'],
                generator_hash=i['generator_hash'], points=c['pts'], points_sha256=i['points_sha256'], instance_hash=i['instance_hash'])


def cmd_instance(a):
    specs = json.load(open(a.spec))
    if isinstance(specs, dict):
        specs = [specs.get('spec', specs)]
    os.makedirs(a.outdir, exist_ok=True)
    for spec in specs:
        t0 = time.time()
        c = build(spec)
        inst = instance_binary(spec, c) if spec['family'] == 'binary' else instance_prime(spec, c)
        path = os.path.join(a.outdir, spec['name'] + '.json')
        with open(path, 'w') as f:
            json.dump(inst, f, separators=(',', ':'))
        log('instance', spec['name'], os.path.getsize(path), 'bytes', round(time.time() - t0, 1), 's')


def cmd_lfd_file(a):
    """prime instance file -> exact LFD ladder (CompactClosure, CPU)."""
    inst = json.load(open(a.instance))
    assert inst['format'] == 'ic-prime-msum-instance/1'
    g = {int(ms): int(cf) for ms, cf in inst['generator']}
    assert S.generator_hash(inst['p'], inst['k'], inst['m'], g) == inst['generator_hash'], 'generator hash mismatch'
    pts = [int(x) for x in inst['points']]
    assert hashlib.sha256(json.dumps(pts).encode()).hexdigest() == inst['points_sha256'], 'points hash mismatch'
    t0 = time.time()
    res = S.prime_ladder(inst['p'], inst['m'], inst['k'], g, pts, inst['g0'], cap=a.cap, engine='compact', sandwich=True, log=log)
    res.update(name=inst['name'], generator_hash=inst['generator_hash'], seconds=round(time.time() - t0, 1), max_rss_gb=round(rss_gb(), 2))
    if a.out:
        write_json(a.out, res)
    print(json.dumps({k: res[k] for k in ('name', 'status', 'exact_lfd', 'strict_falls', 'seconds')}), flush=True)


# =============================================================================== GPU plan

GPU_BIN = [  # (m, n, k, tag): binary cells with r = 16 that exceed the local column budget, and r = 18..24
    (4, 9, 4, 'r16'), (3, 10, 6, 'r18'), (4, 10, 5, 'r20'), (5, 10, 4, 'r20'), (3, 11, 7, 'r21'),
    (3, 12, 8, 'r24'), (4, 12, 6, 'r24'),
]
GPU_PRIME = [  # (m, k, D list, note)
    (4, 4, [12, 13], 'discriminating cell: d = 12 < m k = 16; (a) 17 vs (b) 13'),
    (5, 3, [15, 16], 'd = 15 = m k: LFD = 16 forced by Proposition 1'),
    (3, 6, [8], 'beyond the recorded D = 7 closure at k = 6 (V_7 fell behind, LFD >= 8)'),
    (3, 7, [7, 8], 'k = 7'),
]


def bin_gpu_specs():
    out = []
    for m, n, k, tag in GPU_BIN:
        for j, (cname, a, b) in enumerate(BIN_CURVES[:3]):
            base = 'invariant' if j == 1 else 'random'
            target = ['order-two', 'random', 'planted'][j]
            out.append(dict(family='binary', m=m, n=n, k=k, curve=cname, a=a, b=b, base=base, target=target, seed=SEED,
                            idx=j, tag=tag, name=f'G-m{m}n{n}k{k}-{cname}-{base[:3]}-{target[:3]}'))
    return out


def prime_gpu_specs():
    out = []
    for m, k, Ds, note in GPU_PRIME:
        for j, (p, planted) in enumerate(PRIME_PRIMES if (m, k) in ((4, 4), (5, 3)) else [(7681, False), (7681, True)]):
            if (p - 1) % (1 << k):
                continue
            out.append(dict(family='prime', m=m, k=k, p=p, planted=planted, seed=PRIME_SEED_BASE + j, Ds=Ds, note=note,
                            name=f'GP-m{m}k{k}-p{p}-{"pl" if planted else "r"}'))
    return out


def prime_cost(n, nvar, tflops64=20.0):
    madds = nvar * n ** 3 / 6.0
    return dict(columns=n, madds=madds, seconds_model=round(2 * madds / (tflops64 * 1e12), 1),
                seconds_pessimistic=round(8 * madds / (tflops64 * 1e12), 1),
                stored_peak_gb_float32=round(n * n * 4 / 4 / 2 ** 30, 2), stored_peak_gb_uint16=round(n * n * 2 / 4 / 2 ** 30, 2))


def cmd_plan(a):
    from csd import ic_bounds_bin_big as EB
    gdir = os.path.join(OUT, 'gpu')
    os.makedirs(os.path.join(gdir, 'specs'), exist_ok=True)
    cells, cost = [], []
    summary = dict(binary=[], prime=[])
    for spec in bin_gpu_specs():
        r = spec['m'] * spec['k']
        info = None
        deg_src = None
        if r <= a.build_max_r:
            c = build_binary(spec)
            info = c['info']
            d_in, zc = info['input_degree'], info['solutions']
            deg_src = 'built'
        elif r <= a.degree_max_r:
            c = build_binary_light(spec)
            d_in, zc = c['d_in'], c['Z']
            deg_src = 'evaluated (equations not built)'
        else:
            d_in = 6 if spec['m'] == 3 else 12     # theoretical upper bound on the Boolean degree (see README); NOT built
            zc = None
            deg_src = 'bound, not computed'
        Dcap = min(r, d_in + 2)
        pl = EB.plan(r, Dcap, 'gpu', batch=1024, work_gb=20.0, tile_mb=1024.0, dev_gb=141.0)
        row = dict(name=spec['name'], m=spec['m'], n=spec['n'], k=spec['k'], r=r, built=info is not None, input_degree=d_in,
                   input_degree_source=deg_src, solutions=zc, system_hash=(info or {}).get('system_hash'),
                   points_sha256=(info or {}).get('points_sha256'), instance_hash=(info or {}).get('instance_hash'),
                   cap_degree=Dcap, columns=pl['columns'], basis_worst_gb=pl['basis_worst_gb'], predicted_device_peak_gb=pl['predicted_peak_gb'],
                   fits_device_141gb=pl['fits_device'], seconds_model=pl['seconds_total_estimate'],
                   seconds_pessimistic=round(4 * pl['seconds_total_estimate']),
                   frozen_prediction_c=(dict(upper=min(r, d_in + 2), point=min(r, d_in + 1)) if info else None))
        summary['binary'].append(row)
        write_json(os.path.join(gdir, 'specs', spec['name'] + '.json'), dict(spec=spec, build=row))
        cells.append(dict(name=spec['name'], family='binary', spec=os.path.join('specs', spec['name'] + '.json'), max_degree=Dcap,
                          est_seconds=row['seconds_model'], time_limit_s=int(max(3600, 6 * row['seconds_model']))))
        log('binary', json.dumps({k: row[k] for k in ('name', 'r', 'input_degree', 'solutions', 'cap_degree', 'columns', 'predicted_device_peak_gb', 'seconds_model')}))
    for spec in prime_gpu_specs():
        m, k = spec['m'], spec['k']
        row = dict(name=spec['name'], m=m, k=k, p=spec['p'], Ds=spec['Ds'], note=spec['note'])
        if m * k <= 20 and (m, k) in ((4, 4), (5, 3)):
            c = build_prime(spec)
            i = c['info']
            row.update(generator_degree=i['generator_degree'], generator_terms=i['generator_terms'], generator_hash=i['generator_hash'],
                       solutions=i['solutions'], points_sha256=i['points_sha256'], instance_hash=i['instance_hash'],
                       prop1_variables=i['carry_prop1_variables'], built=True)
        else:
            row.update(built=False, solutions=None)
        costs = []
        for D in spec['Ds']:
            n = sum(math.comb(m * k, w) for w in range(min(D, m * k) + 1))
            costs.append(dict(D=D, **prime_cost(n, m * k)))
        row['cost_by_D'] = costs
        summary['prime'].append(row)
        write_json(os.path.join(gdir, 'specs', spec['name'] + '.json'), dict(spec=spec, build=row))
        cells.append(dict(name=spec['name'], family='prime', spec=os.path.join('specs', spec['name'] + '.json'), Ds=spec['Ds'],
                          est_seconds=max(c_['seconds_model'] for c_ in costs), time_limit_s=int(max(7200, 6 * max(c_['seconds_model'] for c_ in costs)))))
        log('prime', spec['name'], row.get('generator_degree'), row.get('solutions'), [c_['seconds_model'] for c_ in costs])
    write_json(os.path.join(gdir, 'cost.json'), summary)
    write_json(os.path.join(gdir, 'cells.json'), dict(order=[dict(cell='gpu-msum', purpose='cells not run locally', reserve=False, systems=cells)]))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    ins = sub.add_parser('instance')
    ins.add_argument('--spec', required=True)
    ins.add_argument('--outdir', required=True)
    lf = sub.add_parser('lfd-file')
    lf.add_argument('--instance', required=True)
    lf.add_argument('--cap', type=int, default=None)
    lf.add_argument('--out')
    pl = sub.add_parser('plan')
    pl.add_argument('--build-max-r', type=int, default=18)
    pl.add_argument('--degree-max-r', type=int, default=20)
    v = sub.add_parser('validate')
    v.add_argument('--which', nargs='*')
    s = sub.add_parser('survey')
    s.add_argument('--family')
    s.add_argument('--m', type=int, nargs='+')
    s.add_argument('--out')
    l = sub.add_parser('lfd')
    l.add_argument('--cells', nargs='*')
    l.add_argument('--prefix', nargs='*', default=[])
    l.add_argument('--family')
    l.add_argument('--m', type=int, nargs='+')
    l.add_argument('--k', type=int, nargs='+')
    l.add_argument('--engine', default='compact', choices=['compact', 'dense'])
    l.add_argument('--bin-engine', default=None, choices=['py', 'big'])
    l.add_argument('--sandwich', default='auto')
    l.add_argument('--cap', type=int, default=None)
    l.add_argument('--budget-min', type=float, default=20.0)
    l.add_argument('--max-columns', type=int, default=70000)
    l.add_argument('--force', action='store_true')
    a = ap.parse_args()
    if getattr(a, 'sandwich', None) in ('True', 'False'):
        a.sandwich = a.sandwich == 'True'
    {'survey': cmd_survey, 'lfd': cmd_lfd, 'validate': cmd_validate, 'plan': cmd_plan, 'instance': cmd_instance, 'lfd-file': cmd_lfd_file}[a.cmd](a)


if __name__ == '__main__':
    main()
