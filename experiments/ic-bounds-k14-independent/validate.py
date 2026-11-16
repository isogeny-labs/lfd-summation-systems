#!/usr/bin/env python3
"""Validation of the independent lift (csd/ic_bounds_lift_independent.py) at small depth.

For each case: (1) the engine (scripts/ic_bounds_big_run.py, run as a subprocess) writes the stored
basis of V_D; (2) the independent lift decides 1 in V_{D+1} and dim(V_{D+1} n M_<=D) from that
checkpoint, in several settings; (3) ORACLE A is the direct closure in degree D + 1 by the reference
engine csd/ic_bounds_d5_verify.py (no lift); (4) ORACLE B (small cases) is a naive closure written
here in ordinary coordinates, from g alone; (5) the numbers recorded for the first lift in
experiments/ic-bounds-k14/LIFT.md are read from that file.  The first lift (csd/ic_bounds_lift.py)
is never imported or run.

    python3 experiments/ic-bounds-k14-independent/validate.py --sets quick [--device gpu] --out FILE
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'csd'))
import numpy as np                                   # noqa: E402
import ic_bounds_lift_independent as IL              # noqa: E402
import ic_bounds_d5_verify as REF                    # noqa: E402  (oracle A)

ENGINE = os.path.join(ROOT, 'scripts', 'ic_bounds_big_run.py')
LIFT_MD = os.path.join(ROOT, 'experiments', 'ic-bounds-k14', 'LIFT.md')

CURVES = {                       # name -> (p, a, b, X by depth or fixed)
    'unit-65537': (65537, 1, 0, 4271),
    'unit-12289': (12289, 1, 0, 9362),
    'unit-163841': (163841, 3, 0, 73730),
    'proper-65537-X2': (65537, 4021, 9917, 2),
    'proper-65537-X1424': (65537, 4021, 9917, 1424),
    'proper-15361': (15361, 4021, 9917, {6: 4491, 7: 7776}),
}
SETS = {
    'quick': [('unit-65537', 4, 3), ('unit-65537', 6, 3), ('unit-65537', 6, 4), ('proper-65537-X2', 6, 4)],
    'enter': [('unit-65537', k, 3) for k in (3, 4, 5)] + [('unit-65537', k, 4) for k in (6, 7, 8)]
             + [('unit-12289', 4, 3), ('unit-12289', 5, 3), ('unit-12289', 6, 4), ('unit-12289', 7, 4)]
             + [('unit-163841', 5, 3), ('unit-163841', 6, 4), ('unit-163841', 7, 4)],
    'stay': [('unit-65537', k, 2) for k in (2, 3, 4)] + [('unit-65537', k, 3) for k in (6, 7, 8)]
            + [('unit-12289', 3, 2), ('unit-12289', 6, 3), ('unit-12289', 7, 3)]
            + [('unit-163841', 6, 3), ('unit-163841', 7, 3)],
    'proper': [('proper-65537-X2', k, D) for k in (6, 7, 8) for D in (3, 4)]
              + [('proper-15361', 6, 4), ('proper-15361', 7, 4), ('proper-15361', 7, 3)],
    'k9': [('unit-65537', 9, 4), ('proper-65537-X1424', 9, 4)],
    'k10d5': [('unit-65537', 10, 5), ('proper-65537-X1424', 10, 5)],     # oracle: recorded direct closures
}
K14 = os.path.join(ROOT, 'experiments', 'ic-bounds-k14')


def recorded_oracle(name, k, D):
    """direct closures in degree D + 1 that are too large to recompute here: read from the records."""
    if (name, k, D) == ('unit-65537', 10, 5):
        f = os.path.join(K14, 'runs', 'reference-engine-unit-p65537-k10-D6.json')
        d = json.load(open(f))
        f2 = os.path.join(K14, 'runs', 'unit-p65537-k10-D6.json')
        e = json.load(open(f2))
        assert d['rank'] == e['rank'] and d['profile'] == e['profile'] and d['contains_one'] == e['contains_one']
        return dict(dim=d['rank'], contains_one=d['contains_one'], profile=d['profile'], dim_low=d['profile'][D],
                    dim_low_parts=[pp[D] for pp in d['part_profiles']] if 'part_profiles' in d else None,
                    source=[os.path.relpath(f, ROOT), os.path.relpath(f2, ROOT)],
                    note='recorded direct closures (reference engine on CPU, big engine on the H200); not recomputed here')
    if (name, k, D) == ('proper-65537-X1424', 10, 5):
        f = os.path.join(K14, 'validation', 'lift-k10d5-proper.json')
        d = json.load(open(f))['rows'][0]['direct']
        return dict(dim=d['dim'], contains_one=d['contains_one'], profile=d['profile'], dim_low=d['profile'][D],
                    dim_low_parts=None, source=[os.path.relpath(f, ROOT)],
                    note='recorded direct closure by the big engine on CPU (' + d['source'] + '); not recomputed here; no per-part profile in the record')
    return None
MODES = {
    'plain': dict(),
    'tiny': dict(batch=7, first_chunk=5, table_gb=2e-5, mi=5),        # many chunks, split products
    'low': dict(include_low=True),
    'seed9': dict(seed=9, first_chunk=64),
}
NAIVE_MAX_COLUMNS = 1700


def sha(path):
    return IL._sha256_file(path)


def params(name, k):
    p, a, b, X = CURVES[name]
    return p, a, b, (X[k] if isinstance(X, dict) else X)


def solutions(p, a, b, X, k):
    try:
        return len(REF.solutions(p, a, b, X, k)) if (p - 1) % (1 << k) == 0 else None
    except Exception:
        return None


def engine_basis(work, name, k, D):
    """stored basis of V_D written by the engine; returns (ckpt dir, engine record, instance file)."""
    p, a, b, X = params(name, k)
    g = IL.s3_generator(p, a, b, X, k)
    d = os.path.join(work, f'{name}-k{k}-D{D}')
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    inst = os.path.join(d, 'inst.json')
    json.dump(dict(k=k, prime=p, a=a, b=b, seed=0, target=[X, 0],
                   generator={f'{e[0]},{e[1]}': v for e, v in sorted(g.items())}), open(inst, 'w'))
    r = subprocess.run([sys.executable, ENGINE, 'closure', '--instance', inst, '--D', str(D), '--device', 'cpu',
                        '--threads', '2', '--no-eval', '--ckpt', os.path.join(d, 'ck'), '--out', os.path.join(d, 'rec.json')],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'engine failed ({r.returncode}): {r.stdout[-400:]} {r.stderr[-400:]}')
    return os.path.join(d, 'ck'), json.load(open(os.path.join(d, 'rec.json'))), inst


def oracle_reference(g, k, p, D1):
    t0 = time.time()
    C = REF.Closure(g, k, p, D1, split=True)
    C.run()
    tot, parts = C.profile()
    return dict(dim=C.rank(), contains_one=bool(C.contains_one()), profile=tot,
                dim_low=tot[D1 - 1], dim_low_parts=[pp[D1 - 1] for pp in parts],
                full_sweep_clean=bool(C.sweep_full()['passed']), seconds=round(time.time() - t0, 1))


def oracle_naive(g, k, p, D1):
    """least subspace of span(M_<=D1) containing g and closed under the shifts of its elements of
    weight < D1: ordinary coordinates, one echelon form, repeated sweeps until nothing is new."""
    T = IL.Table(*IL.monomials(k, range(D1, -1, -1)), k)
    E = IL.Echelon(T.n, p, np)
    v = np.zeros((1, T.n))
    for (e1, e2), c in g.items():
        v[0, T.find(np.array([e1]), np.array([e2]))[0]] = c
    E.add(v)
    low = np.flatnonzero(T.w < D1)
    maps = []
    for i in (0, 1):
        for j in range(k):
            b1 = (T.a1[low] + (1 << j)) % (1 << k) if i == 0 else T.a1[low]
            b2 = (T.a2[low] + (1 << j)) % (1 << k) if i == 1 else T.a2[low]
            maps.append(T.find(b1, b2))
    while True:
        rows = np.flatnonzero(T.w[E.piv] < D1)             # pivot of weight < D1 => supported on weight < D1
        M = np.zeros((len(rows), T.n))
        M[np.arange(len(rows)), E.piv[rows]] = 1.0
        M[:, E.free] = E.B[rows]
        assert not M[:, T.w == D1].any()
        new = 0
        for m in maps:
            S = np.zeros_like(M)
            S[:, m] = M[:, low]
            new += E.add(S)
        if new == 0:
            break
    one = int(T.find(np.array([0]), np.array([0]))[0])
    return dict(dim=E.rank, contains_one=bool(E.rowof[one] >= 0), dim_low=int(np.sum(T.w[E.piv] < D1)))


def lift_md_table():
    """rows of the validation table of LIFT.md (first lift), keyed by (instance, k, D)."""
    out = {}
    pat = re.compile(r'^\| \w[\w-]*: ([\w-]+)-k(\d+)-D(\d+)to(\d+) \(p = (\d+)\) \| (\w+) \| (\d+) \| (\d+), (yes|no) \| (\d+), (yes|no) \| '
                     r'yes \(([\d, ]+)\) \| [^|]* \| \[(\d+), (\d+)\], \[(\d+), (\d+)\] \| (\d+) \|')
    for line in open(LIFT_MD):
        m = pat.match(line)
        if m:
            prof = [int(x) for x in m.group(12).split(',')]
            D = int(m.group(3))
            out.setdefault((m.group(1), int(m.group(2)), D), dict(
                p=int(m.group(5)), dim_VD=int(m.group(7)), dim_VD1=int(m.group(10)), contains_one=m.group(11) == 'yes',
                dim_low=prof[D], dim_T=int(m.group(13)) + int(m.group(14)), dim_K=int(m.group(15)) + int(m.group(16))))
    return out


def run_lift(ck, k, D, p, g, xp, hook=None, **kw):
    lift = IL.Lift(k, D, p, ck['parts'], g, xp=xp, **kw)
    if hook:
        hook(lift)
    return lift.run()


KEYS = ('contains_one', 'dim_low', 'dim_low_parts', 'dim_VD1')


def do_case(work, name, k, D, xp, modes, md):
    p, a, b, X = params(name, k)
    g = IL.s3_generator(p, a, b, X, k)
    assert g == REF.chain_generator(p, a, b, X, k), 'generator differs from the reference module'
    ckdir, erec, inst = engine_basis(work, name, k, D)
    ck = IL.read_checkpoint(ckdir)
    case = dict(case=f'{name}-k{k}-D{D}to{D + 1}', instance=name, p=p, a=a, b=b, X=X, k=k, D=D,
                solutions=solutions(p, a, b, X, k),
                engine_VD=dict(dim=erec['rank'], closed=erec['closed'], contains_one=erec['contains_one'],
                               basis_sha256=erec['basis_sha256']))
    ora = recorded_oracle(name, k, D)
    case['oracle_kind'] = 'recorded' if ora else 'recomputed'
    if ora is None:
        ora = oracle_reference(g, k, p, D + 1)
    case['oracle_reference_direct'] = ora
    want = dict(contains_one=ora['contains_one'], dim_low=ora['dim_low'], dim_low_parts=ora['dim_low_parts'], dim_VD1=ora['dim'])
    runs, ok = {}, True
    if erec['columns'] > 800:
        modes = [m for m in modes if m != 'tiny'] or modes
    for mode in modes:
        r = run_lift(ck, k, D, p, g, xp, **MODES[mode])
        if want['dim_low_parts'] is None:                     # the record has no per-part profile
            want['dim_low_parts'] = r['dim_low_parts'] if sum(r['dim_low_parts']) == want['dim_low'] else None
        if mode == modes[0]:                                  # hashes of what was read = the engine's record
            L0 = IL.Lift(k, D, p, ck['parts'], g, xp=np)
            case['canonical_sha_equal_engine'] = [IL.canonical_sha(pv, B, nc, p) for (pv, B), nc in zip(ck['parts'], L0.ncol)] == erec['basis_sha256']
            # two routes to phi agree on a sample (gather of two rows against products with whole parts)
            rng = np.random.RandomState(3)
            pdm = np.flatnonzero(L0.Pd)
            smp = pdm[rng.permutation(len(pdm))[:40]]
            U = np.zeros((len(smp), L0.L.n))
            U[np.arange(len(smp)), smp] = 1.0
            case['phi_two_routes_equal'] = bool(len(smp) == 0 or np.array_equal(IL.neg(L0.red_L(U), p), L0.phi(smp)))
            ok = ok and case['canonical_sha_equal_engine'] and case['phi_two_routes_equal']
        agree = all(r[key] == want[key] for key in KEYS) and r['fixed_point'] and r['sigma_stable'] and r['dim_VD'] == erec['rank']
        ok = ok and agree
        runs[mode] = dict({key: r[key] for key in KEYS + ('fixed_point', 'saturated', 'schedule_complete', 'verify_new_rows',
                                                          'sigma_stable', 'dim_VD', 'n_T', 'rank_K', 'n_single', 'falls_new',
                                                          'n_special_tops', 'n_regular_tops', 'closed_check', 'seconds')},
                          chunks=r['counts']['chunks'], fall_rounds=r['counts']['fall_rounds'], agrees_with_oracle=bool(agree))
    case['runs'] = runs
    r0 = runs[modes[0]]
    cols = r0['n_regular_tops'] + r0['n_special_tops'] + erec['columns']
    if cols <= NAIVE_MAX_COLUMNS:
        nv = oracle_naive(g, k, p, D + 1)
        nv['agrees'] = bool(nv['dim'] == r0['dim_VD1'] and nv['contains_one'] == r0['contains_one'] and nv['dim_low'] == r0['dim_low'])
        case['oracle_naive_ordinary'] = nv
        ok = ok and nv['agrees']
    key = (name, k, D)
    if key in md:
        m = md[key]
        cmp_ = dict(dim_VD=m['dim_VD'] == r0['dim_VD'], dim_VD1=m['dim_VD1'] == r0['dim_VD1'],
                    contains_one=m['contains_one'] == r0['contains_one'], dim_low=m['dim_low'] == r0['dim_low'],
                    dim_T=m['dim_T'] == r0['n_T'] - r0['n_single'], dim_K=m['dim_K'] == r0['rank_K'] - r0['n_single'])
        case['first_lift_recorded'] = dict(m, equal=cmp_, all_equal=all(cmp_.values()))
        ok = ok and all(cmp_.values())
    case['ok'] = bool(ok)
    return case, ck, erec, g, ckdir, inst


def negative_controls(work, name, k, D, xp):
    """each control must end in a refusal or in a result different from the clean one."""
    p, a, b, X = params(name, k)
    g = IL.s3_generator(p, a, b, X, k)
    ckdir, erec, inst = engine_basis(work, name, k, D)
    ck = IL.read_checkpoint(ckdir)
    clean = run_lift(ck, k, D, p, g, xp)
    base = {key: clean[key] for key in KEYS}
    out = dict(case=f'{name}-k{k}-D{D}to{D + 1}', clean=base, controls=[])

    def attempt(label, fn):
        try:
            r = fn()
            res = {key: r[key] for key in KEYS}
            c = dict(control=label, outcome='result', result=res, changed=res != base)
        except IL.Refuse as e:
            c = dict(control=label, outcome='refused', reason=str(e))
        c['detected'] = c['outcome'] == 'refused' or c['changed']
        out['controls'].append(c)

    # 1. generator: one coefficient changed
    for e in ((0, 0), (1, 1)):
        g2 = dict(g)
        g2[e] = (g2.get(e, 0) + 1) % p
        attempt(f'generator coefficient {e} + 1', lambda g2=g2: run_lift(ck, k, D, p, g2, xp))
    # 2. one entry of the stored basis changed in memory (file hashes bypassed)
    Lay = IL.Lift(k, D, p, ck['parts'], g, xp=np)
    rng = np.random.RandomState(11)
    for q in (0, 1):
        piv, B = ck['parts'][q]
        wq = Lay.L.w[np.flatnonzero(Lay.L.a1 <= Lay.L.a2 if q == 0 else Lay.L.a1 < Lay.L.a2)]
        free = np.setdiff1d(np.arange(len(wq)), piv)
        for kind, rows in (('low', np.flatnonzero(wq[piv] < D)), ('top', np.flatnonzero(wq[piv] == D))):
            for trial in range(3):
                if len(rows) == 0 or B.shape[1] == 0:
                    continue
                r_ = int(rows[rng.randint(len(rows))])
                cands = np.flatnonzero(wq[free] <= wq[piv[r_]])       # keep the row's weight bound
                c_ = int(cands[rng.randint(len(cands))])
                parts2 = [(pv, Bq.copy()) for pv, Bq in ck['parts']]
                parts2[q][1][r_, c_] = (parts2[q][1][r_, c_] + 1) % p
                ck2 = dict(ck, parts=parts2)
                attempt(f'basis part {q}, {kind} row {r_}, free column {c_} (weight {int(wq[free][c_])}) + 1',
                        lambda ck2=ck2: run_lift(ck2, k, D, p, g, xp))
                c = out['controls'][-1]
                L2 = IL.Lift(k, D, p, parts2, g, xp=np)
                c['canonical_sha_differs_from_engine_record'] = \
                    [IL.canonical_sha(pv, Bq, nc, p) for (pv, Bq), nc in zip(parts2, L2.ncol)] != erec['basis_sha256']
    # 3. one byte of a tile file changed on disk: the reader must refuse
    d2 = ckdir + '-corrupt'
    shutil.rmtree(d2, ignore_errors=True)
    shutil.copytree(ckdir, d2)
    slot = sorted(x for x in os.listdir(d2) if x.startswith('ck-'))[-1]
    tile = os.path.join(d2, slot, 'p0_t00000.npy')
    raw = bytearray(open(tile, 'rb').read())
    raw[-3] ^= 0x40
    open(tile, 'wb').write(bytes(raw))
    attempt('one byte of p0_t00000.npy flipped on disk', lambda: run_lift(IL.read_checkpoint(d2), k, D, p, g, xp))
    # 4. algorithm: the fall rows are never shifted (hook in this harness only)
    def no_falls(lift):
        lift.falls = lambda verify=False: 0
    attempt('TEST HOOK: fall rows never shifted', lambda: run_lift(ck, k, D, p, g, xp, hook=no_falls))
    # 5. algorithm: one shift index left out of the generator pass (hook)
    def drop_shift(lift):
        orig_chunk = lift.chunk
        lift.chunk = lambda tops, direct: orig_chunk(tops, direct[lift.dir_u[direct] != 0])
    attempt('TEST HOOK: direct generators of shift u_{1,0} dropped', lambda: run_lift(ck, k, D, p, g, xp, hook=drop_shift))
    # 6. algorithm: no direct generator at all (hook)
    def drop_direct(lift):
        orig_chunk = lift.chunk
        lift.chunk = lambda tops, direct: orig_chunk(tops, direct[:0])
    attempt('TEST HOOK: all direct generators dropped', lambda: run_lift(ck, k, D, p, g, xp, hook=drop_direct))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sets', default='quick')
    ap.add_argument('--modes', default='plain,tiny,low,seed9')
    ap.add_argument('--device', default='cpu', choices=['cpu', 'gpu'])
    ap.add_argument('--negative', default='unit-65537:7:4,unit-65537:7:3,proper-65537-X2:7:4,unit-163841:6:4')
    ap.add_argument('--work', default=os.path.join(tempfile.gettempdir(), 'lfd-lift-independent'))
    ap.add_argument('--out')
    a = ap.parse_args()
    xp = np
    if a.device == 'gpu':
        import cupy as xp
    os.makedirs(a.work, exist_ok=True)
    md = lift_md_table()
    t0 = time.time()
    rec = dict(harness='experiments/ic-bounds-k14-independent/validate.py', argv=sys.argv[1:], device=a.device,
               numpy=np.__version__, lift_md_rows_parsed=len(md),
               source_sha256=dict(lift=sha(IL.__file__), harness=sha(__file__), oracle_reference=sha(REF.__file__),
                                  engine=sha(os.path.join(ROOT, 'csd', 'ic_bounds_big.py')), engine_driver=sha(ENGINE),
                                  lift_md=sha(LIFT_MD)),
               cases=[], negative=[])
    for s in a.sets.split(','):
        for name, k, D in SETS[s]:
            case, *_ = do_case(a.work, name, k, D, xp, a.modes.split(','), md)
            case['set'] = s
            rec['cases'].append(case)
            r0 = case['runs'][a.modes.split(',')[0]]
            print(json.dumps(dict(set=s, case=case['case'], ok=case['ok'], one=r0['contains_one'], dim_low=r0['dim_low'],
                                  parts=r0['dim_low_parts'], dim=r0['dim_VD1'], saturated=r0['saturated'],
                                  oracle_s=case['oracle_reference_direct'].get('seconds'),
                                  lift_s=[case['runs'][m]['seconds'] for m in case['runs']],
                                  first_lift=case.get('first_lift_recorded', {}).get('all_equal'),
                                  naive=case.get('oracle_naive_ordinary', {}).get('agrees'))), flush=True)
    for spec in [x for x in a.negative.split(',') if x]:
        name, k, D = spec.split(':')
        neg = negative_controls(a.work, name, int(k), int(D), xp)
        rec['negative'].append(neg)
        for c in neg['controls']:
            print(json.dumps(dict(negative=neg['case'], control=c['control'], outcome=c['outcome'],
                                  detected=c['detected'], reason=c.get('reason'), result=c.get('result'))), flush=True)
    cs = rec['cases']
    rec['totals'] = dict(
        cases=len(cs), lift_runs=sum(len(c['runs']) for c in cs), all_ok=all(c['ok'] for c in cs),
        mismatches_with_reference_oracle=sum(not r['agrees_with_oracle'] for c in cs for r in c['runs'].values()),
        one_enters=sum(c['oracle_reference_direct']['contains_one'] for c in cs),
        one_stays_out=sum(not c['oracle_reference_direct']['contains_one'] for c in cs),
        naive_oracle_cases=sum('oracle_naive_ordinary' in c for c in cs),
        naive_oracle_agree=sum(c.get('oracle_naive_ordinary', {}).get('agrees', False) for c in cs),
        first_lift_rows_compared=sum('first_lift_recorded' in c for c in cs),
        first_lift_rows_equal=sum(c.get('first_lift_recorded', {}).get('all_equal', False) for c in cs),
        negative_controls=sum(len(n['controls']) for n in rec['negative']),
        negative_detected=sum(c['detected'] for n in rec['negative'] for c in n['controls']),
        seconds=round(time.time() - t0, 1))
    print(json.dumps(rec['totals']), flush=True)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(rec, f, indent=1)
            f.write('\n')
    sys.exit(0 if rec['totals']['all_ok'] else 1)


if __name__ == '__main__':
    main()
