"""Local validation harness for csd/ic_bounds_lift.py.  Oracle: the DIRECT closure V_{D+1} of the engine
csd/ic_bounds_big.py (dimension, profile, contains_one, and the canonical hash of V_{D+1} n M_<=D per part).

  cases     --set NAME [--devices cpu,fake-gpu] [--modes plain,tiny,low]   direct vs lift, one row per run
  negative  negative controls (dropped shift, falls not shifted) on a list of cases
  claim     independent dense check of "low rows of the direct closure = reduced echelon basis of V n M_<=D"
  ckpt      lift from an engine checkpoint (device and host basis), lift checkpoint + resume after SIGTERM
  expected  values for the BUILTIN table of scripts/ic_bounds_lift_run.py

Run from the repository root with ``ulimit -v 1500000`` and one BLAS thread.
"""
import argparse
import gc
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'csd'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402
import ic_bounds_big as B  # noqa: E402
import ic_bounds_lift as LF  # noqa: E402

HARNESS_SHA256 = B.file_sha256(__file__)
VAL = ROOT / 'experiments' / 'ic-bounds-k14' / 'validation'

# instances: (p, a, b, X); the generator is the chain normal form of S_3(x1, x2, X) (engine chain_generator)
UNIT65537 = (65537, 1, 0, 4271)          # the unit instance of experiments/ic-bounds-k14/instances (no solutions)
UNIT12289 = (12289, 1, 0, 9362)
UNIT163841 = (163841, 3, 0, None)        # X read from the instance file
PROP65537_HI = (65537, 4021, 9917, 1424)  # proper-p65537-s12 (|Z| = 2 at k = 9, 14 at k = 10)
PROP65537_LO = (65537, 4021, 9917, 2)     # same curve, X = 2: |Z| = 2, 2, 4, 4 at k = 6, 7, 8, 9
PROP15361 = {6: (15361, 4021, 9917, 4491), 7: (15361, 4021, 9917, 7776)}

SETS = {
    'smoke': [('unit-12289', UNIT12289, k, D) for k, D in [(3, 2), (4, 2), (3, 3), (4, 3), (5, 3), (6, 3), (6, 4)]],
    'enter': [('unit-65537', UNIT65537, k, D) for k, D in [(3, 3), (4, 3), (5, 3), (6, 4), (7, 4), (8, 4), (9, 4)]],
    'stay': [('unit-65537', UNIT65537, k, D) for k, D in [(2, 2), (3, 2), (4, 2), (6, 3), (7, 3), (8, 3)]],
    'proper': [('proper-65537-X2', PROP65537_LO, k, D) for k, D in [(6, 3), (6, 4), (7, 3), (7, 4), (8, 3), (8, 4)]]
              + [('proper-65537-X1424', PROP65537_HI, 9, 4)]
              + [('proper-15361', PROP15361[6], 6, 4), ('proper-15361', PROP15361[7], 7, 4), ('proper-15361', PROP15361[7], 7, 3)],
    'tiny': [('unit-65537', UNIT65537, k, D) for k, D in [(4, 3), (5, 3), (6, 3), (6, 4), (7, 3), (7, 4)]]
            + [('proper-65537-X2', PROP65537_LO, 6, 4), ('proper-65537-X2', PROP65537_LO, 7, 4), ('proper-15361', PROP15361[7], 7, 4)],
    'bigp': [('p16777213-a1-X5', (16777213, 1, 0, 5), k, D) for k, D in [(5, 3), (6, 3), (6, 4)]],
    'k10d4': [('unit-65537', UNIT65537, 10, 4)],
    'k10d4-proper': [('proper-65537-X1424', PROP65537_HI, 10, 4)],
    'k11d4': [('unit-65537', UNIT65537, 11, 4)],
    'k9d5': [('unit-65537', UNIT65537, 9, 5)],
    'k10d5': [('unit-65537', UNIT65537, 10, 5)],
    'k10d5-proper': [('proper-65537-X1424', PROP65537_HI, 10, 5)],   # 14 solutions: V_6 must stay inside the ideal
}

MODES = {
    'plain': dict(cache_gb=0.25),
    'tiny': dict(tile_rows=7, batch_units=40, direct_rows=5, cache_gb=0.004, col_block=33, fall_every=2),
    'low': dict(include_low=True, fall_every=10 ** 9, cache_gb=0.25),
}


def gen(inst, k):
    p, a, b, X = inst
    return B.chain_generator(p, a, b, X, k)


def n_solutions(inst, k):
    p, a, b, X = inst
    try:
        return len(B.solutions(p, a, b, X, k))
    except AssertionError:
        return None


def lift_once(inst, k, D, device, mode, work_gb=0.4, **over):
    kw = dict(MODES[mode])
    kw.update(over)
    eng = LF.compute_vd(gen(inst, k), k, inst[0], D, device=device, work_gb=work_gb, tile_rows=11 if mode == 'tiny' else None)
    vd = dict(rank=eng.rank(), seconds=round(eng.seconds(), 1), madds=eng.total_madds())
    B._Live.peak = 0
    L = LF.Lift(eng, log=False, work_gb=work_gb, **kw)
    st = L.run()
    r = L.record(st)
    r['vd_closure'] = vd
    del L, eng
    gc.collect()
    return r


def atomic(path, obj):
    path = str(path)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1)
        f.write('\n')
    os.replace(tmp, path)


def cmd_cases(a):
    rows = []
    allok = True
    devices = a.devices.split(',')
    modes = a.modes.split(',')
    for label, inst, k, D in SETS[a.set]:
        t0 = time.time()
        if a.known_full:
            L1 = B.Layout(k, D + 1, True)
            LD = B.Layout(k, D, True)
            from math import comb
            direct = dict(dim=L1.c, profile=[sum(comb(2 * k, w) for w in range(j + 1)) for j in range(D + 2)], contains_one=True,
                          low_basis_sha256=[LF.full_space_sha(P['c']) for P in LD.parts], source=a.known_full)
        else:
            direct = LF.direct_record(gen(inst, k), k, inst[0], D + 1, device='cpu', work_gb=0.4)
        td = time.time() - t0
        for dev in devices:
            for mode in modes:
                r = lift_once(inst, k, D, dev, mode, **({'cache_gb': a.cache_gb} if a.cache_gb else {}))
                c = LF.compare(r, direct)
                ok = c['all'] and r['fixed_point'] and r['status'] == 'done'
                allok &= ok
                row = dict(case=f'{label}-k{k}-D{D}to{D + 1}', p=inst[0], a=inst[1], b=inst[2], X=inst[3], k=k, D=D, device=dev, mode=mode,
                           solutions=n_solutions(inst, k), ok=bool(ok), compare=c,
                           direct=dict(dim=direct['dim'], profile=direct['profile'], contains_one=direct['contains_one'],
                                       low_sha=[h[:16] for h in direct['low_basis_sha256']], seconds=round(td, 1),
                                       madds=direct.get('madds'), source=direct.get('source', 'engine closure, closed, full sweep clean')),
                           lift=dict(dim=r['dim_VD1'], profile=r['profile'], contains_one=r['contains_one'],
                                     low_sha=[h[:16] for h in r['low_basis_sha256']], dim_VD=r['dim_VD'], free=r['free_columns_VD'],
                                     regular_top=r['regular_top_orbits'], special_top=r['special_top_orbits'], dim_T=r['dim_T'],
                                     dim_K=r['dim_K'], falls=r['falls'], falls_new_per_round=r['falls_new_per_round'], rounds=r['rounds'],
                                     reductions=r['reductions'], counts=r['counts'], madds=r['madds'], seconds=r['seconds'],
                                     vd_closure=r['vd_closure'], dev_peak_live_mb=r['dev_peak_live_mb'], maxrss_mb=r['maxrss_mb'], settings=r['settings']))
                rows.append(row)
                print(json.dumps(dict(case=row['case'], device=dev, mode=mode, ok=row['ok'], dim=[r['dim_VD1'], direct['dim']],
                                      one=[r['contains_one'], direct['contains_one']], dim_K=r['dim_K'], dim_T=r['dim_T'],
                                      rounds=r['rounds'], lift_s=r['seconds'], direct_s=round(td, 1))), flush=True)
    out = dict(set=a.set, devices=devices, modes=modes, all_ok=bool(allok), rows=rows,
               source_sha256=dict(lift=LF.LIFT_SHA256, engine=B.ENGINE_SHA256, harness=HARNESS_SHA256))
    atomic(VAL / f'lift-{a.tag or a.set}.json', out)
    print('ALL OK' if allok else 'MISMATCH')
    sys.exit(0 if allok else 1)


def cmd_negative(a):
    rows = []
    cases = SETS['enter'][3:6] + SETS['stay'][3:] + [SETS['proper'][1], SETS['proper'][3]] + (SETS['k10d4'] if a.big else [])
    for label, inst, k, D in cases:
        direct = LF.direct_record(gen(inst, k), k, inst[0], D + 1, device='cpu', work_gb=0.4)
        ctrl = [('no_fall_shifts', dict(test_no_fall_shifts=True))] + [(f'drop_shift_{j}', dict(test_drop_shift=j)) for j in (0, k - 1)]
        for name, kw in ctrl:
            r = lift_once(inst, k, D, 'cpu', 'plain', **kw)
            c = LF.compare(r, direct)
            row = dict(case=f'{label}-k{k}-D{D}to{D + 1}', control=name, answer_changed=not c['all'], compare=c,
                       direct=dict(dim=direct['dim'], contains_one=direct['contains_one']),
                       lift=dict(dim=r['dim_VD1'], contains_one=r['contains_one'], status=r['status'], fixed_point=r['fixed_point']))
            rows.append(row)
            print(json.dumps(row), flush=True)
    out = dict(rows=rows, changed=sum(r['answer_changed'] for r in rows), total=len(rows),
               contains_one_changed=sum(not r['compare']['contains_one'] for r in rows),
               source_sha256=dict(lift=LF.LIFT_SHA256, engine=B.ENGINE_SHA256, harness=HARNESS_SHA256))
    atomic(VAL / 'lift-negative.json', out)
    sys.exit(0 if out['changed'] > 0 else 1)


def rref_modp(M, p):
    """plain reduced row echelon form over F_p (int64, small matrices)."""
    M = M.copy() % p
    r = 0
    piv = []
    for c in range(M.shape[1]):
        nz = np.flatnonzero(M[r:, c])
        if len(nz) == 0:
            continue
        i = r + nz[0]
        M[[r, i]] = M[[i, r]]
        M[r] = M[r] * pow(int(M[r, c]), -1, p) % p
        for j in np.flatnonzero(M[:, c]):
            if j != r:
                M[j] = (M[j] - M[j, c] * M[r]) % p
        piv.append(c)
        r += 1
        if r == M.shape[0]:
            break
    return M[:r], piv


def cmd_claim(a):
    """V = direct closure of degree D+1 (part q), as a dense matrix.  V n L computed independently as
    {x V : x in the left kernel of V[:, top columns]}, then its reduced echelon form on the low columns and the
    Part.sha-style hash; compared with LF.low_sha."""
    import hashlib
    rows = []
    for label, inst, k, D in [('unit-12289', UNIT12289, 5, 3), ('unit-65537', UNIT65537, 6, 3), ('proper-65537-X2', PROP65537_LO, 6, 4)]:
        p = inst[0]
        eng = B.Engine(gen(inst, k), k, p, D + 1, device='cpu', log=False, work_gb=0.3)
        assert eng.run() == 'done' and eng.closed
        for q, P in enumerate(eng.parts):
            sha, nlow = LF.low_sha(P, P.s0)
            V = P.rows_host(np.arange(P.rank)).astype(np.int64)          # scrambled on purpose below
            rng = np.random.default_rng(5 + q)
            S = rng.integers(0, p, size=(P.rank, P.rank))
            V = (S.astype(object) @ V.astype(object)) % p                 # another basis (S is invertible w.h.p.; rank checked)
            V = np.array(V, dtype=np.int64)
            s0 = P.s0
            # left kernel of the top block: rref of [Vtop | I]
            aug = np.concatenate([V[:, :s0], np.eye(P.rank, dtype=np.int64)], axis=1)
            Rr, piv = rref_modp(aug, p)
            kern = np.array([Rr[i, s0:] for i in range(len(piv)) if piv[i] >= s0], dtype=np.int64).reshape(-1, P.rank)
            inter = np.array((kern.astype(object) @ V.astype(object)) % p, dtype=np.int64)
            assert not inter[:, :s0].any()
            Rl, pl = rref_modp(inter[:, s0:], p)
            free = np.setdiff1d(np.arange(P.c - s0), pl)
            h = hashlib.sha256()
            h.update(np.array(pl, dtype=np.int64).tobytes())
            h.update(free.astype(np.int64).tobytes())
            h.update(np.ascontiguousarray(Rl[:, free].astype(np.uint16 if p < 65536 else np.uint32)).tobytes())
            row = dict(case=f'{label}-k{k}-D{D + 1}', part=q, rank=P.rank, dim_intersection_independent=len(pl), low_rows=nlow,
                       sha_independent=h.hexdigest()[:16], sha_low_rows=sha[:16], equal=h.hexdigest() == sha and len(pl) == nlow)
            rows.append(row)
            print(json.dumps(row), flush=True)
    ok = all(r['equal'] for r in rows)
    atomic(VAL / 'lift-claim-lowrows.json', dict(rows=rows, all_equal=ok, source_sha256=dict(lift=LF.LIFT_SHA256, engine=B.ENGINE_SHA256, harness=HARNESS_SHA256)))
    sys.exit(0 if ok else 1)


def cmd_ckpt(a):
    """engine checkpoint of V_D written by the engine driver -> lift driver with --ckpt (device and host basis),
    with a tight device budget, then SIGTERM + --resume of the lift itself."""
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    k, D = a.k, a.D
    inst = dict(k=k, prime=65537, a=1, b=0, seed=10, target=[4271, 11720], points=0,
                generator={f'{e1},{e2}': v for (e1, e2), v in sorted(gen(UNIT65537, k).items())})
    ifile = work / f'unit-k{k}.json'
    atomic(ifile, inst)
    ck = work / f'ck-k{k}-D{D}'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    run = lambda cmd, **kw: subprocess.run([sys.executable] + cmd, cwd=ROOT, env=env, capture_output=True, text=True, **kw)
    r = run(['scripts/ic_bounds_big_run.py', 'closure', '--instance', str(ifile), '--D', str(D), '--ckpt', str(ck), '--tile-rows', '64',
             '--work-gb', '0.4', '--out', str(work / 'closure.json')])
    assert r.returncode == 0, r.stderr[-2000:]
    direct = LF.direct_record(gen(UNIT65537, k), k, 65537, D + 1, device='cpu', work_gb=0.4)
    rows = []
    variants = [('cpu', []), ('fake-gpu-device-basis', ['--device', 'fake-gpu']),
                ('fake-gpu-host-basis-budget', ['--device', 'fake-gpu', '--basis-on', 'host', '--dev-gb', a.dev_gb, '--work-gb', '0.05', '--tile-rows', '64'])]
    for name, extra in variants:
        out = work / f'lift-{name}.json'
        r = run(['scripts/ic_bounds_lift_run.py', 'lift', '--ckpt', str(ck), '--instance', str(ifile), '--out', str(out)] + extra)
        rec = json.load(open(out)) if r.returncode == 0 else None
        c = LF.compare(rec, direct) if rec else None
        rows.append(dict(variant=name, exit=r.returncode, compare=c, fixed_point=rec and rec['fixed_point'],
                         evictions=rec and rec['counts']['follower_evictions'], dev_peak_live_mb=rec and rec['dev_peak_live_mb'],
                         vd=rec and {x: rec['vd'][x] for x in ('source', 'ckpt_reason', 'closure_complete', 'basis_on')}))
        print(json.dumps(rows[-1]), flush=True)
    # SIGTERM + resume of the lift (checkpoint every batch)
    import signal
    lck = work / 'lift-ck'
    cmd = ['scripts/ic_bounds_lift_run.py', 'lift', '--ckpt', str(ck), '--out', str(work / 'lift-resumed.json'), '--lift-ckpt', str(lck),
           '--ckpt-minutes', '0.0001', '--batch-units', '60', '--resume', '--log-seconds', '0']
    codes = []
    for delay in [float(x) for x in a.delays.split(',')]:
        pr = subprocess.Popen([sys.executable] + cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        time.sleep(delay)
        pr.send_signal(signal.SIGTERM)
        pr.communicate()
        codes.append(pr.returncode)
    r = run(cmd)
    rec = json.load(open(work / 'lift-resumed.json')) if r.returncode == 0 else None
    c = LF.compare(rec, direct) if rec else None
    rows.append(dict(variant='sigterm-resume', interrupted_exit_codes=codes, final_exit=r.returncode, compare=c,
                     resumes=rec and rec['resumes'], fixed_point=rec and rec['fixed_point']))
    print(json.dumps(rows[-1]), flush=True)
    ok = all(x['compare'] and x['compare']['all'] for x in rows)
    atomic(VAL / f'lift-ckpt-k{k}-D{D}.json', dict(k=k, D=D, direct=direct, rows=rows, all_ok=ok,
                                                source_sha256=dict(lift=LF.LIFT_SHA256, engine=B.ENGINE_SHA256, harness=HARNESS_SHA256)))
    sys.exit(0 if ok else 1)


def cmd_expected(a):
    specs = [('unit-12289-k5-D3', UNIT12289, 5, 3), ('unit-12289-k6-D3', UNIT12289, 6, 3), ('unit-12289-k6-D4', UNIT12289, 6, 4),
             ('unit-65537-k7-D4', UNIT65537, 7, 4), ('unit-65537-k7-D3', UNIT65537, 7, 3), ('proper-15361-k6-D4', PROP15361[6], 6, 4),
             ('proper-65537-k7-D4', PROP65537_LO, 7, 4), ('unit-65537-k8-D4', UNIT65537, 8, 4)]
    out = {}
    for name, inst, k, D in specs:
        d = LF.direct_record(gen(inst, k), k, inst[0], D + 1, device='cpu', work_gb=0.4)
        out[name] = dict(p=inst[0], a=inst[1], b=inst[2], X=inst[3], k=k, D=D, dim=d['dim'], profile=d['profile'], one=d['contains_one'],
                         sha=[h[:16] for h in d['low_basis_sha256']])
        print(f'    "{name}": {json.dumps(out[name])},', flush=True)
    atomic(VAL / 'lift-selftest-expected.json', dict(cases=out, source='direct closure, csd/ic_bounds_big.py', engine_sha256=B.ENGINE_SHA256,
                                                     lift_sha256=LF.LIFT_SHA256, harness_sha256=HARNESS_SHA256))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['cases', 'negative', 'claim', 'ckpt', 'expected'])
    ap.add_argument('--set', default='smoke', choices=sorted(SETS))
    ap.add_argument('--devices', default='cpu')
    ap.add_argument('--modes', default='plain')
    ap.add_argument('--tag')
    ap.add_argument('--known-full', help='do not run the direct closure: V_{D+1} is the full space by the named record')
    ap.add_argument('--big', action='store_true')
    ap.add_argument('--cache-gb', type=float, default=None)
    ap.add_argument('--work')
    ap.add_argument('--k', type=int, default=8)
    ap.add_argument('--D', type=int, default=4)
    ap.add_argument('--dev-gb', default='0.08')
    ap.add_argument('--delays', default='1.5,1.5')
    a = ap.parse_args()
    {'cases': cmd_cases, 'negative': cmd_negative, 'claim': cmd_claim, 'ckpt': cmd_ckpt, 'expected': cmd_expected}[a.cmd](a)


if __name__ == '__main__':
    main()
