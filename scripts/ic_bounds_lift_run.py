"""Driver for csd/ic_bounds_lift.py (degree lift through the quotient: V_D -> V_{D+1}).

Subcommands
  lift      from a finished degree-D closure (``--ckpt DIR``: an engine checkpoint written by
            ic_bounds_big_run.py closure --ckpt DIR; or ``--instance FILE --D D``: V_D is computed first with the
            engine, small cases) decide 1 in V_{D+1}, compute dim V_{D+1}, its profile and the canonical basis
            hash of V_{D+1} n M_<=D per part
  plan      predicted sizes, memory and number of images for k, D and the free-column counts of V_D
  selftest  arithmetic exactness on the device + built-in small cases on cpu and on the device against
            embedded values from the direct closure

The two engine files (ic_bounds_big.py, ic_bounds_lift.py) are found next to this script (remote layout: all in
/workspace/code) or in ../csd/ (repository layout).  No repository imports, no network.

Exit codes: 0 fixed point reached ; 1 selftest/other failure ; 2 stopped because 1 is in V_{D+1} (--stop-on-one;
the record's dimension, profile and hashes are then lower bounds only) ; 4 a test hook was used (not a result) ;
75 SIGTERM (resume with --resume when --lift-ckpt was given).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
for cand in (HERE, HERE.parent / 'csd'):
    if (cand / 'ic_bounds_lift.py').exists() and (cand / 'ic_bounds_big.py').exists():
        sys.path.insert(0, str(cand))
        break
else:  # pragma: no cover
    raise SystemExit('ic_bounds_lift.py / ic_bounds_big.py not found next to this script or in ../csd/')


def _set_blas_threads(argv):
    n = None
    for i, x in enumerate(argv):
        if x == '--threads' and i + 1 < len(argv):
            n = int(argv[i + 1])
        elif x.startswith('--threads='):
            n = int(x.split('=', 1)[1])
    keys = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')
    if n is not None:
        for key in keys:
            os.environ[key] = str(n)
    else:
        try:
            d = len(os.sched_getaffinity(0))
        except Exception:
            d = 1
        try:
            q, per = open('/sys/fs/cgroup/cpu.max').read().split()[:2]
            if q != 'max':
                d = max(1, min(d, -(-int(q) // int(per))))
        except Exception:
            pass
        for key in keys:
            os.environ.setdefault(key, str(d))


_set_blas_threads(sys.argv)
import numpy as np  # noqa: E402
import ic_bounds_big as B  # noqa: E402
import ic_bounds_lift as LF  # noqa: E402

DRIVER_SHA256 = B.file_sha256(__file__)
T0 = time.time()

# built-in selftest cases: values from the DIRECT closure V_{D+1} of the engine (csd/ic_bounds_big.py), see
# experiments/ic-bounds-k14/validation/lift-selftest-expected.json; sha = first 16 hex digits per part
BUILTIN = {
    "unit-12289-k5-D3": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 5, "D": 3, "dim": 386, "profile": [1, 11, 56, 176, 386], "one": True, "sha": ["51994086169266ae", "9a7f9304eba34229"]},
    "unit-12289-k6-D3": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 6, "D": 3, "dim": 605, "profile": [0, 0, 2, 111, 605], "one": False, "sha": ["e92935475f7482e4", "dea12bcefc3502d3"]},
    "unit-12289-k6-D4": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 6, "D": 4, "dim": 1586, "profile": [1, 13, 79, 299, 794, 1586], "one": True, "sha": ["c13c0eb42e6d9979", "3583762bc711f4a7"]},
    "unit-65537-k7-D4": {"p": 65537, "a": 1, "b": 0, "X": 4271, "k": 7, "D": 4, "dim": 3473, "profile": [1, 15, 106, 470, 1471, 3473], "one": True, "sha": ["46eee6288dc831b9", "32a2dadeed2e7d8a"]},
    "unit-65537-k7-D3": {"p": 65537, "a": 1, "b": 0, "X": 4271, "k": 7, "D": 3, "dim": 829, "profile": [0, 0, 2, 115, 829], "one": False, "sha": ["9d59679046053325", "8186e3b7f0a7169f"]},
    "proper-15361-k6-D4": {"p": 15361, "a": 4021, "b": 9917, "X": 4491, "k": 6, "D": 4, "dim": 1584, "profile": [0, 11, 77, 297, 792, 1584], "one": False, "sha": ["559de8b301aa48db", "6949095244de85fd"]},
    "proper-65537-k7-D4": {"p": 65537, "a": 4021, "b": 9917, "X": 2, "k": 7, "D": 4, "dim": 3471, "profile": [0, 13, 104, 468, 1469, 3471], "one": False, "sha": ["4aa5406301f83103", "5e4462fcf36b7fe9"]},
}


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['maxrss_mb'] = B.maxrss_mb()
    print(json.dumps(kw), flush=True)


def atomic_write_json(path, obj):
    path = str(path)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def parse_gen(d):
    return {tuple(int(x) for x in key.split(',')): int(v) for key, v in d.items()}


def lift_kwargs(a):
    return dict(work_gb=a.work_gb, cache_gb=a.cache_gb, tile_mb=a.tile_mb, tile_rows=a.tile_rows, col_block=a.col_block,
                batch_units=a.batch_units, direct_rows=a.direct_rows, stop_on_one=a.stop_on_one and not a.to_fixed_point,
                include_low=a.include_low, fall_every=a.fall_every, status=a.status, log_seconds=a.log_seconds, ckpt=a.lift_ckpt,
                ckpt_minutes=a.ckpt_minutes, driver_sha=DRIVER_SHA256, test_drop_shift=a.test_drop_shift,
                test_no_fall_shifts=a.test_no_fall_shifts)


def cmd_lift(a):
    inst = json.load(open(a.instance)) if a.instance else None
    inst_sha = B.file_sha256(a.instance) if a.instance else None
    if a.ckpt:
        log(event='loading_checkpoint', ckpt=a.ckpt, basis_on=a.basis_on)
        eng, info = LF.load_vd(a.ckpt, device=a.device, basis_on=a.basis_on, work_gb=a.work_gb, force=a.force_ckpt)
        if inst and (inst['k'], inst['prime']) != (eng.k, eng.p):
            raise SystemExit('instance file does not match the checkpoint (k, p)')
    else:
        if not inst or a.D is None:
            raise SystemExit('lift needs --ckpt DIR, or --instance FILE --D D')
        g = parse_gen(inst['generator'])
        log(event='computing_VD', k=inst['k'], D=a.D, p=inst['prime'])
        eng = LF.compute_vd(g, inst['k'], inst['prime'], a.D, device=a.device, work_gb=a.work_gb, tile_mb=a.tile_mb)
        info = dict(source='engine closure in this process', closed=bool(eng.closed), rank=eng.rank(),
                    basis_sha256=[P.sha() for P in eng.parts], seconds=round(eng.seconds(), 1),
                    basis_id=','.join(P.sha() for P in eng.parts))
    L = LF.Lift(eng, vd_info=info, **lift_kwargs(a))
    if a.dev_gb:
        gb = 1 << 30
        basis_dev = sum(P.dev_bytes() for P in eng.parts) if a.device != 'cpu' else 0
        kworst = max(L.Tdim) ** 2
        slack = int(1.5 * gb) if a.device == 'gpu' else 0
        cache = int(a.dev_gb * gb) - L.work_bytes - kworst - slack - basis_dev
        if a.cache_gb is None:
            if cache <= 0:
                raise SystemExit(f'--dev-gb too small: nothing left for the follower cache ({cache / gb:.2f} GB); use --basis-on host')
            L.cache_bytes = cache
        if a.device == 'fake-gpu':
            B._Live.limit = int(a.dev_gb * gb)
    resumed = False
    if a.resume:
        if not a.lift_ckpt:
            raise SystemExit('--resume needs --lift-ckpt DIR')
        resumed = L.load_checkpoint()
        if not resumed:
            log(event='no_lift_checkpoint_found_starting_fresh')
    L.install_sigterm()
    log(event='start', k=L.k, D=L.D, p=L.p, device=a.device, dim_VD=sum(L.rankD), free=L.nfr, regular_weight_D_orbits=L.nO,
        regular_top=L.nReg, special_top=L.nSp, dim_T=L.Tdim, units=len(L.U_o), direct_rows=[len(x) for x in L.direct],
        batches=[len(L.plan_batches(q)) for q in (0, 1)], schedule_items=len(L.schedule()), cache_gb=None if L.cache_bytes is None else round(L.cache_bytes / (1 << 30), 2),
        resumed=resumed, lift_sha256=LF.LIFT_SHA256, engine_sha256=B.ENGINE_SHA256, driver_sha256=DRIVER_SHA256)
    status = L.run()
    if status == 'sigterm':
        log(event='exit_sigterm')
        sys.exit(LF.EXIT_SIGTERM)
    B._Live.limit = None
    rec = L.record(status)
    rec['resumed_from_checkpoint'] = resumed
    if inst:
        rec.update(instance_file=os.path.basename(a.instance), instance_sha256=inst_sha, a=inst.get('a'), b=inst.get('b'),
                   seed=inst.get('seed'), want=inst.get('want'), unit_nonresidue=inst.get('unit_nonresidue'),
                   instance_target=inst.get('target'), instance_points=inst.get('points'))
    log(**{key: v for key, v in rec.items() if key not in ('vd', 'k_madds_detail', 'part_profiles_low')}, event='result')
    if a.out:
        atomic_write_json(a.out, rec)
        log(event='saved', file=a.out)
    hook = a.test_drop_shift is not None or a.test_no_fall_shifts
    sys.exit(4 if hook else {'done': 0, 'one': 2}.get(status, 1))


def cmd_plan(a):
    free = [int(x) for x in a.free.split(',')]
    ftop = [int(x) for x in a.free_top.split(',')]
    print(json.dumps(LF.plan(a.k, a.D, free, ftop, a.work_gb, dev_gb=a.dev_gb, basis_on=a.basis_on), indent=1))


def run_case(spec, device, **kw):
    g = B.chain_generator(spec['p'], spec['a'], spec['b'], spec['X'], spec['k'])
    eng = LF.compute_vd(g, spec['k'], spec['p'], spec['D'], device=device, work_gb=0.5, batch=64, tile_rows=kw.pop('eng_tile_rows', None))
    L = LF.Lift(eng, log=False, **kw)
    st = L.run()
    return L.record(st)


def cmd_selftest(a):
    dev = a.device
    ok = True
    ar = B.selftest_arith(dev)
    log(event='arith_device', **ar)
    ok &= ar['passed']
    # addmul with the inner-dimension split active (p near 2^24: max_inner = 4), against Python integers
    B.use_device(dev)
    rng = np.random.default_rng(3)
    am_fail = 0
    for p in (16777213, 65537, 163841):
        mi = B.max_inner(p)
        A = rng.integers(0, p, size=(3, 41))
        M = rng.integers(0, p, size=(41, 5))
        acc = rng.integers(0, p, size=(3, 5))
        for sign in (1, -1):
            want = (acc.astype(object) + sign * (A.astype(object) @ M.astype(object))) % p
            got = LF.addmul(B.DV.to_dev(acc.astype(np.float64)), B.DV.to_dev(A.astype(np.float64)), B.DV.to_dev(M.astype(np.float64)), p, min(mi, 7), sign)
            am_fail += int(np.any(B.DV.to_host(got).astype(np.int64).astype(object) != want))
    log(event='addmul', failures=am_fail)
    ok &= am_fail == 0
    rows = []
    for name, spec in BUILTIN.items():
        if a.quick and spec.get('slow'):
            continue
        rc = run_case(spec, 'cpu', work_gb=0.5)
        rd = run_case(spec, dev, work_gb=0.5)
        rt = run_case(spec, dev, work_gb=0.5, tile_rows=7, batch_units=40, direct_rows=5, cache_gb=spec.get('cache_gb', 0.01),
                      col_block=33, eng_tile_rows=9)
        good = True
        for r in (rc, rd, rt):
            good &= (r['status'] == 'done' and r['fixed_point'] and r['dim_VD1'] == spec['dim'] and r['profile'] == spec['profile']
                     and r['contains_one'] == spec['one'] and [h[:16] for h in r['low_basis_sha256']] == spec['sha'])
        row = dict(case=name, dim=rd['dim_VD1'], profile=rd['profile'], contains_one=rd['contains_one'],
                   sha=[h[:16] for h in rd['low_basis_sha256']], dim_K=rd['dim_K'], rounds=rd['rounds'],
                   cpu_device_tiny_all_match_embedded=bool(good), evictions_tiny=rt['counts']['follower_evictions'],
                   seconds=[rc['seconds'], rd['seconds'], rt['seconds']])
        log(event='case', **row)
        rows.append(row)
        ok &= good
    rep = dict(device=dev, arith=ar, addmul_failures=am_fail, cases=rows, passed=bool(ok), guard_hits=B.GUARD_HITS[0],
               source_sha256=dict(lift=LF.LIFT_SHA256, engine=B.ENGINE_SHA256, driver=DRIVER_SHA256))
    log(event='selftest', passed=bool(ok))
    if a.out:
        atomic_write_json(a.out, rep)
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['lift', 'plan', 'selftest'])
    ap.add_argument('--ckpt', help='engine checkpoint directory of the finished degree-D closure')
    ap.add_argument('--instance', help='instance JSON (with --D: compute V_D first; with --ckpt: only recorded)')
    ap.add_argument('--D', type=int, default=None, help='degree of the given closure (the lift computes D+1)')
    ap.add_argument('--k', type=int)
    ap.add_argument('--out')
    ap.add_argument('--status')
    ap.add_argument('--device', default=None, choices=['cpu', 'gpu', 'fake-gpu'])
    ap.add_argument('--work-gb', type=float, default=2.0, help='budget for batch buffers (device memory on gpu)')
    ap.add_argument('--cache-gb', type=float, default=None, help='byte budget of the follower cache (default: no limit, or derived from --dev-gb)')
    ap.add_argument('--dev-gb', type=float, default=None, help='total device budget: cache = dev - work - K worst case - 1.5 GB - basis on the device')
    ap.add_argument('--basis-on', default='device', choices=['device', 'host'], help='where the V_D tiles live (host: batch rows are transferred)')
    ap.add_argument('--tile-mb', type=float, default=256.0)
    ap.add_argument('--tile-rows', type=int, default=None)
    ap.add_argument('--col-block', type=int, default=None)
    ap.add_argument('--batch-units', type=int, default=None, help='max units (regular orbit, shift) per batch (default: from --work-gb)')
    ap.add_argument('--direct-rows', type=int, default=None)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--stop-on-one', action='store_true', help='stop as soon as 1 is in the space found (dimension/profile are then lower bounds)')
    ap.add_argument('--to-fixed-point', action='store_true', help='continue to the fixed point even with --stop-on-one (the default behaviour)')
    ap.add_argument('--include-low', action='store_true', help='also process the low rows of V_D (not needed when V_D is closed; makes the result independent of that)')
    ap.add_argument('--fall-every', type=int, default=8, help='shift the pending falls after this many unit batches (any value gives the same result)')
    ap.add_argument('--lift-ckpt', help='directory for checkpoints of the lift (K_0, K_1, position)')
    ap.add_argument('--ckpt-minutes', type=float, default=0.0)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--force-ckpt', action='store_true', help='accept an engine checkpoint whose closure is not marked complete')
    ap.add_argument('--log-seconds', type=float, default=20.0)
    ap.add_argument('--test-drop-shift', type=int, default=None, help='TEST HOOK (negative control): drop the rows of one shift index')
    ap.add_argument('--test-no-fall-shifts', action='store_true', help='TEST HOOK (negative control): do not shift the falls')
    ap.add_argument('--free', help='plan: free columns of V_D per part, e.g. 25136,25120')
    ap.add_argument('--free-top', help='plan: free columns of weight D per part, e.g. 703,702')
    ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    a.device = a.device or ('gpu' if a.cmd == 'selftest' else 'cpu')
    if a.cmd == 'lift':
        cmd_lift(a)
    elif a.cmd == 'plan':
        if a.k is None or a.D is None or not a.free or not a.free_top:
            ap.error('plan needs --k --D --free a,b --free-top a,b')
        cmd_plan(a)
    else:
        cmd_selftest(a)


if __name__ == '__main__':
    main()
