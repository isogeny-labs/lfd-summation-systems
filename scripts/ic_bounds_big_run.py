"""Driver for csd/ic_bounds_big.py (closure engine for V_D, CPU numpy / GPU cupy).

Subcommands
  closure   run (or resume) one closure from an instance JSON; no repository imports
  plan      predicted peak memory for k, D, device
  selftest  arithmetic exactness + CPU/device agreement on built-in small instances + device report
  instance  (LOCAL ONLY, imports csd.ic_dreg and csd.ic_curve_sweep like the reference driver) draw a target,
            enumerate the solutions Z, write the instance JSON that `closure --instance` consumes

The engine file is found next to this script (remote layout: both files in /workspace/code) or in ../csd/
(repository layout).  Remote layout of paths: inputs /workspace/in, outputs /workspace/out, checkpoints
/workspace/ckpt; every path is a flag.

Exit codes: 0 done (closed, clean sweep) ; 1 selftest/other failure ; 2 stopped because 1 is in V_D
(--stop-on-one) ; 3 sweep failed and --no-repair ; 75 SIGTERM (checkpointed, resume with --resume).
"""
import argparse
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
for cand in (HERE / 'ic_bounds_big.py', HERE.parent / 'csd' / 'ic_bounds_big.py'):
    if cand.exists():
        ENGINE_PATH = cand
        sys.path.insert(0, str(cand.parent))
        break
else:  # pragma: no cover
    raise SystemExit('ic_bounds_big.py not found next to this script or in ../csd/')


def _set_blas_threads(argv):
    """BLAS thread count must be fixed before numpy is imported: --threads overrides the environment;
    otherwise the environment (OPENBLAS_NUM_THREADS etc.) wins, and without it the default is the number of
    CPUs this process may run on (len(os.sched_getaffinity(0)), not os.cpu_count())."""
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
        try:        # cgroup CPU quota (v2, then v1): the affinity count can be the host's cores in a limited pod
            q, per = open('/sys/fs/cgroup/cpu.max').read().split()[:2]
            if q != 'max':
                d = max(1, min(d, -(-int(q) // int(per))))
        except Exception:
            try:
                q = int(open('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read())
                per = int(open('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read())
                if q > 0:
                    d = max(1, min(d, -(-q // per)))
            except Exception:
                pass
        for key in keys:
            os.environ.setdefault(key, str(d))


_set_blas_threads(sys.argv)
import numpy as np  # noqa: E402
import ic_bounds_big as B  # noqa: E402

DRIVER_SHA256 = B.file_sha256(__file__)
ROOT = HERE.parent

# reference values (csd/ic_bounds_d5_verify.py, gens='all', split=True) for the built-in selftest instances;
# sha = first 16 hex digits of each part's basis_sha256
BUILTIN = {
    "unit-k6-D5": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 6, "D": 5, "rank": 1586, "one": True, "sha": ["126f88ebc88184be", "0527370c25f10d1f"]},
    "unit-k6-D6": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 6, "D": 6, "rank": 2510, "one": True, "sha": ["8250bccce8cfeb13", "783b2f9e75e71522"]},
    "unit-k7-D5": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 7, "D": 5, "rank": 3473, "one": True, "sha": ["c2dbabeb09304e51", "d98d9bee18090030"]},
    "unit-k7-D6": {"p": 12289, "a": 1, "b": 0, "X": 9362, "k": 7, "D": 6, "rank": 6476, "one": True, "sha": ["54a5d2ccfd9ce966", "0835a149b84b2d34"]},
    "proper-k6-D5": {"p": 15361, "a": 4021, "b": 9917, "X": 4491, "k": 6, "D": 5, "rank": 1584, "one": False, "sha": ["653448a619f1dd2c", "8a2ba8a7fbf9e19d"]},
    "proper-k6-D6": {"p": 15361, "a": 4021, "b": 9917, "X": 4491, "k": 6, "D": 6, "rank": 2508, "one": False, "sha": ["4111530e47a83425", "83b0db79f53e4fea"]},
    "proper-k7-D5": {"p": 15361, "a": 4021, "b": 9917, "X": 7776, "k": 7, "D": 5, "rank": 3469, "one": False, "sha": ["23d4523919d476f3", "19adca9ad6b4846d"]},
    "proper-k7-D6": {"p": 15361, "a": 4021, "b": 9917, "X": 7776, "k": 7, "D": 6, "rank": 6472, "one": False, "sha": ["6ffba6a7ec1807f1", "3c2359db5c622fd9"]},
}

T0 = time.time()


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
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


def engine_kwargs(a):
    return dict(device=a.device, split=not a.nosplit, gens=a.gens, work_gb=a.work_gb, tile_mb=a.tile_mb,
                tile_rows=a.tile_rows, col_block=a.col_block, batch=a.batch, threads=a.threads,
                compact_frac=a.compact_frac, stop_on_one=a.stop_on_one, t_random=a.t_random,
                sweep_chunk=a.sweep_chunk, dev_gb=a.dev_gb, tile_budget_mb=a.tile_budget_mb, ckpt=a.ckpt, ckpt_minutes=a.ckpt_minutes, ckpt_slots=a.ckpt_slots,
                status=a.status, no_repair=a.no_repair, driver_sha=DRIVER_SHA256, log_seconds=a.log_seconds,
                withhold=tuple(int(x) for x in a.withhold.split(':')) if a.withhold else None)


# ----------------------------------------------------------------------------- closure

def evaluation_fields(eng, Z, D):
    """mirror of the reference driver: points, target, soundness on Z, evaluation rank, deficient."""
    k, p = eng.k, eng.p
    B._Live.limit = None       # the fake-device budget applies to the closure, not to these host-side checks
    rec = dict(points=len(Z), target=eng.L.c - len(Z))
    if Z:
        pts = [tuple(z) for z in Z]
        X = eng.random_elements_full(4, seed=7)
        ev = B.evaluate_full(eng.L, X, pts, p, k)
        rec['soundness_vanish_on_Z'] = bool(np.all(ev == 0))
        Ll = B.Layout(k, D - 1, False) if D >= 1 else None
        rec['eval_rank_below_D'] = B.eval_rank(Ll, pts, p, k, D - 1,
                                               dict(tile_bytes=int(eng.tile_bytes), work_bytes=int(eng.work_bytes),
                                                    tile_rows=eng.tile_rows_opt))
    return rec


def cmd_closure(a):
    inst = json.load(open(a.instance))
    k, p = inst['k'], inst['prime']
    g = parse_gen(inst['generator'])
    D = a.D
    inst_sha = B.file_sha256(a.instance)
    source_id = hashlib.sha256(json.dumps([k, p, D, not a.nosplit, a.gens, sorted(g.items()), inst_sha]).encode()).hexdigest()
    kw = engine_kwargs(a)
    eng = B.Engine(g, k, p, D, source_id=source_id, **kw)
    resumed = False
    if a.resume:
        if not a.ckpt:
            raise SystemExit('--resume needs --ckpt DIR')
        resumed = eng.load_checkpoint(force=a.force_resume)
        if not resumed:
            log(event='no_checkpoint_found_starting_fresh')
    eng.install_sigterm()
    log(event='start', k=k, D=D, p=p, device=a.device, columns=eng.L.c, part_columns=[P.c for P in eng.parts],
        instance=os.path.basename(a.instance), resumed=resumed, engine_sha256=B.ENGINE_SHA256,
        driver_sha256=DRIVER_SHA256, threads=eng.threads)
    status = eng.run()
    if status == 'sigterm':
        log(event='exit_sigterm_checkpointed')
        sys.exit(B.EXIT_SIGTERM)
    rec = eng.record()
    rec['status'] = status
    rec.update(instance_file=os.path.basename(a.instance), instance_sha256=inst_sha,
               a=inst['a'], b=inst['b'], seed=inst['seed'], want=inst.get('want'),
               unit_nonresidue=inst.get('unit_nonresidue'), instance_target=inst['target'],
               attempts=inst.get('attempts'), instance_k=k, resumed_from_checkpoint=resumed,
               source_sha256=dict(engine=B.ENGINE_SHA256, driver=DRIVER_SHA256, instance=inst_sha))
    if status == 'one':
        rec['contains_one'] = True
        rec['closure_complete'] = False
    if status == 'done' and 'Z' in inst and not a.no_eval:
        rec.update(evaluation_fields(eng, inst['Z'], D))
        rec['deficient'] = rec['rank'] < rec['target']
    rec['seconds'] = round(eng.seconds(), 1)
    log(**{key: v for key, v in rec.items() if key not in ('basis_sha256', 'part_profiles')}, event='result')
    if a.out:
        atomic_write_json(a.out, rec)
        log(event='saved', file=a.out)
    sys.exit({'done': 0, 'one': 2, 'failed_sweep': 3}.get(status, 1))


# ----------------------------------------------------------------------------- plan

def cmd_plan(a):
    r = B.plan(a.k, a.D, a.device, a.work_gb, a.tile_mb, split=not a.nosplit, compact_frac=a.compact_frac, batch=a.batch, dev_gb=a.dev_gb)
    print(json.dumps(r, indent=1))


# ----------------------------------------------------------------------------- selftest

def device_report(device):
    out = dict(device=device, numpy=np.__version__, python=sys.version.split()[0])
    if device == 'gpu':
        import cupy as cp
        out['cupy'] = cp.__version__
        rt = cp.cuda.runtime
        out['cuda_runtime'] = rt.runtimeGetVersion()
        out['cuda_driver'] = rt.driverGetVersion()
        props = rt.getDeviceProperties(cp.cuda.Device().id)
        name = props['name']
        out['name'] = name.decode() if isinstance(name, bytes) else str(name)
        out['total_memory_gb'] = round(props['totalGlobalMem'] / (1 << 30), 1)
        free, total = rt.memGetInfo()
        out['free_memory_gb'] = round(free / (1 << 30), 1)
        out['compute_capability'] = f"{props['major']}.{props['minor']}"
        try:
            out['cublas_version'] = int(cp.cuda.cublas.getVersion(cp.cuda.device.get_cublas_handle()))
        except Exception as e:  # pragma: no cover
            out['cublas_version'] = f'unavailable: {e}'
    return out


def run_small(name, spec, device, **extra):
    g = B.chain_generator(spec['p'], spec['a'], spec['b'], spec['X'], spec['k'])
    eng = B.Engine(g, spec['k'], spec['p'], spec['D'], device=device, log=False, work_gb=extra.get('work_gb', 1.0),
                   batch=extra.get('batch', 64), tile_rows=extra.get('tile_rows'), source_id=name,
                   tile_budget_mb=extra.get('tile_budget_mb'))
    st = eng.run()
    r = eng.record()
    return st, r


def cmd_selftest(a):
    ok = True
    dev = a.device
    rep = dict(report=device_report(dev))
    log(event='device', **rep['report'])
    if dev == 'gpu':
        B.use_device('gpu')
        d = B.DV
        used, total = d.mem_info()
        log(event='memory', used_mb=used >> 20, total_mb=total >> 20)
    ar = B.selftest_arith(dev)
    log(event='arith_device', **{k: v for k, v in ar.items()})
    ar_cpu = B.selftest_arith('cpu')
    log(event='arith_cpu', **{k: v for k, v in ar_cpu.items()})
    ok &= ar['passed'] and ar_cpu['passed']
    rows = []
    names = [n for n in BUILTIN if not (a.quick and n.endswith('k7-D6'))]
    for name in names:
        spec = BUILTIN[name]
        stc, rc = run_small(name, spec, 'cpu')
        std, rd = run_small(name, spec, dev, tile_rows=7, batch=5) if a.tiny else run_small(name, spec, dev)
        same = (rc['rank'] == rd['rank'] and rc['profile'] == rd['profile'] and rc['basis_sha256'] == rd['basis_sha256']
                and rc['part_ranks'] == rd['part_ranks'] and rc['contains_one'] == rd['contains_one'])
        ref_ok = (rc['rank'] == spec['rank'] and rc['contains_one'] == spec['one']
                  and [h[:16] for h in rc['basis_sha256']] == spec['sha'])
        closed = bool(rc['closed'] and rd['closed'])
        row = dict(case=name, rank=rd['rank'], profile=rd['profile'], contains_one=rd['contains_one'], cpu_eq_device=same,
                   matches_reference=ref_ok, closed_cpu=rc['closed'], closed_device=rd['closed'],
                   sha=[h[:16] for h in rd['basis_sha256']], seconds_cpu=rc['seconds'], seconds_device=rd['seconds'])
        log(event='case', **row)
        rows.append(row)
        ok &= same and ref_ok and closed
    # forced spilling: tiles of 64 rows and a tile budget of 13 MB / 5 MB (one part needs at most 10 MB / 3.4 MB, both 20 / 7),
    # device run against the cpu run; the spill counters must be nonzero
    for name, tbm in (('proper-k7-D6', 13), ('unit-k7-D5', 5)):
        spec = BUILTIN[name]
        stc, rc = run_small(name, spec, 'cpu')
        std, rd = run_small(name, spec, dev, tile_rows=64, tile_budget_mb=tbm, work_gb=0.5, batch=48)
        sp = rd['spill']
        same = (rc['rank'] == rd['rank'] and rc['profile'] == rd['profile'] and rc['basis_sha256'] == rd['basis_sha256']
                and rc['part_ranks'] == rd['part_ranks'] and rc['contains_one'] == rd['contains_one'])
        ref_ok = (rc['rank'] == spec['rank'] and [h[:16] for h in rc['basis_sha256']] == spec['sha'])
        row = dict(case='spill-' + name, rank=rd['rank'], cpu_eq_device=same, matches_reference=ref_ok, closed_device=rd['closed'],
                   spills=sp['spills'], fetches=sp['fetches'], d2h_mb=sp['d2h_bytes'] >> 20, h2d_mb=sp['h2d_bytes'] >> 20,
                   seconds_device=rd['seconds'])
        log(event='case', **row)
        rows.append(row)
        ok &= same and ref_ok and rd['closed'] and sp['spills'] > 0 and sp['fetches'] > 0
    rep.update(arith_device=ar, arith_cpu=ar_cpu, cases=rows, passed=bool(ok))
    log(event='selftest', passed=bool(ok))
    if a.out:
        atomic_write_json(a.out, rep)
    sys.exit(0 if ok else 1)


# ----------------------------------------------------------------------------- instance (local)

def cmd_instance(a):
    sys.path.insert(0, str(ROOT))
    import random
    from csd import ic_dreg as I
    from csd import ic_curve_sweep as S
    p, k = a.prime, a.k
    rng = random.Random(a.seed)
    E = I.PrimeCurve(p, a.a, a.b)
    attempts = 0
    lim = 5000 if a.unit_nonresidue else 200
    for attempt in range(lim):
        attempts = attempt + 1
        R = S.random_rational_target(E, rng, p)
        if a.unit_nonresidue and pow(R[0], (p - 1) // 2, p) != p - 1:
            continue
        Z = B.solutions(p, a.a % p, a.b % p, R[0], k)
        if Z or not a.want:
            break
    g = B.chain_generator(p, a.a % p, a.b % p, R[0], k)
    assert B.s3_coeffs(p, a.a % p, a.b % p, R[0]) == I.s3_prime_coeffs(E, R[0])
    if a.unit_nonresidue:
        assert pow(R[0], (p - 1) // 2, p) == p - 1
    rec = dict(k=k, prime=p, a=a.a, b=a.b, seed=a.seed, want=a.want, unit_nonresidue=a.unit_nonresidue,
               target=list(R), target_x_is_square=pow(R[0], (p - 1) // 2, p) == 1, attempts=attempts,
               points=len(Z), Z=[list(z) for z in Z],
               columns={str(D): B.count_columns(k, D) for D in (4, 5, 6)}, Z_bound_2_2k=2 << k,
               generator={f'{e1},{e2}': v for (e1, e2), v in sorted(g.items())})
    if a.eval_rank and not Z:
        rec['eval_rank_below_D'], rec['onto'], rec['eval_rank_D'] = 0, True, a.D
    if a.eval_rank and Z:
        B.use_device('cpu')
        L = B.Layout(k, a.D - 1, False)
        er = B.eval_rank(L, [tuple(z) for z in Z], p, k, a.D - 1, {})
        rec['eval_rank_below_D'] = er
        rec['onto'] = er == len(Z)
        rec['eval_rank_D'] = a.D
    rec['source_sha256'] = dict(engine=B.ENGINE_SHA256, driver=DRIVER_SHA256)
    log(**{key: v for key, v in rec.items() if key not in ('generator', 'Z')}, Z_head=rec['Z'][:3])
    atomic_write_json(a.out, rec)
    log(event='saved', file=a.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['closure', 'plan', 'selftest', 'instance'])
    ap.add_argument('--instance')
    ap.add_argument('--out')
    ap.add_argument('--k', type=int)
    ap.add_argument('--D', type=int, default=6)
    ap.add_argument('--device', default=None, choices=['cpu', 'gpu', 'fake-gpu'])
    ap.add_argument('--work-gb', type=float, default=4.0, help='budget for working buffers (device memory on gpu)')
    ap.add_argument('--tile-mb', type=float, default=256.0)
    ap.add_argument('--tile-rows', type=int, default=None)
    ap.add_argument('--col-block', type=int, default=None)
    ap.add_argument('--batch', type=int, default=512, help='max rows per (source part, shift) batch')
    ap.add_argument('--threads', type=int, default=None, help='default: len(os.sched_getaffinity(0))')
    ap.add_argument('--compact-frac', type=float, default=0.05)
    ap.add_argument('--dev-gb', type=float, default=None, help='total device budget: tiles of the inactive part spill to host memory so that tiles + work + slack fit (default: no limit)')
    ap.add_argument('--tile-budget-mb', type=float, default=None, help='test: explicit tile budget (overrides --dev-gb)')
    ap.add_argument('--nosplit', action='store_true')
    ap.add_argument('--gens', default='all', choices=['all', 'g'])
    ap.add_argument('--stop-on-one', action='store_true')
    ap.add_argument('--t-random', type=int, default=6)
    ap.add_argument('--sweep-chunk', type=int, default=256)
    ap.add_argument('--ckpt')
    ap.add_argument('--ckpt-minutes', type=float, default=0.0)
    ap.add_argument('--ckpt-slots', type=int, default=1)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--force-resume', action='store_true')
    ap.add_argument('--status')
    ap.add_argument('--log-seconds', type=float, default=20.0)
    ap.add_argument('--no-repair', action='store_true', help='sweep failure ends the run (exit 3) instead of repairing')
    ap.add_argument('--withhold', help='TEST HOOK part:index - skip the shifts of one low row (negative control)')
    ap.add_argument('--no-eval', action='store_true', help='skip the Z evaluation checks')
    # instance
    ap.add_argument('--prime', type=int, default=65537)
    ap.add_argument('--a', type=int, default=1)
    ap.add_argument('--b', type=int, default=0)
    ap.add_argument('--seed', type=int, default=10)
    ap.add_argument('--want', action='store_true')
    ap.add_argument('--unit-nonresidue', action='store_true')
    ap.add_argument('--eval-rank', action='store_true')
    # selftest
    ap.add_argument('--quick', action='store_true', help='selftest: skip the k=7, D=6 case')
    ap.add_argument('--tiny', action='store_true', help='selftest: device runs use tiny tiles and batches')
    a = ap.parse_args()
    a.device = a.device or ('gpu' if a.cmd == 'selftest' else 'cpu')
    if a.cmd == 'closure':
        if not a.instance:
            ap.error('closure needs --instance')
        cmd_closure(a)
    elif a.cmd == 'plan':
        cmd_plan(a)
    elif a.cmd == 'selftest':
        cmd_selftest(a)
    else:
        if a.k is None or not a.out:
            ap.error('instance needs --k and --out')
        cmd_instance(a)


if __name__ == '__main__':
    main()
