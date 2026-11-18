#!/usr/bin/env python3
"""Pod job: split-presentation closures V_D (m >= 4 summands) on a GPU with CuPy.

One command does, in this order, and writes everything under --out:
  env.json        python / numpy / cupy / CUDA / device name and memory, nvidia-smi, source hashes
  selftest.json   exact-arithmetic self-test of the device layer (csd.ic_bounds_big.selftest_arith)
  gate.json       equality gate on the device against the CPU reference, with negative controls
                  (scripts/ic_split_gpu_gate.py); the job stops here if it fails
  closures.jsonl  one record per closure of the target cells (timings per phase, device memory
                  high-water mark, counted multiply-adds, status), flushed as soon as it exists
  job.log         every line printed, with a UTC time stamp

    python3 scripts/ic_split_gpu_pod.py --device gpu --out /workspace/gpu-split-m \
        --cells '6,2:3-7;5,3:3-8;7,2:3-8' --device-gb 130 --budget-s 14400

A closure stops itself when its stored basis passes the memory budget (--device-gb minus the working
set) or when it passes --budget-s; it is then recorded as 'stopped: ...' (with a checkpoint if
--checkpoint-gb allows its size) and the ladder of that cell ends.  The upper bound printed before
each closure is the worst case (1 enters, about n_q^2/4 cells per part); a closure without 1 uses
rank x (n_q - rank) cells, usually far less.  The batch is cut so that the working set stays below
15 % of the planned device memory.  Any exception is recorded with its traceback and the
job goes on with the next cell.  --device cpu and --device fake-gpu run the same code with numpy.

Instances: the certified instance of experiments/ic-bounds-split-m for cells that have one; else a
seeded unit instance by the rule of scripts/run_ic_split_m.py (prime --prime, seed --seed, target
redrawn until S_{m+1} has no zero on the grid).  The ladder of a cell stops when 1 enters.
"""
import argparse, datetime, json, math, os, pathlib, platform, subprocess, sys, time, traceback

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

LOG = None


def log(msg):
    line = f'{datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M:%S")} {msg}'
    print(line, flush=True)
    if LOG:
        LOG.write(line + '\n')
        LOG.flush()


def columns(m, k, D):
    return sum(math.comb(m * k, w) * math.comb(D - w + m - 2, m - 2) for w in range(min(m * k, D) + 1))


def predicted_device_gib(n, batch, split=True):
    """Upper bound: basis n_q^2/4 cells per part (parts about 0.56 n and 0.44 n with the split) plus
    working set of about six batch x n_q arrays and one tile of temporary."""
    parts = [0.56 * n, 0.44 * n] if split else [n]
    basis = sum(q * q / 4 for q in parts) * 8
    work = 6 * batch * parts[0] * 8 + 4096 * parts[0] * 8 * 2
    return (basis + work) / 2 ** 30


def environment(device):
    env = dict(python=sys.version, platform=platform.platform(), argv=sys.argv, cwd=os.getcwd())
    import numpy as np
    env['numpy'] = np.__version__
    try:
        env['git'] = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    except Exception as e:
        env['git'] = repr(e)
    if device == 'gpu':
        import cupy as cp
        env['cupy'] = cp.__version__
        dev = cp.cuda.Device()
        free, total = cp.cuda.runtime.memGetInfo()
        env['cuda_runtime'] = cp.cuda.runtime.runtimeGetVersion()
        env['device_id'] = dev.id
        env['device_total_gib'] = round(total / 2 ** 30, 2)
        env['device_free_gib'] = round(free / 2 ** 30, 2)
        try:
            env['nvidia_smi'] = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total,driver_version',
                                                '--format=csv,noheader'], capture_output=True, text=True).stdout.strip()
        except Exception as e:
            env['nvidia_smi'] = repr(e)
    return env


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument('--device', default='gpu', choices=['gpu', 'cpu', 'fake-gpu'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--cells', default='6,2:3-7;5,3:3-8;7,2:3-8', help="'m,k:Dmin-Dmax;...'")
    ap.add_argument('--gate', default='4,2:3-6;5,2:3-6;4,3:3-6', help='cells of the on-device equality gate; empty to skip')
    ap.add_argument('--gate-recorded', default='4,3:7', help='gate cells compared with committed reference records')
    ap.add_argument('--device-gb', type=float, default=None, help='device memory the job may plan for (default: 92 %% of free)')
    ap.add_argument('--budget-s', type=float, default=14400, help='time budget per closure')
    ap.add_argument('--batch', type=int, default=4096)
    ap.add_argument('--tile-rows', type=int, default=4096)
    ap.add_argument('--sub', type=int, default=64)
    ap.add_argument('--checkpoint-gb', type=float, default=40, help='write checkpoints only for closures predicted below this size')
    ap.add_argument('--checkpoint-every-s', type=float, default=1800)
    ap.add_argument('--prime', type=int, default=114689)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--no-split', action='store_true')
    ap.add_argument('--resume', action='store_true', help='resume closures from checkpoints in --out and skip recorded ones')
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    LOG = open(out / 'job.log', 'a')
    log(f'start device={a.device} cells={a.cells} out={out}')
    status = 0
    try:
        env = environment(a.device)
    except Exception:
        log('FATAL environment: ' + traceback.format_exc())
        return 2
    from csd import ic_bounds_big as B
    from csd import ic_bounds_split_gpu as G
    import ic_split_gpu_gate as GATE
    import run_ic_split_m as RUN
    env['engine_sha256'] = G.ENGINE_SHA256
    env['device_layer_sha256'] = B.ENGINE_SHA256
    (out / 'env.json').write_text(json.dumps(env, indent=1) + '\n')
    log('env ' + json.dumps({k: v for k, v in env.items() if k not in ('argv', 'python')}))
    if a.device_gb is None:
        a.device_gb = 0.92 * env['device_free_gib'] if a.device == 'gpu' else 8.0
    log(f'device memory planned for: {a.device_gb:.1f} GiB')

    # 1. arithmetic self-test on the device
    try:
        st = B.selftest_arith(a.device)
        (out / 'selftest.json').write_text(json.dumps(st, indent=1, default=str) + '\n')
        log('selftest ' + json.dumps(st, default=str)[:600])
        if not st.get('passed'):
            log('FATAL selftest failed')
            return 3
    except Exception:
        log('FATAL selftest raised: ' + traceback.format_exc())
        return 3

    # 2. equality gate on the device (reference on the host CPU)
    if a.gate:
        arms = f'{a.device}:split,{a.device}:plain'
        cmd = [sys.executable, str(ROOT / 'scripts/ic_split_gpu_gate.py'), '--live', a.gate, '--recorded', a.gate_recorded,
               '--arms', arms, '--controls', '--batch', '1024', '--out', str(out / 'gate.json')]
        log('gate: ' + ' '.join(cmd))
        r = subprocess.run(cmd, capture_output=True, text=True)
        for line in (r.stderr + r.stdout).splitlines():
            log('gate| ' + line)
        if r.returncode:
            log(f'FATAL equality gate failed (exit {r.returncode}); no target cell is run')
            return 4
        log('gate passed')

    # 3. target cells
    done = set()
    rec_path = out / 'closures.jsonl'
    if a.resume and rec_path.exists():
        for line in rec_path.read_text().splitlines():
            r = json.loads(line)
            if r.get('status') in ('one', 'closed'):
                done.add((r['m'], r['k'], r['D']))
    for item in filter(None, a.cells.split(';')):
        mk, dr = item.split(':')
        m, k = map(int, mk.split(','))
        lo, hi = (map(int, dr.split('-')) if '-' in dr else (int(dr), int(dr)))
        try:
            t0 = time.time()
            if (GATE.CELLS / f'cell-m{m}-k{k}.jsonl').exists() or (m, k) == (5, 3):
                p, gens, _ = GATE.instance(m, k)
                source = 'certified instance of experiments/ic-bounds-split-m'
            else:
                from csd import ic_dreg as I, ic_bounds_split_m as S
                E, inst = RUN.unit_instance(a.prime, m, k, a.seed)
                p = a.prime
                gens, _ = S.split_generators_m(E, inst['t'], k, m)
                source = f'unit instance p={p} seed={a.seed} a={inst["a"]} b={inst["b"]} t={inst["t"]} attempt={inst["attempt"]}'
            log(f'cell m={m} k={k}: {source} ({time.time() - t0:.1f}s)')
        except Exception:
            log(f'ERROR instance m={m} k={k}: ' + traceback.format_exc())
            status = 1
            continue
        for D in range(lo, hi + 1):
            if (m, k, D) in done:
                log(f'm={m} k={k} D={D}: already recorded, skipped')
                continue
            n = columns(m, k, D)
            need = predicted_device_gib(n, min(a.batch, 2048), not a.no_split)
            base = dict(m=m, k=k, p=p, D=D, monomials_le_D=n, instance=source, predicted_device_gib_upper=round(need, 2))
            log(f'm={m} k={k} D={D}: {n} columns, upper bound of device memory {need:.1f} GiB '
                f'(a closure without 1 normally uses far less; the run stops itself at the budget)')
            ck = str(out / f'ckpt-m{m}-k{k}-D{D}.npz') if need <= a.checkpoint_gb else None
            batch = max(64, min(a.batch, int(0.15 * a.device_gb * 2 ** 30 / (6 * 0.56 * n * 8))))
            work_gib = 6 * batch * 0.56 * n * 8 / 2 ** 30
            log(f'  batch {batch} (working set about {work_gib:.1f} GiB), checkpoint: {ck}')
            try:
                t0 = time.time()
                rec, eng = G.closure(gens, m, k, p, D, device=a.device, split=not a.no_split, batch=batch,
                                     tile_rows=a.tile_rows, sub=a.sub, budget_s=a.budget_s,
                                     mem_budget_gb=max(0.5, a.device_gb - work_gib - 2 * a.tile_rows * 0.56 * n * 8 / 2 ** 30),
                                     checkpoint=ck, checkpoint_every_s=a.checkpoint_every_s, resume=a.resume, log=log)
                rec.update(base, wall_seconds=round(time.time() - t0, 2))
                del eng
            except Exception:
                rec = dict(base, status='error', traceback=traceback.format_exc())
                log(f'ERROR m={m} k={k} D={D}: ' + rec['traceback'])
                status = 1
            try:
                B.DV.free_pool()
            except Exception:
                pass
            with open(rec_path, 'a') as f:
                f.write(json.dumps(rec) + '\n')
            st = rec.get('stats', {})
            log(f'RESULT m={m} k={k} D={D} status={rec.get("status")} dim={rec.get("dim_V_D")} '
                f'one={rec.get("contains_one")} verified={rec.get("verified_closed")} t={rec.get("seconds")}s '
                f'peak_device={rec.get("peak_device_gib")}GiB reduce={st.get("s_reduce")}s update={st.get("s_update")}s '
                f'rref={st.get("s_rref")}s products={st.get("s_products")}s verify={st.get("s_verify")}s')
            if 'ops' in rec:
                print(json.dumps(dict(rec['ops'], cell=[m, k, D])), flush=True)
            if rec.get('status') not in ('one', 'closed'):
                log(f'm={m} k={k}: ladder stopped at D={D} ({rec.get("status")})')
                if rec.get('status') != 'error':
                    status = max(status, 1)
                break
            if rec.get('contains_one'):
                log(f'm={m} k={k}: 1 enters at D={D}: exact value {D} if every smaller closure above is verified closed')
                break
    log(f'end status={status}')
    return status


if __name__ == '__main__':
    sys.exit(main())
