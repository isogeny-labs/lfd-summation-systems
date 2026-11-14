#!/usr/bin/env python3
"""Pod-side chain for the binary deep runs (experiments/ic-bounds-binary-deep/PLAN.md).

Runs in /workspace/code on the pod and needs no connection to the launching machine:
  1. the self-test of the engine on the device (through nopool_run.py = relaunch/nopool_dev_only.py: the CuPy DEVICE
     pool disabled, the pinned pool left alone); the chain stops if it fails;
  2. the cells of cells.json in their order, one system after another: the ladder D = 0, 1, 2, ... up to the
     certificate degree (fall profile, LFD, local certificate) for each instance in /workspace/in.
Every step logs to /workspace/logs/chain-<step>.log; every record is written to /workspace/out as soon as it exists
(per degree: <system>.json.D<d>.json, per system: <system>.json, plus one line in results.jsonl); the progress is in
/workspace/out/chain-status.json after every step.  A system that exceeds its time limit is killed and recorded as
`timeout`; a cell whose estimate no longer fits before --deadline is recorded as `skipped: budget` (never silently).
WORKSPACE overrides /workspace (tests).  --dry prints the commands.  --device fake-gpu / cpu are for local tests.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import time

ap = argparse.ArgumentParser()
ap.add_argument('--cells', default=None, help='cells.json (default: <WORKSPACE>/in/cells.json)')
ap.add_argument('--device', default='gpu', choices=['gpu', 'fake-gpu', 'cpu'])
ap.add_argument('--deadline', type=float, default=1e12, help='seconds after the start of the chain after which no new system is started')
ap.add_argument('--no-nopool', action='store_true', help='do not use the device-pool-off wrapper (local tests; or a pod where the pool is wanted)')
ap.add_argument('--dry', action='store_true')
ap.add_argument('--only', nargs='*', default=None, help='cell names to run (default: all in order)')
ap.add_argument('--batch', type=int, default=1024)
ap.add_argument('--work-gb', type=float, default=20.0)
ap.add_argument('--tile-mb', type=float, default=1024.0)
ap.add_argument('--tiles-on', default='device')
ap.add_argument('--tiny', action='store_true', help='tiny tiles and batches (local tests with the simulated device)')
ap.add_argument('--skip-selftest', action='store_true')
ap.add_argument('--code-dir', default=None, help='directory holding the driver (default: the directory of this file; pod layout)')
args = ap.parse_args()

WS = os.environ.get('WORKSPACE', '/workspace')
OUT, IN, LOGS = (os.path.join(WS, d) for d in ('out', 'in', 'logs'))
CODE = args.code_dir or os.path.dirname(os.path.abspath(__file__))
STATUS = os.path.join(OUT, 'chain-status.json')
T0 = time.time()
state = dict(started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), device=args.device, steps=[], cells={}, finished=False)
CUR = [None]


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save():
    if args.dry:
        return
    state['elapsed_s'] = round(time.time() - T0)
    tmp = STATUS + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, STATUS)


def on_term(sig, frm):
    state['aborted'] = 'signal %d' % sig
    if CUR[0] is not None:
        try:
            CUR[0].terminate()
        except Exception:   # noqa: BLE001
            pass
    state.update(finished=True, ended=now())
    save()
    sys.exit(75)


def run(name, argv, limit=None):
    cmd = ['python3'] + (['nopool_run.py'] if (args.device == 'gpu' and not args.no_nopool) else []) + argv
    if args.dry:
        print(name + ': ' + ' '.join(cmd))
        return 0
    step = dict(name=name, cmd=cmd, start=now(), rc=None)
    state['steps'].append(step)
    save()
    with open(os.path.join(LOGS, 'chain-%s.log' % name), 'ab') as log:
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=CODE)
        CUR[0] = p
        try:
            rc = p.wait(timeout=limit)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
            rc = 'timeout'
        CUR[0] = None
    step.update(rc=rc, end=now())
    save()
    return rc


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:   # noqa: BLE001
        return None


def main():
    signal.signal(signal.SIGTERM, on_term)
    cells_path = args.cells or os.path.join(IN, 'cells.json')
    cells = load(cells_path)
    if not args.dry:
        for d in (OUT, LOGS):
            os.makedirs(d, exist_ok=True)
    if not args.skip_selftest:
        st = ['ic_bounds_bin_big_run.py', 'selftest', '--device', args.device, '--out', os.path.join(OUT, 'selftest-%s.json' % args.device)]
        if args.tiny:
            st.append('--tiny')
        rc = run('selftest', st)
        if rc not in (0,) and not args.dry:
            state['aborted'] = 'self-test failed (rc %s); nothing else was run' % rc
            state.update(finished=True, ended=now())
            save()
            return 1
    for cell in cells['order']:
        if args.only and cell['cell'] not in args.only:
            continue
        cs = state['cells'].setdefault(cell['cell'], dict(purpose=cell.get('purpose'), systems={}))
        for s in cell['systems']:
            name = s['name']
            left = args.deadline - (time.time() - T0)
            if left < s.get('est_seconds', 0):
                cs['systems'][name] = dict(status='skipped: budget', est_seconds=s.get('est_seconds'), seconds_left=round(left))
                save()
                continue
            inst = os.path.join(IN, s['instance'])
            out = os.path.join(OUT, name + '.json')
            argv = ['ic_bounds_bin_big_run.py', 'ladder', '--instance', inst, '--max-degree', str(s['max_degree']), '--device', args.device,
                    '--batch', str(args.batch), '--work-gb', str(args.work_gb), '--tile-mb', str(args.tile_mb),
                    '--tiles-on', args.tiles_on, '--status', os.path.join(OUT, name + '.status.json'), '--out', out]
            if args.tiny:
                argv += ['--batch', '5', '--tile-rows', '8', '--col-block', '64']
            cs['systems'][name] = dict(status='running', start=now())
            save()
            rc = run(name, argv, limit=s.get('time_limit_s'))
            rec = load(out) if not args.dry else dict(status='dry', exact_lfd=None)
            row = dict(status=('timeout' if rc == 'timeout' else 'done' if rc == 0 else 'failed rc %s' % rc), rc=rc,
                       exact_lfd=(rec or {}).get('exact_lfd'), cert=(rec or {}).get('certificate_degree'),
                       engine_status=(rec or {}).get('status'), n_points=(rec or {}).get('n_points'),
                       seconds=(rec or {}).get('seconds'), quotient_ranks=(rec or {}).get('quotient_ranks'),
                       strict_falls=(rec or {}).get('strict_falls'), end=now())
            cs['systems'][name] = row
            if not args.dry:
                with open(os.path.join(OUT, 'results.jsonl'), 'a') as f:
                    f.write(json.dumps(dict(name=name, cell=cell['cell'], system_hash=(rec or {}).get('system_hash'), **row)) + '\n')
            save()
    state.update(finished=True, ended=now())
    save()
    return 0


if __name__ == '__main__':
    sys.exit(main())
