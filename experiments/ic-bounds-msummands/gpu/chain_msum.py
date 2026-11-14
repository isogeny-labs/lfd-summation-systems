#!/usr/bin/env python3
"""Pod-side chain for the m-summand cells that were NOT run locally (experiments/ic-bounds-msummands/gpu/PLAN.md).
Modelled on experiments/ic-bounds-binary-deep/pod/chain_bin.py; same conventions (per-system log, status file, results.jsonl,
`timeout` and `skipped: budget` recorded, never silent).  NOTHING HERE HAS BEEN RUN: it is the frozen plan.

Layout on the pod (WORKSPACE, default /workspace):
  code/   csd/ (whole package, the repository's), scripts/ic_bounds_msum_run.py, scripts/ic_bounds_bin_big_run.py, scripts/ic_bounds_m3_depth.py,
          chain_msum.py, gpu/cells.json, gpu/specs/*.json
  in/     instance files (built in step 1 by `ic_bounds_msum_run.py instance`; hashes compared with the frozen ones in cost.json)
  out/    records
Steps: (1) engine self-test on the device (binary: ic_bounds_bin_big_run.py selftest; stops the chain on failure);
       (2) for every cell in order: build the instance (binary: ic-bin-instance/1 for the existing GPU engine; prime: ic-prime-msum-instance/1),
           compare its system/generator hash with the frozen one, run the ladder:
           binary -> scripts/ic_bounds_bin_big_run.py ladder (bit-packed GF(2), cupy);
           prime  -> scripts/ic_bounds_msum_run.py lfd-file (CompactClosure, CPU, float64; the m-block GPU engine does not exist, see PLAN.md).
Exit 75 on SIGTERM after recording the state.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import time

ap = argparse.ArgumentParser()
ap.add_argument('--cells', default=None)
ap.add_argument('--device', default='gpu', choices=['gpu', 'fake-gpu', 'cpu'])
ap.add_argument('--deadline', type=float, default=1e12)
ap.add_argument('--dry', action='store_true')
ap.add_argument('--only', nargs='*', default=None)
ap.add_argument('--family', default=None, choices=[None, 'binary', 'prime'])
ap.add_argument('--batch', type=int, default=1024)
ap.add_argument('--work-gb', type=float, default=20.0)
ap.add_argument('--tile-mb', type=float, default=1024.0)
ap.add_argument('--skip-selftest', action='store_true')
ap.add_argument('--code-dir', default=None)
args = ap.parse_args()

WS = os.environ.get('WORKSPACE', '/workspace')
OUT, IN, LOGS = (os.path.join(WS, d) for d in ('out', 'in', 'logs'))
CODE = args.code_dir or os.path.join(WS, 'code')
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
        except Exception:
            pass
    state.update(finished=True, ended=now())
    save()
    sys.exit(75)


def run(name, argv, limit=None, env_extra=None):
    cmd = ['python3'] + argv
    if args.dry:
        print(name + ': ' + ' '.join(cmd))
        return 0
    step = dict(name=name, cmd=cmd, start=now(), rc=None)
    state['steps'].append(step)
    save()
    env = dict(os.environ, **(env_extra or {}))
    with open(os.path.join(LOGS, 'chain-%s.log' % name), 'ab') as log:
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=os.path.join(CODE, 'scripts'), env=env)
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
    except Exception:
        return None


def main():
    signal.signal(signal.SIGTERM, on_term)
    cells_path = args.cells or os.path.join(CODE, 'gpu', 'cells.json')
    cells = load(cells_path)
    frozen = load(os.path.join(CODE, 'gpu', 'cost.json')) or {}
    frozen_hash = {r['name']: r.get('system_hash') or r.get('generator_hash') for fam in ('binary', 'prime') for r in frozen.get(fam, [])}
    if not args.dry:
        for d in (OUT, LOGS, IN):
            os.makedirs(d, exist_ok=True)
    if not args.skip_selftest and args.family in (None, 'binary'):
        rc = run('selftest', ['ic_bounds_bin_big_run.py', 'selftest', '--device', args.device, '--out', os.path.join(OUT, 'selftest-%s.json' % args.device)],
                 env_extra={})
        if rc != 0 and not args.dry:
            state['aborted'] = 'self-test failed (rc %s); nothing else was run' % rc
            state.update(finished=True, ended=now())
            save()
            return 1
    for cell in cells['order']:
        cs = state['cells'].setdefault(cell['cell'], dict(purpose=cell.get('purpose'), systems={}))
        for s in cell['systems']:
            name = s['name']
            if args.only and name not in args.only:
                continue
            if args.family and s['family'] != args.family:
                continue
            left = args.deadline - (time.time() - T0)
            if left < s.get('est_seconds', 0):
                cs['systems'][name] = dict(status='skipped: budget', est_seconds=s.get('est_seconds'), seconds_left=round(left))
                save()
                continue
            spec = os.path.join(CODE, 'gpu', s['spec'])
            rc = run('build-' + name, ['ic_bounds_msum_run.py', 'instance', '--spec', spec, '--outdir', IN], limit=7200)
            inst_path = os.path.join(IN, name + '.json')
            inst = load(inst_path) if not args.dry else dict(system_hash=frozen_hash.get(name), generator_hash=frozen_hash.get(name), msum_system_hash=frozen_hash.get(name))
            got = (inst or {}).get('msum_system_hash') or (inst or {}).get('generator_hash')
            want = frozen_hash.get(name)
            if rc != 0 or (want is not None and got != want):
                cs['systems'][name] = dict(status='instance mismatch or build failure', rc=rc, got=got, want=want)
                save()
                continue
            out = os.path.join(OUT, name + '.json')
            if s['family'] == 'binary':
                argv = ['ic_bounds_bin_big_run.py', 'ladder', '--instance', inst_path, '--max-degree', str(s['max_degree']), '--device', args.device,
                        '--batch', str(args.batch), '--work-gb', str(args.work_gb), '--tile-mb', str(args.tile_mb), '--tiles-on', 'device',
                        '--status', os.path.join(OUT, name + '.status.json'), '--out', out]
            else:
                argv = ['ic_bounds_msum_run.py', 'lfd-file', '--instance', inst_path, '--out', out]
            cs['systems'][name] = dict(status='running', start=now())
            save()
            rc = run(name, argv, limit=s.get('time_limit_s'))
            rec = load(out) if not args.dry else dict(status='dry')
            row = dict(status=('timeout' if rc == 'timeout' else 'done' if rc == 0 else 'failed rc %s' % rc), rc=rc,
                       exact_lfd=(rec or {}).get('exact_lfd'), engine_status=(rec or {}).get('status'),
                       strict_falls=(rec or {}).get('strict_falls'), seconds=(rec or {}).get('seconds'), end=now())
            cs['systems'][name] = row
            if not args.dry:
                with open(os.path.join(OUT, 'results.jsonl'), 'a') as f:
                    f.write(json.dumps(dict(name=name, cell=cell['cell'], **row)) + '\n')
            save()
    state.update(finished=True, ended=now())
    save()
    return 0


if __name__ == '__main__':
    sys.exit(main())
