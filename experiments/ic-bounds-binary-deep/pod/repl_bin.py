#!/usr/bin/env python3
"""Pod-side wrapper for the replication of the binary deep runs (../REPLICATION.md).

Starts two chains of chain_bin.py side by side, each in its own workspace so that no file of one can be taken for
a file of the other:
  r2  every system of cells.json on the GPU, with batch, tile and work sizes different from the first run;
  r3  the systems of cells_cpu.json on the CPU device (numpy), threads limited by --b-threads.
Every --every seconds, and at the end, the out/ and logs/ directories of each run are packed into
<WORKSPACE>/out/repl-r2.tgz and repl-r3.tgz (written to a temporary name, then renamed), and the progress is written
to <WORKSPACE>/out/chain-status.json.  The top-level key "finished" becomes true only when both chains have exited;
no other key of that file is named "finished" (the copy loop greps for it).
"""
import argparse
import json
import os
import shlex
import signal
import subprocess
import sys
import tarfile
import time

ap = argparse.ArgumentParser()
ap.add_argument('--deadline', type=float, default=1e12)
ap.add_argument('--a-args', default='--batch 768 --tile-mb 512 --work-gb 12')
ap.add_argument('--b-args', default='--device cpu')
ap.add_argument('--a-cells', default='cells.json')
ap.add_argument('--b-cells', default='cells_cpu.json')
ap.add_argument('--a-threads', type=int, default=6)
ap.add_argument('--b-threads', type=int, default=16)
ap.add_argument('--every', type=float, default=120.0)
args = ap.parse_args()

WS = os.environ.get('WORKSPACE', '/workspace')
CODE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WS, 'out')
STATUS = os.path.join(OUT, 'chain-status.json')
RUNS = [('r2', args.a_args, args.a_cells, args.a_threads), ('r3', args.b_args, args.b_cells, args.b_threads)]
T0 = time.time()
procs = {}


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def pack(name):
    d = os.path.join(WS, name)
    tmp = os.path.join(OUT, 'repl-%s.tgz.tmp' % name)
    with tarfile.open(tmp, 'w:gz') as tar:
        for sub in ('out', 'logs'):
            p = os.path.join(d, sub)
            if os.path.isdir(p):
                tar.add(p, arcname='%s/%s' % (name, sub))
    os.replace(tmp, os.path.join(OUT, 'repl-%s.tgz' % name))


def progress(name):
    try:
        with open(os.path.join(WS, name, 'out', 'chain-status.json')) as f:
            st = json.load(f)
    except Exception:   # noqa: BLE001
        return None
    systems = [s for c in st.get('cells', {}).values() for s in c.get('systems', {}).values()]
    return dict(steps=len(st.get('steps', [])), systems_seen=len(systems),
                systems_done=sum(1 for s in systems if s.get('status') == 'done'),
                not_done=sorted(set(s.get('status') for s in systems if s.get('status') not in ('done', 'running'))),
                chain_ended=st.get('ended'), chain_aborted=st.get('aborted'))


def save(finished, aborted=None):
    state = dict(started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(T0)), time=now(), elapsed_s=round(time.time() - T0),
                 runs={n: dict(rc=procs[n].poll() if n in procs else None, args=a, cells=c, threads=t, progress=progress(n))
                       for n, a, c, t in RUNS}, aborted=aborted, finished=finished)
    tmp = STATUS + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, STATUS)


def on_term(sig, frm):
    for p in procs.values():
        if p.poll() is None:
            p.terminate()
    time.sleep(5)
    for n, _a, _c, _t in RUNS:
        pack(n)
    save(True, aborted='signal %d' % sig)
    sys.exit(75)


def main():
    os.makedirs(OUT, exist_ok=True)
    signal.signal(signal.SIGTERM, on_term)
    for name, a, cells, threads in RUNS:
        d = os.path.join(WS, name)
        for sub in ('out', 'logs'):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        link = os.path.join(d, 'in')
        if not os.path.exists(link):
            os.symlink(os.path.join(WS, 'in'), link)
        env = os.environ.copy()
        env['WORKSPACE'] = d
        for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            env[v] = str(threads)
        cmd = ['python3', 'chain_bin.py', '--cells', os.path.join(WS, 'in', cells), '--deadline', str(args.deadline)] + shlex.split(a)
        log = open(os.path.join(d, 'logs', 'chain.out'), 'ab')
        procs[name] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=CODE, env=env)
    save(False)
    while any(p.poll() is None for p in procs.values()):
        time.sleep(args.every)
        for name, _a, _c, _t in RUNS:
            pack(name)
        save(False)
    for name, _a, _c, _t in RUNS:
        pack(name)
    save(True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
