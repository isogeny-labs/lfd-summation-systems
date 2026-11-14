#!/usr/bin/env python3
"""Pod-side wrapper for the m-summand runs (gpu/RUN.md): unpacks the code bundle, then runs several chain_msum.py
chains side by side, each in its own workspace (<WS>/<chain>/{in,out,logs}) with its own thread count, and packs each
chain's out/ and logs/ into <WS>/out/msum-<chain>.tgz every --every seconds and at the end (temporary name, then
rename).  Progress in <WS>/out/chain-status.json; its top-level "finished" becomes true only when every chain has exited.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import tarfile
import time

ap = argparse.ArgumentParser()
ap.add_argument('--plan', default='chains.json')
ap.add_argument('--deadline', type=float, default=1e12)
ap.add_argument('--every', type=float, default=120.0)
args = ap.parse_args()
WS = os.environ.get('WORKSPACE', '/workspace')
OUT = os.path.join(WS, 'out')
CODE = os.path.join(WS, 'code')
STATUS = os.path.join(OUT, 'chain-status.json')
T0 = time.time()
procs = {}


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def pack(name):
    d = os.path.join(WS, name)
    tmp = os.path.join(OUT, 'msum-%s.tgz.tmp' % name)
    with tarfile.open(tmp, 'w:gz') as tar:
        for sub in ('out', 'logs'):
            p = os.path.join(d, sub)
            if os.path.isdir(p):
                tar.add(p, arcname='%s/%s' % (name, sub))
    os.replace(tmp, os.path.join(OUT, 'msum-%s.tgz' % name))


def progress(name):
    try:
        with open(os.path.join(WS, name, 'out', 'chain-status.json')) as f:
            st = json.load(f)
    except Exception:   # noqa: BLE001
        return None
    systems = [s for c in st.get('cells', {}).values() for s in c.get('systems', {}).values()]
    return dict(systems_seen=len(systems), systems_done=sum(1 for s in systems if s.get('status') == 'done'),
                other=sorted(set(str(s.get('status')) for s in systems if s.get('status') not in ('done', 'running'))),
                chain_ended=st.get('ended'), chain_aborted=st.get('aborted'))


def save(chains, finished, aborted=None):
    state = dict(started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(T0)), time=now(), elapsed_s=round(time.time() - T0),
                 chains={c['name']: dict(rc=procs[c['name']].poll() if c['name'] in procs else None, progress=progress(c['name'])) for c in chains},
                 aborted=aborted, finished=finished)
    tmp = STATUS + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, STATUS)


def main():
    os.makedirs(OUT, exist_ok=True)
    bundle = os.path.join(WS, 'in', 'msum-bundle.tgz')
    with tarfile.open(bundle) as tar:
        tar.extractall(CODE, filter='data')
    plan = json.load(open(os.path.join(CODE, 'gpu', args.plan)))
    chains = plan['chains']

    def on_term(sig, frm):
        for p in procs.values():
            if p.poll() is None:
                p.terminate()
        time.sleep(5)
        for c in chains:
            pack(c['name'])
        save(chains, True, aborted='signal %d' % sig)
        sys.exit(75)
    signal.signal(signal.SIGTERM, on_term)
    for c in chains:
        d = os.path.join(WS, c['name'])
        for sub in ('out', 'logs', 'in'):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        env = os.environ.copy()
        env['WORKSPACE'] = d
        for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            env[v] = str(c['threads'])
        cmd = ['python3', os.path.join(CODE, 'gpu', 'chain_msum.py'), '--code-dir', CODE, '--deadline', str(args.deadline)] + c['args']
        log = open(os.path.join(d, 'logs', 'chain.out'), 'ab')
        procs[c['name']] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=CODE, env=env)
    save(chains, False)
    while any(p.poll() is None for p in procs.values()):
        time.sleep(args.every)
        for c in chains:
            pack(c['name'])
        save(chains, False)
    for c in chains:
        pack(c['name'])
    save(chains, True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
