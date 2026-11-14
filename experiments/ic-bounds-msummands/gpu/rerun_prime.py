#!/usr/bin/env python3
"""Second job on the m-summand pod (gpu/RUN.md, "Prime rerun"): the prime chains of the first job failed at once
because the code bundle lacked csd/_elim.c.  This adds the missing C sources to /workspace/code/csd (never
overwriting an existing file), then runs the prime chains of --plan in /workspace/p2/<chain> with the same
chain_msum.py, packs their records into /workspace/p2/out/msum-<chain>.tgz, and writes /workspace/p2/out/chain-status.json.
Its top-level "finished" becomes true only when its own chains have ended AND the first job's status
(/workspace/out/chain-status.json) says finished, so that the copy loop (which watches this file) terminates the pod
only after both jobs are done."""
import argparse
import json
import os
import subprocess
import sys
import tarfile
import time

ap = argparse.ArgumentParser()
ap.add_argument('--plan', default='chains-prime.json')
ap.add_argument('--deadline', type=float, default=1e12)
ap.add_argument('--every', type=float, default=120.0)
args = ap.parse_args()
WS0 = '/workspace'
CODE = os.path.join(WS0, 'code')
WS = os.path.join(WS0, 'p2')
OUT = os.path.join(WS, 'out')
STATUS = os.path.join(OUT, 'chain-status.json')
T0 = time.time()
procs = {}


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def first_finished():
    try:
        with open(os.path.join(WS0, 'out', 'chain-status.json')) as f:
            return bool(json.load(f).get('finished'))
    except Exception:   # noqa: BLE001
        return False


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
                chain_ended=st.get('ended'))


def save(chains, finished, own_done):
    state = dict(started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(T0)), time=now(), elapsed_s=round(time.time() - T0),
                 chains={c['name']: dict(rc=procs[c['name']].poll() if c['name'] in procs else None, progress=progress(c['name'])) for c in chains},
                 own_chains_done=own_done, first_job_finished=first_finished(), finished=finished)
    tmp = STATUS + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, STATUS)


def main():
    os.makedirs(OUT, exist_ok=True)
    with tarfile.open(os.path.join(WS0, 'in', 'msum-bundle2.tgz')) as tar:
        members = [m for m in tar.getmembers() if m.name.startswith('csd/') and m.name.endswith('.c')]
        for m in members:
            if not os.path.exists(os.path.join(CODE, m.name)):
                tar.extract(m, CODE, filter='data')
        for name in ('gpu/' + args.plan,):
            tar.extract(tar.getmember(name), CODE, filter='data')
    chains = json.load(open(os.path.join(CODE, 'gpu', args.plan)))['chains']
    for c in chains:
        d = os.path.join(WS, c['name'])
        for sub in ('out', 'logs', 'in'):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        env = os.environ.copy()
        env['WORKSPACE'] = d
        for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            env[v] = str(c['threads'])
        cmd = ['python3', os.path.join(CODE, 'gpu', 'chain_msum.py'), '--code-dir', CODE, '--deadline', str(args.deadline)] + c['args']
        procs[c['name']] = subprocess.Popen(cmd, stdout=open(os.path.join(d, 'logs', 'chain.out'), 'ab'), stderr=subprocess.STDOUT,
                                            stdin=subprocess.DEVNULL, cwd=CODE, env=env)
    save(chains, False, False)
    while True:
        own = all(p.poll() is not None for p in procs.values())
        for c in chains:
            pack(c['name'])
        if own and first_finished():
            save(chains, True, True)
            return 0
        save(chains, False, own)
        time.sleep(args.every)


if __name__ == '__main__':
    sys.exit(main())
