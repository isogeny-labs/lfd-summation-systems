#!/usr/bin/env python3
"""Pod-side chain for the k = 14 runs lost on 2026-11-11 (see ../RELAUNCH.md).

Runs in /workspace/code on the pod and needs no connection to the launching machine:
  self-tests, then for each instance: closure D = 6 with the basis kept, lift D = 6 -> 7,
  independent checker on the kept basis.
Every step logs to /workspace/logs/chain-<step>.log; progress is in /workspace/out/chain-status.json.
`--dry` prints the commands and runs nothing.  WORKSPACE overrides /workspace (tests).
"""
import json
import os
import shutil
import subprocess
import sys
import time

WS = os.environ.get('WORKSPACE', '/workspace')
OUT, CK, IN, LOGS = (os.path.join(WS, d) for d in ('out', 'ckpt', 'in', 'logs'))
DRY = '--dry' in sys.argv
STATUS = os.path.join(OUT, 'chain-status.json')
FULL7 = 1683218          # |M_<=7| at k = 14
INSTANCES = [
    # tag, instance file, solutions, what the lift must satisfy
    dict(tag='proper-p114689-k14', inst='proper-p114689-s22-k14.json', points=2370),
    dict(tag='unit-p163841-k14', inst='unit-p163841-a3-s21-k14.json', points=0),
]
state = dict(started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), steps=[], readings={}, finished=False)


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save():
    if DRY:
        return
    tmp = STATUS + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, STATUS)


def run(name, argv, nopool=True):
    cmd = ['python3'] + (['nopool_run.py'] if nopool else []) + argv
    if DRY:
        print(name + ': ' + ' '.join(cmd))
        return 0
    step = dict(name=name, cmd=cmd, start=now(), rc=None)
    state['steps'].append(step)
    save()
    with open(os.path.join(LOGS, 'chain-%s.log' % name), 'ab') as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    step.update(rc=rc, end=now())
    save()
    return rc


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:   # noqa: BLE001
        return None


def one(I):
    tag, inst = I['tag'], os.path.join(IN, I['inst'])
    ck = os.path.join(CK, tag + '-D6')
    rec = os.path.join(OUT, tag + '-D6.json')
    R = state['readings'].setdefault(tag, {})
    rc = run(tag + '-D6', ['ic_bounds_big_run.py', 'closure', '--instance', inst, '--D', '6', '--device', 'gpu',
                           '--dev-gb', '105', '--work-gb', '20', '--tile-mb', '1024', '--batch', '1024', '--threads', '8',
                           '--ckpt', ck, '--ckpt-minutes', '30', '--resume',
                           '--status', os.path.join(OUT, tag + '-D6.status.json'), '--out', rec])
    if DRY:
        d = dict(status='done', closed=True, contains_one=False, rank=448922)
    else:
        d = load(rec) or {}
    R['closure'] = {k: d.get(k) for k in ('status', 'closed', 'contains_one', 'rank', 'profile', 'part_ranks', 'repairs', 'seconds')}
    R['closure']['rc'] = rc
    save()
    if not (d.get('status') == 'done' and d.get('closed') and d.get('contains_one') is False):
        R['lift'] = 'skipped: no closed degree-6 closure without 1'
        save()
        return
    lrec = os.path.join(OUT, 'lift-%s-D6to7.json' % tag)
    for attempt, work in enumerate(('16', '8')):
        rc = run('lift-%s-try%d' % (tag, attempt + 1),
                 ['ic_bounds_lift_run.py', 'lift', '--ckpt', ck, '--instance', inst, '--device', 'gpu', '--basis-on', 'host',
                  '--dev-gb', '130', '--work-gb', work, '--lift-ckpt', os.path.join(CK, '%s-D6-lift%d' % (tag, attempt + 1)),
                  '--ckpt-minutes', '20', '--resume',
                  '--status', os.path.join(OUT, 'lift-%s-D6to7.status.json' % tag), '--out', lrec])
        if rc == 0 or DRY:
            break
    L = {} if DRY else (load(lrec) or {})
    bound = FULL7 - I['points']
    dim = L.get('dim_VD1')
    R['lift'] = dict(rc=rc, status=L.get('status'), fixed_point=L.get('fixed_point'), dim_VD1=dim, contains_one=L.get('contains_one'),
                     rounds=L.get('rounds'), guard_hits=L.get('guard_hits'), seconds=L.get('seconds'), ideal_bound=bound,
                     within_ideal=(dim is not None and dim <= bound), equals_ideal_part=(dim == bound))
    save()
    rc = run('certcheck-' + tag, ['ic_bounds_cert_check.py', '--device', 'gpu', '--shifts', 'all', '--instance', inst, '--ckpt', ck,
                                  '--D', '6', '--engine-record', rec, '--work-gb', '45', '--row-chunk', '4096',
                                  '--progress-seconds', '30', '--status', os.path.join(OUT, 'certcheck-%s-D6.status.json' % tag),
                                  '--out', os.path.join(OUT, 'certcheck-%s-D6.json' % tag)])
    C = {} if DRY else (load(os.path.join(OUT, 'certcheck-%s-D6.json' % tag)) or {})
    R['certcheck'] = dict(rc=rc, all_true=C.get('all_true'), certificate_valid=C.get('certificate_valid'), one_excluded=C.get('one_excluded'))
    save()
    if not DRY and L.get('status') == 'done' and L.get('fixed_point'):
        for d_ in os.listdir(CK):          # free the disk for the next instance; kept if the lift did not finish
            if d_.startswith(tag):
                shutil.rmtree(os.path.join(CK, d_), ignore_errors=True)


def main():
    if not DRY:
        for d in (OUT, CK, LOGS):
            os.makedirs(d, exist_ok=True)
    a = run('selftest-engine', ['ic_bounds_big_run.py', 'selftest', '--device', 'gpu', '--out', os.path.join(OUT, 'selftest-gpu.json')])
    b = run('selftest-lift', ['ic_bounds_lift_run.py', 'selftest', '--device', 'gpu', '--out', os.path.join(OUT, 'lift-selftest-gpu.json')])
    if a or b:
        state['aborted'] = 'a self-test failed (engine rc %s, lift rc %s); nothing else was run' % (a, b)
    else:
        for I in INSTANCES:
            one(I)
    state.update(finished=True, ended=now())
    save()


if __name__ == '__main__':
    main()
