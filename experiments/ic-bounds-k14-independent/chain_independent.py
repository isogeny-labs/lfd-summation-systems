#!/usr/bin/env python3
"""Pod-side chain for PLAN.md: the independent lift at k = 14 (NOT RUN at the target; see PLAN.md).

All files sit in one directory (the pod's /workspace/code): ic_bounds_big.py, ic_bounds_big_run.py,
ic_bounds_cert_check.py, ic_bounds_lift_independent.py, run_ic_lift_independent.py, this file, and
optionally ic_bounds_lift.py + ic_bounds_lift_run.py (the first lift, for a same-basis comparison).
Instances are read from ROOT/in, results go to ROOT/out, checkpoints to ROOT/ckpt, logs to ROOT/logs.

  python3 chain_independent.py --stage small                 # arithmetic + k = 10 checks on the device
  python3 chain_independent.py --stage k14 --tags unit-p65537 [--first-lift]
  CHAIN_DRY=1 python3 chain_independent.py --stage all --tags unit-p65537,unit-p163841   # print the commands

Order per k = 14 instance: closure D = 6 (engine; one checkpoint, at the end of the closure phase),
independent checker on that checkpoint, independent lift, optionally the first lift.  The
checkpoint is never deleted by this script.  A step that fails stops the chain.
"""
import argparse
import json
import os
import subprocess
import sys
import time

DRY = bool(os.environ.get('CHAIN_DRY'))
HERE = os.path.dirname(os.path.abspath(__file__))

SMALL = {     # k = 10, D = 5 -> 6: values recorded in experiments/ic-bounds-k14-independent/validation.json
    'unit-p65537-a1-s10-k10': dict(contains_one=True, dim_low=21700, dim_low_parts=[10878, 10822], dim_VD1=60460, dim_VD=18569),
    'proper-p65537-s12-k10': dict(contains_one=False, dim_low=21686, dim_low_parts=[10871, 10815], dim_VD1=60446, dim_VD=18569),
}
K14 = {       # expected at k = 14: closure record of the earlier runs, and the first lift's records
    'unit-p65537': dict(inst='unit-p65537-a1-s10-k14.json',
                        basis_sha256=['eeac52ef422214ffe21d411fa9b6da488de2b63e20871875dbd94e81cb3d8c38',
                                      'a58c6c7d982c6f3e1144d8a79223e23aeba34bd965f940ffe8a9182defd0fb24'],
                        lift=dict(contains_one=True, dim_low=499178, dim_low_parts=[249824, 249354], dim_VD1=1683218)),
    'unit-p163841': dict(inst='unit-p163841-a3-s21-k14.json',
                         basis_sha256=['e0c4afd72517d56c3b0f78369ce9aaee4a7872b087cdaf383b4c682e5676c170',
                                       '44a89923b4cee5ef2a831892a83a24f0df98a325fbcbafbc70ec0bb77efde905'],
                         lift=dict(contains_one=True, dim_low=499178, dim_low_parts=[249824, 249354], dim_VD1=1683218)),
    'proper-p114689': dict(inst='proper-p114689-s22-k14.json', basis_sha256=None,
                           lift=dict(contains_one=False, dim_low=496808, dim_low_parts=[248639, 248169], dim_VD1=1680848)),
}
COMMON14 = dict(dim_VD=448922, n_special_tops=11038, n_regular_tops=1173002, n_N_weight_D=1406)


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


class Chain:
    def __init__(self, a):
        self.a = a
        self.IN, self.OUT, self.CK, self.LOGS = (os.path.join(a.root, d) for d in ('in', 'out', 'ckpt', 'logs'))
        self.state = dict(started=now(), device=a.device, steps=[], readings={}, finished=False)
        if not DRY:
            for d in (self.OUT, self.CK, self.LOGS):
                os.makedirs(d, exist_ok=True)

    def save(self):
        if DRY:
            return
        tmp = os.path.join(self.OUT, 'chain-independent-status.json.tmp')
        with open(tmp, 'w') as f:
            json.dump(self.state, f, indent=1)
        os.replace(tmp, tmp[:-4])

    def run(self, name, argv, ok=(0,)):
        cmd = [sys.executable] + argv
        if DRY:
            print(name + ': ' + ' '.join(cmd))
            return 0
        step = dict(name=name, cmd=cmd, start=now(), rc=None)
        self.state['steps'].append(step)
        self.save()
        with open(os.path.join(self.LOGS, 'chain-%s.log' % name), 'ab') as log:
            rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=HERE)
        step.update(rc=rc, end=now())
        self.save()
        if rc not in ok:
            self.stop('%s exited with %s' % (name, rc))
        return rc

    def stop(self, why):
        self.state.update(aborted=why, finished=True, ended=now())
        self.save()
        print('ABORT: ' + why, flush=True)
        sys.exit(1)

    def load(self, path):
        if DRY:
            return {}
        with open(path) as f:
            return json.load(f)

    def gpu(self, argv):
        """engine and checker scripts run with the cupy DEVICE pool disabled (as the recorded runs)."""
        w = os.path.join(HERE, 'nopool_run.py')
        return ([w] if self.a.device == 'gpu' and os.path.exists(w) else []) + argv

    def lift(self, name, ck, inst, crec, out, device, extra=()):
        argv = ['run_ic_lift_independent.py', 'lift', '--ckpt', ck, '--instance', inst, '--closure-record', crec,
                '--device', device, '--out', out] + list(extra)
        if device == 'gpu':
            argv.append('--no-pool')
        self.run(name, argv, ok=(0, 3))
        return self.load(out)

    def compare(self, tag, rec, want, what):
        got = {k: rec.get(k) for k in want}
        self.state['readings'].setdefault(tag, {})[what] = dict(got=got, want=want, equal=got == want,
                                                              fixed_point=rec.get('fixed_point'), saturated=rec.get('saturated'),
                                                              seconds=rec.get('seconds'), status=rec.get('status'))
        self.save()
        return DRY or (got == want and rec.get('fixed_point') is True)

    # ---- stage: small
    def small(self):
        dev = self.a.device
        self.run('arith', ['run_ic_lift_independent.py', 'arith', '--device', dev, '--out', os.path.join(self.OUT, 'indep-arith-%s.json' % dev)])
        for tag, want in SMALL.items():
            inst = os.path.join(self.IN, tag + '.json')
            ck, crec = os.path.join(self.CK, tag + '-D5'), os.path.join(self.OUT, tag + '-D5.json')
            self.run('closure-' + tag, self.gpu(['ic_bounds_big_run.py', 'closure', '--instance', inst, '--D', '5', '--device', dev,
                                                 '--work-gb', '4', '--ckpt', ck, '--out', crec]))
            r = self.lift('indep-' + tag, ck, inst, crec, os.path.join(self.OUT, 'indep-lift-%s-D5to6-%s.json' % (tag, dev)), dev)
            if not self.compare(tag, r, want, 'lift_' + dev):
                self.stop('small check failed on %s (%s)' % (tag, dev))
            if dev == 'gpu':       # the same checkpoint with the basis on the device, and on the CPU
                r2 = self.lift('indep-dev-' + tag, ck, inst, crec, os.path.join(self.OUT, 'indep-lift-%s-D5to6-gpu-basisdev.json' % tag),
                               'gpu', ['--basis-on', 'device', '--first-chunk', '256', '--seed', '5'])
                r3 = self.lift('indep-cpu-' + tag, ck, inst, crec, os.path.join(self.OUT, 'indep-lift-%s-D5to6-cpu.json' % tag), 'cpu')
                if not (self.compare(tag, r2, want, 'lift_gpu_basis_on_device') and self.compare(tag, r3, want, 'lift_cpu')):
                    self.stop('small check failed on %s (gpu basis on device / cpu)' % tag)

    # ---- stage: k14
    def k14(self, tag):
        a, I = self.a, K14[tag]
        dev = a.device
        inst = os.path.join(self.IN, I['inst'])
        ck, crec = os.path.join(self.CK, tag + '-k14-D6'), os.path.join(self.OUT, tag + '-k14-D6.json')
        R = self.state['readings'].setdefault(tag, {})
        if not (a.reuse_closure and os.path.exists(crec)):
            # the recorded settings of the second k = 14 run; --ckpt-minutes 100000: one checkpoint, at the end
            self.run('closure-' + tag, self.gpu(['ic_bounds_big_run.py', 'closure', '--instance', inst, '--D', '6', '--device', dev,
                                                 '--dev-gb', '105', '--work-gb', '20', '--tile-mb', '1024', '--batch', '1024', '--threads', '8',
                                                 '--ckpt', ck, '--ckpt-minutes', '100000', '--resume',
                                                 '--status', os.path.join(self.OUT, tag + '-k14-D6.status.json'), '--out', crec]))
        d = self.load(crec)
        R['closure'] = {k: d.get(k) for k in ('status', 'closed', 'contains_one', 'rank', 'profile', 'part_ranks', 'repairs', 'seconds', 'basis_sha256')}
        R['closure']['hashes_equal_recorded'] = None if I['basis_sha256'] is None else d.get('basis_sha256') == I['basis_sha256']
        self.save()
        if not DRY and not (d.get('status') == 'done' and d.get('closed') and d.get('contains_one') is False and d.get('rank') == 448922):
            self.stop('closure of %s is not the closed degree-6 closure of dimension 448,922 without 1' % tag)
        self.run('certcheck-' + tag, self.gpu(['ic_bounds_cert_check.py', '--device', dev, '--shifts', 'all', '--instance', inst, '--ckpt', ck,
                                               '--D', '6', '--engine-record', crec, '--work-gb', '45', '--row-chunk', '4096', '--progress-seconds', '30',
                                               '--status', os.path.join(self.OUT, 'certcheck-%s-k14-D6.status.json' % tag),
                                               '--out', os.path.join(self.OUT, 'certcheck-%s-k14-D6.json' % tag)]))
        c = self.load(os.path.join(self.OUT, 'certcheck-%s-k14-D6.json' % tag))
        R['certcheck'] = {k: c.get(k) for k in ('all_true', 'certificate_valid', 'one_excluded')}
        self.save()
        if not DRY and not c.get('all_true'):
            self.stop('independent checker did not confirm the basis of %s' % tag)
        r = self.lift('indep-lift-' + tag, ck, inst, crec, os.path.join(self.OUT, 'indep-lift-%s-k14-D6to7.json' % tag), dev,
                      ['--max-hours', str(a.max_hours), '--table-gb', str(a.table_gb), '--batch', '2048', '--verbose'])
        want = dict(I['lift'], **COMMON14)
        R['verdict'] = 'confirm' if self.compare(tag, r, want, 'independent_lift') else \
            ('refute' if r.get('fixed_point') else 'undecided')
        if not DRY and r.get('status') == 'refused':
            R['verdict'] = 'undecided (refused: %s)' % r.get('reason')
        self.save()
        if a.first_lift:
            self.run('first-lift-' + tag, self.gpu(['ic_bounds_lift_run.py', 'lift', '--ckpt', ck, '--instance', inst, '--device', dev,
                                                    '--basis-on', 'host', '--dev-gb', '130', '--work-gb', '16',
                                                    '--out', os.path.join(self.OUT, 'lift-%s-k14-D6to7-again.json' % tag)]), ok=(0, 1, 2, 75))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', default='all', choices=['small', 'k14', 'all'])
    ap.add_argument('--tags', default='unit-p65537')
    ap.add_argument('--device', default='gpu', choices=['gpu', 'cpu'])
    ap.add_argument('--root', default='/workspace')
    ap.add_argument('--max-hours', type=float, default=6.0, help='deadline of the independent lift')
    ap.add_argument('--table-gb', type=float, default=16.0)
    ap.add_argument('--first-lift', action='store_true', help='also rerun the first lift on the same checkpoint')
    ap.add_argument('--reuse-closure', action='store_true', help='skip the closure if its record exists')
    a = ap.parse_args()
    C = Chain(a)
    if a.stage in ('small', 'all'):
        C.small()
    if a.stage in ('k14', 'all'):
        for tag in a.tags.split(','):
            C.k14(tag)
    C.state.update(finished=True, ended=now())
    C.save()
    print(json.dumps(C.state['readings'], indent=1))


if __name__ == '__main__':
    main()
