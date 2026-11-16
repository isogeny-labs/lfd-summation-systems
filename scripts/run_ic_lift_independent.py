#!/usr/bin/env python3
"""Driver of the independent degree lift (csd/ic_bounds_lift_independent.py).

  lift  --ckpt DIR --instance FILE [--closure-record FILE] [--device cpu|gpu] --out FILE
        reads the stored basis of V_D (engine checkpoint), rebuilds the generator from (p, a, b, X, k),
        and decides whether 1 is in V_{D+1}; reports dim(V_{D+1} n M_<=D) in total and per sigma-part.
  plan  --n-free A,B --n-special S --regular R --ranks A,B      memory estimate from counts
  arith [--device cpu|gpu]                                      exact-arithmetic self-test

Exit codes: 0 result written (fixed point reached), 3 result written but not a fixed point,
4 refused (an input check failed; the reason is in the record).
The only repository module imported is csd/ic_bounds_lift_independent.py (found next to this file
or in ../csd).
"""
import argparse
import hashlib
import json
import os
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, '..', 'csd')):
    if os.path.exists(os.path.join(d, 'ic_bounds_lift_independent.py')):
        sys.path.insert(0, d)
        break
import numpy as np                                    # noqa: E402
import ic_bounds_lift_independent as IL               # noqa: E402


def sha(path):
    return IL._sha256_file(path)


def get_xp(device):
    if device == 'gpu':
        import cupy
        return cupy
    return np


def log(**kw):
    print(json.dumps(kw), flush=True)


def cmd_lift(a):
    xp = get_xp(a.device)
    t0 = time.time()
    rec = dict(command='lift', ckpt=a.ckpt, instance_file=os.path.basename(a.instance), device=a.device,
               source_sha256=dict(lift=sha(IL.__file__), driver=sha(__file__), instance=sha(a.instance)))
    try:
        inst = json.load(open(a.instance))
        k, p, X = inst['k'], inst['prime'], inst['target'][0]
        g = IL.s3_generator(p, inst['a'], inst['b'], X, k)
        stored = {tuple(int(x) for x in key.split(',')): int(v) for key, v in inst['generator'].items()}
        if stored != g:
            raise IL.Refuse('the generator in the instance file is not S_3 of its (p, a, b, X, k)')
        ck = IL.read_checkpoint(a.ckpt, verify_files=not a.no_verify_files)
        if (ck['k'], ck['p']) != (k, p):
            raise IL.Refuse(f'checkpoint is for k = {ck["k"]}, p = {ck["p"]}')
        rec.update(basis=dict(path=ck['path'], files_verified=ck['files_verified'],
                              manifest_files_id=ck['manifest_files_id'], engine_sha256=ck['engine_sha256']))
        log(event='loaded', k=k, D=ck['D'], p=p, ranks=[len(pv) for pv, _ in ck['parts']], sec=round(time.time() - t0, 1))
        lift = IL.Lift(k, ck['D'], p, ck['parts'], g, xp=xp, include_low=a.include_low, batch=a.batch,
                       table_gb=a.table_gb, first_chunk=a.first_chunk, seed=a.seed, log=log)
        if a.closure_record:
            cr = json.load(open(a.closure_record))
            mine = [IL.canonical_sha(pv, B, nc, p) for (pv, B), nc in zip(ck['parts'], lift.ncol)]
            rec['basis']['canonical_sha256'] = mine
            rec['basis']['closure_record'] = dict(file=os.path.basename(a.closure_record), closed=cr.get('closed'),
                                                  basis_sha256=cr.get('basis_sha256'))
            if cr.get('basis_sha256') != mine:
                raise IL.Refuse('canonical hashes of the stored basis differ from the closure record')
            if not cr.get('closed'):
                raise IL.Refuse('the closure record does not say closed: true')
        log(event='tables', n_N=lift.nN, n_special=lift.nSp, n_regular=len(lift.regular), n_T=lift.nT,
            plan=IL.plan(lift.nN, lift.nSp, len(lift.regular), [len(pv) for pv, _ in ck['parts']],
                         [B.shape[1] for _, B in ck['parts']], a.batch, a.table_gb), sec=round(time.time() - t0, 1))
        rec.update(lift.run(closed_check=a.closed_check))
        rec['status'] = 'done'
        code = 0 if rec['fixed_point'] else 3
    except IL.Refuse as e:
        rec.update(status='refused', reason=str(e))
        code = 4
    rec['seconds_total'] = round(time.time() - t0, 1)
    rec['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    if a.device == 'gpu':
        rec['device_peak_mb'] = int(xp.get_default_memory_pool().total_bytes() // 2 ** 20)
    log(**{key: v for key, v in rec.items() if key != 'counts'}, event='result')
    if a.out:
        tmp = a.out + '.part'
        with open(tmp, 'w') as f:
            json.dump(rec, f, indent=1)
            f.write('\n')
        os.replace(tmp, a.out)
    sys.exit(code)


def cmd_plan(a):
    fr = [int(x) for x in a.n_free.split(',')]
    rk = [int(x) for x in a.ranks.split(',')]
    print(json.dumps(IL.plan(sum(fr) + a.n_single, a.n_special, a.regular, rk, fr, a.batch, a.table_gb), indent=1))


def cmd_arith(a):
    """modp and mm against Python integers, at the bounds and across the inner-dimension split."""
    xp = get_xp(a.device)
    rng = np.random.RandomState(7)
    out = dict(device=a.device, cases=[])
    ok = True
    for p in (12289, 65537, 163841, 16777213):
        lim = 2 ** 53 - 2 * p
        vals = [0, 1, p - 1, p, p + 1, lim, lim - 1, -lim, -lim + 1, -1, -p, -p - 1]
        vals += [q * p + d for q in (lim // p - 1, lim // p // 2, 3) for d in (-1, 0, 1) if abs(q * p + d) <= lim]
        vals += [int(v) for v in rng.randint(-2 ** 52, 2 ** 52, size=2000)]
        X = xp.asarray(np.array(vals, dtype=np.float64))
        good = np.array_equal(IL.to_host(IL.modp(X, p, xp)), np.array([v % p for v in vals], dtype=np.float64))
        mi = IL.max_inner(p)
        A = np.full((3, min(mi, 5000) + 7), p - 1, dtype=np.float64)
        B = np.full((A.shape[1], 2), p - 1, dtype=np.float64)
        A[1] = rng.randint(0, p, size=A.shape[1])
        want = np.array([[sum(int(x) * int(y) for x, y in zip(A[i], B[:, j])) % p for j in range(2)] for i in range(3)])
        res = []
        for m in (None, 5, 7, A.shape[1] - 3):
            M = IL.mm(xp.asarray(A), xp.asarray(B), p, xp, mi=m)
            res.append(bool(np.array_equal(IL.to_host(M), want)))
        ok = ok and good and all(res)
        out['cases'].append(dict(p=p, max_inner=mi, modp_ok=bool(good), mm_ok=res))
    out['passed'] = bool(ok)
    print(json.dumps(out, indent=1))
    if a.out:
        json.dump(out, open(a.out, 'w'), indent=1)
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['lift', 'plan', 'arith'])
    ap.add_argument('--ckpt')
    ap.add_argument('--instance')
    ap.add_argument('--closure-record', help='engine record of the closure that wrote the checkpoint: its basis_sha256 must match')
    ap.add_argument('--out')
    ap.add_argument('--device', default='cpu', choices=['cpu', 'gpu'])
    ap.add_argument('--no-verify-files', action='store_true')
    ap.add_argument('--include-low', action='store_true', help='also shift the rows of pivot weight < D')
    ap.add_argument('--closed-check', type=int, default=2, help='random combinations in the closedness test (0: off)')
    ap.add_argument('--batch', type=int, default=2048)
    ap.add_argument('--table-gb', type=float, default=2.0)
    ap.add_argument('--first-chunk', type=int, default=2048)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--n-free'); ap.add_argument('--ranks'); ap.add_argument('--n-single', type=int, default=0)
    ap.add_argument('--n-special', type=int); ap.add_argument('--regular', type=int)
    a = ap.parse_args()
    {'lift': cmd_lift, 'plan': cmd_plan, 'arith': cmd_arith}[a.cmd](a)


if __name__ == '__main__':
    main()
