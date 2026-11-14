#!/usr/bin/env python3
"""Driver for csd/ic_bounds_bin_big.py: exact closure V_D of binary (characteristic 2) Weil-descent systems.

Subcommands
  closure   one closure at degree D from an instance JSON (optional checkpoint at the end, optional resume =
            load the checkpoint and run the certificate sweep)
  ladder    D = 0, 1, 2, ... up to the certificate degree (or unit, or --max-degree): fall profile, LFD, local
            certificate; one record per degree is written as soon as it exists
  plan      predicted peak memory and operation estimate for (2k variables, D)
  selftest  F_2 kernels against integer arithmetic + the embedded small systems on cpu and on the device
  instance  (LOCAL: imports the repository builders) build an instance JSON from a system spec
  zcheck    re-enumerate Z from the instance description with the driver's own enumerator (pod-side check)

The engine file is found next to this script (pod layout) or in ../csd/ (repository layout).
Exit codes: 0 ok; 1 failure; 3 sweep failed with --no-repair (negative control).
"""
import argparse
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
for cand in (HERE / 'ic_bounds_bin_big.py', HERE.parent / 'csd' / 'ic_bounds_bin_big.py'):
    if cand.exists():
        sys.path.insert(0, str(cand.parent))
        break
else:  # pragma: no cover
    raise SystemExit('ic_bounds_bin_big.py not found next to this script or in ../csd/')
for cand in (HERE / 'ic_bounds_big.py', HERE.parent / 'csd' / 'ic_bounds_big.py'):
    if cand.exists():
        sys.path.insert(0, str(cand.parent))
        break

import numpy as np  # noqa: E402
import ic_bounds_bin_big as EB  # noqa: E402


# =============================================================================== instance building (independent of the repository)

def system_hash_bin(r, equations):
    """same definition as csd.ic_curve_sweep.prepare_binary / csd.ic_bounds_binary.system_hash('binary', ...)"""
    eqs = sorted(set(tuple(sorted(e)) for e in equations if e))
    system = dict(kind='binary', variables=r, equations=eqs)
    return hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


class GF:
    """F_(2^n) = F_2[x]/(modulus), carry-less arithmetic on Python ints and, vectorised, on uint64 arrays."""

    def __init__(self, n, modulus):
        self.n, self.mod = n, modulus

    def mul(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a >> self.n & 1:
                a ^= self.mod
        return out

    def sq(self, a):
        return self.mul(a, a)

    def vmul(self, A, B):
        A = A.astype(np.uint64).copy()
        B = np.asarray(B, dtype=np.uint64)
        out = np.zeros_like(A)
        one = np.uint64(1)
        for i in range(self.n):
            out ^= np.where((B >> np.uint64(i)) & one, A, np.uint64(0))
            A = A << one
            A = np.where((A >> np.uint64(self.n)) & one, A ^ np.uint64(self.mod), A)
        return out

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.sq(a)
            e >>= 1
        return r

    def inv(self, a):
        return self.pow(a, (1 << self.n) - 2)

    def sqrt(self, a):
        for _ in range(self.n - 1):
            a = self.sq(a)
        return a

    def trace(self, a):
        t, x = 0, a
        for _ in range(self.n):
            t ^= x
            x = self.sq(x)
        return t


def span_list(basis):
    out = [0]
    for v in basis:
        out += [u ^ v for u in out]
    return out


def _solve_affine(vecs, rhs):
    """all c (bit mask over vecs) with XOR of vecs[j] for j in c == rhs: (particular, kernel basis) or None."""
    piv = {}
    kern = []
    for j, v in enumerate(vecs):
        comb = 1 << j
        while v:
            p = v.bit_length() - 1
            if p in piv:
                pv, pc = piv[p]
                v ^= pv
                comb ^= pc
            else:
                piv[p] = (v, comb)
                break
        if not v:
            kern.append(comb)
    comb, v = 0, rhs
    while v:
        p = v.bit_length() - 1
        if p not in piv:
            return None
        pv, pc = piv[p]
        v ^= pv
        comb ^= pc
    return comb, kern


def enumerate_Z(n, modulus, basis, b, xR):
    """Complete zero set Z of S_3(x1, x2, xR) = (x1x2 + x1xR + x2xR)^2 + x1x2xR + b on V x V, V = span(basis),
    as the integers c1 | c2 << k (c = coordinate vector).  Pure field arithmetic: S_3 = alpha x2^2 + beta x2 + gamma
    with alpha = (x1 + xR)^2, beta = x1 xR, gamma = beta^2 + b is F_2-affine in x2; solved on the coordinates of x2
    for each x1.  Uses nothing from the descended Boolean equations and nothing from the repository."""
    F = GF(n, modulus)
    k = len(basis)
    V = np.array(span_list(basis), dtype=np.uint64)
    xr = np.full(len(V), xR, dtype=np.uint64)
    s = V ^ np.uint64(xR)
    alpha = F.vmul(s, s)
    beta = F.vmul(V, xr)
    gamma = F.vmul(beta, beta) ^ np.uint64(b)
    imgs = []
    for v in basis:
        imgs.append(F.vmul(alpha, np.full(len(V), F.sq(v), dtype=np.uint64)) ^ F.vmul(beta, np.full(len(V), v, dtype=np.uint64)))
    imgs = np.stack(imgs)
    out = []
    for c1 in range(len(V)):
        sol = _solve_affine([int(imgs[j, c1]) for j in range(k)], int(gamma[c1]))
        if sol is None:
            continue
        part, kern = sol
        pts = [part]
        for kv in kern:
            pts += [q ^ kv for q in pts]
        out.extend(c1 | (c << k) for c in pts)
    return sorted(out)


def structure_of(equations, r):
    """rank, m_eff, r_eff, const of the span of the equations over F_2 (same as experiments/ic-bounds-binary-law/
    structure.py::analyse) and e = r_eff - m_eff."""
    mons = sorted({m for e in equations for m in e}, key=lambda m: (m.bit_count(), m))
    pos = {m: i for i, m in enumerate(mons)}

    def rank_rows(rows):
        basis = {}
        for v in rows:
            while v:
                p = v.bit_length() - 1
                if p in basis:
                    v ^= basis[p]
                else:
                    basis[p] = v
                    break
        return len(basis), basis
    vecs = []
    for e in equations:
        v = 0
        for m in e:
            v ^= 1 << pos[m]
        vecs.append(v)
    rank, basis = rank_rows(vecs)
    low_mask = sum(1 << pos[m] for m in mons if m.bit_count() <= 1)
    low = [v for p, v in basis.items() if (1 << p) & low_mask]
    lin = []
    for v in low:
        w = 0
        for m in mons:
            if (v >> pos[m]) & 1 and m.bit_count() == 1:
                w ^= m
        lin.append(w)
    lin_rank, _ = rank_rows(lin)
    const = False
    if 0 in pos:
        red = 1 << pos[0]
        for p in sorted(basis, reverse=True):
            if (red >> p) & 1:
                red ^= basis[p]
        const = red == 0
    k = r // 2
    A = (1 << k) - 1
    bil = all(m.bit_count() < 2 or (m.bit_count() == 2 and (m & A).bit_count() == 1) for m in mons)
    m_eff = rank - len(low)
    return dict(rank=rank, m_eff=m_eff, L=len(low), lin_rank=lin_rank, const=const, r_eff=r - lin_rank,
                bilinear=bil, equations=len(equations), e=(r - lin_rank) - m_eff)


def _seed(*parts):
    return int.from_bytes(hashlib.sha256(':'.join(map(str, parts)).encode()).digest()[:8], 'big')


def build_instance(spec, with_z=True):
    """LOCAL: imports the repository builders.  spec: name, n, k, base ('random' | 'invariant' | 'trace-zero' |
    'subfield' | 'given'), curve ('koblitz0' | 'koblitz1' | 'ordinary' | 'given'), seed, target ('order-two' |
    'random-rational' | 'given'), optional index (which invariant space), basis/a/b/target for 'given'."""
    root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(root))
    from csd import ic_dreg as I
    from csd import ic_curve_sweep as S
    from csd.ic_arbitrary_bases import canonical_basis
    from csd.ic_research_campaign import field
    n, k = spec['n'], spec['k']
    F = field(n)
    rng = random.Random(_seed(spec.get('seed', 0), spec['name']))
    cv = spec['curve']
    if cv == 'koblitz0':
        a, b = 0, 1
    elif cv == 'koblitz1':
        a, b = 1, 1
    elif cv == 'ordinary':
        a, b = rng.randrange(1 << n), rng.randrange(2, 1 << n)
    else:
        a, b = spec['a'], spec['b']
    E = I.BinaryCurve(F, a, b)
    bk = spec['base']
    if bk == 'random':
        basis = canonical_basis(S.random_nonstable_subspace(F, k, rng, []))
    elif bk == 'invariant':
        cands = S.invariant_candidates(F, k, spec.get('limit', 4), rng)
        basis = canonical_basis(cands[spec.get('index', 0) % len(cands)])
    elif bk == 'trace-zero':
        assert k == n - 1
        basis = canonical_basis([F.sq(1 << i) ^ (1 << i) for i in range(n)])
        assert len(basis) == k
    elif bk == 'subfield':
        assert n == 2 * k
        imgs = [F.frob(1 << i, k) ^ (1 << i) for i in range(n)]
        piv, kern = {}, []
        for i in range(n):
            v, pre = imgs[i], 1 << i
            while v:
                p = v.bit_length() - 1
                if p in piv:
                    v ^= piv[p][0]
                    pre ^= piv[p][1]
                else:
                    piv[p] = (v, pre)
                    break
            if not v:
                kern.append(pre)
        basis = canonical_basis(kern)
    else:
        basis = canonical_basis(spec['basis'])
    assert len(basis) == k
    tk = spec['target']
    if tk == 'order-two':
        target = (0, F.sqrt(b))
    elif tk == 'random-rational':
        target = S.random_rational_target(E, rng, 1 << n, exclude_x={0})
    else:
        target = tuple(spec['target_point'])
    assert E.on_curve(target)
    eqs = sorted(set(tuple(sorted(row)) for row in I.binary_descent(E, [basis, basis], target[0]) if row))
    r = 2 * k
    st = structure_of([list(e) for e in eqs], r)
    inst = dict(format='ic-bin-instance/1', name=spec['name'], spec=spec, n=n, modulus=F.modulus, k=k, variables=r,
                basis=basis, curve=dict(a=a, b=b, family='koblitz' if a in (0, 1) and b == 1 else 'ordinary'),
                target=list(target), frobenius_invariant=bool(I.is_frobenius_invariant(F, basis)),
                trace_zero=all(F.trace(v) == 0 for v in basis),
                scaled_subfield=bool(spec['base'] == 'subfield'), structure=st,
                equations=[list(e) for e in eqs], system_hash=system_hash_bin(r, eqs))
    inst['class'] = 'const' if st['const'] else ('order-two' if target[0] == 0 else 'affine')
    if with_z:
        Z = enumerate_Z(n, F.modulus, basis, b, target[0])
        inst['points'] = Z
        inst['n_points'] = len(Z)
        inst['points_sha256'] = hashlib.sha256(json.dumps(Z).encode()).hexdigest()
        inst['unit_class'] = bool(len(Z) == 0 and not st['const'])
    return inst


# =============================================================================== common

def load_instance(path):
    with open(path) as f:
        d = json.load(f)
    r = d['variables']
    d['equations'] = [[int(m) for m in e] for e in d['equations']]
    d['points'] = [int(p) for p in d.get('points', [])]
    assert system_hash_bin(r, d['equations']) == d['system_hash'], 'system hash mismatch'
    if 'points_sha256' in d:
        assert hashlib.sha256(json.dumps(d['points']).encode()).hexdigest() == d['points_sha256'], 'points hash mismatch'
    return d


def check_points(inst):
    """every point vanishes on every equation (vectorised) and the points are distinct."""
    pts = np.array(inst['points'], dtype=np.uint64)
    assert len(set(inst['points'])) == len(pts)
    ok = True
    for e in inst['equations']:
        acc = np.zeros(len(pts), dtype=np.uint8)
        for m in e:
            mm = np.uint64(m)
            acc ^= ((pts & mm) == mm).astype(np.uint8)
        ok &= not acc.any()
    return bool(ok)


def write_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=str)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def engine_kw(a):
    return dict(device=a.device, batch=a.batch, tile_rows=a.tile_rows, tile_mb=a.tile_mb, work_gb=a.work_gb,
                col_block=a.col_block, compact_frac=a.compact_frac, tiles_on=a.tiles_on,
                withhold=a.withhold, repair=not a.no_repair, status=a.status, quiet=False)


def cmd_closure(a):
    inst = load_instance(a.instance)
    r = inst['variables']
    eng = EB.Engine(inst['equations'], r, a.D, points=inst['points'], **engine_kw(a))
    if a.resume and a.ckpt and os.path.exists(os.path.join(a.ckpt, 'manifest.json')):
        rec = eng.resume_certificate(a.ckpt)
        rec['resumed'] = True
    else:
        rec = eng.run()
        if a.ckpt and rec['closed']:
            eng.save_checkpoint(a.ckpt, extra=dict(system_hash=inst['system_hash']))
    rec.update(system_hash=inst['system_hash'], n_points=len(inst['points']), instance=os.path.basename(a.instance))
    if a.out:
        write_json(a.out, rec)
    print(json.dumps({k: rec[k] for k in ('D', 'columns', 'rank', 'codim', 'closed', 'unit', 'repairs', 'sweeps', 'seconds',
                                        'rref_sha256', 'peak_rss_mb')}), flush=True)
    return 0 if rec['closed'] else 3


def cmd_ladder(a):
    inst = load_instance(a.instance)
    r = inst['variables']
    t0 = time.time()
    ok_pts = check_points(inst) if inst['points'] else True
    assert ok_pts, 'a listed point does not vanish on the equations'
    recs = []

    def log(rec):
        recs.append(rec)
        if a.out:
            write_json('%s.D%d.json' % (a.out, rec['D']), dict(rec, system_hash=inst['system_hash']))
        print(json.dumps({k: rec[k] for k in ('D', 'columns', 'rank', 'closed', 'unit', 'repairs', 'seconds')}), flush=True)
    prof = EB.local_profile_bin(inst['equations'], r, inst['points'], a.max_degree, engine_kw(a), log=log)
    out = dict(name=inst['name'], system_hash=inst['system_hash'], n=inst['n'], k=inst['k'], n_points=len(inst['points']),
               points_sha256=inst.get('points_sha256'), cls=inst.get('class'), engine=EB.ENGINE_VERSION,
               engine_sha256=EB.ENGINE_SHA256, device=a.device, status=prof['status'], exact_lfd=prof['exact_lfd'],
               certificate_degree=prof['certificate_degree'], strict_falls=prof['strict_falls'],
               quotient_ranks=prof['quotient_ranks'], local_checks=prof['local_checks'],
               profiles=[r_['profile'] for r_ in prof['per_degree']], rref_sha256=[r_['rref_sha256'] for r_ in prof['per_degree']],
               columns=[r_['columns'] for r_ in prof['per_degree']], per_degree_seconds=[r_['seconds'] for r_ in prof['per_degree']],
               per_degree_flops=[r_['flops_update'] + r_['flops_residual'] for r_ in prof['per_degree']],
               closed_all=all(r_['closed'] for r_ in prof['per_degree']), repairs=sum(r_['repairs'] for r_ in prof['per_degree']),
               seconds=round(time.time() - t0, 2), peak_rss_mb=EB.maxrss_mb(),
               device_peak_bytes=max((r_['device_peak_bytes'] for r_ in prof['per_degree']), default=0))
    if a.out:
        write_json(a.out, out)
    print(json.dumps({k: out[k] for k in ('name', 'status', 'exact_lfd', 'certificate_degree', 'strict_falls', 'seconds')}), flush=True)
    return 0 if out['status'] in ('exact-unit', 'exact-proper-local') and out['closed_all'] else 3


def cmd_plan(a):
    out = []
    for k in a.k:
        for D in a.D_list:
            out.append(EB.plan(2 * k, D, a.device or 'gpu', batch=a.batch, work_gb=a.work_gb, tile_mb=a.tile_mb,
                               compact_frac=a.compact_frac, dev_gb=a.dev_gb))
    print(json.dumps(out, indent=1))
    if a.out:
        write_json(a.out, out)
    return 0


def cmd_selftest(a):
    dev = a.device or 'cpu'
    res = dict(device=dev, started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), engine=EB.ENGINE_VERSION,
               engine_sha256=EB.ENGINE_SHA256)
    res['kernels'] = EB.selftest_kernels(dev)
    ok = all(res['kernels'].values())
    cases = EB.SELFTEST_CASES
    base = EB.selftest_systems('cpu', cases)
    res['cpu_cases'] = base
    ok &= all(c['ok'] for c in base)
    if dev != 'cpu':
        dd = EB.selftest_systems(dev, cases, tiny=a.tiny)
        res['device_cases'] = dd
        ok &= all(c['ok'] for c in dd)
        ok &= all(x['hashes'] == y['hashes'] and x['dims'] == y['dims'] for x, y in zip(base, dd))
        EB.use_device(dev)
        res['device_memory'] = dict(zip(('used_bytes', 'total_bytes'), EB.DV.mem_info()))
    res['ok'] = bool(ok)
    if a.out:
        write_json(a.out, res)
    print(json.dumps({k: res[k] for k in ('device', 'ok')} | dict(cases=len(cases))), flush=True)
    return 0 if ok else 1


def cmd_instance(a):
    specs = json.load(open(a.spec))
    if isinstance(specs, dict):
        specs = [specs]
    os.makedirs(a.outdir, exist_ok=True)
    for s in specs:
        t0 = time.time()
        inst = build_instance(s)
        assert check_points(inst) if inst['points'] else True
        path = os.path.join(a.outdir, s['name'] + '.json')
        with open(path, 'w') as f:
            json.dump(inst, f, separators=(',', ':'))
        print('%s k=%d n=%d class=%s |Z|=%d e=%d bytes=%d %.1fs' % (s['name'], inst['k'], inst['n'], inst['class'], inst['n_points'],
              inst['structure']['e'], os.path.getsize(path), time.time() - t0), flush=True)
    return 0


def cmd_zcheck(a):
    inst = load_instance(a.instance)
    t0 = time.time()
    Z = enumerate_Z(inst['n'], inst['modulus'], inst['basis'], inst['curve']['b'], inst['target'][0])
    ok = Z == sorted(inst['points']) and check_points(inst)
    print(json.dumps(dict(name=inst['name'], z_equal=ok, n_points=len(Z), seconds=round(time.time() - t0, 1))))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['closure', 'ladder', 'plan', 'selftest', 'instance', 'zcheck'])
    ap.add_argument('--instance')
    ap.add_argument('--spec')
    ap.add_argument('--outdir')
    ap.add_argument('--out')
    ap.add_argument('--k', type=int, nargs='+')
    ap.add_argument('--D', type=int, default=5)
    ap.add_argument('--D-list', dest='D_list', type=int, nargs='+', default=[5])
    ap.add_argument('--max-degree', type=int, default=7)
    ap.add_argument('--device', default=None, choices=['cpu', 'gpu', 'fake-gpu'])
    ap.add_argument('--batch', type=int, default=1024)
    ap.add_argument('--tile-rows', type=int, default=None)
    ap.add_argument('--tile-mb', type=float, default=1024.0)
    ap.add_argument('--work-gb', type=float, default=20.0)
    ap.add_argument('--col-block', type=int, default=None)
    ap.add_argument('--compact-frac', type=float, default=0.25)
    ap.add_argument('--tiles-on', default='device', choices=['device', 'host'])
    ap.add_argument('--dev-gb', type=float, default=141.0)
    ap.add_argument('--withhold', type=int, default=None, help='TEST HOOK: row index whose shifts are skipped in the closure (negative control)')
    ap.add_argument('--no-repair', action='store_true')
    ap.add_argument('--ckpt')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--status')
    ap.add_argument('--tiny', action='store_true')
    ap.add_argument('--threads', type=int, default=None)
    a = ap.parse_args()
    if a.device is None and a.cmd != 'plan':
        a.device = 'cpu'
    return dict(closure=cmd_closure, ladder=cmd_ladder, plan=cmd_plan, selftest=cmd_selftest, instance=cmd_instance,
                zcheck=cmd_zcheck)[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
