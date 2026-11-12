#!/usr/bin/env python3
"""Test harness for scripts/ic_bounds_cert_check.py (the independent certificate checker).

This file MAY import the engines (it uses them only to PRODUCE test inputs); the checker itself imports none.
Groups (``--groups``, default all):
  arith      modp / product exactness against python integers, near the bound, on cpu and fake-gpu
  spec       the checker's layout, shifts, projection and S_3 expansion against the reference engine; S_3 vanishes
             at x(P1 +- P2); W cap M_<D = span of the low rows (dense computation on a small case)
  positive   checkpoints written by scripts/ic_bounds_big_run.py closure (end of the closure phase)
  reference  bases of csd/ic_bounds_d5_verify.py converted to the checkpoint format
  negative   corrupted inputs, each must be rejected (and for the right reason), and a harmless corruption
  devices    cpu vs fake-gpu, chunked inner dimension, tiny blocks: identical verdicts
Light limits: every subprocess runs under RLIMIT_AS 1.5 GB and one BLAS thread.  Scratch under --scratch,
deleted at the end unless --keep.  Records go to experiments/ic-bounds-k14/validation/certcheck-<group>.json.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import random
import resource
import shutil
import signal
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
CHECK = os.path.join(HERE, 'ic_bounds_cert_check.py')
spec = importlib.util.spec_from_file_location('ic_bounds_cert_check', CHECK)
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
from csd import ic_bounds_d5_verify as REF  # noqa: E402  (test inputs and cross-checks only)

DRIVER = os.path.join(HERE, 'ic_bounds_big_run.py')
VALDIR = os.path.join(ROOT, 'experiments', 'ic-bounds-k14', 'validation')
INSTDIR = os.path.join(ROOT, 'experiments', 'ic-bounds-k14', 'instances')
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
SC = None


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (1_500_000_000, 1_500_000_000))


def sh(cmd, timeout=1800):
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV, preexec_fn=limit, timeout=timeout)
    return r.returncode, r.stdout, r.stderr


def jl(path):
    return json.load(open(path))


def make_instance(name, k, p, a, b, seed, unit):
    out = os.path.join(SC, f'{name}.json')
    if not os.path.exists(out):
        cmd = [sys.executable, DRIVER, 'instance', '--k', str(k), '--prime', str(p), '--a', str(a), '--b', str(b), '--seed', str(seed), '--out', out]
        if unit:
            cmd.append('--unit-nonresidue')
        rc, so, se = sh(cmd)
        assert rc == 0, se[-500:]
    return out


def run_engine(inst, D, ck, rec, extra=()):
    cmd = [sys.executable, DRIVER, 'closure', '--instance', inst, '--D', str(D), '--ckpt', ck, '--out', rec, '--work-gb', '0.25',
           '--tile-mb', '4', '--batch', '128', '--no-eval', *extra]
    rc, so, se = sh(cmd)
    return rc, so, se


def run_check(inst, ck, out, extra=()):
    cmd = [sys.executable, CHECK, '--instance', inst, '--ckpt', ck, '--out', out, '--progress-seconds', '1e9', '--work-gb', '0.1', *extra]
    t = time.time()
    rc, so, se = sh(cmd)
    V = jl(out) if os.path.exists(out) else None
    return rc, V, round(time.time() - t, 2), se[-400:]


def summ(V):
    if V is None:
        return None
    cl = V.get('closure') or {}
    return dict(format_ok=V.get('format_ok'), g_in_W=V.get('g_in_W'), closed=V.get('closed'), one_excluded=V.get('one_excluded'),
                certificate_valid=V.get('certificate_valid'), all_true=V.get('all_true'), dimension=(V.get('dimension') or {}).get('W'),
                part_ranks=(V.get('dimension') or {}).get('parts'), profile=(V.get('profile') or {}).get('pivots_weight_le_j'),
                basis_sha256=V.get('basis_sha256'), rows_checked=cl.get('rows_checked'), triples_checked=cl.get('triples_checked'),
                failures=cl.get('failures'), failures_by_direction=cl.get('failures_by_direction'), shift_mode=cl.get('shift_mode'),
                generator_multiples_in_W=V.get('generator_multiples_in_W'), rref_leading_column_ok=V.get('rref_leading_column_ok'),
                reasons=V.get('reasons'), engine_record_equal=(V.get('engine_record_comparison') or {}).get('all_equal'),
                seconds=(V.get('resources') or {}).get('seconds'), maxrss_mb=(V.get('resources') or {}).get('maxrss_mb'))


def save_group(name, rec):
    os.makedirs(VALDIR, exist_ok=True)
    rec['checker_sha256'] = C.file_sha256(CHECK)
    rec['harness_sha256'] = C.file_sha256(__file__)
    path = os.path.join(VALDIR, f'certcheck-{name}.json')
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
        f.write('\n')
    print(f'[{name}] passed={rec["passed"]} -> {os.path.relpath(path, ROOT)}', flush=True)
    return rec['passed']


# ------------------------------------------------------------------------- checkpoint manipulation
def load_parts(ck):
    ck = C.locate_ckpt(ck)
    man = jl(os.path.join(ck, 'manifest.json'))
    parts = []
    for q, nt in enumerate(man['ntiles']):
        piv = np.load(os.path.join(ck, f'p{q}_host.npz'))['piv'].astype(np.int64)
        B = np.concatenate([np.load(os.path.join(ck, f'p{q}_t{ti:05d}.npy')) for ti in range(nt)]).astype(np.int64)
        parts.append((piv, B))
    return man, parts


def save_parts(root, man, parts, cs=None, name='ck-00000001', state=None):
    """rebuild a checkpoint dir from (piv, B) parts; columns of the part from the layout in man['config']."""
    cfg = man['config']
    L = cs or C.Layout(cfg['k'], cfg['D'], cfg['split'])
    full = []
    for q, (piv, B) in enumerate(parts):
        c = L.parts[q].c
        ispiv = np.zeros(c, dtype=bool)
        ispiv[piv] = True
        free = np.flatnonzero(~ispiv)
        rows = np.zeros((len(piv), c), dtype=np.int64)
        rows[np.arange(len(piv)), piv] = 1
        rows[:, free] = B
        full.append((piv, rows))
    st = dict(man['state'])
    st.update(state or {})
    return C.write_checkpoint(root, full, cfg, tile_rows=64, state=st, name=name)


def copy_tree(src, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


# ------------------------------------------------------------------------- group: arith
def group_arith():
    rec = dict(group='arith', cases=[])
    ok = True
    for dev in ('cpu', 'fake-gpu'):
        r = C.selftest_arith(C.Dev(dev))
        rec['cases'].append(dict(device=dev, **r))
        ok &= r['passed']
    # independent: python-int reference for modp over every residue class boundary near +-(2^53 - 2p), random
    # entries of the product blocks near the bound, and the chunked path with a tiny forced inner dimension
    rng = np.random.default_rng(99)
    extra = []
    for p in (12289, 65537, 114689, 16777213):
        lim = (1 << 53) - 2 * p
        vals = [lim - i for i in range(0, 40)] + [-(lim - i) for i in range(0, 40)] + [int(x) for x in rng.integers(-lim, lim, 20000)]
        X = np.array(vals, dtype=np.float64)
        assert all(int(x) == v for x, v in zip(X, vals))
        got = C.modp(C.Dev('cpu'), X, p).astype(np.int64)
        e1 = all(int(g) == v % p for g, v in zip(got, vals))
        mi = C.max_inner(p)
        # the largest sum a block can produce: (p-1)^2 mi + (p-1) <= lim
        e2 = (p - 1) ** 2 * mi + (p - 1) <= lim and (p - 1) ** 2 * (mi + 1) + (p - 1) > lim
        for forced in (None, 5):
            C.MI_OVERRIDE = forced
            m = 37
            A = rng.integers(0, p, size=(4, m), dtype=np.int64)
            Bm = rng.integers(0, p, size=(m, 6), dtype=np.int64)
            A[0, :] = p - 1
            Bm[:, 0] = p - 1
            g = C.gemm_mod(C.Dev('cpu'), A.astype(np.float64), Bm.astype(np.float64), p).astype(np.int64)
            want = np.array([[sum(int(A[i, t]) * int(Bm[t, j]) for t in range(m)) % p for j in range(6)] for i in range(4)])
            e3 = bool(np.array_equal(g, want))
            extra.append(dict(p=p, forced_inner=forced, max_inner=mi, modp_exact=e1, bound_tight=bool(e2), gemm_exact=e3))
            ok &= e1 and e3 and bool(e2)
        C.MI_OVERRIDE = None
    rec['cases'].append(dict(kind='python-int cross-check', rows=extra))
    rec['passed'] = bool(ok)
    return save_group('arith', rec)


# ------------------------------------------------------------------------- group: spec
def group_spec():
    rec = dict(group='spec', cases=[])
    ok = True
    # S_3: own expansion = reference expansion; vanishes at x(P1 + P2) and x(P1 - P2)
    for (p, a, b, X) in ((12289, 1, 0, 9362), (15361, 4021, 9917, 4491), (65537, 4021, 9917, 1424)):
        e1 = C.s3_expand(p, a, b, X) == REF.s3_coeffs(p, a, b, X)
        rng = random.Random(7)
        pts = []
        while len(pts) < 60:
            x = rng.randrange(p)
            rhs = (x ** 3 + a * x + b) % p
            if pow(rhs, (p - 1) // 2, p) == 1 and rhs:
                y = next(t for t in range(1, p) if t * t % p == rhs) if p < 20000 else None
                if y is None:
                    # tonelli not needed: p = 65537 = 1 mod 4 -> use pow trick with a non residue search
                    y = _sqrt(rhs, p)
                pts.append((x, y))

        def add(P1, P2):
            (x1, y1), (x2, y2) = P1, P2
            if x1 == x2:
                if (y1 + y2) % p == 0:
                    return None
                lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
            else:
                lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
            x3 = (lam * lam - x1 - x2) % p
            return (x3, (lam * (x1 - x3) - y1) % p)

        bad = 0
        n = 0
        for i in range(0, 60, 2):
            P1, P2 = pts[i], pts[i + 1]
            for Q in (add(P1, P2), add(P1, (P2[0], -P2[1] % p))):
                if Q is None:
                    continue
                n += 1
                S = C.s3_expand(p, a, b, Q[0])
                val = sum(c * pow(P1[0], e1_, p) * pow(P2[0], e2_, p) for (e1_, e2_), c in S.items()) % p
                bad += val != 0
        rec['cases'].append(dict(kind='S_3', p=p, a=a, b=b, X=X, equals_reference=bool(e1), vanishing_tests=n, nonvanishing=int(bad)))
        ok &= bool(e1) and bad == 0 and n > 20
    # layout equality with the reference
    for (k, D) in ((6, 4), (7, 5), (8, 4)):
        Lr = REF.Layout(k, D, True)
        Lc = C.Layout(k, D, True)
        same = (np.array_equal(Lr.a1, Lc.a1) and np.array_equal(Lr.a2, Lc.a2) and np.array_equal(Lr.w, Lc.w) and Lr.lowc == Lc.lowc)
        same &= all(np.array_equal(Lr.shift[ij], Lc.shift[ij]) for ij in Lr.shift)
        same &= all(np.array_equal(Lr.parts[i]['reps'], Lc.parts[i].reps) for i in (0, 1))
        rng = np.random.default_rng(3)
        H = rng.integers(0, 12289, size=(5, Lc.c)).astype(np.float64)
        pr = Lr.project(H.copy(), 12289)
        pc = Lc.project(H.astype(np.int64), 12289)
        proj_ok = all(np.array_equal(np.mod(a_, 12289), b_.astype(np.float64)) for a_, b_ in zip(pr, pc))
        # lift o project o lift is 2 lift on the symmetric part, 2 lift on the anti part (up to the engine's conventions)
        rec['cases'].append(dict(kind='layout', k=k, D=D, columns=Lc.c, layout_equal=bool(same), projection_equal=bool(proj_ok)))
        ok &= bool(same and proj_ok)
    # dense check: W cap M_<D has dimension = number of low rows (k=6, D=4)
    p, k, D = 12289, 6, 4
    g = REF.chain_generator(p, 1, 0, 9362, k)
    cl = REF.Closure(g, k, p, D)
    cl.run()
    Lc = C.Layout(k, D, False)
    # full-coordinate basis of W from the reference parts
    rows = []
    for pi, P in enumerate(cl.parts):
        rows.append(cl.L.lift(pi, P.rows(np.arange(P.rank))))
    full = np.mod(np.vstack(rows), p).astype(np.int64)
    high = np.arange(Lc.low0)
    # dim(W cap M_<D) = rank(W) - rank(W restricted to the weight >= D columns)
    r_all = len(C.rref_rows(full, p)[0])
    r_high = len(C.rref_rows(full[:, high], p)[0])
    nlow = sum(len(cl.low_rows(i)) for i in range(2))
    rec['cases'].append(dict(kind='W cap M_<D is spanned by the low rows', k=k, D=D, rank_W=r_all, rank_on_high_columns=r_high,
                             dim_W_cap_M_lt_D=r_all - r_high, low_rows=nlow, equal=(r_all - r_high == nlow)))
    ok &= (r_all - r_high == nlow)
    # dense, full-coordinate closure test (no parts, no projection, no layers): independent of the checker's machinery
    def dense_closed(parts_pb, Lfull):
        rowsf = []
        for pi, (piv, B) in enumerate(parts_pb):
            c_ = cl.L.parts[pi]['c']
            ispiv = np.zeros(c_, dtype=bool)
            ispiv[piv] = True
            fr = np.flatnonzero(~ispiv)
            Rw = np.zeros((len(piv), c_), dtype=np.int64)
            Rw[np.arange(len(piv)), piv] = 1
            Rw[:, fr] = B
            rowsf.append(np.mod(cl.L.lift(pi, Rw.astype(np.float64)), p).astype(np.int64))
        full_ = np.vstack(rowsf)
        pv_, rr_ = C.rref_rows(full_, p)
        low_ = [i for i, ld in enumerate(pv_) if Lfull.w[ld] < D]
        bad_ = 0
        for i in low_:
            for ij in Lfull.all_shifts:
                v = np.zeros(Lfull.c, dtype=np.int64)
                v[Lfull.shift[ij]] = rr_[i][Lfull.low0:]
                for ld, row in zip(pv_, rr_):
                    if v[ld]:
                        v = (v - v[ld] * row) % p
                bad_ += bool(v.any())
        return bad_, len(low_)

    import io as _io
    ref_parts = [(P.piv.astype(np.int64), P.B.astype(np.int64)) for P in cl.parts]
    Lf = C.Layout(k, D, False)
    bad0, nlow0 = dense_closed(ref_parts, Lf)
    c0 = cl.parts[0]
    lowi = [i for i in range(c0.rank) if c0.w[c0.piv[i]] < D]
    i0 = lowi[len(lowi) // 2]
    jj = [j for j, f in enumerate(c0.free) if f > c0.piv[i0] and c0.w[f] < D][0]
    corrupt = [(P.piv.astype(np.int64), P.B.astype(np.int64).copy()) for P in cl.parts]
    corrupt[0][1][i0, jj] = (corrupt[0][1][i0, jj] + 1) % p
    bad1, nlow1 = dense_closed(corrupt, Lf)
    verdicts = []
    for label, pb in (('untouched', ref_parts), ('low-coefficient +1', corrupt)):
        root = os.path.join(SC, 'dense-' + label.split()[0])
        if os.path.exists(root):
            shutil.rmtree(root)
        Ls = C.Layout(k, D, True)
        full_rows = []
        for q_, (piv, B) in enumerate(pb):
            c_ = Ls.parts[q_].c
            ispiv = np.zeros(c_, dtype=bool)
            ispiv[piv] = True
            Rw = np.zeros((len(piv), c_), dtype=np.int64)
            Rw[np.arange(len(piv)), piv] = 1
            Rw[:, np.flatnonzero(~ispiv)] = B
            full_rows.append((piv, Rw))
        C.write_checkpoint(root, full_rows, dict(k=k, p=p, D=D, split=True, gens='all', source_id='dense'), tile_rows=50)
        gi = REF.chain_generator(p, 1, 0, 9362, k)
        ip = os.path.join(SC, 'dense-inst.json')
        json.dump(dict(k=k, prime=p, a=1, b=0, target=[9362, 0], generator={f'{e1},{e2}': v for (e1, e2), v in gi.items()}), open(ip, 'w'))
        ns = argparse.Namespace(instance=ip, ckpt=root, D=None, device='cpu', shifts='all', out=None, status=None, engine_record=None,
                                gen_mult_weight=1, col_block=None, row_chunk=64, work_gb=0.05, progress_seconds=1e9, fail_fast=False)
        old = sys.stdout
        sys.stdout = _io.StringIO()
        try:
            Vd = C.verify(ns)
        finally:
            sys.stdout = old
        verdicts.append(dict(case=label, checker_closed=Vd['closed'], checker_failures=Vd['closure']['failures']))
        shutil.rmtree(root)
    agree = (bad0 == 0 and verdicts[0]['checker_closed'] is True and bad1 > 0 and verdicts[1]['checker_closed'] is False)
    rec['cases'].append(dict(kind='dense full-coordinate closure test agrees with the checker', k=k, D=D, low_rows_dense=nlow0,
                             dense_failing_row_shift_pairs=dict(untouched=bad0, corrupted=bad1), checker=verdicts, agree=bool(agree)))
    ok &= bool(agree)
    rec['passed'] = bool(ok)
    return save_group('spec', rec)


def _sqrt(a, p):
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


# ------------------------------------------------------------------------- group: positive
POS = [
    # name, k, p, a, b, seed, unit, D, expect_one_excluded
    ('unit-p12289-k7-D4', 7, 12289, 1, 0, 10, True, 4, True),
    ('proper-p15361-k7-D4', 7, 15361, 4021, 9917, 12, False, 4, True),
    ('unit-p12289-k8-D4', 8, 12289, 1, 0, 10, True, 4, True),
    ('proper-p65537-k8-D4', 8, 65537, 4021, 9917, 12, False, 4, True),
    ('unit-p65537-k8-D4', 8, 65537, 1, 0, 10, True, 4, True),
    ('proper-p15361-k7-D5', 7, 15361, 4021, 9917, 12, False, 5, True),
    ('unit-p12289-k6-D5', 6, 12289, 1, 0, 10, True, 5, False),          # 1 in V_D
    ('unit-p65537-k10-D5', 10, 65537, 1, 0, 10, True, 5, True),
    ('proper-p65537-k10-D5', 10, 65537, 4021, 9917, 12, False, 5, True),
]


def produce(name, k, p, a, b, seed, unit, D, extra=()):
    inst = make_instance(f'{name.rsplit("-D", 1)[0]}', k, p, a, b, seed, unit)
    ck = os.path.join(SC, f'ck-{name}')
    recp = os.path.join(SC, f'rec-{name}.json')
    if os.path.exists(ck):
        shutil.rmtree(ck)
    t = time.time()
    rc, so, se = run_engine(inst, D, ck, recp, extra)
    return inst, ck, recp, rc, round(time.time() - t, 1), se[-300:]


def group_positive(only=None):
    rec = dict(group='positive', cases=[])
    ok = True
    for (name, k, p, a, b, seed, unit, D, exp1) in POS:
        if only and not any(o in name for o in only):
            continue
        inst, ck, recp, rc, tsec, err = produce(name, k, p, a, b, seed, unit, D)
        case = dict(name=name, k=k, p=p, D=D, engine_rc=rc, engine_seconds=tsec)
        if rc not in (0, 2) or not os.path.exists(recp):
            case.update(error=err, passed=False)
            rec['cases'].append(case)
            ok = False
            continue
        er = jl(recp)
        case['engine'] = dict(rank=er['rank'], profile=er['profile'], part_ranks=er['part_ranks'], contains_one=er['contains_one'],
                              basis_sha256=er['basis_sha256'], closed=er['closed'])
        res = {}
        for mode in ('all', 'block1'):
            outp = os.path.join(SC, f'v-{name}-{mode}.json')
            rc2, V, tc, se = run_check(inst, ck, outp, ['--engine-record', recp, '--shifts', mode, '--gen-mult-weight', '1'])
            res[mode] = dict(rc=rc2, wall_seconds=tc, **(summ(V) or {}))
        good = True
        for mode, r in res.items():
            good &= bool(r.get('certificate_valid') and r.get('format_ok') and r.get('closed') and r.get('g_in_W'))
            good &= (r.get('one_excluded') is exp1)
            good &= bool(r.get('engine_record_equal'))
            good &= (r.get('dimension') == er['rank'] and r.get('profile') == er['profile'] and r.get('basis_sha256') == er['basis_sha256'])
            good &= bool(r.get('generator_multiples_in_W'))
        case.update(checker=res, expect_one_excluded=exp1, passed=bool(good))
        rec['cases'].append(case)
        ok &= good
        print('  positive', name, 'ok' if good else 'FAIL', flush=True)
        if k >= 10:
            shutil.rmtree(ck, ignore_errors=True)
    rec['passed'] = bool(ok)
    return save_group('positive', rec)


# ------------------------------------------------------------------------- group: reference
def group_reference():
    rec = dict(group='reference', cases=[])
    ok = True
    for (name, k, p, a, b, seed, unit, D, exp1) in [x for x in POS if x[7] <= 5 and x[1] <= 7]:
        inst = make_instance(name.rsplit('-D', 1)[0], k, p, a, b, seed, unit)
        I = jl(inst)
        g = REF.chain_generator(p, I['a'] % p, I['b'] % p, I['target'][0], k)
        cl = REF.Closure(g, k, p, D, split=True, gens='all')
        cl.run()
        parts = []
        for P in cl.parts:
            parts.append((P.piv.astype(np.int64), P.rows(np.arange(P.rank)).astype(np.int64)))
        ck = os.path.join(SC, f'refck-{name}')
        if os.path.exists(ck):
            shutil.rmtree(ck)
        cfg = dict(k=k, p=p, D=D, split=True, gens='all', source_id='reference-engine')
        C.write_checkpoint(ck, parts, cfg, tile_rows=50)
        refsha = [P.sha() for P in cl.parts]
        outp = os.path.join(SC, f'v-ref-{name}.json')
        fake = os.path.join(SC, f'refrec-{name}.json')
        prof, _ = cl.profile()
        json.dump(dict(rank=cl.rank(), profile=prof, part_ranks=[P.rank for P in cl.parts], contains_one=cl.contains_one(),
                       basis_sha256=refsha, columns=cl.L.c), open(fake, 'w'))
        rc, V, tc, se = run_check(inst, ck, outp, ['--engine-record', fake])
        s = summ(V)
        good = bool(s and s['certificate_valid'] and s['engine_record_equal'] and s['basis_sha256'] == refsha and s['one_excluded'] is exp1)
        rec['cases'].append(dict(name=name, reference_rank=cl.rank(), reference_profile=prof, reference_sha256=refsha, checker=s,
                                 reference_sweep_full=cl.sweep_full(), passed=good))
        ok &= good
        print('  reference', name, 'ok' if good else 'FAIL', flush=True)
    rec['passed'] = bool(ok)
    return save_group('reference', rec)


# ------------------------------------------------------------------------- group: negative
def group_negative():
    rec = dict(group='negative', cases=[])
    ok = True
    name, k, p, a, b, seed, unit, D = 'unit-p12289-k8-D4', 8, 12289, 1, 0, 10, True, 4
    inst, ck, recp, rc, tsec, err = produce(name, k, p, a, b, seed, unit, D)
    man, parts = load_parts(ck)
    L = C.Layout(k, D, True)
    cfg = man['config']

    def check_parts(label, parts2, instp=inst, extra=(), state=None, expect=None):
        root = os.path.join(SC, f'neg-{label}')
        if os.path.exists(root):
            shutil.rmtree(root)
        save_parts(root, man, parts2, cs=L, state=state)
        res = {}
        for mode in ('all', 'block1'):
            rc2, V, tc, se = run_check(instp, root, os.path.join(SC, f'vneg-{label}-{mode}.json'), ['--shifts', mode, *extra])
            res[mode] = summ(V)
        shutil.rmtree(root)
        return res

    def lowrows(q):
        piv, B = parts[q]
        return [i for i in range(len(piv)) if L.parts[q].w[piv[i]] < D]

    def highrows(q):
        piv, B = parts[q]
        return [i for i in range(len(piv)) if L.parts[q].w[piv[i]] >= D]

    def freecols(q):
        piv, _ = parts[q]
        ispiv = np.zeros(L.parts[q].c, dtype=bool)
        ispiv[piv] = True
        return np.flatnonzero(~ispiv)

    def case(label, res, want, note):
        nonlocal ok
        good = all(want(r) for r in res.values())
        rec['cases'].append(dict(name=label, note=note, results=res, passed=bool(good)))
        ok &= good
        print('  negative', label, 'ok' if good else 'FAIL', flush=True)

    # sanity: the untouched basis passes
    res = check_parts('control', parts)
    case('untouched-control', res, lambda r: r['certificate_valid'] is True, 'unmodified basis from the engine is accepted')
    # 1. one coefficient changed in a low row
    q = 0
    piv, B = parts[q]
    fc = freecols(q)
    lows = lowrows(q)
    i = lows[len(lows) // 2]
    jcol = [j for j, f in enumerate(fc) if f > piv[i] and L.parts[q].w[f] < D][0]
    B2 = B.copy()
    B2[i, jcol] = (B2[i, jcol] + 1) % p
    res = check_parts('low-coefficient', [(piv, B2), parts[1]])
    case('low-row-coefficient-changed', res, lambda r: r['certificate_valid'] is False and (r['closed'] is False or r['g_in_W'] is False),
         f'part 0, row {i} (pivot weight {int(L.parts[q].w[piv[i]])}), free column {int(fc[jcol])}: +1 mod p')
    # 2a. one coefficient changed in a high row hit by a low shift: closure is the only fact affected
    tried = []
    found = None
    for r_ in highrows(0)[:60]:
        B3 = B.copy()
        B3[r_, 0] = (B3[r_, 0] + 1) % p
        root = os.path.join(SC, 'neg-hi-try')
        if os.path.exists(root):
            shutil.rmtree(root)
        save_parts(root, man, [(piv, B3), parts[1]], cs=L)
        rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-hi-try.json'), ['--fail-fast'])
        shutil.rmtree(root)
        tried.append(dict(row=int(r_), pivot=int(piv[r_]), closed=V['closed'], g_in_W=V['g_in_W'], certificate_valid=V['certificate_valid']))
        if V['closed'] is False and V['g_in_W'] is True:
            found = r_
            break
    B3 = B.copy()
    B3[found, 0] = (B3[found, 0] + 1) % p
    res = check_parts('high-coefficient', [(piv, B3), parts[1]])
    case('high-row-coefficient-changed-closure-only', res,
         lambda r: r['certificate_valid'] is False and r['closed'] is False and r['g_in_W'] is True and r['format_ok'] is True,
         f'part 0, row {found} (pivot column {int(piv[found])}, weight {int(L.parts[0].w[piv[found]])}), free column {int(fc[0])}: +1 mod p; '
         f'g and m*g only involve weights <= 3 < D so they stay in W; only the shift of some low row lands on this pivot. '
         f'rows tried before this one: {len(tried) - 1}')
    # 2b. the same corruption on a high row that no shifted low row ever touches is harmless: W' is still a closed containing space
    harmless = None
    htried = 0
    for r_ in highrows(0)[:200]:
        B3 = B.copy()
        B3[r_, 0] = (B3[r_, 0] + 1) % p
        root = os.path.join(SC, 'neg-hi-try')
        if os.path.exists(root):
            shutil.rmtree(root)
        save_parts(root, man, [(piv, B3), parts[1]], cs=L)
        rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-hi-try.json'), ['--fail-fast'])
        shutil.rmtree(root)
        htried += 1
        if V['certificate_valid']:
            harmless = r_
            break
    if harmless is None:
        rec['cases'].append(dict(name='high-row-coefficient-harmless', note=f'no harmless high row among the first {htried}: every high row of part 0 is hit', passed=True))
    else:
        B3 = B.copy()
        B3[harmless, 0] = (B3[harmless, 0] + 1) % p
        res = check_parts('high-harmless', [(piv, B3), parts[1]])
        case('high-row-coefficient-harmless', res, lambda r: r['certificate_valid'] is True,
             f'part 0, row {harmless} (pivot column {int(piv[harmless])}): changing it yields another space W\' that still contains g and is closed, '
             'so a correct checker accepts it (a different, still valid certificate); it is not a soundness failure')
    # 3. a low row removed (pivot moved to the free columns, the zero column inserted: the format stays valid)
    det = []
    first_closure_only = None
    for q_ in (0, 1):
        pv, Bq = parts[q_]
        for i_ in lowrows(q_)[:80]:
            keep = np.arange(len(pv)) != i_
            fcq = freecols(q_)
            newfree = np.sort(np.concatenate([fcq, [pv[i_]]]))
            pos = int(np.searchsorted(newfree, pv[i_]))
            Bn = np.insert(Bq[keep], pos, 0, axis=1)       # other rows are 0 on a pivot column (RREF)
            newp = [(pv[keep], Bn) if q2 == q_ else parts[q2] for q2 in (0, 1)]
            root = os.path.join(SC, 'neg-rm-try')
            if os.path.exists(root):
                shutil.rmtree(root)
            save_parts(root, man, newp, cs=L)
            rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-rm-try.json'), ['--fail-fast'])
            shutil.rmtree(root)
            det.append((q_, i_, V['format_ok'], V['closed'], V['g_in_W'], V['certificate_valid']))
            if V['format_ok'] and V['closed'] is False and V['g_in_W'] is True and first_closure_only is None:
                first_closure_only = (q_, i_, newp)
    n_rej = sum(1 for d in det if d[5] is False)
    rec['cases'].append(dict(name='low-row-removed-all-tested', tested=len(det), rejected=n_rej, format_invalid=sum(1 for d in det if d[2] is False),
                             rejected_by_closure_only=sum(1 for d in det if d[3] is False and d[4] is True),
                             rejected_by_g_only=sum(1 for d in det if d[3] is True and d[4] is False),
                             rejected_by_both=sum(1 for d in det if d[3] is False and d[4] is False),
                             accepted=[dict(part=d[0], row=d[1]) for d in det if d[5]],
                             note='each low row of the first 80 of each part removed in turn; an accepted removal means the removed row was redundant '
                                  '(the smaller space is still closed and contains g), which is a valid smaller certificate',
                             passed=bool(det and n_rej > 0)))
    ok &= bool(det and n_rej > 0)
    if first_closure_only:
        q_, i_, newp = first_closure_only
        res = check_parts('removed-low', newp)
        case('low-row-removed-closure-only', res, lambda r: r['certificate_valid'] is False and r['closed'] is False and r['g_in_W'] is True and r['format_ok'] is True,
             f'part {q_}, row {i_} removed; format valid, g still in W, rejected by the closure check')
    # 4. a basis from an interrupted closure
    inst2 = make_instance('unit-p12289-k8', k, p, a, b, seed, unit)
    ckint = os.path.join(SC, 'ck-interrupted')
    if os.path.exists(ckint):
        shutil.rmtree(ckint)
    rc_, so, se = run_engine(inst2, D, ckint, os.path.join(SC, 'rec-int.json'), ['--ckpt-minutes', '0.0002', '--ckpt-slots', '200', '--batch', '2'])
    cks = sorted(d for d in os.listdir(ckint) if d.startswith('ck-'))
    sel = []
    for d in cks:
        m = jl(os.path.join(ckint, d, 'manifest.json'))
        sel.append((d, m['reason'], m['state']['phase'], m['state']['closure_complete']))
    mid = [d for d, r_, ph, cc in sel if ph == 'closure' and not cc]
    ires = []
    allres = []
    finalrec = jl(os.path.join(SC, 'rec-int.json')) if os.path.exists(os.path.join(SC, 'rec-int.json')) else None
    for d in mid:
        rc2, V, tc, se = run_check(inst2, os.path.join(ckint, d), os.path.join(SC, 'v-int-one.json'))
        allres.append((d, summ(V)))
    for d in ((mid[0], mid[len(mid) // 2], mid[-1]) if mid else ()):
        ires.append(dict(checkpoint=d, state=[x for x in sel if x[0] == d][0][1:], **dict(allres)[d]))
    rejected = [d for d, r_ in allres if r_['certificate_valid'] is False]
    accepted = [d for d, r_ in allres if r_['certificate_valid'] is True]
    # a basis taken inside the closure phase may already be the final space (the last pending rows were redundant): accepted only then
    acc_equal_final = all(dict(allres)[d]['basis_sha256'] == finalrec['basis_sha256'] for d in accepted) if finalrec else None
    good = bool(ires) and all(r['format_ok'] for r in ires) and ires[0]['certificate_valid'] is False and ires[1]['certificate_valid'] is False \
        and len(rejected) > 0 and bool(acc_equal_final is not False) and all(r_['format_ok'] for _, r_ in allres)
    rec['cases'].append(dict(name='interrupted-closure', note=f'{len(cks)} checkpoints kept ({len(mid)} taken inside the closure phase); first, middle, last shown; all {len(mid)} checked',
                             inside_closure_phase=len(mid), rejected=len(rejected), accepted=len(accepted),
                             accepted_checkpoints_equal_final_basis=acc_equal_final, results=ires, passed=good))
    ok &= good
    print('  negative interrupted', 'ok' if good else 'FAIL', flush=True)
    # 4b. the same, taken at the very end of the closure (closure_done) with one low row's shifts withheld by the engine's own test hook
    ckw = os.path.join(SC, 'ck-withheld')
    if os.path.exists(ckw):
        shutil.rmtree(ckw)
    rc_, so, se = run_engine(inst2, D, ckw, os.path.join(SC, 'rec-wh.json'), ['--gens', 'g', '--withhold', '0:3', '--no-repair'])
    wres = None
    if os.path.exists(os.path.join(SC, 'rec-wh.json')) or os.path.isdir(ckw):
        rc2, V, tc, se = run_check(inst2, ckw, os.path.join(SC, 'v-wh.json'))
        wres = summ(V)
    rec['cases'].append(dict(name='engine-withheld-row-closure-done', note='engine test hook --withhold 0:3 --gens g --no-repair; the checkpoint is the closure_done one',
                             engine_rc=rc_, checker=wres, passed=bool(wres and (wres['closed'] is False or wres['certificate_valid'] is True))))
    # 5. g replaced by a different generator in the instance file
    Ij = jl(inst)
    for label, mod in (('g-coefficient-changed', 'coef'), ('g-other-target', 'target')):
        I2 = json.loads(json.dumps(Ij))
        if mod == 'coef':
            I2['generator']['1,1'] = (I2['generator']['1,1'] + 5) % p
        else:
            X2 = (Ij['target'][0] + 17) % p
            g2 = C.chain_nf(C.s3_expand(p, Ij['a'], Ij['b'], X2), k, p)
            I2['generator'] = {f'{e1},{e2}': v for (e1, e2), v in g2.items()}
        ip = os.path.join(SC, f'inst-{label}.json')
        json.dump(I2, open(ip, 'w'))
        res = {}
        for mode in ('all',):
            rc2, V, tc, se = run_check(ip, ck, os.path.join(SC, f'v-{label}.json'), ['--shifts', mode])
            res[mode] = summ(V)
            res[mode]['generator_equals_formula'] = V['generator']['instance_equals_formula']
        case(label, res, lambda r: r['certificate_valid'] is False and r['g_in_W'] is False and r['generator_equals_formula'] is False,
             'instance generator dict replaced; the basis was computed for the original g, so g is not in W and the dict disagrees with S_3(a, b, target[0])')
    # 6. format corruption: a tile byte changed (manifest hash), a pivot duplicated, a coefficient outside [0, p)
    root = copy_tree(ck, os.path.join(SC, 'neg-fmt'))
    tp = os.path.join(C.locate_ckpt(root), 'p0_t00000.npy')
    data = bytearray(open(tp, 'rb').read())
    data[-1] ^= 1
    open(tp, 'wb').write(bytes(data))
    rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-fmt1.json'))
    case('tile-byte-flipped', {'all': summ(V)}, lambda r: r['format_ok'] is False and r['certificate_valid'] is False, 'sha256 of a tile differs from the manifest')
    shutil.rmtree(root)
    piv, B = parts[0]
    piv2 = piv.copy()
    piv2[1] = piv2[0]
    root = os.path.join(SC, 'neg-fmt2')
    if os.path.exists(root):
        shutil.rmtree(root)
    C.write_checkpoint(root, [(piv2, np.zeros((len(piv2), L.parts[0].c), dtype=np.int64)), (parts[1][0], np.zeros((len(parts[1][0]), L.parts[1].c), dtype=np.int64))], cfg, tile_rows=64)
    # files of the manifest are consistent with the (garbage) content; only the pivot list is invalid
    rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-fmt2.json'))
    case('duplicate-pivot', {'all': summ(V)}, lambda r: r['format_ok'] is False and r['certificate_valid'] is False, 'two rows with the same pivot column')
    shutil.rmtree(root)
    B4 = B.copy()
    B4[3, 2] = p
    root = os.path.join(SC, 'neg-fmt3')
    if os.path.exists(root):
        shutil.rmtree(root)
    save_parts(root, man, [(piv, B4), parts[1]], cs=L)
    rc2, V, tc, se = run_check(inst, root, os.path.join(SC, 'v-fmt3.json'))
    case('coefficient-equal-p', {'all': summ(V)}, lambda r: r['format_ok'] is False and r['certificate_valid'] is False, 'an entry equal to p (outside [0, p))')
    shutil.rmtree(root)
    rec['passed'] = bool(ok)
    return save_group('negative', rec)


# ------------------------------------------------------------------------- group: devices
def group_devices():
    rec = dict(group='devices', cases=[])
    ok = True
    for name in ('proper-p15361-k7-D4', 'unit-p12289-k8-D4'):
        row = [x for x in POS if x[0] == name][0]
        inst, ck, recp, rc, tsec, err = produce(*row[:8])
        runs = {}
        for label, extra in (('cpu', []), ('fake-gpu', ['--device', 'fake-gpu']),
                             ('cpu-tiny-blocks', ['--col-block', '5', '--row-chunk', '7']),
                             ('fake-gpu-tiny-blocks', ['--device', 'fake-gpu', '--col-block', '3', '--row-chunk', '4'])):
            rc2, V, tc, se = run_check(inst, ck, os.path.join(SC, f'v-dev-{name}-{label}.json'), extra)
            runs[label] = summ(V)
        # chunked inner dimension: run in-process with a tiny forced max_inner
        C.MI_OVERRIDE = 7
        try:
            ns = argparse.Namespace(instance=inst, ckpt=ck, D=None, device='cpu', shifts='all', out=None, status=None, engine_record=recp,
                                    gen_mult_weight=1, col_block=None, row_chunk=64, work_gb=0.05, progress_seconds=1e9, fail_fast=False)
            import io
            old = sys.stdout
            sys.stdout = io.StringIO()
            try:
                V = C.verify(ns)
            finally:
                sys.stdout = old
        finally:
            C.MI_OVERRIDE = None
        runs['cpu-forced-inner-7'] = summ(V)
        keys = ('closed', 'g_in_W', 'one_excluded', 'dimension', 'profile', 'basis_sha256', 'rows_checked', 'triples_checked', 'failures')
        base = {x: runs['cpu'][x] for x in keys}
        same = all({x: r[x] for x in keys} == base for r in runs.values())
        rec['cases'].append(dict(name=name, runs=runs, identical_verdicts=bool(same), passed=bool(same and runs['cpu']['certificate_valid'])))
        ok &= bool(same and runs['cpu']['certificate_valid'])
        print('  devices', name, 'ok' if same else 'FAIL', flush=True)
    # the built-in self-test with fake-gpu through the CLI
    rc, so, se = sh([sys.executable, CHECK, '--selftest', '--device', 'fake-gpu'])
    rec['cases'].append(dict(name='cli --selftest --device fake-gpu', rc=rc, passed=rc == 0))
    ok &= rc == 0
    rec['passed'] = bool(ok)
    return save_group('devices', rec)


def main():
    global SC
    ap = argparse.ArgumentParser()
    ap.add_argument('--scratch', default='/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/certcheck')
    ap.add_argument('--groups', default='arith,spec,positive,reference,negative,devices')
    ap.add_argument('--only', default=None, help='positive group: comma separated name fragments')
    ap.add_argument('--keep', action='store_true')
    a = ap.parse_args()
    SC = a.scratch
    os.makedirs(SC, exist_ok=True)
    allok = True
    for g in a.groups.split(','):
        t = time.time()
        if g == 'positive':
            allok &= group_positive(a.only.split(',') if a.only else None)
        else:
            allok &= {'arith': group_arith, 'spec': group_spec, 'reference': group_reference, 'negative': group_negative,
                      'devices': group_devices}[g]()
        print(f'  ({g}: {time.time() - t:.1f}s)', flush=True)
    if not a.keep:
        shutil.rmtree(SC, ignore_errors=True)
    print('ALL PASSED' if allok else 'SOME FAILED')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
