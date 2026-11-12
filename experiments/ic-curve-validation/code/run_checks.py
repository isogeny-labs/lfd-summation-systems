#!/usr/bin/env python3
"""Driver: independent LFD recomputation for a deterministic sample of curve-sweep systems.

usage: run_checks.py OUT.jsonl SELECTOR [SELECTOR ...]
SELECTOR = prime2 | prime3 | bin4 | bin6 | bin8 | hash:<system_hash>
Only data files + the lfd_engine.py in this directory are used (no repository closure engine).
"""
import json
import os
import random
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lfd_engine as L

ROOT = '/home/user/ml-cryptanalysis'
SW = ROOT + '/experiments/ic-curve-sweep/'


def load():
    S, R, C = {}, {}, {}
    for l in open(SW + 'systems.jsonl'):
        s = json.loads(l)
        S[s['system_hash']] = s
    for l in open(SW + 'results.jsonl'):
        r = json.loads(l)
        R[r['system_hash']] = r
    for l in open(SW + 'cases.jsonl'):
        c = json.loads(l)
        C.setdefault(c['system_hash'], []).append(c)
    return S, R, C


# ----- own S3 evaluation (input-integrity cross-check of the stored normalised generator) -----

def s3_prime(p, a, b, x1, x2, X):
    t = (x1 - x2) ** 2 * X * X - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * X + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)
    return t % p


def gf2_mul(a, b, mod):
    n = mod.bit_length() - 1
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> n & 1:
            a ^= mod
    return r


def s3_bin(mod, bconst, x1, x2, x3):
    m = lambda u, v: gf2_mul(u, v, mod)
    t = m(x1, x2) ^ m(x1, x3) ^ m(x2, x3)
    return m(t, t) ^ m(m(x1, x2), x3) ^ bconst


def integrity_prime(s, case, G, pts):
    """G (stored normalised generator, evaluated through the chain) must equal S3(c*u1, c*u2, xR) on mu_N^2,
    and the zero set must equal the stored coordinate_solutions (decoded through generator_value)."""
    p, k = s['prime'], s['k']
    N = 1 << k
    a, b, c, xR = case['a'], case['b'], case['base_scale'], case['target'][0]
    mus = L.mu_N(p, N)
    bad = 0
    for x1 in mus:
        for x2 in mus:
            u = [pow(x1, 1 << j, p) for j in range(k)] + [pow(x2, 1 << j, p) for j in range(k)]
            v = 0
            for e, co in G.items():
                t = co
                for j, ej in enumerate(e):
                    if ej:
                        t = t * pow(u[j], ej, p) % p
                v = (v + t) % p
            if v != s3_prime(p, a, b, c * x1 % p, c * x2 % p, xR):
                bad += 1
    g0 = s['generator_value']
    dl = {pow(g0, e, p): e for e in range(N)}
    mine = sorted(dl[pt[0]] | (dl[pt[k]] << k) for pt in pts)
    return dict(generator_mismatches_on_mu2=bad, zero_set_equals_stored=mine == sorted(s['coordinate_solutions']),
                curve=[p, a, b], scale=c, xR=xR)


def integrity_binary(s, case, pts, r):
    k = r // 2
    mod, bc, basis, xR = case['modulus'], case['b'], case['basis'], case['target'][0]
    vals = []
    for c in range(1 << k):
        x = 0
        for j in range(k):
            if c >> j & 1:
                x ^= basis[j]
        vals.append(x)
    bad = 0
    mine = set()
    for bits in range(1 << r):
        x, y = vals[bits & ((1 << k) - 1)], vals[bits >> k]
        zero = s3_bin(mod, bc, x, y, xR) == 0
        if zero:
            mine.add(bits)
        eq_zero = all((sum(int((bits & m) == m) for m in eq) & 1) == 0 for eq in s['equations'])
        if zero != eq_zero:
            bad += 1
    pset = set(sum(v << b for b, v in enumerate(pt)) for pt in pts)
    return dict(assignments_where_S3_and_equations_disagree=bad, zero_set_equals_S3_zero_set=pset == mine,
                zero_set_equals_stored=sorted(pset) == sorted(s['coordinate_solutions']))


def run_one(h, S, R, C, out, Dcap_extra):
    s, rep, case = S[h], R[h], C[h][0]
    t0 = time.time()
    if s['kind'] == 'prime':
        k, p = s['k'], s['prime']
        r = 2 * k
        gen = [(int(m), int(c)) for m, c in s['generators'][0]]
        inputs, G = L.prime_inputs(gen, k, p)
        pts = L.prime_points(G, k, p)
        integ = integrity_prime(s, case, G, pts)
        Dmax = int(os.environ['LFD_DMAX']) if os.environ.get('LFD_DMAX') else 2 * k + 1 + Dcap_extra
        label = 'prime k=%d p=%d' % (k, p)
    else:
        r = s['variables']
        p = 2
        inputs = L.binary_inputs(s['equations'], r)
        pts = L.binary_points(s['equations'], r)
        integ = integrity_binary(s, case, pts, r)
        k = r // 2
        Dmax = k + 2 + Dcap_extra
        label = 'binary(%s) vars=%d' % (case['kind'], r)
    res = L.exact_lfd(inputs, r, p, pts, Dmax)
    mine_status = 'exact-unit' if len(pts) == 0 else 'exact-proper-local'
    rec = dict(system_hash=h, label=label, case_kind=case['kind'], vars=r,
               repo=dict(exact_lfd=rep['exact_lfd'], status=rep['status'], point_count=rep['point_count'],
                         strict_falls=rep['strict_falls'], quotient_ranks=rep['quotient_ranks']),
               mine=dict(lfd=res['lfd'], lfd_largest_fall=res.get('lfd_largest_fall'), falls=res['falls'],
                         points=len(pts), status=mine_status, certified_at=res['certified_at'],
                         per_degree=[{k2: v for k2, v in d.items()} for d in res['per_degree']]),
               integrity=integ, seconds=round(time.time() - t0, 1),
               maxrss_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    rec['match'] = dict(lfd=res['lfd'] == rep['exact_lfd'], points=len(pts) == rep['point_count'],
                        status=mine_status == rep['status'], strict_falls=res['falls'] == rep['strict_falls'])
    out.write(json.dumps(rec) + '\n')
    out.flush()
    print(h[:12], label, 'mine', res['lfd'], 'repo', rep['exact_lfd'], 'pts', len(pts), rep['point_count'],
          'falls', res['falls'], rep['strict_falls'], 'cert', res['certified_at'], '%.0fs' % rec['seconds'],
          rec['maxrss_mb'], 'MB', flush=True)


def select(S, R, C, sel):
    rng = random.Random(20261105)
    if sel.startswith('hash:'):
        return [sel[5:]]
    if sel.startswith('sample:'):
        _, ks, ns = sel.split(':')
        k, n = int(ks), int(ns)
        pool = sorted(h for h, s in S.items() if s['kind'] == 'prime' and s['k'] == k)
        rng.shuffle(pool)
        strata = {}
        for h in pool:
            strata.setdefault((R[h]['status'], R[h]['exact_lfd']), []).append(h)
        chosen, keys = [], sorted(strata)
        while len(chosen) < n and any(strata.values()):
            for key in keys:
                if strata[key] and len(chosen) < n:
                    chosen.append(strata[key].pop())
        return chosen
    if sel.startswith('all'):
        what = sel[3:]
        if what.startswith('prime'):
            k = int(what[5:])
            return sorted(h for h, s in S.items() if s['kind'] == 'prime' and s['k'] == k)
        v = int(what[3:])
        return sorted(h for h, s in S.items() if s['kind'] == 'binary' and s['variables'] == v)
    if sel in ('prime2', 'prime3'):
        k = int(sel[-1])
        pool = [h for h, s in S.items() if s['kind'] == 'prime' and s['k'] == k]
        pool.sort()
        chosen = []
        # strata: rarely-seen LFD, then proper vs unit across primes and base kinds
        odd = [h for h in pool if R[h]['exact_lfd'] != 4]
        chosen += odd
        for status, quota in (('exact-proper-local', 7 if k == 2 else 9), ('exact-unit', 7 if k == 2 else 9)):
            sub = [h for h in pool if R[h]['status'] == status and R[h]['exact_lfd'] == 4]
            rng.shuffle(sub)
            by = {}
            for h in sub:
                key = (S[h]['prime'], C[h][0]['base_kind'])
                by.setdefault(key, []).append(h)
            picked = []
            keys = sorted(by)
            while len(picked) < quota and any(by.values()):
                for key in keys:
                    if by[key] and len(picked) < quota:
                        picked.append(by[key].pop())
            chosen += picked
        return chosen
    if sel.startswith('bin'):
        v = int(sel[3:])
        pool = sorted(h for h, s in S.items() if s['kind'] == 'binary' and s['variables'] == v)
        quota = {4: 16, 6: 16, 8: 14}[v]
        rng.shuffle(pool)
        strata = {}
        for h in pool:
            key = (C[h][0]['kind'], R[h]['status'], R[h]['exact_lfd'])
            strata.setdefault(key, []).append(h)
        chosen = []
        keys = sorted(strata)
        while len(chosen) < quota and any(strata.values()):
            for key in keys:
                if strata[key] and len(chosen) < quota:
                    chosen.append(strata[key].pop())
        return chosen
    raise SystemExit('bad selector ' + sel)


if __name__ == '__main__':
    outp = sys.argv[1]
    S, R, C = load()
    hs = []
    for sel in sys.argv[2:]:
        hs += select(S, R, C, sel)
    done = set()
    if os.path.exists(outp):
        done = set(json.loads(l)['system_hash'] for l in open(outp))
    with open(outp, 'a') as out:
        for h in hs:
            if h in done:
                continue
            run_one(h, S, R, C, out, 0)
