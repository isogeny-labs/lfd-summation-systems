"""Validation checks for csd/ic_bounds_msum.py (called by ``scripts/ic_bounds_msum_run.py validate``).

Each check returns a JSON-able dict with ``ok`` and the exact statement of what was compared.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import time

import numpy as np

from . import ic_dreg as I
from . import ic_bounds_m3 as M
from . import ic_bounds_msum as S
from . import ic_curve_sweep as CS
from .ic_arbitrary_bases import canonical_basis
from .ic_research_campaign import field

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
EXP = os.path.join(ROOT, 'experiments')


def _rows(path):
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]


# ----------------------------------------------------------------------------------------------- summation values

def check_summation_prime(seed=1):
    """V1: the evaluation recursion against s4_value (m = 3) at random points, degenerate ones included; permutation
    symmetry and split independence for S_5, S_6 (exact equality); degree 2^(n-2) in every variable; vanishing at the
    x-coordinates of random rational decompositions P_1 + ... + P_m = R (m = 2..5)."""
    rng = random.Random(seed)
    out = dict(name='summation-prime')
    ratios = {}
    degenerate = 0
    total = 0
    for p, (a, b) in [(17, (3, 5)), (97, (44, 26)), (1000003, (123457, 99991))]:
        E = I.PrimeCurve(p, a, b)
        cur = S.prime_curve_ops(p, a, b)
        fast, slow = S.Summation(cur, True), S.Summation(cur, False)
        bad = 0
        vals = set()
        for it in range(400):
            xs = [rng.randrange(p) for _ in range(4)]
            if it % 5 == 1:
                xs[1] = xs[0]; degenerate += 1
            if it % 5 == 2:
                xs[3] = xs[2]; degenerate += 1
            if it % 7 == 3:
                xs[rng.randrange(4)] = 0; degenerate += 1
            if it % 11 == 4:
                xs[0] = xs[1] = xs[2] = xs[3]; degenerate += 1
            ref = M.s4_value(E, *xs)
            v1, v2 = fast.value(xs), slow.value(xs)
            total += 1
            if v1 != v2:
                bad += 1
            if ref != 0:
                vals.add(v1 * pow(ref, -1, p) % p)
            elif v1 != 0:
                bad += 1
        ratios[p] = sorted(vals)
        out[f'p{p}_fast_vs_sylvester_mismatches'] = bad
    out['ratio_to_s4_value_by_prime'] = ratios       # one global constant per prime expected: {1}
    out['points_tested'] = total
    out['degenerate_points'] = degenerate
    # symmetry, split independence, degrees
    p = 97
    cur = S.prime_curve_ops(p, 44, 26)
    Sm = S.Summation(cur)
    sym = split = 0
    for n in (4, 5, 6, 7):
        for _ in range(60):
            ys = [rng.randrange(p) for _ in range(n)]
            if n == 7 and _ > 10:
                break
            v = Sm.value(ys)
            pr = ys[:]; rng.shuffle(pr)
            sym += Sm.value(pr) != v
            for j in range(2, n - 1):
                split += Sm.value(ys, split=j) != v
    out['permutation_mismatches_n4to7'] = int(sym)
    out['split_mismatches_n4to7'] = int(split)
    degs = {}
    for n in (3, 4, 5, 6):
        ds = []
        for pos in range(n):
            deg = -1
            for _trial in range(4):      # the top coefficient may vanish at one random choice of the other variables
                ys = [rng.randrange(p) for _ in range(n)]
                vals = []
                for x in range(2 ** (n - 2) + 6):
                    y = ys[:]; y[pos] = x; vals.append(Sm.value(y))
                d = vals; dd = -1
                for t in range(len(vals)):
                    if any(v % p for v in d):
                        dd = t
                    d = [(d[i + 1] - d[i]) % p for i in range(len(d) - 1)]
                    if not d:
                        break
                deg = max(deg, dd)
            ds.append(deg)
        degs[n] = ds
    out['degree_in_each_variable_by_n'] = degs
    # decompositions
    van = {}
    for m in (2, 3, 4, 5):
        tested = bad = 0
        for p in (97, 1009, 7681):
            a, b = rng.randrange(1, p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            E = I.PrimeCurve(p, a, b)
            Sm = S.Summation(S.prime_curve_ops(p, a, b))
            for _ in range(80):
                Ps = [E.random_point(rng) for _ in range(m)]
                if m >= 2 and rng.random() < .15:
                    Ps[1] = Ps[0]
                R = None
                for P in Ps:
                    R = E.add(R, P)
                if R is None:
                    continue
                tested += 1
                if Sm.value([P[0] for P in Ps] + [R[0]]) != 0:
                    bad += 1
        van[m] = dict(tested=tested, nonzero=bad)
    out['vanishing_at_decompositions'] = van
    out['ok'] = (all(v == 0 for k, v in out.items() if k.endswith('mismatches')) and sym == 0 and split == 0
                 and all(len(r) <= 1 and (not r or r == [1]) for r in ratios.values())
                 and all(d == [2 ** (n - 2)] * n for n, d in degs.items())
                 and all(v['nonzero'] == 0 and v['tested'] > 150 for v in van.values()))
    return out


def check_summation_binary(seed=2):
    out = dict(name='summation-binary')
    rng = random.Random(seed)
    res = {}
    ok = True
    for n in (5, 6, 8):
        F = field(n)
        for (a, b) in [(0, 1), (1, 1), (rng.randrange(1 << n), rng.randrange(2, 1 << n))]:
            E = I.BinaryCurve(F, a, b)
            cur = S.binary_curve_ops(F, a, b)
            fast, slow = S.Summation(cur, True), S.Summation(cur, False)
            mism = 0
            for it in range(300):
                xs = [rng.randrange(1 << n) for _ in range(4)]
                if it % 5 == 1:
                    xs[1] = xs[0]
                if it % 7 == 3:
                    xs[2] = 0
                ref = I.s4_binary_value(F, b, *xs)
                mism += (fast.value(xs) != ref) + (slow.value(xs) != ref)
            sym = 0
            for it in range(100):
                for nn in (5, 6):
                    ys = [rng.randrange(1 << n) for _ in range(nn)]
                    v = fast.value(ys)
                    pr = ys[:]; rng.shuffle(pr)
                    sym += fast.value(pr) != v
                    sym += sum(fast.value(ys, split=j) != v for j in range(2, nn - 1))
            vanish = tested = 0
            for m in (2, 3, 4, 5):
                for _ in range(40):
                    Ps = [E.random_point(rng) for _ in range(m)]
                    R = None
                    for P in Ps:
                        R = E.add(R, P)
                    if R is None:
                        continue
                    tested += 1
                    vanish += fast.value([P[0] for P in Ps] + [R[0]]) != 0
            res[f'n{n}_a{a}_b{b}'] = dict(s4_mismatches=int(mism), sym_split_mismatches=int(sym), decomposition_tested=tested,
                                        decomposition_nonzero=int(vanish))
            ok &= mism == 0 and sym == 0 and vanish == 0
    out['rows'] = res
    out['statement'] = ('exact equality (characteristic 2, no constant) of the recursion with csd.ic_dreg.s4_binary_value at 900 random points '
                        'incl. equal coordinates and zeros; permutation and split independence of S_5, S_6; vanishing at x-coordinates of random decompositions m = 2..5')
    out['ok'] = bool(ok)
    return out


# ----------------------------------------------------------------------------------------------- prime generators

def _ratio(new, old, p):
    if set(new) != set(old):
        return None
    ratios = {new[m] * pow(old[m], -1, p) % p for m in old}
    return sorted(ratios)


def check_prime_generators():
    """V2: generators by evaluation (grid + inverse DFT) against the existing builders: m = 2 prime_coset_generator (k = 2, 3, 4, 5),
    m = 3 direct_generator (k = 2, 3, 4); same support and ratio 1.  Plus the recorded generator_hash of the depth study."""
    out = dict(name='prime-generators', m2=[], m3=[], depth_hash=[])
    ok = True
    rng = random.Random(7)
    for k in (2, 3, 4, 5):
        count = 0
        for p in (97, 193, 257, 641, 7681, 449, 577):
            if (p - 1) % (1 << k):
                continue
            for _ in range(2):
                a, b = rng.randrange(1, p), rng.randrange(1, p)
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                E = I.PrimeCurve(p, a, b)
                t = rng.randrange(p)
                g, grid, g0 = S.prime_generator(E, t, 2, k)
                old = CS.prime_coset_generator(E, t, k, 1, 1)
                ratio = _ratio(g, old, p)
                out['m2'].append(dict(p=p, a=a, b=b, t=t, k=k, ratio=ratio, terms=len(g)))
                ok &= ratio == [1]
                count += 1
                if count >= 5:
                    break
            if count >= 5:
                break
    for k in (2, 3, 4):
        count = 0
        for p in (97, 193, 257, 641, 7681, 449):
            for _ in range(2):
                a, b = rng.randrange(1, p), rng.randrange(1, p)
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                E = I.PrimeCurve(p, a, b)
                t = rng.randrange(p)
                g, grid, g0 = S.prime_generator(E, t, 3, k)
                old = M.direct_generator(E, t, k)
                ratio = _ratio(g, old, p)
                out['m3'].append(dict(p=p, a=a, b=b, t=t, k=k, ratio=ratio, terms=len(g), degree=S.generator_degree(g)))
                ok &= ratio == [1]
                count += 1
                if count >= 6:
                    break
            if count >= 6:
                break
    # recorded generator hashes of experiments/ic-bounds-m3-depth (instances p97-r ... p641-pl, k = 2..5)
    import importlib.util
    depth = S.load_depth_module()
    rec = {}
    for r in _rows(os.path.join(EXP, 'ic-bounds-m3-depth', 'closures.jsonl')):
        rec[(r['instance'], r['k'])] = r['generator_hash']
    for name, spec in depth.INSTANCES.items():
        inst = depth.make_instance(*spec)
        E = I.PrimeCurve(inst['p'], inst['a'], inst['b'])
        for k in (2, 3, 4, 5):
            if (name, k) not in rec:
                continue
            g, grid, g0 = S.prime_generator(E, inst['t'], 3, k)
            h = hashlib.sha256(json.dumps(sorted(g.items())).encode()).hexdigest()[:16]
            out['depth_hash'].append(dict(instance=name, k=k, recorded=rec[(name, k)], ours=h, equal=h == rec[(name, k)]))
            ok &= h == rec[(name, k)]
    out['statement'] = ('same support and every coefficient equal (ratio {1}) between the evaluation-based normal form and the existing '
                        'symbolic builders: m = 2 (csd.ic_curve_sweep.prime_coset_generator, scale 1), m = 3 (csd.ic_bounds_m3.direct_generator); '
                        'sha256 of the m = 3 generators equal to the recorded generator_hash of experiments/ic-bounds-m3-depth')
    out['counts'] = dict(m2=len(out['m2']), m3=len(out['m3']), depth_hash=len(out['depth_hash']))
    out['ok'] = bool(ok)
    return out


# ----------------------------------------------------------------------------------------------- binary generators

def check_binary_generators():
    out = dict(name='binary-generators', m2=[], m3=[])
    ok = True
    rng = random.Random(11)
    # m = 2 and m = 3: exact equality of the canonical equation lists (the Boolean normal form is unique)
    for n in (5, 6, 7, 8, 9):
        F = field(n)
        for (a, b) in [(0, 1), (1, 1), (rng.randrange(1 << n), rng.randrange(2, 1 << n))]:
            E = I.BinaryCurve(F, a, b)
            for k in (2, 3, 4):
                if k >= n:
                    continue
                basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
                for tx in (0, rng.randrange(1, 1 << n)):
                    for m in (2, 3):
                        ref = tuple(sorted(set(tuple(sorted(e)) for e in I.binary_descent(E, [basis] * m, tx) if e)))
                        eq, truth = S.binary_generator(E, basis, m, tx)
                        same = list(eq) == list(ref)
                        sols = I.binary_solutions(E, [basis] * m, tx)
                        zs = sorted(sum(c << (i * k) for i, c in enumerate(t)) for t in sols)
                        zs2 = S.binary_zero_set(truth)
                        row = dict(n=n, a=a, b=b, k=k, m=m, tx=tx, equations_equal=same, Z_equal=zs == zs2, Z=len(zs2))
                        out['m2' if m == 2 else 'm3'].append(row)
                        ok &= same and zs == zs2
    # m = 3: against ic_bounds_binary.m3_direct_system (including its system_hash) on recorded rows
    import csd.ic_bounds_binary as B
    recs = _rows(os.path.join(EXP, 'ic-bounds-binary', 'results-m3-direct.jsonl'))
    sel = []
    for k in (2, 3, 4, 5):
        sel += [r for r in recs if r['k'] == k][:6]
    out['recorded_hash'] = []
    for r in sel:
        F = I.GF2n(r['n'], r['modulus'])
        E = I.BinaryCurve(F, r['a'], r['b'])
        eq, truth = S.binary_generator(E, r['basis'], 3, r['target'][0])
        h = S.binary_system_hash('binary-m3-direct', 3 * r['k'], eq)
        pts = S.binary_zero_set(truth)
        out['recorded_hash'].append(dict(case_id=r['case_id'], k=r['k'], n=r['n'], hash_equal=h == r['system_hash'],
                                         Z=len(pts), Z_recorded=r['solutions'], Z_equal=len(pts) == r['solutions'],
                                         points_sha_prefix_equal=hashlib.sha256(json.dumps(pts).encode()).hexdigest()[:16] == r['points_sha256']))
        ok &= h == r['system_hash'] and len(pts) == r['solutions']
    out['statement'] = ('exact equality of the canonical Boolean equation lists (sorted tuples of sorted monomial-mask lists; this is the unique '
                        'Boolean normal form, so the polynomial sets are equal, not only their spans) between the Moebius transform of the evaluated '
                        'S_{m+1} truth table and csd.ic_dreg.binary_descent for m = 2 and 3; Z equal to csd.ic_dreg.binary_solutions; '
                        'system_hash equal to the recorded one of experiments/ic-bounds-binary/results-m3-direct.jsonl on %d rows' % len(sel))
    out['counts'] = dict(m2=len(out['m2']), m3=len(out['m3']), recorded=len(out['recorded_hash']))
    out['ok'] = bool(ok)
    return out


def check_binary_m45_consistency():
    """m = 4, 5 have no existing builder: the evaluated truth tables for different splits and permuted block order are equal
    (so the Boolean functions are the symmetric summation polynomial's), planted decompositions are in Z, the Z found by the group law
    (rational decompositions) is a subset of the zero set."""
    out = dict(name='binary-m45-consistency', rows=[])
    ok = True
    rng = random.Random(13)
    for (m, n, k) in [(4, 6, 2), (4, 7, 2), (4, 8, 3), (5, 6, 2), (5, 7, 2)]:
        F = field(n)
        for (a, b) in [(0, 1), (rng.randrange(1 << n), rng.randrange(2, 1 << n))]:
            E = I.BinaryCurve(F, a, b)
            basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
            V = I.span(basis)
            pl = S.planted_target_m(E, V, m, rng)
            if pl is None:
                continue
            R, Ps = pl
            eq1, tr1 = S.binary_generator(E, basis, m, R[0])
            alt = 2 if m == 5 else 3
            eq2, tr2 = S.binary_generator(E, basis, m, R[0], split=alt)
            ops = S.BinOps(F)
            eq3, tr3 = S.binary_generator(E, basis, m, R[0], fast=False, ops=ops)
            z = S.binary_zero_set(tr1)
            pt = sum(V.index(P[0]) << (i * k) for i, P in enumerate(Ps))
            rat = S.rational_decompositions_binary(E, V, R, m)
            idx = {x: i for i, x in enumerate(V)}
            ratz = {tuple(sorted(t)) for t in rat}
            zsets = {tuple(sorted(V[(c >> (i * k)) & ((1 << k) - 1)] for i in range(m))) for c in z}
            row = dict(m=m, n=n, k=k, a=a, b=b, split_equal=bool(np.array_equal(tr1, tr2)),
                       fast_vs_sylvester_equal=bool(np.array_equal(tr1, tr3)), planted_in_Z=pt in set(z), Z=len(z),
                       rational_multisets=len(rat), rational_subset_of_Z=ratz <= zsets)
            out['rows'].append(row)
            ok &= row['split_equal'] and row['fast_vs_sylvester_equal'] and row['planted_in_Z'] and row['rational_subset_of_Z']
    out['ok'] = bool(ok)
    return out


# ----------------------------------------------------------------------------------------------- closures: positive controls

def _prime_ladder_for(p, m, k, g, pts, g0, engine='auto', cap=None, sandwich='auto'):
    n = sum(1 for _ in range(1))  # placeholder to keep flake quiet
    eng = engine
    if engine == 'auto':
        eng = 'dense' if (1 << (m * k)) <= 1100 else 'compact'
    return S.prime_ladder(p, m, k, g, pts, g0, cap=cap, engine=eng, sandwich=sandwich)


def control_m2_prime(per_k=6):
    """m = 2 prime subgroup systems of experiments/ic-curve-sweep (k = 2..5): generator hash (the recorded system_hash covers the
    generator list), exact LFD, strict falls, status, certificate degree, all equal to the recorded values."""
    cases = _rows(os.path.join(EXP, 'ic-curve-sweep', 'cases.jsonl'))
    res = {}
    for r in _rows(os.path.join(EXP, 'ic-curve-sweep', 'results.jsonl')):
        res[r['system_hash']] = r
    rows = []
    ok = True
    for k in (2, 3, 4, 5):
        sel = [c for c in cases if c['kind'] == 'prime' and c['base_kind'] == 'subgroup' and c['k'] == k and c['system_hash'] in res]
        # mix of unit and proper, different target classes
        sel.sort(key=lambda c: (res[c['system_hash']]['status'], c['target_class']))
        pick = []
        seen = set()
        for c in sel:
            key = (res[c['system_hash']]['status'], c['target_class'])
            if key not in seen:
                pick.append(c); seen.add(key)
        for c in sel[::max(1, len(sel) // per_k)]:
            if c not in pick:
                pick.append(c)
        for c in pick[:per_k + 2]:
            p = c['prime']
            E = I.PrimeCurve(p, c['a'], c['b'])
            t = c['target'][0]
            g, grid, g0 = S.prime_generator(E, t, 2, k)
            sysh = hashlib.sha256(json.dumps(dict(kind='prime', prime=p, m=2, k=k, generators=[sorted(g.items())]),
                                             sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            pts = S.prime_zero_set(grid, 2, k)
            old_pts = CS.prime_points(E, k, g0, 1, 1, t)
            out = S.prime_ladder(p, 2, k, g, pts, g0, cap=2 * k + 1, engine='dense', sandwich=False)
            rec = res[c['system_hash']]
            row = dict(k=k, p=p, a=c['a'], b=c['b'], t=t, target_class=c['target_class'], hash_equal=sysh == c['system_hash'],
                       Z_equal=sorted(pts) == sorted(old_pts), new_lfd=out['exact_lfd'], recorded_lfd=rec['exact_lfd'],
                       new_status=out['status'], recorded_status=rec['status'], new_falls=out['strict_falls'],
                       recorded_falls=rec['strict_falls'])
            row['equal'] = (row['hash_equal'] and row['Z_equal'] and row['new_lfd'] == row['recorded_lfd']
                            and row['new_status'] == row['recorded_status'] and row['new_falls'] == row['recorded_falls'])
            ok &= row['equal']
            rows.append(row)
    return dict(name='control-m2-prime', ok=bool(ok), n=len(rows), rows=rows,
                statement='generator-list hash (= recorded system_hash), Z, exact LFD, status and strict falls equal to experiments/ic-curve-sweep rows')


def control_m2_binary(per_k=5):
    cases = _rows(os.path.join(EXP, 'ic-curve-sweep', 'cases.jsonl'))
    res = {r['system_hash']: r for r in _rows(os.path.join(EXP, 'ic-curve-sweep', 'results.jsonl'))}
    rows = []
    ok = True
    for k in (4, 5, 6, 7):
        sel = [c for c in cases if c['kind'] in ('binary', 'koblitz') and c['k'] == k and c['system_hash'] in res]
        sel.sort(key=lambda c: (res[c['system_hash']]['status'], c['kind'], c['base_kind'], c['target_class']))
        pick, seen = [], set()
        for c in sel:
            key = (res[c['system_hash']]['status'], c['kind'], c['target_class'])
            if key not in seen:
                pick.append(c); seen.add(key)
        for c in pick[:per_k + 3]:
            F = I.GF2n(c['n'], c['modulus'])
            E = I.BinaryCurve(F, c['a'], c['b'])
            eq, truth = S.binary_generator(E, c['basis'], 2, c['target'][0])
            h = S.binary_system_hash('binary', 2 * k, eq)
            pts = S.binary_zero_set(truth)
            out = S.binary_ladder(eq, 2 * k, pts, cap=2 * k + 1, engine='py')
            rec = res[c['system_hash']]
            row = dict(k=k, n=c['n'], kind=c['kind'], base=c['base_kind'], target_class=c['target_class'], hash_equal=h == c['system_hash'],
                       Z=len(pts), Z_recorded=rec.get('point_count'), new_lfd=out['exact_lfd'], recorded_lfd=rec['exact_lfd'],
                       new_status=out['status'], recorded_status=rec['status'], new_falls=out['strict_falls'],
                       recorded_falls=rec['strict_falls'])
            row['equal'] = (row['hash_equal'] and row['Z'] == row['Z_recorded'] and row['new_lfd'] == row['recorded_lfd']
                            and row['new_status'] == row['recorded_status'] and row['new_falls'] == row['recorded_falls'])
            ok &= row['equal']
            rows.append(row)
    return dict(name='control-m2-binary', ok=bool(ok), n=len(rows), rows=rows,
                statement='system_hash, |Z|, exact LFD, status, strict falls equal to experiments/ic-curve-sweep rows')


def control_m3_prime(ks=(2, 3, 4)):
    """m = 3 prime: the 12 depth-study instances at k = 2, 3, 4: new generator hash, solution count, profile p_j of V_D at every recorded D,
    exact LFD 7 (README of experiments/ic-bounds-m3-depth)."""
    depth = S.load_depth_module()
    recorded = {}
    for r in _rows(os.path.join(EXP, 'ic-bounds-m3-depth', 'closures.jsonl')):
        recorded[(r['instance'], r['k'], r['D'])] = r
    rows = []
    ok = True
    for name, spec in depth.INSTANCES.items():
        inst = depth.make_instance(*spec)
        E = I.PrimeCurve(inst['p'], inst['a'], inst['b'])
        for k in ks:
            if (name, k, 7) not in recorded:
                continue
            t0 = time.time()
            g, grid, g0 = S.prime_generator(E, inst['t'], 3, k)
            pts = S.prime_zero_set(grid, 3, k)
            out = S.prime_ladder(inst['p'], 3, k, g, pts, g0, engine='dense' if k <= 2 else 'compact', sandwich='auto')
            cmp_rows = []
            for row in out['profile']:
                rec = recorded.get((name, k, row['D']))
                if rec is not None and row['D'] >= 6:
                    cmp_rows.append(dict(D=row['D'], profile_equal=row['profile'] == rec['profile'], rank_equal=row['quotient_rank'] == rec['rank']))
            r = dict(instance=name, k=k, p=inst['p'], Z=len(pts), Z_recorded=recorded[(name, k, 7)]['solutions'],
                     lfd=out['exact_lfd'], status=out['status'], falls=out['strict_falls'], compared=cmp_rows,
                     seconds=round(time.time() - t0, 1))
            r['equal'] = (r['lfd'] == 7 and r['Z'] == r['Z_recorded'] and cmp_rows and all(c['profile_equal'] and c['rank_equal'] for c in cmp_rows))
            ok &= bool(r['equal'])
            rows.append(r)
    return dict(name='control-m3-prime', ok=bool(ok), n=len(rows), rows=rows,
                statement='LFD 7, |Z|, V_D profile and rank at every recorded D >= 6 equal to experiments/ic-bounds-m3-depth/closures.jsonl')


def control_m3_binary(per_k=7):
    recs = _rows(os.path.join(EXP, 'ic-bounds-binary', 'results-m3-direct.jsonl'))
    rows = []
    ok = True
    rng = random.Random(3)
    for k in (2, 3, 4, 5):
        pool = [r for r in recs if r['k'] == k]
        # mix unit / proper, input degree 3 / 6
        pool.sort(key=lambda r: (r['status'], r['input_degree'], r['n'], r['case_id']))
        pick, seen = [], set()
        for r in pool:
            key = (r['status'], r['input_degree'], r['exact_lfd'])
            if key not in seen:
                pick.append(r); seen.add(key)
        rng.shuffle(pool)
        for r in pool:
            if len(pick) >= per_k:
                break
            if r not in pick:
                pick.append(r)
        for r in pick[:per_k]:
            F = I.GF2n(r['n'], r['modulus'])
            E = I.BinaryCurve(F, r['a'], r['b'])
            eq, truth = S.binary_generator(E, r['basis'], 3, r['target'][0])
            h = S.binary_system_hash('binary-m3-direct', 3 * k, eq)
            pts = S.binary_zero_set(truth)
            t0 = time.time()
            out = S.binary_ladder(eq, 3 * k, pts, cap=3 * k + 1, engine='py' if k <= 4 else 'big')
            row = dict(k=k, n=r['n'], case_id=r['case_id'], hash_equal=h == r['system_hash'], Z=len(pts), Z_recorded=r['solutions'],
                       new_lfd=out['exact_lfd'], recorded_lfd=r['exact_lfd'], new_status=out['status'], recorded_status=r['status'],
                       new_falls=out['strict_falls'], recorded_falls=r['strict_falls'], new_cert=out.get('certificate_degree'),
                       recorded_cert=r.get('certificate_degree'), seconds=round(time.time() - t0, 1))
            row['equal'] = (row['hash_equal'] and row['Z'] == row['Z_recorded'] and row['new_lfd'] == row['recorded_lfd']
                            and row['new_status'] == row['recorded_status'] and row['new_falls'] == row['recorded_falls']
                            and row['new_cert'] == row['recorded_cert'])
            ok &= row['equal']
            rows.append(row)
    return dict(name='control-m3-binary', ok=bool(ok), n=len(rows), rows=rows,
                statement='system_hash, |Z|, exact LFD, status, strict falls, certificate degree equal to experiments/ic-bounds-binary/results-m3-direct.jsonl')


# ----------------------------------------------------------------------------------------------- engine cross-checks

def engines_prime():
    """dense chain_degree_closure vs CompactClosure (with and without the dimension sandwich): quotient rank and low profile equal at every D."""
    rows = []
    ok = True
    rng = random.Random(5)
    for (m, k, p) in [(2, 3, 97), (3, 2, 97), (3, 2, 7681), (3, 3, 257), (4, 2, 97), (4, 2, 257), (5, 2, 641)]:
        for planted in (True, False):
            inst = S.make_prime_instance(p, m, k, 3, planted)
            E = I.PrimeCurve(p, inst['a'], inst['b'])
            g, grid, g0 = S.prime_generator(E, inst['t'], m, k)
            pts = S.prime_zero_set(grid, m, k)
            res = {}
            for eng, sand in (('dense', False), ('compact', False), ('compact', True)):
                if eng == 'dense' and (1 << (m * k)) > 1100:
                    continue
                res[(eng, sand)] = S.prime_ladder(p, m, k, g, pts, g0, engine=eng, sandwich=sand)
            keys = list(res)
            base = res[keys[0]]
            same = all(res[kk]['exact_lfd'] == base['exact_lfd'] and res[kk]['strict_falls'] == base['strict_falls']
                       and [(x['quotient_rank'], x['profile']) for x in res[kk]['profile']] == [(x['quotient_rank'], x['profile']) for x in base['profile']]
                       for kk in keys)
            rows.append(dict(m=m, k=k, p=p, planted=planted, Z=len(pts), engines=[f'{a}/sandwich={b}' for a, b in keys],
                             lfd=base['exact_lfd'], status=base['status'], equal=same))
            ok &= same
    return dict(name='engines-prime', ok=bool(ok), n=len(rows), rows=rows,
                statement='per-D quotient rank and low profile identical between chain_degree_closure (dense), CompactClosure without and with the sandwich')


def engines_binary():
    rows = []
    ok = True
    rng = random.Random(6)
    for (m, n, k) in [(3, 7, 2), (3, 8, 3), (4, 6, 2), (4, 8, 3), (5, 6, 2)]:
        F = field(n)
        for (a, b) in [(0, 1), (rng.randrange(1 << n), rng.randrange(2, 1 << n))]:
            E = I.BinaryCurve(F, a, b)
            basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
            V = I.span(basis)
            pl = S.planted_target_m(E, V, m, rng)
            tx = pl[0][0] if pl else 0
            eq, truth = S.binary_generator(E, basis, m, tx)
            pts = S.binary_zero_set(truth)
            r = m * k
            a1 = S.binary_ladder(eq, r, pts, engine='py')
            a2 = S.binary_ladder(eq, r, pts, engine='big')
            same = (a1['exact_lfd'] == a2['exact_lfd'] and a1['strict_falls'] == a2['strict_falls'] and a1['quotient_ranks'] == a2['quotient_ranks']
                    and a1['status'] == a2['status'])
            rows.append(dict(m=m, n=n, k=k, a=a, b=b, Z=len(pts), lfd_py=a1['exact_lfd'], lfd_big=a2['exact_lfd'], equal=same))
            ok &= same
    return dict(name='engines-binary', ok=bool(ok), n=len(rows), rows=rows,
                statement='csd.ic_local_certificate.local_profile (python) vs csd.ic_bounds_bin_big.local_profile_bin (bit-packed, cpu): LFD, falls, quotient ranks equal')


# ----------------------------------------------------------------------------------------------- negative controls

def negative_prime():
    """A deliberately wrong generator: one grid value altered.
    (i) a point of Z gets a nonzero value: g' with the TRUE Z must fail the certificate (V_D not inside the true ideal, or no local certificate);
    (ii) a nonzero value is set to zero (extra zero): same test against the true Z;
    (iii) the consistent altered system (g', Z') is reported next to the original (does the LFD change?)."""
    rows = []
    for (m, k, p, planted, seed) in [(3, 2, 97, True, 3), (3, 3, 257, True, 3), (4, 2, 257, True, 3), (4, 2, 97, True, 4), (5, 2, 641, True, 3), (3, 2, 193, False, 3)]:
        inst = S.make_prime_instance(p, m, k, seed, planted)
        E = I.PrimeCurve(p, inst['a'], inst['b'])
        g, grid, g0 = S.prime_generator(E, inst['t'], m, k)
        pts = S.prime_zero_set(grid, m, k)
        N = 1 << k
        flat = grid.reshape(-1)
        base = S.prime_ladder(p, m, k, g, pts, g0, engine='dense' if (1 << (m * k)) <= 1100 else 'compact', sandwich='auto')
        zero_idx = [i for i in range(len(flat)) if flat[i] == 0]
        nz_idx = [i for i in range(len(flat)) if flat[i] != 0]
        trials = []
        if zero_idx:
            trials.append(('zero-to-nonzero', zero_idx[0], 1))
        trials.append(('nonzero-to-zero', nz_idx[0], 0))
        trials.append(('nonzero-shifted', nz_idx[1], (int(flat[nz_idx[1]]) + 1) % p))
        for kind, idx, newv in trials:
            alt = grid.copy().reshape(-1)
            alt[idx] = newv
            alt = alt.reshape(grid.shape)
            coef = S.inverse_dft(alt, p, g0, k)
            g2 = S.coef_to_masks(coef, k)
            pts2 = S.prime_zero_set(alt, m, k)
            eng = 'dense' if (1 << (m * k)) <= 1100 else 'compact'
            try:
                wrong = S.prime_ladder(p, m, k, g2, pts, g0, engine=eng, sandwich=False if eng == 'dense' else 'auto', true_pts=pts)
                outcome = dict(status=wrong['status'], lfd=wrong['exact_lfd'], inclusion_failures=wrong.get('inclusion_failures_vs_true_Z'))
                detected = bool(wrong.get('inclusion_failures_vs_true_Z')) or wrong['status'] != 'exact-proper-local' or wrong['exact_lfd'] != base['exact_lfd']
            except AssertionError as e:
                outcome = dict(assertion=str(e)[:100]); detected = True
            cons = S.prime_ladder(p, m, k, g2, pts2, g0, engine=eng, sandwich=False if eng == 'dense' else 'auto')
            trials_out = dict(kind=kind, grid_index=int(idx), original_lfd=base['exact_lfd'], wrong_vs_true_Z=outcome, detected=detected,
                              consistent_altered=dict(Z=len(pts2), status=cons['status'], lfd=cons['exact_lfd']),
                              consistent_changes_result=cons['exact_lfd'] != base['exact_lfd'] or len(pts2) != len(pts))
            rows.append(dict(m=m, k=k, p=p, Z=len(pts), **trials_out))
    return dict(name='negative-prime', rows=rows, n=len(rows),
                detected_wrong_vs_true_Z=sum(r['detected'] for r in rows),
                zero_to_nonzero_detected=[r['detected'] for r in rows if r['kind'] == 'zero-to-nonzero'],
                ok=all(r['detected'] for r in rows if r['kind'] == 'zero-to-nonzero'))


def negative_binary():
    rows = []
    rng = random.Random(9)
    for (m, n, k) in [(3, 7, 2), (3, 8, 3), (4, 6, 2), (5, 6, 2)]:
        F = field(n)
        E = I.BinaryCurve(F, rng.randrange(1 << n), rng.randrange(2, 1 << n))
        basis = canonical_basis(CS.random_nonstable_subspace(F, k, rng, []))
        V = I.span(basis)
        pl = S.planted_target_m(E, V, m, rng)
        tx = pl[0][0] if pl else 0
        eq, truth = S.binary_generator(E, basis, m, tx)
        pts = S.binary_zero_set(truth)
        r = m * k
        base = S.binary_ladder(eq, r, pts, engine='py')
        zidx = pts[0] if pts else None
        for kind, idx in (('flip-at-zero', zidx), ('flip-at-nonzero', int(np.flatnonzero(truth != 0)[0]))):
            if idx is None:
                continue
            alt = truth.copy()
            alt[idx] ^= np.uint64(1)
            T = alt.copy()
            for i in range(r):
                T = T.reshape(-1, 2, 1 << i)
                T[:, 1, :] ^= T[:, 0, :]
            T = T.reshape(-1)
            eq2 = sorted(set(tuple(sorted(int(x) for x in np.flatnonzero((T >> np.uint64(b)) & np.uint64(1)))) for b in range(n)) - {()})
            try:
                out = S.binary_ladder(eq2, r, pts, engine='py')
                outcome = dict(status=out['status'], lfd=out['exact_lfd'])
                detected = out['status'] != base['status'] or out['exact_lfd'] != base['exact_lfd']
            except ValueError as e:
                outcome = dict(error=str(e)[:80]); detected = True
            pts2 = S.binary_zero_set(alt)
            cons = S.binary_ladder(eq2, r, pts2, engine='py')
            rows.append(dict(m=m, n=n, k=k, kind=kind, Z=len(pts), original=dict(status=base['status'], lfd=base['exact_lfd']),
                             wrong_vs_true_Z=outcome, detected=bool(detected),
                             consistent_altered=dict(Z=len(pts2), status=cons['status'], lfd=cons['exact_lfd'])))
    return dict(name='negative-binary', rows=rows, n=len(rows), detected=sum(r['detected'] for r in rows),
                ok=all(r['detected'] for r in rows if r['kind'] == 'flip-at-zero'))


def prime_z_vs_group_law():
    """Z (zero set of the evaluated S_{m+1} on H^m) against the group law: every ordered tuple in H^m with signed lifts P_i summing to R
    (rational decomposition, found by meet in the middle) lies in Z.  Reported: counts of both (Z is larger: S_{m+1} also vanishes at
    decompositions with y-coordinates outside F_p)."""
    from .ic_bounds_msum import rational_decompositions_prime, make_prime_instance, prime_generator, prime_zero_set
    rows = []
    ok = True
    for (m, k) in [(3, 2), (3, 3), (4, 2), (4, 3), (5, 2)]:
        for (p, planted, seed) in [(97, True, 10), (193, False, 11), (257, True, 12), (449, False, 13), (641, True, 14), (7681, False, 15)]:
            inst = make_prime_instance(p, m, k, seed, planted)
            E = I.PrimeCurve(p, inst['a'], inst['b'])
            g, grid, g0 = prime_generator(E, inst['t'], m, k)
            Z = set(prime_zero_set(grid, m, k))
            H = M.subgroup(p, k)[1]
            R = (inst['t'], inst['y'])
            rat = rational_decompositions_prime(E, H, R, m)
            ratz = {sum(H.index(x) << (i * k) for i, x in enumerate(tup)) for tup in rat}
            row = dict(m=m, k=k, p=p, planted=planted, Z=len(Z), rational_ordered=len(ratz), subset=ratz <= Z)
            ok &= row['subset']
            rows.append(row)
    return dict(name='prime-Z-vs-group-law', ok=bool(ok), n=len(rows), rows=rows,
                statement='for every instance of the m = 3, 4, 5 prime cells: rational decompositions (group law, signed lifts) subset of Z')
