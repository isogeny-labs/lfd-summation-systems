"""Frozen, setup-charged two-summand torsion crossover campaign."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.field import primitive_root_of_unity, sqrt_mod  # noqa: E402
from csd.ic_dreg import PrimeCurve  # noqa: E402

HERE = Path(__file__).resolve().parent
P, R = 37889, 2389
E = PrimeCurve(P, 1, 0)


def ticks():
    return perf_counter(), process_time()


def elapsed(start):
    return perf_counter() - start[0], process_time() - start[1]


def canonical(a, b):
    return (a, b) if a <= b else (b, a)


def points_at(x):
    p = E.lift_x(x)
    if p is None:
        return ()
    q = E.neg(p)
    return (p,) if q == p else (p, q)


def fixture(n):
    z = primitive_root_of_unity(n, P)
    assert z is not None and pow(z, n // 2, P) != 1
    roots = sorted({pow(z, j, P) for j in range(1, n)})
    assert len(roots) == n - 1 and 1 not in roots
    xs = [x for x in roots if points_at(x)]
    fb = sorted(p for x in xs for p in points_at(x))
    return xs, fb


def generator():
    for x in range(P):
        for p in points_at(x):
            g = E.mul(16, p)
            if g is not None:
                assert E.mul(R, g) is None
                return g
    raise AssertionError('missing subgroup generator')


G = generator()
Q = E.mul(113, G)


def targets(n, fb):
    rng = random.Random(20261125 + n)
    rows = []
    for _ in range(256):
        u, v = rng.randrange(R), rng.randrange(R)
        rows.append({'kind': 'random', 'u': u, 'v': v,
                     'point': E.add(E.mul(u, G), E.mul(v, Q))})
    rng = random.Random(20262125 + n)
    for _ in range(32):
        a, b = rng.choice(fb), rng.choice(fb)
        rows.append({'kind': 'planted', 'point': E.add(a, b),
                     'witness': canonical(a, b)})
    return rows


def direct_scan(fb, panel):
    base = set(fb)
    out = []
    for row in panel:
        t = row['point']
        out.append({canonical(a, b) for a in fb
                    if (b := E.add(t, E.neg(a))) in base})
    return out


def direct_table(fb, panel):
    table = defaultdict(set)
    for i, a in enumerate(fb):
        for b in fb[i:]:
            table[E.add(a, b)].add((a, b))
    return [table.get(row['point'], set()) for row in panel]


def roots(q, r):
    # P(q,r,t) = (q-4r)t^2 - 2(r^2-1)t + q.
    a, b, c = (q - 4*r) % P, (-2*(r*r-1)) % P, q % P
    if a == 0:
        if b == 0:
            return None if c == 0 else ()
        return ((-c * pow(b, -1, P)) % P,)
    d = (b*b - 4*a*c) % P
    s = sqrt_mod(d, P)
    if s is None:
        return ()
    z = pow(2*a, -1, P)
    return tuple({((-b + s)*z) % P, ((-b - s)*z) % P})


def quotient_table(xs, panel, bound):
    start = ticks()
    inv = {x: pow(x-1, -1, P) for x in xs}
    u = {x: ((x+1)*inv[x]) % P for x in xs}
    orbits = defaultdict(set)
    for i, x in enumerate(xs):
        ux = u[x]
        for y in xs[i:]:
            uy = u[y]
            orbits[((ux+uy)**2 % P, ux*uy % P)].add((x, y))
        if i % 8 == 0 and dominated(start, bound):
            return None, 'dominated_censored', 0, len(orbits), elapsed(start)
    table = defaultdict(set)
    universal = set()
    for j, ((q, r), pairs) in enumerate(orbits.items()):
        rt = roots(q, r)
        if rt is None:
            universal.update(pairs)
        else:
            for t in rt:
                table[t].update(pairs)
        if j % 32 == 0 and dominated(start, bound):
            return None, 'dominated_censored', 0, len(orbits), elapsed(start)
    out = []
    for j, row in enumerate(panel):
        target = row['point']
        candidates = table.get(target[0], set()) | universal if target is not None else set()
        # O is outside the affine S3 equation; its inverse-pair locus is exact.
        if target is None:
            candidates = {(x, x) for x in xs}
        found = set()
        for x, y in candidates:
            for a in points_at(x):
                for b in points_at(y):
                    if E.add(a, b) == target:
                        found.add(canonical(a, b))
        out.append(found)
        if j % 8 == 0 and dominated(start, bound):
            return out, 'dominated_censored', len(out), len(orbits), elapsed(start)
    return out, 'complete', len(out), len(orbits), elapsed(start)


def dominated(start, bound):
    if bound is None:
        return False
    w, c = elapsed(start)
    return w > bound[0] and c > bound[1]


def rank(rows, modulus=R):
    pivots = {}
    for raw in rows:
        v = {k: x % modulus for k, x in raw.items() if x % modulus}
        while v:
            k = min(v)
            if k not in pivots:
                z = pow(v[k], -1, modulus)
                pivots[k] = {j: x*z % modulus for j, x in v.items()}
                break
            f = v[k]
            for j, x in pivots[k].items():
                y = (v.get(j, 0) - f*x) % modulus
                if y:
                    v[j] = y
                else:
                    v.pop(j, None)
    return len(pivots)


def relation_rank(fb, panel, answers):
    # Work in the prime-order quotient. A column represents
    # (1/16) log_G([16]P); P and -P share one column with opposite signs.
    reps = sorted({min(p, E.neg(p)) for p in fb
                   if E.mul(16, p) is not None})
    cols = {p: i for i, p in enumerate(reps)}
    rows = []
    for row, pairs in zip(panel[:256], answers[:256]):
        assert row['kind'] == 'random'
        for a, b in pairs:
            v = {}
            for p in (a, b):
                rep = min(p, E.neg(p))
                if rep in cols:
                    j = cols[rep]
                    v[j] = v.get(j, 0) + (1 if p == rep else -1)
            # Q=113G, so the challenge has one known right-hand column.
            v[len(reps)] = -(row['u'] + 113*row['v'])
            rows.append(v)
    return rank(rows)


def run_cell(n, rep):
    setup_start = ticks()
    xs, fb = fixture(n)
    panel = targets(n, fb)
    setup_wall, setup_cpu = elapsed(setup_start)
    arms = {}
    answers = {}
    order = ('scan', 'table') if rep % 2 == 0 else ('table', 'scan')
    for name in order:
        start = ticks()
        ans = direct_scan(fb, panel) if name == 'scan' else direct_table(fb, panel)
        w, c = elapsed(start)
        arms[name] = {'status': 'complete', 'wall_seconds': w, 'cpu_seconds': c}
        answers[name] = ans
    assert answers['scan'] == answers['table'], 'direct arms disagree'
    assert all(row['witness'] in ans for row, ans in zip(panel[256:], answers['scan'][256:]))
    best = (min(arms[x]['wall_seconds'] for x in order),
            min(arms[x]['cpu_seconds'] for x in order))
    qans, state, prefix, orbits, (w, c) = quotient_table(xs, panel, best)
    arms['torsion'] = {'status': state, 'wall_seconds': w, 'cpu_seconds': c,
                       'processed_targets': prefix, 'orbits': orbits}
    for arm in arms.values():
        arm['total_wall_seconds'] = setup_wall + arm['wall_seconds']
        arm['total_cpu_seconds'] = setup_cpu + arm['cpu_seconds']
    if qans is not None:
        assert qans == answers['scan'][:prefix], 'torsion candidate mismatch'
    counts = [len(a) for a in answers['scan']]
    rr = relation_rank(fb, panel, answers['scan'])
    return {'n': n, 'repetition': rep, 'factor_base_points': len(fb),
            'shared_setup_wall_seconds': setup_wall,
            'shared_setup_cpu_seconds': setup_cpu,
            'target_count': len(panel), 'random_relation_pairs': sum(counts[:256]),
            'planted_relation_pairs': sum(counts[256:]),
            'random_targets_with_relation': sum(bool(x) for x in counts[:256]),
            'independent_rank': rr, 'arms': arms,
            'crossover': state == 'complete' and rr > 0 and all(
                w < arms[x]['wall_seconds'] and c < arms[x]['cpu_seconds'] for x in order)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', nargs='*', type=int, default=[8,16,32,64,128,256,512,1024])
    ap.add_argument('--repetitions', type=int, default=3)
    args = ap.parse_args()
    allowed = {8,16,32,64,128,256,512,1024}
    assert set(args.n) <= allowed and args.repetitions in (1,2,3)
    source = Path(__file__).read_bytes()
    protocol = (HERE/'PROTOCOL.md').read_bytes()
    output = HERE/'results.json'
    if output.exists():
        data = json.loads(output.read_text())
        assert data['source_sha256'] == hashlib.sha256(source).hexdigest()
        assert data['protocol_sha256'] == hashlib.sha256(protocol).hexdigest()
    else:
        data = {'source_sha256': hashlib.sha256(source).hexdigest(),
                'protocol_sha256': hashlib.sha256(protocol).hexdigest(),
                'python': sys.version, 'platform': platform.platform(),
                'curve': {'p': P, 'a': 1, 'b': 0, 'order': 38224,
                          'G': G, 'Q': Q}, 'cells': []}
    seen = {(x['n'], x['repetition']) for x in data['cells']}
    for n in args.n:
        for rep in range(args.repetitions):
            if (n, rep) in seen:
                continue
            cell = run_cell(n, rep)
            data['cells'].append(cell)
            output.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
            print(json.dumps(cell, sort_keys=True), flush=True)
    by_n = defaultdict(list)
    for cell in data['cells']:
        by_n[cell['n']].append(cell)
    for n in sorted(by_n):
        cells = by_n[n]
        print(n, 'reps', len(cells), 'torsion_status', [x['arms']['torsion']['status'] for x in cells],
              'median_scan', median(x['arms']['scan']['wall_seconds'] for x in cells),
              'median_table', median(x['arms']['table']['wall_seconds'] for x in cells))


if __name__ == '__main__':
    main()
