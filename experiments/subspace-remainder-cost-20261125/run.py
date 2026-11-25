"""Frozen matched native LRF versus conditional-image cost comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
from pathlib import Path
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402
from csd.ic_rootfind import Solver, _lib  # noqa: E402
from make_manifest import FIXTURES, cell  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE/'manifest.json'
PROTOCOL = HERE/'PROTOCOL.md'
OUT = HERE/'results.json'


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def manifest():
    data = json.loads(MANIFEST.read_text())
    assert data['protocol_sha256'] == sha(PROTOCOL.read_bytes())
    assert data['source_sha256'] == sha((HERE/'make_manifest.py').read_bytes())
    return data


def now():
    return perf_counter(), process_time()


def elapsed(t):
    return perf_counter() - t[0], process_time() - t[1]


def pair_list(rows):
    return [[list(p) for p in pair] for pair in rows]


def worker(n, k, arm, bound_wall, bound_cpu):
    fixture = next(c for c in manifest()['cells'] if (c['n'], c['k']) == (n, k))
    start = now()
    F = GF2n(n, fixture['modulus'])
    E = BinaryCurve(F, 0, fixture['b'])
    S = Solver(E, [1 << i for i in range(k)])
    fb = sorted((int(x), int(y)) for x, y in zip(S.fx, S.fy))
    assert len(fb) == fixture['factor_base_points']
    assert sha(json.dumps(fb, separators=(',', ':')).encode()) == fixture['factor_base_sha256']
    setup_wall, setup_cpu = elapsed(start)
    rows = []
    for idx, target in enumerate(fixture['targets']):
        R = tuple(target['target']) if target['target'] is not None else None
        assert R is not None, 'frozen panel requires affine targets'
        t = now()
        if arm == 'lrf':
            pairs, stats, _ = S.lrf(2, R)
        elif arm == 'brute':
            pairs, _ = S.brute(2, R)
            stats = {}
        elif arm == 'image':
            pairs, stats, _ = S.lrf_image(2, R)
        else:
            raise ValueError(arm)
        assert all(E.add(a, b) == R for a, b in pairs)
        assert all(a in fb and b in fb for a, b in pairs)
        if target['kind'] == 'planted':
            assert target['witness'] in pair_list(pairs)
        wall, cpu = elapsed(t)
        rows.append({'index': idx, 'wall_seconds': wall, 'cpu_seconds': cpu,
                     'pairs': pair_list(pairs), 'stats': stats})
        total_wall, total_cpu = elapsed(start)
        if arm == 'image' and bound_wall is not None and bound_cpu is not None:
            if total_wall > bound_wall and total_cpu > bound_cpu:
                status = 'dominated_censored'
                break
    else:
        status = 'complete'
    total_wall, total_cpu = elapsed(start)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'n': n, 'k': k, 'arm': arm, 'status': status,
            'setup_wall_seconds': setup_wall, 'setup_cpu_seconds': setup_cpu,
            'wall_seconds': total_wall, 'cpu_seconds': total_cpu,
            'peak_rss_raw': rss, 'peak_rss_unit': 'bytes' if sys.platform == 'darwin' else 'KiB',
            'image_table_bytes': getattr(S, 'image_table_bytes', 0),
            'processed_targets': len(rows), 'targets': rows}


def rank_mod(rows, q):
    pivots = {}
    for raw in rows:
        v = {j: x % q for j, x in raw.items() if x % q}
        while v:
            j = min(v)
            if j not in pivots:
                inv = pow(v[j], -1, q)
                pivots[j] = {a: b*inv % q for a, b in v.items()}
                break
            mul = v[j]
            for a, b in pivots[j].items():
                x = (v.get(a, 0) - mul*b) % q
                if x:
                    v[a] = x
                else:
                    v.pop(a, None)
    return len(pivots)


def relation_rank(fixture, reference):
    E = BinaryCurve(GF2n(fixture['n'], fixture['modulus']), 0, fixture['b'])
    h, q = fixture['cofactor'], fixture['q']
    seen = set()
    for row in reference['targets'][:64]:
        for a, b in row['pairs']:
            for p in (tuple(a), tuple(b)):
                if E.mul(h, p) is not None:
                    seen.add(min(p, E.neg(p)))
    columns = {p: j for j, p in enumerate(sorted(seen))}
    vectors = []
    for target, row in zip(fixture['targets'][:64], reference['targets'][:64]):
        assert target['kind'] == 'random'
        for a, b in row['pairs']:
            v = {}
            for p0 in (a, b):
                p = tuple(p0)
                if E.mul(h, p) is None:
                    continue
                rep = min(p, E.neg(p))
                j = columns[rep]
                v[j] = v.get(j, 0) + (1 if p == rep else -1)
            v[len(columns)] = -h*(target['u'] + 113*target['v'])
            vectors.append(v)
    return rank_mod(vectors, q)


def child(n, k, arm, bound=None):
    cmd = [sys.executable, str(Path(__file__).resolve()), '--worker', str(n), str(k), arm]
    if bound is not None:
        cmd.extend([str(bound[0]), str(bound[1])])
    p = subprocess.run(cmd, check=True, text=True, capture_output=True, cwd=ROOT)
    return json.loads(p.stdout)


def run_main():
    frozen = manifest()
    _lib()  # Common software build before timed workers.
    source_hashes = {name: sha((ROOT/name).read_bytes()) for name in
                     ('csd/_rootfind.c', 'csd/ic_rootfind.py')}
    source_hashes['experiments/subspace-remainder-cost-20261125/run.py'] = sha(Path(__file__).read_bytes())
    if OUT.exists():
        result = json.loads(OUT.read_text())
        assert result['source_sha256'] == source_hashes
        assert result['manifest_sha256'] == sha(MANIFEST.read_bytes())
    else:
        result = {'protocol_sha256': sha(PROTOCOL.read_bytes()),
                  'manifest_sha256': sha(MANIFEST.read_bytes()),
                  'source_sha256': source_hashes,
                  'host': {'platform': platform.platform(), 'machine': platform.machine(),
                           'python': sys.version, 'compiler': os.popen('cc --version').read().splitlines()[:1]},
                  'cells': []}
    seen = {(x['n'], x['k'], x['repetition']) for x in result['cells']}
    cached_g = {}
    for args in FIXTURES:
        n, k = args[:2]
        fixture = next(c for c in frozen['cells'] if (c['n'], c['k']) == (n, k))
        t = now()
        generated = cell(*args, cached_g)
        assert generated == fixture
        common_wall, common_cpu = elapsed(t)
        for rep in range(3):
            if (n, k, rep) in seen:
                continue
            order = ('lrf', 'brute') if rep % 2 == 0 else ('brute', 'lrf')
            arms = {name: child(n, k, name) for name in order}
            assert arms['lrf']['status'] == arms['brute']['status'] == 'complete'
            assert [x['pairs'] for x in arms['lrf']['targets']] == [x['pairs'] for x in arms['brute']['targets']]
            bound = (min(arms[x]['wall_seconds'] for x in order),
                     min(arms[x]['cpu_seconds'] for x in order))
            arms['image'] = child(n, k, 'image', bound)
            prefix = arms['image']['processed_targets']
            assert [x['pairs'] for x in arms['image']['targets']] == [
                x['pairs'] for x in arms['lrf']['targets'][:prefix]]
            rr = relation_rank(fixture, arms['lrf'])
            for arm in arms.values():
                arm['total_wall_seconds'] = arm['wall_seconds'] + common_wall
                arm['total_cpu_seconds'] = arm['cpu_seconds'] + common_cpu
                arm['cost_per_rank_wall_seconds'] = (
                    arm['total_wall_seconds']/rr if rr and arm['status'] == 'complete' else None)
            image = arms['image']
            crossover = image['status'] == 'complete' and rr > 0 and all(
                image['wall_seconds'] < arms[x]['wall_seconds'] and
                image['cpu_seconds'] < arms[x]['cpu_seconds'] for x in order)
            row = {'n': n, 'k': k, 'repetition': rep, 'common_setup_wall_seconds': common_wall,
                   'common_setup_cpu_seconds': common_cpu, 'independent_rank': rr,
                   'crossover': crossover, 'arms': arms}
            result['cells'].append(row)
            OUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
            print(n, k, rep, 'rank', rr, 'image', image['status'],
                  'crossover', crossover, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', nargs='+')
    args = ap.parse_args()
    if args.worker:
        n, k, arm = int(args.worker[0]), int(args.worker[1]), args.worker[2]
        bound = tuple(map(float, args.worker[3:5])) if len(args.worker) == 5 else (None, None)
        print(json.dumps(worker(n, k, arm, *bound), separators=(',', ':')))
    else:
        run_main()


if __name__ == '__main__':
    main()
