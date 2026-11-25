"""Exact Boolean-CNF lift with bounded all-model enumeration and matched arms."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from time import monotonic, perf_counter, process_time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd import ic_dreg as I  # noqa: E402
from csd import ic_rootfind as rf  # noqa: E402
from csd.gf2_ghs import GF2n  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'manifest.json'
OUT = HERE / 'results.json'
ARMS = ('sat', 'lrf', 'image', 'brute')
CAP = 10.0


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def manifest():
    data = json.loads(MANIFEST.read_text())
    assert data['protocol_sha256'] == sha((HERE / 'PROTOCOL.md').read_bytes())
    assert data['generator_sha256'] == sha((HERE / 'make_manifest.py').read_bytes())
    assert data['parent_sha256'] == sha((HERE.parent / 'subspace-remainder-cost-20261125' / 'manifest.json').read_bytes())
    return data


def canonical(rows):
    return sorted({tuple(sorted(tuple(map(int, p)) for p in pair)) for pair in rows})


def jsonify(rows):
    return [[list(p) for p in pair] for pair in rows]


def add_xor(clauses: list[list[int]], terms: list[int], parity: int, next_id: int) -> int:
    if not terms:
        if parity:
            clauses.append([])
        return next_id
    current = terms[0]
    for var in terms[1:]:
        next_id += 1
        out = next_id
        clauses += [[-current, -var, -out], [current, var, -out],
                    [current, -var, out], [-current, var, out]]
        current = out
    clauses.append([current if parity else -current])
    return next_id


def affine_system(E, S, target):
    F = E.F
    xT, yT = target
    assert xT != 0
    c1 = E.a ^ xT
    c2 = F.mul(F.sq(yT ^ xT) ^ E.b, F.inv(xT))
    Q = lambda x: rf._apply(S.Q, x)
    pi = lambda x: rf._apply(S.pi, x)
    top = S.w - 1
    cols = [pi(F.mul(xT, Q(v) ^ v)) | (F.trace(v) << top) for v in S.basis]
    cols += [pi(xT)]
    rhs = pi(F.mul(xT, Q(c1) ^ c1) ^ c2) | (F.trace(c1) << top)
    solution = I._f2_solve(cols, rhs, S.w)
    if solution is None:
        return None
    part, kernel = solution

    def at(mask):
        s = part
        for j, vector in enumerate(kernel):
            if mask >> j & 1:
                s ^= vector
        e1 = 0
        for i, v in enumerate(S.basis):
            if s >> i & 1:
                e1 ^= v
        lam = Q(e1 ^ c1) ^ ((s >> S.k) & 1)
        e2 = F.mul(xT, F.sq(lam)) ^ c2
        return e1, e2, lam

    zero = at(0)
    directions = [tuple(a ^ b for a, b in zip(at(1 << j), zero))
                  for j in range(len(kernel))]
    return zero, directions, at


def cnf_system(F, basis, zero, directions):
    k, d = len(basis), len(directions)
    next_id = k + d
    terms = [(j + 1, directions[j][1]) for j in range(d)]
    terms += [(d + i + 1, F.sq(basis[i]) ^ F.mul(zero[0], basis[i]))
              for i in range(k)]
    clauses = []
    for j in range(d):
        for i, v in enumerate(basis):
            coefficient = F.mul(directions[j][0], v)
            if coefficient == 0:
                continue
            next_id += 1
            h, g, z = next_id, j + 1, d + i + 1
            clauses += [[-h, g], [-h, z], [h, -g, -z]]
            terms.append((h, coefficient))
    for bit in range(F.N):
        row = [var for var, coeff in terms if coeff >> bit & 1]
        next_id = add_xor(clauses, row, (zero[1] >> bit) & 1, next_id)
    return clauses, next_id, d + k


def kissat_models(clauses, nvars, original_vars, deadline):
    models = []
    with tempfile.TemporaryDirectory(prefix='quad-lift-') as temp:
        cnf = Path(temp) / 'system.cnf'
        while True:
            left = deadline - monotonic()
            if left <= 0:
                return models, 'censored'
            with cnf.open('w') as stream:
                stream.write(f'p cnf {nvars} {len(clauses)}\n')
                for clause in clauses:
                    stream.write(' '.join(map(str, clause)) + ' 0\n')
            try:
                run = subprocess.run(['kissat', '-q', str(cnf)], text=True,
                                     capture_output=True, timeout=left)
            except subprocess.TimeoutExpired:
                return models, 'censored'
            if run.returncode == 20:
                return models, 'complete'
            if run.returncode != 10:
                raise RuntimeError(f'kissat exit {run.returncode}: {run.stderr[:500]}')
            assignment = {}
            for line in run.stdout.splitlines():
                if line.startswith('v '):
                    for literal in line.split()[1:]:
                        value = int(literal)
                        if value:
                            assignment[abs(value)] = value > 0
            if any(i not in assignment for i in range(1, original_vars + 1)):
                raise AssertionError('incomplete SAT model')
            bits = tuple(assignment[i] for i in range(1, original_vars + 1))
            models.append(bits)
            clauses.append([(-i if bits[i-1] else i) for i in range(1, original_vars + 1)])


def solve_sat(E, S, target, deadline):
    system = affine_system(E, S, target)
    if system is None:
        return [], 'complete', {'linear_status': 'inconsistent', 'models': 0}
    zero, directions, at = system
    clauses, nvars, original_vars = cnf_system(E.F, S.basis, zero, directions)
    models, status = kissat_models(clauses, nvars, original_vars, deadline)
    fb = {tuple(map(int, p)) for p in zip(S.fx, S.fy)}
    rows = []
    for bits in models:
        g = sum(1 << j for j, bit in enumerate(bits[:len(directions)]) if bit)
        x1 = 0
        for i, bit in enumerate(bits[len(directions):]):
            if bit:
                x1 ^= S.basis[i]
        e1, e2, lam = at(g)
        assert E.F.sq(x1) ^ E.F.mul(e1, x1) ^ e2 == 0
        x2 = x1 ^ e1
        mu = E.F.mul(lam, target[0]) ^ target[1] ^ target[0]
        p1 = (x1, E.F.mul(lam, x1) ^ mu)
        p2 = (x2, E.F.mul(lam, x2) ^ mu)
        assert E.on_curve(p1) and E.on_curve(p2)
        assert E.add(p1, p2) == target
        assert p1 in fb and p2 in fb
        rows.append((p1, p2))
    return canonical(rows), status, {'linear_status': 'solvable',
                                    'kernel_dim': len(directions), 'models': len(models),
                                    'cnf_vars': nvars, 'cnf_clauses_initial': len(clauses)-len(models),
                                    'input_degree_bound': 2,
                                    'max_solver_degree': None}


def worker(cell, index, arm):
    start_wall, start_cpu = perf_counter(), process_time()
    child_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    F = GF2n(cell['n'], cell['modulus'])
    E = I.BinaryCurve(F, 0, cell['b'])
    S = rf.Solver(E, [1 << i for i in range(cell['k'])])
    fb = sorted((int(x), int(y)) for x, y in zip(S.fx, S.fy))
    assert sha(json.dumps(fb, separators=(',', ':')).encode()) == cell['factor_base_sha256']
    setup_wall = perf_counter() - start_wall
    target_row = cell['targets'][index]
    target = tuple(target_row['target'])
    assert target[0] != 0
    query_start = perf_counter()
    if arm == 'sat':
        pairs, status, stats = solve_sat(E, S, target, monotonic() + CAP)
    elif arm == 'lrf':
        pairs, stats, _ = S.lrf(2, target)
        status = 'complete'
    elif arm == 'image':
        pairs, stats, _ = S.lrf_image(2, target)
        status = 'complete'
    else:
        pairs, _ = S.brute(2, target)
        stats, status = {}, 'complete'
    pairs = canonical(pairs)
    if status == 'complete':
        assert all(E.add(a, b) == target and a in fb and b in fb for a, b in pairs)
        if 'witness' in target_row:
            assert canonical([target_row['witness']])[0] in pairs
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    return {'cell': [cell['n'], cell['k']], 'index': index, 'kind': target_row['kind'],
            'arm': arm, 'status': status, 'pairs': jsonify(pairs), 'stats': stats,
            'setup_wall_seconds': setup_wall, 'query_wall_seconds': perf_counter()-query_start,
            'wall_seconds_internal': perf_counter()-start_wall,
            'cpu_seconds_internal': process_time()-start_cpu + child.ru_utime-child_start.ru_utime
                                    + child.ru_stime-child_start.ru_stime,
            'peak_rss_raw': max(self_usage.ru_maxrss, child.ru_maxrss),
            'peak_rss_unit': 'bytes' if sys.platform == 'darwin' else 'KiB'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', nargs=4, metavar=('N', 'K', 'INDEX', 'ARM'))
    ap.add_argument('--pilot', action='store_true')
    args = ap.parse_args()
    data = manifest()
    if args.worker:
        n, k, index = map(int, args.worker[:3])
        arm = args.worker[3]
        assert arm in ARMS
        cell = next(c for c in data['cells'] if (c['n'],c['k']) == (n,k))
        print(json.dumps(worker(cell, index, arm), sort_keys=True), flush=True)
        return
    cases = [(c, i, arm) for c in data['cells'] for i in range(len(c['targets']))
             for arm in ARMS]
    if args.pilot:
        cases = [(data['cells'][0], 4, arm) for arm in ARMS]
    receipts = []
    for c, i, arm in cases:
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   str(c['n']), str(c['k']), str(i), arm]
        start = perf_counter()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        try:
            process = subprocess.run(command, capture_output=True, text=True,
                                     timeout=CAP+5)
            if process.returncode:
                row = {'cell':[c['n'],c['k']], 'index':i, 'arm':arm,
                       'status':'error', 'stderr':process.stderr[-1000:],
                       'stdout':process.stdout[-1000:]}
            else:
                row = json.loads(process.stdout)
        except subprocess.TimeoutExpired:
            row = {'cell':[c['n'],c['k']], 'index':i, 'arm':arm, 'status':'censored_worker'}
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        row['charged_wall_seconds'] = perf_counter()-start
        row['charged_cpu_seconds'] = (after.ru_utime-before.ru_utime
                                      + after.ru_stime-before.ru_stime)
        receipts.append(row)
        print(c['n'],c['k'],i,arm,row['status'],round(row['charged_wall_seconds'],3), flush=True)
    if not args.pilot:
        kissat = shutil.which('kissat')
        assert kissat is not None
        version = subprocess.run([kissat, '--version'], capture_output=True,
                                 text=True, check=True).stdout.strip()
        environment = {'python':sys.version.split()[0], 'platform':platform.platform(),
                       'machine':platform.machine(), 'kissat_path':kissat,
                       'kissat_version':version,
                       'ic_rootfind_sha256':sha((ROOT/'csd'/'ic_rootfind.py').read_bytes()),
                       'rootfind_c_sha256':sha((ROOT/'csd'/'_rootfind.c').read_bytes()),
                       'ic_dreg_sha256':sha((ROOT/'csd'/'ic_dreg.py').read_bytes()),
                       'gf2_ghs_sha256':sha((ROOT/'csd'/'gf2_ghs.py').read_bytes())}
        OUT.write_text(json.dumps({'manifest_sha256':sha(MANIFEST.read_bytes()),
                                   'runner_sha256':sha(Path(__file__).read_bytes()),
                                   'environment':environment,
                                   'receipts':receipts}, sort_keys=True, indent=2)+'\n')
    else:
        print(json.dumps(receipts, sort_keys=True))


if __name__ == '__main__':
    main()
