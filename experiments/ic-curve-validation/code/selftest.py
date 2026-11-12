#!/usr/bin/env python3
"""Self-tests of the independent engine (lfd_engine.py): linear algebra vs exact Python
integers, variable-wise closure vs the literal Lean multiplier rule, and tiny known cases."""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lfd_engine as L

ROOT = '/home/user/ml-cryptanalysis'


def rank_py(M, p):
    M = [[x % p for x in row] for row in M]
    rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        pr = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if pr is None:
            continue
        M[rk], M[pr] = M[pr], M[rk]
        inv = pow(M[rk][c], -1, p)
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def test_rank():
    rng = random.Random(1)
    for p in (2, 3, 97, 257, 12289):
        for _ in range(12):
            n, m = rng.randint(1, 60), rng.randint(1, 60)
            r0 = rng.randint(0, min(n, m))
            A = np.array([[rng.randrange(p) for _ in range(m)] for _ in range(r0)] or [[0] * m])
            B = np.array([[rng.randrange(p) for _ in range(r0 or 1)] for _ in range(n)])
            M = (B @ A) % p if r0 else np.zeros((n, m), dtype=np.int64)
            assert L.rank_mod(M, p) == rank_py(M.tolist(), p), (p, n, m)
    print('rank_mod == pure-Python exact rank on random low-rank matrices: OK')


def test_tiny():
    # F = {x^2 - x}: principal ideal, LFD 0 (every ideal element is q*f with deg q + deg f = deg h)
    out = L.exact_lfd([{(2,): 1, (1,): 2}], 1, 3, [(0,), (1,)], 4)
    assert out['lfd'] == 0, out
    # F = {x^2-1, x^2-x}: ideal (x-1), zero set {1}; one-variable.
    out = L.exact_lfd([{(2,): 1, (0,): 2}, {(2,): 1, (1,): 2}], 1, 3, [(1,)], 4)
    print('toy {x^2-1,x^2-x} p=3: lfd', out['lfd'], 'falls', out['falls'], [ (r['D'], r['dimV'], r['dimI']) for r in out['per_degree']])
    # By hand: x-1 = (x^2-1)-(x^2-x)... = (x^2 - 1) - (x^2 - x) = x - 1, so x-1 in V_2 (deg 1 element of V_2),
    # V_1 = 0 while I_1 = span{x-1}: dim 1 > 0 -> fall at 2 -> LFD 2.
    assert out['lfd'] == 2 and out['falls'] == [2], out
    # unit ideal {x, x-1}... x-(x-1) = 1 in V_1: LFD 1
    out = L.exact_lfd([{(1,): 1}, {(1,): 1, (0,): 2}], 1, 3, [], 3)
    assert out['lfd'] == 1, out
    # unit ideal {x^2, x^2-1}: 1 in V_2: LFD 2
    out = L.exact_lfd([{(2,): 1}, {(2,): 1, (0,): 2}], 1, 3, [], 3)
    assert out['lfd'] == 2, out
    print('hand-computed toy LFDs (0, 2, 1, 2): OK')


def test_naive(nsys=6):
    S = {}
    for line in open(ROOT + '/experiments/ic-curve-sweep/systems.jsonl'):
        s = json.loads(line)
        S[s['system_hash']] = s
    pr = [s for s in S.values() if s['kind'] == 'prime' and s['k'] == 2]
    bi = [s for s in S.values() if s['kind'] == 'binary' and s['variables'] == 4]
    rng = random.Random(7)
    rng.shuffle(pr)
    rng.shuffle(bi)
    n = 0
    for s in pr[:nsys] + bi[:nsys]:
        if s['kind'] == 'prime':
            inputs, G = L.prime_inputs([(int(m), int(c)) for m, c in s['generators'][0]], s['k'], s['prime'])
            r, p = 2 * s['k'], s['prime']
        else:
            r, p = s['variables'], 2
            inputs = L.binary_inputs(s['equations'], r)
        for D in range(0, 6):
            a = L.closure(inputs, r, D, p)
            b = L.closure_naive(inputs, r, D, p)
            assert (a['dim'], a['dim_low'], a['contains_one']) == (b['dim'], b['dim_low'], b['contains_one']), (s['system_hash'], D, a['dim'], b['dim'])
        n += 1
    print('variable-wise closure == literal monomial-multiplier closure (D=0..5, r=4) on', n, 'systems: OK')


if __name__ == '__main__':
    test_rank()
    test_tiny()
    test_naive()
