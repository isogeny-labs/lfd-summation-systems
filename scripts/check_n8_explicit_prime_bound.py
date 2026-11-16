#!/usr/bin/env python3
"""Exact finite regressions for GrowthExplicitPrimeBound; no large-depth closure run."""
import hashlib
import json
import sys
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / 'experiments/explicit-prime-bound-20261115/regression.json'


def compute():
    x = s.symbols('x')
    rows, eigenchecks = [], 0
    for k in range(5):
        n = 2**k
        shift = s.Matrix(n, n, lambda i, j: int((i+1) % n == j))
        identity = s.eye(n)
        mat = (shift+8*identity)**n - (shift-identity)**n
        det = int(mat.det(method='domain-ge'))
        bound = 2**n * 9**(4**k)
        assert mat.applyfunc(lambda a: a % 2) == identity
        assert det % 2 == 1 and 0 < abs(det) <= bound
        assert all(sum(abs(mat[i, j]) for i in range(n)) <= 9**n+2**n for j in range(n))
        eliminant = (x+8)**n-(x-1)**n
        assert det == int(s.resultant(x**n-1, eliminant, x))
        # Deliberately incorrect constant: the mod-two unit certificate must fail.
        mutant = (shift+7*identity)**n - (shift-identity)**n
        assert mutant.applyfunc(lambda a: a % 2) != identity
        cells = []
        for p in (17, 41, 97, 193, 257):
            if (p-1) % n:
                continue
            roots = [a for a in range(1, p) if pow(a, n, p) == 1]
            pairs = []
            for a in roots:
                vec = s.Matrix([pow(a, i, p) for i in range(n)])
                scalar = (pow(a+8, n, p)-pow(a-1, n, p)) % p
                assert (mat*vec-scalar*vec).applyfunc(lambda b: b % p) == s.zeros(n, 1)
                eigenchecks += 1
                for b in roots:
                    if ((a-1)*b-a-8) % p == 0:
                        pairs.append([a, b])
            assert bool(pairs) == (det % p == 0)
            cells.append(dict(p=p, solutions=len(pairs)))
        rows.append(dict(k=k, N=n, determinant=str(det), bound=str(bound), cells=cells))

    # A real positive instance of p>B_k at depth 2; primality is deterministic
    # in SymPy's <2^64 range. This checks the unit criterion, not the k>=24 window.
    n, bound = 4, 2**4 * 9**16
    p = bound+1 + (1-(bound+1)) % 8
    while not s.isprime(p):
        p += 8
    assert bound < p < 2**64 and (p-1) % 8 == 0
    nonresidue = next(a for a in range(2, p) if pow(a, (p-1)//2, p) == p-1)
    g = pow(nonresidue, (p-1)//n, p)
    roots = {pow(g, i, p) for i in range(n)}
    assert len(roots) == n
    assert not [(a, b) for a in roots for b in roots if ((a-1)*b-a-8) % p == 0]
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                rows=rows, eigenvector_checks=eigenchecks,
                above_threshold_unit_example=dict(k=2, p=p, bound=str(bound)),
                scope='exact finite matrix/unit regressions; the all-depth growth theorem is checked in Lean')


def main():
    result = compute()
    if '--write' in sys.argv[1:]:
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECORD.read_text()), 'recomputed record differs'
    print('PASS explicit prime bound: five matrices,', result['eigenvector_checks'],
          'eigenvectors, five negative controls, above-threshold unit at p =',
          result['above_threshold_unit_example']['p'])


if __name__ == '__main__':
    main()
