#!/usr/bin/env python3
"""Replay sharp-cutoff checks and bounded rational/exceptional-prime probes."""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_curve_sweep as sweep

RECORD = ROOT / 'experiments/sharp-prime-bound-20261115/regression.json'


def compute():
    matrices = []
    for k in range(5):
        n = 2**k
        shift = sp.Matrix(n, n, lambda i, j: int((i+1) % n == j))
        mat = (shift+8*sp.eye(n))**n-(shift-sp.eye(n))**n
        assert all(a >= 0 for a in mat)
        assert all(sum(mat[i, j] for i in range(n)) == 9**n for j in range(n))
        det = int(mat.det(method='domain-ge'))
        sharp, old = 9**(n*n), 2**n*9**(n*n)
        assert 0 < abs(det) <= sharp < old
        matrices.append(dict(k=k, determinant=str(det), sharp=str(sharp), old=str(old)))

    # Denominator-cleared cyclic eliminants for a=m/q. These are numerical
    # regressions, not an all-parameter nonvanishing or LFD theorem.
    rational = []
    parameters = ['-3', '-1', '1', '-3/4', '-5/2', '-2', '-1/2', '0', '1/2', '2', '5/2', '6']
    x = sp.symbols('x')
    for text in parameters:
        a = Fraction(text)
        m, q = a.numerator, a.denominator
        for k in (2, 3):
            n = 2**k
            det = int(sp.resultant(x**n-1, (q*x+m+2*q)**n-(q*(x-1))**n, x))
            assert (det == 0) == (a in (-3, -1, 1))
            tested = 0
            for p in (17, 41, 97, 193, 257):
                if (p-1) % n or q % p == 0:
                    continue
                roots = [r for r in range(1, p) if pow(r, n, p) == 1]
                pairs = [(u, v) for u in roots for v in roots
                         if (q*(u-1)*v-q*u-m-2*q) % p == 0]
                assert bool(pairs) == (det % p == 0)
                tested += 1
            rational.append(dict(a=text, k=k, determinant=str(det), prime_cells=tested))

    exceptional = []
    for k, p in [(2, 17), (3, 257), (4, 65537), (5, 641), (5, 6700417), (6, 274177)]:
        assert sp.isprime(p) and pow(2, 2**k, p) == p-1
        curve = sweep.I.PrimeCurve(p, 6, -7)
        generator = sweep.prime_coset_generator(curve, 1, k, 1, 1)
        gv = sweep.subgroup_generator_value(p, k)
        points = sweep.prime_points(curve, k, gv, 1, 1, 1)
        profile = sweep.prime_local_profile(generator, k, p, points, gv,
                    max_degree=5, max_columns=16384, max_cells=30_000_000)
        assert profile['status'] == 'exact-proper-local'
        assert profile['containing_spaces_checked']
        exceptional.append(dict(k=k, p=p, profile=profile))
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                matrices=matrices, rational=rational, exceptional=exceptional,
                scope='sharp bound is proved in Lean; rational and exceptional profiles are finite regressions only')


def main():
    result = compute()
    if '--write' in sys.argv:
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECORD.read_text()), 'record differs from replay'
    print('PASS sharp cutoff:', len(result['matrices']), 'matrices;', len(result['rational']),
          'rational-depth cells;', len(result['exceptional']), 'exact exceptional profiles')


if __name__ == '__main__':
    main()
