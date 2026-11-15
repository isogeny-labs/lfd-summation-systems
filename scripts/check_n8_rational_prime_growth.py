#!/usr/bin/env python3
"""Exact finite checks of denominator-cleared rational prime-field growth inputs."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys
from types import SimpleNamespace

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd.ic_dreg import s3_prime_value

RECORD = ROOT / 'experiments/rational-prime-growth-20261115/regression.json'


def cutoff(k, m, q):
    n = 2**k
    height = 3*q+abs(m+2*q)
    return 2**n * height**(n*n) + 2*n + q + abs(m+3*q) + abs(4*m+3*q)


def roots(p, n):
    u = next(a for a in range(2, p) if pow(a, (p-1)//2, p) == p-1)
    g = pow(u, (p-1)//n, p)
    h = sorted({pow(g, i, p) for i in range(n)})
    assert len(h) == n and all(pow(x, n, p) == 1 for x in h)
    return h


def compute():
    rows, positives = [], []
    grid_checks = denominator_rejections = wrong_constant_disagreements = 0
    x = sp.symbols('x')
    parameters = ['-5/2', '-2', '-1/2', '0', '1/2', '2', '6', '1/17', '-3', '-1', '1', '-3/4']
    for text in parameters:
        a = Fraction(text)
        m, q = a.numerator, a.denominator
        for k in (2, 3):
            n = 2**k
            shift = sp.Matrix(n, n, lambda i, j: int((i+1) % n == j))
            matrix = (q*shift+(m+2*q)*sp.eye(n))**n-(q*(shift-sp.eye(n)))**n
            det = int(matrix.det(method='domain-ge'))
            assert det == int(sp.resultant(x**n-1, (q*x+m+2*q)**n-(q*(x-1))**n, x))
            assert (det == 0) == (a in (-3, -1, 1))
            height = 3*q+abs(m+2*q)
            assert abs(det) <= 2**n*height**(n*n) < cutoff(k, m, q)
            cells = []
            for p in (17, 41, 97, 193, 257):
                if (p-1) % n:
                    continue
                if q % p == 0:
                    denominator_rejections += 1
                    try:
                        pow(q, -1, p)
                    except ValueError:
                        continue
                    raise AssertionError('noninvertible denominator accepted')
                h = roots(p, n)
                ap = m*pow(q, -1, p) % p
                # Include singular reductions as polynomial controls. They are
                # excluded by the theorem's discriminant guards, not curves.
                curve = SimpleNamespace(p=p, a=ap, b=-1-ap)
                solutions = 0
                for u in h:
                    for v in h:
                        raw = s3_prime_value(curve, u, v, 1)
                        line = (q*(u-1)*v-q*u-m-2*q) % p
                        assert q*q*raw % p == line*line % p
                        solutions += raw == 0
                        mutant = (q*(u-1)*v-q*u-m-q) % p
                        wrong_constant_disagreements += mutant*mutant % p != q*q*raw % p
                        grid_checks += 1
                assert bool(solutions) == (det % p == 0)
                cells.append(dict(p=p, solutions=solutions,
                    nonsingular=(4*ap**3+27*(-1-ap)**2) % p != 0))
            rows.append(dict(a=text, k=k, determinant=str(det), cutoff=str(cutoff(k,m,q)), cells=cells))
        if a not in (-3, Fraction(-3,4), -1, 1) and q <= 2:
            bound = cutoff(2, m, q)
            p = bound+1+(1-(bound+1)) % 8
            while not sp.isprime(p):
                p += 8
            assert bound < p < 2**64 and (p-1) % 8 == 0
            assert q % p and (m+3*q) % p and (4*m+3*q) % p
            h = roots(p, 4)
            assert all((q*(u-1)*v-q*u-m-2*q) % p for u in h for v in h)
            positives.append(dict(a=text, k=2, p=p, cutoff=str(bound)))
    assert denominator_rejections and wrong_constant_disagreements
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        rows=rows, above_cutoff_units=positives, grid_checks=grid_checks,
        denominator_rejections=denominator_rejections,
        wrong_constant_disagreements=wrong_constant_disagreements,
        scope='finite eliminant, denominator and actual S3 checks; k>=24 growth is symbolic Lean evidence')


def main():
    result = compute()
    if '--write' in sys.argv:
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECORD.read_text()), 'record differs from replay'
    print('PASS rational prime growth:', len(result['rows']), 'parameter-depth cells;',
          result['grid_checks'], 'S3 checks;', len(result['above_cutoff_units']), 'above-cutoff units;',
          result['denominator_rejections'], 'denominator rejections')


if __name__ == '__main__':
    main()
