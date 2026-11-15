#!/usr/bin/env python3
"""Independent exact checks of fixed-curve units and rational parameters.

Run with --write to save the deterministic record; otherwise recompute and compare.
No closure degree, prime-density claim, or fresh Lean proof is inferred.
The subgroup scan imports no project generator or closure implementation.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / 'experiments/fixed-curve-validation-20261114/results.json'
X = sp.symbols('X')


def subgroup(p, n):
    # For n a power of two, g^(n/2) != 1 proves order exactly n.
    for a in range(2, p):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) != 1:
            h = {pow(g, j, p) for j in range(n)}
            assert len(h) == n and all(pow(x, n, p) == 1 for x in h)
            return h
    raise AssertionError((p, n))


def graph_pairs(p, h, offset=8):
    pairs = []
    for x in sorted(h):
        if x == 1:
            continue  # h(1,z)=-9; all scanned primes exceed 8.
        z = (x + offset) * pow(x - 1, -1, p) % p
        if z in h:
            pairs.append([x, z])
    return pairs


def raw_s3(x, z, p):
    return ((x-z)**2 - 2*((x+z)*(x*z+6)-14)
            + (x*z-6)**2 + 28*(x+z)) % p


def compute():
    primes = list(sp.primerange(9, 100000))
    rows, raw_checks, mutant_disagreements = [], 0, 0
    for k in range(2, 11):
        n, cells = 2**k, []
        eligible = [p for p in primes if p > 2*n and (p-1) % (2*n) == 0]
        for j, p in enumerate(eligible):
            h = subgroup(p, n)
            pairs = graph_pairs(p, h)
            for x, z in pairs:
                assert raw_s3(x, z, p) == 0
            if k <= 6 and j < 3:
                brute = [[x, z] for x in sorted(h) for z in sorted(h)
                         if raw_s3(x, z, p) == 0]
                assert pairs == brute, (k, p)
                raw_checks += n*n
                mutant_disagreements += graph_pairs(p, h, offset=7) != brute
            diagonal = sum(pow(a, n, p) == 1 for a in (4, p-2))
            assert sum(x == z for x, z in pairs) == diagonal
            assert (len(pairs)-diagonal) % 2 == 0
            cells.append(dict(p=int(p), solutions=len(pairs),
                diagonal_solutions=diagonal, offdiagonal_solutions=len(pairs)-diagonal,
                witness=pairs[:1]))
        bad = [c for c in cells if c['solutions']]
        good = [c for c in cells if not c['solutions']]
        rows.append(dict(k=k, N=n, cells=cells, checked=len(cells),
                         unit=len(good), proper=len(bad),
                         smallest_unit=good[0]['p'] if good else None))
    assert sum(r['checked'] for r in rows) == 4759
    assert sum(r['unit'] for r in rows) == 4636
    assert raw_checks == 16368 and mutant_disagreements > 0

    resultants = []
    for k in range(2, 5):
        n = 2**k
        f = sp.Poly(X**n-1, X)
        g = sp.Poly((X+8)**n-(X-1)**n, X)
        r = int(sp.resultant(f, g))
        assert r != 0 and abs(r) <= 2**n * 9**(n*n)
        diagonal_factor = (4**n-1)*(2**n-1)
        assert r % diagonal_factor == 0
        offdiagonal_factor = math.isqrt(r // diagonal_factor)
        assert r == diagonal_factor * offdiagonal_factor**2
        remainder, factors = abs(r), {}
        # Bounded trial division: never launch an unbounded factorization.
        for p in sp.primerange(2, 20001):
            e = 0
            while remainder % p == 0:
                remainder //= p
                e += 1
            if e:
                factors[str(p)] = e
        for cell in rows[k-2]['cells']:
            assert (r % cell['p'] == 0) == bool(cell['solutions'])
            if cell['p'] < 20000:
                p = cell['p']
                degree = int(sp.degree(sp.gcd(sp.Poly(f, X, modulus=p),
                                             sp.Poly(g, X, modulus=p))))
                assert degree == cell['solutions'], (k, p, degree)
        resultants.append(dict(k=k, resultant=str(r), small_factors=factors,
                               diagonal_factor=str(diagonal_factor),
                               offdiagonal_factor=str(offdiagonal_factor),
                               unfactored_cofactor=str(remainder),
                               factorization_complete=remainder == 1))
    assert resultants[0]['factorization_complete']
    assert resultants[1]['factorization_complete']

    # An independent integral multiplication matrix, rather than elimination.
    certificates = []
    for k in range(2, 5):
        n = 2**k
        m = sp.zeros(n*n)
        for i in range(n):
            for j in range(n):
                for u, w, c in [(1, 1, 1), (1, 0, -1), (0, 1, -1), (0, 0, -8)]:
                    m[((i+u) % n)*n + (j+w) % n, i*n+j] += c
        r = abs(int(resultants[k-2]['resultant']))
        if k == 2:
            assert int(m.det()) == r
        rhs = sp.zeros(n*n, 1)
        rhs[0] = 1
        sol, denominator = DomainMatrix.from_Matrix(m).solve_den(DomainMatrix.from_Matrix(rhs))
        v = sol.to_Matrix()
        common = int(sp.igcd(int(denominator), sp.igcd(*list(v))))
        d = int(denominator) // common
        v = v / common
        if d < 0:
            d, v = -d, -v
        assert m*v == d*rhs
        # Both divisibilities prove identical prime support without factorization.
        assert r % d == 0 and pow(d, n*n, r) == 0
        vv = [0]*(n*n)
        for i in range(n*n):
            for j in range(n*n):
                vv[((i//n+j//n) % n)*n + (i % n+j % n) % n] += int(v[i])*int(v[j])
        assert m*m*sp.Matrix(vv) == d*d*rhs
        certificates.append(dict(k=k, scalar=str(d),
            coefficients=[str(int(c)) for c in v], basis='x^i*z^j, index=N*i+j',
            same_prime_support_as_resultant=True))

    parameters = sorted({Fraction(a, b) for a in range(-12, 13) for b in range(1, 5)})
    rational_rows = []
    for k in range(2, 6):
        n = 2**k
        for a in parameters:
            aa = sp.Rational(a.numerator, a.denominator)
            degree = int(sp.degree(sp.gcd(sp.Poly(X**n-1, X),
                sp.Poly((X+2+aa)**n-(X-1)**n, X))))
            assert bool(degree) == (a in (Fraction(-3), Fraction(-1), Fraction(1)))
            rational_rows.append(dict(k=k, a=str(a), gcd_degree=degree))
    return dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Exact finite-instance validation; general arguments are written in README.md, not Lean',
                scan_bound_exclusive=100000, scan=rows, raw_pair_checks=raw_checks,
                deliberately_wrong_graph_disagreements=mutant_disagreements,
                resultants=resultants, integral_inverse_certificates=certificates,
                rational_parameters=rational_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = compute()
    if args.write:
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(RECORD.read_text()), 'record differs from recomputation'
    print(json.dumps(dict(status='passed', depth_prime_cells=sum(r['checked'] for r in result['scan']),
        unit_cells=sum(r['unit'] for r in result['scan']), raw_pair_checks=result['raw_pair_checks'],
        rational_checks=len(result['rational_parameters']),
        negative_control_disagreements=result['deliberately_wrong_graph_disagreements']), indent=2))


if __name__ == '__main__':
    main()
