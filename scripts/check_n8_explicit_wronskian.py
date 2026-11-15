#!/usr/bin/env python3
"""Exact Q[X] regression for N8's elementary polynomial-Wronskian route.

Checks the recurrence, degree ceiling, and all cofactor identities
for the actual 4-term rational relation x*z - x - z - 5/2 = 0.
This is a finite test, NOT a proof of GrowthSpread.SpreadBound.
"""
from fractions import Fraction as F
from itertools import permutations


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def poly(*a):
    return trim([F(x) for x in a])


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def scale(a, c):
    return trim([x * c for x in a])


def diff(a):
    return trim([F(i) * a[i] for i in range(1, len(a))]) if len(a) > 1 else poly(0)


def power(a, n):
    out = poly(1)
    for _ in range(n):
        out = mul(out, a)
    return out


def det(matrix):
    n = len(matrix)
    result = poly(0)
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j] for i in range(n)
                         for j in range(i + 1, n))
        term = poly(1)
        for i, j in enumerate(perm):
            term = mul(term, matrix[i][j])
        result = add(result, scale(term, (-1) ** inversions))
    return result


def main():
    c = F(5, 2)
    x = poly(0, 1)
    one = poly(1)
    xm1 = poly(-1, 1)
    xc = poly(c, 1)
    D = mul(mul(x, xm1), xc)
    # Relation: x*z - x - z - 5/2 = 0; z=(x+5/2)/(x-1).
    exponents = [(1, 1), (1, 0), (0, 1), (0, 0)]
    coefficients = [F(1), F(-1), F(-1), -c]
    r = len(exponents) - 1
    bound = r * (r - 1)
    ps = []
    for u, v in exponents:
        Q = add(sub(scale(mul(xm1, xc), u),
                    scale(mul(x, xc), v)),
                scale(mul(x, xm1), v))
        sequence = [one]
        for j in range(r - 1):
            sequence.append(add(mul(D, diff(sequence[-1])),
                                mul(sub(Q, scale(diff(D), j)),
                                    sequence[-1])))
        assert all(len(sequence[j]) - 1 <= 2 * j for j in range(r))
        ps.append(sequence)
    cofactors = []
    for omitted in range(len(exponents)):
        matrix = [[ps[i][j] for j in range(r)]
                  for i in range(len(exponents)) if i != omitted]
        result = det(matrix)
        assert result != poly(0) and len(result) - 1 <= bound
        cofactors.append(result)

    # W_j/W_0 = (-1)^j*coefficient_j/coefficient_0, using the
    # determinant cofactor identity for a minimal four-term relation.
    for j, (u, v) in enumerate(exponents):
        u0, v0 = exponents[0]
        fj_num = mul(power(x, u), power(xc, v))
        fj_den = power(xm1, v)
        f0_num = mul(power(x, u0), power(xc, v0))
        f0_den = power(xm1, v0)
        lhs = scale(mul(mul(f0_num, fj_den), cofactors[j]),
                    coefficients[0])
        rhs = scale(mul(mul(fj_num, f0_den), cofactors[0]),
                    (-1) ** j * coefficients[j])
        assert lhs == rhs, j

    assert [len(A) - 1 for A in cofactors] == [5, 5, 4, 4]
    print("PASS: Q[X] derivative numerator recurrence and degrees <= 6")
    print("PASS: all four exact Wronskian cofactor identities")


if __name__ == "__main__":
    main()
