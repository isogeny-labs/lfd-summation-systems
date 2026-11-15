#!/usr/bin/env python3
"""Independent exact regression for the N8 rational two-torsion coordinate.

Tests the Möbius involution and its simple-pole identity over prime fields,
and the characteristic-p derivative boundary on p-th-power coefficients.
This is evidence for the algebraic lemmas, not a proof of the Wronskian
height bound or CurveModel.
"""
from math import factorial


def test_mobius():
    primes = (5, 7, 11, 13, 17, 31, 101)
    count = 0
    for p in primes:
        for a in range(p):
            if (a + 3) % p == 0:
                continue
            for x in range(p):
                if x == 1:
                    continue
                z = ((x + 2 + a) * pow(x - 1, -1, p)) % p
                assert z != 1
                assert (x - 1) * (z - 1) % p == (a + 3) % p
                assert ((z + 2 + a) * pow(z - 1, -1, p)) % p == x
                count += 1
    assert count == 11452, count
    return count


def test_frobenius_coefficients():
    # The derivative of h^p vanishes in char p even for *nonconstant*
    # h=(X+2)/(X-1).  A p-circuit argument must use F^p, not just constants.
    for p in (3, 5, 7, 11, 17):
        # numerator/denominator of h^p; their derivatives each vanish mod p.
        numerator = [0] * (p + 1)
        denominator = [0] * (p + 1)
        for j in range(p + 1):
            numerator[j] = (factorial(p) // (factorial(j) * factorial(p - j))
                            * 2 ** (p - j)) % p
            denominator[j] = (factorial(p) // (factorial(j) * factorial(p - j))
                              * (-1) ** (p - j)) % p
        assert numerator[0] == 2 % p and numerator[p] == 1
        assert denominator[0] == (-1) % p and denominator[p] == 1
        for j in range(1, p + 1):
            assert (j * numerator[j]) % p == 0
            assert (j * denominator[j]) % p == 0
        # The p successive derivatives of 1,X,...,X^(p-1)
        # have nonzero triangular Wronskian determinant.
        det = 1
        for j in range(p):
            det = det * factorial(j) % p
        assert det != 0
    return 5


if __name__ == "__main__":
    print("N8 Möbius tests:", test_mobius(), "exact field tuples")
    print("N8 Frobenius/Wronskian boundary:", test_frobenius_coefficients(), "primes")
