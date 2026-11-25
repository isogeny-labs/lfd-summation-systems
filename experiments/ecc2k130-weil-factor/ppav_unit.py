"""An explicit CM unit correcting all signs of the trial polarization.

All sign checks use integers or rational Taylor enclosures, not floats.
"""
from fractions import Fraction as Q
from math import factorial, gcd

import flint
import sympy as s


def negative_m_embeddings():
    """Signs of -1/2+sqrt(7)/2*cot(2*pi*i/131), i=1,...,130."""
    def atan_bounds(inv, terms=45):
        val = sum((Q((-1)**k, (2*k+1)*inv**(2*k+1))
                   for k in range(terms)), Q(0))
        following = Q((-1)**terms, (2*terms+1)*inv**(2*terms+1))
        return min(val, val+following), max(val, val+following)

    a, b = atan_bounds(5)
    c, d = atan_bounds(239)
    # Machin's identity gives a rigorous enclosure of pi.
    lo = Q("3.14159265358979323846264338327950288419716939937510")
    hi = lo+Q(1, 10**50)
    assert lo < 16*a-4*d < 16*b-4*c < hi
    rlo = Q(264575131106459059050, 10**20)
    rhi = rlo+Q(1, 10**20)
    assert rlo*rlo < 7 < rhi*rhi
    negative = []
    for i in range(1, 131):
        j = i if i <= 65 else i-131
        center = (lo+hi)*j/131
        scaled = center*10**20
        q = Q(scaled.numerator//scaled.denominator, 10**20)
        error = Q(1, 10**18)
        input_error = abs(q-center)+2*(hi-lo)
        assert input_error+abs(q)**44/factorial(44) < error
        assert input_error+abs(q)**45/factorial(45) < error
        sine = sum(((-1)**k*q**(2*k+1)/factorial(2*k+1)
                    for k in range(22)), Q(0))
        cosine = sum(((-1)**k*q**(2*k)/factorial(2*k)
                      for k in range(23)), Q(0))
        slo, shi = sine-error, sine+error
        clo, chi = cosine-error, cosine+error
        nlo = min(rlo*clo, rlo*chi, rhi*clo, rhi*chi)-shi
        nhi = max(rlo*clo, rlo*chi, rhi*clo, rhi*chi)-slo
        assert slo > Q(1, 1000) or shi < -Q(1, 1000)
        assert nlo > Q(1, 1000) or nhi < -Q(1, 1000)
        if (slo > 0) != (nlo > 0):
            negative.append(i)
    return negative


def build_unit():
    conductor = 917

    def crt(r, i):
        return next(b for b in range(i, conductor, 131) if b % 7 == r)

    embeddings = [crt(1, i) for i in range(1, 131)]
    subgroup = [crt(h, 1) for h in (1, 2, 4)]

    def signature(a):
        out = 0
        for i, b in enumerate(embeddings):
            sign = 0
            for h in subgroup:
                sign ^= ((a*b*h) % conductor > conductor//2)
                sign ^= ((b*h) % conductor > conductor//2)
            out |= sign << i
        return out

    # Norm of (theta^a-theta^-a)/(theta-theta^-1) from the real
    # cyclotomic field to K+. Its three-factor signatures are exact.
    basis, units = {}, []
    for a in range(1, conductor//2+1):
        if gcd(a, conductor) != 1:
            continue
        val, combination = signature(a), 1 << len(units)
        units.append(a)
        while val:
            pivot = val.bit_length()-1
            if pivot not in basis:
                basis[pivot] = val, combination
                break
            val ^= basis[pivot][0]
            combination ^= basis[pivot][1]
    assert len(basis) == 129
    negative = negative_m_embeddings()
    assert len(negative) % 2 == 0
    target = sum(1 << (i-1) for i in negative)
    val, combination = target, 0
    while val:
        pivot = val.bit_length()-1
        assert pivot in basis
        val ^= basis[pivot][0]
        combination ^= basis[pivot][1]
    chosen = [a for j, a in enumerate(units) if combination >> j & 1]
    check = 0
    for a in chosen:
        check ^= signature(a)
    assert check == target

    x = s.Symbol("x")
    phi = flint.fmpz_poly(list(map(int, reversed(
        s.cyclotomic_poly(conductor, x, polys=True).all_coeffs()))))
    unit = flint.fmpz_poly([1])
    for a in chosen:
        for h in subgroup:
            eta = [0]*conductor
            for j in range(a):
                eta[((a-1-2*j)*h) % conductor] += 1
            unit = (unit*flint.fmpz_poly(eta)) % phi

    # theta = zeta_7^3*zeta_131^75. Convert the degree-780
    # representation to the 1,alpha basis over Z[zeta_131].
    converted = [[0]*130 for _ in range(6)]
    for k, coefficient in enumerate(unit):
        coefficient = int(coefficient)
        if not coefficient:
            continue
        r, j = 3*k % 7, 75*k % 131
        rterms = [(z, -1) for z in range(6)] if r == 6 else [(r, 1)]
        jterms = [(z, -1) for z in range(130)] if j == 130 else [(j, 1)]
        for rr, rs in rterms:
            for jj, js in jterms:
                converted[rr][jj] += rs*js*coefficient
    # alpha=zeta_7+zeta_7^2+zeta_7^4.
    assert converted[1] == converted[2] == converted[4]
    assert not any(converted[3]) and not any(converted[5])
    coeffs = [converted[0], converted[1]]
    return coeffs, {
        "cyclotomic_conductor": conductor,
        "norm_subgroup_exponents": subgroup,
        "chosen_circular_units": chosen,
        "signature_rank": len(basis),
        "negative_embeddings_of_m": negative,
        "unit_norm_Kplus_Q": 1,
        "positivity_certificate": "rational Taylor enclosures and exact circular-unit signs",
        "unit_coefficients_1_alpha_over_zeta": coeffs,
    }
