#!/usr/bin/env python3
"""Exact arithmetic certificate for docs/ecc2k130-integral-symmetry.md.

Standard library only. Arithmetic checks support the stated mathematical
arguments; this is not a formal verification of Rosati, Torelli, or Howe.
"""
import json
from pathlib import Path


def prime64(n):
    """Deterministic Miller--Rabin on unsigned 64-bit integers."""
    assert 0 <= n < 2**64
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, shifts = n - 1, 0
    while d % 2 == 0:
        d //= 2
        shifts += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(shifts - 1):
            x = x*x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def certificate():
    n, genus = 131, 130
    assert prime64(n)
    order2 = next(k for k in range(1, n) if pow(2, k, n) == 1)
    assert order2 == 130
    # Tame Riemann--Hurwitz, including every possible quotient genus.
    possibilities = []
    for h in range(genus + 1):
        numerator = 2*genus - 2 - n*(2*h - 2)
        if numerator >= 0 and numerator % (n - 1) == 0:
            possibilities.append([h, numerator // (n - 1)])
    assert possibilities == [[0, 4]]
    assert not any(r % order2 == 0 for _, r in possibilities)

    # alpha^131 = a*alpha + b, alpha^2 + alpha + 2 = 0.
    a, b = 0, 1
    for _ in range(n):
        a, b = b - a, -2*a
    ell = 146505763881528721
    assert prime64(263) and prime64(ell)
    assert a == 263*ell == 38531015900842053623
    assert b == 8123678690673902378
    assert 2*b - a == -22283658519494248867
    assert ell % n == n - 1
    assert pow(ell % n, 2, n) == 1
    assert pow((-7) % ell, (ell - 1)//2, ell) == ell - 1

    # Work in F_ell[alpha]/(alpha^2+alpha+2), an actual quadratic field
    # because -7 is a nonsquare. Coordinates are (constant, alpha).
    def mul(x, y):
        c, d = x
        e, f = y
        return ((c*e - 2*d*f) % ell, (c*f + d*e - d*f) % ell)

    def power(x, k):
        out = (1, 0)
        while k:
            if k & 1:
                out = mul(out, x)
            x = mul(x, x)
            k >>= 1
        return out

    alpha, beta, one = (0, 1), (ell - 1, ell - 1), (1, 0)
    assert mul(alpha, beta) == (2, 0)
    ratio = mul(beta, power(alpha, ell*ell - 2))
    assert ratio != one and power(ratio, n) == one
    zeta = power(ratio, pow(2, -1, n))
    assert zeta != one and power(zeta, n) == one
    assert mul(zeta, zeta) == ratio
    assert power(alpha, ell) == beta
    assert power(zeta, ell) == power(zeta, n - 1)
    pi = mul(alpha, zeta)
    verschiebung = mul(beta, power(zeta, n - 1))
    assert pi == verschiebung == power(pi, ell)
    assert pi[1] == 0 and mul(pi, pi) == (2, 0)

    # Two reductions of O_K to F_263 that coincide on F,V but disagree
    # on zeta prove zeta is not in the universal ring R = Z[F,V].
    p, alpha1, alpha2 = 263, 139, 123
    assert all((x*x + x + 2) % p == 0 for x in (alpha1, alpha2))
    ratio263 = alpha1*pow(alpha2, -1, p) % p
    assert ratio263 != 1 and pow(ratio263, n, p) == 1
    zeta1 = next(x for x in range(2, p)
                 if pow(x, n, p) == 1 and ratio263*x % p != 1)
    zeta2 = ratio263*zeta1 % p
    assert zeta1 != zeta2 and zeta1 != 1 and zeta2 != 1
    f1, f2 = alpha1*zeta1 % p, alpha2*zeta2 % p
    v1 = alpha2*pow(zeta1, -1, p) % p
    v2 = alpha1*pow(zeta2, -1, p) % p
    assert (f1, v1) == (f2, v2)

    # Finish the finite point/place/effective-divisor checks with rigorous
    # tails. These are consistency conditions, never Jacobian sufficiency.
    here = Path(__file__).parent
    coeffs = json.loads((here / "P_A.json").read_text())["P_A_coeffs_low_to_high"]
    remainder = [0]*263
    remainder[262], remainder[131], remainder[0] = 1, a - 2*b, 2**131
    quotient = [0]*261
    for k in range(262, 1, -1):
        t = remainder[k]
        quotient[k - 2] = t
        for j, c in enumerate((2, 1, 1)):
            remainder[k - 2 + j] -= t*c
    assert remainder == [0]*263 and coeffs == quotient
    traces = [2, -1]
    for k in range(2, 32):
        traces.append(-traces[-1] - 2*traces[-2])
    points = {k: 2**k + 1 - (130 if k % 131 == 0 else -1)*traces[k]
              for k in range(1, 32)}
    # Euler exponents, via recursive subtraction, avoid a second Mobius
    # implementation. Integrality for all degrees follows from the unique
    # Euler product of an integer power series with constant term one.
    places = {}
    for k in range(1, 32):
        numerator = points[k] - sum(d*places[d] for d in places if k % d == 0)
        assert numerator % k == 0
        places[k] = numerator // k
    assert min(points.values()) >= 0 and min(places.values()) >= 0
    # For d>=32 the negative terms in the crude Weil/Mobius tail,
    # divided by 2^d, are bounded by the following expression at d=32.
    # Every term decreases thereafter. We use d rather than floor(d/2)
    # as the number of possible proper divisors, deliberately overcounting.
    from fractions import Fraction
    tail_bound = (Fraction(260 + 32, 2**16) + Fraction(32, 2**32)
                  + Fraction(260*32, 2**24))
    assert tail_bound < 1
    # Thus N_d - sum_{e|d,e<d} N_e > 0 for all d>=32.
    # Also N_k>0 for k>=32 directly from the Weil bound.
    zeta_numerator = list(reversed(coeffs))
    effective = [sum(zeta_numerator[k]*(2**(m - k + 1) - 1)
                     for k in range(min(m, 260) + 1)) for m in range(261)]
    assert min(effective) >= 0
    order = sum(coeffs)
    assert sum(zeta_numerator) == order
    assert sum(c*2**(260-k) for k, c in enumerate(zeta_numerator)) == order*2**130
    assert all(effective[m] == order*(2**(m - 129) - 1) for m in (259, 260))

    return {
        "schema": "ecc2k130.integral-symmetry.v1",
        "order_131_of_2": order2,
        "riemann_hurwitz_quotient_genus_and_branch_count": possibilities,
        "rational_order_131_action_at_genus_130": "excluded by branch orbit lengths",
        "alpha_power_131": {"alpha_coefficient": a, "constant": b,
                            "alpha_coefficient_prime_factors": [263, ell]},
        "howe_theorem_1_1_witness": {
            "ell": ell, "ell_mod_131": ell % n,
            "minus_7_is_nonsquare": True,
            "zeta_coordinates": list(zeta), "pi_equals_V_coordinates": list(pi),
            "interpretation": "inert prime of K+ dividing F-V",
            "conclusion": "the F2-isogeny class contains a principally polarized variety",
        },
        "zeta_not_in_Z_F_V_witness": {
            "modulus": p, "alpha_images": [alpha1, alpha2],
            "zeta_images": [zeta1, zeta2], "common_F_V_images": [f1, v1],
        },
        "all_degree_zeta_consistency": {
            "exact_point_and_place_degrees": [1, 31],
            "tail_start": 32,
            "normalized_negative_tail_upper_bound_at_32": str(tail_bound),
            "effective_divisor_exact_degrees": [0, 258],
            "effective_divisor_tail": "B_m = P_A(1)*(2^(m-129)-1), m>=259",
            "conclusion": "all point, place, and effective-divisor counts are nonnegative integers",
        },
        "jacobian_in_isogeny_class": "unresolved",
    }


if __name__ == "__main__":
    result = certificate()
    path = Path(__file__).with_name("integral_symmetry.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
