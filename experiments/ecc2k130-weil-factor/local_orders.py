#!/usr/bin/env python3
"""Arithmetic inputs for the isogeny-degree obstruction, not a curve search.

Uses SymPy's exact modular polynomial arithmetic; python-flint is optional.
See docs/ecc2k130-isogeny-class-obstruction.md for the written proof.
"""
import json
from pathlib import Path

import sympy as s
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_pow_mod

from integral_symmetry import prime64


def main():
    here = Path(__file__).parent
    coeffs = json.loads((here / "P_A.json").read_text())["P_A_coeffs_low_to_high"]
    n, ell = 131, 146505763881528721
    a, b = 0, 1
    for _ in range(n):
        a, b = b-a, -2*a
    assert a == 263*ell and prime64(263) and prime64(ell)

    # alpha and zeta are both integral in R[1/a], including at 2:
    # alpha=(F^131-b)/a, zeta=V^130-t*F,
    # t=(beta^130-1)/alpha in Z[alpha].
    ua, va = 0, 1
    for _ in range(130):
        ua, va = va-ua, -2*ua
    # beta^130 = -ua*alpha + va-ua; division by alpha is integral.
    c, d = -ua, va-ua
    assert (d-1) % 2 == 0
    t0, t1 = c - (d-1)//2, -(d-1)//2
    # alpha*(t0+t1*alpha) = beta^130-1.
    assert (-2*t1, t0-t1) == (d-1, c)

    x = s.Symbol("x")
    X = lambda p: s.Poly(x, x, modulus=p)

    def frobenius(poly, exponent, p):
        return s.Poly.from_list(gf_pow_mod([1, 0], exponent,
                                list(map(int, poly.all_coeffs())), p, ZZ),
                                x, modulus=p)

    patterns = {}
    for p in (263, ell):
        f = s.Poly.from_list(list(reversed(coeffs)), x, modulus=p)
        repeated = s.gcd(f, f.diff())
        assert repeated.degree() == 129
        assert s.gcd(repeated, repeated.diff()).degree() == 0
        single = f.exquo(repeated*repeated)
        assert single.degree() == 2 and s.gcd(single, repeated).degree() == 0
        linear_repeated = s.gcd(repeated, frobenius(repeated, p, p)-X(p))
        linear_single = s.gcd(single, frobenius(single, p, p)-X(p))
        assert s.rem(frobenius(repeated, p*p, p)-X(p), repeated).is_zero
        assert s.rem(frobenius(single, p*p, p)-X(p), single).is_zero
        if p == 263:
            assert (linear_repeated.degree(), linear_single.degree()) == (129, 2)
        else:
            assert (linear_repeated.degree(), linear_single.degree()) == (1, 0)
        patterns[str(p)] = {
            "linear_factors_multiplicity_1": linear_single.degree(),
            "linear_factors_multiplicity_2": linear_repeated.degree(),
            "quadratic_factors_multiplicity_1": (2-linear_single.degree())//2,
            "quadratic_factors_multiplicity_2": (129-linear_repeated.degree())//2,
        }

    assert s.legendre_symbol(-7, 131) == -1
    assert s.legendre_symbol(-7, 263) == 1
    assert s.legendre_symbol(-7, ell) == -1
    assert 263 % 131 == 1 and ell % 131 == 130
    # A product polarization on the Weil restriction restricts to I+J
    # on E^130; the determinant lemma gives det(I+J)=131.
    assert 1+130 == n

    record = {
        "schema": "ecc2k130.local-orders.v1",
        "field_discriminant_factorization": {"7": 130, "131": 258},
        "power_order_index_factorization": {"2": 8385, "263": 129, str(ell): 129},
        "frobenius_verschiebung_order_index_factorization": {"263": 129, str(ell): 129},
        "zeta_reconstruction_t_coefficients_constant_alpha": [t0, t1],
        "modular_factor_patterns": patterns,
        "K_over_Kplus_at_131": "unramified inert, residue degree 2/1",
        "K_over_Kplus_at_263": "split at every prime",
        "K_over_Kplus_at_ell": "inert at all 130 primes of Kplus, residue degree 2/1",
        "natural_polarization_kernel": "E[131], one simple O_K module of length 1",
        "written_proof_consequence": {
            "for": "every F2-isogeny f:A->B with B principally polarized over F2",
            "condition": f"v_{ell}(deg(f)) is odd",
            "maximal_order_representative_principally_polarizable": False,
            "jacobian_in_isogeny_class": "unresolved",
        },
    }
    (here / "local_orders.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
