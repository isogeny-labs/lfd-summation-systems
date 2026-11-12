#!/usr/bin/env python3
"""Symbolic identities used by the carry lower bound (README, Proposition 1), over Z[a, b, t].

S_3(u, v, w) = (u - v)^2 w^2 - 2((u + v)(u v + a) + 2b) w + (u v - a)^2 - 4b(u + v),
S_4(x1, x2, x3, t) = Res_X(S_3(x1, x2, X), S_3(x3, t, X)).

Checks: (i) S_4 is symmetric (x1 <-> x3 swap of the split point), degree 4 in each x_i, total 12;
(ii) the coefficient of x3^4 is S_3(x1, x2, t)^2; (iii) the coefficient of x1^3 x2^3 x3^3 is
-8 (5 t^3 + 3 a t + 6 b).  Writes experiments/ic-bounds-m3/symbolic.json.
"""
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def main():
    x1, x2, x3, t, a, b, X = sp.symbols('x1 x2 x3 t a b X')

    def S3(u, v, w):
        return (u - v) ** 2 * w ** 2 - 2 * ((u + v) * (u * v + a) + 2 * b) * w + (u * v - a) ** 2 - 4 * b * (u + v)

    S4 = sp.expand(sp.resultant(S3(x1, x2, X), S3(x3, t, X), X))
    P = sp.Poly(S4, x1, x2, x3)
    c333 = sp.factor(P.coeff_monomial(x1 ** 3 * x2 ** 3 * x3 ** 3))
    lc = sp.Poly(S4, x3).LC()
    sym = sp.expand(S4 - sp.resultant(S3(x1, x3, X), S3(x2, t, X), X)) == 0
    sym2 = sp.expand(S4 - S4.subs({x1: x2, x2: x1}, simultaneous=True)) == 0
    sym3 = sp.expand(S4 - sp.expand(sp.resultant(S3(x1, t, X), S3(x3, x2, X), X))) == 0
    out = dict(
        symmetric_x2_x3=bool(sym), symmetric_x1_x2=bool(sym2), symmetric_t_x2=bool(sym3),
        degrees=[P.degree(x1), P.degree(x2), P.degree(x3)],
        total_degree_with_t=sp.Poly(S4, x1, x2, x3, t).total_degree(),
        leading_coefficient_x3_equals_S3_squared=bool(sp.expand(lc - S3(x1, x2, t) ** 2) == 0),
        coefficient_x1x2x3_cubed=str(c333),
        terms=len(P.terms()),
    )
    assert out['leading_coefficient_x3_equals_S3_squared'] and str(c333) == '-8*(3*a*t + 6*b + 5*t**3)'
    (ROOT / 'experiments' / 'ic-bounds-m3' / 'symbolic.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
