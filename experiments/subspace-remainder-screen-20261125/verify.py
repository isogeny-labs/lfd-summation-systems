"""Independent finite-field check of the subspace-remainder correspondence."""
from __future__ import annotations

import json
import sys
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402


def linearized_coefficients(F: GF2n, basis: list[int]) -> list[int]:
    # L_{V+<v>}(Z) = L_V(Z)^2 + L_V(v)L_V(Z), starting with L_{0}(Z)=Z.
    coeff = [1]
    for v in basis:
        h = linearized_eval(F, coeff, v)
        assert h != 0, 'basis must be independent'
        nxt = [F.mul(h, c) for c in coeff] + [0]
        for i, c in enumerate(coeff):
            nxt[i + 1] ^= F.sq(c)
        coeff = nxt
    return coeff


def linearized_eval(F: GF2n, coeff: list[int], x: int) -> int:
    power, out = x, 0
    for c in coeff:
        out ^= F.mul(c, power)
        power = F.sq(power)
    return out


def quadratic_remainder(F: GF2n, coeff: list[int], e1: int, e2: int):
    # Z^(2^j) = A_j Z + B_j modulo Z^2+e1 Z+e2.
    a, b = 1, 0
    ra, rb = 0, 0
    for c in coeff:
        ra ^= F.mul(c, a)
        rb ^= F.mul(c, b)
        aa = F.sq(a)
        a, b = F.mul(e1, aa), F.mul(e2, aa) ^ F.sq(b)
    return ra, rb


def check(n: int, modulus: int, k: int):
    F = GF2n(n, modulus)
    basis = [1 << j for j in range(k)]
    V = {j for j in range(1 << k)}
    coeff = linearized_coefficients(F, basis)
    assert len(coeff) == k + 1 and coeff[0] != 0
    checked = off_diagonal = diagonal = 0
    for e1 in range(1 << n):
        for e2 in range(1 << n):
            a, b = quadratic_remainder(F, coeff, e1, e2)
            roots = [x for x in range(1 << n)
                     if F.sq(x) ^ F.mul(e1, x) ^ e2 == 0]
            assert (a == 0 and b == 0) == (len(roots) == 2 and all(x in V for x in roots))
            # Symmetric image equations retain the repeated-root branch.
            image_sum = F.mul(e1, a)
            image_product = F.mul(e2, F.sq(a)) ^ F.mul(e1, F.mul(a, b)) ^ F.sq(b)
            assert (image_sum == 0 and image_product == 0) == (
                len(roots) > 0 and all(x in V for x in roots))
            if e1:
                assert F.mul(a, e1) == linearized_eval(F, coeff, e1)
                off_diagonal += int(a == 0 and b == 0)
            else:
                x = F.sqrt(e2)
                assert roots == [x]
                assert (x in V) == (linearized_eval(F, coeff, x) == 0)
                diagonal += int(x in V)
            checked += 1
    return {'n': n, 'k': k, 'modulus': modulus, 'checked_quadratics': checked,
            'off_diagonal_member_quadratics': off_diagonal,
            'diagonal_member_quadratics': diagonal, 'status': 'passed'}


def check_full_two_torsion_quotient():
    x, i = sp.symbols('x'), sp.I
    j = (x*x - 1)**2 / (x*(x*x + 1))
    duplicate_x = (3*x*x + 1)**2 / (4*x*(x*x + 1)) - 2*x
    transforms = (1/x, i*(x+i)/(x-i), -i*(x-i)/(x+i))
    assert sp.cancel(j - 4*duplicate_x) == 0
    assert all(sp.cancel(j.subs(x, t) - j) == 0 for t in transforms)
    return {'identity': 'J(x)=4*x([2]P)', 'torsion_translations_checked': 3,
            'status': 'passed'}


if __name__ == '__main__':
    rows = [check(5, (1 << 5) | (1 << 2) | 1, 2),
            check(7, (1 << 7) | (1 << 1) | 1, 3)]
    out = Path(__file__).with_name('verification.json')
    out.write_text(json.dumps({'rows': rows,
                               'full_two_torsion_quotient': check_full_two_torsion_quotient()},
                              indent=2) + '\n')
    for row in rows:
        print(row)
