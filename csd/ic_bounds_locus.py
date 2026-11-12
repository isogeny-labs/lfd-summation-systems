"""Algebraic features of a prime target (a, b, t) for the exceptional-locus map
(bounds program item 2; scripts/ic_bounds_locus*.py).

Every feature is a polynomial condition in (a, b, t) over F_p, or a membership
condition in the 2-power subgroup H = {x : x^(2^k) = 1} of F_p^*.  The PKM
generator NF(S_3(x1, x2, t)) has coefficients

    c22 = 1, c20 = c02 = t^2, c21 = c12 = -2t, c11 = -2(t^2 + a),
    c10 = c01 = -2(a t + 2b), c00 = a^2 - 4 b t,

so the vanishing of each coefficient is one of the conditions below.
"""
from __future__ import annotations


def in_H(x: int, p: int, k: int) -> bool:
    x %= p
    return x != 0 and pow(x, 1 << k, p) == 1


def x_double(p: int, a: int, b: int, t: int):
    """x(2R) for R = (t, y); None when y = 0 (2R = O)."""
    f = (t ** 3 + a * t + b) % p
    if f == 0:
        return None
    return (t ** 4 - 2 * a * t * t - 8 * b * t + a * a) * pow(4 * f, -1, p) % p


def two_torsion_roots(p: int, a: int, b: int):
    return [e for e in range(p) if (e ** 3 + a * e + b) % p == 0]


def polynomial_values(p: int, a: int, b: int, t: int) -> dict[str, int]:
    """Named polynomials in (a, b, t); the feature '<name>=0' is their vanishing."""
    return {
        't': t % p,
        'rhs': (t ** 3 + a * t + b) % p,                    # y = 0: target of order two
        '3t2+a': (3 * t * t + a) % p,                        # rhs'(t)
        't2+a': (t * t + a) % p,                             # c11 = 0
        'at+2b': (a * t + 2 * b) % p,                        # c10 = 0
        'a2-4bt': (a * a - 4 * b * t) % p,                   # c00 = 0
        'psi3': (3 * t ** 4 + 6 * a * t * t + 12 * b * t - a * a) % p,   # R of order 3
        'x2R_num': (t ** 4 - 2 * a * t * t - 8 * b * t + a * a) % p,     # x(2R) = 0
        't2-a': (t * t - a) % p,
        't2+3a': (t * t + 3 * a) % p,
        't3-b': (t ** 3 - b) % p,
        't3+b': (t ** 3 + b) % p,
        '2t3-b': (2 * t ** 3 - b) % p,
        'at-b': (a * t - b) % p,
        'at+b': (a * t + b) % p,
        'at+3b': (a * t + 3 * b) % p,
        '3at+4b': (3 * a * t + 4 * b) % p,
        'a': a % p,
        'b': b % p,
    }


def features(p: int, a: int, b: int, t: int, k: int) -> dict[str, bool]:
    vals = polynomial_values(p, a, b, t)
    out = {f'{name}=0': v == 0 for name, v in vals.items()}
    # membership of the nonzero values in H (class conditions; H contains -1)
    for name in ('t', 'rhs', '3t2+a', 't2+a', 'at+2b', 'a2-4bt', 'a', 'b'):
        v = vals[name]
        out[f'{name} in H'] = in_H(v, p, k)
    x2 = x_double(p, a, b, t)
    out['x2R in H'] = x2 is not None and in_H(x2, p, k)
    out['x2R=t'] = x2 is not None and (x2 - t) % p == 0
    # ratios of coefficients of the generator (orbit invariants of the H-torus scaling)
    c = generator_coefficients(p, a, b, t)
    for name, (num, den) in {
        'c20/c22': (c['20'], c['22']),
        'c00/c22': (c['00'], c['22']),
        'c10/c00': (c['10'], c['00']),
        'c11/c22': (c['11'], c['22']),
        'c21/c20': (c['21'], c['20']),
        'c10/c11': (c['10'], c['11']),
        'c00/c11': (c['00'], c['11']),
        'c20/c00': (c['20'], c['00']),
    }.items():
        out[f'{name} in H'] = bool(num % p and den % p and in_H(num * pow(den, -1, p), p, k))
    roots = [e for e in two_torsion_roots(p, a, b) if (t - e) % p]
    out['some t-e in H'] = any(in_H(t - e, p, k) for e in roots)       # e: other 2-torsion abscissae
    out['all t-e in H'] = bool(roots) and all(in_H(t - e, p, k) for e in roots)
    out['three rational 2-torsion'] = len(two_torsion_roots(p, a, b)) == 3
    out['j=0'] = a % p == 0
    out['j=1728'] = b % p == 0
    return out


def generator_coefficients(p: int, a: int, b: int, t: int) -> dict[str, int]:
    return {'22': 1, '20': t * t % p, '21': -2 * t % p, '11': -2 * (t * t + a) % p,
            '10': -2 * (a * t + 2 * b) % p, '00': (a * a - 4 * b * t) % p}
