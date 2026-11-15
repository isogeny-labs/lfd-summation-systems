#!/usr/bin/env python3
"""Exact count of the primes over which the fixed E6 system is a proper ideal.

The system of y^2 = x^3 + 6x - 7 with target abscissa 1 at depth k has a common zero
over F_p (p = 1 mod 2^k) exactly when some pair of 2^k-th roots of unity (x1, x2)
satisfies (x1 - 1) x2 - x1 - 8 = 0, that is x2 = T(x1) = (x1 + 8)/(x1 - 1): the chain
variables are then the successive squares. Otherwise the ideal is the unit ideal
(the grid is the whole zero set of the chain equations, and it is reduced).

For every depth k = 2, ..., 10 this script enumerates all primes p = 1 mod 2^(k+1)
below 64 * 4^k and decides the question for each by running over the 2^k roots.
It reports the exceptional primes by size class relative to 4^k, the count expected
from the heuristic "each x1 gives a root of unity with probability 2^k/(p-1)", and
checks that every prime factor of the Fermat number 2^(2^k)+1 in range is exceptional
(formal/ICFallDegrees/GrowthFermatDiagonalObstruction.lean proves this for all k).

    python3 scripts/check_n8_prime_grid_exceptions.py            # recompute, compare with the record
    python3 scripts/check_n8_prime_grid_exceptions.py --write    # recompute and write the record

Record: experiments/ic-bounds-growth/prime-grid-exceptions.json. Integers only.
"""
import json
import math
import sys
from pathlib import Path

RECORD = Path(__file__).resolve().parents[1] / 'experiments/ic-bounds-growth/prime-grid-exceptions.json'
DEPTHS = range(2, 11)
CLASSES = [(0, 1), (1, 4), (4, 16), (16, 64)]          # p / 4^k in [lo, hi)
MR_BASES = (2, 3, 5, 7, 11, 13, 17)                    # deterministic below 3.4e14


def is_prime(n):
    if n < 2:
        return False
    for q in MR_BASES:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def roots_of_unity(p, k):
    """The 2^k-th roots of unity of F_p, p = 1 mod 2^k."""
    n = 1 << k
    r = 2
    while pow(r, (p - 1) // 2, p) != p - 1:            # least quadratic non-residue
        r += 1
    g = pow(r, (p - 1) // n, p)                        # order exactly 2^k
    out, x = [], 1
    for _ in range(n):
        out.append(x)
        x = x * g % p
    assert x == 1 and len(set(out)) == n
    return out


def common_zeros(p, k):
    """Pairs (x1, x2) of 2^k-th roots of unity with (x1 - 1) x2 - x1 - 8 = 0 mod p."""
    mu = roots_of_unity(p, k)
    members = set(mu)
    found = []
    for x in mu:
        if x == 1:                                     # then the form is -9, nonzero for p != 3
            continue
        y = (x + 8) * pow(x - 1, -1, p) % p
        if y in members:
            assert ((x - 1) * y - x - 8) % p == 0
            found.append((x, y))
    return found


def fermat_factors_below(k, bound):
    f = (1 << (1 << k)) + 1
    return [p for p in range(1, bound, 1 << (k + 1)) if p > 1 and f % p == 0 and is_prime(p)]


def compute():
    rows = []
    for k in DEPTHS:
        four_k = 4 ** k
        bound = 64 * four_k
        primes = [p for p in range(1 + (1 << (k + 1)), bound, 1 << (k + 1)) if is_prime(p)]
        exceptional = {p for p in primes if common_zeros(p, k)}
        fermat = fermat_factors_below(k, bound)
        assert all(p in exceptional for p in fermat), (k, fermat)
        classes = []
        for lo, hi in CLASSES:
            inside = [p for p in primes if lo * four_k <= p < hi * four_k]
            bad = [p for p in inside if p in exceptional]
            # unordered solutions {x1, x2}: about 4^k / (2p) expected; Poisson heuristic for "none"
            expected = sum(1 - math.exp(-four_k / (2 * p)) for p in inside)
            classes.append(dict(range_over_4k=[lo, hi], primes=len(inside), exceptional=len(bad),
                                heuristic_expected=round(expected, 1)))
        rows.append(dict(k=k, modulus=1 << (k + 1), bound=bound, primes=len(primes),
                         exceptional=len(exceptional),
                         smallest_exceptional=sorted(exceptional)[:5],
                         largest_exceptional=max(exceptional) if exceptional else None,
                         smallest_good=min(p for p in primes if p not in exceptional),
                         fermat_factors_in_range=fermat, classes=classes))
    totals = []
    for i, (lo, hi) in enumerate(CLASSES):
        totals.append(dict(range_over_4k=[lo, hi],
                           primes=sum(r['classes'][i]['primes'] for r in rows),
                           exceptional=sum(r['classes'][i]['exceptional'] for r in rows),
                           heuristic_expected=round(sum(r['classes'][i]['heuristic_expected'] for r in rows), 1)))
    return dict(
        statement='primes p = 1 mod 2^(k+1), p < 64*4^k, over which the depth-k system of '
                  'y^2 = x^3 + 6x - 7 with target abscissa 1 has a common zero (proper ideal)',
        method='exhaustive over the 2^k-th roots of unity x1: is (x1+8)/(x1-1) a 2^k-th root of unity',
        depths=rows, totals_by_class=totals)


def main():
    record = compute()
    text = json.dumps(record, indent=1) + '\n'
    if '--write' in sys.argv[1:]:
        RECORD.write_text(text)
        print(f'wrote {RECORD.relative_to(RECORD.parents[2])}')
    else:
        saved = json.loads(RECORD.read_text())
        assert saved == json.loads(text), 'recomputed counts differ from the record'
    for row in record['depths']:
        print(f"k={row['k']:2d}: {row['exceptional']:4d} exceptional of {row['primes']:4d} primes below 64*4^k; "
              f"smallest good prime {row['smallest_good']}; Fermat factors {row['fermat_factors_in_range']}")
    for t in record['totals_by_class']:
        print(f"p/4^k in [{t['range_over_4k'][0]}, {t['range_over_4k'][1]}): {t['exceptional']} of {t['primes']} "
              f"exceptional (heuristic {t['heuristic_expected']})")
    print('PASS exact prime-grid exceptions, k = 2..10')


if __name__ == '__main__':
    main()
