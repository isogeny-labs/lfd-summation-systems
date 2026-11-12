"""Exact last-fall-degree sweep over prime, binary ordinary and Koblitz curves
with random and symmetric (automorphism-invariant) factor bases.

Factor-base classes:

* prime curves, subgroup chain presentation of depth k (2^k | p-1):
  ``subgroup`` (symmetric: the subgroup H of order 2^k, closed under
  negation and under multiplication by H) and ``random-coset`` (x = c*u with
  u in H and a random c outside H; the coset is still a factor base of
  size 2^k but is not a group);
* binary ordinary and Koblitz curves over F_(2^n), F_2-linear coordinate
  spaces of dimension k: ``random-subspace`` (a uniformly random basis) and
  ``frobenius-invariant`` (symmetric: a subspace stable under the squaring
  Frobenius, which includes every subfield; on a Koblitz curve the Frobenius
  is a curve automorphism, so this base is invariant under it).

Every case stores its system hash, the scalar-S3 solution count, the number
of rational decompositions checked with the group law, and an exact ordinary
iterated LFD from a certified closure profile.  Prime chain profiles use the
same local-to-global rank criterion as the Boolean profiles: evaluation of the
squarefree monomials of degree below D onto the solution set, and the closure
quotient of codimension equal to the number of solutions.
"""
from __future__ import annotations

from collections import Counter
from itertools import product
import hashlib
import json
import math
import random

import numpy as np

from . import elim
from . import ic_dreg as I
from . import ic_fall_certificates as C
from .ic_arbitrary_bases import canonical_basis, scaled_subfield, boolean_points
from .ic_local_certificate import local_profile
from .ic_research_campaign import field, geometry


# ---------------------------------------------------------------------------
# prime curves: subgroup and coset chain presentations
# ---------------------------------------------------------------------------


def two_adic_valuation(n: int) -> int:
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def prime_coset_generator(E: I.PrimeCurve, xR: int, k: int, c1: int, c2: int) -> dict[int, int]:
    """NF(S_3(c1*u1, c2*u2, xR)) in F_p[u1,u2]/(u_i^{2^k}-1) as {chain mask: coefficient}.
    With c1 = c2 = 1 this is the PKM subgroup generator."""
    N = 1 << k
    g: dict[int, int] = {}
    for (e1, e2), c in I.s3_prime_coeffs(E, xR).items():
        scaled = c * pow(c1, e1, E.p) * pow(c2, e2, E.p) % E.p
        mask = (e1 % N) | ((e2 % N) << k)
        g[mask] = (g.get(mask, 0) + scaled) % E.p
    return {mask: c for mask, c in g.items() if c}


def prime_points(E: I.PrimeCurve, k: int, generator_value: int, c1: int, c2: int, xR: int) -> list[int]:
    """All chain assignments (E1 | E2 << k) with S_3(c1 g^E1, c2 g^E2, xR) = 0."""
    N = 1 << k
    powers = [pow(generator_value, e, E.p) for e in range(N)]
    out = []
    for e1 in range(N):
        x1 = c1 * powers[e1] % E.p
        for e2 in range(N):
            x2 = c2 * powers[e2] % E.p
            if I.s3_prime_value(E, x1, x2, xR) == 0:
                out.append(e1 | (e2 << k))
    return out


def prime_evaluation_rank(points: list[int], masks: list[int], k: int, generator_value: int, p: int,
                          max_cells: int) -> int:
    """Rank over F_p of the evaluation matrix of squarefree chain monomials at the points.
    A monomial with block exponents (a1, a2) evaluates to g^(a1 E1 + a2 E2) at (E1, E2)."""
    if len(points) * max(1, len(masks)) > max_cells:
        raise ValueError('Evaluation matrix budget exceeded.')
    N = 1 << k
    order = N
    matrix = np.zeros((len(points), len(masks)), dtype=np.int64)
    for i, bits in enumerate(points):
        E1, E2 = bits & (N - 1), bits >> k
        for j, mask in enumerate(masks):
            a1, a2 = mask & (N - 1), mask >> k
            matrix[i, j] = pow(generator_value, (a1 * E1 + a2 * E2) % order, p)
    return elim.fp_echelon(matrix, p)[0]


def prime_local_profile(generator: dict[int, int], k: int, p: int, points: list[int], generator_value: int, *,
                        max_degree: int = 6, max_columns: int = 16384, max_cells: int = 120_000_000):
    """Exact ordinary LFD of {g} + chain for two summands, certified by ranks.

    unit ideal: the first degree deriving 1.  proper ideal: once evaluation of
    the squarefree monomials of degree < D is onto the solution set and the
    closure quotient at D has codimension #points, every polynomial of degree
    <= D reduces at budget D to a representative of degree < D, so the full
    ideal has LFD <= D (local-to-global theorem) and the exact value is the
    last strict fall.  Each closure is checked against the containing-space
    certificate before it is used.
    """
    points = list(points)
    if len(points) != len(set(points)):
        raise ValueError('Points must be distinct.')
    N = 1 << k
    for bits in points:
        E1, E2 = bits & (N - 1), bits >> k
        value = sum(c * pow(generator_value, ((mask & (N - 1)) * E1 + (mask >> k) * E2) % N, p)
                    for mask, c in generator.items()) % p
        if value:
            raise ValueError('A supplied point does not vanish on the generator.')
    spaces = []
    falls = []
    checks = []
    reason = 'configured degree cap'
    for D in range(max_degree + 1):
        try:
            current = C.chain_degree_closure([generator], 2, k, p, D, max_columns=max_columns, max_cells=max_cells)
        except ValueError as error:
            reason = str(error)
            break
        if spaces and current.low_degree_dimension(D - 1) > spaces[-1].rank:
            falls.append(D)
        spaces.append(current)
        if current.contains_one:
            assert not points
            return dict(status='exact-unit', exact_lfd=D, lfd_lower_bound=D, lfd_upper_bound=D,
                        strict_falls=falls, scanned_through=D,
                        quotient_ranks=[s.quotient_rank for s in spaces], point_count=0,
                        local_checks=checks, containing_spaces_checked=True)
        assert C.check_closed_exclusion_certificate(current, [generator])
        if D > 0 and points:
            try:
                rank = prime_evaluation_rank(points, [m for m in current.masks if m.bit_count() < D],
                                             k, generator_value, p, max_cells)
            except ValueError as error:
                reason = str(error)
                break
            codim = len(current.masks) - current.quotient_rank
            checks.append(dict(degree=D, evaluation_rank=rank, codimension=codim))
            if rank == len(points) == codim:
                last = max(falls, default=0)
                return dict(status='exact-proper-local', exact_lfd=last, lfd_lower_bound=last,
                            lfd_upper_bound=last, strict_falls=falls, scanned_through=D,
                            local_upper_bound=D, point_count=len(points),
                            quotient_ranks=[s.quotient_rank for s in spaces], local_checks=checks,
                            containing_spaces_checked=True)
    return dict(status='capped', exact_lfd=None, lfd_lower_bound=max(falls, default=0),
                lfd_upper_bound=2 * k + 1, strict_falls=falls, scanned_through=len(spaces) - 1,
                quotient_ranks=[s.quotient_rank for s in spaces], point_count=len(points),
                local_checks=checks, cap_reason=reason, containing_spaces_checked=True)


def subgroup_generator_value(p: int, k: int) -> int:
    """The generator g with H = {g^E}, matching ic_dreg.subgroup_2power."""
    values = I.subgroup_2power(p, k)
    return values[1] if k >= 1 else 1


def prepare_prime(E: I.PrimeCurve, k: int, base_kind: str, scale: int, target, target_class, planted=None):
    p = E.p
    N = 1 << k
    g0 = subgroup_generator_value(p, k)
    H = [pow(g0, e, p) for e in range(N)]
    assert len(set(H)) == N
    if base_kind == 'subgroup':
        assert scale == 1
    else:
        assert base_kind == 'random-coset' and scale % p and scale not in H
    values = [scale * h % p for h in H]
    generator = prime_coset_generator(E, target[0], k, scale, scale)
    points = prime_points(E, k, g0, scale, scale, target[0])
    system = dict(kind='prime', prime=p, m=2, k=k, variables=2 * k, generators=[sorted(generator.items())],
                  generator_value=g0, coordinate_solutions=points)
    key = hashlib.sha256(json.dumps(dict(kind='prime', prime=p, m=2, k=k, generators=system['generators']),
                                    sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    geo = geometry(E, values, 2, target, 'direct')
    assert geo['polynomial_endpoint_tuples'] == len(points)
    if planted:
        P, Q = map(tuple, planted)
        assert E.add(P, Q) == target and P[0] in values and Q[0] in values
        assert geo['rational_endpoint_tuples'] > 0
    case = dict(kind='prime', curve_class=('j0' if E.a == 0 else 'j1728' if E.b == 0 else 'random'),
                prime=p, a=E.a, b=E.b, k=k, m=2, base_kind=base_kind, symmetric=base_kind == 'subgroup',
                base_scale=scale, base_size=N, target=list(target), target_class=target_class,
                planted=planted, polynomial_tuples=len(points),
                rational_ordered_tuples=geo['rational_endpoint_tuples'],
                rational_unordered_multisets=geo['rational_unordered_coordinate_multisets'],
                independent_relation_rank=None, system_hash=key)
    case['case_id'] = hashlib.sha256(json.dumps(case, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:20]
    return case, dict(system_hash=key, **system)


def evaluate_prime(system, max_degree=6, max_columns=16384, max_cells=200_000_000):
    generator = dict((int(m), int(c)) for m, c in system['generators'][0])
    result = prime_local_profile(generator, system['k'], system['prime'], system['coordinate_solutions'],
                                 system['generator_value'], max_degree=max_degree, max_columns=max_columns,
                                 max_cells=max_cells)
    if result['exact_lfd'] is not None:
        assert result['exact_lfd'] <= 2 * system['k'] + 1
    return dict(system_hash=system['system_hash'], **result)


# ---------------------------------------------------------------------------
# binary ordinary and Koblitz curves: random and Frobenius-invariant subspaces
# ---------------------------------------------------------------------------


def prepare_binary(E: I.BinaryCurve, basis: list[int], base_kind: str, target, target_class, planted=None):
    F = E.F
    k = len(basis)
    basis = canonical_basis(basis)
    assert len(basis) == k and E.on_curve(target)
    invariant = I.is_frobenius_invariant(F, basis)
    if base_kind == 'frobenius-invariant':
        assert invariant
    equations = sorted(set(tuple(sorted(row)) for row in I.binary_descent(E, [basis, basis], target[0]) if row))
    system = dict(kind='binary', variables=2 * k, equations=equations)
    key = hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    points = boolean_points(equations, k)
    values = I.span(basis)
    rational = 0
    unordered = set()
    for bits in points:
        x, y = values[bits & ((1 << k) - 1)], values[bits >> k]
        assert I.s3_binary_value(F, E.b, x, y, target[0]) == 0
        P = E.lift_x(x)
        if P is None:
            continue
        if any((Q := E.add(target, E.neg(P0))) is not None and Q[0] == y for P0 in {P, E.neg(P)}):
            rational += 1
            unordered.add(tuple(sorted((x, y))))
    if planted:
        P, Q = map(tuple, planted)
        assert E.add(P, Q) == target and P[0] in values and Q[0] in values
        assert rational > 0
    koblitz = E.a in (0, 1) and E.b == 1
    case = dict(kind='koblitz' if koblitz else 'binary', n=F.N, modulus=F.modulus, k=k, m=2, a=E.a, b=E.b,
                basis=basis, base_kind=base_kind, symmetric=base_kind == 'frobenius-invariant',
                frobenius_invariant=invariant, scaled_subfield=scaled_subfield(F, basis),
                subfield=bool(scaled_subfield(F, basis) and 1 in values and all(F.frob(v, k) == v for v in basis)),
                base_size=1 << k, target=list(target), target_class=target_class, planted=planted,
                polynomial_tuples=len(points), rational_ordered_tuples=rational,
                rational_unordered_multisets=len(unordered), independent_relation_rank=None, system_hash=key)
    case['case_id'] = hashlib.sha256(json.dumps(case, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:20]
    return case, dict(system_hash=key, **system, coordinate_solutions=points)


def evaluate_binary(system, max_degree=6, max_columns=16384, max_bytes=384_000_000):
    k = system['variables'] // 2
    result = local_profile(system['equations'], system['variables'], system['coordinate_solutions'],
                           max_degree=max_degree, max_columns=max_columns, max_bytes=max_bytes)
    if result['exact_lfd'] is not None:
        assert result['exact_lfd'] <= k + 2
    return dict(system_hash=system['system_hash'], **result)


def random_nonstable_subspace(F, k: int, rng: random.Random, forbidden: list[list[int]]) -> list[int]:
    for _ in range(10000):
        B = canonical_basis(I.random_subspace(F, k, rng))
        if not I.is_frobenius_invariant(F, B) and B not in forbidden:
            return B
    raise RuntimeError('No non-invariant random subspace found.')


def _f2_mul(a: int, b: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def frobenius_invariant_subspaces(F, d: int) -> list[list[int]]:
    """Bases of every Frobenius-invariant subspace of dimension d, for every n.

    F_(2^n) is a cyclic F_2[sigma]-module isomorphic to F_2[x]/(x^n-1), so the
    invariant subspaces are the images of g(sigma) for the divisors g of x^n-1,
    of dimension n-deg g.  Write n = 2^e * n' with n' odd; x^n-1 = (x^n'-1)^(2^e)
    and the divisors are products of the irreducible factors of x^n'-1 with
    multiplicities up to 2^e.  For odd n this agrees with ic_dreg.invariant_subspaces.
    """
    n = F.N
    e = two_adic_valuation(n)
    odd = n >> e
    factors = [(1 << 1) | 1] if odd == 1 else I._f2_poly_factor_xn1(odd)
    degrees = [f.bit_length() - 1 for f in factors]
    out = []
    for exponents in product(range((1 << e) + 1), repeat=len(factors)):
        if sum(a * deg for a, deg in zip(exponents, degrees)) != n - d:
            continue
        g = 1
        for f, a in zip(factors, exponents):
            for _ in range(a):
                g = _f2_mul(g, f)
        basis = I._f2_basis([I._apply_sigma_poly(F, g, 1 << j) for j in range(n)])
        assert len(basis) == d, (n, d, exponents, len(basis))
        assert I.is_frobenius_invariant(F, basis)
        out.append(basis)
    return out


def invariant_candidates(F, k: int, limit: int, rng: random.Random) -> list[list[int]]:
    spaces = [canonical_basis(B) for B in frobenius_invariant_subspaces(F, k)]
    unique = []
    for B in spaces:
        if B not in unique:
            unique.append(B)
    if len(unique) <= limit:
        return unique
    # Prefer the subfield when it exists, then a seeded sample.
    subfield = [B for B in unique if 1 in I.span(B) and all(F.frob(v, k) == v for v in B)]
    rest = [B for B in unique if B not in subfield]
    chosen = subfield[:1] + rng.sample(rest, limit - len(subfield[:1]))
    return chosen


def random_rational_target(E, rng: random.Random, modulus: int, exclude_x=()):
    for _ in range(5000):
        x = rng.randrange(1, modulus)
        if x in exclude_x:
            continue
        P = E.lift_x(x)
        if P is not None:
            return P
    raise RuntimeError('No rational target found.')


def planted_target(E, values, rng: random.Random):
    points = [P for x in values if x and (P := E.lift_x(x)) is not None]
    for _ in range(400):
        if len(points) < 2:
            return None
        P, Q = rng.sample(points, 2)
        if rng.randrange(2):
            P = E.neg(P)
        if rng.randrange(2):
            Q = E.neg(Q)
        R = E.add(P, Q)
        if R is not None and R[0]:
            return R, [list(P), list(Q)]
    return None


def order_two_prime_target(E: I.PrimeCurve):
    for t in range(E.p):
        if E.rhs(t) % E.p == 0:
            return (t, 0)
    return None


def generate(seed: int = 20261104, prime_list=(97, 193, 257, 641, 7681), binary_degrees=range(5, 13),
             max_k_prime: int = 7, max_k_binary: int = 7, random_bases: int = 2, invariant_limit: int = 2):
    rng = random.Random(seed)
    cases = []
    systems = {}
    skips = []
    # prime curves
    for p in prime_list:
        kmax = min(max_k_prime, two_adic_valuation(p - 1))
        curves = []
        for a, b in [(0, 1), (1, 0)]:
            if (4 * a ** 3 + 27 * b * b) % p:
                curves.append((a, b))
        while len(curves) < 4:
            a, b = rng.randrange(1, p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p and (a, b) not in curves:
                curves.append((a, b))
        for a, b in curves:
            E = I.PrimeCurve(p, a, b)
            for k in range(2, kmax + 1):
                N = 1 << k
                g0 = subgroup_generator_value(p, k)
                H = [pow(g0, e, p) for e in range(N)]
                bases = [('subgroup', 1)]
                if N == p - 1:
                    skips.append(dict(kind='prime', prime=p, a=a, b=b, k=k, base_kind='random-coset',
                                      reason='The subgroup is the whole multiplicative group; no proper coset exists'))
                while N < p - 1 and len(bases) < 1 + random_bases:
                    c = rng.randrange(2, p)
                    if c not in H and ('random-coset', c) not in bases:
                        bases.append(('random-coset', c))
                for base_kind, scale in bases:
                    values = [scale * h % p for h in H]
                    targets = []
                    two = order_two_prime_target(E)
                    if two is not None:
                        targets.append((two, 'order-two', None))
                    for _ in range(2):
                        R = random_rational_target(E, rng, p, exclude_x={t[0][0] for t in targets})
                        targets.append((R, 'random-rational', None))
                    plant = planted_target(E, values, rng)
                    if plant is None:
                        skips.append(dict(kind='prime', prime=p, a=a, b=b, k=k, base_kind=base_kind,
                                          reason='No nonzero planted sum found'))
                    else:
                        targets.append((plant[0], 'planted-rational-decomposition', plant[1]))
                    for target, kind, planted in targets:
                        case, system = prepare_prime(E, k, base_kind, scale, target, kind, planted)
                        cases.append(case)
                        systems.setdefault(system['system_hash'], system)
    # binary ordinary and Koblitz
    for n in binary_degrees:
        F = field(n)
        curves = [(0, 1), (1, 1)]
        while len(curves) < 4:
            a, b = rng.randrange(1 << n), rng.randrange(2, 1 << n)
            if (a, b) not in curves:
                curves.append((a, b))
        for k in range(2, min(max_k_binary, n - 1) + 1):
            bases = []
            for B in invariant_candidates(F, k, invariant_limit, rng):
                bases.append(('frobenius-invariant', B))
            if not any(kind == 'frobenius-invariant' for kind, _ in bases):
                skips.append(dict(kind='binary', n=n, k=k, reason='No Frobenius-invariant subspace of this dimension'))
            forbidden = [B for _, B in bases]
            for _ in range(random_bases):
                B = random_nonstable_subspace(F, k, rng, forbidden)
                forbidden.append(B)
                bases.append(('random-subspace', B))
            for a, b in curves:
                E = I.BinaryCurve(F, a, b)
                for base_kind, B in bases:
                    values = I.span(B)
                    targets = [((0, F.sqrt(b)), 'order-two', None)]
                    for _ in range(2):
                        R = random_rational_target(E, rng, 1 << n, exclude_x={t[0][0] for t in targets})
                        targets.append((R, 'random-rational', None))
                    plant = planted_target(E, values, rng)
                    if plant is None:
                        skips.append(dict(kind='binary', n=n, a=a, b=b, k=k, base_kind=base_kind,
                                          reason='No nonzero planted sum found'))
                    else:
                        targets.append((plant[0], 'planted-rational-decomposition', plant[1]))
                    for target, kind, planted in targets:
                        case, system = prepare_binary(E, B, base_kind, target, kind, planted)
                        cases.append(case)
                        systems.setdefault(system['system_hash'], system)
    return cases, systems, skips


def evaluate(system):
    if system['kind'] == 'prime':
        return evaluate_prime(system)
    return evaluate_binary(system)


def summarize(cases, results):
    def histogram(subset):
        return dict(Counter(str(results[h]['exact_lfd']) for h in sorted({c['system_hash'] for c in subset})))

    groups = []
    for kind in ('prime', 'binary', 'koblitz'):
        for base_kind in sorted({c['base_kind'] for c in cases if c['kind'] == kind}):
            for k in sorted({c['k'] for c in cases if c['kind'] == kind}):
                subset = [c for c in cases if c['kind'] == kind and c['base_kind'] == base_kind and c['k'] == k]
                if not subset:
                    continue
                hashes = sorted({c['system_hash'] for c in subset})
                exact = [results[h]['exact_lfd'] for h in hashes if results[h]['exact_lfd'] is not None]
                groups.append(dict(kind=kind, base_kind=base_kind, symmetric=subset[0]['symmetric'], k=k,
                                   cases=len(subset), systems=len(hashes), exact_systems=len(exact),
                                   capped_systems=len(hashes) - len(exact),
                                   exact_lfd_histogram=histogram(subset),
                                   max_exact_lfd=max(exact, default=None),
                                   mean_exact_lfd=(sum(exact) / len(exact)) if exact else None,
                                   lean_upper_bound=(2 * k + 1),
                                   written_upper_bound=(2 * k + 1 if kind == 'prime' else k + 2)))
    return groups
