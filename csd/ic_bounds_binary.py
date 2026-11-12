"""Item 6 of the bounds programme: larger binary ordinary and Koblitz systems.

Exact ordinary iterated last fall degree (LFD) of Weil-descended summation
polynomial systems over F_(2^n), Boolean presentation (closure modulo the field
equations x_i^2 = x_i, as in ``csd.ic_boolean_closure``), with three families:

* ``m2``: two summands, S_3(x1, x2, xR), x1, x2 in one F_2-linear coordinate
  space of dimension k (2k Boolean variables), built by
  ``csd.ic_curve_sweep.prepare_binary`` (which cross-checks every coordinate
  zero against the scalar S_3 and counts rational decompositions with the
  group law);
* ``m3-direct``: three summands, Weil descent of S_4(x1, x2, x3, xR) built by
  ``csd.ic_dreg.binary_descent`` with three copies of the basis (3k variables);
* ``m3-split``: three summands, S_3(x1, x2, z) and S_3(z, x3, xR) with an
  auxiliary full-field coordinate z (n Boolean variables), built by
  ``csd.ic_research_campaign.binary_split`` (3k + n variables).

Every LFD is certified by ``csd.ic_local_certificate.local_profile``: either
1 in V_D (unit ideal, |Z| = 0) or, at the certificate degree D*, evaluation of
the squarefree monomials of degree < D* onto the complete coordinate zero set
Z and codim V_D* = |Z|; the exact value is then the last strict fall
(local-to-global criterion, formal/ICFallDegrees/LocalReductionBound). Each
closure is re-checked by ``check_containing_space`` (contains the inputs and is
stable under every variable shift). A budget overflow is recorded as
``capped`` with the last strict fall as a lower bound only.

Z is enumerated independently of the descended equations (field arithmetic):
``csd.ic_dreg.binary_solutions`` for m = 2 and m = 3 direct, and the scalar
S_3 chain for m = 3 split; ``local_profile`` re-checks that every point
vanishes on every Boolean equation.

Cost-related quantities recorded for every system (no ECDLP speed claim):
|Z|, input degree d_in, LFD, certificate degree D*, the squarefree
Macaulay column count |M_<=D| and the standard Macaulay row count
sum_e |M_<=D-deg e| at D = D* and at D = max(LFD, d_in), and dim V_D*.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import resource
import time

from . import ic_dreg as I
from . import ic_curve_sweep as S
from .ic_arbitrary_bases import canonical_basis
from .ic_boolean_closure import ClosureBudgetExceeded
from .ic_local_certificate import local_profile
from .ic_research_campaign import field, binary_split

SEED = 20261108


def _seed(*parts) -> int:
    return int.from_bytes(hashlib.sha256(':'.join(map(str, parts)).encode()).digest()[:8], 'big')


def curves_for(n: int, seed: int = SEED):
    """(0,1), (1,1) (Koblitz) and two random ordinary curves (b not in {0,1})."""
    rng = random.Random(_seed(seed, 'curves', n))
    curves = [(0, 1, 'koblitz'), (1, 1, 'koblitz')]
    while len(curves) < 4:
        a, b = rng.randrange(1 << n), rng.randrange(2, 1 << n)
        if all((a, b) != c[:2] for c in curves):
            curves.append((a, b, 'ordinary-random'))
    return curves


def bases_for(F, k: int, n_random: int, n_invariant: int, seed: int = SEED):
    rng = random.Random(_seed(seed, 'bases', F.N, k))
    out = []
    for B in (S.invariant_candidates(F, k, n_invariant, rng) if n_invariant > 0 else []):
        out.append(('frobenius-invariant', canonical_basis(B)))
    forbidden = [B for _, B in out]
    for _ in range(n_random):
        B = S.random_nonstable_subspace(F, k, rng, forbidden)
        forbidden.append(B)
        out.append(('random-subspace', B))
    return out


def targets_for(E, values, seed_parts, n_random: int = 2):
    rng = random.Random(_seed(*seed_parts))
    F = E.F
    targets = [((0, F.sqrt(E.b)), 'order-two', None)]
    for _ in range(n_random):
        R = S.random_rational_target(E, rng, 1 << F.N, exclude_x={t[0][0] for t in targets})
        targets.append((R, 'random-rational', None))
    plant = S.planted_target(E, values, rng)
    if plant is not None:
        targets.append((plant[0], 'planted-rational-decomposition', plant[1]))
    return targets


def squarefree_count(r: int, D: int) -> int:
    return sum(math.comb(r, i) for i in range(0, min(r, D) + 1)) if D >= 0 else 0


def macaulay_size(equations, r: int, D: int):
    """Standard Boolean Macaulay matrix at degree D: rows = every equation times
    every squarefree monomial of degree <= D - deg(equation); columns = M_<=D."""
    rows = sum(squarefree_count(r, D - max(m.bit_count() for m in eq)) for eq in equations if eq)
    return dict(degree=D, rows=rows, columns=squarefree_count(r, D))


def system_hash(kind: str, variables: int, equations) -> str:
    system = dict(kind=kind, variables=variables, equations=sorted(tuple(sorted(e)) for e in equations if e))
    return hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


# ---------------------------------------------------------------------------
# systems
# ---------------------------------------------------------------------------


def m2_system(E, basis, base_kind, target, target_class, planted):
    case, system = S.prepare_binary(E, basis, base_kind, target, target_class, planted)
    k = len(basis)
    # independent field-level enumeration of Z (no descended equations)
    sols = I.binary_solutions(E, [case['basis'], case['basis']], target[0])
    indep = sorted(c1 | (c2 << k) for c1, c2 in sols)
    assert indep == sorted(system['coordinate_solutions']), 'Z mismatch between descent and field enumeration'
    return dict(family='m2', variables=2 * k, equations=system['equations'], points=system['coordinate_solutions'],
                system_hash=system['system_hash'], rational_ordered_tuples=case['rational_ordered_tuples'],
                rational_unordered_multisets=case['rational_unordered_multisets'],
                scaled_subfield=case['scaled_subfield'], subfield=case['subfield'])


def m3_direct_system(E, basis, target):
    k = len(basis)
    bases = [basis, basis, basis]
    equations = sorted(set(tuple(sorted(e)) for e in I.binary_descent(E, bases, target[0]) if e))
    sols = I.binary_solutions(E, bases, target[0])
    points = sorted(c1 | (c2 << k) | (c3 << 2 * k) for c1, c2, c3 in sols)
    rational = I.binary_decompositions(E, bases, target, sols)
    return dict(family='m3-direct', variables=3 * k, equations=[list(e) for e in equations], points=points,
                system_hash=system_hash('binary-m3-direct', 3 * k, equations),
                rational_unordered_multisets=rational)


def m3_split_points(E, basis, target_x):
    """Every (c1, c2, c3, z): S_3(x1, x2, z) = 0 = S_3(z, x3, xR), z in F_(2^n)."""
    F = E.F
    k, n = len(basis), F.N
    V = I.span(basis)
    out = []
    after = {}
    for z in range(1 << n):
        after[z] = [l for l in range(1 << k) if I.s3_binary_value(F, E.b, z, V[l], target_x) == 0]
    for i in range(1 << k):
        for j in range(1 << k):
            for z in range(1 << n):
                if after[z] and I.s3_binary_value(F, E.b, V[i], V[j], z) == 0:
                    for l in after[z]:
                        out.append(i | (j << k) | (l << 2 * k) | (z << 3 * k))
    return sorted(out)


def m3_split_system(E, basis, target):
    k, n = len(basis), E.F.N
    equations, r = binary_split(E, basis, 3, target[0])
    equations = sorted(set(tuple(sorted(e)) for e in equations if e))
    points = m3_split_points(E, basis, target[0])
    return dict(family='m3-split', variables=r, equations=[list(e) for e in equations], points=points,
                system_hash=system_hash('binary-m3-split', r, equations))


# ---------------------------------------------------------------------------
# certification
# ---------------------------------------------------------------------------


def certify(system, max_degree: int, max_columns: int, max_bytes: int):
    r = system['variables']
    eqs = system['equations']
    t0 = time.time()
    try:
        prof = local_profile(eqs, r, system['points'], max_degree=max_degree, max_columns=max_columns,
                             max_bytes=max_bytes)
    except ClosureBudgetExceeded as error:  # defensive; local_profile already catches
        prof = dict(status='capped', exact_lfd=None, cap_reason=str(error))
    seconds = time.time() - t0
    d_in = max(m.bit_count() for e in eqs for m in e)
    out = dict(status=prof['status'], exact_lfd=prof['exact_lfd'], lfd_lower_bound=prof.get('lfd_lower_bound'),
               strict_falls=prof.get('strict_falls'), scanned_through=prof.get('scanned_through'),
               quotient_ranks=prof.get('quotient_ranks'), local_checks=prof.get('local_checks'),
               cap_reason=prof.get('cap_reason'), input_degree=d_in, equations=len(eqs),
               closure_seconds=round(seconds, 3))
    if prof['exact_lfd'] is not None:
        D_cert = prof['exact_lfd'] if prof['status'] == 'exact-unit' else prof['local_upper_bound']
        out['certificate_degree'] = D_cert
        out['macaulay_at_certificate'] = macaulay_size(eqs, r, D_cert)
        D_cost = max(prof['exact_lfd'], d_in)
        out['cost_degree'] = D_cost
        out['macaulay_at_cost_degree'] = macaulay_size(eqs, r, D_cost)
        # dim of V_D* modulo the field equations = |M_<=D*| - codim
        q = prof['quotient_ranks'][D_cert]
        out['dim_closure_at_certificate'] = q
        out['codim_closure_at_certificate'] = squarefree_count(r, D_cert) - q
        assert prof['status'] == 'exact-unit' or out['codim_closure_at_certificate'] == len(system['points'])
    out['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    return out


def bounds(family: str, k: int, n: int, r: int):
    b = dict(kernel_r_plus_1=r + 1)
    if family == 'm2':
        b['kernel_k_plus_l_plus_1'] = 2 * k + 1
        b['written_min_k_l_plus_2'] = k + 2
    elif family == 'm3-direct':
        b['kernel_mk_plus_1_analogue'] = 3 * k + 1  # = r + 1 here
    return b
