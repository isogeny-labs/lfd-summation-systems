"""Small exact certificate checks for the accompanying FFD/LFD proof note.

The arbitrary-parameter results are proved mathematically in docs/IC_FALL_PROOFS.md.
These routines verify symbolic identities and small ordinary-ring closures;
they do not infer a general theorem by fitting or enumerating examples.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from . import elim
from .field import is_prime


def finite_quotient_lfd_bound(standard_degree: int, input_degree: int) -> int:
    """Requires an explicit degree-compatible zero-dimensional GB subset of F."""
    if standard_degree < 0 or input_degree < 0:
        raise ValueError("Degrees must be nonnegative.")
    return max(standard_degree + 1, input_degree)


def pkm_lfd_bounds(m: int, k: int) -> dict:
    """Bounds for direct PKM chains, 2^k | p-1, p>3, with no auxiliary variables.

    normalized is a different generating presentation from raw. The raw bound
    retains the exact input degree m*2^(m-1); preprocessing is not free.
    Homogeneous FFD means syzygies modulo Koszul syzygies in the ordinary ring.
    The truncated finite-field convention agrees here when p exceeds the
    stated degree. The normalized FFD bound requires a nonconstant generator.
    """
    if m < 2 or k < 1:
        raise ValueError("Expected m >= 2 and k >= 1.")
    variables = m*k
    raw_degree_bound = m*2**(m-1)
    return dict(m=m, k=k, variables=variables, quotient_dimension=2**variables,
                normalized_lfd_upper_bound=variables+1,
                raw_input_degree_upper_bound=raw_degree_bound,
                raw_lfd_upper_bound=max(variables+1, raw_degree_bound),
                raw_homogeneous_ffd=raw_degree_bound,
                normalized_homogeneous_ffd_upper_bound=(4 if m == 2 and k >= 2
                                                       else m*min(k, m)+1),
                normalized_ffd_bound_requires_nonconstant_generator=True)


def m2_normalized_homogeneous_ffd(prime: int, a: int, target_x: int, k: int) -> int:
    """Exact homogeneous-syzygy FFD of normalized two-summand PKM chains.

    Requires k>=2, p>3, and an actual subgroup of size 2^k in F_p^*. This is
    not an affine strict-fall measurement or an LFD/solving-degree formula.
    """
    if prime <= 3 or not is_prime(prime) or k < 2 or (prime-1) % (1 << k):
        raise ValueError("Require prime p>3, k>=2, and 2^k dividing p-1.")
    return 3 if (3*target_x*target_x+a) % prime == 0 else 4


def symbolic_identities() -> dict[str, bool]:
    import sympy as s
    x, y, a, c = s.symbols("x y a c")
    Cx, Cy, h = x*x-1, y*y-1, x*y-c
    raw = x*x*y*y-2*a*x*y+a*a
    A = y*h-x*Cy                 # x-c*y, formed within degree 3
    B = x*h-y*Cx                 # y-c*x, formed within degree 3
    normalized = 1+a*a-2*a*x*y
    identities = {
        "A_equals_x_minus_c_y": A-(x-c*y),
        "B_equals_y_minus_c_x": B-(y-c*x),
        "linear_combination": A+c*B-(1-c*c)*x,
        "unit_from_linear_fall": x*((A+c*B)/(1-c*c))-Cx-1,
        "raw_chain_reduction": raw-y*y*Cx-Cy-normalized,
        "normalized_equals_scaled_h": normalized-(-2*a*h.subs(c, (1+a*a)/(2*a))),
        "c_squared_minus_one": ((1+a*a)/(2*a))**2-1-(a*a-1)**2/(4*a*a),
    }
    u, v, w, z, alpha, beta, t = s.symbols("u v w z alpha beta t")
    q = alpha*u*w+beta*(u*z+v*w)+v*z
    multiplication = s.Matrix([[beta, alpha, 0, 0], [1, beta, 0, 0],
                               [0, 0, beta, alpha], [0, 0, 1, beta]])
    identities.update({
        "m2_linear_annihilator_determinant": multiplication.det()-(alpha-beta**2)**2,
        "m2_derivative_determinant": (-2*(t*t+a))-(-2*t)**2+2*(3*t*t+a),
        "m2_degree_four_syzygy": u*v*q-u*u*(alpha*v*w+beta*v*z)-v*v*(beta*u*w+u*z),
        "m2_degenerate_degree_three_syzygy": (v-beta*u)*q.subs(alpha, beta**2)
            -(v*v-beta*beta*u*u)*(z+beta*w),
    })
    return {name: s.cancel(s.expand(value)) == 0 for name, value in identities.items()}


def compositions(total, variables):
    if variables == 1:
        yield (total,)
    else:
        for i in range(total+1):
            for rest in compositions(total-i, variables-1):
                yield (i,) + rest


@dataclass
class DegreeClosure:
    """The exact ordinary-polynomial vector space V_(F,D), not a quotient proxy."""
    variables: tuple
    prime: int
    degree_bound: int
    monomials: tuple
    basis: np.ndarray
    pivots: tuple
    rounds: int

    @property
    def rank(self):
        return len(self.basis)

    def contains(self, expression):
        import sympy as s
        poly = s.Poly(expression, *self.variables, modulus=self.prime)
        if poly.is_zero:
            return True
        if poly.total_degree() > self.degree_bound:
            return False
        index = {m: i for i, m in enumerate(self.monomials)}
        row = np.zeros(len(self.monomials), dtype=np.int64)
        for m, c in poly.terms():
            row[index[m]] = int(c) % self.prime
        stacked = np.vstack([self.basis, row])
        return elim.fp_echelon(stacked, self.prime)[0] == self.rank

    def low_degree_dimension(self, degree):
        # The basis and pivot columns use a descending degree-compatible order.
        return sum(sum(self.monomials[pivot]) <= degree for pivot in self.pivots)


def degree_closure(generators, variables, prime, degree_bound, max_columns=20000):
    """Construct V_(F,D) by exact echelonization and repeated low-degree products.

    Multipliers of newly fallen rows are included. No field/chain reduction is
    performed before checking the ordinary total-degree budget.
    """
    import sympy as s
    variables = tuple(variables)
    if not variables or degree_bound < 0 or not is_prime(prime) or prime >= 2**31:
        raise ValueError("Invalid variables, degree, or finite-field characteristic.")
    n = len(variables)
    if math.comb(n+degree_bound, n) > max_columns:
        raise ValueError("Exact-closure verification column budget exceeded.")
    monomials = tuple(sorted((m for d in range(degree_bound+1)
                             for m in compositions(d, n)), key=lambda m: (sum(m), m), reverse=True))
    index = {m: i for i, m in enumerate(monomials)}
    rows = []
    for expression in generators:
        poly = s.Poly(expression, *variables, modulus=prime)
        if poly.is_zero or poly.total_degree() > degree_bound:
            continue
        row = np.zeros(len(monomials), dtype=np.int64)
        for m, c in poly.terms():
            row[index[m]] = int(c) % prime
        rows.append(row)
    A = np.asarray(rows, dtype=np.int64).reshape((-1, len(monomials)))
    rounds = 0
    while True:
        out = []
        rank, pivots, _ = elim.fp_echelon(A, prime, out)
        basis = out[0]
        products = []
        for row, pivot in zip(basis, pivots):
            if sum(monomials[pivot]) >= degree_bound:
                continue
            for j in range(n):
                result = np.zeros(len(monomials), dtype=np.int64)
                for column in np.flatnonzero(row):
                    m = list(monomials[column])
                    m[j] += 1
                    result[index[tuple(m)]] = row[column]
                products.append(result)
        if not products:
            break
        A = np.vstack([basis, *products])
        next_rank = elim.fp_echelon(A, prime)[0]
        rounds += 1
        if next_rank == rank:
            break
    return DegreeClosure(variables, prime, degree_bound, monomials, basis, tuple(pivots), rounds)


def exact_family(prime, a):
    """Nonsingular E:y^2=x^3+a*x, target (0,0), coordinate subgroup +/-1."""
    import sympy as s
    if prime <= 3 or not is_prime(prime) or a % prime in (0, 1, prime-1):
        raise ValueError("Require prime p>3 and a nonzero with a != +/-1.")
    a %= prime
    x, y = s.symbols("x y")
    c = (1+a*a)*pow(2*a, -1, prime) % prime
    assert (c*c-1) % prime != 0
    chain = [x*x-1, y*y-1]
    return dict(variables=(x, y), prime=prime, a=a, b=0, target=(0, 0), c=c,
                normalized=chain+[x*y-c],
                raw=chain+[x*x*y*y-2*a*x*y+a*a])
