"""Exact degree-bounded closure certificates for PKM chain systems.

The quotient computation is justified by a degree-compatible Groebner-basis
subsystem already present in the input; see docs/IC_FALL_COUNTEREXAMPLES.md.
Products are allowed only when their *unreduced* degree fits the budget.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math

import numpy as np

from . import elim
from .field import is_prime


@dataclass
class ChainDegreeClosure:
    m: int
    k: int
    prime: int
    degree_bound: int
    masks: tuple[int, ...]
    basis: np.ndarray
    pivots: tuple[int, ...]
    rank_history: tuple[int, ...]

    @property
    def quotient_rank(self):
        return len(self.basis)

    @property
    def rank(self):
        """Dimension of V_(F,D) in the ordinary ring, including chain rows."""
        return (math.comb(self.m*self.k+self.degree_bound, self.degree_bound)
                -len(self.masks)+self.quotient_rank)

    @property
    def contains_one(self):
        return len(self.masks)-1 in self.pivots

    def low_degree_dimension(self, degree):
        if degree < 0:
            return 0
        degree = min(degree, self.degree_bound)
        ordinary = math.comb(self.m*self.k+degree, degree)
        standard = sum(mask.bit_count() <= degree for mask in self.masks)
        quotient = sum(self.masks[j].bit_count() <= degree for j in self.pivots)
        return ordinary-standard+quotient

    def contains_normal_form(self, polynomial):
        index = {mask: j for j,mask in enumerate(self.masks)}
        row = np.zeros(len(self.masks), dtype=np.int64)
        for mask,c in polynomial.items():
            c %= self.prime
            if not c:
                continue
            if mask not in index:
                return False
            row[index[mask]] = c
        for source,pivot in zip(self.basis,self.pivots):
            if row[pivot]:
                scale = int(row[pivot])*pow(int(source[pivot]),-1,self.prime) % self.prime
                row = (row-scale*source) % self.prime
        return not np.any(row)


def _product_columns(masks, m, k, degree_bound):
    index = {mask: j for j,mask in enumerate(masks)}
    columns = np.full((m*k,len(masks)),-1,dtype=np.int64)
    N = 1 << k
    for source,mask in enumerate(masks):
        if mask.bit_count() >= degree_bound:
            continue
        for j in range(m*k):
            block,bit = divmod(j,k)
            shift = block*k
            exponent = ((mask >> shift)+(1 << bit)) % N
            product = (mask & ~((N-1) << shift)) | (exponent << shift)
            columns[j,source] = index[product]
    return columns


def _products(basis,pivots,masks,columns,degree_bound,done=None):
    todo = []
    for row,pivot in zip(basis,pivots):
        # Do not reduce a degree-(D+1) multiplication and then charge only D.
        if masks[pivot].bit_count() >= degree_bound:
            continue
        if done is not None:
            key = row.tobytes()
            if key in done:
                continue
            done.add(key)
        todo.append(row)
    result = np.zeros((len(todo)*len(columns),len(masks)),dtype=np.int64)
    for i,row in enumerate(todo):
        support = np.flatnonzero(row)
        for j in range(len(columns)):
            destination = columns[j,support]
            assert np.all(destination >= 0)
            # Multiplication by a chain variable permutes the full quotient
            # monomials (addition modulo 2^k in one exponent block).
            result[i*len(columns)+j,destination] = row[support]
    return result


def chain_degree_closure(generators,m,k,prime,degree_bound,input_degrees=None,
                         max_columns=20000,max_cells=25000000):
    """Compute the exact ordinary-ring V_(F,D), represented modulo chains.

    Generators are *normal forms*, dicts {squarefree mask: coefficient}.
    input_degrees optionally records their raw input degrees: a raw generator
    cannot be used below its original degree. All m*k chain equations are
    inputs. There is no iteration cap; strictly increasing finite ranks prove
    termination. The return value includes the ordinary-ring dimension.
    """
    r = m*k
    if (m < 1 or k < 1 or r > 60 or degree_bound < 0 or prime <= 3
            or prime >= 2**31 or not is_prime(prime) or (prime-1) % (1 << k)):
        raise ValueError("Require valid PKM parameters, prime p>3, and 2^k | p-1.")
    count = sum(math.comb(r,d) for d in range(min(r,degree_bound)+1))
    if count > max_columns:
        raise ValueError("Exact chain-closure column budget exceeded.")
    masks = tuple(sorted((sum(1 << j for j in choice)
                          for d in range(min(r,degree_bound)+1)
                          for choice in combinations(range(r),d)),
                         key=lambda mask:(mask.bit_count(),mask),reverse=True))
    index = {mask:j for j,mask in enumerate(masks)}
    cleaned = []
    degrees = []
    for generator in generators:
        if any(not isinstance(mask,int) or mask < 0 or mask >= 1 << r for mask in generator):
            raise ValueError("Normal-form masks must fit the chain variables.")
        normal = {mask:int(c) % prime for mask,c in generator.items() if c % prime}
        cleaned.append(normal)
        degrees.append(max((mask.bit_count() for mask in normal),default=0))
    if input_degrees is not None:
        if len(input_degrees) != len(generators) or any(
                d < normal for d,normal in zip(input_degrees,degrees)):
            raise ValueError("Input degrees must cover the normal-form degrees.")
        degrees = list(input_degrees)
    rows = []
    for generator,degree in zip(cleaned,degrees):
        if degree > degree_bound or not generator:
            continue
        row = np.zeros(count,dtype=np.int64)
        for mask,c in generator.items():
            row[index[mask]] = c
        rows.append(row)
    matrix = np.asarray(rows,dtype=np.int64).reshape((-1,count))
    columns = _product_columns(masks,m,k,degree_bound)
    done = set()
    history = []
    while True:
        if matrix.size > max_cells:
            raise ValueError("Exact chain-closure matrix budget exceeded.")
        out = []
        rank,pivots,_ = elim.fp_echelon(matrix,prime,out)
        basis = out[0]
        history.append(rank)
        if count-1 in pivots:
            # A derived 1 generates every monomial within degree D. This is
            # a justified closure step, not a quotient-only shortcut.
            basis = np.eye(count,dtype=np.int64)
            pivots = list(range(count))
            break
        if len(history) > 1 and history[-1] == history[-2]:
            break
        low_count = sum(masks[j].bit_count() < degree_bound for j in pivots)
        if (rank+low_count*r)*count > max_cells:
            raise ValueError("Exact chain-closure product budget exceeded.")
        products = _products(basis,pivots,masks,columns,degree_bound,done)
        if not len(products):
            break
        matrix = np.vstack([basis,products])
    return ChainDegreeClosure(m,k,prime,degree_bound,masks,basis,tuple(pivots),tuple(history))


def check_closed_exclusion_certificate(closure,generators,input_degrees=None):
    """Check a containing, multiplication-stable subspace that excludes 1.

    This verifies the *lower bound* without trusting how the supplied basis was
    generated. Upper membership is checked by rerunning chain_degree_closure.
    """
    r = closure.m*closure.k
    degree = closure.degree_bound
    if (closure.m < 1 or closure.k < 1 or r > 60 or degree < 0
            or closure.prime <= 3 or not is_prime(closure.prime)
            or (closure.prime-1) % (1 << closure.k)):
        return False
    count = sum(math.comb(r,d) for d in range(min(r,degree)+1))
    if (len(closure.masks) != count or len(set(closure.masks)) != count
            or any(mask < 0 or mask >= 1 << r or mask.bit_count() > degree
                   for mask in closure.masks)
            or closure.masks != tuple(sorted(closure.masks,
                key=lambda mask:(mask.bit_count(),mask),reverse=True))
            or closure.basis.ndim != 2 or closure.basis.shape[1] != count):
        return False
    # Recompute an echelon basis: the degree filtration, independence, and
    # membership checks must not trust externally supplied pivot metadata.
    out = []
    rank,pivots,_ = elim.fp_echelon(closure.basis,closure.prime,out)
    if rank != closure.quotient_rank or count-1 in pivots:
        return False
    canonical = ChainDegreeClosure(closure.m,closure.k,closure.prime,degree,
                                   closure.masks,out[0],tuple(pivots),())
    if input_degrees is None:
        input_degrees = [max((mask.bit_count() for mask,c in g.items() if c % closure.prime),
                            default=0) for g in generators]
    if len(input_degrees) != len(generators) or any(
            d < max((mask.bit_count() for mask,c in g.items() if c % closure.prime),default=0)
            for g,d in zip(generators,input_degrees)):
        return False
    if any(d <= degree and not canonical.contains_normal_form(g)
           for g,d in zip(generators,input_degrees)):
        return False
    columns = _product_columns(closure.masks,closure.m,closure.k,degree)
    products = _products(out[0],pivots,closure.masks,columns,degree)
    return elim.fp_echelon(np.vstack([out[0],products]),closure.prime)[0] == rank


def exclusion_dual(closure):
    """An explicit linear functional annihilating the basis but sending 1 to 1."""
    out = []
    _,pivots,_ = elim.fp_echelon(closure.basis,closure.prime,out)
    if len(closure.masks)-1 in pivots:
        raise ValueError("Cannot exclude 1 from a space that contains it.")
    dual = np.zeros(len(closure.masks),dtype=np.int64)
    dual[-1] = 1
    for row,pivot in reversed(list(zip(out[0],pivots))):
        total = sum(int(row[j])*int(dual[j]) for j in np.flatnonzero(row[pivot+1:])+pivot+1)
        dual[pivot] = -total*pow(int(row[pivot]),-1,closure.prime) % closure.prime
    assert np.all((closure.basis.astype(object) @ dual.astype(object)) % closure.prime == 0)
    return dual
