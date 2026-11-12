"""Exact iterated ordinary-degree closure modulo explicit Boolean equations.

Normal-form multiplication is charged BEFORE Boolean reduction. The rows
represent the image of V_(F,D) under the degree-nonincreasing normal form
for {x_i^2-x_i}. This is not a static Macaulay or linear-solving stopping rule.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math

import numpy as np

from . import elim


class ClosureBudgetExceeded(ValueError):
    """An exact computation was not completed within the explicit resource cap."""


def monomial_masks(variables, degree, max_columns=16384):
    if not 1 <= variables <= 60 or degree < 0:
        raise ValueError('Require 1<=variables<=60 and a nonnegative degree.')
    count = sum(math.comb(variables,d) for d in range(min(variables,degree)+1))
    if count > max_columns:
        raise ClosureBudgetExceeded(f'{count} normal-form columns exceed cap {max_columns}')
    return tuple(sorted((sum(1<<i for i in c)
                         for d in range(min(variables,degree)+1)
                         for c in combinations(range(variables),d)),
                        key=lambda mask:(mask.bit_count(),mask),reverse=True))


def _clean(equations, variables):
    result = []
    for equation in equations:
        parity = set()
        for mask in equation:
            if not isinstance(mask,int) or not 0<=mask<1<<variables:
                raise ValueError('A monomial mask lies outside the variable set.')
            parity.symmetric_difference_update([mask])
        result.append(tuple(sorted(parity)))
    return result


def _support(row,columns):
    return np.flatnonzero(np.unpackbits(row.view(np.uint8),bitorder='little')[:columns])


@dataclass
class BooleanClosure:
    variables: int
    degree_bound: int
    masks: tuple[int,...]
    basis: np.ndarray
    pivots: tuple[int,...]
    rank_history: tuple[int,...]

    @property
    def quotient_rank(self):
        return len(self.basis)

    @property
    def rank(self):
        return math.comb(self.variables+self.degree_bound,self.variables)-len(self.masks)+self.quotient_rank

    @property
    def contains_one(self):
        return len(self.masks)-1 in self.pivots

    def low_degree_dimension(self,degree):
        if degree<0:
            return 0
        degree=min(degree,self.degree_bound)
        ordinary=math.comb(self.variables+degree,self.variables)
        standard=sum(mask.bit_count()<=degree for mask in self.masks)
        quotient=sum(self.masks[p].bit_count()<=degree for p in self.pivots)
        return ordinary-standard+quotient

    def contains(self,equation):
        equation=_clean([equation],self.variables)[0]
        index={mask:j for j,mask in enumerate(self.masks)}
        if any(mask not in index for mask in equation):
            return False
        row=elim.pack_rows([[index[mask] for mask in equation]],len(index))[0]
        for source,pivot in zip(self.basis,self.pivots):
            if (int(row[pivot//64])>>(pivot%64))&1:
                row ^= source
        return not np.any(row)


def _products(space,max_bytes,done=None):
    index={mask:j for j,mask in enumerate(space.masks)}
    todo=[]
    for row,pivot in zip(space.basis,space.pivots):
        if space.masks[pivot].bit_count()>=space.degree_bound:
            continue
        key=row.tobytes()
        if done is not None:
            if key in done:
                continue
            done.add(key)
        todo.append(row)
    shape=(len(todo)*space.variables,space.basis.shape[1])
    if (shape[0]+space.quotient_rank)*shape[1]*8>max_bytes:
        raise ClosureBudgetExceeded('Boolean closure product matrix exceeds byte cap')
    products=np.zeros(shape,dtype=np.uint64)
    for i,row in enumerate(todo):
        support=[space.masks[j] for j in _support(row,len(space.masks))]
        for variable in range(space.variables):
            # Several monomials can collide under idempotent multiplication:
            # accumulate with XOR rather than assigning a coefficient.
            destination=np.asarray([index[mask | (1<<variable)] for mask in support],dtype=np.int64)
            np.bitwise_xor.at(products[i*space.variables+variable],destination>>6,
                np.left_shift(np.uint64(1),(destination&63).astype(np.uint64)))
    return products


def boolean_degree_closure(equations,variables,degree_bound,*,max_columns=16384,max_bytes=256_000_000):
    equations=_clean(equations,variables)
    masks=monomial_masks(variables,degree_bound,max_columns)
    index={mask:j for j,mask in enumerate(masks)}
    inputs=[[index[mask] for mask in equation] for equation in equations
            if equation and max(mask.bit_count() for mask in equation)<=degree_bound]
    matrix=elim.pack_rows(inputs,len(masks))
    history=[]
    done=set()
    while True:
        if matrix.nbytes>max_bytes:
            raise ClosureBudgetExceeded('Boolean closure input matrix exceeds byte cap')
        output=[]
        rank,pivots,_=elim.gf2_echelon(matrix,len(masks),output)
        history.append(rank)
        space=BooleanClosure(variables,degree_bound,masks,output[0],tuple(pivots),tuple(history))
        if space.contains_one:
            # A unit at D generates all monomials of degree <=D at that budget.
            if len(masks)*max(1,(len(masks)+63)//64)*8>max_bytes:
                raise ClosureBudgetExceeded('Unit closure basis exceeds byte cap')
            diagonal=elim.pack_rows([[i] for i in range(len(masks))],len(masks))
            return BooleanClosure(variables,degree_bound,masks,diagonal,tuple(range(len(masks))),tuple(history))
        products=_products(space,max_bytes,done)
        if not len(products) or len(history)>1 and history[-1]==history[-2]:
            return space
        matrix=np.vstack([space.basis,products])


def check_containing_space(space,equations,*,max_bytes=256_000_000):
    """Independent inclusion/stability check; no trust in the row-generation path."""
    equations=_clean(equations,space.variables)
    expected=monomial_masks(space.variables,space.degree_bound,max(len(space.masks),1))
    if space.masks!=expected or space.basis.dtype!=np.uint64 or space.basis.ndim!=2:
        return False
    if space.basis.shape[1]!=max(1,(len(expected)+63)//64):
        return False
    if len(expected)%64 and np.any(space.basis[:,-1] >> np.uint64(len(expected)%64)):
        return False
    output=[]
    rank,pivots,_=elim.gf2_echelon(space.basis,len(expected),output)
    if rank!=len(space.basis):
        return False
    canonical=BooleanClosure(space.variables,space.degree_bound,expected,output[0],tuple(pivots),())
    if any(equation and max(m.bit_count() for m in equation)<=space.degree_bound
           and not canonical.contains(equation) for equation in equations):
        return False
    products=_products(canonical,max_bytes)
    return elim.gf2_echelon(np.vstack([canonical.basis,products]),len(expected))[0]==rank


def closure_profile(equations,variables,*,max_degree=None,max_columns=16384,max_bytes=256_000_000):
    """Exact LFD only if the proved r+1 bound is scanned, or a first unit is found.

    For proper ideals, the last strict fall is used only after reaching r+1.
    A resource cap returns an explicit partial profile, never a fabricated label.
    """
    equations=_clean(equations,variables)
    upper=variables+1
    limit=upper if max_degree is None else min(max_degree,upper)
    spaces=[]
    falls=[]
    cap=None
    for degree in range(limit+1):
        try:
            current=boolean_degree_closure(equations,variables,degree,
                max_columns=max_columns,max_bytes=max_bytes)
        except ClosureBudgetExceeded as error:
            cap=str(error)
            break
        if spaces and current.low_degree_dimension(degree-1)>spaces[-1].rank:
            falls.append(degree)
        spaces.append(current)
        if current.contains_one:
            return dict(status='exact-unit',exact_lfd=degree,lfd_lower_bound=degree,
                        lfd_upper_bound=degree,strict_falls=falls,
                        first_affine_strict_fall=min(falls,default=None),
                        quotient_ranks=[s.quotient_rank for s in spaces],
                        ordinary_ranks=[s.rank for s in spaces],
                        scanned_through=degree,cap_reason=None),spaces
    complete=len(spaces)==upper+1
    last=max(falls,default=0)
    return dict(status='exact-proper' if complete else 'capped',exact_lfd=last if complete else None,
                lfd_lower_bound=last,lfd_upper_bound=last if complete else upper,
                strict_falls=falls,first_affine_strict_fall=min(falls,default=None),
                quotient_ranks=[s.quotient_rank for s in spaces],ordinary_ranks=[s.rank for s in spaces],
                scanned_through=len(spaces)-1,cap_reason=cap or (None if complete else 'configured degree cap')) ,spaces
