"""Exact row-span affine elimination in the Boolean ring over GF(2).

A polynomial is a set of squarefree monomial bitmasks (0 is the constant).
No oracle, assignment search or Gröbner calculation is used in reduction.
"""
from __future__ import annotations

from dataclasses import dataclass


def normalize(terms):
    out = set()
    for m in terms:
        if m in out:
            out.remove(m)
        else:
            out.add(m)
    return frozenset(out)


def substitute(poly, pivot, rhs):
    """Replace the singleton monomial pivot by the affine polynomial rhs."""
    out = set()
    for monomial in poly:
        terms = ((monomial ^ pivot) | term for term in rhs) if monomial & pivot else (monomial,)
        for term in terms:
            if term in out:
                out.remove(term)
            else:
                out.add(term)
    return frozenset(out)


def row_basis(equations):
    """Span-preserving echelon form, higher degree monomials pivoted first."""
    monomials = sorted(set().union(*equations) if equations else set(),
                       key=lambda m: (m.bit_count(), m))
    positions = {m:i for i,m in enumerate(monomials)}
    pivots, xors = {}, 0
    for p in equations:
        row = sum(1 << positions[m] for m in p)
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
            xors += 1
    rows=[]
    for pivot in sorted(pivots, reverse=True):
        row=pivots[pivot]
        terms=[]
        while row:
            low=row & -row
            terms.append(monomials[low.bit_length()-1])
            row ^= low
        rows.append(frozenset(terms))
    return rows, xors


@dataclass
class Reduction:
    equations: list[frozenset[int]]
    variables: int
    recipes: list[frozenset[int]]
    inconsistent: bool
    statistics: dict

    def recover(self, mask):
        """Lift one residual Boolean assignment to the original variables."""
        out=0
        for i,p in enumerate(self.recipes):
            if sum((mask & m)==m for m in p) & 1:
                out |= 1 << i
        return out


def reduce_affine(equations, variables):
    equations=[normalize(p) for p in equations]
    if any(m < 0 or m >= (1 << variables) for p in equations for m in p):
        raise ValueError('monomial outside the variable range')
    recipes=[frozenset((1 << i,)) for i in range(variables)]
    active=list(range(variables))
    stats=dict(input_variables=variables, eliminated=0, row_rounds=0,
               row_xors=0, substitutions=0, row_ranks=[], row_columns=[],
               row_input_rows=[], affine_rows=[])
    while True:
        stats['row_input_rows'].append(len(equations))
        stats['row_columns'].append(len(set().union(*equations)) if equations else 0)
        equations,xors=row_basis(equations)
        stats['row_rounds']+=1
        stats['row_xors']+=xors
        stats['row_ranks'].append(len(equations))
        if frozenset((0,)) in equations:
            stats['residual_variables']=0
            return Reduction([frozenset((0,))],0,[],True,stats)
        affine=[p for p in equations if all(m.bit_count()<=1 for m in p)]
        stats['affine_rows'].append(len(affine))
        if not affine:
            break
        relation=max(affine,key=lambda p: (max(p),tuple(sorted(p))))
        pivot=max(relation)
        rhs=relation ^ {pivot}
        equations=[substitute(p,pivot,rhs) for p in equations]
        recipes=[substitute(p,pivot,rhs) for p in recipes]
        active.remove(pivot.bit_length()-1)
        stats['eliminated']+=1
        stats['substitutions']+=len(equations)+len(recipes)

    positions={old:new for new,old in enumerate(active)}
    def compact(poly):
        return frozenset(sum(1 << positions[j] for j in active if m >> j & 1) for m in poly)
    stats['residual_variables']=len(active)
    return Reduction([compact(p) for p in equations],len(active),
                     [compact(p) for p in recipes],False,stats)


def trace_constraints(E,basis,r,formulation):
    """Valid only on the rational-point PDP domain (and its finite-u split)."""
    ell=len(basis)
    def trace(blocks,constant):
        p=set()
        if E.F.trace(constant):
            p.add(0)
        for offset,B in blocks:
            for j,v in enumerate(B):
                if E.F.trace(v):
                    p.add(1 << (offset+j))
        return sorted(p)
    blocks=[(i*ell,basis) for i in range(3)]
    if formulation=='original':
        return [trace(blocks,r)]
    if formulation=='split':
        u=(3*ell,[1 << j for j in range(E.F.N)])
        return [trace(blocks[:2]+[u],E.a),trace([blocks[2],u],r ^ E.a)]
    raise ValueError(formulation)
