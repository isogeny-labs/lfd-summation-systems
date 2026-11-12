"""Homogeneous FFD in F2[x]/(x_i^2), modulo Koszul AND Frobenius rows.

Inputs are the distinct nonzero normalized descended coordinate polynomials.
Boolean field polynomials are represented by the quotient, not appended as
zero top forms. This is an explicitly chosen convention, not a claim that it
is interchangeable with every filtered definition in the literature.
"""
from __future__ import annotations

from itertools import combinations
import math

from . import elim
from .ic_boolean_closure import ClosureBudgetExceeded, _clean


def _degree_masks(r,d):
    if d<0 or d>r:
        return []
    return [sum(1<<i for i in c) for c in combinations(range(r),d)]


def _multiply(top,monomial):
    # In the associated graded/truncated algebra a repeated variable KILLS
    # the monomial. This differs from Boolean idempotent multiplication.
    return [mask|monomial for mask in top if not mask&monomial]


def truncated_ffd(equations,variables,*,max_columns=16384,max_bytes=128_000_000):
    equations=sorted(set(row for row in _clean(equations,variables) if row))
    if any(row==(0,) for row in equations):
        return dict(status='constant-input',exact_ffd=None,convention='F2 square-zero Koszul/Frobenius',degrees=[])
    degrees=[max(m.bit_count() for m in row) for row in equations]
    tops=[[m for m in row if m.bit_count()==d] for row,d in zip(equations,degrees)]
    records=[]
    for D in range(min(degrees,default=0),variables+max(degrees,default=0)+1):
        components=[_degree_masks(variables,D-d) for d in degrees]
        columns=[(i,mask) for i,monomials in enumerate(components) for mask in monomials]
        if len(columns)>max_columns:
            return dict(status='capped',exact_ffd=None,degrees=degrees,scanned=records,
                        cap_reason=f'{len(columns)} syzygy columns exceed cap {max_columns}')
        col_index={key:j for j,key in enumerate(columns)}
        targets=_degree_masks(variables,D)
        target_index={mask:j for j,mask in enumerate(targets)}
        images=[sum(1<<target_index[t] for t in _multiply(tops[i],mask)) for i,mask in columns]
        image_rows=[[target_index[t] for t in _multiply(tops[i],mask)] for i,mask in columns]
        image_matrix=elim.pack_rows(image_rows,len(targets))
        if image_matrix.nbytes>max_bytes:
            return dict(status='capped',exact_ffd=None,degrees=degrees,scanned=records,cap_reason='syzygy matrix byte cap')
        rank=elim.gf2_echelon(image_matrix,len(targets))[0]
        kernel=len(columns)-rank
        trivial=[]
        for i in range(len(tops)):
            for monomial in _degree_masks(variables,D-2*degrees[i]):
                trivial.append([col_index[(i,m)] for m in _multiply(tops[i],monomial)])
            for j in range(i+1,len(tops)):
                for monomial in _degree_masks(variables,D-degrees[i]-degrees[j]):
                    trivial.append([col_index[(i,m)] for m in _multiply(tops[j],monomial)]+
                                   [col_index[(j,m)] for m in _multiply(tops[i],monomial)])
        # Directly verify that every proposed trivial row lies in the kernel.
        for row in trivial:
            value=0
            for column in row:value ^= images[column]
            assert value==0,'A purported trivial row is not a syzygy'
        if len(trivial)*max(1,(len(columns)+63)//64)*8>max_bytes:
            return dict(status='capped',exact_ffd=None,degrees=degrees,scanned=records,cap_reason='trivial-module byte cap')
        trivial_matrix=elim.pack_rows(trivial,len(columns))
        trivial_rank=elim.gf2_echelon(trivial_matrix,len(columns))[0]
        assert trivial_rank<=kernel
        record=dict(degree=D,coefficient_columns=len(columns),image_rank=rank,
                    syzygy_dimension=kernel,trivial_dimension=trivial_rank,
                    nontrivial_dimension=kernel-trivial_rank)
        records.append(record)
        if kernel>trivial_rank:
            return dict(status='exact',exact_ffd=D,degrees=degrees,scanned=records)
    return dict(status='infinite',exact_ffd=None,degrees=degrees,scanned=records)
