"""Complete small-field PDP solving with optional exact affine preprocessing.

The ordinary SymPy GB backend is unchanged. The new preprocessing runs in
the Boolean ring, then uses a bijective affine map to recover all variables.
"""
from __future__ import annotations

import signal
import time

from . import ic_dreg as I
from . import semaev_matched as M
from .boolean_affine import Reduction, reduce_affine, trace_constraints


def lifted_points(E, values):
    points=set()
    for x in values:
        p=E.lift_x(x)
        if p is not None:
            points.update((p,E.neg(p)))
    return sorted(points)


def build_system(E,basis,r,formulation):
    """Same equations as the audited builder; membership scans V, not k²."""
    F,ell=E.F,len(basis)
    if I._f2_rank(basis)!=ell:
        raise ValueError('basis must be independent')
    if formulation=='original':
        eqs=I.binary_descent(E,[basis]*3,r)
        nv=3*ell
    elif formulation=='split':
        add=I._bp_add
        mul=lambda p,q:I._bp_mul(F,p,q)
        sq=lambda p:I._bp_sq(F,p)
        xs=[I._bp_block(basis,i*ell) for i in range(3)]
        u=I._bp_block([1<<j for j in range(F.N)],3*ell)
        def s3(x,y,z):
            return add(add(sq(add(add(mul(x,y),mul(x,z)),mul(y,z))),
                           mul(mul(x,y),z)),I._bp_const(E.b))
        polys=[s3(xs[0],xs[1],u),s3(u,xs[2],I._bp_const(r))]
        eqs=[sorted(m for m,c in p.items() if c>>j&1) for p in polys for j in range(F.N)]
        eqs=[p for p in eqs if p]
        nv=3*ell+F.N
    else:
        raise ValueError(formulation)
    V=I.span(basis)
    anf=[int(E.lift_x(x) is None) for x in V]
    for j in range(ell):
        for mask in range(1<<ell):
            if mask>>j&1:
                anf[mask]^=anf[mask^(1<<j)]
    for block in range(3):
        member=[m<<(block*ell) for m,c in enumerate(anf) if c]
        if member:
            eqs.append(member)
    return eqs,nv


def infinity_solutions(E,basis,r):
    V=I.span(basis)
    return {(x,x,r) for x in V if E.lift_x(x) is not None} if r in V else set()


class GBTimeout(Exception):
    pass


def solve(E,basis,R,formulation='original',precondition=True,gb_timeout=10):
    """Return (all ordered x triples, charged phase metrics); no oracle input."""
    import sympy as sp
    start=time.perf_counter()
    eqs,nv=build_system(E,basis,R[0],formulation)
    input_seconds=time.perf_counter()-start
    input_stats=M.polynomial_stats(eqs,nv)
    start=time.perf_counter()
    if precondition:
        augmented=eqs+trace_constraints(E,basis,R[0],formulation)
        reduction=reduce_affine(augmented,nv)
    else:
        reduction=Reduction([frozenset(p) for p in eqs],nv,
                            [frozenset((1<<j,)) for j in range(nv)],False,{})
    preprocess=time.perf_counter()-start
    residual,nr=reduction.equations,reduction.variables
    start=time.perf_counter()
    variables=sp.symbols(f'z0:{nr}')
    do_gb=not reduction.inconsistent and bool(residual) and nr>0
    if do_gb:
        exprs=[sp.Add(*(sp.Mul(*(variables[j] for j in range(nr) if m>>j&1))
                         for m in p)) for p in residual]
        exprs += [x*x-x for x in variables]
    encoding=time.perf_counter()-start
    G=None
    start=time.perf_counter()
    if do_gb:
        def alarm(_signum,_frame):
            raise GBTimeout()
        oldhandler=signal.signal(signal.SIGALRM,alarm)
        signal.setitimer(signal.ITIMER_REAL,gb_timeout)
        try:
            G=sp.groebner(exprs,*variables,modulus=2,order='grevlex',method='buchberger')
        finally:
            signal.setitimer(signal.ITIMER_REAL,0)
            signal.signal(signal.SIGALRM,oldhandler)
    gb=time.perf_counter()-start
    start=time.perf_counter()
    assignments=[]
    if reduction.inconsistent:
        pass
    elif G is None:
        if not any(0 in p for p in residual):
            assignments=range(1<<nr)
    else:
        lex=G if list(G)==[1] else G.fglm(order='lex')
        constraints=[]
        for p in lex.polys:
            terms=[sum(1<<j for j,e in enumerate(m) if e) for m,c in p.terms() if int(c)&1]
            support=0
            for m in terms:
                support |= m
            constraints.append((support,terms))
        def visit(j,mask,assigned):
            if any((support&assigned)==support and M.eval_boolean(p,mask)
                   for support,p in constraints):
                return
            if j<0:
                assignments.append(mask)
                return
            visit(j-1,mask,assigned|(1<<j))
            visit(j-1,mask|(1<<j),assigned|(1<<j))
        visit(nr-1,0,0)
    tuples={M.decode(reduction.recover(mask) if precondition else mask,basis) for mask in assignments}
    if formulation=='split':
        tuples |= infinity_solutions(E,basis,R[0])
    recovery=time.perf_counter()-start
    metrics=dict(setup_seconds=input_seconds+encoding,input_seconds=input_seconds,
                 encoding_seconds=encoding,preprocess_seconds=preprocess,
                 gb_seconds=gb,recovery_seconds=recovery,
                 total_seconds=input_seconds+encoding+preprocess+gb+recovery,
                 input=input_stats,residual_variables=nr,
                 residual_equations=len(residual),preprocessing_solved=not do_gb,
                 preprocessing_inconsistent=reduction.inconsistent,
                 preprocessing=reduction.statistics,
                 gb_degree=max((p.total_degree() for p in G.polys),default=0) if G is not None else None,
                 projected_solutions=len(tuples))
    return tuples,metrics
