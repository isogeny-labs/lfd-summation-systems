#!/usr/bin/env python3
"""Independent exact checks of the N8 omitted-column Wronskian identity.

No Lean imports: SymPy computes the rational derivatives and separately
implements the polynomial numerator recurrence.  Each tested dependence
is the genuine identity x^(u+1) z^v = x^u z^v + x^(u+1) z^(v-1)
+ c x^u z^(v-1), where z=(x+c)/(x-1).
"""
import sympy as sp

x = sp.Symbol("x")
checks = 0

def mono(c, u, v):
    return x**u * ((x+c)/(x-1))**v

def cofactor(c, uv):
    d = x*(x+c)*(x-1)
    columns = []
    for u, v in uv:
        q = u*(x+c)*(x-1) - v*(c+1)*x
        p = sp.Integer(1)
        col = []
        for s in range(len(uv)):
            col.append(p)
            p = sp.expand(d*sp.diff(p,x) + (q-s*sp.diff(d,x))*p)
        columns.append(col)
    return sp.expand(sp.det(sp.Matrix(columns).T))

for c in (sp.Integer(6), sp.Integer(8), sp.Integer(-3)):
    for u in (0, 1):
        for v in (1, 2):
            exponent_family = [(u, v), (u+1, v-1), (u, v-1)]
            replacement = (u+1, v)
            coefficients = [sp.Integer(1), sp.Integer(1), c]
            fs = [mono(c,*uv) for uv in exponent_family]
            g = mono(c,*replacement)
            assert sp.cancel(g-sum(a*f for a,f in zip(coefficients,fs))) == 0
            A0 = cofactor(c, exponent_family)
            W0 = sp.det(sp.Matrix([
                [sp.diff(f,x,s) for f in fs] for s in range(3)]))
            for i in range(3):
                replaced = exponent_family.copy()
                replaced[i] = replacement
                gs = fs.copy()
                gs[i] = g
                Ai = cofactor(c,replaced)
                Wi = sp.det(sp.Matrix([
                    [sp.diff(f,x,s) for f in gs] for s in range(3)]))
                t = sp.Rational(13,7)
                sub = lambda p: sp.cancel(p).subs(x,t)
                assert sub(Wi-coefficients[i]*W0) == 0
                d=x*(x+c)*(x-1)
                normalized=lambda ws, ap, uv: sub(
                    d**3*ws - ap*sp.prod(mono(c,*ej) for ej in uv))
                assert normalized(W0,A0,exponent_family) == 0
                assert normalized(Wi,Ai,replaced) == 0
                assert sub(Ai*g-coefficients[i]*A0*fs[i]) == 0
                checks += 1

print(f"N8 independent genuine omitted-minor checks passed: {checks}")
