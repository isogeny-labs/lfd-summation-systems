#!/usr/bin/env python3
"""Independent exact Wronskian/cofactor and infinity/local-order regression.

This is an exact SymPy regression of one explicit four-monomial circuit,
not a replacement for the arbitrary-circuit Lean theorem.
"""
from sympy import Matrix, Poly, Symbol, cancel, diff, factor, prod, S

X = Symbol("X")
tested = 0
for c in (2, 5, 8, -3):
    z = (X + c)/(X - 1)
    functions = [X*z, z, X, S.One]
    exponents = [(1,1), (0,1), (1,0), (0,0)]
    D = X*(X+c)*(X-1)
    assert cancel(functions[0]-functions[1]-functions[2]-c) == 0
    cofactors = []
    for omit in range(4):
        ff = [functions[j] for j in range(4) if j != omit]
        W = Matrix([[diff(f, X, k) for f in ff] for k in range(3)]).det()
        A = factor(cancel(D**3 * W / prod(ff)))
        p = Poly(A,X)
        assert p != 0
        assert p.degree() <= 6
        cofactors.append(p)
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            A, B = cofactors[i], cofactors[j]
            gamma = cancel(functions[i]/functions[j] * B.as_expr()/A.as_expr())
            assert gamma != 0 and not gamma.has(X), (c,i,j,gamma)
            assert exponents[i][0]-exponents[j][0] == A.degree()-B.degree()
            def ord_one(p):
                order = 0
                while p.eval(1) == 0:
                    p = p.quo(Poly(X-1,X))
                    order += 1
                return order
            assert -(exponents[i][1]-exponents[j][1]) == ord_one(A)-ord_one(B)
            tested += 1
print(f"PASS: {tested} exact cofactor ratios and both pole comparisons over four parameters")
