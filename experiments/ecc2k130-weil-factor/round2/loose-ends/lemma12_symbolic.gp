\\ Symbolic (generic s, q, a, b, c, d) check of Lemma 1 and Lemma 2 of lead2_index/PROOF.md.
\\ NOT a formal proof (Lean was unavailable); a computer-algebra identity check in Z[s,q,a,b,c,d][X]/(X^2 - sX + q).
m = X^2 - s*X + q;
cj(u) = Mod(subst(lift(u), X, s - X), m);          \\ alpha -> beta = s - alpha
al = Mod(X, m); be = s - al;
Dnum(u, v) = u*cj(v) - cj(u)*v;                     \\ = D(u,v)*(alpha - beta)
x1 = a + b*al; y1 = c + d*al;
printf("Lemma 1: x*ybar - xbar*y == -(ad-bc)(alpha-beta) identically: %d\n", Dnum(x1, y1) == -(a*d - b*c)*(al - be));
printf("Lemma 1: D(1,alpha) = -1: %d ; D alternating: %d\n", Dnum(1, al) == -(al - be), Dnum(x1, x1) == 0);
U = vector(62); U[1] = 0; U[2] = 1; for (n = 2, 61, U[n+1] = s*U[n] - q*U[n-1]);   \\ U[n+1] = U_n
ok = 1; for (i = 0, 30, for (j = 0, 30, if (Dnum(al^i, be^j) != U[i+j+1]*(al - be), ok = 0)));
printf("Lemma 2: D(alpha^i, beta^j) = U_{i+j} for 0<=i,j<=30 (generic s,q): %d\n", ok);
ok2 = 1; for (n = 0, 60, if (Dnum(1, al^n) != -U[n+1]*(al - be), ok2 = 0));
printf("Lemma 2: D(1, alpha^n) = -U_n for n<=60: %d\n", ok2);
\\ negative control: a wrong sign / wrong index must fail
printf("negative control D(alpha^2,beta^3) = U_4 (wrong)? %d (expect 0)\n", Dnum(al^2, be^3) == U[5]*(al - be));
\\ Lemma 1 index formula, numerically: [O : Zx+Zy] = |ad-bc| via SNF over 200 random cases (s=-1,q=2 not needed: lattice in Z^2)
ok3 = 1; for (t = 1, 200, my(M = matrix(2,2,i,j,random(41)-20)); if (matdet(M) != 0 && abs(matdet(M)) != prod(k=1,2,matsnf(M)[k]), ok3 = 0));
printf("index = |det| on 200 random lattices: %d\n", ok3);
