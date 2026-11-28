\\ U_ell = (alpha^ell - beta^ell)/(alpha-beta), alpha^2+alpha+2=0; two independent computations
lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }
lucasU2(n) = { my(al = Mod(x, x^2+x+2)); my(t = (al^n - (-1-al)^n)/(2*al+1)); t = lift(t); if (poldegree(t) > 0, error("not rational")); polcoeff(t,0); }
{
for (i = 1, 12, my(ell = prime(i+1)); if (ell == 7, next);
  my(U = lucasU(ell), U2 = lucasU2(ell)); if (U != U2, error("mismatch"));
  my(F = factor(abs(U)), eps = kronecker(-7, ell));
  my(types = vector(#F~, j, my(p = F[j,1]);
     if (p % ell == ell-1 && kronecker(-7,p) == -1, "inert",
     if (p % ell == 1 && kronecker(-7,p) == 1, "split", Str("OTHER(p mod ell=", p%ell, ",(-7/p)=", kronecker(-7,p), ")")))));
  printf("ell=%d eps=%d U=%d factors=%s exps=%s types=%s\n", ell, eps, U, F[,1]~, F[,2]~, types));
}
