\\ K+ = Q(c, s), c = z+1/z, s = sqrt(7(4-c^2)) = sqrt(-7)*(z-1/z) = (2*alpha+1)*(z-1/z).
default(parisizemax, 4000000000);
default(realprecision, 300);
\\ minimal polynomial of c = z + 1/z, in variable y
cpoly(ell) = { my(R = polresultant(polcyclo(ell, 'z), 'z^2 - 'y*'z + 1, 'z)); my(f = factor(R)[1,1]); f/pollead(f); }
\\ returns [P (abs pol of K+ in x, reduced), mapping data]: c and s as polmods mod P
setupKp(ell) = {
  my(fc = cpoly(ell), Kc = nfinit(fc));
  my(eq = rnfequation(Kc, 'x^2 - 7*(4 - 'y^2), 1));  \\ eq[1](x): x_abs = s + k*c ; eq[2] = c as polmod
  my(c0 = eq[2], s0 = Mod('x, eq[1]) - eq[3]*c0);
  my(red = polredbest(eq[1], 1));             \\ red[2]: old root x expressed in new variable
  my(P = red[1], xo = red[2]);
  my(c1 = Mod(subst(lift(c0), 'x, lift(xo)), P), s1 = Mod(subst(lift(s0), 'x, lift(xo)), P));
  if (subst(fc, 'y, c1) != 0, error("c bad")); if (s1^2 != 7*(4 - c1^2), error("s bad"));
  [P, c1, s1];
}

\\ Full setup with maps K+ -> K.  K = Q(zeta, alpha) as absolute field polK (polcompositum), ZE, AL polmods.
\\ Returns [P, c1, s1, xK] where xK is the image in K of the K+ generator x (root of P).
setupBoth(ell) = {
  my(fc = cpoly(ell), Kc = nfinit(fc));
  my(eq = rnfequation(Kc, 'x^2 - 7*(4 - 'y^2), 1));
  my(k = eq[3], red = polredbest(eq[1], 1), P = red[1]);
  my(xrev = modreverse(red[2]));            \\ new root x as polmod in old root X (mod eq[1])
  my(c1 = Mod(subst(lift(eq[2]), 'x, lift(red[2])), P));
  my(s1 = Mod(lift(red[2]), P) - k*c1);
  if (subst(fc, 'y, c1) != 0 || s1^2 != 7*(4 - c1^2), error("K+ data bad"));
  my(C = polcompositum(polcyclo(ell, 'x), 'x^2 + 'x + 2, 1)[1]);
  my(polK = C[1], ZE = C[2], AL = C[3]);
  if (polcyclo(ell, 'x) != 0 && subst(polcyclo(ell,'x),'x,ZE) != 0, error("zeta"));
  if (AL^2 + AL + 2 != 0, error("alpha"));
  my(cK = ZE + 1/ZE, sK = (2*AL + 1)*(ZE - 1/ZE));
  my(XK = sK + k*cK);
  my(xK = subst(lift(xrev), 'x, XK));
  if (subst(P, 'x, xK) != 0, error("embedding K+ -> K failed"));
  if (subst(lift(c1), 'x, xK) != cK || subst(lift(s1), 'x, xK) != sK, error("c,s images wrong"));
  [P, c1, s1, xK, polK, ZE, AL];
}
toK(f, xK) = subst(lift(f), 'x, xK);
