\\ Red-team (redteam-lead1): reproduce lead 1's decide() UNCHANGED, then the same routine with ONE fix:
\\ GP evaluates  [0,...,0]~ == 0  as TRUE, so the original tests  `if (sol == 0, next)`  and
\\ `if (x == 0, next)`  also discard the case where matsolvemod returns the ZERO SOLUTION
\\ (g already antisymmetric, resp. lambda already Phi-positive).  The fix tests type(..) == "t_INT"
\\ (matsolvemod returns the integer 0 only when there is no solution).
\\ usage: gp -q < bugcheck_decide.gp     (reads ../../lead1_polarization/decide_ppav.gp read-only)
\r ../../lead1_polarization/decide_ppav.gp
print("GP check: [0,0]~ == 0 evaluates to ", [0,0]~ == 0);

decide_fixed(J, Phi) = {
  my(cyc = BNF.cyc, gens = BNF.gen, w = BNF.tu[1], tu = BNF.tu[2], fu = BNF.fu, r = #fu);
  my(hexp = vector(r), out = List(), nclass = 0, nprinc = 0, nanti = 0, zsol = 0, zx = 0);
  for (i = 1, r, hexp[i] = dlog(nfeltdiv(NF, cj(fu[i]), fu[i])));
  my(A = Mat(concat([dlog(nfeltdiv(NF, cj(tu), tu))], hexp)));
  my(ker = matsolvemod(A, w, [0]~, 1)[2]);
  my(kunits = vector(#ker[1,], k, my(c = ker[,k], u = nfeltpow(NF, tu, c[1])); for (i = 1, r, u = nfeltmul(NF, u, nfeltpow(NF, fu[i], c[i+1]))); u));
  my(Ms = matrix(#Phi, #kunits, i, k, 0));
  for (k = 1, #kunits, my(s = signrow(kunits[k], Phi)); for (i = 1, #Phi, Ms[i,k] = s[i]));
  forvec(e = vector(#cyc, i, [0, cyc[i]-1]),
    nclass++;
    my(S = idealhnf(NF, 1)); for (i = 1, #cyc, S = idealmul(NF, S, idealpow(NF, gens[i], e[i])));
    my(N = idealmul(NF, J, idealinv(NF, idealmul(NF, S, cjI(S)))));
    my(pr = bnfisprincipal(BNF, N));
    if (pr[1] != 0*pr[1], next);
    nprinc++;
    my(g = pr[2], tt = dlog(nfeltdiv(NF, nfeltmul(NF, -1, g), cj(g))));
    my(sol = matsolvemod(A, w, [tt]~));
    if (type(sol) == "t_INT", next);          \\ FIX (was: sol == 0)
    if (sol == 0, zsol++);
    nanti++;
    my(u0 = nfeltpow(NF, tu, sol[1])); for (i = 1, r, u0 = nfeltmul(NF, u0, nfeltpow(NF, fu[i], sol[i+1])));
    my(lam = nfeltmul(NF, g, u0));
    if (cj(lam) != nfeltmul(NF, -1, lam), error("antisymmetry failed"));
    my(tgt = vector(#Phi, j, imag(emb(lam, Phi[j])) < 0)~);
    my(x = matsolvemod(Ms, 2, tgt));
    if (type(x) == "t_INT", next);            \\ FIX (was: x == 0)
    if (x == 0, zx++);
    my(u = nfalgtobasis(NF, 1)); for (k = 1, #kunits, if (x[k] % 2, u = nfeltmul(NF, u, kunits[k])));
    lam = nfeltmul(NF, lam, u);
    for (j = 1, #Phi, if (imag(emb(lam, Phi[j])) <= 0, error("positivity failed")));
    listput(out, [S, lam]));
  [nclass, nprinc, nanti, #out, out, zsol, zx];
}
cmpdec(tag, J, Phi) = {
  my(a = decide(J, Phi), b = decide_fixed(J, Phi));
  printf("%s: ORIGINAL classes=%d princ=%d antisym=%d ppav=%d | FIXED antisym=%d ppav=%d (zero-vector solutions wrongly discarded by original: antisym step %d, sign step %d)\n",
    tag, a[1], a[2], a[3], a[4], b[3], b[4], b[6], b[7]);
  \\ verify every FIXED solution by an exact Gram determinant on the O_K-lattice S (only valid when J = D^-1)
  b;
}
{
  initell(11);
  my(Pp = cmtype(1), J0 = idealinv(NF, NF.diff));
  my(b = cmpdec("ell=11 End=O_K", J0, Pp));
  for (k = 1, #b[5], my(B = idealhnf(NF, b[5][k][1]), G = gram(B, b[5][k][2])); printf("   class %d: det Gram on S = %s, integral %d\n", k, matdet(G), denominator(G) == 1));
}
{
  initell(17);
  my(Pp = cmtype(1), J0 = idealinv(NF, NF.diff), p = 271, L = idealprimedec(NF, p));
  my(mp = vector(#L, i, minpoly(resid(L[i], PI))));
  my(v0 = select(i -> poldegree(mp[i]) == 1, [1..#L]), pairs = List());
  for (i = 1, #L, for (j = i+1, #L, if (poldegree(mp[i]) == 2 && mp[i] == mp[j], listput(pairs, [i, j]))));
  my(Q0 = idealhnf(NF, L[v0[1]]), P1 = idealmul(NF, L[pairs[1][1]], L[pairs[1][2]]), P2 = idealmul(NF, L[pairs[2][1]], L[pairs[2][2]]));
  cmpdec("ell=17 R: v0 only", idealmul(NF, J0, idealinv(NF, Q0)), Pp);
  cmpdec("ell=17 R: one glued pair", idealmul(NF, J0, idealinv(NF, P1)), Pp);
  cmpdec("ell=17 R: v0 + one pair", idealmul(NF, J0, idealinv(NF, idealmul(NF, Q0, P1))), Pp);
  cmpdec("ell=17 R: two pairs", idealmul(NF, J0, idealinv(NF, idealmul(NF, P1, P2))), Pp);
  cmpdec("ell=17 R = O_K", J0, Pp);
}
quit;
