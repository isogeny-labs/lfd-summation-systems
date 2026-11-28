\\ rt_gor.gp -- red-team (redteam-lead1): lattice-level checks of the local claims used in Theorems A/B.
\\ For R_c := Z[pi,pibar] + c  (c = conductor ideal: product of primes of K at the non-maximal components):
\\   (G1) R_c is Gorenstein:  R^dual * (R : R^dual) = R   (so every lattice with multiplier ring R is invertible);
\\   (G2) R^dual * O_K = (D_K c)^-1, i.e. the O_K-hull of R^dual is the inverse different times c^-1
\\        (this is the identity behind "necessity": lambda T Tbar = R^dual implies lambda S Sbar = D^-1 c^-1);
\\   (G3) the multiplier ring of R^dual is R, and [O_K : R] = p^(#v0 in c) * p^(2 #glued inert pairs) * p^(#split pairs);
\\   (G4) Z[pi,pibar] + c is a ring (closed under multiplication).
\\ usage: echo 'ELLS=[11,17];' | cat - rt_gor.gp | gp -q
\r rt_lib.gp

rt_hnfq(M) = { my(d = denominator(M)); mathnf(d*M)/d; }
rt_trmat() = { my(n = NDEG, E = matid(n)); matrix(n, n, i, j, nfelttrace(NFK, nfeltmul(NFK, E[,i], E[,j]))); }
rt_dual(B) = rt_hnfq((B~ * TRM)^(-1));
rt_prodlat(B1, B2) = rt_hnfq(matconcat(concat(vector(#B1, i, vector(#B2, j, nfeltmul(NFK, B1[,i], B2[,j]))))));
rt_zdual(B) = (B^(-1))~;                                      \\ dual for the standard coordinate pairing
rt_intersect(Ls) = { my(S = rt_hnfq(matconcat(vector(#Ls, k, rt_zdual(Ls[k]))))); rt_hnfq(rt_zdual(S)); }
rt_colon(L, M) = rt_intersect(vector(#M, i, rt_hnfq(matconcat(vector(#L, j, nfeltdiv(NFK, L[,j], M[,i]))))));
rt_eq(B1, B2) = rt_hnfq(B1) == rt_hnfq(B2);

rt_gorcase(tag, cond, expect_index) = {
  my(R = rt_hnfq(concat(rt_Zpp(), idealhnf(NFK, cond))));
  my(ring = rt_eq(rt_prodlat(R, R), R));
  my(Rd = rt_dual(R), C = rt_colon(R, Rd), gor = rt_eq(rt_prodlat(Rd, C), R));
  my(hull = rt_prodlat(Rd, matid(NDEG)), target = idealinv(NFK, idealmul(NFK, NFK.diff, cond)));
  my(mult = rt_colon(Rd, Rd));
  my(idx = abs(matdet(R)));
  printf("  %-30s ring=%d  Gorenstein (R^v (R:R^v) = R)=%d  O_K-hull of R^v == (D c)^-1: %d  End(R^v) == R: %d  [O_K:R] = %s (expected %s: %d)\n",
    tag, ring, gor, rt_eq(hull, idealhnf(NFK, target)), rt_eq(mult, R), factor(idx), factor(expect_index), idx == expect_index);
  [ring, gor, rt_eq(hull, idealhnf(NFK, target)), rt_eq(mult, R), idx == expect_index];
}

rt_gorrun(ell) = {
  rt_init(ell, 0);
  TRM = rt_trmat();
  my(U = rt_U(ell), FU = factor(abs(U)), out = List());
  for (k = 1, #FU[,1], my(p = FU[k,1], sp = rt_special(p), L = sp[1], v0 = sp[2], pairs = sp[3], cpart = sp[4]);
    if (#v0,
      my(Q0 = idealhnf(NFK, L[v0[1]]), P1 = idealmul(NFK, L[pairs[1][1]], L[pairs[1][2]]), P2 = idealmul(NFK, L[pairs[2][1]], L[pairs[2][2]]));
      listput(out, rt_gorcase(Str("v0 above ", p), Q0, p));
      listput(out, rt_gorcase(Str("one glued pair above ", p), P1, p^2));
      listput(out, rt_gorcase(Str("v0 + two pairs above ", p), idealmul(NFK, Q0, idealmul(NFK, P1, P2)), p^5)),
      my(sc = select(pq -> cpart[pq[1]] == pq[2], pairs), q = select(pq -> cpart[pq[1]] != pq[2], pairs)[1]);
      listput(out, rt_gorcase(Str("self-conj split pair at ", p), idealmul(NFK, L[sc[1][1]], L[sc[1][2]]), p));
      listput(out, rt_gorcase(Str("split pair + conj at ", p), idealmul(NFK, idealmul(NFK, L[q[1]], L[q[2]]), idealmul(NFK, L[cpart[q[1]]], L[cpart[q[2]]])), p^2)));
    \\ the minimal order Z[pi,pibar] itself: conductor = all non-maximal components; index |U|^(ell-2)
    my(Rmin = rt_Zpp(), Rd = rt_dual(Rmin), C = rt_colon(Rmin, Rd));
    printf("  Z[pi,pibar] itself: [O_K:Z[pi,pibar]] = %s (|U|^(ell-2) = %s), Gorenstein: %d\n", factor(abs(matdet(Rmin))), factor(abs(U)^(ell-2)), rt_eq(rt_prodlat(Rd, C), Rmin)));
  printf("  ell=%d: all lattice checks pass: %d\n\n", ell, prod(i = 1, #out, vecprod(out[i])));
}
{ for (i = 1, #ELLS, rt_gorrun(ELLS[i])); }
