\\ h(K+) for ell = 5, 13, 17 via the totally real subfield of degree ell-1 (GRH-conditional bnfinit)
default(parisizemax, 4000000000);
{
foreach([5,13,17], ell,
  my(pK = polredbest(polcompositum(polcyclo(ell), 'x^2+7)[1]), S = nfsubfields(pK, ell-1), R = []);
  foreach(S, s, if(polsturm(s[1]) == ell-1, R = concat(R, [s[1]])));
  if(#R != 1, error("expected one totally real subfield"));
  my(bKp = bnfinit(polredbest(R[1]), 1), bK = bnfinit(pK, 1));
  printf("ell=%d  h_K=%d  h_K+=%d  h_K/h_K+=%s\n", ell, bK.no, bKp.no, bK.no/bKp.no));
}
