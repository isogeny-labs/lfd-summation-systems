\r kplus_setup.gp
{
for (i = 1, 4, my(ell = [11,17,29,37][i], t0 = getabstime());
  my(S = setupKp(ell), P = S[1]);
  my(bnf = bnfinit(nfinit([P, [2,3,5,7,ell]]), 1));
  my(nar = bnfnarrow(bnf));
  my(sgn = bnfsignunit(bnf));
  printf("ell=%d deg=%d bnfinit %d ms h=%d cyc=%s narrow=%s  units %d  sign-matrix rank F2 (incl -1) %d\n", ell, poldegree(P), getabstime()-t0,
     bnf.no, bnf.cyc, nar[1..2], #bnf.fu, matrank(Mod(concat(vectorv(poldegree(P),j,1), lift(sgn)), 2)));
  writebin(Str("bnf_kplus_", ell, ".bin"), bnf));
}
