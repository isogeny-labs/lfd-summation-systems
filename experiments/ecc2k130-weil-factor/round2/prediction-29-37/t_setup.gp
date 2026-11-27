\r kplus_setup.gp
{
for (i = 1, 4, my(ell = [11,17,29,37][i], t0 = getabstime());
  my(S = setupKp(ell), P = S[1]);
  my(dg = 7^((ell-1)/2) * ell^(ell-2));
  my(nf = nfinit([P, [2,3,5,7,ell]]));
  printf("ell=%d deg=%d nfinit %d ms  disc ok %d  certify-list=%s\n", ell, poldegree(P), getabstime()-t0,
     abs(nf.disc) == dg, nfcertify(nf)));
}
