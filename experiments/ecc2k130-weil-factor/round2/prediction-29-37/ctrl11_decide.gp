\\ Control: rerun lead-1 decide() at ell=11 and list which ideal classes S of K carry a ppav with End=O_K.
\r ../../lead1_polarization/decide_ppav.gp
initell(11);
Pp = cmtype(1); J0 = idealinv(NF, NF.diff);
r = decide(J0, Pp);
report("ell=11 End=O_K", r);
{
forvec(e = vector(#BNF.cyc, i, [0, BNF.cyc[i]-1]),
  my(S = idealhnf(NF, 1)); for (i = 1, #BNF.cyc, S = idealmul(NF, S, idealpow(NF, BNF.gen[i], e[i])));
  my(ok = 0); for (k = 1, #r[5], if (bnfisprincipal(BNF, idealmul(NF, S, idealinv(NF, r[5][k][1])), 0) == 0*BNF.cyc~, ok = 1));
  printf("  class e=%s : ppav found %d\n", e, ok));
}
