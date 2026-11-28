\r ../../lead1_polarization/decide_ppav.gp
initell(11);
J0 = idealinv(NF, NF.diff);
w = BNF.tu[1]; tu = BNF.tu[2]; fu = BNF.fu; r = #fu;
hexp = vector(r, i, dlog(nfeltdiv(NF, cj(fu[i]), fu[i])));
A = Mat(concat([dlog(nfeltdiv(NF, cj(tu), tu))], hexp));
print("A = ", A, "  image of h mod w: ", matimagemod(A, w));
{
for (e = 0, 4,
  my(S = idealpow(NF, BNF.gen[1], e));
  my(N = idealmul(NF, J0, idealinv(NF, idealmul(NF, S, cjI(S)))));
  my(pr = bnfisprincipal(BNF, N), g = pr[2]);
  my(tt = dlog(nfeltdiv(NF, nfeltmul(NF, -1, g), cj(g))));
  my(sol = matsolvemod(A, w, [tt]~));
  printf("e=%d class=%s tt=%d sol=%s  ; N(S) principal in K+ check: S*Sbar == cj? %d\n", e, pr[1], tt, sol, idealhnf(NF,cjI(idealmul(NF,S,cjI(S))))==idealhnf(NF,idealmul(NF,S,cjI(S)))));
}
