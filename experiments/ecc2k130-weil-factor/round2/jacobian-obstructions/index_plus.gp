\\ Real Weil polynomial h of A_ell, disc(h), and the plus/minus split of [O_K : Z[pi,pibar]].
\\ Independent of lead 2's HNF computation: uses discriminants only.
default(parisizemax, 4000000000);
realweil(P, g, q) = {
  my(R = P, h = 0);
  forstep(k = g, 0, -1, my(c = polcoeff(R, g+k, 'x)); h += c*'y^k; R -= c*'x^(g-k)*('x^2+q)^k);
  if(R != 0, error("not a q-Weil polynomial"));
  h;
}
PA(ell) = {
  my(f = polresultant(polcyclo(ell,'z), 'x^2 + 'x*'z + 2*'z^2, 'z)); f / pollead(f);
}
Uell(ell) = {my(w = quadgen(-7), al = w - 1, be = -w); (al^ell - be^ell)/(al - be)};
run(ell) = {
  my(P = PA(ell), g = ell-1, h, dh, dKp, dK, ip2, NT, iK2, u = Uell(ell), ip, iK);
  h = realweil(P, g, 2);
  dh = poldisc(h);
  dKp = ell^(ell-2) * 7^((ell-1)/2);
  dK = ell^(2*(ell-2)) * 7^(ell-1);
  ip2 = abs(dh) / dKp;
  NT = abs(polresultant(h, 'y^2 - 8));
  iK2 = dh^2 * NT / dK;
  printf("ell=%d  U=%s  log2|disc h|=%.3f  log2 N(t^2-8)=%.3f\n", ell, u, log(abs(dh))/log(2), log(NT)/log(2));
  printf("   [O_K+:Z[t]]^2 integral %d square %d ; [O_K:Z[pi,pibar]]^2 integral %d square %d\n", denominator(ip2)==1, issquare(ip2), denominator(iK2)==1, issquare(iK2));
  ip = sqrtint(ip2); iK = sqrtint(iK2);
  printf("   [O_K:R0] == |U|^(ell-2) : %d\n", iK == abs(u)^(ell-2));
  printf("   [O_K+:R0+] factored = %s\n", if(ip == 1, "1", Str(factor(ip))));
  printf("   log2[O_K:R0]=%.3f  log2[O_K+:R0+]=%.3f  log2 minus-index=%.3f\n", log(iK)/log(2), if(ip==1,0,log(ip)/log(2)), log(iK/ip)/log(2));
  printf("   log2 sqrt(|disc h|*N(t^2-8)) = %.3f ; (g(g+1)/4) = %.2f\n", (log(abs(dh)) + log(NT))/2/log(2), g*(g+1)/4);
  [P, h];
}
{
foreach([3,5,11,13,17,19], ell, run(ell));
my(r = run(131), ref = readvec("PA131.txt"));
printf("P_A(131) from resultant == P_A.json : %d\n", Vecrev(r[1]) == ref);
}
