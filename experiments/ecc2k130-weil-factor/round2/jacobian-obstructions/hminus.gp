\\ Relative class number h^- of K_ell = Q(zeta_ell, sqrt(-7)) by the analytic formula
\\ h^- = Q * w * prod_{chi odd} (-B_{1,chi}/2)  (Washington, Thm 4.17), Q in {1,2}; independent of bnfinit.
\\ Controls: (a) Q(zeta_p) vs the known table; (b) K_ell, ell<=19, vs bnfinit h_K and h_{K+} (GRH-conditional).
default(parisizemax, 4000000000);
default(realprecision, 400);
B1(ell, i, e) = {
  my(n = ell-1, g = znprimroot(ell), f = if(i % n == 0, 1, ell) * if(e, 7, 1), s = 0, z = exp(2*Pi*I/n), tab);
  tab = vector(ell-1); for(k = 0, n-1, tab[lift(g^k)] = z^((i*k) % n));
  for(a = 1, f-1,
    if(gcd(a, f) > 1, next);
    my(v = if(i % n, tab[a % ell], 1));
    if(e, v *= kronecker(-7, a));
    s += v*a);
  s / f;
}
logprodK(ell) = {
  my(n = ell-1, L = 0., P = 1.);
  for(i = 0, n-1, for(e = 0, 1,
     if((-1)^i * (-1)^e == -1, my(b = -B1(ell, i, e)/2); P *= b)));
  P;
}
logprodCyc(p) = {my(P = 1.); forstep(i = 1, p-2, 2, P *= -B1(p, i, 0)/2); P};
{
print("== control (a): w*prod for Q(zeta_p) (Q=1 for cyclotomic); known h^-: 23:3 29:8 31:9 37:37 41:121 43:211 47:695 59:41241 ==");
foreach([23,29,31,37,41,43,47,59], p, my(v = 2*p*logprodCyc(p)); printf("p=%d  %.8g  (imag %.1e)\n", p, real(v), abs(imag(v))));
print("== control (b): K_ell: analytic w*prod vs bnfinit(K).no / bnfinit(K+).no ==");
foreach([3,5,11,13,17,19], ell,
  my(v = 2*ell*logprodK(ell), pK = polredbest(polcompositum(polcyclo(ell), 'x^2+7)[1]), bK, sub, pKp, bKp);
  bK = bnfinit(pK, 1);
  \\ K+ = Q(zeta+zeta^-1, sqrt(7*ell*(-1)^((ell-1)/2) * -1)): even characters; build it directly
  pKp = polredbest(polcompositum(polsubcyclo(ell, (ell-1)/2), 'x^2 - (-7)*(-1)^((ell-1)/2)*ell)[1]);
  if(ell == 3, pKp = polredbest('x^2 - 21));
  bKp = bnfinit(pKp, 1);
  printf("ell=%d  analytic w*prod=%.10g (imag %.1e)   h_K=%d  h_K+=%d  h_K/h_K+=%s  deg K+=%d totally real=%d\n",
     ell, real(v), abs(imag(v)), bK.no, bKp.no, bK.no/bKp.no, poldegree(pKp), polsturm(pKp)==poldegree(pKp)));
my(v = 2*131*logprodK(131));
printf("ell=131 analytic w*prod = %.8e  imag %.1e  log2 = %.4f\n", real(v), abs(imag(v)), log(abs(real(v)))/log(2));
printf("ell=131 frac(w*prod) = %.3e   frac(2*w*prod) = %.3e\n", abs(real(v) - round(real(v))), abs(2*real(v) - round(2*real(v))));
printf("ell=131 w*prod rounded = %d\n", round(real(v)));
my(c = 2*131*logprodCyc(131)); printf("control: h^-(Q(zeta_131)) analytic = %d (frac %.1e)\n", round(real(c)), abs(real(c)-round(real(c))));
}
