\\ Pic of the index-p orders O = {x in O_K : chi1(x) = chi2(x) mod p} (lead 3) at small ell.
\\ K = Q(zeta_ell, sqrt(-7)), pi = alpha*zeta, pibar = beta/zeta.  GRH (bnfinit) for h_K and the units.
\\ Exact sequence: 1 -> O^* -> O_K^* -> (O_K/f)^*/(O/f)^* = F_p^* -> Pic(O) -> Pic(O_K) -> 1,  f = P1*P2,
\\ so #Pic(O) = h_K * (p-1) / #im,  im = {chi1(u)/chi2(u) : u in O_K^*}  (a subgroup of F_p^*).
\\ Since O is Gorenstein (lead 3), lattices with End = O are invertible O-ideals, so by Deligne
\\ #{iso classes of A' over F_2 with End(A') = O} = #Pic(O).
default(parisizemax, 4000000000);
run(ell, p) = {
  my(t0 = getabstime());
  my(C = polcompositum(polcyclo(ell), x^2+x+2, 1)[1]);
  my(g = C[1], z = C[2], a = C[3]);
  my(R = polredbest(g, 1)); g = R[1];
  z = Mod(subst(lift(z), x, R[2]), g); a = Mod(subst(lift(a), x, R[2]), g);
  if (z^ell != 1 || a^2+a+2 != 0 || z == 1, error("embedding"));
  my(bnf = bnfinit(g, 1), nf = bnf.nf);
  my(h = bnf.no, tu = bnf.tu, fu = bnf.fu);
  printf("ell=%d p=%d  deg=%d  h_K=%d (GRH)  cyc=%s  w=%d  #fu=%d  reg=%.6g  time=%ds\n",
         ell, p, poldegree(g), h, bnf.cyc, tu[1], #fu, bnf.reg, (getabstime()-t0)\1000);
  \\ sanity: units have norm +-1
  for (j = 1, #fu, if (abs(norm(fu[j])) != 1, error("fu not unit")));
  \\ Hasse unit index Q = [E : mu E^+]: Q = 2 iff some unit has u/ubar not in mu^2.
  \\ complex conjugation on K: zeta -> 1/zeta, alpha -> -1-alpha; find it among nfgaloisconj.
  my(G = nfgaloisconj(nf), cc = 0);
  for (i = 1, #G, my(s = Mod(G[i], g)); if (subst(lift(z), x, s) == z^-1 && subst(lift(a), x, s) == -1-a, cc = s));
  if (cc == 0, error("no conj"));
  my(cj = u -> subst(lift(u), x, cc));
  my(zw = tu[2], Qidx = 1);   \\ zw generates mu_K (order w)
  for (j = 1, #fu, my(r = fu[j] / cj(fu[j]));
     my(k = -1); for (e = 0, tu[1]-1, if (zw^e == r, k = e; break));
     if (k < 0, error("u/ubar not a root of unity"));
     if (k % 2 == 1, Qidx = 2));
  printf("   Hasse unit index Q = %d\n", Qidx);
  \\ characters O_K -> F_p (p splits completely)
  my(D = idealprimedec(nf, p));
  if (#D != poldegree(g) || vecmax(apply(P->P.f, D)) != 1, error("p not totally split"));
  my(ch = vector(#D, i, my(m = nfmodprinit(nf, D[i])); [m, nfmodpr(nf, z, m), nfmodpr(nf, a, m)]));
  my(gen = znprimroot(p));
  my(lg = v -> znlog(Mod(simplify(v.pol), p), gen));
  my(units = concat([tu[2]], fu));
  my(res = List(), npairs = 0);
  for (i = 1, #D, for (j = i+1, #D,
    my(wi = ch[i][2], ai = ch[i][3], wj = ch[j][2], aj = ch[j][3]);
    if (ai*wi == aj*wj && (-1-ai)/wi == (-1-aj)/wj,
      npairs++;
      \\ is the pair conjugation-stable?  chi_j = chi_i o conj  <=>  w_j = 1/w_i and a_j = -1-a_i
      my(cs = (wj == 1/wi && aj == -1-ai));
      my(L = vector(#units, k, lg(nfmodpr(nf, units[k], ch[i][1]) / nfmodpr(nf, units[k], ch[j][1]))));
      my(d = gcd(concat(L, [p-1])), im = (p-1)/d);
      \\ independent lower bound: image of mu_K alone (zeta -> w_i/w_j)
      my(imz = znorder(Mod(simplify((wi/wj).pol), p)));
      my(pic = h*(p-1)/im);
      listput(res, [i, j, cs, im, imz, pic]);
      printf("   pair (%d,%d) conj-stable=%d  #im(O_K^*)=%d  #im(mu_K)=%d  -> #Pic(O)=h_K*(p-1)/#im=%d  ker(Pic O->Pic O_K)=%d\n",
             i, j, cs, im, imz, pic, (p-1)/im))));
  printf("   #index-p orders containing R = %d ; total iso classes over all of them = %d ; with End=O_K: %d\n",
         npairs, sum(k=1,#res,res[k][6]), h);
  printf("   predicted for the conj-stable pair (#im = ell*Q, since chi1(u)/chi2(u) = chi1(u/ubar) and {u/ubar} = mu_K^2 or mu_K): #Pic = h_K*(p-1)/(ell*Q) = %d\n", h*(p-1)/(ell*Qidx));
  \\ negative control: a prime q not dividing U_ell gives no merged pair (no index-q order containing R)
  my(Ul = (polcoef(lift(Mod(x,x^2+x+2)^ell),1)));
  printf("   U_ell=%d ; time %ds\n", Ul, (getabstime()-t0)\1000);
  [ell, p, h, Qidx, npairs, Vec(res)]
};
negctrl(ell, q) = {
  my(C = polcompositum(polcyclo(ell), x^2+x+2, 1)[1], g = C[1], z = C[2], a = C[3]);
  my(nf = nfinit(g), D = idealprimedec(nf, q), n = 0);
  my(ch = vector(#D, i, my(m = nfmodprinit(nf, D[i])); [nfmodpr(nf, z, m), nfmodpr(nf, a, m), D[i].f]));
  for (i = 1, #D, for (j = i+1, #D, if (ch[i][3] == 1 && ch[j][3] == 1 && ch[i][1]*ch[i][2] == ch[j][1]*ch[j][2] && (-1-ch[i][2])/ch[i][1] == (-1-ch[j][2])/ch[j][1], n++)));
  printf("negative control ell=%d q=%d: merged pairs = %d (expect 0)\n", ell, q, n);
};

\\ Self-glued inert type (the p2 analogue): p = -1 mod ell, (-7/p) = -1, p | U_ell.
\\ O_P = {x : x mod P in F_p}, P the unique prime of degree 2 with chi(pi) in F_p.
\\ (O_K/P)^*/(O_P/P)^* = F_{p^2}^*/F_p^*, cyclic of order p+1; u -> u^(p-1) identifies it with the norm-1 group.
run_inert(ell, p) = {
  my(t0 = getabstime());
  my(C = polcompositum(polcyclo(ell), x^2+x+2, 1)[1]);
  my(g = C[1], z = C[2], a = C[3]);
  my(R = polredbest(g, 1)); g = R[1];
  z = Mod(subst(lift(z), x, R[2]), g); a = Mod(subst(lift(a), x, R[2]), g);
  my(bnf = bnfinit(g, 1), nf = bnf.nf, h = bnf.no, tu = bnf.tu, fu = bnf.fu);
  my(G = nfgaloisconj(nf), cc = 0);
  for (i = 1, #G, my(s = Mod(G[i], g)); if (subst(lift(z), x, s) == z^-1 && subst(lift(a), x, s) == -1-a, cc = s));
  my(cj = u -> subst(lift(u), x, cc), zw = tu[2], Qidx = 1);
  for (j = 1, #fu, my(r = fu[j]/cj(fu[j]), k = -1); for (e = 0, tu[1]-1, if (zw^e == r, k = e; break)); if (k < 0, error("bad")); if (k%2, Qidx = 2));
  printf("ell=%d p=%d (inert type)  deg=%d h_K=%d (GRH) cyc=%s w=%d Q=%d\n", ell, p, poldegree(g), h, bnf.cyc, tu[1], Qidx);
  my(D = idealprimedec(nf, p), found = 0);
  for (i = 1, #D, my(P = D[i]); if (P.f != 2, next);
    my(m = nfmodprinit(nf, P), w = nfmodpr(nf, z, m), al = nfmodpr(nf, a, m));
    my(lam = al*w, inFp = (lam^p == lam), lamb = (-1-al)/w);
    if (!inFp, next);
    found++;
    my(selfg = (idealval(nf, nfbasistoalg(nf, P.gen[2]) , P) >= 0) && idealval(nf, idealhnf(nf, [p, nfalgtobasis(nf, cj(nfbasistoalg(nf, P.gen[2])))]), P) >= 1);
    my(units = concat([zw], fu), ords = vector(#units, k, my(v = nfmodpr(nf, units[k], m)^(p-1)); fforder(v)));
    my(im = lcm(ords));  \\ cyclic group: subgroup generated has order lcm of orders
    printf("   P#%d: lambda=pi mod P in F_p: %d (lambda^2=2: %d), pibar in F_p: %d, conj(P)=P: %d\n", i, inFp, lam^2==2, lamb^p==lamb, selfg);
    printf("   orders of unit images in F_{p^2}^*/F_p^* (mu_K gen first): %s\n", ords);
    printf("   #im(O_K^*) = %d (predicted ell*Q = %d) -> #Pic(O_P) = h_K*(p+1)/#im = %d ; ker = %d\n", im, ell*Qidx, h*(p+1)/im, (p+1)/im));
  printf("   primes P of degree 2 with pi mod P in F_p: %d (expect 1); time %ds\n", found, (getabstime()-t0)\1000);
};
