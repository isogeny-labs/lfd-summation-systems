\\ rt_run.gp -- driver for rt_lib.gp.  usage:  echo 'ELLS=[3,5,11,13]; \r rt_run.gp' | gp -q
\r rt_lib.gp

\\ conductor ideal of K for a list of prime ideals (product, each once)
rt_prod(L) = { my(c = idealhnf(NFK, 1)); for (i = 1, #L, c = idealmul(NFK, c, L[i])); c; }

rt_case(tag, cond, nv0, pred) = {
  my(t0 = getabstime(), J = idealmul(NFK, idealinv(NFK, NFK.diff), idealinv(NFK, cond)));
  my(a = rt_A(J), b = rt_B(cond));
  my(existsA = (a[4] > 0), existsB = b[2]);
  my(cons = (existsA == existsB) && (b[2] == (b[3] == 1)) && (a[4] == 0 || a[4] == b[5]));
  printf("  %-34s n(v0)=%d pred=%-3s | A: classes=%d princ=%d antisym=%d W=%d | B: [Cl+:N Cl]=%d inN=%d Artin=%+d |kerN|=%d | consistent=%d agree-with-rule=%d [%d ms]\n",
    tag, nv0, if (pred, "yes", "no"), a[1], a[2], a[3], a[4], b[1], b[2], b[3], b[5], cons, existsA == pred, getabstime() - t0);
  if (a[4] > 0,
    my(S = a[5][1][1], lam = a[5][1][2], rr = a[5][1][3], bt = rt_build(S, lam, cond));
    if (type(bt[1]) == "t_INT", printf("    !! no unimodular lattice found after %d twists\n", bt[2]),
      printf("    explicit T found after %d local-generator tries\n", bt[2]);
      my(v = rt_verify(bt[1], lam, tag));
      \\ negative control: the same lambda on the O_K-lattice S has |det| = N(cond) and is not integral if cond != 1
      if (idealhnf(NFK, cond) != matid(NDEG), my(G0 = rt_gram(idealhnf(NFK, S), lam)); printf("    neg. control (same lambda on S): det=%s integral=%d\n", matdet(G0), denominator(G0) == 1));
      my(h = rt_hilbert(rr));
      printf("    Hilbert symbols (r,-7)_v of this ppav: product=%d, #real places with r<0=%d, places with symbol -1: %s\n", h[1], h[2], h[3])));
  [a[4], b[2], b[3], pred, cons];
}

rt_run(ell) = {
  my(T0 = getabstime());
  rt_init(ell);
  rt_structure();
  rt_cmtype();
  my(U = rt_U(ell), FU = if (abs(U) == 1, matrix(0, 2), factor(abs(U))));
  printf("  U_ell = %d = %s ; all prime factors proven prime: %d ; squarefree: %d\n", U, FU, vecprod(apply(isprime, FU[,1]~)), issquarefree(abs(U)));
  \\ natural (restricted product) polarization lambda_nat = (2 - z - 1/z)/(ell sqrt-7) on T = O_K
  my(zi = nfeltpow(NFK, ZE, -1), m = nfalgtobasis(NFK, 2) - ZE - zi, lnat = nfeltdiv(NFK, m, nfeltmul(NFK, ELL, S7)));
  my(Gn = rt_gram(matid(NDEG), lnat), rn = rt_down(nfeltdiv(NFK, lnat, S7)), sn = nfeltsign(NFKP, rn));
  my(hn = rt_hilbert(rn));
  printf("  natural polarization on O_K: det=%s (ell^2=%d) integral=%d alternating=%d ; r_nat signs: #pos=%d #neg=%d (totally definite) ; Hilbert product=%d, places with -1: %s\n",
    matdet(Gn), ell^2, denominator(Gn) == 1, Gn + Gn~ == 0, #select(t -> t > 0, sn), #select(t -> t < 0, sn), hn[1], hn[3]);
  \\ sanity control of rt_hilbert: product formula on random elements of K+ (small coefficients so the norm is easy to factor)
  my(okprod = 1); for (k = 1, 3, my(x = Mod(random(41) - 20 + (random(5) - 2)*'y + (random(3) - 1)*'y^2, NFKP.pol)); if (x != 0, okprod = okprod && (rt_hilbert(x)[1] == 1)));
  printf("  control: Hilbert product formula holds on 3 random elements of K+: %d\n", okprod);
  my(res = List());
  listput(res, rt_case("End = O_K", idealhnf(NFK, 1), 0, EPS == 1));
  for (k = 1, #FU[,1], my(p = FU[k,1], sp = rt_special(p), L = sp[1], v0 = sp[2], pairs = sp[3], cpart = sp[4]);
    my(typ = if (p % ell == ell - 1 && kronecker(-7, p) == -1, "inert", if (p % ell == 1 && kronecker(-7, p) == 1, "split", "OTHER")));
    printf("  p=%d (%s type): %d primes of K, f=%s, v0=%s, glued pairs=%d\n", p, typ, #L, Set(vector(#L, i, L[i].f)), v0, #pairs);
    if (typ == "inert",
      if (#v0 != 1, error("expected exactly one self-glued place"));
      my(Q0 = L[v0[1]], P1 = rt_prod([L[pairs[1][1]], L[pairs[1][2]]]));
      \\ inert type: every place conj-stable; check
      printf("    every prime above p conj-stable: %d\n", prod(i = 1, #L, cpart[i] == i));
      listput(res, rt_case(Str("R: v0 above ", p), idealhnf(NFK, Q0), 1, EPS == -1));
      listput(res, rt_case(Str("R: one glued pair above ", p), P1, 0, EPS == 1));
      listput(res, rt_case(Str("R: v0 + one pair above ", p), idealmul(NFK, Q0, P1), 1, EPS == -1));
      if (#pairs >= 2, listput(res, rt_case(Str("R: two pairs above ", p), idealmul(NFK, P1, rt_prod([L[pairs[2][1]], L[pairs[2][2]]])), 0, EPS == 1))));
    if (typ == "split",
      my(sc = select(pq -> cpart[pq[1]] == pq[2], pairs), nsc = select(pq -> cpart[pq[1]] != pq[2], pairs));
      printf("    self-conjugate glued pairs: %d ; others: %d\n", #sc, #nsc);
      if (#sc, listput(res, rt_case(Str("R: self-conj split pair at ", p), rt_prod([L[sc[1][1]], L[sc[1][2]]]), 0, EPS == 1)));
      if (#nsc, my(q = nsc[1], cq = [cpart[q[1]], cpart[q[2]]]);
        listput(res, rt_case(Str("R: split pair + conj pair at ", p), rt_prod([L[q[1]], L[q[2]], L[cq[1]], L[cq[2]]]), 0, EPS == 1))));
    if (typ == "OTHER", error("prime of unexpected type")));
  printf("  SUMMARY ell=%d: %d cases, all methods consistent: %d, rule eps*(-1)^n=+1 <=> ppav agrees in all: %d [%d s]\n\n", ell, #res,
    prod(i = 1, #res, res[i][5]), prod(i = 1, #res, (res[i][1] > 0) == res[i][4]), (getabstime() - T0)\1000);
  Vec(res);
}

{ for (i = 1, #ELLS, rt_run(ELLS[i])); }
