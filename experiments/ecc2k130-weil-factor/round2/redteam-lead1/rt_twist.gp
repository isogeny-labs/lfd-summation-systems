\\ rt_twist.gp -- red-team (redteam-lead1): Theorem B at ell = 17 beyond the cases lead 1 ran.
\\ (a) Method B (narrow class group, rt_lib) on ALL 2^8 = 256 orders R = Z[pi,pibar] + c between Z[pi,pibar] and O_K
\\     (c = any product of {v0, the 7 glued place-pairs above 271}); compare with the rule eps*(-1)^n = +1.
\\ (b) Method A (class-by-class search, no CFT) on a sample of these orders.
\\ (c) Explicit unimodular lattices for orders with SEVERAL glued pairs: the local generator is twisted one
\\     component at a time; the local condition at a glued pair c1 is isolated by the mixed Gram b(M, T) with
\\     M = T  cap  (c / c1) S, which is automatically integral at every other component.
\\ usage: gp -q < rt_twist.gp
\r rt_lib.gp

rt_hnfq(M) = { my(d = denominator(M)); mathnf(d*M)/d; }
rt_zdual(B) = (B^(-1))~;
rt_intersect2(A, B) = rt_hnfq(rt_zdual(rt_hnfq(matconcat([rt_zdual(A), rt_zdual(B)]))));
rt_mixgram(M, T, lam) = matrix(#M, #T, i, j, nfelttrace(NFK, nfeltmul(NFK, lam, nfeltmul(NFK, M[,i], rt_cj(T[,j])))));

\\ comps: list of components, each a vector of prime ideals (v0: [Q0]; glued pair: [Pa, Pb])
rt_build_twist(S, lam, comps, maxtries = 4000) = {
  my(allP = concat(comps), cond = idealhnf(NFK, 1));
  for (i = 1, #allP, cond = idealmul(NFK, cond, allP[i]));
  my(t = rt_locgen(S, allP), totaltries = 0);
  for (k = 1, #comps, if (#comps[k] == 1, next);           \\ v0: no twist needed (local unit condition automatic)
    my(c1 = idealmul(NFK, comps[k][1], comps[k][2]), rest = idealmul(NFK, cond, idealinv(NFK, c1)), ok = 0);
    my(fx = matconcat([Col(allP), vectorv(#allP, i, 1)]));
    for (tr = 1, maxtries, totaltries++;
      my(T = rt_T(S, cond, t), M = rt_intersect2(T, idealhnf(NFK, idealmul(NFK, rest, S))));
      if (denominator(rt_mixgram(M, T, lam)) == 1, ok = 1; break);
      \\ new twist: eta = eps mod comps[k][1], 1 mod every other prime of c
      my(eps = vectorv(NDEG, i, random(1000)));
      while (idealval(NFK, eps, comps[k][1]) > 0, eps = vectorv(NDEG, i, random(1000)));
      my(yv = vector(#allP, i, if (idealhnf(NFK, allP[i]) == idealhnf(NFK, comps[k][1]), eps, nfalgtobasis(NFK, 1))));
      my(eta = idealchinese(NFK, fx, yv));
      t = nfeltmul(NFK, t, eta));
    if (!ok, return([0, totaltries])));
  [rt_T(S, cond, t), totaltries];
}

{
  my(t0 = getabstime());
  rt_init(17);
  my(p = 271, sp = rt_special(p), L = sp[1], v0 = sp[2], pairs = sp[3]);
  printf("  p=%d: v0=%s, %d glued pairs\n", p, v0, #pairs);
  my(Q0 = L[v0[1]], PP = vector(#pairs, k, [L[pairs[k][1]], L[pairs[k][2]]]));
  \\ (a) all 256 orders by method B
  my(nagree = 0, ntot = 0, tab = matrix(2, 2));
  forvec(e = vector(1 + #PP, i, [0, 1]),
    my(cond = idealhnf(NFK, 1), n = e[1]);
    if (e[1], cond = idealmul(NFK, cond, Q0));
    for (k = 1, #PP, if (e[k + 1], cond = idealmul(NFK, cond, idealmul(NFK, PP[k][1], PP[k][2]))));
    my(b = rt_B(cond), pred = (EPS * (-1)^n == 1));
    ntot++; nagree += (b[2] == pred); tab[1 + n, 1 + b[2]]++);
  printf("  (a) method B over all %d orders: agree with rule %d/%d ; table rows n(v0)=0/1, cols ppav no/yes: %s [%d s]\n", ntot, nagree, ntot, tab, (getabstime() - t0)\1000);
  \\ (b)+(c) method A and explicit lattices on orders with several pairs
  my(cases = [[[Q0], PP[1], PP[2]], [PP[1], PP[2]], [[Q0], PP[1], PP[2], PP[3]], [PP[1], PP[2], PP[3]]]);
  for (c = 1, #cases, my(comps = cases[c], allP = concat(comps), cond = idealhnf(NFK, 1), n = (#comps[1] == 1));
    for (i = 1, #allP, cond = idealmul(NFK, cond, allP[i]));
    my(J = idealmul(NFK, idealinv(NFK, NFK.diff), idealinv(NFK, cond)), a = rt_A(J));
    printf("  (b) order with v0:%d and %d glued pairs: method A W=%d of %d classes (rule predicts ppav: %d)\n", n, #comps - n, a[4], a[1], EPS * (-1)^n == 1);
    if (a[4] > 0,
      my(S = a[5][1][1], lam = a[5][1][2], rr = a[5][1][3], bt = rt_build_twist(S, lam, comps));
      if (type(bt[1]) == "t_INT", printf("      !! twist search failed after %d tries\n", bt[2]),
        printf("      (c) explicit lattice after %d twist tries\n", bt[2]);
        rt_verify(bt[1], lam, Str("v0:", n, " pairs:", #comps - n));
        my(h = rt_hilbert(rr)); printf("      Hilbert symbols: product=%d, places with -1: %s\n", h[1], h[3]))));
  printf("  total %d s\n", (getabstime() - t0)\1000);
}
