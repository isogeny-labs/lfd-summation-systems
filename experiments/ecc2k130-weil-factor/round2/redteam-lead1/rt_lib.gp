\\ rt_lib.gp -- red-team (slug redteam-lead1) INDEPENDENT re-implementation for reviewing lead 1
\\ (experiments/ecc2k130-weil-factor/lead1_polarization).  Nothing here is copied from decide_ppav.gp.
\\
\\ Independence from lead 1's code:
\\  * K is built as a relative quadratic extension K = K+(sqrt-7) of K+ (rnfinit), with K+ defined from
\\    the minimal polynomial of th = (z+1/z) + sqrt-7 (z-1/z) obtained by a resultant (no polcompositum);
\\  * zeta is found as a root of Phi_ell in K (nfroots), alpha = (-1+sqrt-7)/2;
\\  * Method A: for every class S of Cl(K), solve lambda S Sbar = J with lambda antisymmetric, and decide
\\    positivity with the units and real signatures of K+ (bnfinit of K+), not of K;
\\  * Method B: NO unit or generator search in K at all. Compute the narrow class group Cl+(K+) (bnrinit),
\\    the image of the norm map Cl(K) -> Cl+(K+), and test whether the class of b+ (b+ O_K = sqrt-7 D_K c)
\\    lies in it.  Existence <=> [b+] in N(Cl(K)) (see the report for the reduction).
\\  * Every "exists" answer is certified by an explicit lattice T and an exact integer Gram determinant.
\\  * Hilbert symbols (r,-7)_v over K+ are computed with nfhilbert at every place that can be nontrivial.
\\ All class groups / units are under GRH (bnfinit) unless stated otherwise.

default(parisizemax, 4000000000);
default(realprecision, 120);

\\ Lucas number from alpha^n - beta^n in Q(sqrt-7) (closed form, not a recurrence)
rt_U(n) = { my(al = Mod('w, 'w^2 + 'w + 2), be = -1 - al, q = lift((al^n - be^n)/(al - be)));
  if (poldegree(q) > 0, error("U not rational")); polcoef(q, 0); }

rt_cj(x) = nfgaloisapply(NFK, SIG, x);
rt_cjI(A) = idealhnf(NFK, nfgaloisapply(NFK, SIG, A));
rt_up(xp) = nfalgtobasis(NFK, rnfeltup(RNF, xp));
rt_down(x) = { my(d = rnfeltdown(RNF, nfbasistoalg(NFK, x))); d; }
rt_sgnvec(xp) = { my(s = nfeltsign(NFKP, xp)); vector(#s, i, s[i] < 0)~; }

rt_init(ell, wantbnfK = 1) = {
  my(t0 = getabstime());
  ELL = ell; BNR = 0;
  my(Pz = polcyclo(ell, 'z), zz = Mod('z, Pz), a = lift(zz + 1/zz), b = lift(zz - 1/zz));
  my(Rs = polresultant(Pz, ('y - a)^2 + 7*b^2, 'z), F = factor(Rs));
  if (#F[,1] != 1 || poldegree(F[1,1]) != ell - 1 || F[1,2] != 2, error("theta+ not primitive: ", F));
  my(Pp = polredbest(F[1,1]));
  NFKP = nfinit(Pp); BNFKP = bnfinit(NFKP, 1);
  RNF = rnfinit(NFKP, 'x^2 + 7, 1);
  NFK = nfinit(RNF);
  if (wantbnfK, BNFK = bnfinit(NFK, 1));
  S7 = nfalgtobasis(NFK, rnfeltreltoabs(RNF, 'x));
  my(YK = rt_up('y), G = nfgaloisconj(NFK));
  SIG = 0;
  for (i = 1, #G, if (nfgaloisapply(NFK, G[i], YK) == YK && nfgaloisapply(NFK, G[i], S7) == -S7, SIG = G[i]));
  if (type(SIG) == "t_INT", error("no complex conjugation"));
  my(TT = varhigher("TT"), R = nfroots(NFK, polcyclo(ell, TT)));
  if (#R != ell - 1, error("zeta not in K"));
  ZE = nfalgtobasis(NFK, R[1]);
  AL = nfalgtobasis(NFK, (-1 + nfbasistoalg(NFK, S7))/2);
  PI = nfeltmul(NFK, AL, ZE); PIB = rt_cj(PI);
  if (nfeltmul(NFK, PI, PIB) != nfalgtobasis(NFK, 2), error("pi*pibar != 2"));
  if (rt_cj(ZE) != nfeltpow(NFK, ZE, -1), error("conj(zeta) != zeta^-1"));
  if (nfeltmul(NFK, AL, AL) + AL + nfalgtobasis(NFK, 2) != 0*AL, error("alpha"));
  EPS = kronecker(-7, ell);
  NDEG = poldegree(NFK.pol);
  printf("[init] ell=%d  deg K=%d  deg K+=%d  h(K+)=%d (cyc %s)  h(K)=%s (cyc %s)  w(K)=%s  #fu(K+)=%d  eps=(-7/ell)=%d  [%d ms]\n",
    ell, NDEG, poldegree(NFKP.pol), BNFKP.no, BNFKP.cyc, if (wantbnfK, BNFK.no, "-"), if (wantbnfK, BNFK.cyc, "-"),
    if (wantbnfK, BNFK.tu[1], "-"), #BNFKP.fu, EPS, getabstime() - t0);
}

\\ ---------------- structural checks ----------------
rt_structure() = {
  my(rd = rnfdisc(NFKP, 'x^2 + 7));
  printf("  K/K+ relative discriminant ideal trivial (unramified at all finite places): %d\n", idealhnf(NFKP, rd[1]) == matid(poldegree(NFKP.pol)));
  my(Pl = idealprimedec(NFKP, ELL));
  printf("  places of K+ above ell: %d, e = %s, f = %s ; v_ell(D_K+) = %s (ell-2 = %d, odd)\n", #Pl,
    vector(#Pl, i, Pl[i].e), vector(#Pl, i, Pl[i].f), vector(#Pl, i, idealval(NFKP, NFKP.diff, Pl[i])), ELL - 2);
  my(up = idealfactor(NFK, rnfidealup(RNF, Pl[1], 1)));
  printf("  place above ell is %s in K/K+ (eps = %d)\n", if (#up[,1] == 1 && up[1,2] == 1, "INERT", if (#up[,1] == 2, "SPLIT", "RAMIFIED?")), EPS);
  my(P7 = idealprimedec(NFKP, 7));
  printf("  places of K+ above 7: %d, e = %s, v(D_K+) = %s, inert in K: %s\n", #P7, Set(vector(#P7, i, P7[i].e)),
    Set(vector(#P7, i, idealval(NFKP, NFKP.diff, P7[i]))), Set(vector(#P7, i, my(u = idealfactor(NFK, rnfidealup(RNF, P7[i], 1))); #u[,1] == 1)));
  my(P2 = idealprimedec(NFKP, 2));
  printf("  places of K+ above 2: %d, all split in K: %d\n", #P2, prod(i = 1, #P2, #idealfactor(NFK, rnfidealup(RNF, P2[i], 1))[,1] == 2));
}

\\ CM type of the canonical lift: Phi = { phi0 o s : s(pi) in P } for a prime P of K above 2.
\\ Claim under test: this set equals { s : s(sqrt-7) = s0 sqrt-7 } for one fixed s0.
rt_cmtype() = {
  my(P2 = idealprimedec(NFK, 2), G = nfgaloisconj(NFK), allok = 1);
  for (k = 1, #P2, my(cnt = matrix(2, 2));
    for (i = 1, #G, my(s = G[i], v = (idealval(NFK, nfgaloisapply(NFK, s, PI), P2[k]) > 0), sg = (nfgaloisapply(NFK, s, S7) == S7));
      cnt[1 + v, 1 + sg]++);
    my(ok = (cnt[1,1] == 0 && cnt[2,2] == 0) || (cnt[1,2] == 0 && cnt[2,1] == 0));
    allok = allok && ok;
    if (k <= 2, printf("  CM type via prime %d above 2: rows v(s(pi))>0 no/yes, cols s(sqrt-7)=-/+ : %s  induced from Q(sqrt-7): %d\n", k, cnt, ok)));
  printf("  CM type induced from Q(sqrt-7) for all %d primes above 2: %d\n", #P2, allok);
  allok;
}

\\ ---------------- Gram / lattice utilities ----------------
rt_gram(B, lam) = matrix(#B, #B, i, j, nfelttrace(NFK, nfeltmul(NFK, lam, nfeltmul(NFK, B[,i], rt_cj(B[,j])))));
rt_inlat(B, x) = denominator(matsolve(B, x)) == 1;
rt_stable(B, x) = { for (i = 1, #B, if (!rt_inlat(B, nfeltmul(NFK, x, B[,i])), return(0))); 1; }
rt_Zpp() = { if (type(ZPP) == "t_MAT" && #ZPP == NDEG && ZPPELL == ELL, return(ZPP));
  my(c = List()); for (i = 0, NDEG - 1, listput(c, nfeltpow(NFK, PI, i))); for (j = 1, NDEG - 1, listput(c, nfeltpow(NFK, PIB, j)));
  ZPPCOLS = Vec(c); ZPPELL = ELL; ZPP = mathnf(matconcat(ZPPCOLS)); ZPP; }

\\ ---------------- primes above p | U: v0, glued pairs, conj partners ----------------
rt_resval(pr, x) = { my(m = nfmodprinit(NFK, pr)); nfmodpr(NFK, x, m); }
rt_special(p) = {
  my(L = idealprimedec(NFK, p), n = #L, mp = vector(n), out);
  for (i = 1, n, mp[i] = minpoly(rt_resval(L[i], PI)));
  my(v0 = [], pairs = List(), cpart = vector(n));
  for (i = 1, n, my(ci = rt_cjI(L[i])); for (j = 1, n, if (idealhnf(NFK, L[j]) == ci, cpart[i] = j)));
  for (i = 1, n, if (L[i].f == 2 && poldegree(mp[i]) == 1, v0 = concat(v0, i)));
  for (i = 1, n, for (j = i + 1, n, if (mp[i] == mp[j] && (L[i].f == 1 || poldegree(mp[i]) == 2), listput(pairs, [i, j]))));
  [L, v0, Vec(pairs), cpart, mp];
}

\\ ---------------- Method A ----------------
rt_dlog(x, tu, w) = { my(y = nfalgtobasis(NFK, 1)); for (t = 0, w - 1, if (y == x, return(t)); y = nfeltmul(NFK, y, tu)); error("rt_dlog: not a root of unity"); }

\\ J: target ideal (lambda S Sbar = J).  Returns [#classes, #principal, #antisym, #W, list of [S, lambda, r]]
rt_A(J) = {
  my(cyc = BNFK.cyc, gens = BNFK.gen, w = BNFK.tu[1], tu = nfalgtobasis(NFK, BNFK.tu[2]));
  my(fu = apply(u -> nfalgtobasis(NFK, u), BNFK.fu), r = #fu);
  my(Hm = Mat(concat([rt_dlog(nfeltdiv(NFK, rt_cj(tu), tu), tu, w)], vector(r, i, rt_dlog(nfeltdiv(NFK, rt_cj(fu[i]), fu[i]), tu, w)))));
  my(fup = BNFKP.fu, Ms = matconcat(concat([rt_sgnvec(-1)], vector(#fup, i, rt_sgnvec(fup[i])))));
  my(ncl = 0, npr = 0, nan = 0, out = List());
  forvec(e = vector(#cyc, i, [0, cyc[i] - 1]),
    ncl++;
    my(S = idealhnf(NFK, 1)); for (i = 1, #cyc, S = idealmul(NFK, S, idealpow(NFK, gens[i], e[i])));
    my(N = idealmul(NFK, J, idealinv(NFK, idealmul(NFK, S, rt_cjI(S)))));
    my(pr = bnfisprincipal(BNFK, N, 1));
    if (pr[1] != 0*pr[1], next);
    npr++;
    my(g = nfalgtobasis(NFK, pr[2]));
    my(om = nfeltdiv(NFK, nfeltmul(NFK, -1, g), rt_cj(g)));
    my(sol = matsolvemod(Hm, w, [rt_dlog(om, tu, w)]~));
    if (type(sol) == "t_INT", next);
    nan++;
    my(u = nfeltpow(NFK, tu, sol[1])); for (i = 1, r, u = nfeltmul(NFK, u, nfeltpow(NFK, fu[i], sol[i + 1])));
    my(lam = nfeltmul(NFK, g, u));
    if (rt_cj(lam) != -lam, error("antisymmetry"));
    my(rr = rt_down(nfeltdiv(NFK, lam, S7)));
    my(x = matsolvemod(Ms, 2, rt_sgnvec(rr)));
    if (type(x) == "t_INT", next);
    my(adj = Mod((-1)^x[1], NFKP.pol)); for (i = 1, #fup, if (x[i + 1] % 2, adj *= fup[i]));
    rr = rr * adj;
    if (rt_sgnvec(rr) != 0*rt_sgnvec(rr), error("sign fix failed"));
    lam = nfeltmul(NFK, S7, rt_up(rr));
    if (idealhnf(NFK, idealmul(NFK, lam, idealmul(NFK, S, rt_cjI(S)))) != idealhnf(NFK, J), error("ideal equation"));
    listput(out, [S, lam, rr]));
  [ncl, npr, nan, #out, out];
}

\\ ---------------- Method B ----------------
rt_normdown(A) = rnfidealdown(RNF, idealmul(NFK, A, rt_cjI(A)));
rt_B(cond) = {
  if (type(BNR) == "t_INT", BNR = bnrinit(BNFKP, [1, vector(NFKP.r1, i, 1)]));
  my(cyc = BNR.cyc, ng = #BNFK.gen, V = matrix(#cyc, ng));
  for (j = 1, ng, V[,j] = bnrisprincipal(BNR, rt_normdown(BNFK.gen[j]), 0));
  my(idx = if (#cyc, matdet(mathnf(matconcat([V, matdiagonal(cyc)]))), 1));
  my(big = idealmul(NFK, S7, idealmul(NFK, NFK.diff, cond)));
  my(bp = rnfidealdown(RNF, big));
  if (idealhnf(NFK, rnfidealup(RNF, bp, 1)) != idealhnf(NFK, big), error("b+ not extended from K+"));
  my(c = bnrisprincipal(BNR, bp, 0));
  my(inH = if (#cyc, type(matsolvemod(V, cyc~, c)) != "t_INT", 1));
  \\ Artin symbol of b+ in Gal(K/K+) computed from the factorization (parity of inert exponents)
  my(fa = idealfactor(NFKP, bp), art = 1);
  for (k = 1, #fa[,1], my(u = idealfactor(NFK, rnfidealup(RNF, fa[k,1], 1)));
    if (#u[,1] == 1 && u[1,2] == 1 && (fa[k,2] % 2), art = -art));
  my(hplus = vecprod(cyc));
  [idx, inH, art, hplus, BNFK.no / (hplus / idx)];
}

\\ ---------------- Hilbert symbols (r,-7)_v over K+ at all places that can be nontrivial ----------------
rt_hilbert(rr, extraP = []) = {
  my(nm = idealnorm(NFKP, rr), P = Set(concat(concat([2, 7, ELL], extraP), concat(factor(numerator(abs(nm)))[,1]~, factor(denominator(nm))[,1]~))));
  my(neg = List(), prodf = 1);
  for (k = 1, #P, my(p = P[k], L = idealprimedec(NFKP, p));
    for (i = 1, #L, my(h = nfhilbert(NFKP, rr, -7, L[i]));
      if (h == -1, my(u = idealfactor(NFK, rnfidealup(RNF, L[i], 1)));
        listput(neg, [p, i, L[i].f, if (#u[,1] == 1, "inert", "split"), idealval(NFKP, rr, L[i])]));
      prodf *= h));
  my(sg = nfeltsign(NFKP, rr), nreal = #select(t -> t < 0, sg));
  prodf *= (-1)^nreal;
  [prodf, nreal, Vec(neg)];
}

\\ ---------------- explicit Deligne module T = Z[pi,pibar] t + c S, and full verification ----------------
rt_locgen(S, Fc, tries = 100000) = {
  my(H = idealhnf(NFK, S));
  for (k = 1, tries, my(t = H * vectorv(#H, i, random(41) - 20));
    my(ok = (t != 0*t)); for (j = 1, #Fc, if (!ok, break); if (idealval(NFK, t, Fc[j]) != idealval(NFK, S, Fc[j]), ok = 0));
    if (ok, return(t)));
  error("no local generator");
}
rt_T(S, cond, t) = { rt_Zpp(); mathnf(concat(matconcat(vector(#ZPPCOLS, k, nfeltmul(NFK, ZPPCOLS[k], t))), idealhnf(NFK, idealmul(NFK, cond, S)))); }
rt_verify(B, lam, label) = {
  my(G = rt_gram(B, lam), d = matdet(G));
  my(rr = rt_down(nfeltdiv(NFK, lam, S7)));
  my(res = [d, denominator(G) == 1, G + G~ == 0, rt_cj(lam) == -lam, rt_stable(B, PI), rt_stable(B, PIB),
            rt_sgnvec(rr) == 0*rt_sgnvec(rr), rt_stable(B, ZE)]);
  printf("    [%s] det Gram=%s integral=%d alternating=%d lam antisym=%d pi,pibar-stable=%d,%d r>>0=%d zeta in End=%d\n",
    label, res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8]);
  res;
}
\\ search for a ppav lattice with End = Z[pi,pibar] + cond: twist the local generator until unimodular
rt_build(S, lam, cond, maxtries = 5000) = {
  my(Fc = if (cond == 1 || idealhnf(NFK, cond) == matid(NDEG), [], idealfactor(NFK, cond)[,1]));
  for (k = 1, maxtries, my(t = rt_locgen(S, Fc), B = rt_T(S, cond, t), G = rt_gram(B, lam));
    if (denominator(G) == 1 && abs(matdet(G)) == 1, return([B, k])));
  [0, maxtries];
}
