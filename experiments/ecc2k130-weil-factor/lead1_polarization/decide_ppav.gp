\\ Lead 1, method 2 (PARI/GP, direct enumeration; class groups/units under GRH unless certified).
\\ K = Q(zeta_ell, sqrt-7), pi = alpha*zeta, pibar = beta/zeta, alpha^2+alpha+2 = 0.
\\ Criterion used (Deligne 1969; Howe 1995, Trans. AMS 347, sect. 4, re-derived in REPORT.md sect. 2):
\\   a Deligne module T (a Z[pi,pibar]-stable lattice in K) carries a principal polarization iff
\\   there is lambda in K^* with conj(lambda) = -lambda, Im phi(lambda) > 0 for phi in Phi, and
\\   lambda * conj(T) = T^dual  (trace dual), i.e. b(x,y) = Tr(lambda x conj(y)) is unimodular on T.
\\ Phi = {phi : phi(sqrt-7) = s0*i*sqrt7}; both s0 = +-1 are run (the canonical-lift type is one of them).
default(parisizemax, 4000000000);
default(realprecision, 250);

setupK(ell) = {
  my(P = polcompositum(polcyclo(ell, 'x), 'x^2 + 'x + 2, 1)[1]);
  my(Q = polredbest(P[1], 1), a = Q[2], pol = Q[1]);
  my(ze = Mod(subst(lift(P[2]), 'x, a), pol), al = Mod(subst(lift(P[3]), 'x, a), pol));
  [pol, ze, al];
}

\\ globals per ell: BNF, NF, SIG (conjugation), ZE, AL, SQ (=sqrt-7), PI, PIB
initell(ell) = {
  my(S = setupK(ell), t0 = getabstime());
  NF = nfinit(S[1]); ZE = nfalgtobasis(NF, S[2]); AL = nfalgtobasis(NF, S[3]);
  SQ = nfalgtobasis(NF, 2*S[3] + 1);
  PI = nfeltmul(NF, AL, ZE); PIB = nfeltdiv(NF, 2, PI);
  BNF = bnfinit(S[1], 1);
  my(G = nfgaloisconj(NF)); SIG = 0;
  for (i = 1, #G, my(s = G[i]);
    if (nfgaloisapply(NF, s, ZE) == nfeltpow(NF, ZE, -1) && nfgaloisapply(NF, s, AL) == nfalgtobasis(NF, -1 - nfbasistoalg(NF, AL)), SIG = s));
  if (SIG == 0, error("no complex conjugation found"));
  printf("ell=%d deg=%d bnfinit %d ms, h=%d cyc=%s, #fu=%d, w=%d\n", ell, poldegree(S[1]), getabstime()-t0, BNF.no, BNF.cyc, #BNF.fu, BNF.tu[1]);
}

cj(x) = nfgaloisapply(NF, SIG, x);
dlog(x) = { my(w = BNF.tu[1], tu = BNF.tu[2]); for (t = 0, w-1, if (nfeltpow(NF, tu, t) == x, return(t))); error("not a root of unity"); }
signrow(x, Phi) = vector(#Phi, j, my(v = emb(x, Phi[j])); if (abs(imag(v)) > 1e-50*abs(v), error("unit not real")); real(v) < 0);
locgen(t, S, F) = { if (t == 0, return(0)); for (k = 1, #F, if (idealval(NF, t, F[k]) != idealval(NF, S, F[k]), return(0))); 1; }
cjI(I) = nfgaloisapply(NF, SIG, I);
emb(x, r) = subst(lift(nfbasistoalg(NF, x)), variable(NF.pol), r);

\\ CM type: list of complex roots r (one per conjugate pair), oriented so s0*Im(emb(sqrt-7)) > 0.
cmtype(s0) = {
  my(R = NF.roots, out = List());
  for (j = NF.r1 + 1, #R, my(r = R[j]); if (s0 * imag(emb(SQ, r)) < 0, r = conj(r)); listput(out, r));
  Vec(out);
}

\\ Decide, for a target ideal J (conj-stable), whether exist ideal S and lambda with lambda S Sbar = J,
\\ conj(lambda) = -lambda, Im phi(lambda) > 0 on Phi. Returns list of [S, lambda] (one per class that works).
decide(J, Phi) = {
  my(cyc = BNF.cyc, gens = BNF.gen, w = BNF.tu[1], tu = BNF.tu[2], fu = BNF.fu, r = #fu);
  my(hexp = vector(r), out = List(), nclass = 0, nprinc = 0, nanti = 0);
  \\ h(u) = conj(u)/u as a power of tu
  for (i = 1, r, hexp[i] = dlog(nfeltdiv(NF, cj(fu[i]), fu[i])));
  my(A = Mat(concat([dlog(nfeltdiv(NF, cj(tu), tu))], hexp)));
  my(ker = matsolvemod(A, w, [0]~, 1)[2]);
  \\ real units from kernel basis; their signs at Phi
  my(kunits = vector(#ker[1,], k, my(c = ker[,k], u = nfeltpow(NF, tu, c[1])); for (i = 1, r, u = nfeltmul(NF, u, nfeltpow(NF, fu[i], c[i+1]))); u));
  my(Ms = matrix(#Phi, #kunits, i, k, 0));
  for (k = 1, #kunits, my(s = signrow(kunits[k], Phi)); for (i = 1, #Phi, Ms[i,k] = s[i]));
  forvec(e = vector(#cyc, i, [0, cyc[i]-1]),
    nclass++;
    my(S = idealhnf(NF, 1)); for (i = 1, #cyc, S = idealmul(NF, S, idealpow(NF, gens[i], e[i])));
    my(N = idealmul(NF, J, idealinv(NF, idealmul(NF, S, cjI(S)))));
    my(pr = bnfisprincipal(BNF, N));
    if (pr[1] != 0*pr[1], next);
    nprinc++;
    my(g = pr[2], tt = dlog(nfeltdiv(NF, nfeltmul(NF, -1, g), cj(g))));
    my(sol = matsolvemod(A, w, [tt]~));
    if (sol == 0, next);
    nanti++;
    my(u0 = nfeltpow(NF, tu, sol[1])); for (i = 1, r, u0 = nfeltmul(NF, u0, nfeltpow(NF, fu[i], sol[i+1])));
    my(lam = nfeltmul(NF, g, u0));
    if (cj(lam) != nfeltmul(NF, -1, lam), error("antisymmetry failed"));
    \\ target: need real unit u with sign(phi(u)) = sign(Im phi(lam)) ; vector of 1 where Im<0
    my(tgt = vector(#Phi, j, imag(emb(lam, Phi[j])) < 0)~);
    my(x = matsolvemod(Ms, 2, tgt));
    if (x == 0, next);
    my(u = nfalgtobasis(NF, 1)); for (k = 1, #kunits, if (x[k] % 2, u = nfeltmul(NF, u, kunits[k])));
    lam = nfeltmul(NF, lam, u);
    for (j = 1, #Phi, if (imag(emb(lam, Phi[j])) <= 0, error("positivity failed")));
    listput(out, [S, lam]));
  [nclass, nprinc, nanti, #out, out];
}

\\ Z-basis (columns, nf-basis coordinates) of a lattice given as HNF; Gram of b(x,y)=Tr(lam x conj y)
gram(B, lam) = matrix(#B, #B, i, j, nfelttrace(NF, nfeltmul(NF, lam, nfeltmul(NF, B[,i], cj(B[,j])))));
inlat(B, x) = denominator(matsolve(B, x)) == 1;
stable(B, x) = { for (i = 1, #B, if (!inlat(B, nfeltmul(NF, x, B[,i])), return(0))); 1; }

\\ Full verification of a Deligne module T (Z-basis B, nf coordinates) with polarization lam
verify(B, lam, Phi, label) = {
  my(G = gram(B, lam), d = matdet(G));
  my(okanti = (G + G~ == 0), okint = (denominator(G) == 1));
  my(okF = stable(B, PI), okV = stable(B, PIB), zst = stable(B, ZE));
  my(okpos = 1); for (j = 1, #Phi, if (imag(emb(lam, Phi[j])) <= 0, okpos = 0));
  my(okl = (cj(lam) == nfeltmul(NF, -1, lam)));
  printf("  [%s] det Gram=%s integral=%d alternating=%d lam antisym=%d F,V-stable=%d,%d Phi-positive=%d zeta in End=%d\n",
    label, d, okint, okanti, okl, okF, okV, okpos, zst);
  [d, okint, okanti, okl, okF, okV, okpos, zst];
}

\\ Primes of K above p, with residue data of pi, pibar
primesabove(p) = idealprimedec(NF, p);
resid(pr, x) = nfmodpr(NF, x, nfmodprinit(NF, pr));

\\ Lattice T = Z t + Q*S, Q a product of primes at which S is locally generated by t
sublattice(S, Q) = {
  my(QS = idealmul(NF, Q, S), H = idealhnf(NF, S), t = 0);
  for (i = 1, #H, if (!inlat(QS, H[,i]), t = H[,i]; break));
  \\ need t to generate S locally at every prime of Q: test t + combinations until valuations match
  my(F = idealfactor(NF, Q)[,1]);
  my(tries = 0); while (!locgen(t, S, F), tries++; t = H * vectorv(#H, i, random(7) - 3); if (tries > 10000, error("no local generator")));
  mathnf(concat(Mat(t), idealhnf(NF, QS)));
}

lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }

\\ Report one decision
report(tag, res) = printf("  %s: classes=%d, principal J*(S Sbar)^-1: %d, with antisymmetric generator: %d, Phi-positive (=> ppav): %d\n", tag, res[1], res[2], res[3], res[4]);

run(ell) = {
  initell(ell);
  my(D = NF.diff, n = poldegree(NF.pol));
  my(dgen = nfeltmul(NF, SQ, nfeltpow(NF, nfeltadd(NF, ZE, nfeltmul(NF, -1, nfeltpow(NF, ZE, -1))), ell - 2)));
  printf("  different check: N(diff)=|disc| %d ; diff = (sqrt-7 (z-1/z)^(ell-2)) %d ; generator real (conj=self) %d\n",
    idealnorm(NF, D) == abs(NF.disc), idealhnf(NF, dgen) == idealhnf(NF, D), cj(dgen) == dgen);
  my(Pp = cmtype(1), Pm = cmtype(-1));
  my(J0 = idealinv(NF, D));
  my(rp = decide(J0, Pp), rm = decide(J0, Pm));
  report("End = O_K, Phi_+", rp); report("End = O_K, Phi_-", rm);
  for (k = 1, min(2, #rp[5]), verify(idealhnf(NF, rp[5][k][1]), rp[5][k][2], Pp, "O_K ppav"));
  \\ primes above ell: count and inertness in K/K+
  my(Pl = idealprimedec(NF, ell));
  printf("  primes of K above ell: %d, each conj-stable: %s, f=%s\n", #Pl, vector(#Pl, i, idealhnf(NF, cjI(Pl[i])) == idealhnf(NF, Pl[i])), vector(#Pl, i, Pl[i].f));
  my(U = lucasU(ell), F = factor(abs(U))[,1]);
  printf("  U_ell = %d, primes %s\n", U, F~);
  for (a = 1, #F, my(p = F[a], L = primesabove(p));
    my(data = vector(#L, i, [resid(L[i], PI), resid(L[i], PIB), L[i].f, idealhnf(NF, cjI(L[i])) == idealhnf(NF, L[i])]));
    printf("  p=%d: %d primes of K, residue degrees %s, #conj-stable %d\n", p, #L, Set(vector(#L, i, L[i].f)), sum(i = 1, #L, data[i][4]));
    \\ gluing classes: primes with equal (pi,pibar) residues (degree-1) or self-glued (residue in F_p, degree 2)
    my(conds = List());
    for (i = 1, #L,
      if (L[i].f == 1,
        for (j = i+1, #L, if (L[j].f == 1 && data[i][1] == data[j][1] && data[i][2] == data[j][2],
          my(sc = (idealhnf(NF, cjI(L[i])) == idealhnf(NF, L[j])));
          printf("    glued pair (%d,%d), self-conjugate pair: %d\n", i, j, sc);
          if (sc, listput(conds, [idealmul(NF, L[i], L[j]), Str("split self-conj pair at ", p)])))),
        my(r1 = data[i][1]);
        if (r1^p == r1 && data[i][2] == r1,
          printf("    self-glued inert place %d: pi == pibar mod Q and lies in F_p; conj-stable %d\n", i, data[i][4]);
          listput(conds, [idealhnf(NF, L[i]), Str("self-glued inert place at ", p)]))));
    for (c = 1, #conds,
      my(Q = conds[c][1], J = idealmul(NF, J0, idealinv(NF, Q)));
      my(rq = decide(J, Pp), rq2 = decide(J, Pm));
      report(Str("End = R (conductor: ", conds[c][2], "), Phi_+"), rq);
      report(Str("End = R (conductor: ", conds[c][2], "), Phi_-"), rq2);
      for (k = 1, min(2, #rq[5]), my(S = rq[5][k][1], lam = rq[5][k][2]);
        my(B = sublattice(S, Q)); verify(B, lam, Pp, "R ppav  T = Z t + Q S"));
      \\ negative control: the same lambda on the O_K-lattice S is NOT unimodular (degree N(Q))
      if (#rq[5], verify(idealhnf(NF, rq[5][1][1]), rq[5][1][2], Pp, "neg. control: same lambda on S=O_K T"))));
}
