\\ Narrow-class analysis in K+ for ell in {11,17,29,37}; GRH-conditional (bnfinit not certified).
\r kplus_setup.gp
lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }
{
for (i = 1, 4, my(ell = [11,17,29,37][i]);
  my(S = setupKp(ell), P = S[1], c = S[2], s = S[3]);
  my(bnf = bnfinit(nfinit([P, [2,3,5,7,ell]]), 1), nf = bnf.nf, n = poldegree(P));
  my(bnr = bnrinit(bnf, [1, vector(n, j, 1)], 1));
  my(N = rnfnormgroup(bnr, varhigher("t")^2 + 7));   \\ subgroup of Cl^+ that are norms from K=K+(sqrt-7)
  my(cyc = bnr.cyc);
  my(sg = bnfsignunit(bnf));
  my(rk = matrank(Mod(concat(Mat(vectorv(n,j,1)), apply(t -> (1 - t)/2, sg)), 2)));
  printf("ell=%d n=%d h=%d Cl+ cyc=%s  [Cl+:norms]=%d  rank F2 of sign map on {-1,fu}=%d (=> h+/h = 2^%d)\n",
     ell, n, bnf.no, cyc, matdet(N), rk, n - rk);
  my(innorm(e) = denominator(matsolve(N, e)) == 1);
  \\ q above ell
  my(Lq = idealprimedec(nf, ell)); if (#Lq != 1 || Lq[1].e != n, error("q not totally ramified"));
  my(q = Lq[1], eq = bnrisprincipal(bnr, q, 0));
  my(epsK = kronecker(-7, ell));
  printf("  q: e=%d f=%d; narrow class %s; in norm group %d ; eps=%d\n", q.e, q.f, eq, innorm(eq), epsK);
  \\ v0 places for each inert-type p | U_ell
  my(U = lucasU(ell), F = factor(abs(U))[,1]);
  my(d = (s + 7*c)/2);  \\ = sqrt-7 * (pi - pibar), up to sign
  for (a = 1, #F, my(p = F[a], L = idealprimedec(nf, p));
    my(hits = select(pr -> idealval(nf, lift(d), pr) > 0, L));
    printf("  p=%d (p mod ell=%d, (-7/p)=%d): %d places in K+, f=%s; places with v(sqrt-7*(pi-pibar))>0: %d\n",
       p, p % ell, kronecker(-7,p), #L, Set(vector(#L,j,L[j].f)), #hits);
    for (h = 1, #hits, my(v = hits[h], ev = bnrisprincipal(bnr, v, 0));
       printf("    v0: f=%d val=%d narrow class %s in norm group %d\n", v.f, idealval(nf, lift(d), v), ev, innorm(ev)));
    );

);
}
