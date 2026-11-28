\\ (a) Local Hilbert symbols (r,-7)_v and sign parity for the polarization candidates r (product formula check)
\\ (b) Unconditional proof that q is not narrowly principal at ell=29,37: a quadratic L/K+ unramified at all
\\     finite primes (checked with rnfdisc) in which q does not split.
\r construct.gp
nnegs(x) = #select(t -> t < 0, nfeltsign(NF, nfalgtobasis(NF, x)));
hilb(r, pr) = nfhilbert(NF, nfalgtobasis(NF, r), -7, pr);
\\ all finite places where the symbol can be -1: primes dividing 2*7*num/den of N(r)
symbols(r, label) = {
  my(Nr = abs(norm(r)), S = Set(concat([2, 7], concat(factor(numerator(Nr))[,1]~, factor(denominator(Nr))[,1]~))));
  my(out = List(), prod = 1);
  for (i = 1, #S, my(L = idealprimedec(NF, S[i])); for (j = 1, #L, my(h = hilb(r, L[j])); prod *= h; if (h == -1, listput(out, [S[i], L[j].f, idealval(NF, nfalgtobasis(NF, r), L[j])]))));
  my(ng = nnegs(r));
  printf("  [%s] #negative real embeddings=%d ; finite places with symbol -1 (p, f, v(r)) = %s ; product finite=%d * infinite=%d = %d\n",
     label, ng, Vec(out), prod, (-1)^ng, prod*(-1)^ng);
}
{
for (ii = 1, 3, my(ell = [11, 29, 37][ii]);
  init(ell);
  my(q = qideal(), PR = factor(abs(lucasU(ell)))[,1], inert = select(p -> p % ell == ell-1, PR));
  \\ T = O_K candidate (generator of q, possibly non-definite)
  my(g1 = Mod(lift(nfbasistoalg(NF, bnfisprincipal(BNF, q)[2])), P)); symbols(lamfrom(g1)[1], Str("ell=", ell, " r for T=O_K (any generator of q)"));
  for (a = 1, #inert, my(v0 = v0ideal(inert[a]), I = idealmul(NF, q, v0));
     my(g = Mod(lift(nfbasistoalg(NF, bnfisprincipal(BNF, I)[2])), P)); symbols(lamfrom(g)[1], Str("ell=", ell, " r for R(v0 ", inert[a], "), n=1")));
  if (#inert == 2, my(I = idealmul(NF, q, idealmul(NF, v0ideal(inert[1]), v0ideal(inert[2]))), g = posgen(I));
     symbols(lamfrom(g)[1], Str("ell=", ell, " r for R(v0,v0'), n=2 (the constructed one)")));
  \\ (b) unramified quadratic extensions of K+ (narrow 2-class field), q's splitting
  my(Py = subst(P, 'x, 'y), bnfy = bnfinit(nfinit([Py, [2,3,5,7,ell]]), 1), nfy = bnfy.nf);
  my(bnry = bnrinit(bnfy, [1, vector(poldegree(Py), j, 1)], 1), qy = idealprimedec(nfy, ell));
  if (#qy != 1, error("q in y-model"));
  my(pols = bnrclassfield(bnry, 2), qinert = 0);
  printf("  narrow Cl+ = %s ; %d quadratic generators of the 2-elementary narrow class field\n", bnry.cyc, #pols);
  for (k = 1, #pols, my(f = pols[k]);
    my(d = rnfdisc(nfy, f)[1], unr = (idealnorm(nfy, d) == 1));
    my(rnf = rnfinit(nfy, f), nq = #rnfidealprimedec(rnf, qy[1]));
    printf("   L_%d: rel disc norm %d (unramified at finite places: %d); #primes of L over q = %d -> q %s\n", k, idealnorm(nfy, d), unr, nq, if (nq == 2, "splits", "does NOT split"));
    if (unr && nq == 1, qinert = 1));
  printf("  => unconditional: q not narrowly principal: %d\n\n", qinert));
}
