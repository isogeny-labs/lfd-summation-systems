\\ Independent unconditional check that q (ramified prime above ell in K+) has no totally
\\ definite generator, via explicit units u of K+ with K+(sqrt u)/K+ unramified at all finite
\\ primes and u a non-square mod q. No bnfinit/bnrclassfield used.
default(parisizemax, 4000000000);
\\ Chebyshev U_{a-1}(c) = (z^a - z^-a)/(z - z^-1)
W = varhigher("W");
ff(t) = { if (type(t) == "t_FFELT", t = t.pol); t = lift(t); if (type(t)=="t_POL", t = polcoeff(t,0)); t; }
cheb(a, c) = { my(u0 = 0, u1 = 1); for (i = 2, a, [u0, u1] = [u1, c*u1 - u0]); u1; }
run(ell) = {
  my(R = polresultant(polcyclo(ell,'z), 'z^2 - 'y*'z + 1, 'z), fc = factor(R)[1,1]);
  my(Kc = nfinit(fc));
  my(eq = rnfequation(Kc, 'x^2 - 7*(4 - 'y^2), 1));
  my(P = eq[1], c = eq[2], k = eq[3], s = Mod('x, P) - k*c);
  my(nf = nfinit(P), n = poldegree(P));
  printf("ell=%d deg K+=%d  disc K+ = %s\n", ell, n, factor(nf.disc));
  \\ K/K+ = K+(sqrt-7) unramified at finite primes:
  my(dK = rnfdisc(nf, W^2 + 7)); printf("  rnfdisc(K+(sqrt-7)) finite part norm = %s\n", idealnorm(nf, dK[1]));
  my(Lq = idealprimedec(nf, ell)); if (#Lq != 1 || Lq[1].e != n, error("q"));
  my(q = Lq[1]);
  \\ candidate units: -1, cyclotomic units U_{a-1}(c), a=2..(ell-1)/2 (in Q(zeta)^+ subset K+), plus
  \\ units x - y with x,y such that... just cyclotomic ones and their K+-conjugates via s? keep simple.
  my(U = List([Mod(-1,P)]));
  for (a = 2, (ell-1)/2, listput(U, Mod(lift(subst(cheb(a, 'c), 'c, c)), P)));
  \\ also units from other K+ structure: (s + 7 c)/2 type elements are not units; skip.
  if (USEBNF, my(bnf = bnfinit(nf, 1)); for (i = 1, #bnf.fu, listput(U, Mod(lift(bnf.fu[i]), P))));
  U = Vec(U);
  for (i = 1, #U, if (abs(norm(U[i])) != 1, error("not a unit")));
  \\ residues mod q (F_ell) and squareness: Legendre
  my(rq = vector(#U, i, kronecker(ff(nfmodpr(nf, nfalgtobasis(nf, U[i]), nfmodprinit(nf, q))), ell)));
  \\ structure of (O/4O)^*
  my(bid = idealstar(nf, 4, 2), cyc = bid.cyc);
  my(Lg = matconcat(vector(#U, i, ideallog(nf, nfalgtobasis(nf, U[i]), bid))));
  \\ square in (O/4)^* iff log coords even at even cyc entries
  my(ev = select(j -> cyc[j] % 2 == 0, [1..#cyc]));
  my(A = Mod(matrix(#ev, #U, a, b, Lg[ev[a], b]), 2));
  my(sq = Mod(vector(#U, i, (1 - rq[i])/2), 2));  \\ 1 = nonsquare mod q
  my(Kr = matker(A));
  printf("  #units %d; (O/4)^* cyc len %d (even %d); kernel dim (unramified at 2 subgroup) %d\n", #U, #cyc, #ev, #Kr);
  \\ find kernel vector with odd non-square parity at q
  my(found = 0);
  for (j = 1, #Kr, if (sq * Kr[,j] != 0, found = Kr[,j]; break));
  if (!found, printf("  no cyclotomic-unit witness found\n"); return(0));
  my(e = lift(found), u = prod(i = 1, #U, U[i]^e[i]));
  \\ certify: rnfdisc trivial at finite places, u nonsquare mod q
  my(dL = rnfdisc(nf, W^2 - lift(u)));
  my(t = ff(nfmodpr(nf, nfalgtobasis(nf, u), nfmodprinit(nf, q))));
  my(sg = nfeltsign(nf, nfalgtobasis(nf, u)));
  my(dec = rnfidealprimedec(rnfinit(nf, W^2 - lift(u)), q));
  printf("  witness exps %s\n  rnfdisc finite norm = %s ; u mod q = %d, Legendre = %d ; #primes above q in L = %d ; #neg embeddings of u = %d\n",
    e~, idealnorm(nf, dL[1]), t, kronecker(t, ell), #dec, #select(x->x<0, sg));
  \\ Hilbert-reciprocity consequence for the ideal needed by T=O_K: I = 7^-1 (c^2-4)^-k q^-1
  \\ check: for ANY r generating I, (# real places where r<0 and u<0) is odd, so r is not
  \\ totally positive, and -r neither (#neg(u) even).
  1;
}
USEBNF = 0; run(11); run(29); run(37);
print("--- with additional units found by bnfinit (GRH used only to FIND them; each verified to be a unit) ---");
USEBNF = 1; run(11); run(29); run(37);
