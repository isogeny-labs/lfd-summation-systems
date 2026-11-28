\r construct.gp
\\ find a split place v (p = 1 mod ell, (-7/p)=1) with q*v narrowly principal
findsplit(q) = {
  forprime(p = 2, 10^7, if (p % ELL != 1 || kronecker(-7, p) != 1, next);
    my(L = idealprimedec(NF, p));
    for (i = 1, #L, my(I = idealmul(NF, q, L[i]));
      if (bnrisprincipal(BNR, I, 0) == 0*BNR.cyc~, return([p, L[i]]))));
}
signvec(x) = apply(t -> (1 - t)/2, nfeltsign(NF, x));
\\ GRH-dependent: can the sign vector of x be made all-0 or all-1 by a unit of K+?
unitfix(x) = {
  my(n = poldegree(P), Su = apply(t -> (1 - t)/2, bnfsignunit(BNF)));
  my(A = concat(Mat(vectorv(n, j, 1)), Su), sv = signvec(x)~);
  my(s1 = matsolvemod(A, 2, sv), s2 = matsolvemod(A, 2, (sv + vectorv(n,j,1)) % 2));
  [type(s1) != "t_INT", type(s2) != "t_INT", vecsum(sv)];
}
{
for (ii = 1, 3, my(ell = [11, 29, 37][ii], t0 = getabstime());
  init(ell);
  my(q = qideal(), U = lucasU(ell), PR = factor(abs(U))[,1]);
  my(eps = kronecker(-7, ell));
  printf(" eps=%d, U=%d, primes %s\n", eps, U, PR~);
  \\ (1) A_ell itself: T = O_K, needs q narrowly principal
  my(g = posgen(q));
  printf(" (1) T = O_K (A_ell itself): q narrow class %s -> %s\n", bnrisprincipal(BNR, q, 0)~, if (g, "narrowly principal: construct", "NOT narrowly principal"));
  if (g, my(L = lamfrom(g)); verifyT(matid(N), L[2], L[1], Str("ell=", ell, " T=O_K")));
  if (!g, my(g1 = nfbasistoalg(NF, bnfisprincipal(BNF, q)[2]), L = lamfrom(Mod(lift(g1), P)));
     printf("     q = (g1) with g1 of #negative signs %d; unit correction possible (GRH units)? to + : %d, to - : %d\n", unitfix(nfalgtobasis(NF, g1))[3], unitfix(nfalgtobasis(NF, g1))[1], unitfix(nfalgtobasis(NF, g1))[2]);
     verifyT(matid(N), L[2], L[1], Str("ell=", ell, " NEG CONTROL T=O_K with lambda from non-definite generator")));
  \\ (2) some T with End = O_K: T = P above a split place v with [q v] = 0 in Cl+
  if (!g, my(sv = findsplit(q), p = sv[1], v = sv[2]);
     my(gv = posgen(idealmul(NF, q, v)), L = lamfrom(gv));
     my(M = kerlat(splitfunctional(p, v), p));
     printf(" (2) End = O_K, T = P | v, v split place above p=%d with [q v]=0 in Cl+\n", p);
     verifyT(M, L[2], L[1], Str("ell=", ell, " T=P (p=", p, ")")));
  \\ (3) orders non-maximal at self-glued places v0 (inert-type primes)
  my(inert = select(p -> p % ell == ell-1, PR));
  for (a = 1, #inert, my(p = inert[a], v0 = v0ideal(p), f = v0functional(p, v0), M = kerlat(f, p));
     my(I = idealmul(NF, q, v0), gg = posgen(I));
     printf(" (3) R non-max at v0 above %d only (n=1): [q v0] in Cl+ = %s; narrowly principal: %d\n", p, bnrisprincipal(BNR, I, 0)~, gg != 0);
     if (gg, my(L = lamfrom(gg)); verifyT(M, L[2], L[1], Str("ell=", ell, " R(v0 ", p, ")")),
        my(g1 = nfbasistoalg(NF, bnfisprincipal(BNF, I)[2]), L = lamfrom(Mod(lift(g1), P)), uf = unitfix(nfalgtobasis(NF, g1)));
        printf("     unit correction to definite possible (GRH units)? +:%d -:%d\n", uf[1], uf[2]);
        verifyT(M, L[2], L[1], Str("ell=", ell, " NEG R(v0 ", p, ") with non-definite lambda"))));
  if (#inert == 2, my(p1 = inert[1], p2 = inert[2], v1 = v0ideal(p1), v2 = v0ideal(p2));
     my(M1 = kerlat(v0functional(p1, v1), p1), f2 = v0functional(p2, v2) * M1, M2 = kerlat(apply(t -> t % p2, f2), p2), M = M1*M2);
     my(I = idealmul(NF, q, idealmul(NF, v1, v2)), gg = posgen(I));
     printf(" (4) R non-max at both v0 (n=2): [q v0 v0'] in Cl+ = %s; narrowly principal: %d\n", bnrisprincipal(BNR, I, 0)~, gg != 0);
     if (gg, my(L = lamfrom(gg)); verifyT(M, L[2], L[1], Str("ell=", ell, " R(v0,v0')"))));
  printf(" [ell=%d done in %d ms]\n\n", ell, getabstime() - t0));
}
