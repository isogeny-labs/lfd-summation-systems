\\ Method A (PARI finite fields): for each odd prime p | U_ell, list the R = Z[pi,pibar]-residue data at the
\\ places of K = Q(zeta_ell, sqrt-7) above p, find self-glued places (pi = pibar in F_p), and the order of
\\ lambda = pi mod v0 (lambda^2 = 2).  No number field of degree 2(ell-1) is built.
default(parisizemax, 4000000000);
lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }
run(ell) = {
  my(U = abs(lucasU(ell)), F = factor(U)[,1]);
  printf("ell=%d |U|=%s primes=%s eps=(-7/ell)=%d\n", ell, U, F~, kronecker(-7, ell));
  for (i = 1, #F, my(p = F[i]); if (p == 2, next);
    my(f = 1); while ((p^f - 1) % ell, f++);
    my(t = ffgen(ffinit(p, lcm(2, f)), 't), roots_a = polrootsmod('x^2 + 'x + 2, t), z, out = List());
    \\ primitive ell-th root of unity in F_{p^lcm(2,f)}
    my(q = p^lcm(2, f)); z = random(t)^((q - 1)/ell); while (z == 1, z = random(t)^((q - 1)/ell));
    my(pairs = List());
    for (ia = 1, #roots_a, my(a = roots_a[ia], b = 2/a);
      for (j = 1, ell - 1, listput(pairs, [a*z^j, b*z^(-j)])));
    \\ group homs by the R-maximal ideal: two homs give the same m iff the pair (pi,pibar) agree up to p-Frobenius
    my(keys = vector(#pairs, k, my(P = pairs[k], o = List()); for (s = 0, 2*f, listput(o, [P[1]^(p^s), P[2]^(p^s)])); vecsort(apply(x -> Str(x), Vec(o)))[1]));
    my(nm = #Set(keys));
    my(sg = [k | k <- [1..#pairs], pairs[k][1] == pairs[k][2]]);
    printf("  p=%s  p mod ell=%d (-7/p)=%d  type=%s  homs=%d  self-glued homs (pi==pibar)=%d  distinct R-maximal ideals=%d\n", p, p % ell, kronecker(-7, p),
      if (p % ell == 1 && kronecker(-7,p) == 1, "split", if ((p+1) % ell == 0 && kronecker(-7,p) == -1, "inert", "other")), #pairs, #sg, nm);
    if (#sg,
      my(lam = pairs[sg[1]][1], ok = (lam^2 == 2), inFp = (lam^p == lam));
      my(L = Set(apply(k -> pairs[k][1], sg)));
      printf("    distinct lambda values at self-glued homs: %d ; lambda^2==2: %d ; lambda in F_p: %d\n", #L, ok, inFp);
      my(lz = lift(minpoly(lam)));  \\ degree 1: lambda = -coeff
      my(lv = Mod(-polcoef(lz, 0), p));
      my(o1 = znorder(lv), o2 = znorder(-lv), o2_ = znorder(Mod(2, p)));
      printf("    lambda = %s ; ord(lambda) = %s = (p-1)/%s ; ord(-lambda) = %s ; ord(2) = %s ; p-1 = %s\n", lift(lv), o1, (p-1)/o1, o2, o2_, factor(p-1));
    );
  );
}
run(11); run(17); run(19); run(131);
\\ negative controls: primes not dividing U_ell give no self-glued places (pi == pibar mod v forces p | U_ell)
