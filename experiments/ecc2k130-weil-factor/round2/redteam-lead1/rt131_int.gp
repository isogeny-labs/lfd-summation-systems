\\ red-team (redteam-lead1): integer facts at ell = 131 from the closed form U = (alpha^131 - beta^131)/(alpha - beta)
rt_U(n) = { my(al = Mod('w, 'w^2 + 'w + 2), be = -1 - al, q = lift((al^n - be^n)/(al - be))); polcoef(q, 0); }
{
  my(U = rt_U(131), F = factor(abs(U)), p1 = 263, p2 = 146505763881528721);
  printf("U_131 = %d\nfactor = %s ; squarefree %d ; isprime (proven) %s\n", U, F, issquarefree(U), apply(isprime, F[,1]~));
  printf("(-7/131) = %d ; 263 mod 131 = %d, (-7/263) = %d ; p2 mod 131 = %d, (-7/p2) = %d\n", kronecker(-7,131), p1 % 131, kronecker(-7,p1), p2 % 131, kronecker(-7,p2));
  printf("v_131(D_K+) = e-1 = %d (odd) ; 131 mod 7 = %d\n", 129, 131 % 7);
  \\ order of 7 mod 131 (parity decides whether places above 7 are inert in K/K+; contribution is +1 either way)
  printf("ord(7 mod 131) = %d ; ord(2 mod 131) = %d\n", znorder(Mod(7,131)), znorder(Mod(2,131)));
  \\ every prime p | U_ell, ell <= 199, is of split or inert type; report squarefreeness of U_ell (hypothesis of the general rule)
  my(bad = List());
  forprime(l = 3, 113, if (l == 7, next); my(u = rt_U(l), f = factor(abs(u), 10^7));
    for (k = 1, #f[,1], my(p = f[k,1]); if (!ispseudoprime(p), next);
      my(t = if (p % l == 1 && kronecker(-7,p) == 1, "S", if (p % l == l-1 && kronecker(-7,p) == -1, "I", "X")));
      if (t == "X" || f[k,2] > 1, listput(bad, [l, p, t, f[k,2]]))));
  printf("ell <= 113: primes p | U_ell (fully factored part) of type other than split/inert, or with p^2 | U_ell: %s\n", Vec(bad));
}
