B1(ell, i, e) = {
  my(n = ell-1, g = znprimroot(ell), f = if(i % n == 0, 1, ell) * if(e, 7, 1), s = 0, z = exp(2*Pi*I/n), tab);
  tab = vector(ell-1); for(k = 0, n-1, tab[lift(g^k)] = z^((i*k) % n));
  for(a = 1, f-1,
    if(gcd(a, f) > 1, next);
    my(v = if(i % n, tab[a % ell], 1));
    if(e, v *= kronecker(-7, a));
    s += v*a);
  s / f;
}
logprodK(ell) = {
  my(n = ell-1, L = 0., P = 1.);
  for(i = 0, n-1, for(e = 0, 1,
     if((-1)^i * (-1)^e == -1, my(b = -B1(ell, i, e)/2); P *= b)));
  P;
}
logprodCyc(p) = {my(P = 1.); forstep(i = 1, p-2, 2, P *= -B1(p, i, 0)/2); P};
