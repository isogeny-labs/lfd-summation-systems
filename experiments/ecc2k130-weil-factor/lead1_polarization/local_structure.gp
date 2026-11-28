\\ Lead 1, method 2b (PARI/GP): local structure of Z[pi,pibar] at the primes dividing U_ell.
\\ usage: gp -q < local_structure.gp
\\ For each ell and each p | U_ell:
\\  (a) p-part of [O_K : Z[pi,pibar]] from the HNF of the monomials mod p^2 together with the rank mod p;
\\  (b) residue points: values (pi, pibar) at the 2(ell-1) embeddings K -> Fpbar, the gluing
\\      (i,a) ~ (i+k,b), and the conj-orbits (= places of K+); the unique self-glued place v0.
default(parisizemax, 4000000000)
default(threadsizemax, 2000000000)
lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }

\\ O_K basis: z^j * w^e, j<ell-1, e<2, alpha = w (root of w^2+w+2). Coordinates of an element mod m.
coords(e, ell) = { my(d = ell-1, v = vector(2*d), p = lift(lift(e)));
  for (k = 0, 1, my(c = polcoef(p, k, 'w)); for (i = 0, d-1, v[k*d+i+1] = polcoef(c, i, 'z))); v~; }

pindex(ell, p) = {
  my(d = ell-1, m = p^2, Pz = polcyclo(ell, 'z), Pw = 'w^2 + 'w + 2);
  my(z = Mod(Mod(1, m)*'z, Pz), one = Mod(Mod(1,m)*'z^0, Pz));
  my(w = Mod(one*'w, Pw), Z = Mod(z*'w^0, Pw));
  my(pi = w*Z, pb = (-1 - w)*Z^(ell-1), cols = List(), x = Mod(one*'w^0, Pw));
  for (i = 0, 2*d-1, listput(cols, coords(x, ell)); x *= pi);
  x = pb; for (j = 1, 2*d-1, listput(cols, coords(x, ell)); x *= pb);
  my(M = lift(matconcat(Vec(cols))));
  my(H = mathnfmodid(M, m), v2 = valuation(prod(i = 1, #H, H[i,i]), p));
  my(r1 = matrank(M*Mod(1,p)));
  [2*d - r1, v2];
}

residue(ell, p) = {
  my(d = ell - 1, T = ffinit(p, 2, 't), g = ffgen(T, 't));
  my(R = polrootsmod(polcyclo(ell, 'x), [T, p]), om1 = subst(lift(R[1]), 't, g));
  my(ab = polrootsmod('x^2 + 'x + 2, [T, p]), a = subst(lift(ab[1]), 't, g), b = -1 - a);
  \\ embedding (i,s): zeta -> om1^i, alpha -> a (s=1) or b (s=2)
  my(val(i, s) = my(al = if (s == 1, a, b), be = if (s == 1, b, a)); [al*om1^i, be*om1^(-i)]);
  my(pts = Map(), keys = List());
  for (i = 1, d, for (s = 1, 2, my(v = val(i, s)); my(k = Str(v)); if (mapisdefined(pts, k), my(o = mapget(pts, k)); mapput(pts, k, concat(o, [[i,s]])), mapput(pts, k, [[i,s]]); listput(keys, k))));
  my(sizes = vector(#keys, j, #mapget(pts, keys[j])));
  \\ conj orbit of (i,1) is {(i,1),(-i,2)} as embeddings into C; over Fpbar with p = -1 mod ell the
  \\ Frobenius orbit is {(i,a),(i,a)^p}. Report: number of residue points, and the glued classes.
  my(glued = select(j -> sizes[j] == 2, [1..#keys]));
  printf("  ell=%d p=%d: residue points of Z[pi,pibar] at embeddings: %d (expected ell=%d); classes of size 2: %d, size 1: %d\n",
    ell, p, #keys, ell, #glued, #keys - #glued);
  \\ self-glued test: a glued class {(i,1),(j,2)} with j = -i mod ell (i.e. the class is a conj-orbit)
  my(selfg = 0, selfi = 0);
  for (t = 1, #glued, my(c = mapget(pts, keys[glued[t]])); my(e1 = c[1], e2 = c[2]);
    if (e1[2] != e2[2] && Mod(e1[1] + e2[1], ell) == 0, selfg++; selfi = if (e1[2] == 1, e1[1], e2[1])));
  my(frobfix = 0);
  if (selfg, my(v = val(selfi, 1)); frobfix = (v[1]^p == v[1]) && (v[1] == v[2]) );
  printf("    self-glued conj-orbits: %d ; at it pi == pibar mod p and pi mod p in F_p: %d ; zeta mod p in F_p there: %d\n",
    selfg, frobfix, if (selfg, (om1^selfi)^p == om1^selfi, -1));
  [#keys, #glued, selfg, frobfix];
}

{
  foreach([11, 17, 19, 29, 131], ell,
    my(U = lucasU(ell), F = factor(abs(U)));
    printf("ell=%d U=%d factor=%s ; all factors proven prime (isprime): %d ; (-7/ell)=%d\n", ell, U, F, vecprod(apply(isprime, F[,1]~)), kronecker(-7, ell));
    for (k = 1, #F~, my(p = F[k,1], t = getabstime(), r = pindex(ell, p));
      printf("  p=%d: p mod ell=%d, (-7/p)=%d ; #elementary divisors divisible by p = %d ; sum min(v,2) = %d  => v_p(index) = %s (expected v_p(U)*(ell-2) = %d) [%d ms]\n",
        p, p % ell, kronecker(-7, p), r[1], r[2], if (r[1] == r[2], r[2], "?"), F[k,2]*(ell-2), getabstime()-t);
      residue(ell, p)));
}
