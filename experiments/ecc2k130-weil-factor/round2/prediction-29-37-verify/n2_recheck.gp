\\ Independent re-construction and re-check of "End = O_K ppav exists" at ell = 29, 37 (and 11 control).
\\ Own code: K+ = Q(c,s) via rnfequation (no polredbest); K as tower Q(alpha)(zeta) with manual traces.
\\ GRH (bnfinit) is used ONLY to find generators / a suitable split prime; all checks are exact.
default(parisizemax, 4000000000);
ZV = varhigher("ZV"); AV = varlower("AV");
ff(t) = { if (type(t) == "t_FFELT", t = t.pol); t = lift(t); if (type(t)=="t_POL", t = polcoeff(t,0)); t; }
setup(ell) = {
  ELL = ell;
  my(R = polresultant(polcyclo(ell,'z), 'z^2 - 'y*'z + 1, 'z), fc = factor(R)[1,1]);
  my(eq = rnfequation(nfinit(fc), 'x^2 - 7*(4 - 'y^2), 1));
  P = eq[1]; Cc = eq[2]; KK = eq[3]; Ss = Mod('x, P) - KK*Cc;
  NF = nfinit(P); BNF = bnfinit(NF, 1); NN = poldegree(P);
  Q = idealprimedec(NF, ell)[1];
  PHI = polcyclo(ell, ZV); AL = Mod(AV, AV^2 + AV + 2); ZE = Mod(ZV, PHI);
  N = 2*(ell-1);
  BAS = vector(N, i, my(j = (i-1) % (ell-1), e = (i-1) \ (ell-1)); ZE^j * AL^e);
  CBAS = vector(N, i, my(j = (i-1) % (ell-1), e = (i-1) \ (ell-1)); ZE^(-j) * (-1-AL)^e);
  SU = apply(t -> (1-t)/2, bnfsignunit(BNF));
}
\\ K+ element (polmod mod P in x) -> tower, using x = s + KK c
toT(r) = { my(cT = ZE + 1/ZE, sT = (2*AL+1)*(ZE - 1/ZE)); subst(lift(r), 'x, sT + KK*cT); }
trT(y) = { my(f = lift(y), t = 0);
  for (j = 0, ELL-2, my(cj = polcoeff(f, j, ZV)); if (cj == 0, next);
     my(cl = lift(cj), u = polcoeff(cl, 0, AV), v = polcoeff(cl, 1, AV));
     t += (2*u - v) * if (j == 0, ELL-1, -1));
  t; }
coordsT(y) = { my(f = lift(y)); vectorv(N, i, my(j = (i-1) % (ELL-1), e = (i-1) \ (ELL-1)); polcoeff(lift(polcoeff(f, j, ZV)), e, AV)); }
\\ definite generator of ideal I: use bnfisprincipal + unit sign correction (GRH only for finding)
defgen(I) = {
  my(t = bnfisprincipal(BNF, I)); if (t[1] != 0*t[1], error("not principal"));
  my(g = t[2], sv = apply(x -> (1-x)/2, nfeltsign(NF, g))~);
  my(A = concat(Mat(vectorv(NN, j, 1)), SU), sol = matsolvemod(A, 2, sv));
  if (type(sol) == "t_INT", return(0));
  my(u = (-1)^sol[1] * prod(i = 1, #BNF.fu, nfbasistoalg(NF, BNF.fu[i])^sol[i+1]));
  my(G = Mod(lift(nfbasistoalg(NF, g) * u), P));
  if (#Set(nfeltsign(NF, nfalgtobasis(NF, G))) != 1, error("sign fix failed"));
  G; }
\\ r with r O = (7 (c^2-4)^k q I)^-1, lambda = sqrt-7 * r
mkr(G) = { my(k = (ELL-3)/2); 1/(7 * (Cc^2-4)^k * G); }
checkall(r, M, extra, label) = {
  my(k = (ELL-3)/2);
  \\ (i) ideal criterion in K+
  my(target = idealinv(NF, idealmul(NF, idealmul(NF, 7*(Cc^2-4)^k, Q), extra)));
  my(idok = idealhnf(NF, idealdiv(NF, r, target)) == matid(NN)); if (!idok, print("   r: ", idealfactor(NF, r)[,2]~, " target: ", idealfactor(NF, target)[,2]~));
  my(sg = nfeltsign(NF, nfalgtobasis(NF, r)));
  \\ (ii) Gram in own tower arithmetic
  my(lam = (2*AL+1) * toT(r));
  my(LB = vector(N, i, lam*BAS[i]));
  my(G0 = matrix(N, N, i, j, trT(LB[i]*CBAS[j])));
  my(G = M~ * G0 * M);
  \\ stability of lattice M under pi, pibar
  my(Mi = M^-1, PI = AL*ZE, PIB = (-1-AL)/ZE);
  my(stab(y) = denominator(Mi * matconcat(vector(N, i, coordsT(y * (BAS*M[,i])))) ) == 1);
  printf("  [%s] ideal criterion=%d  signs: #neg=%d of %d  Gram: integral=%d alternating=%d det=%s  pi/pibar-stable=%d/%d\n",
     label, idok, #select(x->x<0, sg), #sg, denominator(G)==1, G+G~==0, matdet(G), stab(PI), stab(PIB));
}
run(ell) = {
  setup(ell);
  printf("ell=%d\n", ell);
  \\ T = O_K: needs definite generator of q
  my(g = defgen(Q));
  printf("  T=O_K: definite generator of q found (GRH units)? %d\n", g != 0);
  if (g, checkall(mkr(g), matid(N), 1, "T=O_K"),
     my(g1 = Mod(lift(nfbasistoalg(NF, bnfisprincipal(BNF, Q)[2])), P)); checkall(mkr(g1), matid(N), 1, "NEG T=O_K nondefinite"));
  if (g, return);
  \\ split prime search: choose residue data (z0,a0) in F_p, define v = ker on K+, test [q v] narrowly trivial
  forprime(p = 2, 10^6, if (p % ell != 1 || kronecker(-7, p) != 1, next);
    my(z0 = polrootsmod(polcyclo(ell,'z), p)[1], a0 = polrootsmod('a^2+'a+2, p)[1]);
    my(c0 = z0 + 1/z0, s0 = (2*a0+1)*(z0 - 1/z0), x0 = lift(s0 + KK*c0));
    my(v = idealadd(NF, p, lift(Mod('x,P)) - x0));
    if (idealnorm(NF, v) != p, next);
    my(G = defgen(idealmul(NF, Q, v))); if (!G, next);
    \\ T = kernel of O_K -> F_p, z->z0, alpha->a0
    my(f = vector(N, i, my(j = (i-1) % (ell-1), e = (i-1) \ (ell-1)); lift(z0^j * a0^e)));
    my(M = matkermod(Mat(f), p));  \\ not a lattice; build properly below
    my(L = matconcat([matrix(N,N,i,j, if(i==j,p,0)), matrix(N,1)]));
    my(k0 = 1, inv = lift(1/Mod(f[1], p)), Mt = matid(N));
    for (j = 2, N, Mt[1, j] = -lift(Mod(f[j]*inv, p))); Mt[1,1] = p;
    printf("  split p=%d with [q v] narrowly trivial (GRH search); index of T = %d\n", p, abs(matdet(Mt)));
    checkall(mkr(G), Mt, v, Str("T=P above ", p));
    break);
}
\\ n=2 at ell=37: order non-maximal at both self-glued places (p = 73, 2663); T = R itself.
selfglued(p) = {
  my(t = ffgen(p^2, 't), a0 = polrootsmod('a^2+'a+2, t)[1], b0 = -1 - a0);
  my(zs = polrootsmod(polcyclo(ELL, 'z), t), z0 = 0);
  for (i = 1, #zs, if (a0*zs[i] == b0/zs[i], z0 = zs[i]; break));
  if (z0 == 0, error("no self-glued z0"));
  my(c0 = z0 + 1/z0, s0 = (2*a0+1)*(z0 - 1/z0), x0 = s0 + KK*c0);
  if (x0^p != x0, error("x0 not in F_p"));
  my(v = idealadd(NF, p, lift(Mod('x,P)) - ff(x0)));
  if (idealnorm(NF, v) != p, error("v0 not prime of norm p"));
  \\ functional: coefficient of t in image (image in F_p iff 0)
  my(f = vector(N, i, my(j = (i-1) % (ELL-1), e = (i-1) \ (ELL-1), w = (z0^j * a0^e).pol); lift(polcoeff(w, 1))));
  [v, f];
}
kerlat(f, p) = { my(k0 = 0); for (k = 1, #f, if (f[k] % p, k0 = k; break));
  my(inv = lift(1/Mod(f[k0], p)), M = matid(#f));
  for (j = 1, #f, if (j == k0, M[k0, j] = p, M[k0, j] = -lift(Mod(f[j]*inv, p)))); M; }
{
setup(37); printf("ell=37 n=2 recheck\n");
my(A = selfglued(73), B = selfglued(2663));
my(M1 = kerlat(A[2], 73), f2 = (B[2] * M1) % 2663, M = M1 * kerlat(f2, 2663));
my(G = defgen(idealmul(NF, Q, idealmul(NF, A[1], B[1]))));
printf("  definite generator for q v0 v0' found: %d ; index = %d\n", G != 0, abs(matdet(M)));
if (G, checkall(mkr(G), M, idealmul(NF, A[1], B[1]), "R(v0,v0')"));
my(Mz = M^-1 * matconcat(vector(N, i, coordsT(ZE * (BAS*M[,i])))));
printf("  zeta in End(R)? %d\n", denominator(Mz) == 1);
\\ negative control: v0 over 73 only
my(G1 = defgen(idealmul(NF, Q, A[1])));
printf("  n=1 (73 only): definite generator found (GRH units): %d\n", G1 != 0);
}
