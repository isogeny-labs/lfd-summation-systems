\\ Constructive principal polarizations at ell = 11 (control), 29, 37.
\\ Deligne modules T in K = Q(zeta_ell, sqrt-7), pi = alpha*zeta, pibar = beta/zeta.
\\ Polarization b(x,y) = Tr_{K/Q}(lambda x conj(y)), lambda = sqrt-7 * r, r in K+ totally positive
\\ (positive on the CM type Phi = {phi : Im phi(sqrt-7) > 0}).
\\ Everything that is CHECKED below (ideal equalities, Gram integrality / alternation / det, stability,
\\ signs) is unconditional; GRH (bnfinit) is used only to FIND generators.
\r kplus_setup.gp
lucasU(n) = { my(a = 0, b = 1); for (k = 1, n-1, [a, b] = [b, -b - 2*a]); b; }

\\ ---------- per-ell globals ----------
init(ell) = {
  ELL = ell; BB = setupBoth(ell);
  P = BB[1]; C1 = BB[2]; S1 = BB[3]; XK = BB[4]; POLK = BB[5]; ZE = BB[6]; AL = BB[7];
  BNF = bnfinit(nfinit([P, [2,3,5,7,ell]]), 1); NF = BNF.nf;
  if (NF.pol != P, error("bnf polynomial mismatch"));
  BNR = bnrinit(BNF, [1, vector(poldegree(P), j, 1)], 1);
  N = 2*(ell-1);
  BAS = vector(N, i, my(j = (i-1) % (ell-1), e = (i-1) \ (ell-1)); ZE^j * AL^e);
  CBAS = vector(N, i, my(j = (i-1) % (ell-1), e = (i-1) \ (ell-1)); ZE^(-j) * (-1-AL)^e);
  BMAT = matrix(N, N, a, i, polcoeff(lift(BAS[i]), a-1));
  BMATI = BMAT^-1;
  PI = AL*ZE; PIB = (-1-AL)/ZE;
  printf("ell=%d: deg K=%d, deg K+=%d\n", ell, N, poldegree(P));
}
ff2int(u) = { if (type(u) == "t_FFELT", u = u.pol); if (type(u) == "t_POL", if (poldegree(u) > 0, error("not in prime field")); u = polcoeff(u, 0)); lift(u); }
coords(y) = BMATI * vectorv(N, a, polcoeff(lift(y), a-1));
multmat(y) = matconcat(vector(N, i, coords(y*BAS[i])));
gram0(lam) = { my(LB = vector(N, i, lam*BAS[i])); matrix(N, N, i, j, trace(LB[i]*CBAS[j])); }

\\ lattice {x in Z^N : f.x = 0 mod p} as basis matrix (columns), f integer row vector
kerlat(f, p) = {
  my(k0 = 0); for (k = 1, #f, if (f[k] % p, k0 = k; break)); if (!k0, error("zero functional"));
  my(inv = lift(1/Mod(f[k0], p)), M = matrix(#f, #f));
  for (j = 1, #f, if (j == k0, M[k0, j] = p, M[j, j] = 1; M[k0, j] = -lift(Mod(f[j]*inv, p))));
  M;
}

\\ real-embedding signs of r in K+ (exact routine) and numeric Im(phi(lambda)) on Phi in K
signsKp(r) = nfeltsign(NF, nfalgtobasis(NF, r));
imOnPhi(lamK) = {
  my(R = polroots(POLK), out = List());
  for (j = 1, #R, my(sq = subst(lift(2*AL+1), 'x, R[j]));
     if (imag(sq) > 0, listput(out, imag(subst(lift(lamK), 'x, R[j])))));
  Vec(out);
}

\\ verify a polarized Deligne module (T = lattice M in O_K coordinates, lambda)
verifyT(M, lamK, rKp, label) = {
  my(G = M~ * gram0(lamK) * M, d = matdet(G));
  my(Mi = M^-1, stab(y) = denominator(Mi * multmat(y) * M) == 1);
  my(sg = signsKp(rKp), im = imOnPhi(lamK));
  my(res = [d, denominator(G) == 1, G + G~ == 0, conj_ok(lamK), stab(PI), stab(PIB), stab(ZE),
            #select(t -> t > 0, sg) == #sg, #im == N/2 && vecmin(im) > 0]);
  printf("  [%s] index[O_K:T]=%s det=%s integral=%d alternating=%d conj(lam)=-lam %d  pi,pibar-stable=%d,%d zeta-stable=%d  r totally positive (nfeltsign)=%d  Im phi(lam)>0 on all %d phi (numeric, min %.3e)=%d\n",
    label, abs(matdet(M)), d, res[2], res[3], res[4], res[5], res[6], res[7], res[8], #im, vecmin(im), res[9]);
  res;
}
conj_ok(lamK) = { my(cl = subst(lift(lamK), 'x, 1/ZE)); \\ apply conj: zeta->1/zeta, alpha->beta
  \\ write lam = A(z) + alpha B(z) using coords, then conj = A(1/z) + beta B(1/z)
  my(v = coords(lamK), cv = sum(i = 1, N, v[i]*CBAS[i])); cv == -lamK; }

\\ lambda from a totally positive generator g of the K+ ideal I (I = q * extra): r = (-1)^k / (7 (c^2-4)^k g)
lamfrom(g) = {
  my(k = (ELL-3)/2, r = (-1)^k / (7 * (C1^2 - 4)^k * g));
  [r, (2*AL + 1) * toK(r, XK)];
}
\\ find totally positive generator of ideal I of K+; returns 0 if not narrowly principal
posgen(I) = {
  my(t = bnrisprincipal(BNR, I, 1));
  if (t[1] != 0*t[1], return(0));
  my(g = nfbasistoalg(NF, t[2]));
  if (idealhnf(NF, g) != idealhnf(NF, I), error("generator mismatch"));
  if (vecmin(nfeltsign(NF, t[2])) < 1, error("generator not totally positive"));
  Mod(lift(g), P);
}
qideal() = { my(L = idealprimedec(NF, ELL)); if (#L != 1, error("q")); L[1]; }
\\ v0 above p (self-glued): place of K+ where (s+7c)/2 vanishes
v0ideal(p) = { my(L = select(pr -> idealval(NF, lift((S1 + 7*C1)/2), pr) > 0, idealprimedec(NF, p))); if (#L != 1, error("v0 count")); L[1]; }
\\ self-glued residue data at p: (z0,a0) in F_{p^2} with pi = pibar mod V, and its im-functional
v0functional(p, v0) = {
  my(t = ffgen(p^2, 't), a0 = polrootsmod('a^2 + 'a + 2, t)[1]);  \\ a0 in F_{p^2}\F_p
  my(b0 = -1 - a0, w = b0/a0, z0 = w^((ELL+1)/2));
  if (z0^ELL != 1 || z0 == 1, error("z0")); if (a0*z0 != b0/z0, error("not self-glued"));
  \\ consistency with the K+ place v0: residues of c and s
  my(modpr = nfmodprinit(NF, v0), cres = nfmodpr(NF, nfalgtobasis(NF, C1), modpr), sres = nfmodpr(NF, nfalgtobasis(NF, S1), modpr));
  my(cz = z0 + 1/z0, sz = (2*a0+1)*(z0 - 1/z0));
  if (cz^p != cz || sz^p != sz, error("c,s not in F_p"));
  my(okc = (ff2int(cres) == ff2int(cz)) , oks = (ff2int(sres) == ff2int(sz)));
  if (!(okc && oks), error("residue mismatch with v0"));
  my(rv = vector(N, i, my(j = (i-1) % (ELL-1), e = (i-1) \ (ELL-1)); z0^j * a0^e));
  vector(N, i, lift(polcoeff(rv[i].pol, 1)));   \\ coefficient of t: residue in F_p iff 0
}
\\ split prime: kernel functional of P above split place v of K+ (degree 1 in K)
splitfunctional(p, v) = {
  my(modpr = nfmodprinit(NF, v), cres = Mod(ff2int(nfmodpr(NF, nfalgtobasis(NF, C1), modpr)), p), sres = Mod(ff2int(nfmodpr(NF, nfalgtobasis(NF, S1), modpr)), p));
  my(zs = polrootsmod('z^2 - lift(cres)*'z + 1, p), as = polrootsmod('a^2 + 'a + 2, p));
  for (i = 1, #zs, for (j = 1, #as, my(z0 = zs[i], a0 = as[j]);
    if ((2*a0+1)*(z0 - 1/z0) == sres,
      return(vector(N, k, my(jj = (k-1) % (ELL-1), e = (k-1) \ (ELL-1)); lift(z0^jj * a0^e))))));
  error("no residue data");
}
