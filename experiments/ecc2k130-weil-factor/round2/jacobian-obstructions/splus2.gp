\\ |S+(O_K)| = [U+_{K+} : N_{K/K+} U_K] for K_ell (GRH, bnfinit): the polarization-unit factor 2^u.
default(parisizemax, 4000000000);
for(t = 1, 3, ell = [11, 17, 19][t]; \
  pK = polredbest(polcompositum(polcyclo(ell), 'x^2+7)[1]); S = nfsubfields(pK, ell-1); \
  for(j = 1, #S, if(polsturm(S[j][1]) == ell-1, pKp = polredbest(S[j][1]))); \
  pKp = subst(pKp, 'x, 'y); bKp = bnfinit(pKp, 1); n = ell-1; nfp = nfinit(pKp); rnf = rnfinit(nfp, 'x^2+7); \
  bK = bnfinit(rnfequation(nfp, 'x^2+7), 1); \
  Up = concat([-1], bKp.fu); \
  Mtp = matrix(n, #Up, k, j, Mod(nfeltsign(bKp, Up[j])[k] < 0, 2)); \
  kerS = matker(Mtp); Ku = concat([bK.tu[2]], bK.fu); \
  Mi = matrix(#Up, #Ku, a, b, 0); \
  for(b = 1, #Ku, nu = rnfeltnorm(rnf, rnfeltabstorel(rnf, Ku[b])); ex = bnfisunit(bKp, nu); \
     Mi[1, b] = Mod(ex[#ex], 2); for(k = 1, #ex-1, Mi[k+1, b] = Mod(ex[k], 2))); \
  printf("ell=%d: h(K+)=%d h(K)=%d; dim U/U^2=%d, dim U+/U^2=%d, dim N(U_K)/U^2=%d, norms tot.pos=%d => log2|S+(O_K)|=%d\n", \
     ell, bKp.no, bK.no, #Up, #kerS, matrank(Mi), matrank(Mtp*Mi)==0, #kerS - matrank(Mi)));
