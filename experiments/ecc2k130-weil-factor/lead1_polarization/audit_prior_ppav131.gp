\\ usage: python3 dump_prior_ppav131.py && gp -q < audit_prior_ppav131.gp
\\ Independent audit (PARI, own arithmetic) of the prior producer's explicit ell=131 module
\\ (../construct_ppav.py, ../explicit_ppav.json). Checks done here from the raw matrices only:
\\  (1) the form on O_K is b(x,y) = Tr_{K/Q}(lam x conj(y)) for an explicit lam, recovered from row 1, verified on all entries;
\\  (2) conj(lam) = -lam;  (3) sign of Im phi(lam) at all 130 embeddings of each s-constant CM type;
\\  (4) the lattice I = O_K + Z v/DEG: Gram = PRINC, exact det(PRINC) = 1, pi- and pibar-stability, zeta-stability;
\\  (5) where End(I) is non-maximal: is the line Z/DEG * v stable under zeta mod 263, mod p2?
default(parisizemax, 4000000000)
default(threadsizemax, 4000000000)
\r prior_ppav131_dump.txt
ell = 131; d = ell - 1; n = 2*d; P2 = 146505763881528721;
PHI = polcyclo(ell, 'y); AP = 'x^2 + 'x + 2;
el(c) = Mod(Mod(sum(j = 0, d-1, c[j+1]*'y^j), PHI) + Mod(sum(j = 0, d-1, c[d+j+1]*'y^j), PHI)*'x, AP);
co(e) = { my(p = lift(e), c0 = lift(polcoef(p, 0, 'x)), c1 = lift(polcoef(p, 1, 'x)));
  concat(vector(d, j, polcoef(c0, j-1, 'y)), vector(d, j, polcoef(c1, j-1, 'y)))~; }
trQ(a) = { my(q = lift(a)); ell*polcoef(q, 0, 'y) - sum(j = 0, d-1, polcoef(q, j, 'y)); }
tr(e) = { my(p = lift(e)); trQ(2*Mod(polcoef(p, 0, 'x), PHI) - Mod(polcoef(p, 1, 'x), PHI)); }
cjQ(a) = Mod(subst(lift(a), 'y, Mod('y, PHI)^(ell-1)), PHI);
cj(e) = { my(p = lift(e)); Mod(cjQ(Mod(polcoef(p, 0, 'x), PHI)) + cjQ(Mod(polcoef(p, 1, 'x), PHI))*(-1 - 'x), AP); }
E = vector(n, k, el(vector(n, i, i == k)));
t0 = getabstime();
TR = matrix(n, n, k, j, tr(E[k]*cj(E[j])));           \\ TR[k,j] = Tr(e_k conj(e_j))
\\ (1) recover lam: OMEGA[1,j] = b(e_1, e_j) = Tr(lam conj e_j) = sum_k c_k TR[k,j]
c = matsolve(TR~, OMEGA[1,]~); lam = el(c);
C = matrix(n, n, i, k, co(lam*E[i])[k]);
ok1 = (C*TR == OMEGA);
printf("(1) OMEGA == [Tr(lam e_i conj e_j)] on all %d entries: %d   (denominator of lam coords: %d) [%d ms]\n", n^2, ok1, denominator(c), getabstime()-t0);
printf("(2) conj(lam) == -lam: %d\n", cj(lam) == -lam);
\\ (3) signs at the embeddings phi_k,s : zeta -> exp(2 pi i k/131), sqrt-7 -> s*i*sqrt7  (alpha = (-1+sqrt-7)/2)
default(realprecision, 3000);
sg(s) = { my(al = (-1 + s*I*sqrt(7))/2, p = lift(lam), out = vector(d));
  for (k = 1, d, my(zz = exp(2*Pi*I*k/ell), v = subst(lift(polcoef(p, 0, 'x)), 'y, zz) + subst(lift(polcoef(p, 1, 'x)), 'y, zz)*al);
    out[k] = sign(imag(v)); if (abs(imag(v)) < 10^-100, out[k] = 0)); out; }
sp = sg(1); sm = sg(-1);
{printf("(3) Phi_+ (sqrt-7 -> +i sqrt7): #Im>0 = %d, #Im<0 = %d, #~0 = %d ; Phi_-: #Im>0 = %d, #Im<0 = %d\n",
  #select(t -> t > 0, sp), #select(t -> t < 0, sp), #select(t -> t == 0, sp), #select(t -> t > 0, sm), #select(t -> t < 0, sm));}
printf("    mixed-s CM types are not the canonical-lift type (v_2(phi(pi)) depends only on phi|Q(sqrt-7))\n");
\\ (4) the lattice I
B = matid(n); B[, PIVOT + 1] = VNUM / DEG;
G = B~ * OMEGA * B;
printf("(4) Gram of b on I == PRINC (up to transpose convention): %d ; exact det(PRINC) = %d ; integral %d\n", G == PRINC || G~ == PRINC, matdet(PRINC), denominator(PRINC) == 1);
vI = el(VNUM / DEG); pie = el(concat(vector(d, j, j == 2), vector(d, j, 0))) * 0 + Mod(Mod('y, PHI)*'x, AP); pibe = 2 / pie;
inI(e) = denominator(matsolve(B, co(e))) == 1;
stab(m) = { for (k = 1, n, if (!inI(m*el(B[,k])), return(0))); 1; }
zet = Mod(Mod('y, PHI), AP);
printf("    I stable under pi: %d, under pibar=2/pi: %d, under zeta: %d\n", stab(pie), stab(pibe), stab(zet));
\\ (5) local position of the extra generator
vv = VNUM; zv = co(zet * el(vv));
prop(p) = { my(a = vv*Mod(1, p), b = zv*Mod(1, p)); matrank(matconcat([a, b])) == 1; }
printf("(5) zeta*v proportional to v mod 263: %d ; mod p2=%d: %d  (so End(I) is non-maximal exactly at p2)\n", prop(263), P2, prop(P2));
pv = co(pie * el(vv)); pbv = co(pibe * el(vv));
propx(w, p) = matrank(matconcat([vv*Mod(1, p), w*Mod(1, p)])) == 1;
printf("    pi*v ~ v and pibar*v ~ v mod p2: %d %d (line is the self-glued residue point: pi == pibar in F_p)\n", propx(pv, P2), propx(pbv, P2));
