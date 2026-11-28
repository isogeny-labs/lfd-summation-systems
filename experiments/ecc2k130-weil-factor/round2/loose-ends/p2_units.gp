\\ Check of the unit-image argument at ell=131, p = p2, in F_{p^2} (independent of bnfinit; no GRH).
\\ P = the prime of K above p2 with pi mod P in F_p: alpha -> a, zeta -> w with w^2 = b/a.
p = 146505763881528721; l = 131;
default(parisizemax, 4000000000);
if (kronecker(-7, p) != -1, error("x^2+x+2 not irreducible mod p"));
a = ffgen(Mod(1, p)*(y^2+y+2), 'a); b = -1 - a;   \\ F_{p^2} = F_p(alpha) since (-7/p2) = -1
printf("a in F_p? %d (expect 0: (-7/p2) = %d)\n", a^p == a, kronecker(-7, p));
w = (b/a)^((l+1)/2);
printf("w^131 = 1: %d, w != 1: %d, w^2 = b/a: %d\n", w^l == 1, w != 1, w^2 == b/a);
lam = a*w; printf("pi mod P = a*w in F_p: %d ; (a w)^2 = 2: %d ; pibar = b/w in F_p: %d\n", lam^p == lam, lam^2 == 2, (b/w)^p == b/w);
printf("Frobenius at P = complex conjugation: w^p = w^-1: %d, a^p = b: %d\n", w^p == w^-1, a^p == b);
\\ images in F_{p^2}^*/F_p^* via u -> u^(p-1)
img(u) = fforder(u^(p-1));
printf("zeta: order %d (expect 131); -zeta: order %d\n", img(w), img(-w));
\\ real cyclotomic units xi_k = zeta^((1-k)/2) (1-zeta^k)/(1-zeta) (k=2..65): must map to 1
bad = 0; for (k = 2, 65, my(e = lift(Mod((1-k)/2, l))); my(xi = w^e*(1-w^k)/(1-w)); if (img(xi) != 1, bad++));
printf("real cyclotomic units xi_k, k=2..65, with nontrivial image: %d (expect 0)\n", bad);
\\ negative control: the non-real cyclotomic units (1-zeta^k)/(1-zeta) do map nontrivially
nn = 0; for (k = 2, 65, if (img((1-w^k)/(1-w)) != 1, nn++)); printf("non-real units (1-zeta^k)/(1-zeta) with nontrivial image: %d of 64 (expect 64)\n", nn);
\\ fundamental unit of Q(sqrt(917)), sqrt(917) = sqrt(-7)*sqrt(-131) -> (2a+1)*G, G = Gauss sum, G^2 = -131
G = sum(k = 1, l-1, kronecker(k, l) * w^k); printf("G^2 = -131: %d\n", G^2 == -131);
s917 = (2*a+1)*G; printf("s917^2 = 917: %d\n", s917^2 == 917);
q = quadunit(917); u = real(q) + imag(q)*(1 + s917)/2;   \\ disc 917 = 1 mod 4: w = (1+sqrt917)/2
printf("fundamental unit of Q(sqrt917): norm %d, image order %d (expect 1)\n", norm(q), img(u));
printf("(p2+1)/262 = %d = %s\n", (p+1)/262, factor((p+1)/262));
