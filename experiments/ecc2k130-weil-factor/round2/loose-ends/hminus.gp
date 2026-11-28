\\ Relative class number h^- of CM abelian fields via h^-/Q = w * prod_{chi odd} (-B_{1,chi}/2)  (Washington Thm 4.17).
\\ Unconditional (no GRH).  Fields: Q(zeta_l) and K_l = Q(zeta_l, sqrt(-7)).
default(realprecision, 600);
B1(chi, f) = sum(a = 1, f, chi(a) * a) / f;   \\ chi primitive mod f, chi != 1
\\ psi_j mod l: psi_j(g^k) = e(jk/(l-1)), g primitive root; eps = Legendre (./7) (odd, conductor 7)
hm(l, withsqrt7) = {
  my(g = znprimroot(l), n = l-1, logt = vector(l-1), P = 1., cnt = 0);
  for (k = 0, n-1, logt[lift(g^k)] = k);
  my(E = exp(2*Pi*I/n));
  for (j = 0, n-1, for (e = 0, withsqrt7,
     if ((j + e) % 2 == 0, next);          \\ odd iff psi(-1)eps(-1)^e = (-1)^(j+e) = -1
     my(f = if (j == 0, 7, l) * if (e && j, 7, 1));
     my(chi = a -> if (gcd(a, f) > 1, 0, if (j, E^(j*logt[a % l]), 1) * if (e, kronecker(a, 7), 1)));
     P *= -B1(chi, f)/2; cnt++));
  my(w = if (withsqrt7, 2*l, 2*l));  \\ mu(K) = <-zeta_l> in both cases (l != 7, sqrt(-7)*... gives no extra roots of unity)
  my(v = w * P);
  if (abs(imag(v)) > 1e-50, error("not real"));
  \\ return 2*w*prod (an integer whenever h^-/Q is an integer or half-integer)
  my(r = round(2*real(v)), err = abs(2*real(v) - r));
  [cnt, r, err]
};
\\ controls: Q(zeta_p), Q=1: known h^- (Washington, Tables Sec. 4): 23:3, 31:9, 37:37, 41:121, 43:211, 47:695, 59:41241, 67:853513, 71:3882809
known = [[23,3],[31,9],[37,37],[41,121],[43,211],[47,695],[59,41241],[67,853513],[71,3882809]];
for (i = 1, #known, my(r = hm(known[i][1], 0)); printf("Q(zeta_%d): h^- = %d  (known %d)  match=%d  roundoff=%.3g\n", known[i][1], r[2]/2, known[i][2], r[2]==2*known[i][2], r[3]));
\\ Q(sqrt-7): h = 1 (w=2): direct
printf("Q(sqrt-7): h = %s (expect 1)\n", 2 * (-B1(a -> kronecker(a,7), 7)/2));
foreach([11, 17, 19], l, my(r = hm(l, 1)); printf("K_%d = Q(zeta_%d, sqrt-7): 2*h^-/Q = %d  roundoff=%.3g  (#odd chars %d)  [half-integer h^-/Q forces Q=2]\n", l, l, r[2], r[3], r[1]));
r = hm(131, 0); printf("Q(zeta_131): h^- = %d (%d digits) roundoff=%.3g\n  factor: %s\n", r[2]/2, #Str(r[2]/2), r[3], factor(r[2]/2, 10^7));
r = hm(131, 1); printf("K_131: 2*h^-/Q = %d  odd? %d  roundoff=%.3g  #odd chars=%d\n", r[2], r[2]%2, r[3], r[1]);
hm131 = if (r[2] % 2, r[2], r[2]/2);  \\ if 2h^-/Q is odd then Q = 2 and h^- = 2h^-/Q
printf("K_131: Q = %d, h^-(K_131) = %d (%d digits)\n  factor (trial to 1e7): %s\n  isprime(cofactor)? %s\n", if (r[2]%2, 2, 1), hm131, #Str(hm131), factor(hm131, 10^7), my(F=factor(hm131,10^7)); [F[#F~,1], ispseudoprime(F[#F~,1])]);
printf("h^-(Q(zeta_131)) divides h^-(K_131)? %d\n", hm131 % 28496379729272136525 == 0);
