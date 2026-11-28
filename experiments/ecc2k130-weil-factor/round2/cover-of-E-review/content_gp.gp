\\ cover-of-E review: gluing-degree constraint and ray-class numbers (PARI/GP)
default(parisizemax,4000000000);
\\ E: Frobenius char poly x^2+x+2; alpha = Mod(x, x^2+x+2)
al = Mod(x, x^2+x+2);
tr(l) = trace(al^l);           \\ s_l = alpha^l + beta^l
PA(l) = { my(PWl = x^(2*l) - tr(l)*x^l + 2^l, [q,r] = divrem(PWl, x^2+x+2)); if(r!=0, error("not divisible")); q };
content_at_alpha(P) = { my(v = lift(subst(P, x, al))); gcd(polcoeff(v,0), polcoeff(v,1)) };
U(n) = { my(a=0,b=1,t); for(i=2,n, t = -b - 2*a; a=b; b=t); if(n==0,0,b) };
for(l=3, 37, if(isprime(l) && l!=7, my(P=PA(l), c=content_at_alpha(P)); printf("l=%3d content(P_A(alpha))=%s   l*|U_l|=%s  match=%d  trace(A_l)=%d\n", l, c, l*abs(U(l)), c==l*abs(U(l)), -polcoeff(P, poldegree(P)-1))));
P = PA(131); c = content_at_alpha(P);
print("l=131: content(P_A(alpha)) = ", factor(c));
print("        131*|U_131|       = ", factor(131*abs(U(131))), "  equal: ", c==131*abs(U(131)));
print("trace(pi | A_131) = ", -polcoeff(P,259), "   => #C(F_2) for Jac(C)~ExA = 2+1-(-1)-trace = ", 2+1+1-(-polcoeff(P,259)));
print("P_A in Z[x^m] for some m|130, m>1? ", sum(i=1,#divisors(130)-1, my(m=divisors(130)[i+1]); prod(j=0,260, if(j%m!=0 && polcoeff(P,j)!=0, 0, 1))));
print("ord_131(2) = ", znorder(Mod(2,131)));
print("kronecker(-7,131) = ", kronecker(-7,131), "   kronecker(-7,263) = ", kronecker(-7,263), "  263 mod 131 = ", 263%131);
NE(m) = 2^m + 1 - tr(m);
print("#E(F_2^m) mod 131 for m | 130: ", vector(#divisors(130), i, [divisors(130)[i], NE(divisors(130)[i])%131]));
print("#E(F_2)=", NE(1), "  #E(F_4)=", NE(2));
print("2^d-1 mod 131, d=1..4: ", vector(4,d,(2^d-1)%131), "   4^d-1 mod 131, d=1..4: ", vector(4,d,(4^d-1)%131));
\\ alpha mod 131 in F_{131^2}: order
a131 = Mod(Mod(1,131)*x, Mod(1,131)*(x^2+x+2)); o = 1; t = a131; while(t != 1, t*=a131; o++); print("ord(alpha mod 131) in F_131^2^* = ", o);
quit
