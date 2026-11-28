default(parisizemax,4000000000);
p2=146505763881528721;
print("isprime p2 (APRCL) ", isprime(p2,2), " isprime 263 ", isprime(263,2));
print(factor(p2-1)); print(factor(p2+1));
print("p2 mod 131 = ", p2%131, "  kron(-7,p2)=", kronecker(-7,p2), " kron(-7,263)=",kronecker(-7,263)," 263%131=",263%131);
\\ Lucas U_n for x^2+x+2: U_0=0,U_1=1,U_{n+1}=-U_n-2U_{n-1}
U(n)=my(a=0,b=1);for(i=1,n,[a,b]=[b,-b-2*a]);a;
print("U_131 = ",U(131)," = ",factor(U(131)));
\\ independent h^-: odd characters of (Z/917)^* trivial on subgroup fixing K
\\ K = Q(zeta131,sqrt-7) ; chars chi = psi^a * lam^b, lam = (-7/.) = (./7)
hm(l)={my(g=znprimroot(l),f=7*l,P=1.,cnt=0, prec=0);
 for(a=0,l-2, for(b=0,1,
   par = (a%2) + b; if(par%2==0, next); \\ psi^a odd iff a odd; lam odd
   my(cond = if(a==0,1,l)*if(b,7,1), s=0.);
   for(x=1,cond, if(gcd(x,cond)>1,next);
     my(v=1.); if(a, v*=exp(2*Pi*I*a*znlog(x%l,g)/(l-1))); if(b, v*=kronecker(x,7));
     s+=v*x);
   P*= -s/cond/2; cnt++));
 [P*2*l, cnt]};
\\ h^- = Q*w*prod ; w=2l
default(realprecision,300);
foreach([11,17,19],l, r=hm(l); print(l," w*prod=",round(real(r[1])*2)/2," im ",abs(imag(r[1]))," #",r[2]));
r=hm(131); v=real(r[1]); print("131: 2*w*prod=",round(2*v), " err ", abs(2*v-round(2*v)), " imag ",abs(imag(r[1])));
h=round(2*v); print("h^-(K131) with Q=2: ",h, "  ok=", h==1128247071505617895582943540743356840380563482194685535482598575);
print("(p2+1)/262=",(p2+1)/262, "  prod=",h*(p2+1)/262);
\\ bnf check small ell
foreach([11,17],l, K=bnfinit(polcompositum(polcyclo(l),x^2+x+2)[1],1); print(l," h_K=",K.no," cyc ",K.cyc));
quit
