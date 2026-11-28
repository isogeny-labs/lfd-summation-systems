default(parisizemax,2000000000);
\\ U_ell = (a^ell-b^ell)/(a-b), a,b roots of x^2+x+2  (Lucas U with P=-1,Q=2)
U(n)={my(u0=0,u1=1,t);for(i=2,n,t=-u1-2*u0;u0=u1;u1=t);if(n==0,0,u1)};
\\ second route: resultant / norm in Q(sqrt-7)
U2(n)={my(K=nfinit(y^2+y+2),a=Mod(y,y^2+y+2),b=-1-a); lift((a^n-b^n)/(a-b))};
l=131; u=U(l); print("U131=",u," check2=",u==U2(l));
f=factor(u); print(f);
for(i=1,#f~,p=f[i,1]; print("p=",p," isprime(proof)=",isprime(p,2)," (-7/p)=",kronecker(-7,p)," p mod 131=",p%131, " ord(p) mod 131=",znorder(Mod(p,131))));
print("eps=(-7/131)=",kronecker(-7,131));
\\ local structure over F_p: places w of K=Q(zeta,sqrt-7) above p, images of pi, pibar
loc(p,l)={my(T=ffinit(p,4*l),F,res=[]);
  \\ factor Phi_l(x) over F_p, then y^2+y+2 over each residue field
  F=factormod(polcyclo(l),p)[,1]; my(cnt=0,info=[]);
  for(i=1,#F, my(g=F[i],d=poldegree(g),q=ffgen(ffinit(p,2*d),'t),z=polrootsmod(lift(g),p));
    0);
  F};
