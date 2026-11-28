\\ independent h^- for K=Q(zeta_l, sqrt(-7)) and control h^-(Q(zeta_l)) 
default(realprecision,500);
hm(l, withq7)={my(g=znprimroot(l), lg=vector(l-1), prodv=1., f, S, n=0, r);
 for(k=0,l-2, lg[lift(g^k)]=k);
 for(j=0,l-2,
   \\ psi_j odd iff j odd
   if(j%2==1,
     S=sum(a=1,l-1, exp(2*Pi*I*j*lg[a]/(l-1))*a)/l; prodv*= -S/2; n++);
   if(withq7 && j%2==0,
     if(j==0, S=-1, f=7*l; S=sum(a=1,f, if(a%7 && a%l, kronecker(a,7)*exp(2*Pi*I*j*lg[a%l]/(l-1))*a, 0))/f);
     prodv*= -S/2; n++));
 [n, prodv]};
{
\\ controls: h^-(Q(zeta_p)) = Q w prod, Q=1, w=2p
foreach([23,29,31,37,41,43,47,59,67,71,131], p, my(r=hm(p,0)); print("Q(zeta_",p,") h^- = ", round(real(r[2])*2*p,&e), "  err e=",e, " imag ", exponent(imag(r[2]))));
foreach([3,5,11,13,17,19,131], p, my(r=hm(p,1), v=real(r[2])*2*p);
  print("K_",p," #oddchars=",r[1]," Q=1 value=", v, "\n    Q=2 value=", round(2*v,&e), " err=",e, " frac(Q=1) ", v-round(v)));
}
quit
