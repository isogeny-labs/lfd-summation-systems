default(parisizemax,4000000000);
mk(l)=my(R=polresultant(polcyclo(l,y), x^2+x*y+2*y^2, y)); R/pollead(R);
realpart(l)={my(P=mk(l),g=l-1,h=0,Q=P); forstep(k=g,0,-1, my(c=polcoeff(Q,g+k)); h+=c*t^k; Q-=c*x^(g-k)*(x^2+2)^k); if(Q!=0,error("bad")); h};
U=263*146505763881528721;
{
for(i=1,1, l=131; P=mk(l); h=realpart(l);
 dP=abs(poldisc(P)); dK=131^258*7^130; q=dP/dK; print("[O_K:Z[pi]]^2 integral? ",denominator(q)==1, "  issquare ",issquare(q));
 dh=abs(poldisc(h)); dKp=131^129*7^65; qp=dh/dKp; print("[O_K+:Z[t]]^2 = U^128 ? ", qp==U^128);
 Nt=abs(polresultant(h,t^2-8)); dR0=dh^2*Nt; q0=dR0/dK; print("[O_K:Z[pi,pibar]]^2 == U^258 ? ", q0==U^258);
 print("N(t^2-8) = ", factor(Nt));
)}
\\ control ell=11, U=23: expect [O_K:R0]=23^9
{ l=11; P=mk(l); h=realpart(l); dh=abs(poldisc(h)); Nt=abs(polresultant(h,t^2-8));
  dK=11^18*7^10; print("ell=11 [O_K:R0]^2 = ", factor((dh^2*Nt)/dK)); }
quit
