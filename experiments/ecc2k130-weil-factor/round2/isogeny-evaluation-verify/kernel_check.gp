\\ Independent PARI recomputation of the pi action on I/O_K for the PR #57 module.
\\ Element = [A(z), B(z)] meaning A + alpha*B, alpha^2 = -alpha - 2, z = zeta_131.
default(parisizemax, 4000000000);
read("v.gp");
p1=263; p2=146505763881528721; if(p1*p2!=DEG, error("deg"));
Phi=polcyclo(131,'z);
toel(W)={[Pol(Vecrev(W[1..130]),'z), Pol(Vecrev(W[131..260]),'z)]};
mulel(x,y,p)={my(ac=x[1]*y[1], bd=x[2]*y[2], ad=x[1]*y[2]+x[2]*y[1]);
  [ ((ac-2*bd)*Mod(1,p))%Phi, ((ad-bd)*Mod(1,p))%Phi ]};
coeffs(x,p)={concat(Vecrev(lift(x[1]*Mod(1,p)),130), Vecrev(lift(x[2]*Mod(1,p)),130))};
propc(w,u,p)={my(W=coeffs(w,p),U=coeffs(u,p),k=1); while(U[k]==0,k++); my(c=Mod(W[k],p)/U[k]); if(Mod(W,p)==c*Mod(U,p), lift(c), -1)};
v=toel(V);
PI=[0,'z]; PIB=[-'z^130,-'z^130]; ALPHA=[0,1]; ZETA=['z,0];
tt=mulel(PI,PIB,10^9); print("pi*pibar = ", lift(tt[1]), " + alpha*(", lift(tt[2]),")");
{for(i=1,2, my(p=[p1,p2][i]);
  my(c=propc(mulel(PI,v,p),v,p), cb=propc(mulel(PIB,v,p),v,p), ca=propc(mulel(ALPHA,v,p),v,p), cz=propc(mulel(ZETA,v,p),v,p));
  printf("p=%s: pi->%s pibar->%s alpha->%s zeta->%s (-1 = not proportional); JS mod p=%s\n", p, c, cb, ca, cz, JS%p);
  if(c>=0, printf("   c^2-2 mod p = %s, ord(c)=%s\n", (c^2-2)%p, znorder(Mod(c,p)))));}
printf("gcd(content(V),DEG)=%s\n", gcd(content(V),DEG));
w=V; w[2]+=1; wv=toel(w);
printf("negative control (perturbed v) pi-proportional mod p2: %s\n", propc(mulel(PI,wv,p2),wv,p2)>=0);
printf("lcm of orders = %s\n", lcm(znorder(Mod(110,263)), znorder(Mod(14975904886999252,p2))));
quit;
