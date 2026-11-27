default(parisizemax, 4000000000);
default(threadsizemax, 2000000000);
\\ Is zeta_ell in Z[pi, pibar] for the Weil factor A_ell of Res_{F_2^ell/F_2}(E), E: y^2+xy=x^3+1?
\\ O_K = Z[z]/(Phi_ell) (x) Z[w]/(w^2-w+2), w=(1+sqrt-7)/2. alpha=w-1, beta=-w. pi=alpha*z, pibar=beta*z^-1.
coords(e, ell) = {my(d=ell-1, v=vector(2*d), p=lift(lift(e)));
  for(k=0,1, my(c=polcoef(p,k,'w)); for(i=0,d-1, v[k*d+i+1]=polcoef(c,i,'z))); v~};
run(ell) = {my(d=ell-1, Pz=polcyclo(ell,'z), Pw='w^2-'w+2, z=Mod(Mod('z,Pz),Pw), w=Mod(Mod(1,Pz)*'w,Pw));
  my(pi=(w-1)*z, pb=(-w)*z^(ell-1), cols=List());
  my(x=z^0); for(i=0,2*d-1, listput(cols,coords(x,ell)); x*=pi);
  x=pb; for(j=1,2*d-1, listput(cols,coords(x,ell)); x*=pb);
  my(M=matconcat(Vec(cols)), H=mathnf(M), idx=abs(prod(k=1,matsize(H)[1],H[k,k])));
  my(sol=matsolve(H, coords(z,ell)), inZ=(denominator(sol)==1));
  my(sola=matsolve(H, coords(w,ell)));
  printf("ell=%d rank=%d index[O_K:Z[pi,pibar]]=%s factored=%s  zeta in order? %d  omega in order? %d  denom(zeta coords)=%s\n",
    ell, matsize(H)[2], idx, factor(idx), inZ, denominator(sola)==1, denominator(sol));};
foreach([5,11,13,19,131], l, run(l));
