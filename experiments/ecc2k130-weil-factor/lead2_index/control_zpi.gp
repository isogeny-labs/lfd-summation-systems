\\ PARI side of control_zpi.py: [O_F (x) Z[zeta] : Z[pi]] = |disc(charpoly pi)| / |disc(Lambda)| ^ (1/2)  via HNF.
coords(e, ell) = {my(d=ell-1, v=vector(2*d), p=lift(lift(e)));
  for(k=0,1, my(c=polcoef(p,k,'w)); for(i=0,d-1, v[k*d+i+1]=polcoef(c,i,'z))); v~};
{foreach([5,11,13,19,31], ell, my(d=ell-1, Pz=polcyclo(ell,'z), Pw='w^2-'w+2, z=Mod(Mod('z,Pz),Pw), w=Mod(Mod(1,Pz)*'w,Pw));
  my(pi=(w-1)*z, x=z^0, M=matrix(2*d,2*d)); for(i=1,2*d, M[,i]=coords(x,ell); x*=pi);
  print(ell, " ", abs(matdet(M))));}
quit
