\\ Lead 3 (Deligne modules), METHOD A: pure F_p linear algebra, no character theory.
\\ O_K = Z[z]/Phi_ell (x) Z[w]/(w^2-w+2), basis z^i w^j, coordinate index i + d*j (d=ell-1).
\\ alpha=w-1, beta=-w, pi=alpha*z, pibar=beta*z^-1.  R = Z[pi,pibar].
\\ For a prime p:  B = O_K/pO_K,  S = image of R in B.
\\  * index-p R-stable SUB-lattices T (pO_K c T c O_K) <-> lines of f in B^* with f*Mpi = l f, f*Mpb = (2/l) f
\\  * index-p R-stable SUPER-lattices (O_K c T' c p^-1 O_K) <-> lines of v in B with Mpi v = l v, Mpb v = (2/l) v
\\  O_K-stable ones: additionally eigen for Mz (then for w too, since alpha = pi z^-1 and w = alpha+1).
\\ Counts are exact: sum over eigenvalues l in F_p of (p^dim - 1)/(p-1).
default(parisizemax, 4000000000);
mats(ell) = {my(d=ell-1, C=matcompanion(polcyclo(ell)), Z=matrix(2*d,2*d), W=matrix(2*d,2*d), I=matid(d));
  for(a=1,d, for(b=1,d, Z[a,b]=C[a,b]; Z[d+a,d+b]=C[a,b]));
  for(a=1,d, W[d+a,a]=1; W[a,d+a]=-2; W[d+a,d+a]=1);   \\ w*(x0 + x1 w) = -2 x1 + (x0+x1) w
  my(A=W-matid(2*d), Bm=-W, Mpi=A*Z, Mpb=Bm*Z^-1);
  [Z, W, Mpi, Mpb]};
\\ rank of R mod p (dim_Fp of S = (R+pO_K)/pO_K)
rankS(ell, p, M) = {my(d=ell-1, n=2*d, Mpi=M[3]*Mod(1,p), Mpb=M[4]*Mod(1,p), e=vectorv(n,k,k==1)*Mod(1,p), cols=List(), x=e);
  for(i=0,n-1, listput(cols,x); x=Mpi*x); x=Mpb*e; for(j=1,n-1, listput(cols,x); x=Mpb*x);
  matrank(matconcat(Vec(cols)))};
count(ell, p, M, side) = {my(d=ell-1, n=2*d, Mpi=M[3]*Mod(1,p), Mpb=M[4]*Mod(1,p), Mz=M[1]*Mod(1,p));
  if(side=="sub", Mpi=Mpi~; Mpb=Mpb~; Mz=Mz~);   \\ left eigenvectors = right eigenvectors of transposes
  my(ls=polrootsmod(charpoly(M[3]), p), zs=polrootsmod(polcyclo(ell), p), I=matid(n), tot=0, okst=0, res=List());
  foreach(ls, l, my(K=matker(matconcat([Mpi-l*I; Mpb-(2/l)*I])), dl=#K, ok=0);
     if(dl==0, next);
     \\ W is z-stable (z commutes with pi, pibar); O_K-stable lines in W = z-eigenlines of z|W
     my(Kz=Mz*K);
     my(Zw2=matinverseimage(K, Kz));
     if(dl==1, ok=1, my(cz=charpoly(Zw2), rz=polrootsmod(lift(cz), p)); ok=sum(j=1,#rz, my(e=#matker(Zw2-rz[j]*matid(dl))); (p^e-1)/(p-1)));
     tot += (p^dl-1)/(p-1); okst += ok; listput(res, [lift(l), dl, ok]));
  [tot, okst, tot-okst, Vec(res)]};
\\ dim_Fp of End(T)/pO_K for the sub-lattice T = ker f (+pO_K):  {x : x.f in F_p f}
enddim(ell, p, M, f) = {my(d=ell-1, n=2*d, Mz=M[1]*Mod(1,p), Mw=M[2]*Mod(1,p), cols=vector(n), g=f);
  for(i=0,d-1, cols[i+1]=g~; cols[d+i+1]=(g*Mw)~; g=g*Mz);
  n - (matrank(matconcat([matconcat(cols), f~]))-1)};
report(ell, p) = {my(M=mats(ell), t=getabstime(), s=rankS(ell,p,M), cs=count(ell,p,M,"sub"), cp=count(ell,p,M,"super"));
  printf("ell=%d p=%d  dim S=%d (codim %d)  #idx-p R-stable SUB-lattices=%s  O_K-stable=%s  NOT O_K-stable=%s | SUPER: total=%s OK-stable=%s NOT=%s\n",
     ell, p, s, 2*(ell-1)-s, cs[1], cs[2], cs[3], cp[1], cp[2], cp[3]);
  my(v=cs[4], h=Map()); foreach(v, r, my(key=[r[2],r[3]]); mapput(h,key, if(mapisdefined(h,key),mapget(h,key),0)+1));
  printf("   eigenvalue-type histogram [dim W_l, #O_K-stable lines] -> count: %s\n", Mat(h));
  \\ End(T) for one non-O_K-stable sub-lattice: take any f in a W_l of dim>=2 which is not z-eigen
  my(Mpi=M[3]~*Mod(1,p), Mpb=M[4]~*Mod(1,p), Mz=M[1]~*Mod(1,p), n=2*(ell-1), I=matid(n), found=0);
  foreach(v, r, if(r[2]>=2 && !found, my(l=Mod(r[1],p), K=matker(matconcat([Mpi-l*I; Mpb-(2/l)*I])));
     my(f=(K[,1]+K[,2])~);  \\ generic element; check not z-eigen
     my(fz=f*M[1]*Mod(1,p)); if(matrank(matconcat([f~, fz~]))==2,
        my(ed=enddim(ell,p,M,f)); my(zin=vector(ell-1,k, my(g=f*(M[1]^k)*Mod(1,p)); matrank(matconcat([f~,g~]))==1));
        printf("   sample sub-lattice (eigenvalue pi->%d): dim End(T)/pO_K = %d  => [O_K:End(T)] = p^%d ; zeta^k in End(T) for k in %s\n",
          r[1], ed, n-ed, select(x->x, zin, 1));
        found=1)));
  printf("   (%.1fs)\n", (getabstime()-t)/1000.)};
\\ controls: primes q not dividing U_ell (R_q = O_K,q) must give NOT=0; ell=5,13 have |U|=1.
report(5, 11); report(5, 31); report(13, 53); report(13, 79);
report(11, 29); report(11, 67);          \\ negative controls at ell=11 (q not | 23)
report(11, 23);                          \\ positive: U_11 = 23
report(19, 191); report(19, 457);        \\ 191 control, 457 | U_19
report(131, 787);
report(131, 263);
report(131, 146505763881528721);

