\\ Lead 3, METHOD A cross-checks (PARI, global lattices; independent of deligne_chars.py / group_structure.py
\\ except that the pi-eigenvalue lambda of the chosen merged pair is passed in and RE-VERIFIED here by dim W_lambda = 2).
\\ For each (ell, p, lambda): build f in W_lambda (left common eigenvector of pi, pibar mod p, not a zeta-eigenvector),
\\   T  = {x in O_K : f.x = 0 mod p}  (sub-lattice, index p)   and   T' = O_K + Z v/p (v in right eigenspace; super-lattice)
\\ and compute: End(T) mod p, conductor rank (Gorenstein test), conj-stability of End(T), zeta^k membership,
\\   p-part of T/(pi^k-1)T, T'/(pi^k-1)T', O_K/(pi^k-1) by Smith normal form (full SNF for ell<=19, HNF mod p^6 for 131),
\\   v_p(#A(F_2^k)) from the resultant Res(P_A, x^k-1) where P_A = charpoly(pi).
\\ ell=11,19 also: GLOBAL Gorenstein test (trace dual O^v of O=End(T) is invertible: O^v (O:O^v) = O).
default(parisizemax, 4000000000);
default(threadsizemax, 2000000000);
hdet(H) = abs(prod(i=1,matsize(H)[1],H[i,i]));
mats(ell) = {my(d=ell-1, C=matcompanion(polcyclo(ell)), Z=matrix(2*d,2*d), W=matrix(2*d,2*d));
  for(a=1,d, for(b=1,d, Z[a,b]=C[a,b]; Z[d+a,d+b]=C[a,b]));
  for(a=1,d, W[d+a,a]=1; W[a,d+a]=-2; W[d+a,d+a]=1);
  my(Mpi=(W-matid(2*d))*Z, Mpb=(-W)*Z^-1);
  my(n=2*d, Mb=vector(n), Zi=Z^-1, Cj=matrix(n,n), e0=vectorv(n,k,k==1), ew=vectorv(n,k,k==d+1));
  my(Zp=matid(n), ci=e0, cw=e0-ew);
  for(i=0,d-1, Mb[i+1]=Zp; Mb[d+i+1]=Zp*W; Cj[,i+1]=ci; Cj[,d+i+1]=cw; Zp=Zp*Z; ci=Zi*ci; cw=Zi*cw);
  [Z, W, Mpi, Mpb, Mb, Cj]};
mulmat(M, x) = sum(k=1, #x, x[k]*M[5][k]);
setupT(ell, p, lam, M) = {my(n=2*(ell-1), I=matid(n), l=Mod(lam,p));
  my(K=matker(matconcat([(M[3]~-l*I); (M[4]~-(2/l)*I)])*Mod(1,p)), KR=matker(matconcat([(M[3]-l*I); (M[4]-(2/l)*I)])*Mod(1,p)));
  my(f, v, okf=0); if(#K!=2 || #KR!=2, error("dim W_lambda != 2"));
  for(c=1,5, f=(K[,1]+c*K[,2])~; if(matrank(matconcat([f~, (f*M[1])~]))==2, okf=1; break)); if(!okf, error("no f"));
  okf=0; for(c=1,5, v=KR[,1]+c*KR[,2]; if(matrank(matconcat([v, M[1]*v]))==2, okf=1; break)); if(!okf, error("no v"));
  my(BT=mathnf(matconcat([lift(matker(Mat(f))), p*I])));
  my(BTs=mathnf(matconcat([p*I, lift(v)]))/p);
  [f, v, BT, BTs]};
enddata(ell, p, M, f) = {my(n=2*(ell-1), N=matrix(n,n));
  for(k=1,n, N[,k]=(f*M[5][k]*Mod(1,p))~);
  my(KE=matker(matconcat([N, f~]))); my(E=matrix(n,#KE,i,j,KE[i,j]));  \\ x-part
  my(Eb=matimage(E), g=matker(Eb~));   \\ g: functional with ker g = E
  if(#g!=1, error("codim E != 1")); g=g[,1]~;
  my(Fr=matconcat(vector(n,k, g*M[5][k]*Mod(1,p))~)); my(r=matrank(Fr));   \\ rows g*M_{b_k}; f/pO_K = ker
  my(conjstab = matrank(matconcat([g~, (g*M[6]*Mod(1,p))~]))==1);
  my(zk=vector(2*ell, k, g*(M[1]^k)[,1]*Mod(1,p)==0));
  [#Eb, r, conjstab, select(x->x, zk, 1), g]};
\\ p-part exponents of the cokernel of an integer matrix C
ppart(C, p, full) = {my(D = if(full, matsnf(C), matsnf(mathnfmodid(C, p^6)))); my(v=List()); foreach(D, x, if(x && valuation(x,p)>0, listput(v, valuation(x,p)))); vecsort(Vec(v))};
pp_str(v,p) = if(#v==0, "0", my(s=""); foreach(v, e, s=concat(s, Str(if(s=="","","+"), "Z/", p, if(e>1, Str("^",e), "")))); s);
gorglobal(ell, p, M, g) = {my(n=2*(ell-1), I=matid(n));
  my(BO=mathnf(matconcat([lift(matker(Mat(g))), p*I])));
  my(G=matrix(n,n,i,j, trace(M[5][i]*M[5][j])));
  my(D=BO*(BO~*G*BO)^-1);   \\ trace dual of O, coordinates in O_K basis
  \\ (O : D) = {x : x D in O} = dual (std pairing) of the Z-span of rows of A = [BO^-1 * M_{D_j}] (as functionals of x)
  my(BOi=BO^-1, rows=List());
  for(j=1,n, my(Mj=mulmat(M, D[,j])); my(Aj=BOi*Mj); for(i=1,n, listput(rows, Aj[i,])));
  my(A=matconcat(Vec(rows)~)); my(den=denominator(A)); my(H=mathnf((den*A)~)/den); my(J=(H^-1)~);
  \\ product D * J
  my(prods=List()); for(i=1,n, my(Mi=mulmat(M, D[,i])); for(j=1,n, listput(prods, Mi*J[,j])));
  my(P=matconcat(Vec(prods))); my(dn=denominator(P)); my(HP=mathnf(dn*P)/dn);
  [HP==mathnf(BO), abs(matdet(BO))]};
run(ell, p, lam, tag, kl, full) = {my(M=mats(ell), t=getabstime(), n=2*(ell-1), S=setupT(ell,p,lam,M), f=S[1], v=S[2], BT=S[3], BTs=S[4]);
  printf("== ell=%d p=%d pair=%s lambda=%d :  dim W_lambda=2 verified; [O_K:T]=%d  [T':O_K]=%d\n", ell, p, tag, lam, hdet(BT), 1/hdet(BTs));

  printf("   pi,pibar preserve T: %d %d ; zeta preserves T: %d ; pi,pibar preserve T': %d %d ; zeta preserves T': %d\n",
     denominator(BT^-1*M[3]*BT)==1, denominator(BT^-1*M[4]*BT)==1, denominator(BT^-1*M[1]*BT)==1,
     denominator(BTs^-1*M[3]*BTs)==1, denominator(BTs^-1*M[4]*BTs)==1, denominator(BTs^-1*M[1]*BTs)==1);
  my(ed=enddata(ell,p,M,f));
  printf("   End(T): dim E/pO_K=%d -> [O_K:End(T)]=p^%d ; conductor: [O_K:f]=p^%d, [O:f]=p^%d -> Gorenstein (len(O_K/f)=2 len(O/f)): %d ; End(T) conj-stable: %d ; k<=2ell with zeta^k in End(T): %s\n",
     ed[1], n-ed[1], ed[2], ed[2]-1, ed[2]==2*(ed[2]-1), ed[3], ed[4]);

  if(full, my(gg=gorglobal(ell,p,M,lift(ed[5]))); printf("   GLOBAL Gorenstein test: trace dual of End(T) invertible: %d  ([O_K:End(T)]=%d)\n", gg[1], gg[2]));
  my(P=charpoly(M[3]), Y=matid(n));
  foreach(kl, k, my(Mk=M[3]^k-matid(n)); my(vA=valuation(polresultant(P, 'x^k-1), p));
     if(vA==0, next);
     my(a=ppart(Mk,p,full), b=ppart(BT^-1*Mk*BT,p,full), c=ppart(BTs^-1*Mk*BTs,p,full));
     printf("   k=%-4d v_p#A(F_2^k)=%d  A: %s   A'(T): %s   A'(T'): %s  %s\n", k, vA, pp_str(a,p), pp_str(b,p), pp_str(c,p), if(a==b&&a==c,"same","DIFFERS")));
  printf("   (%.1fs)\n", (getabstime()-t)/1000.)};
run(11, 23, 18, "conj", [1..60], 1);
run(11, 23, 1, "loword", [1..25], 1);
run(19, 457, 377, "conj", [8,16,152,304], 1);
run(19, 457, 250, "loword", [8,16,152], 1);
run(131, 263, 110, "conj", [1,2,4,131,262], 0);
run(131, 263, 262, "loword", [1,2,4,262], 0);
run(131, 146505763881528721, 14975904886999252, "conj", [1,2,3,4,5,6,7,8], 0);
