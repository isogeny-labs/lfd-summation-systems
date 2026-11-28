\\ Lead 2: exact HNF index [O_F (x) Z[zeta_ell] : Z[pi,pibar]] for all primes 3 <= ell <= LMAX.
\\ O_F (x) Z[zeta] = Z[z]/Phi_ell (x) Z[w]/(w^2-w+2), w=(1+sqrt-7)/2, alpha=w-1, beta=-w (roots of x^2+x+2).
\\ pi = alpha*z, pibar = beta*z^{-1}. Z[pi,pibar] = Z-span{pi^i (0<=i<2d), pibar^j (1<=j<2d)}, d=ell-1
\\ (pi*pibar = 2 and pi, pibar satisfy monic integer polynomials of degree 2d).
\\ For ell != 7 the ambient ring is O_K, K=Q(zeta_ell,sqrt-7). For ell=7 it is the (non-domain) tensor algebra;
\\ see field7() below for the honest field case.
default(parisizemax, 4000000000);
default(threadsizemax, 2000000000);
LMAX = 113;
coords(e, ell) = {my(d=ell-1, v=vector(2*d), p=lift(lift(e)));
  for(k=0,1, my(c=polcoef(p,k,'w)); for(i=0,d-1, v[k*d+i+1]=polcoef(c,i,'z))); v~};
U(n) = {my(a=0,b=1); if(n==0,return(0)); for(k=2,n, [a,b]=[b,-b-2*a]); b};
run(ell) = {my(d=ell-1, Pz=polcyclo(ell,'z), Pw='w^2-'w+2, z=Mod(Mod('z,Pz),Pw), w=Mod(Mod(1,Pz)*'w,Pw));
  my(pi=(w-1)*z, pb=(-w)*z^(ell-1), cols=List());
  my(x=z^0); for(i=0,2*d-1, listput(cols,coords(x,ell)); x*=pi);
  x=pb; for(j=1,2*d-1, listput(cols,coords(x,ell)); x*=pb);
  my(M=matconcat(Vec(cols)), H=mathnf(M));
  if(matsize(H)[2]!=2*d, error("rank deficient at ell=",ell));
  my(idx=abs(matdet(H)), u=U(ell), f=abs(u));
  my(snf=matsnf(H), expsnf=concat(vector(ell-2,k,f), vector(ell,k,1)));
  my(snfsorted=vecsort(snf,,4), ok1=(idx==f^(ell-2)), ok2=(snfsorted==expsnf));
  my(inZ=(denominator(matsolve(H, coords(z,ell)))==1));
  my(fOK = (denominator(matsolve(H, f*matid(2*d)))==1)); \\ f*O_K inside Z[pi,pibar]?
  my(o2=znorder(Mod(2,ell)));
  printf("%d|%d|%s|%s|%s|%d|%d|%d|%d|%d|%d\n", ell, o2, u, if(f==1,"1",Str(factor(f))), idx,
     ok1, ok2, inZ, (f==1), fOK, issquarefree(f));
  [ell,o2,u,idx,ok1,ok2,inZ,fOK]};
print("# ell|ord_ell(2)|U_ell|factor|U_ell||index|index==|U|^(ell-2)|SNF==(Z/|U|)^(ell-2)|zeta_in_R|U=+-1|f*O_K in R|U squarefree");
forprime(ell=3, LMAX, run(ell));
\\ ell = 7 as an honest field: K = Q(zeta_7) contains sqrt(-7) = Gauss sum g. Both identifications of alpha.
field7() = {my(Pz=polcyclo(7,'z), z=Mod('z,Pz), g=z+z^2+z^4-z^3-z^5-z^6);
  if(g^2!=-7, error("gauss sum"));
  for(s=-1,1, if(s==0,next); my(sq=s*g, al=(-1+sq)/2, be=(-1-sq)/2);
    if(al^2+al+2!=0, error("alpha"));
    my(pi=al*z, pb=be*z^6, cols=List(), x=z^0);
    for(i=0,5, listput(cols, Colrev(lift(x),6)); x*=pi);
    x=pb; for(j=1,5, listput(cols, Colrev(lift(x),6)); x*=pb);
    my(H=mathnf(matconcat(Vec(cols))));
    printf("ell=7 field Q(zeta_7), sqrt-7 = %d*g: rank %d, index [Z[zeta_7]:Z[pi,pibar]] = %s, SNF %s, zeta in R: %d\n",
      s, matsize(H)[2], abs(matdet(H)), matsnf(H), denominator(matsolve(H,Colrev(lift(z),6)))==1)) };
field7();
quit
