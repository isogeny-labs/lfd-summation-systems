\\ Independent check: which modules T=I (End=O_K) of K=Q(zeta_l,sqrt-7) carry a principal polarization
\\ for the CM type induced from Q(sqrt-7).  Written from scratch (does not reuse decide_ppav.gp or rt_lib.gp).
default(parisizemax,4000000000);
isvec0(v)=type(v)!="t_INT";  \\ matsolvemod returns 0 (t_INT) iff no solution
run(l)={
  my(C=polcompositum(polcyclo(l,x), x^2+x+2, 1)[1], RB=polredbest(C[1],1), R=RB[1], mp=RB[2], Z=subst(lift(C[2]),x,mp), AL=subst(lift(C[3]),x,mp));
  my(bnf=bnfinit(R,1), nf=bnf.nf, n=poldegree(R));
  my(z=nfalgtobasis(nf,Z), al=nfalgtobasis(nf,AL), s7=nfalgtobasis(nf,2*AL+1));
  \\ complex conjugation: z->1/z, al->-1-al
  my(G=nfgaloisconj(nf), sig=0);
  for(i=1,#G, my(g=G[i]); if(nfgaloisapply(nf,g,z)==nfeltdiv(nf,1,z) && nfgaloisapply(nf,g,al)==nfalgtobasis(nf,-1-AL), sig=g));
  if(sig==0,error("no conj"));
  my(cj(a)=nfgaloisapply(nf,sig,a));
  \\ check conj is complex conj numerically at first embedding
  my(tu=bnf.tu, w=tu[1], t=nfalgtobasis(nf,tu[2]), fu=apply(u->nfalgtobasis(nf,u),bnf.fu), r=#fu);
  if(w!=2*l,error("torsion"));
  my(tpow=vector(w,k,nfeltpow(nf,t,k-1)));
  my(dlog(m)=for(k=1,w,if(tpow[k]==m,return(k-1))); error("not a root of unity"));
  \\ map U -> mu_w : u -> conj(u)/u
  my(A=matrix(1,r+1)); A[1,1]=dlog(nfeltdiv(nf,cj(t),t)); for(i=1,r,A[1,i+1]=dlog(nfeltdiv(nf,cj(fu[i]),fu[i])));
  my(uexp(c)=my(u=nfeltpow(nf,t,c[1]));for(i=1,r,u=nfeltmul(nf,u,nfeltpow(nf,fu[i],c[i+1])));u);
  \\ real units: kernel of A mod w (lattice), basis via matkermod-like: use HNF of [A; w]
  my(K=matkerint(concat(A,Mat(w)))); \\ integer kernel of [A | w]
  my(rb=vector(#K,j,uexp(vector(r+1,i,K[i,j]))));
  rb=concat(rb,[nfalgtobasis(nf,-1)]);
  foreach(rb,u,if(cj(u)!=u,error("real unit failed")));
  my(sg(a)=apply(v->(real(v)<0), nfeltembed(nf,a)));  \\ 1 where negative (r is real)
  my(M=matrix(n/2,#rb,i,j,sg(rb[j])[i]));
  print("l=",l," degree=",n," h=",bnf.no," cyc=",bnf.cyc," rank(real-unit signs mod 2)=",matrank(Mod(M,2)), " of ",n/2);
  my(Dinv=idealinv(nf,nf.diff), cyc=bnf.cyc, gens=bnf.gen, nc=0,np=0,nanti=0,nok=0);
  forvec(e=vector(#cyc,i,[0,cyc[i]-1]),
    nc++;
    my(I=idealhnf(nf,1)); for(i=1,#cyc,I=idealmul(nf,I,idealpow(nf,gens[i],e[i])));
    my(Ib=idealhnf(nf,matconcat(vector(n,k,cj(I[,k])))));
    my(J=idealmul(nf,Dinv,idealinv(nf,idealmul(nf,I,Ib))));
    my(pr=bnfisprincipal(bnf,J)); if(pr[1]!=0*pr[1],next); np++;
    my(g=pr[2]);
    \\ need u with cj(u)/u = -g/cj(g)
    my(tgt=dlog(nfeltdiv(nf,nfeltmul(nf,-1,g),cj(g))), sol=matsolvemod(A,w,[tgt]~));
    if(!isvec0(sol),next); nanti++;
    my(lam=nfeltmul(nf,g,uexp(sol))); if(cj(lam)!=nfeltmul(nf,-1,lam),error("anti"));
    if(idealhnf(nf,lam)!=J,error("gen"));
    my(rr=nfeltdiv(nf,lam,s7)); if(cj(rr)!=rr,error("r not real"));
    my(x);
    x=matsolvemod(lift(Mod(M,2)),2,sg(rr)~);
    if(!isvec0(x),next);
    my(lam2=lam); for(j=1,#rb,if(x[j]%2,lam2=nfeltmul(nf,lam2,rb[j])));
    my(rr2=nfeltdiv(nf,lam2,s7)); if(vecsum(sg(rr2))!=0,error("pos fail"));
    \\ verify Gram of Tr(lam2 x ybar) on I is unimodular & alternating
    my(Gm=matrix(n,n,i,j,nfelttrace(nf,nfeltmul(nf,lam2,nfeltmul(nf,I[,i],cj(I[,j]))))));
    if(abs(matdet(Gm))!=1 || Gm+Gm~!=0, error("gram"));
    nok++);
  print("  classes=",nc," principal J=",np," antisym gen=",nanti," ppav=",nok,"   eps=",kronecker(-7,l));
};
for(i=1,#L,run(L[i]));
