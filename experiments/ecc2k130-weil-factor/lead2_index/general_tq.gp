\\ Stress test of the GENERAL statement proved in PROOF.md (Theorem 2):
\\ for alpha,beta roots of x^2 - t x + q (t,q in Z, t^2-4q not a square), O = Z[alpha], ell odd prime,
\\ U_ell = (alpha^ell-beta^ell)/(alpha-beta) != 0:
\\   [O (x) Z[zeta_ell] : Z[alpha z, beta z^-1]] = |U_ell|^(ell-2),  quotient = (Z/|U_ell|)^(ell-2).
\\ The point: this family contains many U_ell with SQUARE factors (prime-power case), which the
\\ x^2+x+2 sequence does not provide for ell <= 113.
default(parisizemax, 4000000000);
default(threadsizemax, 2000000000);
coords(e, ell) = {my(d=ell-1, v=vector(2*d), p=lift(lift(e)));
  for(k=0,1, my(c=polcoef(p,k,'w)); for(i=0,d-1, v[k*d+i+1]=polcoef(c,i,'z))); v~};
Ut(n,t,q) = {my(a=0,b=1); if(n==0,return(0)); for(k=2,n, [a,b]=[b,t*b-q*a]); b};
test(t,q,ell) = {my(d=ell-1, Pz=polcyclo(ell,'z), Pw='w^2-t*'w+q, z=Mod(Mod('z,Pz),Pw), al=Mod(Mod(1,Pz)*'w,Pw), be=t-al);
  my(pi=al*z, pb=be*z^(ell-1), cols=List(), x=z^0);
  for(i=0,2*d-1, listput(cols,coords(x,ell)); x*=pi);
  x=pb; for(j=1,2*d-1, listput(cols,coords(x,ell)); x*=pb);
  my(H=mathnf(matconcat(Vec(cols))), f=abs(Ut(ell,t,q)));
  if(matsize(H)[2]!=2*d, return([0,f,0]));
  my(snf=vecsort(matsnf(H),,4), exp=concat(vector(ell-2,k,f), vector(ell,k,1)));
  [abs(matdet(H))==f^(ell-2), f, abs(matdet(H)), snf==exp, snf]};
\\ Theorem 2(a) (index) is claimed for all (t,q); Theorem 2(b) (quotient (Z/f)^(ell-2)) only for gcd(t,q)=1.
ntest=0; nidxfail=0; nsq=0; nsqcop=0; nsnfcopfail=0; nsnfnoncop=0; nnoncop=0;
{forprime(ell=3,11, for(t=-5,5, for(q=-12,12,
   if(q==0 || issquare(t^2-4*q), next);
   my(f=abs(Ut(ell,t,q))); if(f==0 || f > 10^6, next);
   my(r=test(t,q,ell), cop=(gcd(t,q)==1)); ntest++;
   my(sq=(f>1 && !issquarefree(f))); nsq+=sq; nsqcop+=(sq&&cop); nnoncop+=!cop;
   if(!r[1], nidxfail++; printf("INDEX FAIL t=%d q=%d ell=%d U=%d index=%s\n",t,q,ell,f,r[3]));
   if(cop && !r[4], nsnfcopfail++; printf("SNF FAIL (coprime) t=%d q=%d ell=%d U=%d snf=%s\n",t,q,ell,f,r[5]));
   if(!cop && !r[4], nsnfnoncop++; printf("snf differs (gcd(t,q)=%d, outside Thm 2(b)) t=%d q=%d ell=%d |U|=%s snf=%s\n",gcd(t,q),t,q,ell,factor(f),r[5]));
   if(sq, printf("square-factor case: index_ok=%d snf_ok=%d gcd(t,q)=%d t=%d q=%d ell=%d |U|=%s index=%s\n", r[1], r[4], gcd(t,q), t,q,ell,factor(f),r[3])))));}
printf("SUMMARY: tested %d (t,q,ell) cases (ell<=11, |U|<=10^6); %d with non-squarefree |U_ell| (%d of them with gcd(t,q)=1); index failures: %d; SNF failures with gcd(t,q)=1: %d; non-coprime cases: %d, of which SNF differs: %d\n", ntest, nsq, nsqcop, nidxfail, nsnfcopfail, nnoncop, nsnfnoncop);
quit
