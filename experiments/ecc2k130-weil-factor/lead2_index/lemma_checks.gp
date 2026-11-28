\\ Numerical sanity checks of the individual lemmas in PROOF.md for alpha^2+alpha+2=0, all primes 3<=ell<=400.
\\ (The lemmas are proved in PROOF.md; this guards against algebra slips.)
default(parisizemax, 4000000000);
F = nfinit(x^2+x+2); \\ O_F = Z[alpha], alpha = x
U(n) = {my(a=0,b=1); if(n==0,return(0)); for(k=2,n, [a,b]=[b,-b-2*a]); b};
\\ coordinates of alpha^n in basis (1,alpha): alpha^n = -2 U_{n-1} + U_n alpha
co(n) = {my(v=lift(Mod(x,x^2+x+2)^n)); [polcoef(v,0), polcoef(v,1)]};
bad=0;
{forprime(ell=3,400, my(f=abs(U(ell)));
  \\ Lemma 3: alpha^ell = a + U_ell*alpha (so A = Z + fO)
  if(co(ell)[2]!=U(ell), bad++; print("L3 fail ",ell));
  \\ Lemma 4: [O : Z alpha^i + Z beta^(ell-i)] = f, beta = -1-alpha
  for(i=1,ell-1, my(a=co(i), b=lift(Mod(-1-x,x^2+x+2)^(ell-i)), bb=[polcoef(b,0),polcoef(b,1)]);
     if(abs(a[1]*bb[2]-a[2]*bb[1])!=f, bad++; print("L4 fail ",ell," ",i)));
  \\ Remarks: f odd; gcd(f,7*ell)=1 for ell!=7; every prime p|f has p = (-7/p) = +-1 mod ell
  if(f%2==0, bad++; print("odd fail ",ell));
  if(ell!=7 && gcd(f,7*ell)!=1, bad++; print("gcd fail ",ell));
  if(ell<=200, my(P=factor(f)[,1]); for(k=1,#P, if((P[k]-kronecker(-7,P[k]))%ell!=0, bad++; print("p=+-1 fail ",ell," ",P[k]))));
  \\ U_ell = (-7/ell) mod ell
  if(ell!=7 && (U(ell)-kronecker(-7,ell))%ell!=0, bad++; print("U_ell mod ell fail ",ell)));}
print("lemma checks for primes 3..400: failures = ", bad);
\\ squarefreeness of U_ell (prime-power case) for ell up to 200 (factoring may be partial beyond; use issquarefree)
{my(L=List()); forprime(ell=3,200, if(!issquarefree(U(ell)), listput(L,ell)));
 print("primes ell<=200 with U_ell NOT squarefree: ", L);}
quit
