default(parisizemax,4000000000); default(realprecision,400);
\\ Weil poly as resultant, my own construction
mk(l)=my(R=polresultant(polcyclo(l,y), x^2+x*y+2*y^2, y)); R/pollead(R);
P=mk(131);
print("deg ",poldegree(P), " P(1)=",subst(P,x,1));
print("P(1)==ECC2K130 N: ", subst(P,x,1)==680564733841876926932320129493409985129);
\\ point counts: s_k power sums
s(Q,k)=my(r=polroots(Q)); round(real(sum(i=1,#r,r[i]^k)));
\\ exact Newton: coefficients
c=Vec(P); print("coef x^259 = ",c[2]," coef x^258 = ",c[3]);
s1=-c[2]; s2=s1^2-2*c[3]; print("s1=",s1," s2=",s2," N1=",3-s1," N2=",5-s2);
Pm=subst(P,x,-x); cm=Vec(Pm); t1=-cm[2]; t2=t1^2-2*cm[3]; print("twist N1=",3-t1," N2=",5-t2);
\\ x^k structure
print("nonzero coeff positions gcd: ", gcd(vector(#c,i,if(c[i],261-i,0))));
\\ mod 2
print("P mod 2 == x^130*Phi131 mod2: ", Mod(1,2)*P==Mod(1,2)*x^130*polcyclo(131));
print("znorder(2 mod 131)=",znorder(Mod(2,131)), " polisirreducible Phi131 mod2: ", polisirreducible(Mod(1,2)*polcyclo(131)));
\\ negative control ell=7: 2 has order 3 mod 7
P7=mk(7); print("ell=7 P mod 2 = ", factor(Mod(1,2)*P7));
quit
