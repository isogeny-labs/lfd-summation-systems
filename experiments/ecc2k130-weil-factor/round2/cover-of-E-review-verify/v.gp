default(parisizemax,4000000000);
\\ E: t=-1, alpha root of x^2+x+2
ell=131; PE=x^2+x+2;
U(n)=my(a=Mod(y,y^2+y+2)); lift((a^n-(-1-a)^n)/(2*a+1)); \\ (a^n-b^n)/(a-b), a-b=2a+1
u=U(ell); print("U131=",u," factor=",factor(u));
\\ P_W(x)=x^(2l)-t_l x^l+2^l, t_l=a^l+b^l
a=Mod(y,y^2+y+2); tl=lift(a^ell+(-1-a)^ell); tl=polcoeff(tl,0);
PW=x^(2*ell)-tl*x^ell+2^ell; [PA,r]=divrem(PW,PE); print("rem=",r);
print("trace A=",-polcoeff(PA,poldegree(PA)-1));
v=subst(lift(PA),x,a); v=lift(v); c=content(v); print("content P_A(alpha)=",factor(c)," eq=",c==ell*abs(u));
\\ alternative: P_A(alpha)= ell*alpha^(ell-1)*U
print("identity check: ", lift(subst(PA,x,a)-ell*a^(ell-1)*u)==0);
print("ord_131(2)=",znorder(Mod(2,131)));
print("#E(F4)=",4+1-lift(a^2+(-1-a)^2));
print("2^d-1 d<=4:",vector(4,d,(2^d-1)%131)," 4^d-1:",vector(4,d,(4^d-1)%131));
\\ divisors of 131*263*p2 <=390
print(select(d->d<=390,divisors(131*abs(u))));
\\ negative control ell=5,7 content
for(l=3,13,if(isprime(l),tl=polcoeff(lift(a^l+(-1-a)^l),0);P=(x^(2*l)-tl*x^l+2^l)/PE;w=lift(lift(subst(P,x,a)));print(l," content=",content(w)," l*|U|=",l*abs(U(l)))));
