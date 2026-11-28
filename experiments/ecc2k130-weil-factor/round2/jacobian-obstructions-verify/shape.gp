default(parisizemax,4000000000);
mk(l)=my(R=polresultant(polcyclo(l,y), x^2+x*y+2*y^2, y)); R/pollead(R);
realpart(l)={my(P=mk(l),g=l-1,h=0,Q=P); forstep(k=g,0,-1, my(c=polcoeff(Q,g+k)); h+=c*t^k; Q-=c*x^(g-k)*(x^2+2)^k); if(Q!=0,error("bad")); h};
h=realpart(131);
foreach([263,146505763881528721],p, my(F=factormod(h,p)); my(sh=Map()); for(i=1,#F~, my(k=[poldegree(F[i,1]),F[i,2]]); mapput(sh,k,if(mapisdefined(sh,k,&v),v+1,1))); print(p," (deg,mult)->count: ",Mat(sh)));
P=mk(131); foreach([263,146505763881528721],p, my(F=factormod(P,p)); my(sh=Map()); for(i=1,#F~, my(k=[poldegree(F[i,1]),F[i,2]]); mapput(sh,k,if(mapisdefined(sh,k,&v),v+1,1))); print("P_A mod ",p," (deg,mult)->count: ",Mat(sh)));
quit
