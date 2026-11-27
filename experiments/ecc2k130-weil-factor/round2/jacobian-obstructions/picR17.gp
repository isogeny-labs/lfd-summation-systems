\\ ell=17: |Pic R| for R = O_K + Q (glued at a conj-stable degree-2 prime Q above 271) via the ray class group
\\ mod Q (GRH, bnfinit), vs the heuristic h^-*(p+1)/(2*ell) used for the minus local factor.
default(parisizemax, 4000000000);
{
my(ell = 17, p = 271, pK = polredbest(polcompositum(polcyclo(ell), 'x^2+7)[1]), bK = bnfinit(pK, 1), Q, bnr, m);
foreach(idealprimedec(bK, p), P, if(P.f == 2, Q = P; break));
bnr = bnrinit(bK, Q);
m = (p^2-1)*bK.no/bnr.no/(p-1);
printf("ell=17: h_K=%d |Cl_Q|=%d => image of O_K^* in (O/Q)^* has order (p-1)*%d\n", bK.no, bnr.no, m);
printf("   |Pic R| = h*(p+1)/m = %d ; heuristic h^-*(p+1)/(2*ell) = %d\n", bK.no*(p+1)/m, 45*(p+1)/(2*ell));
}
