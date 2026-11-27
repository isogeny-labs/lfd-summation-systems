default(parisizemax, 4000000000);
realweil(P, g, q) = {my(R = P, h = 0); forstep(k = g, 0, -1, my(c = polcoeff(R, g+k, 'x)); h += c*'y^k; R -= c*'x^(g-k)*('x^2+q)^k); h;}
{
my(ell=131, P = Pol(Vecrev(readvec("PA131.txt")),'x), h = realweil(P,130,2), dh = poldisc(h), NT = abs(polresultant(h,'y^2-8)));
my(iK2 = dh^2*NT/(ell^258*7^130), u = 263*146505763881528721, ip2 = abs(dh)/(ell^129*7^65));
print("U*: ", iK2 == u^258, "  plus index == (263*p2)^64 : ", ip2 == u^128, "  ratio iK2/u^258 = ", iK2/u^258);
print("N(t^2-8) = 7^a * ... : ", factor(NT, 10^6)[1..min(5,#factor(NT,10^6)~),]);
}
