p=146505763881528721;
print("p2 < 2^64: ", p < 2^64);
c=primecert(p); print("primecert(p2) = ", c, "  valid: ", primecertisvalid(c));
c1=primecert(p,1); print("primecert(p2,1) [Pratt flag] = ", c1, "  valid: ", primecertisvalid(c1));
c2=primecert(263); print("primecert(263) = ", c2, " valid: ", primecertisvalid(c2));
print("negative: primecertisvalid(p2+2) = ", primecertisvalid(p+2), " ; primecertisvalid(561) = ", primecertisvalid(561));
N=nextprime(2^80); c3=primecert(N); print("control 2^80-range prime: ECPP cert length ", #c3, " valid: ", primecertisvalid(c3));
