\\ Interpolating an isogeny on A[N] would need the points of A[N]; the Weil-dual line to A(F_2) has Frobenius
\\ eigenvalue 2, so A[N] is not rational over F_{2^k} unless ord_N(2) | k (embedding degree of N).
default(parisizemax, 4000000000);
s = 2; t = -1; a = [2, -1]; for (k = 2, 131, a = concat(a, -a[k] - 2*a[k-1]));
tr131 = a[132]; nE = 2^131 + 1 - tr131; N = nE / 4;
printf("#E(F_2^131) = %s ; N = %s ; N prime: %d ; bits %d\n", nE, N, isprime(N), #binary(N));
o = znorder(Mod(2, N));
printf("ord_N(2) = %s ; (N-1)/ord = %s\n", o, (N-1)/o);
