default(parisizemax, 4000000000);
default(realprecision, 400);
read("hminus_defs.gp");
{ my(v = 2*131*logprodK(131), h); print(precision(v), " ", default(realprecision)); h = round(2*real(v)); printf("h^-(K_131) = Q*w*prod with Q=2: %d\nfrac = %.2e, log2 = %.4f\nfactor (partial, <1e7): %s\n", h, abs(2*real(v)-h), log(h)/log(2), factor(h, 10^7)); }
