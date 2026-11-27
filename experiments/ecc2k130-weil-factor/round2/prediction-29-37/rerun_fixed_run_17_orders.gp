\\ usage: gp -q < run_17_orders.gp
\\ ell = 17 (small analogue of 131). Test the predicted rule on several orders R between Z[pi,pibar] and O_K:
\\   ppav with End = R exists  <=>  R is conj-stable and non-maximal at the self-glued inert place v0.
\\ Orders are R = Z[pi,pibar] + f for conj-stable conductors f supported on the primes above 271.
\\ T = Z[pi,pibar] t + f S with t a local generator of S at the primes of f.
\r decide_ppav_fixed.gp
initell(17);
p = 271; L = idealprimedec(NF, p);
mp = vector(#L, i, minpoly(resid(L[i], PI)));
v0 = select(i -> poldegree(mp[i]) == 1, [1..#L]);
pairs = List(); for (i = 1, #L, for (j = i+1, #L, if (poldegree(mp[i]) == 2 && mp[i] == mp[j], listput(pairs, [i, j]))));
printf("primes above 271: %d; self-glued v0 = %s; glued place-pairs = %d; distinct residue classes = %d\n", #L, v0, #pairs, #Set(mp));
Pp = cmtype(1); J0 = idealinv(NF, NF.diff);
monos(B0) = { my(out = List()); for (a = 0, 2*poldegree(NF.pol), for (b = 0, 1, listput(out, nfeltmul(NF, nfeltpow(NF, PI, a), nfeltpow(NF, PIB, b))))); Vec(out); }
MON = monos();
\\ at a glued place-pair the local generator t must be twisted by a unit w (REPORT.md sect. 4.3):
\\ we search t at random among local generators until the form is integral (expected ~ p-1 tries).
lat(S, f, lam) = {
  my(Fp = idealfactor(NF, f)[,1], H = idealhnf(NF, S), t = H[,1], tries = 0, B);
  while (1, while (!locgen(t, S, Fp), t = H * vectorv(#H, i, random(201) - 100));
    B = latt(S, f, t); tries++;
    if (denominator(gram(B, lam)) == 1, printf("    (local generator found after %d tries)\n", tries); return(B));
    t = 0; if (tries > 20000, error("no twist found")));
}
latt(S, f, t) = {
  my(gens = matconcat(vector(#MON, k, nfeltmul(NF, MON[k], t))));
  mathnf(concat(gens, idealhnf(NF, idealmul(NF, f, S))));
}
test(f, label, expect) = {
  my(r = decide(idealmul(NF, J0, idealinv(NF, f)), Pp));
  printf("%s: classes=%d principal=%d antisym=%d ppav-classes=%d  (predicted ppav: %d)\n", label, r[1], r[2], r[3], r[4], expect);
  if (#r[5], my(B = lat(r[5][1][1], f, r[5][1][2])); verify(B, r[5][1][2], Pp, Str(label, " explicit T")));
  r[4] > 0;
}
Q0 = idealhnf(NF, L[v0[1]]);
P1 = idealmul(NF, L[pairs[1][1]], L[pairs[1][2]]);
P2 = idealmul(NF, L[pairs[2][1]], L[pairs[2][2]]);
res = [test(Q0, "R: v0 only", 1), test(P1, "R: one glued place-pair only", 0), test(idealmul(NF, Q0, P1), "R: v0 + one pair", 1), test(idealmul(NF, P1, P2), "R: two pairs, v0 maximal", 0), test(idealhnf(NF, 1), "R = O_K", 0)];
printf("agreement with prediction: %d\n", res == [1, 0, 1, 0, 0]);
