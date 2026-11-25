# Frozen matched torsion crossover protocol (2026-11-25)

Question: for a fixed public curve, identical targets and identical factor-base points, does the two-torsion symmetrized S3 relation finder return independently useful relations faster in total than the best direct method? A negative result applies only to the tested implementation and range.

Fixture: E/F_37889: y^2=x^3+x, #E=38224=16*2389 (Sage checked before freezing). T=(0,0). For N in {8,16,32,64,128,256,512,1024}, H=mu_N and V=H minus {1}; the factor base is every rational point with x in V, both signs. Thus x -> 1/x preserves membership. The coordinate u=(x+1)/(x-1), q=(u+v)^2, r=uv gives P(q,r,t)=(1+t^2)q-4t^2r-2t(r^2-1)=0 for t=x(R).

Target panel per N: 256 random same-target subgroup points R=uG+vQ, G the first nonzero [16]P in ascending rational-point search, Q=[113]G, (u,v) drawn from a frozen Python PRNG seed 20261125+N. Add 32 planted sums of factor-base points using a separate seed 20262125+N. Preserve every target including duplicates and failed attempts. Planted targets check solution recovery; only random subgroup-target relations enter the independent-rank calculation over 2389.

Arms:
1. Direct per-target group scan: build a point lookup once, then for each target and each P in the factor base test whether R-P is also in the base.
2. Direct pair-sum table: enumerate unordered factor-base point pairs once, index their group sum, and query each target.
3. Torsion-quotient S3 table: enumerate unordered factor-base abscissa pairs, map them to (q,r) orbits, solve P(q,r,t)=0 for t, index candidate abscissa pairs by t, then recover and group-check signed point pairs for each target.

All three arms must return the same canonical point-pair set per target. The group check is mandatory for the polynomial arm; a polynomial root by itself is not a relation. All setup, every target, all failures, root recovery, and validation are timed with wall and process CPU time. Each arm runs three times, with rotated order. Report medians and raw per-arm/per-rung receipts. Memory and timeout/censored states are retained. The comparison uses the fastest measured direct arm at each N and the complete batch time divided by the same independent relation rank; if rank is zero, compare batch completion time and report no cost-per-relation ratio. The first crossover needs the torsion arm to beat both direct arms on wall and CPU in all three repeats at the same N, with identical relation sets and positive rank. If none does across all eight N, record no crossover in this frozen range. No extrapolation to other curves, m, or ECC2K-130.
