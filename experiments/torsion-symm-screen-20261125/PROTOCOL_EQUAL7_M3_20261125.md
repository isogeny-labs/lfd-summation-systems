# Exploratory equal-7 three-summand Koblitz panel — frozen 2026-11-25

This follow-up reduces factor-base dimension from the equal-11 panel
because that direct arm reached full base-column rank and admitted all 24
targets. It is an adaptive new candidate, not a replication of that panel.
No target relation or degree result for this base was inspected before this
protocol was committed.

Use the same `K_1/F_(2^9)` curve, modulus 515, group order `518=14·37`,
generator `Q=(135,123)` and first 24 targets `Q,...,24Q` as the equal-11
protocol. The base-size-only rule is the first linearly independent basis
`[1,a,b]` with `2<=a<b<512`, equal but distinct rational `x∈V` and
`u=1/(x+1)∈V` bases of at least four abscissae, and at least three
nonzero cofactor-projected columns per base. It selects
`V=span_F2[1,2,4]`. Each frame has seven rational abscissae, 13 signed
points and six nonzero projected columns; the abscissa sets intersect in
one value. The `u` base is stable under order-two translation because
`1∈V`. No Frobenius orbit compression is used.

Use the exact same direct `S4` and native 18-term symmetrized `S4`
constructions as `PROTOCOL_EQUAL11_M3_20261125.md`, with nine and seven
unknown bits respectively. On every target count all ordered signed
triple attempts, group-checked decompositions, distinct relation rows,
admitted targets and factor-base-column rank modulo 37. Known projected
target multiples are right-hand sides. Verify saved rows independently
against subgroup logs; retain failures and torsion/identity boundaries.

For targets 1–4, check complete Boolean solution sets against finite-field
evaluation and group-law lifts. Compute exact iterated LFD with degree cap
8, at most 512 normal-form columns, 512 MiB per matrix and 60 seconds per
arm/target. Record input degree separately, preserve censored cells, and
compare LFD only on matched target classes. Charge base construction,
equation construction, every attempted triple, lifting/checking, rank
updates and closure work in wall and CPU. A positive joint lead would
still require independent fresh-instance replication, a complete
factor-base solve/descent and matched cold same-target rho before any
end-to-end cost claim.
