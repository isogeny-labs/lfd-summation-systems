# Equal-11 Koblitz three-summand torsion-frame panel — frozen 2026-11-25

## Preflight and target stream

Use `K_1: y²+xy=x³+x²+1` over `F_(2^9)` with polynomial modulus `515`.
The group order is `518=14·37`. The first rational point in ascending
abscissa whose cofactor-14 projection is nonidentity is `(2,47)`; fix
`Q=[14](2,47)=(135,123)` of order 37. Targets are the first 24 nonzero
multiples `Q,...,24Q`, in order, with every failed target retained.

The base was selected using only the curve, rational base sizes and
nontrivial projected-column counts, before any target relation or degree
measurement: enumerate linearly independent integer bases `[1,a,b,c]` with
`2<=a<b<c<512` in lexicographic order, require both `x∈V` and
`u=1/(x+1)∈V` to have equal rational abscissa counts of at least six,
different abscissa sets, and at least four nonzero columns after projection
by 14. The first is `V=span_F2[1,2,4,56]`. Both bases have 11 rational
abscissae, 21 signed points and 10 nonzero projected columns; their
abscissa sets intersect in one value. The `u` base is stable under
translation by `(0,1)` because `1∈V`. Frobenius invariance is not required
for the three-summand polynomial or this factor base; no orbit compression
is used.

## Frozen comparison and gates

On each target, compare the repository's complete direct Boolean `S4`
descent with three `x_i∈V` (12 unknown bits) to the native 18-term
symmetrized `S4` over `(w1,w2,w3,wR,s)` with
`w_i=u_i²+u_i`, `u_i=eps_i+Σ c_ij b_j`, and one parity bit
`eps=eps1+eps2+eps3` (10 unknown bits). Freeze the exact native term list
from `src/cryptanalysis/koblitz_symmetrised.rs::symmetrised_terms(3)` before
measuring. The target is fixed, not an unknown. Record input degree and
complete Boolean solution sets, checked against independent field
evaluation and group-law lifting including all identity/torsion boundaries.

Count admitted targets, distinct checked relation rows and **factor-base
column rank** modulo 37 after each target. Place the known cofactor-projected
target multiple on the right-hand side; do not add one formal target column
per target when deciding useful equation rank. A yield gain requires a
strict rank or admitted-target increase at equal total base size. Check
rows against independent discrete logs in the order-37 subgroup.

Compute exact iterated Boolean LFD on the first four targets with degree
cap 8, at most 4,096 squarefree normal-form columns, 512 MiB per closure
and 60 seconds per arm/target. Preserve capped or failed cells. A degree
gain requires strictly lower exact LFD on a matched target class. Report
native solver degree separately if measured; input degree alone is not an
LFD or total-solving-degree result.

Charge both base constructions, equation construction, all 24 attempted
targets, failed decompositions, root/sign lifting, relation verification,
rank updates and closure work in wall and CPU. A complete factor-base
linear solve, descent, and matched cold rho would require a separate
end-to-end protocol; do not infer them from this panel. Keep this panel
separate from the equal-67 two-summand screen and the frozen n=9
same-base Boolean screen.
