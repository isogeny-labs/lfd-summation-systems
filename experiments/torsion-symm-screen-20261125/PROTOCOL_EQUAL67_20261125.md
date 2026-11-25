# Equal-67 Koblitz torsion-frame panel — 2026-11-25

## Preflight and frozen workload

Use `K_1: y²+xy=x³+x²+1` over `F_(2^15)` with polynomial modulus `0x8003`.
The curve order is 32,494 = 154·211. Choose the first lexicographic
rational point whose multiplication by 154 is nonidentity; the resulting
order-211 point is `Q=(15115,1293)`. Use the first 24 nonzero multiples
`Q,2Q,...,24Q` for the same-target degree panel. A separate full-subgroup
relation inventory may extend through `210Q` in this fixed order, retaining
every failed target.

The first non-subfield, proper Frobenius-invariant dimension-seven subspace
containing 1 that has equal rational base sizes in the two native frames is
`V=span_F2[1,892,1694,3678,5682,8344,22848]`. This was selected by
base-size-only preflight, before relation or degree outcomes. Define

- `F_x`: rational points with `x ∈ V`;
- `F_u`: rational points with `u=1/(x+1) ∈ V`, omitting the pole and
  retaining the rational order-two point at `x=0`.

Both have 67 rational abscissae, 133 signed rational points, and 60
nonzero columns after projection by 154 to the order-211 subgroup. Their
abscissa sets intersect in seven values. Record the exact sets and hashes.
The `u` base is stable under translation by `(0,1)` away from the identity
boundary; the `x` base is the same-size direct-coordinate control.

## Equations and yield

For each target, the direct arm uses the repository's complete Boolean
Weil descent of `S3(x1,x2,xR)` with `x_i ∈ V` (14 variables). The
symmetrized arm uses the existing native Koblitz formula
`w1*w2*wR+w1+w2+wR+s=0`, where `w_i=u_i²+u_i`, `u_i=eps_i+Σ c_ij b_j`,
`s=eps+Σ c_ij b_j+uR`, and `eps=eps1+eps2` (13 variables). The target
is omitted as an unknown in both arms. If `uR` is undefined or zero,
retain a censored symmetrized cell. Check the solution set of each system
against independent finite base-pair enumeration and group-law lifting.

Count admitted targets, distinct verified relation rows and rank modulo
211 after each target. A yield increase requires a strict rank or admitted-
target increase at equal total base size, not duplicate sign or torsion
images. Do not use target outcomes to choose the base.

Compute exact iterated LFD for the first four targets in the fixed stream,
with degree cap 8, at most 8192 normal-form columns, 512 MiB per closure,
and 60 seconds per arm/target. Preserve censored cells. The degree gate
requires a strict reduction on the *same* target class with complete
solution-set checks. Report input degree separately from LFD.

Charge base construction, equation construction, all attempted targets,
failed decompositions, root/sign lifting, relation verification and rank
updates in wall and process CPU. If a complete factor-base linear solve
and matched cold same-target rho comparison are run, report them separately
and include setup; do not infer a full ECDLP cost from relation-collection
components. Keep this panel separate from the n=9 frozen-base test and the
archived native `n=17,23,31` pipeline comparison.
