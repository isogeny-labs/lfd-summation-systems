# Frozen binary Boolean torsion follow-up — 2026-11-25

## Fixed case and question

Use the existing `ic-curve-sweep` Koblitz case
`7c65f50a88d2c5e980b2`: `n=9`, modulus `0x203`,
`E: y^2+xy=x^3+x^2+1`, Frobenius-invariant basis `[302,252,1]`,
and its planted target `(467,0)`. Also retain the other two existing
`a=1`, `n=9`, `k=3` targets for the same base, in their case-file order:
`1bac99e1e73e08b00999` and `a5eb99758876c7f2fee6`.
Do not select targets after looking at the degree outcomes.

The rational two-torsion point is `(0,1)`. For nonzero `x`, translation by
it sends `x` to `1/x`. The frozen base includes `x=0`; its seven nonzero
elements are inverse-closed. Keep the full frozen eight-abscissa base for
the direct comparison, and explicitly record that the torsion action takes
`x=0` outside the affine base.

## Equivalent equations on the same Boolean variables

For two summands put `u=x1+x2`, `v=x1*x2`, and `z=x(R)` in `F_(2^9)`.
The direct equation is the repository's `S3=(v+zu)^2+zv+1=0`,
Weil-descended to six Boolean variables. For `v != 0` define the
simultaneous-inversion invariants

`s=v+1/v`, `w=u^2/v=x1/x2+x2/x1`.

Then `S3/v=s+z^2*w+z`. Use `inverse(0)=0` only as a total Boolean function
for generating coordinate polynomials. Both arms use the identical six
Boolean variables, field, curve, factor base and target.
Check the identity for every pair of the eight base elements and verify
the two equation systems have exactly the same 64-point Boolean truth table
before computing degree.

## Pre-degree boundary amendment

The first finite truth-table preflight found two direct solutions at `v=0`
for the planted target: `(x1,x2)=(0,302),(302,0)`. The proposed nonzero guard
would delete these valid frozen-base decompositions, so it cannot be used for
the matched comparison. This was discovered before computing any LFD.

Set `delta=1+v^(2^9-1)`, which is 1 at `v=0` and 0 otherwise. The corrected
invariant presentation uses the nine coordinate equations of

`H=delta*S3+(1+delta)*(s+z^2*w+z)=0`.

It uses the invariant quotient equation for `v!=0` and the original Semaev
equation on the boundary `v=0`. This is a patched invariant presentation,
not a globally torsion-invariant coordinate system on the frozen affine
base. The truth-table identity and exact solution-set check remain required.

## Measurements and decision gates

For each target, preserve the exact direct and invariant Boolean generator
lists and hashes, solution masks, input degree, and exact iterated LFD from
the repository closure engine. Use a degree cap of 8 and record a censored
result if the bound is not reached. Independently count rational checked
decompositions on the same target/base, including failures; record the
resulting relation rows and rank over any odd prime subgroup order when
applicable. The same base and target imply identical checked relation sets
for the two presentations; the degree comparison alone cannot pass the
two-part yield-plus-degree gate.

Charge polynomial construction and closure wall/CPU, plus the complete
enumeration, relation checking and rank work. These measurements are
component costs; no complete same-target rho comparison is authorized by
this protocol. Preserve all negative and censored outcomes. Keep this
follow-up separate from the earlier binary Edwards scalar-coordinate
screen and from the prime-field PR #60.
