# Frozen characteristic-two matched screen

## Construction and workload

Use the binary Edwards model `E_(d,d)` with `d=1` over the repository's
polynomial-basis fields `F_(2^n)` for `(n,l)=(5,2),(7,3),(9,4)`, with moduli
`0x25, 0x83, 0x203`. Here `t(P)=x(P)+y(P)`, the rational order-four point is
`T4=(1,0)`, and translation by `T4` sends `t` to `t+1`. Verify the complete
Edwards addition law and the published third summation polynomial by finite
field enumeration before measuring. Retain all group-law and formula failures.

Set `V_j=span_(F2){1,theta,...,theta^(j-1)}`. The baseline permits both
summands from `W=V_(l+1)`; the colored arm permits one from `A=V_l` and one
from `B=theta^l+V_l`. Thus `W=A disjoint-union B` and the **total** available
factor-base set is identical. Both constructions use the same canonical
representative of each order-four torsion orbit, with signs and explicit
torsion correction when checking a relation. Select the largest odd prime
divisor `r` of `#E`; choose the first lexicographic point whose multiplication
by `#E/r` is nonidentity as prime-order generator `Q`. Targets are
`Q,2Q,...,min(24,r-1)Q`, in that order, for every arm. Count every target,
including failed attempts.

The direct scalar equation is the paper's `f3(t1,t2,t(R))`. The invariant
equation uses `e=t1+t2`, `s=(t1^2+t1)(t2^2+t2)` and
`f3=e^2+s(t(R)^2+t(R))+t(R)^2` for `d=1`. Verify equality for all field
values used. For each base, derive the exact finite image ideal of allowed
`(e,s)` pairs by elimination from the two base polynomials; retain and charge
this setup. Compare four presentations: baseline direct/invariant and colored
direct/invariant. All comparisons use the same targets within a cell.

## Measurements and gates

Enumerate all candidate point decompositions independently of the polynomial
solver. Record the number of target points admitting at least one checked
relation, distinct checked relation rows, and matrix rank over `F_r` after
each target, for both arms. The colored arm may use only cross-partition
pairs. A claimed yield increase requires a strict increase in independent
useful relation rank or admitted targets at the same total base size; duplicate
symmetry images do not count.

For `n=5,7`, measure exact unit-ideal LFD for the first four targets in the
fixed stream without a checked decomposition, separately for all four
presentations. Include **all** factor-base and field equations in the closure.
Use a degree cap of 12, at most 2,000 monomial columns, and a 60-second cap
per closure; retain censored cells. Also record Gröbner solving degree and
wall/CPU for the same systems when available. The `n=9` cell is yield and
construction-cost only. Record setup, every attempted target, failed solves,
root lifting and relation checks, and rank work; do not promote a sparse
equation or a component timing to a total-cost claim.

The positive result sought is **both** higher independent useful-relation
yield at equal total base size and lower exact LFD or total solving degree on
the same targets, with a favorable total charged cost. Otherwise record the
specific failed gate. This screen is separate from the merged prime-field
PR #60 and does not extrapolate to ECC2K-130.

## Amendment after the initial yield enumeration

The `n=5` and `n=7` cells have no non-torsion column in the first colored
coset, so their colored ideals are the unit ideal before the summation
equation. Their degree values are retained as degenerate controls. Extend the
same exact-LFD protocol to the first four baseline-failed `n=9` targets,
where both colors have nonempty columns (two and four). This is a labelled
follow-up selected after inspecting the initial yield data, not part of the
original confirmatory panel. The degree/column/time caps are unchanged.

The `n=9` baseline-direct arm then reached the predeclared degree cap 12 on
all four targets, while the other three presentations completed. Preserve
that censored result and run a separate extension at degree cap 18 with the
same 2,000-column and 60-second limits. This extension is exploratory and
must not overwrite the original JSON.
