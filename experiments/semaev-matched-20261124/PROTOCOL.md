# Matched binary Semaev comparison

Frozen before measured runs, 24 November 2026. Local CPU only; no rented compute.

## Correctness before timing

Audit raw S4 zeros, projections of the affine two-S3 split, and actual ordered
x-coordinate triples of rational-point decompositions. Enumerate all a, nonzero b
and affine rational targets over GF(4) and GF(8); also a in {0,1}, b in {1,3}
over GF(16). Use the whole field for V in this audit. Enumerate rational points
by scanning both coordinates and use the curve group law, independently of
Semaev evaluation. Record missing/extra tuples and first witnesses. Classify
degree drops, quadratic gcds and the infinity branch x1=x2, x3=x(R).

For the actual PDP comparison, both formulations receive the same constraints
that each xi lifts to a rational point. The split solve is the union of the
finite-u solve and the explicitly recovered infinity branch. Compare each
result with the independently enumerated rational-point decompositions.

## Frozen grid and selection

Grid (n,ell): (3,1), (3,2), (4,2), (5,2), (5,3), (7,2).
Seeds: 0,1,2. Curve a is one seeded bit; b is seeded uniform nonzero.
V is the polynomial-basis span of 1,w,...,w^(ell-1). Seeded shuffle of all
affine targets, select the first with a PDP solution and the first without;
record a missing category rather than replacing it with a favorable target.
Keep the same selected inputs, field basis, coordinate ordering and constraints
for both formulations. Repeated summands and targets in the factor base are
allowed; output is the set of ordered x-coordinate triples, not signed point
triples or an index-calculus relation count.

## Solver and cost accounting

SymPy 1.14.0, groebner(method='buchberger', order='grevlex', modulus=2),
ordinary polynomial ring with all Boolean equations. Original variables:
x1 block, x2 block, x3 block. Split: same blocks followed by u block.
This is a small-field Buchberger experiment, not a reproduction of Magma F4.

Three fresh-process repetitions per formulation and instance; alternate which
formulation runs first. Primary metric: time inside the GB call. Also report
input construction, FGLM conversion plus complete Boolean solution recovery,
infinity-branch recovery, verification, process wall time and peak process RSS
(runtime baseline included). Charge membership-polynomial construction to both.
Use a 5-second GB-call limit and 20-second process limit. Keep timeouts/errors
in records; no numerical speedup is assigned to an unfinished pair. Cap total
benchmark wall time at 1200 seconds and explicitly report any unrun rows.

A separate diagnostic run wraps SymPy's S-polynomial routine, recording its
call count, largest returned S-polynomial degree and term count. Instrumented
runs do not enter timing summaries. No Macaulay matrices are constructed by
this solver: actual matrix dimensions, matrix ranks and matrix sparsity are
unavailable, recorded as null. Observed intermediate degree is not called a
solving degree or degree of regularity.

Report per-case medians and every repetition, paired wins/losses/ties, solution
strata, parameter strata, failures, and a descriptive crossover if the grid
supports one. Do not infer an asymptotic advantage or a large-field crossover.
