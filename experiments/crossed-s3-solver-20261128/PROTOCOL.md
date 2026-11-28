# Matched crossed-S3 solver-cost panel

Frozen before generating solver inputs or collecting measurements, 28 November
2026. This is a bounded solver and relation-query comparison, not a complete
discrete-log measurement. The binary crossed-witness correctness screen is
`../binary-crossed-s3-witness-20261128/`.

## Fixed instances

Fields are GF(2^3), GF(2^4), GF(2^5) in polynomial bases with respective
moduli `0b1011`, `0b10011`, `0b100101`. The curves are
`y^2+xy=x^3+a*x^2+1`, with `a=0,1` in each field. The factor-base abscissas
are the rationally liftable elements of the span of the first `k` basis
vectors, with `(n,k)=(3,2),(4,3),(5,3)`. Both signs of each abscissa are
admissible. For each curve, inspect affine rational target abscissas in
ascending integer order. Select the first whose four-summand abscissa
relation set is nonempty and the first whose set is empty; record a missing
category if necessary. The target point's sign is irrelevant to an abscissa
query but full-point verification must check both possible target signs.

## Three presentations on identical queries

1. One balanced S3 tree, pairing `(x1,x2)` and `(x3,x4)`.
2. Crossed S3 trees, adding pairings `(x1,x3)` and `(x2,x4)` and a second
   root equation. The two trees have distinct auxiliary coordinates.
3. A direct Boolean eliminant for the exact four-point relation predicate on
   the same factor base. It is generated from exhaustive independent group
   addition, and is an empirical small-field baseline, **not** raw S5.

The first two arms encode each auxiliary point as an abscissa plus an infinity
flag. Flags are canonical: the abscissa bits are zero when the flag is set.
Each pair node has the affine `S3` equation when finite and the equality
`x_i=x_j` when infinity. The root node uses the affine `S3` equation when both
children are finite; when exactly one is infinity the other must equal the
target abscissa; two infinity children do not match a finite target. All
constraints are descended to Boolean ANF over GF(2), with Boolean equations
for every variable. Rational liftability of every leaf is imposed in all
arms. No particular sign assignment is fixed inside the x-only equations.
For the direct eliminant, the zero set is exactly the ordered abscissa tuples
for which some rational point lifts sum to a point with the target abscissa.

Before timing, verify every generated ANF equation against its field equation
on a deterministic evaluation sample and verify that the three presentations
have the same projected relation set by exhaustive field-level witness
enumeration on the fixed instances. If this equivalence fails, stop timing
that cell, record the counterexample, and correct the model only under a
versioned protocol amendment before rerun. Independently evaluate full-point
witnesses for each positive cell.

## Measurement and stopping rule

Use `msolve` 0.11.0, characteristic 2, one thread, fixed random seed 12345,
grevlex F4, verbose level 2. Run each arm in a fresh process, with three
repetitions per query, rotating arm order. Cap each process at 20 seconds and
2 GiB virtual memory. Record all timeouts, errors, and skipped cells. Charge
input generation, solver wall and child CPU, solution extraction, O-branch
handling, full-point verification, and process overhead. Record maximum
logged F4 degree, number of rounds, cumulative and largest matrix cells,
rows reduced, zero reductions, ideal/output class, and peak memory where
available. Preserve exact `.ms` inputs, raw stdout/stderr, hashes, versions,
curve/target/base identity, and per-run JSON records. Compare only completed
matched cells and report censoring separately; lower F4 degree alone is not
a solver-cost gain. The direct finite enumeration is timed as a separate
matched baseline, with its setup charged.

This panel can falsify a solver advantage on the fixed small fields. Any
positive result would require repetition at larger fields and a cold,
same-target relation collection and linear algebra comparison before a total
ECDLP cost claim.
