# Crossed S3 witness screen

Frozen before data: 2026-11-28. This is an exhaustive finite-cell correctness
check, not a solver benchmark or a claim of improved relation-finding cost.

## Question

For rational points on a nonsingular short Weierstrass curve, does a conjunction
of two crossed, O-aware three-node S3 trees accept exactly the same four-leaf
x-coordinate tuples and target labels as direct point addition? The first tree
pairs leaves (1,2) and (3,4); the second pairs (1,3) and (2,4). Target O is a
distinct symbol, not x=0. An affine-only version is a diagnostic for missed
relations when an intermediate sum is O.

## Frozen cells and method

Use p=7,11,13 and (a,b)=(1,1),(1,0), six cells in total. Enumerate every
rational affine point, every attained x-coordinate tuple of length four, and
every attained target x-coordinate plus O. For each pair of leaf coordinates,
independently compare the finite roots of S3(x,y,Z) plus the O branch at x=y
with x-coordinates of direct sums of all rational lifts. Compare each tree and
their conjunction with direct four-point addition in every tuple-target cell.
Keep the first affine-only missed relation, all mismatch rows, and exact counts.

The S3 formula is Semaev's 2004 polynomial for y^2=x^3+ax+b. For a finite
target, an O child matches the other child's target x. For target O, two child
labels match iff they are equal, including O/O. In the ordinary case evaluate
S3 on the three finite labels. No timing, F4 or degree reached is measured.

## Decision and limits

Any direct/S3 pair-label mismatch or tree/direct mismatch invalidates the
construction as implemented and must be preserved. Zero mismatches establishes
only these six exact finite cells. The paper proof obligation remains: on a
general field, every genuine relation supplies both trees' internal points,
including O; conversely each tree's x-only witness must be lifted and checked
at point level before its output is used as a relation. A solver-cost study must
then compare one tree, crossed trees, and monolithic S5 on identical targets and
factor bases, charging setup, failures, extraction, verification and linear
algebra. The extra equations alone do not imply lower solver degree.

## Replay

From the repository root, in a clean copy without `result.json`:

```sh
python3 experiments/crossed-s3-witness-20261128/run.py
```

The script refuses to overwrite the first result and records its source and
protocol SHA-256 digests.
