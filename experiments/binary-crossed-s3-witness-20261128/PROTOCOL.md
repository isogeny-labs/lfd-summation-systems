# Binary crossed S3 witness transfer gate

Frozen before data: 2026-11-28. This is an exhaustive exact-cell correctness
check, not a solver benchmark. It follows the [odd-characteristic screen](../crossed-s3-witness-20261128/REPORT.md), whose formula and counts do not transfer
automatically to characteristic two.

## Fixed cells and relation predicate

Use `GF(2^3)` modulo `z^3+z+1` and `GF(2^4)` modulo `z^4+z+1`. On each field,
test `E: y^2+xy=x^3+a x^2+b` for `(a,b)=(0,1),(1,1),(0,z),(1,z)`;
`b != 0` makes these ordinary curves nonsingular. Enumerate every rational
affine point, every attained x-coordinate quadruple, and each attained target
x-coordinate plus the identity O. All leaf lifts are included; matching a
target x means matching either sign of the corresponding rational point.

The binary polynomial is
`S3(x,y,z)=(xy+xz+yz)^2+xyz+b` in the chosen field. Compute its finite
roots directly for each attained `(x,y)`, then include the separate O label
when `x=y`. Compare that label set with direct two-point sums. Verify field
inverses, curve membership, identity, inverse and associativity of the
independently written group law before the S3 comparison.

Two balanced trees pair leaves `(1,2)|(3,4)` and `(1,3)|(2,4)`. Each tree
matches a finite target using S3 at its root, with separate O-child rules;
for target O, equal child labels match, including O/O. Compare each tree and
their conjunction against direct four-point addition for every tuple-target
cell. Count an affine-only diagnostic that omits the O rules. Save every
pair-root and tree mismatch row, source and protocol hashes, and exact counts
before asserting zero mismatches.

## Decision and cost boundary

Any field/group check failure invalidates the instrument; any pair-root or
tree/direct mismatch invalidates equivalence on that cell and is preserved.
Zero mismatches proves only the enumerated cells. No solver degree, F4 matrix,
relation yield or complete ECDLP cost is measured. A later cost experiment
must use the same target and factor base across one tree, crossed trees and
monolithic S5, and charge failures and point extraction before a gain claim.

Replay from the repository root in a clean copy without `result.json`:

```sh
python3 experiments/binary-crossed-s3-witness-20261128/run.py
```
