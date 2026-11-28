# Crossed S3 witnesses: exact finite-field screen

## Result and scope

The [frozen protocol](PROTOCOL.md) and [v1 source](run-v1.py) were committed at
`d6f03e4` before enumeration. The [v1 raw result](result.json) was committed
separately at `dcc2c53`. [Amendment A1](AMENDMENT.md) and the revised
[v2 source](run-v2.py) were committed at `e1ef8d0` before a complete repeat;
the [v2 raw result](result-v2.json) was committed at `27f44fb`. V2 still
asserted on a pair-root mismatch before writing its record. [Amendment A2](AMENDMENT-A2.md)
and the final [source](run.py) were committed at `cc3c57e` before the complete
v3 rerun; its [raw result](result-v3.json) was committed at `9bae4c7` and
decides this report. All three raw files have identical original per-cell
fields; v3 additionally has empty pair-root mismatch lists. Every source,
protocol and applicable amendment SHA-256 field matches the corresponding
file. V3 writes both kinds of failures before its final assertion.
On six nonsingular prime-field curves, v3 enumerated every
attained rational x-coordinate quadruple and every attained target x-coordinate
plus O: **270,910 tuple-target cells**, including **75,577** with a direct
four-point relation. Each O-aware tree and their conjunction agreed with direct
group addition in every cell. The S3 finite-root-plus-O label set also agreed
with direct two-point sums for every ordered pair of attained x-coordinates.

| p | (a,b) | group order | attained x | tuple-target cells | direct relations | O-aware mismatches | affine-only missed |
|---:|:---:|---:|---:|---:|---:|---:|---:|
| 7 | (1,1) | 5 | 2 | 48 | 48 | 0 | 18 |
| 7 | (1,0) | 8 | 4 | 1,280 | 583 | 0 | 235 |
| 11 | (1,1) | 14 | 7 | 19,208 | 8,089 | 0 | 1,711 |
| 11 | (1,0) | 12 | 6 | 9,072 | 3,911 | 0 | 975 |
| 13 | (1,1) | 18 | 9 | 65,610 | 25,905 | 0 | 4,323 |
| 13 | (1,0) | 20 | 11 | 175,692 | 37,041 | 0 | 5,997 |

The 13,259 affine-only misses include target O and cancelling partial sums;
the frozen record does not separate those causes. They demonstrate that an
affine-only instrument is not complete on this grid. The exact agreement is a
finite-cell correctness check, not a solver-degree or cost measurement.

## Construction and proof obligation

Let `x_i=x(P_i)` for four rational factor-base points and let `r` be the
target x-coordinate. Tree A uses

`S3(x1,x2,u)=S3(x3,x4,v)=S3(u,v,r)=0`;

tree B replaces the pairings by `(x1,x3,w)` and `(x2,x4,z)` and uses
`S3(w,z,r)=0`. An internal variable also has a distinct O value when its
children can cancel. At the root, O plus a finite label matches that label;
two finite labels sum to O exactly when they are equal. This is an explicit
disjunction, not the substitution `x(O)=0`.

Every point relation supplies witnesses for both trees by taking the x of its
partial sums, or O. Conversely, for rationally liftable leaves, a finite S3
edge has a signed point-addition witness. At a tree node, flipping every sign
in one child subtree changes the sign of its sum without changing any x
labels, so the edge witnesses can be oriented consistently from the root
downward. The O cases obey the same induction with the explicit identity
rules. Thus a single complete tree and the conjunction have the same
*projected* signed leaf relation set. A solver output still needs rational
point lifts and a direct point-level check of the target relation. The exact
enumeration checks this argument on the frozen finite cells; it is not a
formal proof over every field.

The conjunction cannot increase relation yield at a fixed target and factor
base. Its only proposed benefit is an algebraic presentation effect: the
shared leaf variables give a more overdetermined system. With four leaf
coordinates in a binary subspace of dimension `k`, one balanced S3 tree has
`4k+2n` Boolean variables and `3n` descended equations; two crossed trees
have `4k+4n` variables and `6n` equations. At `4k=n`, these are `3n` variables
versus `5n` variables, respectively. All equations retain the S3 degree;
the extra variables and equations may raise or lower actual F4 cost.

This is separate from the established single-chain degree reduction in
[Semaev's work](https://eprint.iacr.org/2004/031.pdf) and
[Karabina's auxiliary-variable construction](https://eprint.iacr.org/2015/319.pdf).
The existing binary-tree topology record in the linked research program
explicitly does not cover merge DAGs. That scope gap is a reason to test the
presentation, not evidence of an advantage.

## Deciding cost test

Freeze one dataset of identical curves, targets, rational factor bases and
summand order. Compare monolithic S5, one complete balanced S3 tree, and the
crossed O-aware network. First prove equal extracted point relations and
retain failed or censored solves. Then count degree reached, Macaulay/F4
matrix cells, field operations, retries, extraction and verification for all
arms. Finally charge setup, relation collection, useful rank and linear
algebra on the same fixed target, beside a matched generic-group baseline.
The crossed network is falsified as a cost lead if its added variables and
equations outweigh any solver-degree benefit after those charges.

## Replay

Run `python3 experiments/crossed-s3-witness-20261128/run.py` in a clean copy
of this branch without `result-v3.json`. The program refuses to overwrite the
v3 record, verifies pair-root sets independently by group addition, and writes
both kinds of mismatch if any appears before its final assertion. V1 and v2
replays use `run-v1.py` and `run-v2.py`, respectively, without their matching
result files in separate clean copies.
