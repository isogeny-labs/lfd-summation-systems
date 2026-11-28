# Binary crossed S3 witness screen

## Exact result

The [protocol](PROTOCOL.md) and [source](run.py) were frozen and committed at
`0c50c26` before any data. The [raw result](result.json) was committed
separately at `26a2a82`; its source and protocol SHA-256 fields match the
files. Eight ordinary binary curves over GF(8) and GF(16) were checked. The
enumeration covered every attained rational leaf x-coordinate quadruple and
every attained target x-coordinate plus O: **171,742 tuple-target cells**,
including **70,180** with a direct four-point relation. The finite S3 roots
plus explicit O branch agreed with direct pair sums in every one of 350 pair
checks. Each O-aware tree and their conjunction agreed with direct point
addition in every tuple-target cell. There were **zero pair-root and zero
tree mismatches**.

| Field | a | b | Group order | Attained x | Cells | Direct relations | Affine-only misses |
|:---|---:|---:|---:|---:|---:|---:|---:|
| GF(8) | 0 | 1 | 4 | 2 | 48 | 23 | 17 |
| GF(8) | 1 | 1 | 14 | 7 | 19,208 | 8,089 | 1,711 |
| GF(8) | 0 | z | 8 | 4 | 1,280 | 583 | 235 |
| GF(8) | 1 | z | 10 | 5 | 3,750 | 1,681 | 485 |
| GF(16) | 0 | 1 | 16 | 8 | 36,864 | 14,951 | 2,927 |
| GF(16) | 1 | 1 | 16 | 8 | 36,864 | 14,951 | 2,927 |
| GF(16) | 0 | z | 16 | 8 | 36,864 | 14,951 | 2,927 |
| GF(16) | 1 | z | 16 | 8 | 36,864 | 14,951 | 2,927 |

The affine-only form missed **14,156** direct relations in aggregate; the
counter mixes O targets with cancelling internal sums. It shows why the
identity cannot be encoded as field coordinate zero. The four GF(16) rows
have the same counts and should not be treated as four independent samples
from a population of curves. This is exact enumeration of the specified
finite cells, with no statistical inference.

## Instrument and correctness boundary

The source implements GF(2^n) multiplication independently of the S3
formula. It verifies every nonzero field inverse, every point and sum on the
curve, identity and inverse laws, and **20,704** associativity triples before
the relation check. The S3 evaluator agrees with the existing repository
formula in `experiments/ic-curve-validation/code/run_checks.py`, but point
addition is implemented separately. Both pair-root and tree mismatch lists
are written to the raw file before the final assertion, so a failing cell
remains inspectable.

The odd-characteristic [companion screen](../crossed-s3-witness-20261128/REPORT.md)
checks the same relation predicate with a different curve law. Together they
support the O-aware witness construction on these finite cells. A general
equivalence proof still needs to cover all fields, rational lifts, signs,
partial sums at O, and the exact factor-base membership rule. The crossed
network has the same projected relation set as a complete single tree; its
possible benefit is solver presentation cost, not extra relation yield.

## Next deciding test

Freeze identical binary curves, factor bases, targets and extraction rules
for monolithic S5, one O-aware balanced S3 tree, and the crossed network.
Verify every extracted relation by point addition, count all failed and
censored attempts, and compare degree reached, F4 matrix cells, field
operations, memory, extraction and verification. Carry any solver-stage gain
through relation collection, independent rank, linear algebra and a matched
same-target generic baseline. An equation-count increase or a lower first
fall degree by itself does not meet that cost test.

## Replay

From the repository root in a clean copy without `result.json`:

```sh
python3 experiments/binary-crossed-s3-witness-20261128/run.py
```

The script refuses to overwrite the first result and writes source and
protocol SHA-256 digests into the raw record.
