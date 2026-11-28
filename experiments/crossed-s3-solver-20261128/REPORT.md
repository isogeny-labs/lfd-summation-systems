# Matched crossed-S3 solver panel: status

The [protocol](PROTOCOL.md) and generator were committed before data. The
generator produced 27 inputs for 9 selected finite-target queries on six
binary curves (n=3,4,5; a=0,1). Every generated S3 ANF passed 258
deterministic arithmetic evaluations per field. Exhaustive field-level
projection checks established that each O-aware single tree and crossed tree
has the same ordered factor-base-abscissa relation set as direct four-point
addition on all nine queries. The exact Boolean eliminant was generated from
that independent group-addition predicate. The n=3 curves had no selected
zero-relation target; the manifest records all selected queries.
An [independent Sage point replay](sage-replay.json), with source committed
before its result, matched the full ordered tuple set on all nine queries,
including the n=5 cells. It used Sage finite fields and elliptic-curve point
addition rather than the generator's point arithmetic.
The generated msolve syntax was also evaluated directly on 1,503 positive
tree/crossed/direct witness assignments. All 29,184 Boolean leaf assignments
for the direct eliminants matched the expected relation predicate; see
[`input-verification.json`](input-verification.json).

The original full-solve measurement is incomplete. [A1](AMENDMENT-A1.md) records an
unsupported macOS address-space cap before any msolve run. [A2](AMENDMENT-A2.md)
records two 20-second censored, CPU-contended rows and why the panel was
stopped. The remaining 79 scheduled full-solve rows have no measurement.

The repeated headers came from msolve's default genericity retries during
rational parametrization. [A3](AMENDMENT-A3.md) freezes a separate
`-g 2 -c 0` GB-stage replay, which completed a read-only n=3 diagnostic
without the loop. The approved remote box timed out on SSH. The local host
then had a clear window, and all **81 frozen GB-stage runs** completed (nine
queries, three arms, three fresh-process repetitions). The next unrelated
batch started after the last row. All raw hashes, input hashes and within-arm
F4 geometry were checked by [`analyze_gb.py`](analyze_gb.py); the derived
records are in [`analysis-gb.json`](analysis-gb.json).

## GB-stage results

Solver: msolve 0.11.0, one thread, seed 12345, grevlex, `-g 2 -c 0 -v 2`.
The table gives median process wall time over three runs, maximum logged F4
degree, and crossed/single ratios of cumulative F4 matrix cells. All arms
used identical targets and factor bases. `D` is the direct Boolean eliminant,
whose construction used exhaustive point addition.

| n, a, target x | Relations | Tree s | Crossed s | D s | F4 degree tree / crossed | Crossed/tree matrix cells |
|---|---:|---:|---:|---:|---:|---:|
| 3, 0, 0 | 7 | 0.008 | 0.013 | 0.007 | 4 / 4 | 4.93× |
| 3, 1, 0 | 16 | 0.011 | 0.036 | 0.008 | 4 / 4 | 10.12× |
| 4, 0, 0 | 107 | 0.019 | 0.078 | 0.013 | 6 / 5 | 10.36× |
| 4, 0, 8 | 0 | 0.013 | 0.036 | 0.007 | 4 / 4 | 4.78× |
| 4, 1, 0 | 107 | 0.019 | 0.078 | 0.013 | 6 / 5 | 10.36× |
| 4, 1, 8 | 0 | 0.013 | 0.036 | 0.008 | 4 / 4 | 4.78× |
| 5, 0, 0 | 263 | 0.528 | 15.554 | 0.116 | 5 / 5 | 8.36× |
| 5, 1, 0 | 0 | 0.014 | 0.045 | 0.007 | 4 / 4 | 4.29× |
| 5, 1, 8 | 1 | 0.016 | 0.047 | 0.007 | 4 / 4 | 5.08× |

Crossed trees were slower in **9/9** matched GB queries. The median
crossed/single wall ratio was **3.12×**, median child-CPU ratio **3.35×**,
and median cumulative-matrix-cell ratio **5.08×**. In the two cells where
crossing lowered maximum logged F4 degree from 6 to 5, it increased matrix
cells by **10.36×** and wall time by about **4.2×**. The n=5 positive cell
was especially costly: crossed/single median wall ratio **29.48×**.

This is a concrete negative for the crossed presentation's GB stage on this
frozen range. The direct eliminant was fastest at the GB stage in 9/9, but
its exhaustive construction cost has not been charged to a cold query.
Maximum logged F4 degree is an observed solver metric, not a proof of degree
of regularity. The GB-only command did not extract point relations. Input
construction, exceptional-branch recovery, solution extraction, verification,
relation failures across queries, linear algebra, and an identical-target
rho comparison remain outside this measurement. A total-cost advantage is
therefore unestablished; the crossed candidate has a strong negative stage
result and should not be promoted without a materially different encoding.
