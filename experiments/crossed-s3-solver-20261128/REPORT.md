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

The solver measurement is incomplete. [A1](AMENDMENT-A1.md) records an
unsupported macOS address-space cap before any msolve run. [A2](AMENDMENT-A2.md)
records two 20-second censored, CPU-contended rows and why the panel was
stopped. The remaining 79 scheduled rows have no measurement. The existing
records cannot rank solver cost or show a reduction in solving degree.
The direct enumeration and output-extraction costs have not yet been timed.

The repeated headers came from msolve's default genericity retries during
rational parametrization. [A3](AMENDMENT-A3.md) freezes a separate
`-g 2 -c 0` GB-stage replay, which completed a read-only n=3 diagnostic
without the loop. Its timings remain unrun because other msolve workloads
are active on the local host; the approved remote box timed out on SSH.
The next measurement needs an isolated host, followed by a separate solution
extraction and complete relation-query comparison. Only then can a
larger-field and full relation-collection comparison be considered.
