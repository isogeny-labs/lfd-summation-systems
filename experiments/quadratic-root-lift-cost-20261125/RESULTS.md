# Bounded quadratic root-lift solver result

The frozen panel in `manifest.json` contains 15 two-summand queries: four
parent-campaign targets and one diagonal target at each of (n,k)=(23,11),
(29,13), and (29,14). The exact factor bases and parent manifest are hashed.
`run.py` encodes the quadratic root-membership equations into CNF, calls
Kissat until UNSAT, recovers signed points, and checks the group law. Each
arm and target runs in a fresh process; the outer charged wall and CPU totals
include process startup. The 10 s SAT cap was never reached.

| Cell | Complete signed relation counts across the five targets | SAT total wall / CPU (s) | Best native arm total wall / CPU (s) | Median SAT query wall (ms) | Median native LRF query wall (ms) |
|---|---|---:|---:|---:|---:|
| (23,11) | 2, 0, 1, 1, 1 | 0.727 / 0.690 | image 0.468 / 0.452 | 45.0 | 0.147 |
| (29,13) | 0, 0, 2, 1, 1 | 0.761 / 0.737 | image 0.489 / 0.471 | 34.3 | 0.182 |
| (29,14) | 0, 0, 2, 1, 1 | 0.984 / 0.959 | brute 0.514 / 0.498 | 47.9 | 0.352 |

All **60 receipts completed**. For each target the SAT, native LRF,
conditional-image LRF, and direct native enumeration arms returned the same
complete canonical signed pair set: zero parity mismatches. The SAT arm
enumerated 23 Boolean models, which deduplicate to 13 signed pairs across all
15 targets. The three diagonal queries each returned the repeated point pair.
The SAT arm had greater charged wall and CPU cost on every cell, and its
median query wall was over 34 ms while native LRF was below 0.36 ms. Median
process peak RSS for SAT was 38.0, 40.5 and 42.9 MB across the three cells;
the corresponding native LRF medians were 37.5, 38.8 and 40.3 MB.

**Decision.** The exact quadratic lift passes relation correctness but fails
the matched query-cost gate on this frozen panel. The Boolean input equations
have degree at most two; Kissat exposes no polynomial solving-degree trace,
so these results do not establish a lower maximal solving degree. They also
do not establish a new end-to-end scalar-solving method. The native two-point
query is already fast on these cells, and the m=3 method still requires an
outer factor-base scan. A different algebraic solver or relation family would
need its own complete cost and degree measurements before a larger campaign.

Reproduce with `python3 experiments/quadratic-root-lift-cost-20261125/make_manifest.py`,
`python3 experiments/quadratic-root-lift-cost-20261125/run.py`, and
`python3 experiments/quadratic-root-lift-cost-20261125/verify.py`. The raw
receipts, source hashes, Kissat version, host environment and exact pairs are
in `results.json`. The verifier
reconstructs the curve and checks every signed pair and complete-set parity.
This is one run per case; the observed order is decisive within this panel,
but no large-field scaling claim follows.
