# Paper reconciliation, 16 November 2026

This paper update starts from evidence snapshot `2cb270c` on pending PR #24,
which includes merged `main` at `82173f7`. The final package is prepared for its CI-gated merge and does not change the numerical engines or Lean sources.

## Results reflected in the manuscript

The combined manuscript retains the solving-degree section and the formal
non-two-torsion function-field and pole-order results. The split-presentation
table distinguishes exact completed cells from a capped computation:

| Summands | Depth | Instances | Result |
|---|---|---|---|
| 4 | 2 | 6 | exact LFD 6 |
| 4 | 3 | 6 | exact LFD 7 |
| 4 | 4 | 2 | exact LFD 7 |
| 5 | 2 | 6 | exact LFD 6 |
| 3 | 5 | 2 | exact LFD 7 |
| 6 | 2 | 2 | exact LFD 7 |
| 5 | 3 | 1 | LFD at least 8; degree-8 run exhausted memory |

The completed records contain closed exclusions below the stated exact
degree and a unit witness at that degree. Their source hashes match the
current engine or its recorded historical version at `efac47b9`. The capped
cell is supported by a log-derived record and has no degree-8 upper bound.
These are record checks; the large split closures were not recomputed here.

The second depth-14 implementation on both unit instances, over p=65537
and p=163841, reports the same
predicted fields as the first lift: input dimension 448,922, low dimension
499,178 with parts 249,824 and 249,354, and full dimension 1,683,218.
The archived record says the input files were verified, its canonical
hashes match the recorded degree-6 run, and the independent closed-space
check passed. The recorded lift and driver hashes match the checked source.

The paper qualifies the conclusion at the point where it is reported.
Both lifts depend on the input space being contained in the true degree-6
closure. The independent closed-space check proves the reverse inclusion,
which supports the lower bound. The second lift does not supply a second
degree-6 producer or remove this dependence. The large stored basis is not
available in the archive for a fresh replay.

## Fresh bounded validation

The arithmetic self-test passed on CPU at four primes. A fresh small-depth
replay compared four cases in 16 lift settings against direct reference
closures and a naive oracle. It found no mismatches; two cases derive 1 and
two do not. Of 18 altered-input controls, 16 were detected; the two unchanged
results drop direct generators that are redundant in those cases.

Reproduction:

```sh
python3 scripts/run_ic_lift_independent.py arith --device cpu --out /tmp/arith.json
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  python3 experiments/ic-bounds-k14-independent/validate.py \
  --sets quick --negative unit-65537:7:4 \
  --work /tmp/lfd-paper-replay --out /tmp/lift-replay.json
```

The fresh outputs are in `review-evidence-20261116/`. The archived larger
validation has 37 cases, 138 lift runs and 69 controls; its six source hashes
match the checked files. That larger suite was not rerun for this update.
No cloud compute was launched.

## Editorial and artifact checks

Restored 21 compatible prose paragraphs from the earlier deslop pass and
edited the new material under the same scientific-writing guidance. The
main paper and submission sheet retain `adam@isogenylabs.io` only.

The earlier prose reconciliation left all 41 theorem, lemma, corollary,
proposition, definition and proof blocks and all 2,925 mathematical
expressions unchanged, ignoring whitespace. This finalization adds the
second prime and its run history without changing any theorem or proof. The new text qualifies recorded results
without promoting the capped cell or the depth-14 lift to a stronger claim.

The rebuilt PDF has 49 pages and no LaTeX warnings, unresolved references
or overfull boxes. All 439 publication consistency checks pass. The
submission sheet records the current page count, abstract length and
artifact SHA-256 values, and the data appendix cites the frozen evidence
snapshot `2cb270c`.

## Finalization for PR #24

Both archived depth-14 records were checked against the current lift and
driver source hashes and against all reported closure dimensions. The
p=163841 restart and the controls not repeated after it are disclosed in
the methods appendix. Build provenance now names Tectonic 0.16.9 and
xdvipdfmx (0.1), replacing the incorrect pdfTeX attribution. The frozen
evidence snapshot 2cb270c contains both replication records.

Two consecutive final-source builds with Tectonic 0.16.9 produced the same
PDF SHA-256: `2c8357fdb2a000c33e246cdb36f2e6bdd95c8cd4dae4657f62fcfc9d9ce59b82`.
The PDF has 49 A4 pages; all 439 artifact consistency checks pass. The
title page, replication paragraph and methods page were rendered and
visually inspected. This validates the package; the large closures were
not rerun.
