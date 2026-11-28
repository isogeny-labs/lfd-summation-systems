# Amendment A1: preserve mismatch rows on failure

Filed 2026-11-28 after the v1 result and before the complete v2 rerun.

The frozen v1 protocol says to preserve every mismatch. The v1 source collected
mismatch rows but asserted that all rows were empty *before* writing
`result.json`. The completed v1 cells had zero mismatches, so no observed row
was lost, but a future failing replay would exit without its raw mismatch
record. This is a custody defect in the instrument, not a change to the
mathematical question or the fixed grid.

`run-v1.py` is a byte-for-byte copy of the source that produced `result.json`;
that result's `source_sha256` and `protocol_sha256` remain bound to v1. The
amended `run.py` writes all cell records to a separate `result-v2.json` and
then asserts that the mismatch lists are empty. It refuses to overwrite an
existing v2 result, and records the SHA-256 of this amendment alongside the
unchanged protocol and amended source digests. All six cells are rerun from
scratch. V1 and v2 are not pooled; v2 decides the reported finite-cell result.

No solver-degree, cost, relation-yield or population claim is added by this
amendment.
