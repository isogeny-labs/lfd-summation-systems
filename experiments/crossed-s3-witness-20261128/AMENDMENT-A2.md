# Amendment A2: preserve pair-root mismatches

Filed 2026-11-28 after v2 and before the complete v3 rerun. V2 fixed custody
for tree/direct mismatches, but its pair-root equivalence check still asserted
before writing any row. Its six completed cells had zero pair-root mismatches,
so the observed v2 result is unchanged. On a failing cell, however, that
assertion would have lost the diagnostic pair and the already completed cells.

`run-v2.py` is a byte-for-byte copy of v2's source. The amended `run.py` records
every pair-root mismatch, including both direct and polynomial root labels,
and every tree/direct mismatch in a separate `result-v3.json`. It writes the
whole six-cell result before asserting that both mismatch lists are empty.
The output includes SHA-256 digests of the source, unchanged protocol, A1 and
this amendment. V1 and v2 raw records remain separate; v3 decides the final
finite-cell report. No grid cell, relation predicate, or decision threshold
changes.
