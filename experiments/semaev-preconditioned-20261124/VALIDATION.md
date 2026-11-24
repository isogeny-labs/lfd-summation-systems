# Local validation, 24 November 2026

- `python3 -m pytest -q tests`: **68 passed, 1 skipped**, 20.02 seconds.
  This includes the generic and planted affine-map controls, point/trace
  tests, shared-builder equality, and cold/warm saved-record checks.
- `python3 scripts/check_publication_artifacts.py`: **452 checks passed**.
- `python3 scripts/check_semaev_preconditioned.py`: **45 cases / 675 records
  passed**; independent enumeration of y over V for every oracle, exact
  source/input/protocol hashes, record coverage, phase accounting, complete-pair
  summaries and the secondary interval analysis. No timing replay.
- `python3 scripts/run_semaev_preconditioned_warm.py check`: **108 records
  passed**; held-out inputs, source/protocol hashes, coverage and summaries.
  No timing replay.
- `python3 -m compileall -q csd scripts tests`: passed.
- Plots generated from saved records and visually inspected.

All 135 fresh and 36 warm primary solves assert full equality of their
ordered x-coordinate tuple sets with the exact oracle, not only matching
counts. Both timed baseline formulations use the same fast input builder;
its equation lists are checked against the prior audited construction.

The original cold complete-pair gate remains false because ten split
baseline cases are censored. The separate conservative interval analysis
resolves all 45 best baseline medians and supports the 2x practical target;
the warm protocol independently supports it on every held-out case. These
are finite-grid measured conclusions, not statistical or asymptotic theorems.

No fresh Lean build was run. The proofs are written, not kernel-checked.
No existing manuscript, PDF, Lean source or previous experiment record was
changed by this task.
