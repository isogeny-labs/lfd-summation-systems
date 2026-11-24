# Warm-process confirmation

Added before warm measurements, after observing the fresh-process prior-cohort
gain and some held-out baseline timeouts. This is a confirmation with the
same candidate and inputs, not a new choice of algorithm or favorable cases.

Concern: fresh-process SymPy object initialization can contribute to setup
time, particularly for very small systems. Check that a practically significant
gain also holds in a long-lived process on the challenging held-out cohort.

- Use all 12 frozen held-out cases; no timing-based selection.
- Same primary candidate, shared fast builder, exact output set and backend.
- Methods: original, split, original_preconditioned. Three repetitions each,
  cyclic order so every method occupies every position for a case.
- Import SymPy and warm a trivial Boolean Gröbner solve before timing. Reuse
  one process, while charging each call's actual input construction,
  preprocessing, GB and complete recovery. No case-specific cached answer.
- Same ten-second GB cap. Retain every timeout. Whole-run cap 900 seconds;
  unrun rows remain explicit. No extra subprocess/startup metric or per-method
  RSS claim is meaningful in this shared process.
- Apply the separately documented conservative interval analysis, never
  treating a cap as an exact runtime. Include every case in its geometric
  mean lower bound. Check exact candidate outputs against the frozen oracle.
- Confirmation supports the practical 2x target only if all primary
  repetitions complete correctly and the all-case geometric-mean lower bound
  over the faster GB baseline is at least 2x. Keep the original fresh-process
  complete-pair gate and all failures visible.

Run this after the fresh-process benchmark finishes, so the measurements do
not compete with each other. No implementation tuning on these results.
