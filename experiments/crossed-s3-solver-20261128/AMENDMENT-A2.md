# A2: stop the contended, repetitive-output run

28 November 2026. After the frozen input/source commits and A1, the first
measurement process began. The host simultaneously ran an unrelated msolve
experiment from `graph-transport-recovery-20261125`, so its wall and CPU
numbers are not an isolated comparison. The first two rows (n=3, a=0,
target x=0, tree and crossed, repetition 0) each reached the 20-second cap.
Their verbose stderr repeated input/timing headers 13,082 and 3,077 times,
respectively, producing 51,176,279 and 12,367,017 bytes. The third row's
process was interrupted before completion and has no result record. The
driver process group was terminated; all remaining 79 scheduled rows were
unrun. The two captured stderr files were compressed losslessly with
`gzip -n -9`; `runs.json` retains the original and compressed SHA-256 values.

These two rows are **censored and CPU-contended**, not evidence that one
presentation is faster or slower. The unexpected header loop may indicate
msolve genericity restarts; it needs a separate diagnosis. A later run must
use a clean host and either a diagnosed solver mode or a prospectively
amended output/stopping rule. Do not append a rerun to `runs.json` as if it
were part of this same panel. Preserve the failed attempt and start a
versioned replay record.
