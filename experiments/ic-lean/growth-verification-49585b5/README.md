# Completed verification of the growth theorem

The Linux CI run [43679167384](https://github.com/isogeny-labs/lfd-summation-systems/actions/runs/43679167384)
passed at commit `49585b5ff70add29f0a805fd496f01caf66379e8`, using Lean 4.19.0 and pinned Mathlib:
219 library modules, 82 files printing axiom dependencies, and 848 distinct audited declarations,
with only `propext`, `Classical.choice` and `Quot.sound`. `build.log`, `axioms.log` and `summary.json`
are copied without alteration from the run's `lean-audit` artifact.

The follow-up correction concerns an inference in the prose: a proper ideal prevents the unit-ideal
argument from proving an LFD lower bound, but does not refute that numerical bound.
`reconciliation.json` compares all 270 source hashes from the CI summary against the corrected tree.
Only `GrowthFermatDiagonalObstruction.lean` differs, in comments; its declarations and proofs are unchanged.
The original CI source hashes are retained, rather than rebound to a different commit.

The corrected sources were also built and audited locally by `scripts/verify_lean_fall_degrees.py`:
636 declarations, foundational axioms only. The resulting `../verification.json`, `../axioms.log`
and `../build.log` bind that local audit to the corrected sources. A redundant local all-module rerun
was stopped after the Linux full audit succeeded; it is not used as completion evidence.

This verifies the Lean development and its axiom dependencies. The general growth theorem for
genus-one function fields still rests on the written proof and its reviews.
