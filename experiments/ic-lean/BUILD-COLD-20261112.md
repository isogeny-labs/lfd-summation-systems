# Cold build, 2026-11-12

Fresh container: elan installed, `lake exe cache get` (Mathlib oleans only), then `lake build` of
every project module from scratch (4 cores, 15 GB RAM).

- `build-cold-20261112-pass1.log`: 1815 of 1823 targets built. The process was killed at
  `ICFallDegrees.GeneratedPrimeUpper` with exit code 137 (out of memory from parallel jobs), not
  by a Lean error.
- `build-cold-20261112-pass2.log`: `lake build ICFallDegrees.GeneratedPrimeUpper` alone (5 min 36 s),
  then `lake build` of everything: "Build completed successfully".
- `audit-20261112.log`: `scripts/verify_lean_fall_degrees.py`: 513 audited declarations,
  no project axioms, no placeholders.
