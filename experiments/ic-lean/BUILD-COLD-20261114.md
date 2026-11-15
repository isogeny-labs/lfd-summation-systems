# Cold build of every module, 2026-11-14

Machine: macOS (arm64), Lean 4.19.0 (commit 6caaee842e94), Mathlib `c44e0c8ee63ca166450922a373c7409c5d26b00b`.
`formal/.lake/build` did not exist before the run; only the Mathlib oleans came from the Mathlib cache
(`lake env lean --run .lake/packages/mathlib/Cache/Main.lean get`). Command, from the repository root:

    python3 scripts/lean_check_all.py --out <directory>

Records in `all-modules-20261114/`:

- `build-pass1.log`: `lake build ICFallDegrees.GeneratedPrimeUpper` alone (315 s), then `lake build` of the
  root and all 218 module files then present (797 s). 217 built. One failed:
  `ICFallDegrees.GrowthUnitFrobenius`, with `bad import 'Mathlib.Tactic.Omega'`. That file was a side result
  that nothing imports; its proof also fails once the import is corrected, so it was removed from this
  repository (see `formal/README.md`). No other module reported an error.
- `build-pass2.log`: the same command after the removal and after `ICFallDegrees.lean` and `Audit.lean` were
  extended by the fixed-curve growth theorem: "Build completed successfully" for the root and all 217 modules.
  This pass rebuilt only the two changed files; the other oleans are those of pass 1.
- `axioms.log`: the output of every file that contains `#print axioms` (80 files: `Audit.lean`, the 47 files
  `formal/N8*Audit.lean`, and 32 modules that print their own).
- `summary.json`: 884 answers for 840 distinct declarations. 830 depend on `propext`, `Classical.choice` and
  `Quot.sound`, 8 on `propext` and `Quot.sound`, 2 on `propext`; none on anything else. No `sorryAx`.
  It also lists the SHA-256 of every Lean source file and of the toolchain and manifest files.

The same run then executed `scripts/verify_lean_fall_degrees.py`, which wrote `build.log`, `axioms.log` and
`verification.json` in this directory: 628 declarations of `Audit.lean`, foundational axioms only.

**Addendum, later on 2026-11-14.** Two modules were added after the cold build
(`GrowthTwoTorsionFamily`, `GrowthFermatDiagonalObstruction`), and eight lines to `Audit.lean`. The build step of
the same command was rerun on the same machine on the existing build directory, so it is incremental, not cold:
`all-modules-20261114b/build.log` has "Build completed successfully" for the root and all 219 modules (413 s).
`scripts/verify_lean_fall_degrees.py` then rewrote `build.log`, `axioms.log` and `verification.json` in this
directory: 636 declarations of `Audit.lean`, foundational axioms only; the two new modules print their own axiom
lists with the same result. The audit of every file that prints axioms (`summary.json` above: 840 declarations)
was not regenerated here for the two new modules; the CI run on Linux does it, from an empty build directory.

The build covers the modules added after the cold build of 2026-11-12 (`Growth*`, `Selector*`, `UnitInstance`,
`DepthFourteen`) and the fixed-curve growth modules; for the latter it is the first build in this repository.
The workflow `.github/workflows/ci.yml` repeats it on Linux on every pull request.

## Completed Linux audit and corrected growth-theorem prose

The Linux CI run 43679167384 at `49585b5` passed: all 219 library modules and
848 declarations over 82 files that print axioms, foundational axioms only.
The full logs and original source hashes are retained in `growth-verification-49585b5/`.
The corrected tree differs in one Lean file's comments only, verified against all
270 source hashes from the CI report; its 636-declaration audit was rerun locally.
See that directory's README and `reconciliation.json` for the precise evidence boundary.
