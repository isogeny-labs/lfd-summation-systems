# Fresh project build and axiom audit, 14 November 2026

**Passed.** [GitHub Actions run 43648104843](https://github.com/isogeny-labs/lfd-summation-systems/actions/runs/43648104843)
rebuilt all 120 project modules from an empty project build directory, then
ran the complete verifier and axiom audit of 582 declarations. The checkout
was commit `f13c9ed385f0dfb41a982c8494da3b05064c07fd` on a fresh standard
`ubuntu-24.04` runner. The compilation/audit step ran from 22:29:02 to
23:23:16 UTC; these are build execution times, not cryptanalytic timings.

## What was verified

- Lean 4.19.0, compiler commit `6caaee842e94`, x86-64 Linux release.
- Mathlib `c44e0c8ee63ca166450922a373c7409c5d26b00b`; all nine
  dependency revisions match the committed Lake manifest.
- No existing `.lake/build` project directory before compilation. Mathlib's
  published dependency cache was used; compiled project artifacts were not
  restored or reused from an earlier run.
- All 120 project modules compiled in dependency order, one Lake invocation
  per module. All 120 exits and the final verifier exit were zero.
- The fresh report matches all 124 reviewed source/configuration/verifier
  digests and all 582 saved axiom-dependency entries. Only `propext`,
  `Classical.choice` and `Quot.sound` occur. The verifier rejects project
  axioms, placeholders and prohibited execution overrides.
- The executable sparse-certificate demonstration passed.

The existing Lean linter warnings were retained. “Fresh” describes the
project artifacts and environment; it does not mean warning-free source.
No theorem statement, proof source or raw research record was changed to
obtain the successful build. This closes the clean-build gap identified in
the publication review. It does not remove hypotheses from the growth
theorem or supply the missing independent depth-14 lift certificate.

## Retained evidence

| File | Evidence |
|---|---|
| `receipt.json` | Commit, initial empty-build check, exact dependency revisions, source digests, ordered modules, argv, per-stage exit codes and elapsed times |
| `000-*.log` through `119-*.log` | Fresh module-build logs, including the new selector, growth, unit-instance and depth-fourteen modules |
| `verification.json`, `axioms.log` | Fresh verifier report and complete axiom audit |
| `build.log`, `sparse-demo.log`, `verifier.log` | Final aggregate build, sparse demonstration and verifier output |
| `dependency-setup.log`, `mathlib-cache.log` | Pinned dependency and Mathlib-cache setup |
| `driver.log` | Ordered build progress and final successful summary |
| `artifact-manifest.json` | SHA-256 of all 129 files fetched from the CI artifact, without changing their contents |

Downloaded artifact: `lean-clean-f13c9ed385f0dfb41a982c8494da3b05064c07fd`,
GitHub artifact ID `14137215246`. The files are committed here so their
availability does not depend on the CI artifact's 30-day retention window.
The earlier warm audit and 513-declaration cold-build records are preserved.

## Reproduction and method sources

The workflow is `.github/workflows/lean-clean.yml`. It installs the toolchain
named in `formal/lean-toolchain`, fetches the pinned dependencies, checks that
the manifest is unchanged, and obtains the Mathlib cache. In a fresh checkout
after that setup, run:

```sh
python3 scripts/clean_lean_audit.py --out experiments/ic-lean/clean-20261114
```

The script refuses an existing project build directory and deletes nothing.
Its receipt and module logs distinguish actual fresh compilation from the
final verifier's aggregate up-to-date build. The final verifier writes its
normal outputs in the disposable checkout; the checked-in original records
in this evidence checkout were left intact.

The method follows `formal/README.md` and the earlier
`experiments/ic-lean/BUILD-COLD-20261112.md`, with sequential compilation to
avoid the earlier parallel-memory failure. The primary implementation
references checked were [Mathlib 4.19 cache commands](https://github.com/leanprover-community/mathlib4/blob/v4.19.0/Cache/Main.lean)
and [elan toolchain selection](https://github.com/leanprover/elan/blob/master/README.md).
No new mathematical or performance claim is made.

For a read-only check of the retained evidence and its source binding, use:

```sh
python3 scripts/check_publication_artifacts.py
```
