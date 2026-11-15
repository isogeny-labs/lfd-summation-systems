# Last fall degrees of elliptic curve summation systems: code and data

Companion repository of the paper

> Adam Buran, *Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values*.

It holds the Lean development, the closure engines, the independent checker, the written proofs and every
record behind the numbers in the paper. The paper is `paper/curve-lfd-bounds.pdf`. See [publication status and build instructions](paper/README.md) before using a submission package.

## Layout

| Path | Content |
|---|---|
| `paper/` | LaTeX source, figures and PDF of the paper; `binary-lfd.*` is the companion manuscript it cites for the binary bounds (its own records are not in this repository, as its abstract says); `eprint-curve-lfd-bounds/` is the ePrint submission package |
| `formal/` | Lean 4 project `ICFallDegrees` (Lean 4.19.0, Mathlib pinned in `lake-manifest.json`) |
| `csd/` | closure engines and certificate code (Python) |
| `scripts/` | drivers: one script per experiment, the figure script and the paper build |
| `tests/` | unit tests of the certificate code (`pytest tests`) |
| `docs/proofs/` | written proofs that are not formalized (the depth reduction, genericity, the growth theorem `N8-growth.md`, the block-of-four lemmas `N9-block-four.md`), the plans of the Lean formalizations, and the notes behind the kernel-checked fixed-curve growth theorem (the other `N8-*.md`; start with `N8-two-torsion-rational-route.md`) |
| `docs/reviews-*.md`, `docs/reviews-20261112/`, `docs/proofs/N8-growth-review*.md` | reviews of the written proofs, the red-team reviews of the growth proof and the audit of the paper (see *Working material* below) |
| `docs/IC_FALL_COUNTEREXAMPLES.md`, `docs/IC_FALL_PROOFS.md` | notes on the fall counterexamples, checked by `scripts/verify_ic_fall_counterexamples.py` against `experiments/ic-fall-counterexamples/` |
| `experiments/ic-lean/` | audit of the Lean development (declarations and their axioms) |
| `experiments/ic-curve-sweep/` | the sweep of Sections 5 and 6 |
| `experiments/ic-chain-depth8/`, `-depth9/`, `-depth10/` | two-summand systems at depths 8, 9 and 10 |
| `experiments/ic-depth-step/` | dim V_3 and dim V_4 against the ideal, depths 3 to 10 |
| `experiments/ic-curve-paper/` | record of the paper build |
| `experiments/ic-bounds-d5/`, `ic-bounds-exact10/`, `ic-bounds-fall/` | independent engine at depths 10 to 13, the exact value at depth 10, closure dimensions |
| `experiments/ic-bounds-k14/` | GPU engine, degree lift and independent checker, depths 10 to 16 |
| `experiments/ic-bounds-freeext/` | free extension, and the predictions written before the depth-14 runs |
| `experiments/ic-bounds-m3/`, `ic-bounds-m3-depth/` | three summands |
| `experiments/ic-bounds-msum/` | three to eight summands |
| `experiments/ic-bounds-msummands/` | three to five summands with the second implementation (`csd/ic_bounds_msum.py`), local and pod records, predictions |
| `experiments/ic-bounds-upper/`, `ic-bounds-upper-redteam/` | the upper bound `k+2` for two summands, and its independent red-team checks |
| `experiments/ic-bounds-growth/`, `ic-bounds-block4/` | records behind the growth theorem and the block-of-four lemmas |
| `experiments/ic-bounds-genericity/` | genericity of the depth-7 base |
| `experiments/ic-bounds-locus/` | exceptional locus at depth 6 |
| `experiments/ic-bounds-binary/`, `ic-bounds-binary-law/` | binary and Koblitz curves |
| `experiments/ic-bounds-binary-deep/` | deeper binary and Koblitz descents, with two independent replications |
| `experiments/ic-bounds-cert/`, `ic-bounds-exact-k6/` | data of the certificates checked in Lean |
| `experiments/ic-curve-validation/` | independent implementation of the closure in the ordinary ring |

Each experiment directory has a `README.md` (or `RESULTS.md`, `VERIFY.md`, `NOTES.md`, `REPORT.md`) that names the
command behind every record; `experiments/publication-review-20261114/` holds the check records of
`docs/reviews-publication-20261114.md`. Appendix C of the manuscript identifies the evidence used.

## Requirements

Python 3.10 or later (CI uses 3.11) with `numpy`; `matplotlib` for the figures and `sympy` for a few symbolic checks
(`requirements.txt`). The runs at depths 11 to 14 used a GPU through `cupy`; nothing else needs one. The paper
and `binary-lfd.tex` build with pdfLaTeX or Tectonic. The Lean project needs Lean 4.19.0 (see `formal/README.md`).

## Quick checks

    python3 scripts/ic_bounds_msum.py selftest      # three to eight summands: must print SELFTEST PASSED
    python3 scripts/plot_ic_lfd_figures.py          # regenerates three figures from the records
    python3 scripts/build_ic_curve_paper.py         # regenerates tables, all manuscript figures and PDF
    python3 scripts/check_publication_artifacts.py  # checks saved hashes, audit, sweep and submission abstract
    for s in scripts/check_n8_*.py; do python3 $s; done   # 33 exact regressions behind the fixed-curve growth theorem

## Reproducing the main computations

    python scripts/verify_lean_fall_degrees.py --lake /path/to/lake
    python scripts/run_ic_curve_sweep.py --check-only --workers 4
    python scripts/replay_with_current_sources.py experiments/ic-chain-depth8 -- --workers 3
    python scripts/run_ic_chain_depth8.py --check-only --depths 9 --primes 7681 \
      --out experiments/ic-chain-depth9 --workers 1 --product-cap 12000 \
      --shift-max-cells 450000000
    python scripts/run_ic_chain_depth10_probe.py --tag rerun
    python scripts/run_ic_chain_depth10_degree6.py --combos 3000
    python scripts/run_ic_depth_step.py
    python scripts/analyze_ic_base_classes.py
    python scripts/ic_bounds_msum.py run --max-log2-points 24

The first command builds the Lean project and prints the axiom dependencies of every audited declaration; run
`lake exe cache get` in `formal/` first to use Mathlib's prebuilt files, otherwise Mathlib is compiled from
source. The continuous-integration workflow (`.github/workflows/ci.yml`) runs that build and audit, and the
Python self-tests, on every pull request and every push to `main`. Three sources of the depth-8 run were edited after it
(`experiments/ic-chain-depth8/NOTES.md`), so its record no longer passes the driver's source-hash check;
`replay_with_current_sources.py` keeps the data hashes, reports the changed sources and replays every system
with the current code. `python3 scripts/lean_check_all.py --out <dir>` goes further than the first command: it
builds every module under `formal/ICFallDegrees/`, including those the root file does not import, and audits
every `#print axioms` in the project (`formal/README.md`); the workflow runs it too. The
depth-ten runs take about two hours each and the degree-six closure about five hours. The commands for the GPU
runs and for the degree lift are in `experiments/ic-bounds-k14/README.md` and `LIFT.md`.

## What kind of evidence each result has

The paper uses three terms. *Kernel-checked*: the proof is verified by the Lean kernel. *Written*: there is a
proof on paper (`docs/proofs/`) that is not formalized. *Exact*: the outcome of a computation in exact arithmetic
that ends with a certificate. The value 7 at depth 14 rests on the degree lift of
`experiments/ic-bounds-k14/LIFT.md`, which has a written proof and passed a control on a curve with solutions but
has no independent certificate.

## Working material

Besides the paper, its proofs and its records, the repository keeps the working material of the project, so
that readers can see how the results were checked:

- `docs/reviews-*.md`, `docs/reviews-20261112/` and `docs/proofs/N8-growth-review*.md` are internal reviews: red-team reviews of the written proofs,
  referee-style reports and an audit of the paper. They record findings and how each was resolved. They are not
  part of the paper and have not been externally refereed.
- `experiments/ic-bounds-k14/POD.md`, `RELAUNCH.md` and `relaunch/`, `experiments/ic-bounds-binary-deep/pod/` and
  `scripts/ic_bounds_big_pod.py` are runbooks and launch scripts for the rented GPU machines (Runpod) used at
  depths 11 to 16. They contain no credentials; the control token of a pod lives only on the pod and in the
  local state directory (`$POD_STATE_DIR`).
- `docs/proofs/*-lean-plan.md` are plans of formalizations; the Lean sources are the result.
- `docs/IC_BOUNDS_PROGRAM.md` is the working plan of the bounds programme (dated 2026-11-08, with status
  sections up to 2026-11-12): claims with their labels, run orders and the task statements of the next proofs. The written proofs cite it for their task statements; where it and the
  paper differ, the paper and the proofs are current.

## Notes

- The records are the files produced by the runs and are kept unchanged. Some working notes and the validator
  copies under `experiments/ic-curve-validation/code/` contain the absolute path `/home/user/ml-cryptanalysis` of
  the machine they ran on; replace it with the root of this repository to rerun them.
- Some `verification.json` files hash logs that were not kept (`experiments/ic-bounds-d5/logs/`); the results
  they summarize are in the same directory.

## Citing

Please cite the paper; `CITATION.cff` has the entry for this repository.

## Licence

The software (the Python code, the Lean development and the scripts) is under the MIT licence (`LICENSE`). The
papers, figures, written proofs, notes and data records are under the Creative Commons Attribution 4.0
International licence (`LICENSE-CC-BY-4.0.txt`). In detail: MIT covers `csd/`, `scripts/`, `tests/`, `formal/`
and the code inside `experiments/`; CC BY 4.0 covers `paper/`, `docs/`, the Markdown files in `experiments/` and
the data records (the JSON, JSONL, CSV, NPZ and log files in `experiments/`).

Three Lean files contain portions adapted from Apache-2.0 sources and keep those notices in their headers:
`formal/ICFallDegrees/GrowthComplexRatFuncDerivative.lean` (JacobianTwo; licence text in
`formal/third_party/JacobianTwo-APACHE-2.0-LICENSE`), `GrowthPolynomialWronskianBasis.lean` and
`GrowthPolynomialWronskianDistinctDegree.lean` (ArkLib). The adapted portions remain under Apache-2.0.
