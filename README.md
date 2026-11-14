# Last fall degrees of elliptic curve summation systems: code and data

Companion repository of the paper

> Adam Buran, *Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values*.

It holds the Lean development, the closure engines, the independent checker, the written proofs and every
record behind the numbers in the paper. The paper is `paper/curve-lfd-bounds.pdf`. See [publication status and build instructions](paper/README.md) before using a submission package.

## Layout

| Path | Content |
|---|---|
| `paper/` | LaTeX source, figures and PDF of the paper; `binary-lfd.*` is the companion manuscript it cites for the binary bounds |
| `formal/` | Lean 4 project `ICFallDegrees` (Lean 4.19.0, Mathlib pinned in `lake-manifest.json`) |
| `csd/` | closure engines and certificate code (Python) |
| `scripts/` | drivers: one script per experiment, the figure script and the paper build |
| `docs/proofs/` | written proofs that are not formalized (the depth reduction, genericity, the growth theorem `N8-growth.md`, the block-of-four lemmas `N9-block-four.md`) and the plans of the Lean formalizations |
| `docs/reviews-*.md` | reviews of the written proofs, the red-team reviews of the growth proof and the audit of the paper |
| `experiments/ic-lean/` | audit of the Lean development (declarations and their axioms) |
| `experiments/ic-curve-sweep/` | the sweep of Sections 5 and 6 |
| `experiments/ic-chain-depth8/`, `-depth9/`, `-depth10/` | two-summand systems at depths 8, 9 and 10 |
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

Most experiment directories have a `README.md`, `RESULTS.md`, `VERIFY.md` or `REPORT.md`; others retain
summary and verification records. Appendix C of the manuscript identifies the evidence used.

## Requirements

Python 3.10 or later with `numpy`; `matplotlib` for the figures and `sympy` for a few symbolic checks
(`requirements.txt`). The runs at depths 11 to 14 used a GPU through `cupy`; nothing else needs one. The paper
builds with pdfLaTeX or Tectonic. The Lean project needs Lean 4.19.0 (see `formal/README.md`).

## Quick checks

    python3 scripts/ic_bounds_msum.py selftest      # three to eight summands: must print SELFTEST PASSED
    python3 scripts/plot_ic_lfd_figures.py          # regenerates three figures from the records
    python3 scripts/build_ic_curve_paper.py         # regenerates tables, all manuscript figures and PDF
    python3 scripts/check_publication_artifacts.py  # checks saved hashes, audit, sweep and submission abstract

## Reproducing the main computations

    python scripts/verify_lean_fall_degrees.py --lake /path/to/lake
    python scripts/run_ic_curve_sweep.py --check-only
    python scripts/run_ic_chain_depth8.py --check-only
    python scripts/run_ic_chain_depth8.py --check-only --depths 9 --primes 7681 \
      --out experiments/ic-chain-depth9 --workers 1 --product-cap 12000 \
      --shift-max-cells 450000000
    python scripts/run_ic_chain_depth10_probe.py --tag rerun
    python scripts/run_ic_chain_depth10_degree6.py --combos 3000
    python scripts/run_ic_depth_step.py
    python scripts/analyze_ic_base_classes.py
    python scripts/ic_bounds_msum.py run --max-log2-points 24

The first command builds the Lean project and prints the axiom dependencies of every audited declaration. The
depth-ten runs take about two hours each and the degree-six closure about five hours. The commands for the GPU
runs and for the degree lift are in `experiments/ic-bounds-k14/README.md` and `LIFT.md`.

## What kind of evidence each result has

The paper uses three terms. *Kernel-checked*: the proof is verified by the Lean kernel. *Written*: there is a
proof on paper (`docs/proofs/`) that is not formalized. *Exact*: the outcome of a computation in exact arithmetic
that ends with a certificate. The value 7 at depth 14 rests on the degree lift of
`experiments/ic-bounds-k14/LIFT.md`, which has a written proof and passed a control on a curve with solutions but
has no independent certificate.

## Notes

- The records are the files produced by the runs and are kept unchanged. Some working notes and the validator
  copies under `experiments/ic-curve-validation/code/` contain the absolute path `/home/user/ml-cryptanalysis` of
  the machine they ran on; replace it with the root of this repository to rerun them.
- No licence has been chosen yet.
