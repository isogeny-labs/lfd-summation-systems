# Main manuscript: ePrint preparation

**Status of the 582-declaration candidate, 14 November 2026:** corrected ePrint candidate with a fresh clean Lean build and full axiom audit passed; ready for author submission review. No ePrint
submission or editorial acceptance is recorded by this review. This is a
technical and artifact review, not independent certification of every proof.

Title: *Last fall degrees of elliptic curve summation systems: machine-checked
bounds and exact values*, Adam Buran.

- [Manuscript PDF](curve-lfd-bounds.pdf), with its [LaTeX source](curve-lfd-bounds.tex).
- [Review, corrections and remaining evidence gaps](../docs/reviews-publication-20261114.md).
- [Submission form sheet](eprint-curve-lfd-bounds/EPRINT.md) and
  [plain-text abstract](eprint-curve-lfd-bounds/eprint-abstract.txt).

## Scope and evidence

The manuscript gives bounds on specified polynomial presentations, exact
finite-instance closure results, and a conjectural depth profile. Its chain,
selector and carry bounds have a saved source-bound Lean audit. The growth
theorem has a written proof; its Lean reduction retains the curve-model
hypothesis. At depth 14 the lower bound is certified, and the upper bound
uses a degree lift without an independent certificate. Binary replications
use the same engine. These qualifications are part of the paper's claims.
The manuscript makes no discrete-logarithm cost claim.

The [fresh clean Linux build](../experiments/ic-lean/clean-20261114/README.md)
compiled all 120 project modules from an empty project build directory and
passed the full 582-declaration axiom audit. All 124 reviewed source digests
and all audited axiom dependencies match the saved audit. Lean 4.19.0 and
the committed Mathlib revision were used; only the Mathlib dependency cache
was reused. The earlier 513-declaration cold build and warm audit remain
archived with their original scope. The retained Linux run supplies this
clean-build evidence. The newer 628-declaration integration in PR #8 is a
separate manuscript version and is not covered by this receipt.

## Build and checks

From the repository root, with Python 3.10+, `requirements.txt`, Tectonic or
pdfLaTeX, and Poppler (`pdfinfo`, `pdftotext`):

```sh
python3 scripts/build_ic_curve_paper.py
python3 scripts/check_publication_artifacts.py
```

The build regenerates the sweep tables and all manuscript figures. It rejects
unresolved references, LaTeX warnings and overfull boxes before replacing the
PDF. The second command checks saved digests, audit-log/report agreement,
sweep identities and counts, the submission abstract, PDF metadata, page
format and first-page author/contact. It performs no Lean or closure replay.
If the abstract changes, update the submission text from the LaTeX abstract
and rerun the check. The pinned Lean reproduction is in
[formal/README.md](../formal/README.md); use a separate clean checkout and
retain its logs rather than overwriting the archived audit.

## Before submission

- Completed: the clean build and 582-declaration audit are retained under
  `experiments/ic-lean/clean-20261114/`; the read-only publication checker
  validates their source binding and per-module success.
- Make the reviewed source, PDF and evidence snapshot publicly addressable;
  this review starts from the update in PR #4, not the older public `main`.
- Read the review's remaining limits, particularly the depth-14 lift and the
  written growth proof. An independent depth-14 certificate would strengthen
  the result; the qualified preprint does not claim to have one.
- The author confirms authorship, publication information and the ePrint
  distribution license. The earlier sheet's CC BY entry is a proposal.

The existing `paper/` layout is retained to preserve this companion repository's
published paths; this change adds no separate research manuscript.
