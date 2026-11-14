# Main manuscript: ePrint preparation

**Status, 14 November 2026:** corrected candidate for author review. No ePrint
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

The saved audit covers 582 declarations and matches all 124 recorded source
digests. A clean build of the complete current project has not been performed
in this review. The archived cold build covers the earlier 513 declarations.
The Python CI cited in the upstream PR does not build Lean. Obtain a fresh
clean build and axiom audit before treating the current formalization as
independently rebuilt.

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

- Obtain and retain a clean build/audit of the 582-declaration version, or
  explicitly accept the saved warm-build evidence as the publication limit.
- Make the reviewed source, PDF and evidence snapshot publicly addressable;
  this review starts from the update in PR #4, not the older public `main`.
- Read the review's remaining limits, particularly the depth-14 lift and the
  written growth proof. An independent depth-14 certificate would strengthen
  the result; the qualified preprint does not claim to have one.
- The author confirms authorship, publication information and the ePrint
  distribution license. The earlier sheet's CC BY entry is a proposal.

The existing `paper/` layout is retained to preserve this companion repository's
published paths; this change adds no separate research manuscript.
