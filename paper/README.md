# Main manuscript: ePrint preparation

**Status, 15 November 2026:** corrected candidate for author review. No ePrint
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
hypothesis, except for curves with a target of order two, where the
bound is kernel-checked without it (over the complex numbers for a
one-parameter family of curves, and for every allowed rational parameter over admissible prime fields above
an explicit denominator-sensitive cutoff; E6 retains the sharper cutoff B_k=9^(4^k)). For prime factors of Fermat numbers,
the ideal is proper (also kernel-checked), so the unit-ideal argument does not
establish the LFD lower bound; the numerical bound there remains undecided. At depth 14 the lower bound is certified, and the upper bound
uses a degree lift without an independent certificate. Binary replications
use the same engine. These qualifications are part of the paper's claims.
The manuscript makes no discrete-logarithm cost claim.

The rational prime-field extension passed Linux CI at `1c929ae`:
[run 43777495526](https://github.com/isogeny-labs/lfd-summation-systems/actions/runs/43777495526).
The complete build covers 225 modules and audits 883 distinct declarations
across 87 files; the default audit covers 671 declarations. All 275 recorded
source hashes match, with only the three allowed foundational axioms and no
project axioms or placeholders. The paper update preserves those Lean sources.
The new corollary, its explicit denominator guards, the six exceptional-prime
profiles and the finite rational-reduction checks are mapped to declarations
and records in the appendices.

Earlier build history follows. The saved audit covers 671 declarations of `Audit.lean` and matches every
recorded source digest. After the review of 14 November the complete project
was built from an empty build directory (macOS, Lean 4.19.0, pinned Mathlib):
all 217 modules build, and the 840 declarations over all files that print
axioms use only `propext`, `Classical.choice` and `Quot.sound`
(`../experiments/ic-lean/BUILD-COLD-20261114.md`). Two modules were added
afterwards (219 modules; incremental build, addendum of the same file). The workflow
`.github/workflows/ci.yml` repeats the build and audit on Linux for every
pull request. The archived cold build of 12 November covers the earlier 513
declarations.

The complete Linux audit subsequently passed at commit `49585b5`: 219 library
modules and 848 declarations across 82 audit files. The proper-ideal wording
correction changes no Lean declarations or proofs; the corrected sources also
passed the local 636-declaration audit. Original CI logs and the comparison of
all 270 source hashes are retained in
`../experiments/ic-lean/growth-verification-49585b5/`.

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

- Done on 14 November: a clean build and audit of the current Lean project
  (see above). Check that the same run is green in CI on the commit that is
  submitted.
- Make the reviewed source, PDF and evidence snapshot publicly addressable:
  merge the pull requests into `main`. The data appendix of the manuscript
  cites a fixed snapshot commit; it must stay reachable (merge with a merge
  commit, or update the citation and rebuild).
- Read the review's remaining limits, particularly the depth-14 lift and the
  written growth proof. An independent depth-14 certificate would strengthen
  the result; the qualified preprint does not claim to have one.
- The author confirms authorship, publication information and the ePrint
  distribution license. The earlier sheet's CC BY entry is a proposal.

The existing `paper/` layout is retained to preserve this companion repository's
published paths; this change adds no separate research manuscript.

The earlier complex rational-parameter extension is in
`formal/ICFallDegrees/GrowthRationalParameterClassification.lean`: the exact
complex grid-zero classification is a=-3,-1,1, and the growth window now
also covers every rational a outside {-3,-3/4,-1,1}, including a=0,2,-2.
The local default build and 644-declaration audit pass with compiled
dependencies reused; `experiments/rational-parameter-lean-20261114/` records
the commands and scope. Linux CI on the submitted commit is the independent
build check for this extension.

## Revised-paper synchronization

The prose and author contacts from `ml-cryptanalysis` commit `94290676` are
reconciled into this version. Both `adam@isogenylabs.io` and
`a.buran28@gmail.com` appear on page 1 and in the submission sheet. This
version retains the later companion corrections and fixed-curve growth
theorem. The PDF, figures, verification report and ePrint abstract were
regenerated from these sources; the build and all 295 publication checks
pass. `SYNC-94290676.json` records the artifact hashes and source commits.

After that synchronization the growth section of the manuscript was extended
(the kernel-checked theorem now covers a family of curves over the complex
numbers and proves properness of the ideal at Fermat-factor primes; Theorem "Growth for a
target of order two"). The abstract, the author block and the figures did not
change. `SYNC-94290676.json` is the record of the synchronized version and
its hashes of the manuscript source, PDF and build report describe that
version, not the current one; the current hashes are in
`../experiments/ic-curve-paper/verification.json`.

The prime-field extension is in `GrowthRationalPrimeNative.lean` and
`GrowthRationalPrimeCutoff.lean` under `formal/ICFallDegrees/`. For a=m/q
with q>0 it proves the least-LFD window for every admissible prime above
C=2^N H^(N*N)+2N+q+|m+3q|+|4m+3q|, with N=2^k and
H=3q+|m+2q|. The cutoff proves denominator invertibility and nonsingularity;
no lowest-terms assumption is needed. E6 retains the sharper 9^(4^k) cutoff.
