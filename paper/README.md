# Main manuscript: ePrint preparation

**Status, 19 November 2026:** ready for submission; what remains is the
author's own confirmation of the ePrint form (license, publication
information, authorship) and the submission itself, which only an author can
make. The cleanup of
17 November put the full proof of the growth theorem into Appendix C,
shortened the abstract, added a statement on the use of AI tools (the author
must confirm its wording before submission) and removed wording about
earlier drafts. The revision of 18 November shortened the abstract further,
added a significance paragraph in the introduction, and added the
second-engine corroboration at depth 14 on the unit instance over p = 65537
([update note](latest-update-20261118.md)). The revision of 19 November
(abstract unchanged) adds the remaining corroborations of those runs: a
third, separately written lift of the second engine's depth-14 spaces on
both unit instances (1 in V_7), its depth-14 closure over p = 163841, and degree-6 closures
without 1 at depths 15 and 16 whose profiles follow the recursion predicted
on 10 November; none is a certificate, and LFD <= 7 at depths 15 and 16 is
open ([update note](latest-update-20261119.md)). Also on 19 November the
E6 proper-ideal growth theorem of PR #33 was folded into item (5) of the
target-of-order-two theorem: the lower bound holds at every admissible prime
for k >= 38 without a unit-ideal premise, which also settles the
Fermat-factor caveat of item (4) in that range; the capped local searches of
PR #35 were merged as a record and change no claim. It
combines the rational prime-field growth results, the solving-degree section
and the kernel-checked function field and pole orders for targets that are
not of order two (Proposition "Function field and pole orders"). No ePrint
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
establish the LFD lower bound; for k >= 38 item (5) of that theorem supplies
it from a degree-two witness, and below k = 38 it remains undecided. At depth 14 the lower bound is certified, and the upper bound
uses a degree lift replicated by a second implementation on both unit instances.
The lift still assumes that its input space is contained in the true degree-6
closure; the independent closed-space check proves the reverse inclusion.
There is no independent certificate for that input containment; since 18 November
it also follows if the later version of the second engine, whose degree-6 closure of
the unit instance over p = 65537 has the same dimension (one run, not certified),
adds only admissible products. Since 19 November a third, separately written lift of
that engine's space over p = 65537 also gives 1 in V_7 (one run, not certified). Binary replications
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

- CI on the submitted commit must pass both the default Lean audit and the
  all-module check (`scripts/lean_check_all.py`). The reconciliation commit
  `1f3f55a` (PR #32) passed: 235 modules built, 1,038 declarations across 93
  files audited, and no axioms beyond `propext`, `Classical.choice` and
  `Quot.sound`. The older saved default audit covers 671 declarations and is
  not the full count.
- The data appendix cites the combined evidence commit `2cb270c`. It stays
  reachable only if this branch is merged with a merge commit; after a squash
  or a rebase, update the citation and rebuild.
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
reconciled into this version. At that synchronization two contact addresses
appeared on page 1; since 16 November 2026 the only one, on page 1 and in the
submission sheet, is `adam@isogenylabs.io`. This
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

## Revision of 16 November 2026

The manuscript sources of pull requests 13 (solving degree) and 19 (rational
prime-field growth) were merged; the one conflict, in the list of
contributions, keeps both items. New text: Section "Targets with nonzero
ordinate: the function field and two places", with matching sentences in the
contributions, the status paragraph of the growth theorem, the limitations,
the future directions, the conclusion, the Lean appendix and the data table.
The abstract gained no claim from that section, because it proves no growth
bound. The section on the solving degree and the abstract sentence about it
are those of pull request 13.

The fall-degree measurements of the `crypto` repository (stage B: three to
six summands) are the split presentation of Section "Three to eight
summands", with one free auxiliary per partial sum. On the cells both have
(three summands, depths 2 to 4) they equal the certified values of this
repository on every unit instance (`scripts/check_crossrepo_stage_b.py`).
The paper reports this agreement and distinguishes the completed certified
split-presentation cells from capped or external-only readings
(`../experiments/crossrepo-pkm-stage-b/`). The 2-isogeny arm of that study is
a different presentation and is not in the paper. An earlier version of this
note called the stage B values inconsistent with the bound 7 for three
summands; that compared the split with the direct presentation and was wrong.

Pull requests 1 (review note of 12 November) and 9 (clean build of the
582-declaration snapshot) are merged as archives; the current status notes
were kept where they conflicted.

## Latest paper reconciliation, 16 November 2026

This revision is based on frozen evidence snapshot `2cb270c` from pending
PR #24. The split-presentation table now includes exact values 6 and 7
on the completed cells; the five-summand depth-3 cell remains a lower
bound of 8 after an out-of-memory termination. The second implementation
of the depth-14 lift agrees on both unit instances over p=65537 and p=163841, with the
shared degree-6 producer limitation stated in the text. The reviewed
[update record](latest-update-20261116.md) distinguishes the saved large-run
records from the fresh bounded CPU replay. No cloud run was launched for
this paper update.

The prose pass retains the isogenylabs.io-only contact. The solving-degree
section and the kernel-checked non-two-torsion function field and pole
orders from the combined manuscript are preserved.
