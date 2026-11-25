# IACR ePrint submission: Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values

Reviewed against the [current submission form](https://eprint.iacr.org/submit)
on 14 November 2026, again on 15 November 2026 and again on 22 November 2026
(same fields, categories and license options each time). The abstract, the page count and the digests below are checked
against the manuscript by
`python3 scripts/check_publication_artifacts.py` (from the repository root).
Update this sheet after any change to the title, authors, abstract or PDF. Submit at
https://eprint.iacr.org/submit. **Only an author submits**: the license is
irrevocable, an accepted paper can never be fully removed, and every named
author must agree.

| field | value |
|---|---|
| Title | Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values |
| Authors (in order, UTF-8) | Adam Buran <adam@isogenylabs.io> |
| Affiliation / contact (on page 1) | adam@isogenylabs.io (a contact address; page 1 has no affiliation line, which the archive does not require when a contact address is given) |
| Abstract | paste `eprint-abstract.txt` (1908 chars; `$…$` math renders with MathJax) |
| Category | Attacks and cryptanalysis |
| Keywords | last fall degree, index calculus, summation polynomials, elliptic curve discrete logarithm, Koblitz curves, Groebner bases, Lean formalization |
| Publication information | Published nowhere else (the author confirms; this sheet does not attest it) |
| Content | `../curve-lfd-bounds.pdf` (66 pages) |
| Additional notes (public) | empty for a first version; on a revision, what changed |
| Message to editors (private) | empty unless something needs explaining |
| License | CC BY (form option `CC BY: Creative Commons Attribution`), the licence of `paper/` in this repository (CC BY 4.0, `../../LICENSE-CC-BY-4.0.txt`); the author confirms, and checks it against any later venue's copyright agreement |
| Confirmation | tick after reading the acceptance, license and withdrawal conditions |

Enter each author with an email where known (at least one is required);
ORCID and ROR identifiers are optional.

## Source snapshot

Regenerated on 19 November 2026 (revision below) on top of the combined
evidence snapshot of the manuscript merged in PR #32, which includes the merged manuscript and proof
work, the completed split-presentation records, the capped five-summand run,
both second-implementation depth-14 lift records in PR #24, the full proof of
the growth theorem (Appendix C), and the validation follow-ups of PR #31.
The limits of the input-space certificate remain in the manuscript.
The data appendix cites the frozen snapshot `2cb270c` for the records added
before 16 November; later additions cite their own commits. The files below
are those of the commit carrying this sheet.

Revision of 18 November 2026: shortened abstract (180 words, 1291
characters; previously 432 words), a significance paragraph in the
introduction, the full proof of the growth theorem in Appendix C, a statement
on the use of AI tools, and the second engine's degree-6 closure at k = 14 on
the unit instance over p = 65537 (a corroboration from the `isogeny-labs/crypto`
records at commit `98b7d87`, not a certificate; `../latest-update-20261118.md`).

Revision of 19 November 2026 (abstract unchanged): the remaining results of
the runs `paper-runs-20261118` from the `isogeny-labs/crypto` records at
commit `12308b4`, each a corroboration and not a certificate: a third,
separately written degree lift on the second engine's k = 14 spaces on
both unit instances (1 in V_7, same dimensions as the recorded lifts; the
p = 163841 lift from the records at commit `6fc6ccf`); the second
engine's k = 14 closure on the unit instance over p = 163841; degree-6
closures without 1 at k = 15 (p = 65537) and k = 16 (new unit instance over
p = 786433), whose profiles follow the recursion predicted on 10 November;
the prediction scoring, limitations and future directions updated
accordingly (`../latest-update-20261119.md`).

Revision of 19 November 2026, later the same day (abstract unchanged):
merged the E6 proper-ideal growth theorem (PR #33) into item (5) of
the target-of-order-two growth theorem: for $E_6$ over $\F_p$, every $k\ge38$ and
every prime $p$ with $2^{k+1}\mid p-1$, the lower bound
$\lfd>\lfloor(k-10)/14\rfloor$ holds whether or not the ideal is the
unit ideal, by a degree-two ideal witness rather than the unit-ideal
step. This also decides the Fermat-factor caveat of item (4) for
$k\ge38$. The capped local searches of PR #35 were merged as a record;
they do not change any manuscript claim.

Revision of 19 November 2026, last (abstract unchanged): the records of the
second engine's runs and of the third lift, the lift itself and a standalone
extract of the second engine were moved from `isogeny-labs/crypto` into this
repository (`second-engine/`); the paper now cites those paths. The runs were
made with binaries built from crypto commits `450d5b2` and `43d5c40` (listed
with their hashes in `second-engine/PROVENANCE.md`); the extract reproduces the
recorded closures that fit in 16 GB (`second-engine/gate-20261119/`).

Revision of 20 November 2026 (readability; abstract changed): the abstract
now says what the last fall degree is and states the contrast with the first
fall degree of the two-summand prime systems, $3$ or $4$ at every depth
(kernel-checked, `PrimeM2General.exact_ordinary_ffd`); it makes no claim that
the manuscript does not already prove or certify (1599 characters;
previously 1291). The manuscript gained a framing paragraph and three
conceptual figures in Section 2 (index calculus at a glance; the three
degrees; the degree-fall timeline), roadmaps for Sections 4 to 6 and for the
depth-11-to-16 evidence, and sentence-level splits in the introduction,
limitations, future directions and conclusion; one parenthesis in the future
directions now names both unit instances for the third lift, as Section 5.5
already did. No bound, label, citation or evidence label changed.

Revision of 20 November 2026, later the same day (abstract changed): the
abstract was rewritten for accuracy and readability after an independent
claim-by-claim audit (1278 characters; previously 1599). It now states the
assumption under test (a roughly constant solver degree as the factor base
grows), the kernel-checked bound $\mathrm{LFD}\le k+2$ at every depth
$k\ge2$, the exact values by depth range with the instances that carry them,
the corrected breakpoints of Petit, Kosters and Messeng, that the law
$\lfloor(k+14)/4\rfloor$ was fitted on $k\le13$, and the hypothesis
$p>2^{k+1}$ of the growth theorem; every clause was checked against the
manuscript body and the Lean sources. The text contains one non-ASCII
character, the "ö" of "Gröbner" (UTF-8). The manuscript body is unchanged;
`../proof-status.json` now lists the kernel-checked bounds the paper claims.

Revision of 20 November 2026, second cleanup pass (abstract unchanged): a
sharper prose pass removed remaining mechanical phrasing and tightened the
reporting text throughout, chiefly the depth-11-to-16 corroboration, the
multi-summand and depth-six accounts, the conclusion and the computational
appendices, and moved the not-certified larger split-cell readings into a
compact table. No theorem, lemma, definition, proposition, conjecture body,
proof, table value, figure, label, citation or evidence label changed, and
every numeric value of the previous version is preserved somewhere in the
text. The manuscript is now 59 pages (previously 60). The abstract and the
"Use of AI tools" paragraph are byte-identical to the previous version.

Revision of 20 November 2026, discrete-logarithm consequence (abstract
changed): the blanket sentence "We make no claim about the cost of discrete
logarithms" is replaced by the consequence the theorems already carry. In
the two-summand relation search of Petit, Kosters and Messeng, an attempt on
a target of order greater than two finds no relation exactly when its ideal
is the unit ideal (a short Frobenius argument, Section "What the results
mean"); there the solving degree equals the LFD, which the growth theorem
bounds by at least floor((k-10)/14)+1 for t != 0 and p > 2^(k+1). Each such
attempt therefore costs a Macaulay-matrix solver at least
binom(2k, floor((k-10)/14)+1) = 2^((0.44+o(1))k) columns, and estimates that
assume a bounded solver degree understate the relation-search cost by a
factor exponential in k. The abstract, introduction, contributions, the
section "What the results mean, and what they do not" and the conclusion
state this; the paragraph "What does not follow" now limits the claim to
the relation search (no total discrete-logarithm cost, no comparison with
Pollard rho). No theorem, lemma, definition, proposition, conjecture body,
proof, table value, figure, label or citation changed. The manuscript is now
60 pages (previously 59); the abstract is 1622 characters.

Revision of 20 November 2026, more summands and a cost model (abstract
changed): Section 5.8 now states the carry bound for every number of
summands (Theorem "Carry bound for every number of summands"): for the
direct m-summand presentation, m >= 3, LFD >= m min(k, m-1) + 1; the exact
value mk + 1 for k <= m-1 (kernel-checked for an arbitrary generator,
`CarryAllM.lfd_eq_full`; the grid step is a written proof); the onset
m(m-1) + 1 for k >= m when c_m(a,b,t) != 0 (kernel-checked for m = 3,
`CarryAllM.m3_onset`; for m >= 4 a written proof with a kernel-checked core,
the step outside Lean being the degree bound of S_{m+1} in each variable);
and genericity (at most m excluded targets per curve when char K = 0 or
char K does not divide kappa_m), with a proof sketch citing the validated
note `../notes/carry-allm.tex`. The second engine's depth-growth closures
(`experiments/msum-growth/`) are added as corroborations, not certificates:
three summands 7 at k = 2..5, 8 at k = 6, >= 8 at k = 7; four summands 13 at
k = 4, >= 13 at k = 5. Section "What the results mean" adds the m-summand
consequence of the proved floor and, as a clearly labelled model-based
conclusion, the fully charged cost model (`experiments/ic-endtoend/`,
`experiments/ic-costmodel/`): with the solver charged at the proved degree
floor, index calculus with m <= 8 does not beat Pollard rho between 2^40 and
2^256 (2.6 to 5.0 bits above rho at prices 0.1 and 1 GAE), with its
assumptions and the excluded literal-top-layer convention. "What does not
follow", the limitations, future directions, conclusion, contributions, the
Lean table and the data tables were updated accordingly. The abstract gains
one sentence (1796 characters; previously 1622); the manuscript is 63 pages
(previously 60).

**Revision of 23 November 2026** (PR #49; abstract and `eprint-abstract.txt`
unchanged). New Section 6.3 with Corollary "Growth for every number of
summands" (written proof; reduction kernel-checked with the splitting
property as a hypothesis, `formal/ICFallDegrees/GrowthSlice.lean`, built
locally and by the all-module CI step) and Remark 28 (the PR #45 route
through a point of order two); matching updates in the contributions, the
ECDLP section, limitations, future directions, conclusion and the Lean
table. The manuscript is 66 pages (previously 63). The PDF was rebuilt by
`python3 scripts/build_ic_curve_paper.py` with Tectonic 0.15.0 as the only
TeX engine on the `PATH`; before the edit, the same build reproduced the
committed PDF of 20 November byte for byte. Figures are byte-identical. The
digests below are those of this revision.

**Revision of 24 November 2026** (PR #58): the abstract gains the
owner-approved clause "and a two-summand slice extends the unbounded growth
to unit-ideal instances with any fixed number of summands" (Corollary
`cor:growthm`, mapped in `../proof-status.json`, `abstract_claims`);
`eprint-abstract.txt` is resynced (1908 characters; previously 1796). No
other source change; 66 pages; rebuilt with Tectonic 0.15.0 by
`python3 scripts/build_ic_curve_paper.py`. The digests below are for the 24 November certified-grid revision in PR #59.
The abstract remains unchanged; the manuscript adds the finite FFD/LFD grid,
the two exact split (5,3) values, and explicit depth-14/depth-18 evidence limits.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `59ffec43bbd0fd56038a18db2938955f2a5014e9fcdc255ac390fcaaaa67e048` |
| `paper/curve-lfd-bounds.tex` | `2e5096c14d63368559eb467a5fc037744b1da3ddafb92d9f26d9eadeedf60d69` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `9264a95a970d8cbc4ae87c46dadc09876db04c9affac6a52706bab7677c6d3da` |

The PDF was regenerated on 20 November 2026 with Tectonic 0.17.0
(conda-forge; bundle `tlextras-2022.0r0`) by
`python3 scripts/build_ic_curve_paper.py`, run with only Tectonic on the
`PATH` (the script prefers pdfLaTeX when it is installed); the PDF producer is
xdvipdfmx and the figures are byte-identical to the previous build. The
earlier build of 20 November used pdfTeX 1.40.25; the versions of 19 November
also used Tectonic 0.17.0.
After the citation check of 20 November (three citation locators corrected; see
`../README.md`, "Citation check"), the PDF was rebuilt with Tectonic 0.17.0 built
from crates.io (`cargo install tectonic --version 0.17.0 --locked`), again with
only Tectonic on the `PATH`; the producer is xdvipdfmx and the figures are
byte-identical to the previous build.
The artifact hashes above identify the submitted files; no byte-for-byte
reproducibility claim is made for other TeX engines.

Lean: the saved default audit reports 679 declarations (taken on 20 November 2026 from the CI artifact `lean-audit` of run 44815755702 on commit 13415c6, artifact digest sha256:269350ad…4c58; it adds the eight `CarryAllM` declarations). Final PR CI builds
all modules and checks every printed axiom list with
`python3 scripts/lean_check_all.py`, including the non-two-torsion modules
outside the older default audit, and the modules added by PR #31. The merge
is gated on that final CI run; see PR #32 for the tested commit and run.
This sheet does not treat the older saved audit as evidence for modules
added later.

## Checks

- PDF builds without LaTeX warnings, unresolved references or overfull boxes.
- A4, author email on page 1, PDF metadata title matches the manuscript.
- Abstract agrees with the LaTeX source.
- Form constraints checked on 15 November 2026 and again on the texts of 16, 18, 19 and 20 November (including the revision adding more summands): the abstract is ASCII except the "ö" of
  "Gröbner" in the text of 20 November (later), as in the texts before 16 November, with `$…$` math only (no `\emph`, `\cite` or other text commands; the form asks
  for at least 64 characters), and every keyword phrase is at most 40 characters.
- Mathematical and reproducibility limits: see `../README.md` and
  `../../docs/reviews-publication-20261114.md`.
- Publication information and license require the author’s confirmation at submission;
  this sheet does not attest that no other submission or publication exists.

## Pre-submission check, 22 November 2026

No manuscript source, abstract or PDF changed in this check; the digests
above still match the files.

- **Form.** The submission form and the operations page were re-read on 22
  November 2026: the same fields, the same seven categories and the same six
  licences; the PDF must fit A4 or US Letter and include the author emails;
  the first page must give the title, author names and a contact address or
  affiliation.
- **Generator.** `.agents/skills/iacr-eprint/scripts/eprint_form.py`
  (workspace trunk) run on this manuscript and the committed PDF, with the
  category, keywords and licence above: every hard rule passes. Its abstract
  output differs from `eprint-abstract.txt` only in leaving `Gr\"obner`
  unconverted; paste `eprint-abstract.txt`, which has the UTF-8 "ö".
- **Build.** A separate build from the same sources with pdfLaTeX (TeX Live
  2023, pdfTeX 1.40.25, Ubuntu 24.04) passed the build gate of
  `scripts/build_ic_curve_paper.py`: no LaTeX warnings, no unresolved
  references or citations, no overfull boxes, 63 A4 pages, figures
  byte-identical. Its PDF differs in bytes and in some line breaks from the
  committed Tectonic PDF, as expected across engines; the committed PDF
  remains the submission file. The text layer of the committed PDF extracts
  the floor brackets and "≠" as `b c` and `6=`; the pdfLaTeX build extracts
  them correctly. If that matters for indexing, rebuild with pdfLaTeX and
  update the digests and this sheet before submitting.
- **Artifact checks.** `python3 scripts/check_publication_artifacts.py`:
  452 checks passed (679 audited declarations, 3005 sweep systems).
- **Page 1.** Inspected: title, "Adam Buran", `adam@isogenylabs.io`, the
  date and the abstract. PDF metadata: title as above, author "Adam Buran".
- **Use of AI tools.** The paragraph is present (page 3, end of Section 1);
  it names the assistants and says that the proof reviews of Appendix B.4
  ("Independent checks") were AI sessions, not human experts. The author
  must confirm its wording (open owner question since PR #27).
- **Claims against `../proof-status.json`.** Every clause of the abstract is
  mapped there (`curve_bounds_manuscript.abstract_claims`) to its evidence
  label: one kernel-checked bound (LFD ≤ k+2), exact computations (the values
  4, 5, 6 and LFD ≥ 7 at k = 14), a lifted upper bound at k = 14, an
  empirical law, written proofs (growth, solving degree, relation-search
  cost) and the m-summand bound (written proof, kernel-checked core). The
  abstract names each label itself, except that the relation-search
  sentence does not repeat that it rests on the written growth proof of the
  sentence before it. No independent human review exists
  (`independent_review: not obtained`).
- **CI.** The `ci` workflow (Python self-tests, Lean build and axiom audit
  with `scripts/lean_check_all.py`) passed on `main` at `7a54f8a`, run
  44864341390. CI on the merge of this sheet's pull request is the check
  for the submitted commit.
- **Work not in this manuscript.** `isogeny-labs/ml-cryptanalysis` PR #356
  (merged 22 November 2026) proves sharper growth constants in Lean there,
  not here: LFD ≥ ⌊(k+1)/14⌋+1 for k ≥ 27 and ⌊(k−5)/13⌋+1 for k ≥ 31
  under the same curve-model hypothesis, and sharper E6 windows (six modules
  `Growth*Sharp`, `GrowthSlope`; written parts still unreviewed). This
  manuscript's weaker statements, ⌊(k−10)/14⌋+1 for k ≥ 38, remain correct.
  The author decides whether to port that work before submitting or in a
  revision.
- **Bundle.** This repository has no `make eprint` target; the form takes
  only the PDF, and the sources are this repository at the submitted commit.

## Author's steps (only the author)

1. Read the form sheet above and confirm authorship, publication
   information, licence and the "Use of AI tools" wording; decide whether
   the sharper constants of ml-cryptanalysis PR #356 go into this version.
2. Optionally send `REVIEW-REQUEST.md` to an outside expert first.
3. Submit at https://eprint.iacr.org/submit with the fields above, the PDF
   `../curve-lfd-bounds.pdf` and the abstract `eprint-abstract.txt`.
4. After acceptance, keep the password, and record the ePrint number in
   `CITATION.cff` (`preferred-citation`), in `../README.md` and on the
   website.

## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
