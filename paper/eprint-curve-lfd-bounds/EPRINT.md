# IACR ePrint submission: Last fall degrees of elliptic curve summation systems: machine-checked bounds and exact values

Reviewed against the [current submission form](https://eprint.iacr.org/submit)
on 14 November 2026 and again on 15 November 2026 (same fields, categories and
license options). The abstract, the page count and the digests below are checked
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
| Affiliation / contact (on page 1) | adam@isogenylabs.io |
| Abstract | paste `eprint-abstract.txt` (1622 chars; `$…$` math renders with MathJax) |
| Category | Attacks and cryptanalysis |
| Keywords | last fall degree, index calculus, summation polynomials, elliptic curve discrete logarithm, Koblitz curves, Groebner bases, Lean formalization |
| Publication information | Published nowhere else |
| Content | `../curve-lfd-bounds.pdf` (60 pages) |
| Additional notes (public) | empty for a first version; on a revision, what changed |
| Message to editors (private) | empty unless something needs explaining |
| License | CC BY proposed in the earlier form; author must choose and confirm |
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

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `b0a723ea94ffa97ad76731b391a9d26b73e0e1a393b33e3f05958ad1895e338f` |
| `paper/curve-lfd-bounds.tex` | `7db350ee60f375b01c9a1328b314c5b91c54836350d68dd38811bb809a096a57` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `35d4028ee813c0bad789921ae116e44c79785bcfa1a81aeb28123fd266af70c2` |

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

Lean: the saved default audit reports 671 declarations. Final PR CI builds
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
- Form constraints checked on 15 November 2026 and again on the texts of 16, 18, 19 and 20 November: the abstract is ASCII except the "ö" of
  "Gröbner" in the text of 20 November (later), as in the texts before 16 November, with `$…$` math only (no `\emph`, `\cite` or other text commands; the form asks
  for at least 64 characters), and every keyword phrase is at most 40 characters.
- Mathematical and reproducibility limits: see `../README.md` and
  `../../docs/reviews-publication-20261114.md`.
- Publication information and license require the author’s confirmation at submission;
  this sheet does not attest that no other submission or publication exists.

## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
