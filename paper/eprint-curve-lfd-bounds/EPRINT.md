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
| Abstract | paste `eprint-abstract.txt` (1291 chars; `$…$` math renders with MathJax) |
| Category | Attacks and cryptanalysis |
| Keywords | last fall degree, index calculus, summation polynomials, elliptic curve discrete logarithm, Koblitz curves, Groebner bases, Lean formalization |
| Publication information | Published nowhere else |
| Content | `../curve-lfd-bounds.pdf` (59 pages) |
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
separately written degree lift on the second engine's k = 14 space over
p = 65537 (1 in V_7, same dimensions as the recorded lift); the second
engine's k = 14 closure on the unit instance over p = 163841; degree-6
closures without 1 at k = 15 (p = 65537) and k = 16 (new unit instance over
p = 786433), whose profiles follow the recursion predicted on 10 November;
the prediction scoring, limitations and future directions updated
accordingly. The lift of the p = 163841 space had not finished
(`../latest-update-20261119.md`).

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `b44cdd8c56a94bdf1ca167c7eff89459b2760618c4b7ac7bc85eb440caab24d5` |
| `paper/curve-lfd-bounds.tex` | `e95692753c5d3833db83575b31daf08705ea75a142a11cfea5d65550a0567170` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `081f075e151196e44eb685d80e59293385d1f7c06fdcae68850b145d48ca7c43` |

The PDF and figures were regenerated on 19 November 2026 with Tectonic 0.17.0
(conda-forge; bundle `tlextras-2022.0r0`) by
`python3 scripts/build_ic_curve_paper.py`. The PDF producer is xdvipdfmx
(0.1). This replaces the earlier, incorrect pdfTeX build attribution.
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
- Form constraints checked on 15 November 2026 and again on the texts of 16, 18 and 19 November: the abstract is ASCII (the earlier texts had one "ö", in
  "Gröbner"), with `$…$` math only (no `\emph`, `\cite` or other text commands; the form asks
  for at least 64 characters), and every keyword phrase is at most 40 characters.
- Mathematical and reproducibility limits: see `../README.md` and
  `../../docs/reviews-publication-20261114.md`.
- Publication information and license require the author’s confirmation at submission;
  this sheet does not attest that no other submission or publication exists.

## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
