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
| Abstract | paste `eprint-abstract.txt` (2901 chars; `$…$` math renders with MathJax) |
| Category | Attacks and cryptanalysis |
| Keywords | last fall degree, index calculus, summation polynomials, elliptic curve discrete logarithm, Koblitz curves, Groebner bases, Lean formalization |
| Publication information | Published nowhere else |
| Content | `../curve-lfd-bounds.pdf` (49 pages) |
| Additional notes (public) | empty for a first version; on a revision, what changed |
| Message to editors (private) | empty unless something needs explaining |
| License | CC BY proposed in the earlier form; author must choose and confirm |
| Confirmation | tick after reading the acceptance, license and withdrawal conditions |

Enter each author with an email where known (at least one is required);
ORCID and ROR identifiers are optional.

## Source snapshot

Regenerated on 16 November 2026 from the combined evidence snapshot
[`14b6145`](https://github.com/isogeny-labs/lfd-summation-systems/commit/14b6145).
It includes the merged manuscript and proof work, the completed split-presentation
records, the capped five-summand run, and the first second-implementation depth-14
lift record from pending PR #24. That research PR is not merged by this paper
revision. The limits of the input-space certificate remain in the manuscript.
The data appendix cites this frozen snapshot. The files below are those of
the commit carrying this sheet.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `6fa55466b930308a257a1bb63c3543e1961a85afddfad49c6624a03263ebd078` |
| `paper/curve-lfd-bounds.tex` | `ede29ef22653a3006db6325e1de27780a9d7a137465a31206c290e47b9db21b4` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `b344d65276ff962ce972cf76a47d1fcdf68e217750c728ca895a1f23fd191de5` |

The PDF was built with pdfTeX 1.40.25 (TeX Live 2023, Debian) by
`python3 scripts/build_ic_curve_paper.py`; two consecutive builds gave the
same bytes. Earlier packages were built with Tectonic, so a Tectonic build of
the same source gives a different file with the same content. The figures
were regenerated in the same run.

Lean: the audited sources are those of commit `1c929ae`, whose complete
Linux build and audit passed (225 modules, 671 default declarations). The
eight `GrowthNonTorsion*.lean` modules added since are outside that audit
and are checked by the workflow `non-torsion-bridge.yml`. Check that both
workflows are green on the commit that is submitted.

## Checks

- PDF builds without LaTeX warnings, unresolved references or overfull boxes.
- A4, author email on page 1, PDF metadata title matches the manuscript.
- Abstract agrees with the LaTeX source.
- Form constraints checked on 15 November 2026 and again on the text of 16 November: the abstract is ASCII apart from the "ö" of
  "Gröbner", with `$…$` math only (no `\emph`, `\cite` or other text commands; the form asks
  for at least 64 characters), and every keyword phrase is at most 40 characters.
- Mathematical and reproducibility limits: see `../README.md` and
  `../../docs/reviews-publication-20261114.md`.
- Publication information and license require the author’s confirmation at submission;
  this sheet does not attest that no other submission or publication exists.

## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
