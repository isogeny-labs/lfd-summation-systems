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
| Abstract | paste `eprint-abstract.txt` (2876 chars; `$…$` math renders with MathJax) |
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

Regenerated on 16 November 2026 on the branch that combines the open
manuscript, proof and evidence pull requests (1 and 9 to 20) on
top of `main` at `e16cfb8`. The combined evidence, before the manuscript
edits of this revision, is commit
[`26be75c`](https://github.com/isogeny-labs/lfd-summation-systems/commit/26be75c),
which the data appendix cites. The files below are those of the commit that
carries this sheet.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `8237902dc9f600874662c7f96f065de1f425fbe630571b64cd46a25dd65a1404` |
| `paper/curve-lfd-bounds.tex` | `409452a05f32fdf04ea3946c2f5698a29f6af4fd73bc6135796fdcde1f6d9af1` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `9033676d9d5156f441e9230aa9d0212879b6efb25501c52717e80172ba4ebe66` |

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
