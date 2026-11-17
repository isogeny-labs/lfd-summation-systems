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
| Abstract | paste `eprint-abstract.txt` (2903 chars; `$…$` math renders with MathJax) |
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
[`2cb270c`](https://github.com/isogeny-labs/lfd-summation-systems/commit/2cb270c).
It includes the merged manuscript and proof work, the completed split-presentation
records, the capped five-summand run, and both second-implementation depth-14
lift records in PR #24. This package is prepared for merging that PR after
its final CI checks pass. The limits of the input-space certificate remain in the manuscript.
The data appendix cites this frozen snapshot. The files below are those of
the commit carrying this sheet.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `2c8357fdb2a000c33e246cdb36f2e6bdd95c8cd4dae4657f62fcfc9d9ce59b82` |
| `paper/curve-lfd-bounds.tex` | `9b2764756520a04d08b1dd07447f7fb454cb54ad5499f2cb0ed24eeba9328d0b` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `c9eda86df3d7de4e064999dd4d078fc5aba3b3a7586f906b85360bc9e9eca819` |

The PDF and figures were regenerated with Tectonic 0.16.9 by
`python3 scripts/build_ic_curve_paper.py`. The PDF producer is xdvipdfmx
(0.1). This replaces the earlier, incorrect pdfTeX build attribution.
The artifact hashes above identify the submitted files; no byte-for-byte
reproducibility claim is made for other TeX engines.

Lean: the saved default audit reports 671 declarations. Final PR CI builds
all modules and checks every printed axiom list with
`python3 scripts/lean_check_all.py`, including the non-two-torsion modules
outside the older default audit. The merge is gated on that final CI run;
see PR #24 for the tested commit and run. This sheet does not treat the
older saved audit as evidence for modules added later.

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
