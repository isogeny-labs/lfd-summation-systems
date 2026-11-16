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
| Affiliation / contact (on page 1) | adam@isogenylabs.io; a.buran28@gmail.com |
| Abstract | paste `eprint-abstract.txt` (2446 chars; `$…$` math renders with MathJax) |
| Category | Attacks and cryptanalysis |
| Keywords | last fall degree, index calculus, summation polynomials, elliptic curve discrete logarithm, Koblitz curves, Groebner bases, Lean formalization |
| Publication information | Published nowhere else |
| Content | `../curve-lfd-bounds.pdf` (42 pages) |
| Additional notes (public) | empty for a first version; on a revision, what changed |
| Message to editors (private) | empty unless something needs explaining |
| License | CC BY proposed in the earlier form; author must choose and confirm |
| Confirmation | tick after reading the acceptance, license and withdrawal conditions |

Enter each author with an email where known (at least one is required);
ORCID and ROR identifiers are optional.

## Source snapshot

Regenerated on 15 November 2026 from `main` at commit
[`e16cfb8`](https://github.com/isogeny-labs/lfd-summation-systems/commit/e16cfb896cb05607a7f78c43f9db08b118c0b1dc)
(the merge of pull request 8). `python3 scripts/build_ic_curve_paper.py` on a
checkout of that commit reproduces the committed PDF byte for byte, and the
abstract extracted from the LaTeX source equals `eprint-abstract.txt`. The
commit that adds this section changes this sheet, the checker and two notes
only; the files below are those of `e16cfb8`.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `377ef3ddd21ecf5d75538db5642fcff1c3c33592c43bffd744422aa7a0c2c483` |
| `paper/curve-lfd-bounds.tex` | `e83bebb7020a0e4c19b74ac79ddf661f5f85b17f5e195fd2345c971e7db237ff` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `3e73357c02dc04cb12291c20f85765b40ee3e20b5462a70b499b585f93902011` |

The [CI run on `e16cfb8`](https://github.com/isogeny-labs/lfd-summation-systems/actions/runs/43684708863)
passed: the Lean build of all modules, the axiom audit, the Python checks and
the publication checks as they stood at that commit (the checks of this
section were added afterwards). The data appendix of the manuscript cites the
snapshot commit `5a98043`, which is an ancestor of `main`.

## Checks

- PDF builds without LaTeX warnings, unresolved references or overfull boxes.
- A4, author email on page 1, PDF metadata title matches the manuscript.
- Abstract agrees with the LaTeX source.
- Form constraints checked on 15 November 2026: the abstract is plain ASCII with
  `$…$` math only (no `\emph`, `\cite` or other text commands; the form asks
  for at least 64 characters), and every keyword phrase is at most 40 characters.
- Mathematical and reproducibility limits: see `../README.md` and
  `../../docs/reviews-publication-20261114.md`.
- Publication information and license require the author’s confirmation at submission;
  this sheet does not attest that no other submission or publication exists.

## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
