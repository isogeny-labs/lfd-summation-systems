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
| Abstract | paste `eprint-abstract.txt` (2978 chars; `$…$` math renders with MathJax) |
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

Regenerated on 18 November 2026. The evidence reviewed for the manuscript is the
frozen snapshot [`2cb270c`](https://github.com/isogeny-labs/lfd-summation-systems/commit/2cb270c)
(merged with PR #24 at `ce4d785`). This revision applies the independent review of
the 16 November sections (corrected solving-degree counts, narrowed claims, Lean
scope; see `../README.md`, "Revision of 18 November 2026") and adds the root-finding
item, whose records are at commit
[`a28475a`](https://github.com/isogeny-labs/lfd-summation-systems/commit/a28475a).
Both commits must stay reachable: merge with a merge commit. The limits of the
input-space certificate remain in the manuscript. The files below are those of
the commit carrying this sheet.

| file | SHA-256 |
|---|---|
| `paper/curve-lfd-bounds.pdf` | `4d4b944a6907c525b55f485bb7dc5c5ddcb8e01782c64e1d0a2a6a8c791fb16a` |
| `paper/curve-lfd-bounds.tex` | `2ff21661a24cf1b34d6815bafc213a5ff77a0cfea0b8b1d8aa95caf8eb1dbf44` |
| `paper/eprint-curve-lfd-bounds/eprint-abstract.txt` | `9ff0cf29587bc1c532eaa39a28547e6aa47f22d7f2f06c5dfb73740bc2625627` |

The PDF and figures were regenerated on 18 November 2026 with pdfTeX 1.40.25 (TeX Live 2023/Debian),
matplotlib 3.11.2 and numpy 2.5.3 by `python3 scripts/build_ic_curve_paper.py`; the PDF producer is
pdfTeX-1.40.25. The previous version was built with Tectonic 0.16.9 (producer xdvipdfmx); the figure
script and its input records are unchanged, so the figures differ only by the plotting and TeX toolchain.
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
- Form constraints checked on 15 November 2026 and again on the texts of 16 and 18 November: the abstract is ASCII apart from the "ö" of
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
