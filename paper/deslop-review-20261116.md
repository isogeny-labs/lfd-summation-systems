# Deslop prose review, 16 November 2026

Applied Stephen Turner's [deslop skill](https://github.com/stephenturner/skill-deslop/blob/a906154bef375d9d49ed2ad7da13b2db16f0d3d2/SKILL.md),
including its phrase, structure, trope and example references, to the main
manuscript `curve-lfd-bounds.tex`.

## Changes

Made 38 targeted prose revisions across the abstract, introduction, related
work, computation discussion, growth scope, conclusion and methods. The pass
cuts repetitive setup, replaces passive or indirect wording where the actor
is known, shortens structural announcements, and makes the prediction status
explicit. Technical uses of words such as "strictly" and "independently"
remain where they specify a bound or the provenance of a check.

Examples:

- "The theorem is weaker ... in its range and stronger in its status"
  becomes a direct account of the order-two target restriction and the
  remaining genus-one proof obligation.
- "An independent implementation ... with a deterministic full closure
  sweep" becomes a complete sentence describing the implementation and
  how it computes the closure.
- The conjecture discussion now states that computed values agree wherever
  tested and that one prediction remains untested.

The main manuscript and ePrint submission contact now use
`adam@isogenylabs.io` only. The artifact checker enforces the sole first-page
email and the matching submission contact. The paper date is 16 November 2026.

## Preservation and validation

Compared the edited source with the pre-edit manuscript:

- All 32 theorem, lemma, corollary, proposition, definition and proof blocks
  are byte-identical.
- All 2,513 mathematical expressions are unchanged, ignoring whitespace.
- All 275 source digests from the passing Linux Lean audit at `1c929ae`
  still match. This editorial pass changes no Lean source.
- The regenerated PDF has 42 pages, no unresolved references, no LaTeX
  warnings and no overfull boxes. All 303 publication checks pass.
- The submission abstract matches the manuscript and has 2,418 characters,
  excluding the final newline.
- The title/abstract, genericity discussion and conclusion pages were
  rendered and inspected. The first page lists the isogenylabs.io contact.

The depth-14 lift remains uncertified by a final closure certificate. The
general genus-one theorem retains its stated formalization limits. The
finite exceptional-prime computations retain their scope and engine limits.

## Editorial rubric

These scores are a subjective assessment of the prose, not an authorship
test or a mathematical review.

| Dimension | Before | After |
|---|---:|---:|
| Directness | 7 | 8 |
| Rhythm | 6 | 8 |
| Trust | 9 | 9 |
| Authenticity | 6 | 8 |
| Density | 6 | 8 |
| Total | 34/50 | 41/50 |
