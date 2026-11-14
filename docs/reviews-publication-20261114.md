# Publication review for IACR ePrint, 14 November 2026

> **Note added when this review was combined with two other changes (14 November 2026).** The review below
> is unchanged and describes commit `548ec4e`. Since then:
> - item 1 of "Remaining decisions" is done: the Lean project was built from an empty build directory and
>   audited (`experiments/ic-lean/BUILD-COLD-20261114.md`; 628 declarations of `Audit.lean`, 840 over all
>   files), and CI now repeats that build;
> - item 3 changed in part: for one curve with a target of order two the growth bound is kernel-checked
>   without the curve-model hypothesis (Theorem "Growth on a fixed curve" of the manuscript, added after
>   this review and not covered by it); the general theorem is as described below;
> - item 5: the repository licence is set (`LICENSE`, `LICENSE-CC-BY-4.0.txt`);
> - the manuscript has 39 pages, and the records in `experiments/publication-review-20261114/` are those of
>   the reviewed 38-page version.

## Scope and verdict

Reviewed the manuscript update at commit `548ec4e` (PR #4), rather than the
older public `main` at `368f129`. Corrections are on
`codex/publication-review-20261114`. The original checkout and all raw
experiment records were left intact. The companion binary manuscript was
consulted for scope, but is not cleared for submission by this review.

**Verdict: corrected ePrint candidate prepared for the author.** The review
found and repaired a contradictory theorem hypothesis and an unsupported
universal inference from finite data. It also corrected a formal-build
coverage claim and submission metadata. The revised manuscript preserves
the depth-14 certificate qualification and the written/conditional status
of the growth theorem. This pass supplies no blanket mathematical
certification, no fresh Lean build, and no new large closure computation.
A fresh clean Lean build remains the most useful validation before submission.

## Literature reviewed before corrections

Started with the repository's previous paper audits and proof reviews:
`docs/reviews-paper-audit-20261112.md`,
`docs/reviews-papers-rereview-20261112.md`,
`docs/reviews-curve-lfd-novelty-20261105.md`, and the second growth review
`docs/reviews-N8-growth-redteam2-20261112.md`. These are prior review evidence,
not fresh independent endorsements. The current Lean statements and source
hashes were checked directly where used below.

Primary sources checked in this pass:

- Huang–Kosters–Yang–Yeo, [On the last fall degree of zero-dimensional Weil
  descent systems](https://arxiv.org/html/1505.02532v2), Definitions 2.1 and
  2.5, Proposition 2.6 and Proposition 2.8 with its proof. The multiplier
  rule applies to derived rows. The solution algorithm requires a radical
  ideal and uses both LFD and solution count. A closure fixed point alone
  does not establish an upper bound for the full ideal.
- Caminata–Gorla, [Solving degree, last fall degree, and related
  invariants](https://arxiv.org/html/2112.05579), Section 2, Theorem 2.9 and
  Example 4.3. These support the presentation/degree distinctions and the
  chain-subsystem comparison. Gaps between successive observed falls cannot
  by themselves certify that all later falls are absent.
- Petit–Kosters–Messeng, [Algebraic approaches for the elliptic curve
  discrete logarithm problem over prime fields](https://people.maths.ox.ac.uk/petit/files/16PKC_primeECDLP.pdf),
  the chain presentation and Sections 4.1–4.2, Tables 4–5. Read the actual
  numerical tables and correction rule: the corrected Magma statistic has
  the reported breakpoints, but is a different quantity from a certified
  LFD. Their unreduced generator has degree 4; the manuscript explicitly
  explains transport to its normalized presentation at values at least 4.
- Kosters–Yeo, [Notes on summation polynomials](https://arxiv.org/html/1503.08001v3),
  the introduction and Assumption 3.3/Theorem 3.5. The NP-completeness claim
  for nonsingular curves is conditional; the manuscript retains that
  qualifier and distinguishes the binary first-fall observation.
- Brownawell–Masser, [publisher record and accessible extract](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/vanishing-sums-in-function-fields/34E7F8E89B830D9CBE606324B7E254D9).
  The extract confirms the characteristic-zero setting, but the complete
  original theorem was not retrieved. Voloch's original theorem was not
  retrieved either. The manuscript discloses this gap and relies on the
  self-contained written argument in `docs/proofs/N8-growth.md`; this pass
  does not independently validate every step of that argument or its
  historical attribution. No new novelty priority is asserted here.

The [ePrint submission form](https://eprint.iacr.org/submit) and
[operations/acceptance conditions](https://eprint.iacr.org/operations.html)
were also checked. The archive requires readable, self-contained work with
proofs or convincing arguments, first-page author/contact details, a PDF
on A4 or US letter, and author agreement to a distribution license. Archive
acceptance is not peer review or verification of the claims.

## Findings and repairs

### P1: contradictory carry-theorem hypothesis — repaired

Theorem `thm:carry` required `h - λg ∈ (G)` with **λ a scalar**, while also
requiring reduced `h` and `g` and `h ≠ cg` for every scalar `c`. Reduced-ideal
triviality makes the first condition imply `h = λg`, contradicting the last.
The printed theorem therefore could not justify its intended applications.

`formal/ICFallDegrees/CarryLower.lean`, theorem
`not_hasLFDBound_of_carry`, instead takes `g h a : P r K`: its multiplier
is a **polynomial**. The chain specialization uses a variable multiplier.
The manuscript now says `h - A g ∈ (G)` for a polynomial `A`, and the proof
first establishes `h ∈ (G ∪ {g})`. This matches the existing Lean statement
and the applications `h = NF(x_1 g)`. No Lean proof or numerical result changed.

### P2: finite law refutations promoted to a universal nonexistence claim — repaired

Section 5.8 said that there is no bound in solution count alone. Refuting
the particular preregistered laws does not exclude all functions of solution
count. The text now says exactly which inference the data support: the
tested solution-count law is false, and other bounds are not ruled out.
The later paragraph and the binary figure caption likewise describe
variation on the computed instances rather than a necessity theorem about
all possible tight bounds.

### P3: Lean CI coverage incorrectly described — repaired

`experiments/ic-lean/README.md` said that the 69 added declarations were
covered by the warm audit **and CI**. PR #4 explicitly reports that the
upstream CI runs Python, not Lean. The README now attributes the evidence to
the saved source-bound warm audit alone. The 513-declaration cold build
is kept with its actual scope.

All 124 source digests in `experiments/ic-lean/verification.json` match
this checkout. The saved axiom log, the report and the `Audit.lean` declaration
list agree on 582 declarations using only `propext`, `Classical.choice` and
`Quot.sound`. These consistency checks do not replay those proofs. The local
Lean 4.19.0 toolchain exists, but this isolated checkout has no resolved
Mathlib dependencies or project build cache. No clean rebuild was completed
in this pass; no build success is inferred from source-token or hash checks.

### P4: publication package and build gate — repaired

The PDF metadata still named the older title. It now matches the printed
title; the manuscript date is 14 November 2026. The ePrint sheet's PDF path
pointed to an absent file in its own directory, and its generator path did
not exist in this standalone repository. Both instructions were corrected.
The author must confirm the publication-information and license fields;
the previous CC BY entry is retained as a proposal, not consent.

The build previously returned success even with warnings or overfull boxes.
It now rejects either before replacing the PDF, in addition to the existing
unresolved-reference gate. Added `scripts/check_publication_artifacts.py`
for read-only artifact consistency checks and `paper/README.md` for status,
reproduction and submission requirements. The main manuscript now links to
the fixed public evidence snapshot `548ec4e`; its evidence no longer depends
on the older contents of the default branch.

## Validation and its limits

- Rebuilt the complete 38-page manuscript, generated tables and all five
  external manuscript figures. Zero LaTeX warnings, unresolved references
  or overfull boxes. FontTools emitted timestamp diagnostics while loading
  existing font files; these are not LaTeX layout errors.
- Rendered all 38 pages and inspected contact sheets for pagination, clipping,
  tables and diagrams. Contact-sheet inspection is a layout check, not a
  word-by-word proof review. Inspected the changed theorem and binary
  paragraphs at reading resolution after the final build.
- Read-only publication check: saved source/output hashes, 582-declaration
  audit agreement, identities/counts/histogram for all 3005 sweep systems,
  abstract/source agreement, resolved citation keys, PDF title, A4 size,
  and first-page author/contact. A deliberate temporary abstract mismatch
  is rejected; the original submission text is restored afterwards.
- `python3 scripts/ic_bounds_msum.py selftest`: passed, including scalar
  reference checks and comparison of 344 saved three-summand records.
- Re-scored all 57 saved deep binary systems against the frozen predictions:
  zero identity/certification-field problems. The failed predictions remain
  failed. This reads records; it does not reproduce the closures.
- Compared saved binary replications: r2 agrees on 57/57 systems and 307/307
  degree records/basis hashes (186/186 proper-subspace hashes); r3 agrees
  on 46/46 systems and 232/232 degree records/basis hashes (133/133 proper
  subspaces). These are same-engine replications, not independent engines.
- Raw experiment records were not rewritten. New check outputs are under
  `experiments/publication-review-20261114/`.

The diagram sources checked were the inline pipeline diagram in the main
LaTeX source, `scripts/build_ic_curve_paper.py` and
`scripts/make_paper_figures.py`. No numerical or structural finding changes
the graph data; figures were regenerated with the paper. The binary figure's
interpretive caption was narrowed with the surrounding text. This review
introduces no new mechanism, experiment or performance claim.

## Remaining decisions and evidence gaps

1. **Recommended before submission:** a clean Lean build and full axiom audit
   of the 582-declaration snapshot on a fresh supported environment, retaining
   logs and source hashes. Otherwise the author must accept the saved warm
   audit as the disclosed verification boundary. No extra cold-build evidence
   is claimed by this review.
2. The depth-14 upper bound still lacks an independently replayable derivation
   certificate, and the large input bases were not archived. Existing checker
   and lift records cannot be bound file-for-file solely by identical
   manifest hashes. The qualified preprint states these limits; an independent
   certificate is needed to strengthen the label, not to preserve the current
   qualified wording.
3. The growth result's analytic curve-model input remains written, not proved
   in Lean. Its primary-source attribution was not checked against the full
   Brownawell–Masser/Voloch originals. A specialist review of the complete
   written argument would strengthen the preprint; the conditional Lean
   statement must remain conditional.
4. The 149,122-system selector test retains summaries rather than individual
   rows and has no demonstrated negative control at equal depths. It must
   remain secondary evidence to the formal proof. No test rows or frozen
   experiment protocols were reconstructed after the fact in this review.
5. Commit/publish the reviewed changes before distributing the revised PDF,
   retaining the fixed evidence snapshot. The author then confirms the
   selected ePrint license and publication history and performs submission.
   The repository software/data license is a separate choice, still unset.

A stronger assertion that every theorem has been independently rebuilt or
that every computed value is independently certified is not supported.
