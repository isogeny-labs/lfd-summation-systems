# Re-review of the papers, 2026-11-12

Reviewed at commit `7c66226`. Scope: `paper/curve-lfd-bounds.tex`, `paper/curve-sweep-data.tex` and
`paper/binary-lfd.tex`, read in full, against the first review of the same day
(`docs/reviews-papers-proofs-20261112.md`, at `705279f`).

## Verdict

**None of the 27 findings of the first review has been addressed.** Between `705279f` and `7c66226` the
only changes are the review files. The `.tex` sources and the committed PDFs are unchanged, and the
PDFs match the sources: the phrases "ordinary ordinary", "was not among them", "eight with 2^20 or
2^21" and "424 declarations" are all still in the PDF text. Every finding in Section 2 of the first
review was rechecked against the current text and still applies, at the same line numbers.

This read found no mathematical error in a theorem. It did find four new issues (Section 2): two
inferences in the text that do not follow as written, though the records support their conclusions;
a sharper form of first-review finding 1; and a framing problem in the "bound is attained" result.

Limits: nothing was compiled (there is no `pdflatex` or Lean toolchain here), and no closure was rerun.

## 1. Rechecked independently in this pass

- **Lemma `lem:unit` at `p = 163841`:** `73730` is a non-residue, and `73730^3 + 3·73730 = 28291^2`.
  Because `-3` is also a non-residue, `x^2 + 3` has no rational root. So the summands cannot be
  2-torsion, and the proof carries over as the paper says, even though `α((0,0)) = 3` is a non-square.
- **Proof of `lem:unit`, cases (a)–(c):** correct, including summands at `(±i, 0)` when
  `p = 65537`, `k ≥ 2`.
- **Numbers:** `|M_≤6|`, `|M_≤7|` at `k = 14`; the lift identity
  `C(28,7) = 2·(586,501 + 5519)`; fall totals `494 + 12,260 + 12,382 = 25,136` and
  `… + 12,366 = 25,120` (equal to the free-column counts); `T` dimensions; the `k = 11, 12, 13`
  dimension gaps (`31,347`, `47,263`, `67,298`); the codimensions `10, 189, 3131, 50,256`; the closed
  form for `d_5` at `n = 0..6`; `4k^2+172k−571` against the `D = 4` recursion step `8k+176`;
  `389 ≥ 386` and `605 < 794`; `|M_≤7|`, `|M_≤8|` for 18 variables; `|M_≤7| = 2^14` for 15
  variables.
- **Sweep tables:** class totals give 1260 prime + 1745 binary/Koblitz = 3005. The 455 prime systems
  at `k = 2, 3` and the 985 binary/Koblitz systems with at most 8 Boolean variables (`k ≤ 4`) match
  the validator counts. 12 + 9 systems have value 4 at `k = 6, 7`. Table `tab:msum` gives
  409 of 412.
- **Companion paper:** `thm:uniform` (section and induction), `lem:moore`, the two displayed identities
  of `lem:rewrite`, `lem:rigidity`, `lem:positive` (the carry-run weight `z + 1`), `prop:lower`,
  `cor:scaled` and `cor:envelope` were reread with no gap found. The 702 cases are
  `20 ambient degrees × 4 spaces × 3 curves × 3 target classes − 18`; the 20 arises because
  `k + 3 = 2k` at `k = 3`.

## 2. New findings

N1. **The evaluation-rank inferences are circular as written** (`curve-lfd-bounds.tex`,
    "Two more primes", `k = 13`, `p = 114689`; and the `k = 14` control lift). The text says:
    "Since `V_6` lies in the ideal and the evaluation on `Z` has rank at most 594, the evaluation of
    `M_<6` on `Z` has rank 594". It also says, from `dim(V_7 ∩ M_≤6) = 499,178 − 2370`: "So the
    evaluation of `M_<7` on `Z` has rank 2370". `V ⊆ I` gives only `rank ≤ |M| − dim V = |Z|`, which
    is the trivial bound. A lower bound on the rank, needed for the local-to-global criterion,
    does not follow from the dimensions. Both values were computed directly (`eval_rank_below_D` is
    594 in `experiments/ic-bounds-k14/runs/proper-p114689-k13-D6.json` and 2370 in
    `relaunch/fetched/proper-p114689-k14-D6.json`), so the results stand. The sentences should cite
    that computation instead of deducing the rank. Note that `runs/proper-p114689-k14-D6.json` does not
    carry the field; only the `relaunch/fetched` copy does.

N2. **Finding 1 of the first review also covers Koblitz curves.** The inverse family of `thm:inverse`
    and the witness of `cor:envelope` use `a = 0, b = 1`, which is a Koblitz curve. So the claim
    "Koblitz values are at most 5 for dimensions up to 10" (abstract, Section 7, conclusion,
    `tab:twosided`) is contradicted on a Koblitz curve too: the full field `F_{2^10}`, or the subfield
    `F_{2^10} ⊂ F_{2^20}`, gives 6. The sweep is consistent with the theorem where they overlap:
    `n = 12, k = 6` gives value `max(3, 4) = 4`, which is the Koblitz symmetric maximum at `k = 6`.

N3. **"The bound `mk+1` is attained" is the trivial `δ + 1` fall** (Section 7, "More summands";
    abstract; Figure `fig:variables`). At budget `deg g` the generator cannot be multiplied, so
    `V_{deg g} ⊆ ⟨G⟩ + K·g`, and any second ideal element of degree `≤ deg g` gives `LFD ≥ deg g + 1`.
    On the 409 instances `deg g = mk`, so this lower bound and the upper bound `max(mk+1, δ)` are both
    one more than the input degree. Read that way, the systems have one fall right above their input
    degree. The bound being met says nothing about closure growth. The paper should say so next to
    the claim, so that readers do not compare it with the two-summand profile, where the falls lie
    above `δ + 1`.

N4. **Minor:**
    - The depth-ten paragraph matches PKM "at every depth `k = 2..10`", but the introduction and
      Table `tab:profile` match it through `k = 11`.
    - The companion abstract writes the inverse result in terms of `d` while `thm:inverse` uses `n`.
    - The companion still says "Independent proof … review remain[s] outstanding" and is dated
      30 October.
    - The two papers cite Gorla–Mueller–Petit differently: JSC 109 (2022) in one, no volume and key
      `gmp21` in the other.
    - In the companion, `(a,b) = (2,3)` is undefined notation for field elements.

## 3. Status of the first review's findings

| # | Finding (short) | Status at `7c66226` |
|---|---|---|
| 1 | "≤ 5 through dimension 10" too wide | open (see N2) |
| 2 | selector-test negative control ineffective at equal depths | open |
| 3 | "three predictions … 7 not among them" vs six in the file | open |
| 4 | "every `k ≥ 8`" for `p = 11`, `40961`; `16777213`; `p = 193` | open |
| 5 | binary `⌊k/2⌋+1` row wrong for odd `k`, written-only | open |
| 6 | "exceptional set known completely" at `p = 193` | open |
| 7 | semicontinuity uses conjectured 605 | open |
| 8 | "generic" undefined | open |
| 9 | `thm:selector` sketch invalid as written | open |
| 10 | `thm:reduction` labels, omitted containment and `2 ≠ 0`, `W_k` scope | open |
| 11 | free-extension equality overstated; instance counts | open |
| 12 | direct sum / monotonicity not Lean declarations | open |
| 13 | internal inconsistencies (three vs four, 2^24 count, "ordinary ordinary", caption) | open |
| 14 | small count slips (43 chunks, 150 monomials, 80 systems, 1678/1775, row 13) | open |
| 15 | omitted Lean hypotheses; wrong namespace | open |
| 16 | no cold `lake build` log | open |
| 17 | no curve in Lean | open (disclosed) |
| 18 | lift–checker input binding | open |
| 19 | selector-test data absent | open |
| 20 | companion scripts and data absent; 424 vs 513 | open |
| 21 | `ic-depth-step/`, `ic-bounds-genericity/` absent | open |
| 22 | cited `docs/IC_*.md` absent | open |
| 23–26 | stale `LIFT.md`, `U1`, `N5`, experiment READMEs | open |
| 27 | attribution (HKY15, Petit–Quisquater/Diem, "automorphism") | open |

## 4. Recommended order of fixes

1. Do the text-only corrections. They change no result: findings 1–8, 9, 10, 13–15 and new N1–N4.
2. Commit the missing evidence (findings 19–21), or remove the claims that rest on it. This matters
   most for the 149,122-system test of `thm:selector` and for every number in Section 6 of the
   companion.
3. Commit a cold `lake build` log (16) and the lift `manifest_sha256` (18).
4. Refresh the stale notes (22–26).

## 5. Fixes applied (same day)

All edits are in the texts; no theorem or number changed except where a count was wrong.

| Findings | Fix |
|---|---|
| 1, N2 | "≤ 5 through dimension 10" restricted to the grids computed in the abstract, introduction, Section 7, item 5, conclusion and `tab:twosided`; the inverse family is named as a Koblitz exception with value `⌊k/2⌋+1` |
| 2, 19 | The selector-bound test is stated to have no working negative control, citing the 292-generator check; it is also stated that the per-instance records were not kept |
| 3 | The prediction paragraph now says six predictions were frozen, that 7 at `k = 14` was among them, and that two remain untested; ordering is attributed to the recorded dates |
| 4 | Instances `p = 11`, `40961`, `16777213` and the `p = 193` orbit statement are qualified as polynomial systems, with the range of `k` on which each is a subgroup system |
| 5 | `tab:twosided`: the binary/Koblitz lower entry moved to "Written", stated for the worst case over dimension `≤ k` with the odd-`k` caveat |
| 6, 7, 21 | "Known completely" → "among the targets computed"; the degree-18 statement is qualified; the semicontinuity sentence is conditioned on the conjectured value 605; the degree-480 bound is labelled as having no records in the repository, and 52 CPU-hours is labelled an estimate |
| 8 | "Generic" defined in Section 2 as "found on random instances" |
| 9 | Selector sketch rewritten along U1 §2.5 (representatives of weight ≤ k+1, projector, induction on monomials) |
| 10, 11 | "Written" → kernel-checked reduction applied to a computed base; `W_k` defined for general `k_0, D`; `2 ≠ 0` and the containment (`ShiftWitness.normalForm_mem_map_W`) added; free-extension equality stated as an upper bound proved and equality computed (k = 8..12 on one instance, 8, 9 on two more) |
| 12 | Direct-sum form of `thm:bridge` and monotonicity in `thm:restrict` marked as written corollaries |
| 13, 14, 15 | Counts corrected (four instances; sixteen/four grid sizes; five/three free-extension records; 42 chunks; 161 monomials; 80 cases / 75 systems; grouping fields; row 13 prime); "ordinary ordinary"; `tab:bounds` caption; `0 ∉ G`; `k_0 ≥ 2`; `FiniteQuotient` namespace |
| 16 | Cold `lake build` run in this session; log committed as `experiments/ic-lean/build-cold-20261112.log` (see the PR for its outcome) |
| 17 | The Lean paragraph states that there is no curve in Lean and that the binary/Koblitz bounds are one theorem |
| 18 | The paper discloses that the lift and checker records name the input basis by incompatible hashes and lists the fields that agree. The records cannot be repaired, because the basis was not kept |
| 20 | Companion: Lean count 513; the reproduction section says that only `verify_lean_fall_degrees.py` is public; the Section 6 campaign is labelled as having no public records |
| 21 | `experiments/ic-depth-step/` regenerated with `scripts/run_ic_depth_step.py`: 34 rows; `dim V_3 = 2k+1`, `dim V_4 = 4k^2+172k−571` for `k ≥ 6`, low part `4k+87`, depth coverage 5/5/…/4/2, all as in the paper |
| 22–26 | Dated notes added to `formal/README.md`, `experiments/ic-lean/README.md`, `ic-bounds-binary/VERIFY.md`, N5, N6, U1, `LIFT.md`, the m3, m3-depth, d5, exact10, locus and k14 READMEs; `RESULTS_PLACEHOLDER` replaced; N5 column count corrected to 1471; `NOTES.md` added on the stale hashes in the depth-8 and depth-10 directories |
| 27 | HKY15 credited for introducing LFD; Diem (2011) and Petit–Quisquater (2012) cited; Frobenius called an endomorphism. The two new bibliography entries were written from memory and should be checked against the sources |
| N1 | The rank sentences now cite the direct `eval_rank_below_D` computations |
| N3 | A paragraph says that attaining `mk+1` is the trivial `δ+1` fall |
| N4 | PKM match "through k = 11" made consistent; the companion uses `n`; review status, date and `(a,b)=(z,z+1)` notation fixed; the GMP volume added to the main paper |

Not fixable from here: the missing selector-test records (19), the genericity denominator records (21), the companion's scripts and campaign data (20), and the lift/checker binding records (18). Each is now disclosed in the text.
