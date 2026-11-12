# Review of the papers and proofs, 2026-11-12

Reviewed at commit `705279f`. Scope: `paper/curve-lfd-bounds.tex`, `paper/binary-lfd.tex`, the four
written proofs in `docs/proofs/`, the degree-lift proof in `experiments/ic-bounds-k14/LIFT.md`, the
statements of the Lean development in `formal/`, and the records behind the numbers in the main paper.

## Verdict

No mathematical error was found in any theorem, written or kernel-checked. Every written proof that
was refereed holds as stated. The defects are of four kinds: claims worded more widely than their
evidence, statements and labels in the paper that lag the proofs, evidence the paper cites that is
not in this repository, and stale notes. Section 2 lists them in order of importance.

Two limits on this review:

- **Nothing was compiled.** No Lean toolchain was available. "Kernel-checked" was reviewed by reading
  the Lean statements and the committed audit logs, not by building.
- **No large closure was rerun.** Record values were recomputed from stored rows or read from record
  fields. Small cases were recomputed with independent code (Section 4).

### How to read the evidence labels

Each finding is marked **[checked]** when the consolidating reviewer confirmed it directly against
the files, or **[referee]** when it comes from one of the six referee reports in
`docs/reviews-20261112/` and was not re-derived. The referee reports are included unedited; the
scripts they mention were scratch files and are not in the repository.

## 1. What was confirmed

- **Arithmetic of the main paper [checked].** All binomial counts, dimension sums, codimensions, the
  recursion `p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k)` on every quoted profile, the closed form
  for `d_5`, the identity `4k^2+172k-571` at `k = 6..10`, the `k = 14` lift count
  `448,922 + 2·586,501 + 30,655 + 30,639 = 1,683,218`, the checker's `4,120,872` triples and the six
  sign-test p-values were recomputed and agree.
- **Lemma `lem:unit` [checked].** The proof is correct. Independent enumeration confirms: `4271` is a
  non-residue mod `65537` with `4271^3 + 4271 = 11720^2`; `73730` likewise mod `163841`; the solution
  set is empty at `k = 10, 12, 13` over `65537` and at `k = 13` over `163841`, and non-empty at
  `k = 16`. The quoted solution counts `1028` (`p = 65537`, `k = 13`) and `594` (`p = 114689`,
  `k = 13`) and the empty set for the `p = 1153` instance of `thm:exactsix` were reproduced.
- **`S_4` [checked].** The Lean definition `M3CarryLower.S4` equals
  `Res_X(S_3(x_1,x_2,X), S_3(x_3,t,X))` identically; it has 125 terms, `(3,3,3)` is its only
  weight-6 mask, and that coefficient is `-8(5t^3+3at+6b)`.
- **Theorems `thm:three`, `thm:carry`, `lem:degree`, `cor:allm` [checked].** The proofs in the paper
  are correct as written; `alpha - beta^2 = -2(3t^2+a)` was re-derived.
- **Degree lift (`LIFT.md` §§1–2) [checked].** Lemmas 1–5 and the proof of the theorem were read line
  by line; no gap. The three `k = 14` lift records, the six `k = 14` degree-6 closure records and the
  three checker records agree with the paper in every quoted field.
- **Lean statements [referee].** No `sorry`, `admit`, `axiom`, `native_decide`, `unsafe`,
  `implemented_by`, `opaque` or `extern` in `formal/`. Every declaration the paper cites exists, and
  each kernel-checked theorem states what the paper says, with the caveats in 2.3. Instance constants
  match. The axiom log has 513 entries, all within `propext`, `Classical.choice`, `Quot.sound`.
- **Written proofs [referee].** `U1` (the `k+2` bound) proves `thm:selector` as the paper states it.
  `N5` §§2–6 and `N7` §2 are correct. `N6` proves what it claims as proved; genericity itself is its
  Conjecture 6.2 and no theorem of the paper depends on it. The three written arguments of
  `binary-lfd.tex` hold line by line. Each referee reproduced small cases with independent code.
- **Records [referee].** The sweep counts, the validator's 1440 systems, depths 8 to 10, the lift
  validation totals (170 runs, 37 cases), the three-summand counts, the more-summand table (409 of
  412) and the exceptional-locus counts match their records, apart from the items in 2.4.
- **Self-test [checked].** `python3 scripts/ic_bounds_msum.py selftest` prints `SELFTEST PASSED`
  (344 three-summand records compared, 0 inconsistent).

## 2. Findings

### 2.1 Claims wider than their evidence

1. **"Binary and Koblitz values are at most 5 for dimensions up to 10" [checked].** Abstract,
   introduction (l.557), Section 7 (l.1069, l.1507), conclusion (l.1613), Table `tab:twosided`.
   `thm:inverse` of the companion gives `n/2 + 1 = 6` for the full-field base at `n = 10`. The claim
   is true of the seeded grids, which contain no such base. Each sentence needs "on the grids
   computed" or the inverse family as a stated exception.
2. **Negative control of the selector-bound test [checked].** l.851–853: "random degree-2 generators
   do exceed the bound, so the test could have failed". The earlier review
   (`reviews-U1-upper-bound-redteam-20261110.md`, l.261–265) says its control exceeded the bound only
   at unequal depths and that "equal depths alone would not discriminate". The theorem and the
   149,122-system test are at equal depths. The U1 referee's equal-depth runs (292 random generators)
   found none above `k+2`. As stated, the test has no working negative control.
3. **Pre-registered predictions [checked].** l.376–385 says three predictions were committed and
   that the value 7 at `k = 14` "was not among them".
   `experiments/ic-bounds-freeext/predictions-2026-11-10.json` holds six predictions; `lfd-table`
   lists `"14": 7` and `d6-falls-at-14` ends "(LFD(14) = 7)". The sentence should be corrected; the
   file also predicts depths 17 to 22 and a base-depth statement (`d6-base-15`) that the paper does
   not report. The ordering claim (predictions before runs) rests on a self-attested string in
   `outcomes-2026-11-10.json`; it could not be checked from a single-commit clone.
4. **"LFD ≥ 5 at every depth k ≥ 8" for the two instances [checked].** l.527–529, l.1328–1357.
   `p = 11` has no subgroup of order `2^k` for any `k ≥ 2`, and `t = 7` is not a rational abscissa
   for `(a,b) = (2,5)`. `p = 40961` is a subgroup system only for `k ≤ 13`. The Lean statements are
   about the polynomial systems and are correct; the paper gives the extension-field caveat under
   `thm:restrict` for `12289` and `65537` but not here. The same applies to `p = 16777213` among the
   lift validation primes (`v_2(p-1) = 2`) and to "every `k ≥ 8`" at `p = 193` in item 1 of
   Section 7.
5. **Table `tab:twosided`, binary lower bound [referee].** "`≥ floor(k/2)+1`, inverse family,
   `k ≥ 6`" fails for odd `k`: the companion theorem gives 3 for odd `n`, and
   `ic-bounds-binary-law/table.jsonl` has value 3 at `k = 7, 9, 11`. The entry sits under "Proved
   lower" beside Lean entries although it is written-only.
6. **"The exceptional set is known completely" at `k = 6`, `p = 193` [referee].** l.1430. What
   exists is certified LFD on 1276 targets plus one sheet at `D ≤ 4`. "On no curve of degree at most
   18" carries little information: there are 190 monomials of degree at most 18, so 190 points in
   general position lie on no such curve.
7. **Closure-dimension semicontinuity [referee].** l.1441–1443 uses "generic dimension 605" as a
   fact to conclude failure in both directions; that value is the conjecture. Only "special below
   generic" on the line `b = t = 0` is rigorous.
8. **"Generic" is undefined [referee].** In the abstract, Table `tab:twosided`, the growth conjecture
   and the definition of `d_D(k)` it means "random instance computed".

### 2.2 Paper text that lags the proofs

9. **Proof sketch of `thm:selector` [checked].** l.839–845. "Every element of the ideal supported on
   the fibre is derived within weight `k+2`" is false as written (such elements have weight up to
   `2k`), and "summing over alpha" is not a valid step at budget `k+2`. `U1` §2.5 uses a monomial
   induction with a projector, which is correct. The sketch should follow §2.5.
10. **`thm:reduction` statement [checked for labels, referee for the rest].**
    - l.562 and l.1045 call it a "written reduction theorem"; the theorem is labelled kernel-checked.
    - The stated consequence is `1 ∉ V_4(k)`, but l.346–349 and l.1343–1351 use the containment
      `V_D(k) ⊆ W_k(S)` (`ShiftWitness.normalForm_mem_map_W`), which the statement omits.
    - The hypothesis `2 ≠ 0` of the Lean statement is omitted.
    - `W_k(S)` is defined for depth 7 and degree 4 and then quantified over all `k_0, D, m`.
11. **"`V_4(k)` is the free extension of `V_4(7)`" [referee].** l.1043, l.1463. Only the inclusion is
    proved (for `k ≥ 8`); equality is a computation, and the records give `k = 8..12` on one instance
    and `k = 8, 9` on two more, not "`k = 7..12` on three instances" (l.1047).
12. **Direct-sum form of `thm:bridge`; monotonicity in `thm:restrict` [referee].** Both are correct
    corollaries but are not Lean declarations; only the iff and `lfd_gt_of_le` are.
13. **Internal inconsistencies [checked].**
    - l.1472 says `dim V_6(14) = 448,922` "on three instances"; l.314, l.355 and the records say four.
    - l.1289: "eight with `2^20` or `2^21` grid points, two with `2^24`". Table `tab:msum` gives eight
      at `2^20`, eight at `2^21`, and four at `2^24`.
    - l.341–344: free extension tested on "four" and "two" instances; there are five `D5-*` records.
    - l.719–720: "ordinary ordinary".
    - Caption of `tab:bounds`: "the bound `mk+1` for the two-summand direct encodings".
14. **Smaller count slips [referee].** "43 `decide +kernel` chunks" (42 theorems per file); "about
    150 monomials" (161); "80 new systems" (75 distinct); "1678 of 1775 depend only on `(n,k)` and
    the target class" (the record also groups by two further fields); `tab:profile` row 13 "as
    `k = 11`" inherits `p = 12289` though the `k = 13` instances are over `40961`.
15. **Lean hypotheses the paper omits [referee].** `0 ∉ G` in the general part of `thm:carry`;
    `2 ≤ k_0` in `GrowthSchema.growth_schema`. The paper's
    `QuotientLowerBound.closure_le_quotient_span` is `FiniteQuotient.closure_le_quotient_span`.

### 2.3 What "kernel-checked" covers

16. **No from-scratch build evidence [checked].** `experiments/ic-lean/build.log` is a warm-tree
    replay; the certificate modules do not appear in it. All 110 source hashes in
    `verification.json` match the tree [referee]. A cold `lake build` log should be committed.
17. **No curve in Lean [referee, disclosed by the paper].** There is no group law and no proof that
    `S_3`, `S_4` are summation polynomials; `t` need not be an abscissa. `binary_descent_lfd_bound`
    and `koblitz_descent_lfd_bound` are one generic Boolean theorem (the second is `b := 1`), and
    `BinaryS3Bounds.normalized` is only partly reduced modulo `Z^2 - Z`. The bounds are unaffected.
18. **Lift input binding [checked].** The lift record identifies its input by `vd.basis_id` (a hash
    over name, size and digest of each file); the checker record stores `manifest_sha256` and
    per-file digests without sizes. The two cannot be matched from the records. Checkpoint slot,
    `source_id`, rank and free-column counts do agree. Recording `manifest_sha256` in the lift
    record would close this.

### 2.4 Evidence cited but not in the repository

19. **Selector-bound test [checked].** `experiments/ic-bounds-upper/` and
    `experiments/ic-bounds-upper-redteam/` are absent. The counts 149,122 / 132,303 / 2570 appear
    only in the prose of the earlier review, which says per-instance rows were not kept.
20. **Companion paper [checked].** Eight of the nine scripts in the reproduction section of
    `binary-lfd.tex` are absent (only `verify_lean_fall_degrees.py` exists), as is
    `experiments/ic-binary-inverse/`. No record was found for any number in its Section 6
    [referee]. It reports 424 audited declarations against 513 elsewhere.
21. **Other absent directories [referee].** `experiments/ic-depth-step/` (the five depth-step
    curves, l.1019) and `experiments/ic-bounds-genericity/` (the denominator degree `≥ 480`, called
    "exact" at l.1446; the 52 CPU-hours is an estimate).
22. **Cited guides [checked].** `formal/README.md` and the proof notes cite twelve `docs/IC_*.md`
    files; only `IC_FALL_COUNTEREXAMPLES.md` is present.

### 2.5 Stale notes

23. **`LIFT.md` [checked].** The header and §9 say nothing has run on a real GPU; §9a and the paper
    report three `k = 14` GPU lifts.
24. **`docs/proofs/U1-upper-bounds.md` [referee].** Not edited since the 2026-11-10 review; its
    findings F1, F2, F5, F6, F7, F9, F10, F11 remain open in the note (three were fixed in the paper
    only). None affects the proof of the theorem.
25. **`docs/proofs/N5-reduction.md` [referee].** §9 still says the base hypotheses are not
    kernel-replayable; "829 x 3473" should be 1471 columns; it does not say the primary instance is
    a proper ideal.
26. **Experiment READMEs [referee; first two checked].** `ic-lean/README.md` says 339 audited
    declarations; `ic-bounds-freeext/README.md` ends with `RESULTS_PLACEHOLDER`; `ic-bounds-m3`,
    `ic-bounds-m3-depth`, `ic-bounds-d5`, `ic-bounds-exact10`, `ic-bounds-k14/RELAUNCH.md` and
    `ic-bounds-locus` still list as "not run" runs that the paper reports. Source hashes in
    `ic-chain-depth8/verification.json` and `ic-chain-depth10/verification.json` no longer match
    the tree for three files.

### 2.6 Attribution

27. **[referee, from memory, not checked against the sources.]** Last fall degree was introduced by
    Huang, Kosters and Yeo (CRYPTO 2015), which the main paper cites as `hky15` but credits the
    definition to the 2018 paper. Petit–Quisquater and Diem are not cited for linear-subspace factor
    bases. l.500 calls Frobenius "a curve automorphism"; it is an endomorphism.

## 3. Observations that are not defects

- **`thm:restrict` already gives part of `thm:reduction` [referee].** With `1 ∉ V_4(7)`, restriction
  gives `1 ∉ V_4(k)` for every `k ≥ 7`. What the reduction theorem adds is the containment and the
  dimension bound.
- **The depth-7 certificate is strictly sufficient [referee].** At the C2-class point `p = 257`,
  `(148, 227, 143)`, the certificate fails for every `S`, yet `1 ∉ V_4(8)` at depth 8. The referee
  also argues the certificate holds iff 1 is outside the least depth-7 subspace containing `g` and
  closed under both translations, which would extend N6 Theorems 3.3 and 6.1 to `k ≥ 7`.
- **The lift proof does not need closedness if `--include-low` is used.** The `k = 14` runs did not
  use it, so they rely on the closure records' `closed: true`, which the checker confirms on all
  three inputs.

## 4. What was recomputed independently

| Item | By | Result |
|---|---|---|
| `lem:unit` hypotheses and solution counts at five (p, k) | consolidating reviewer | agree |
| `S_4` against the resultant; `(3,3,3)` coefficient | consolidating reviewer | agree |
| Bridge, 36 small cases; depth reduction, 18 instances; `V_4(6..8)`, `W_8`, `W_9` on both certified instances; N7 witness at E = 30, 46, 47 | N5/N7 referee | agree |
| `thm:selector`: 7,235 + 1,400 systems at `k = 2`, 120 at `k = 3`, 194 in the quotient model; Lemma 2.3 on 7,825 polynomials | U1 referee | no violation |
| `thm:three`: 9,431 systems at `k = 2, 3`, `p ≤ 17` | N6 referee | no `1 ∈ V_3` with `3t^2+a ≠ 0` |
| `thm:inverse` at `n = 1..8`; `cor:scaled` on eight pairs; `thm:uniform` on 260 systems | binary referee | agree |
| Evaluation rank 4068 at `k = 14`, `p = 65537` random curve | numbers referee | agree |

## 5. Not reviewed

The Lean build and the `decide +kernel` certificates; the internals of the `N5Cert*` and `ExactK6*`
modules; closures at `k ≥ 10`; `U1` §§3–7 beyond one finding; the pinpoint citations to HKYY
Prop. 2.8, Caminata–Gorla Thm 2.9 and Petit–Kosters–Messeng Table 4, which could not be checked
offline; novelty.
