# Referee report: `paper/binary-lfd.tex` (companion manuscript), `paper/proof-status.json`, and the binary/Koblitz parts of `paper/curve-lfd-bounds.tex`

Repository: `/home/user/lfd-summation-systems` (commit 705279f), read-only. Referee scripts and logs:
`/tmp/claude-0/-home-user/c25f3f7b-b111-516d-a596-5c5b49ca2dae/scratchpad/ref/` (`lfd.py`, `inverse.py`, `uniform.py`, `scaled.py`, `naf.py`, `rewrite.py`, `repro.py`, `inv10.py`).

## Verdict

**Major revision.** I found no fatal mathematical error. The three load-bearing written arguments
(selector bound `min(k,l)+2`; inverse-family upper bound; even lower obstruction) survive a line-by-line
check, and every exact value I could reach independently agrees (inverse family n = 1..8, scaled
subfields up to F_64 in F_4096, 260 random systems for the uniform bound, nine k = 8 records).
The manuscript is honest that Theorems `thm:uniform` and `thm:inverse` are written, not kernel-checked, and this matches
`proof-status.json`. What blocks acceptance:

1. the finite-evidence section (613 systems etc.) and the reproduction section of `binary-lfd.tex`
   cannot be traced to anything in this repository (data and 8 of 9 scripts are absent);
2. `curve-lfd-bounds.tex` repeatedly states "binary and Koblitz values are at most 5 for dimension up
   to 10" without restriction, which the companion's own theorem contradicts (value 6 at dimension 10);
3. stale or inconsistent audit counts, and an understatement of the open formal obligations;
4. the termination argument of Lemma `lem:rewrite` is a sketch, not a proof;
5. positioning against HKY/HKYY/Petit–Quisquater/Diem is missing or mis-attributed.

Lean was **not rebuilt** (no `lake` on this machine); all statements about the kernel rest on reading the
sources and the stored logs in `experiments/ic-lean/`.

## Findings

### A. Mathematics of `binary-lfd.tex`

1. **[fixable gap] `binary-lfd.tex:236-249`, Lemma `lem:rewrite`, termination.** The case analysis is
   informal ("if the position just left of the pair was occupied, it becomes an opposite-sign adjacency";
   "must stop within n steps: revisiting a starting position would already encounter one of the fixed
   opposite-sign digits"). The lexicographic measure (weight, occupied adjacencies) does not decrease in
   the "propagate" case, so the inner propagation needs its own explicit argument, separately for n even
   (moving digit lands on the fixed minus digit and cancels) and n odd (it becomes a mixed adjacency).
   The existence statement is true for n = 2..10: my exhaustive search (`rewrite.py`) finds, from each of
   the 3^n signed words, a path of the stated rules to a cyclic nonadjacent word that never exceeds the
   starting weight (0 failures). The status file already lists this as open
   (`all-length-polynomial-normalization`). Repair: write the propagation as a lemma with an explicit
   invariant, or replace it by a direct construction (reduce the residue mod 2^n-1 to its cyclic NAF
   digit by digit and bound the weight of each intermediate word).

2. **[unclear] `binary-lfd.tex:329-357`, Proposition `prop:lower`.** The proof is correct but stated as
   an invariant on `V_{F,d}`, which hides that it needs no induction: every element of `V_{F,d}` is a sum
   of inputs and of products `z*h` with `z` a variable and `h` in the ideal of ordinary degree <= d-1; for
   any such `h`, Lemma `lem:positive` kills the coefficient of `U/u_i`, and the preimage analysis shows
   `NF(z*h)` has no `U` term. Say this. Also say explicitly that the normal form of a mixed monomial
   (both `u_i` and `v_i`) is covered by the "cancel" step in the proof of Lemma `lem:positive`.
   Checked by hand; consistent with my computations: `dim I_{<=3} - dim V_3 = 2` at n = 6 and
   `dim I_{<=4} - dim V_4 = 2` at n = 8 (exactly `U-V` and its rotation).

3. **[typo/unclear] `binary-lfd.tex:182-183`.** "The Moore matrices are invertible because their columns
   are F_2-linearly independent" should read: because the basis elements are F_2-linearly independent
   (the Moore criterion); K-independence of columns is the conclusion, not the reason.
   Lean `MooreMatrix.determinant_nonzero` proves the right statement by Dedekind independence.

4. **[typo] `binary-lfd.tex:297`.** "The short case d=1 is the carry relation itself": for n = 2 one needs
   two carry rows (`u_1 v_0 - u_0` and `u_1 v_0 - v_1`), not one. Harmless (budget 3).

5. **[unclear] `binary-lfd.tex:391-399`.** The first-fall-degree paragraph ("FFD three", "retaining
   redundant ambient rows can produce FFD two") is asserted without proof or reference to a record.
   **Unverified** by me.

6. **[unclear] `binary-lfd.tex:496-497`.** "with m = 2,...,9" uses `m` for a dimension, while lines 451,
   543 and 576 use `m` for the number of summands. Fix the notation.

7. **[no defect found] Theorem `thm:uniform` (lines 91-159).** Checked each step: the fibre polynomial
   `(x+t)^2 y^2 + xt y + (xt)^2 + b` is right, is nonzero for b != 0, and is F_2-affine in y; selector
   defects cost k+1; fibre rows cost k+2; the section has degree <= k+1; the monomial induction is
   sound. 260 random instances (n = 3..6, dimensions 1..3, random t and b, common and distinct spaces):
   no violation; the bound is attained at min(k,l) = 1 (LFD 3).

8. **[no defect found] S_3, descent and degrees (lines 46-58).** `S_3 = (xy+xt+yt)^2 + xyt + b` is the
   correct binary third summation polynomial (the coefficient a does not occur); with separate variable
   blocks each reduced coordinate has bidegree <= (1,1) because squaring is F_2-linear. The raw bound
   `max(4, min(k,l)+2)` (lines 161-164) follows as stated.

9. **[no defect found] Lemmas `lem:moore`, `lem:rigidity`, `lem:bridge`, `lem:positive`,
   Corollaries `cor:scaled`, `cor:envelope`.** Proofs checked by hand. Rigidity verified exhaustively for
   n = 2..14, positive-word isolation for n = 1..12 (`naf.py`). `cor:envelope` is correct as worded
   (dimension *at most* k); see finding 13 for how the other manuscript restates it.

### B. Evidence level against `proof-status.json` and the Lean sources

10. **[overclaim] `curve-lfd-bounds.tex:1527-1529`.** "The kernel specialization of the selector theorem
    is *the* open formal item of the companion manuscript." `proof-status.json` lists four open
    obligations (Moore equivalence and input filtration; all-length normalization; even lower bound;
    selector specialization), and `main_theorem_fully_formalized: false`. Repair: "one of four".

11. **[overclaim/unclear] `curve-lfd-bounds.tex:792-809`, Table `tab:twosided`.** The binary lower bound
    `>= floor(k/2)+1` sits in a column headed "Proved lower (family)" next to Lean-checked entries; it is
    a written proof with two open formal obligations. The caption says "written entries are proved in the
    companion manuscript". Use "written, not independently reviewed" as the status file does.

12. **[overclaim, the most serious inconsistency] `curve-lfd-bounds.tex:477-478, 557-558, 1069-1070,
    1507-1508, 1613-1614`.** "Binary and Koblitz values are at most 5 for dimensions up to 10" is stated
    without restriction in the abstract-level summary, the results list, Section `sec:missing` and the
    conclusion. By `thm:inverse`/`cor:scaled` of the companion, the Koblitz curve b = 1 with the space
    F_{2^10} (full field, or the proper subfield of F_{2^20}) and the order-two target has value
    10/2+1 = 6. The sweep excludes k = n and has no scaled subfield of degree 10, which is the only
    reason the grid shows 5. `experiments/ic-bounds-binary-law/README.md` itself says so ("hold only on
    the grid"). Every one of these sentences needs "on the seeded grid, which contains no full-field or
    scaled-subfield base of dimension >= 10". (My own run confirms value 5 at n = 8 and the proof gives 6
    at n = 10; a direct n = 10 degree-5 closure was started but had not finished when this report was
    written: **unverified numerically**.)

13. **[overclaim] `curve-lfd-bounds.tex:801`.** "Binary, spaces of dimension k,l: >= floor(k/2)+1,
    k >= 6". For odd k the witness has dimension k-1; the companion (lines 447-451) is careful that the
    statement is for dimension *at most* k. Restate the table entry.

14. **[unclear] Lean fidelity of `binary_descent_lfd_bound` / `koblitz_descent_lfd_bound`
    (`formal/ICFallDegrees/CurveBounds.lean:52-71`, `BinaryS3Bounds.lean:95-139`).** The statement is
    faithful as a bound but is not a statement about a curve: `constants` and `coefficients` are
    arbitrary, there is no curve, no point and no coordinate space in the hypotheses, and the proof is
    `boolean_lfd_bound` (field equations plus input degree <= 3). The Koblitz theorem is the binary one
    with b := 1 and never uses a in {0,1}. Two differences from the text at
    `curve-lfd-bounds.tex:617-625`: (i) Lean's `normalizedConstraint` only replaces `x^2,y^2,z^2` by
    their reduced affine forms; it is not reduced modulo `Z^2-Z` and has degree <= 3, not <= 2 (it
    coincides with the paper's presentation when the three forms use disjoint variables and the target is
    constant, which is the two-summand case); (ii) the hypothesis `2 <= r` is missing from the status
    file line "binary/Koblitz descended S3 systems LFD <= r+1". The paper does say the bound "uses only
    the Boolean field equations and the input degrees"; it should also say that the two table rows are
    one theorem, and that "kernel-checked for binary and Koblitz curves" means "for any Boolean system of
    this shape".

15. **[typo/stale] audit counts.** `binary-lfd.tex:527` says 424 audited declarations;
    `experiments/ic-lean/verification.json` and `formal/README.md:45` have 513; `Audit.lean` has 513
    `#print axioms`; `experiments/ic-lean/README.md:6` says 339; `scripts/verify_lean_fall_degrees.py:74`
    asserts only `>= 424`. Pick one number and generate it.

16. **[accurate] `binary-lfd.tex:527-543`.** The verification-boundary paragraph matches the status file:
    all 20 `kernel_components` exist as theorems in the named files and appear in the stored axiom audit
    with only `propext`, `Classical.choice`, `Quot.sound`; I found no `sorry`, `axiom` or `native_decide`
    in `formal/ICFallDegrees/*.lean` by text search. `even_alternating_bridge`, `rigidity_of_modEq`,
    `lfd_bound_of_local` and `MooreMatrix.*` state what the manuscript says they state. They are
    conditional lemmas (hypotheses such as `InClosure F 2 (u_i^2 - u_{i+1})`); none is instantiated on the
    inverse system, as the status file says.

### C. Reproducibility and data

17. **[fixable gap, blocks the evidence section] `binary-lfd.tex:453-502` and `545-566`.** None of the
    numbers of Section 6 can be traced: 613 distinct systems, 702 = 240+240+222 cases, 18 skips, the
    histogram 7/1/50/381/174, 167/141, 1,803/1,474, 29 caps, 587, "traces through n = 12". No file in the
    repository contains seed 20261030 data; the directory `experiments/ic-binary-inverse` cited by
    `experiments/ic-bounds-binary-law/README.md:6` does not exist. Of the nine commands listed at lines
    549-565, eight scripts are absent (`reproduce_ic_paper.py`, `build_ic_paper.py`,
    `run_ic_arbitrary_bases.py`, `verify_ic_positive_naf.py`, `verify_ic_inverse_parity.py`,
    `run_ic_inverse_subfields.py`, `verify_ic_binary_m2_uniform.py`, `verify_ic_research_artifact.py`);
    only `verify_lean_fall_degrees.py` exists. The "SHA256 manifests" are absent. Internal arithmetic is
    consistent (20 (k,n) pairs x 4 spaces x 3 curves = 240; 240 x 3 = 720 = 702 + 18; histogram sums to
    613), but that is all I can say. Either ship the data and scripts or delete the section.

18. **[typo] `experiments/ic-bounds-binary-law/README.md:6-7` against `score.json`.** "80 new
    preregistered cases (60 distinct systems)"; `score.json` and my count give 75 distinct hashes (60 is
    the count without the 15 order-two rows). `curve-lfd-bounds.tex:1515` calls them "80 new systems".

### D. Citations

19. **[mis-attribution] `binary-lfd.tex:505-506`, `curve-lfd-bounds.tex:486`.** The last fall degree was
    introduced by Huang, Kosters and Yeo (CRYPTO 2015); the 2018 JSC paper of Huang, Kosters, Yang and
    Yeo is the Weil-descent follow-up. `binary-lfd.tex` does not cite the 2015 paper at all.

20. **[fixable gap] no comparison with the HKY/HKYY descent bounds.** Those results bound the LFD of a
    Weil-descent system independently of n under a zero-dimensionality hypothesis on the system before
    descent. The inverse family (descent of `xy+1`, a positive-dimensional ideal before field equations)
    has LFD n/2+1. The manuscript must state why there is no conflict and what its theorem adds; the
    status file lists this as an open review question. I state the HKYY hypotheses from memory:
    **unverified**.

21. **[omission] Petit–Quisquater (ASIACRYPT 2012) and Diem are not cited in either manuscript.** The
    first-fall-degree assumption for binary summation systems originates with Petit–Quisquater, and the
    F_2-linear-subspace factor base with Diem and FPPR; `curve-lfd-bounds.tex:497-498` credits
    "FPPR; Semaev" for it, and Semaev 2004 has no subspace factor base.

22. **[unverified] `binary-lfd.tex:210-211, 517-518`.** Attribution of cyclic (mod 2^n-1) nonadjacent
    uniqueness, with the even exceptional pair, to Preparata 1971. The reference exists; I could not
    confirm that it contains the cyclic statement. The usual reference for cyclic NAF is Clark and Liang
    (IEEE Trans. Inf. Theory, 1973-74). Check the primary text.

23. **[unverified]** `binary-lfd.tex:509-511` ("Section 6, Theorem 11 and Corollary 1" of
    Gorla–Mueller–Petit); `curve-lfd-bounds.tex:489, 1538` (HKYY Prop. 2.8); `1540` (Caminata–Gorla
    Thm 2.9, Ex. 4.3); `1567-1570` (Kosters–Yeo: FFD 2, NP-completeness). Bibliographic data of HKYY,
    Caminata–Gorla, GMP, FPPR, Semaev, PKM, Kosters–Yeo, FHJRV, FGHR, GGMP, Bosma look right from memory.

24. **[typo] `curve-lfd-bounds.tex:500-501`.** "Koblitz curves ... whose Frobenius is a curve
    automorphism": the 2-power Frobenius is an endomorphism of E defined over F_2; it is a bijection on
    E(F_{2^n}), not an automorphism of the curve.

## What I verified and how

All computations use my own code (`lfd.py`): GF(2^n) by shift-and-reduce, the normalized system built
from the formula for S_3, closure `V_D` in the ordinary polynomial ring over F_2 (all monomials of degree
<= D, field equations as inputs, variable-wise products of rows of degree < D, graded echelon basis),
`I_{<=D}` from the rank of evaluation on the zero set found by scalar S_3. LFD = least d with
`V_e = I_{<=e}` for all e >= d; the tail is certified by: `V_D = I_{<=D}` and evaluation of degree < D
onto the zero set imply `V_e = I_{<=e}` for all e >= D (short induction, my own). No repository code was
imported.

| Claim | Result |
|---|---|
| `thm:inverse`, n = 1..8, Boolean coordinates over F_2 | LFD 3,3,3,3,3,4,3,5; matches 3 / max(3,n/2+1) |
| deficiency of `V_{n/2}` | 2 at n = 6 (390 vs 392) and at n = 8 (4588 vs 4590) |
| `cor:scaled`: (n,d) = (4,2),(6,3),(6,2),(8,4),(8,2),(12,6),(12,4),(12,3), alpha = s = 1 and random | 3,3,3,3,3,4,3,3; all match |
| `thm:uniform`, 260 random systems | 0 violations |
| `lem:rigidity`, n = 2..14 | odd: no collision; even: exactly the two alternating pairs |
| `lem:positive`, n = 1..12 | 0 violations |
| `lem:rewrite` (existence), n = 2..10 | 0 unreachable words |
| degree-3 identities at lines 228-230, bridge path of `lem:bridge` | by hand |
| n = 10 inverse, `V_5` against `I_{<=5}` | not finished; unverified |
| Lean build | not run; sources and stored logs read only |

Numbers traced to records (path, result):

| Quoted | Where | Record | Result |
|---|---|---|---|
| m2 histograms k = 8: 2:3, 3:70, 4:119; k = 9: 2:16, 3:30, 4:146; k = 10: 2:17, 3:4, 4:133, 5:30 | `ic-bounds-binary/README.md` | `experiments/ic-bounds-binary/results-m2.jsonl` (568 rows) | match |
| "28 systems with value 5, all at k = 10" | `curve-lfd-bounds.tex:1508` | same file | match: 30 cases, 28 distinct hashes, n = 11 (20), n = 12 (10) |
| "exactly those with more than 470 solutions" | `:1509` | same | match: range 472..1056; largest count with LFD <= 4 is 352 |
| max 5 through k = 10 on the grid | `:801` | same | match (max 5) |
| written bound never exceeded | — | `bounds.written_min_k_l_plus_2` | match, 0 violations |
| F_{2^8} in F_{2^16}, 255 solutions, value 5 | `:1522-1523` | `ic-bounds-binary-law/subfield-check.jsonl` | match; also 63 solutions, value 4 at (12,6), same deficiencies as my run |
| 80 new systems; 47 discriminating | `:1515-1517` | `results-prereg.jsonl`, `score.json` | 80 rows match; 47 = 11+12+12+12 match; 75 distinct (finding 18) |
| 11 hits / 36 misses and 36 / 11 | `:1517` | `score.json`, README table | match |
| k = 9 trace-zero hyperplane, about 490 solutions, 5 in 11 of 11 | `:1518-1519` | `results-prereg.jsonl` (n = 10, k = 9, 11 rows LFD 5, |Z| 482..499) | match |
| unit ideals reach 4 at k = 7 | `:1520` | law README (three systems at (12,7)) | README only; rows not opened |
| 1678 of 1775 | `:1525` | law README section 2 | README only; not recomputed |
| 985 binary systems, 1440 total | `:399-401` | `experiments/ic-curve-validation/REPORT.md:9, 300` | match (report text) |
| 22 of 30; 31 of 38; 20 of 22; 13.2 vs 6.4; -0.18, -0.28 | `:957-965` | `experiments/ic-curve-sweep/base-class-analysis.json` | present (13.19, 6.35) |
| nine k = 8 rows (LFD 2, 3, 4; n = 11, 12, 14) | — | `results-m2.jsonl` case ids 66338b02…, 9d21fb1e…, 416a9613…, feee4b35…, 904c6a87…, 675cbe09…, afbb3edc…, cb75f68b…, 28e0e986… | recomputed: LFD and |Z| equal in 9 of 9; codimensions at degrees 2..4 equal where comparable |
| 613, 702, 240/240/222, 18, 7/1/50/381/174, 167/141, 1,803/1,474, 29, 587 | `binary-lfd.tex:459-498` | none in repository | **not traceable** (finding 17) |
| 424 audited declarations | `binary-lfd.tex:527` | `experiments/ic-lean/verification.json` | **mismatch**: 513 |
