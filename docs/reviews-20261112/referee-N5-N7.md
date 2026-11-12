# Referee report: N5 (depth reduction) and N7 (proper-ideal lower bound)

Objects: `docs/proofs/N5-reduction.md`, `docs/proofs/N7-proper-ideal-lower-bound.md`,
and `thm:bridge`, `thm:reduction`, `thm:restrict` in `paper/curve-lfd-bounds.tex`
(repository `/home/user/lfd-summation-systems`, commit 705279f, read-only).

## Verdict in one paragraph

I found no fatal gap. The mathematical arguments of N5 Sections 2-6 (bridge, translation
identity, two-level model, `W_closed`) and N7 Section 2 are correct as written, and every
numerical base fact I could recompute independently came out as claimed. The defects are in
the statements and labels: the paper's `thm:reduction` omits a hypothesis and the conclusion
that the later text actually uses; "every depth k >= 8" is not a statement about a subgroup
factor base for p = 11 at any depth nor for p = 40961 beyond k = 13; the paper still calls the
reduction "written" in two places and claims an equality where only an inclusion is proved;
the N5 note has a stale formal-status section, a wrong column count and an incorrect remark.
Lean was NOT compiled (no toolchain in this environment), so "kernel-checked" is unverified by
me; I read the Lean statements and the main proof terms only.

## Findings

Severity scale: fatal gap / fixable gap / unclear exposition / typo.

**F1. Unclear exposition (overclaim of scope). Paper l.1328-1341, l.1352-1357, l.527-529;
N7 Sections 1 and 6.** "lfd >= 5 at every depth k >= 8 for p = 40961 and for p = 11".
Section 2 of the paper defines the presentation for `2^k | p-1`.
- p = 11: `p - 1 = 10`, so no subgroup of order `2^k` exists for any k >= 2. In addition
  `t^3 + a t + b = 10` is a non-residue mod 11 for (a,b,t) = (2,5,7), so t = 7 is not the
  abscissa of an F_11-rational point of the curve. The p = 11 theorem is a statement about a
  formal polynomial system, at every depth.
- p = 40961: `p - 1 = 2^13 * 5`. The statement is about a subgroup system only for
  8 <= k <= 13; for k >= 14 it concerns the same polynomials with solutions in an extension
  field.
The paper makes exactly this caveat for p = 12289 under `thm:restrict` (l.1404-1407) but not
here. Repair: repeat the caveat at `thm:reduction`, in the contributions list and in N7; call
p = 11 a formal test of the certificate pipeline, not an instance of the presentation.
(The Lean statements themselves are fine: they quantify over `PrimeM2General.inputs` and make
no subgroup claim.)

**F2. Fixable gap (statement, not proof). Paper `thm:reduction`, l.1311-1326.**
(a) The "Consequently" clause needs `2 != 0` in K (Lean: `two_ne`); the statement does not
say so (`thm:bridge` does).
(b) The statement concludes only `1 notin V_4(k)`. What the text then uses is the containment
`NF(V_D(k)) subseteq Phi(W_k(S))` (l.346-349 "V_5(k) subseteq W_k(V_5(11))"; l.1343-1351 "an
element of the low box outside S stays outside W_k(S)", which is useless without the
containment). That containment is `ShiftWitness.normalForm_mem_map_W` and N5 Theorem 7.2; it
is not part of the stated theorem. Add it as a conclusion.
(c) `W_k(S)` is defined just before the theorem only for base depth 7 and degree 4, and the
theorem then quantifies over "all k_0, D, m". Define `W_k(S)` with k_0 and D.

**F3. Unclear exposition (label mismatch written vs kernel-checked).**
- Paper l.562 "explained for k >= 7 by a written reduction theorem" and l.1045 "by a written
  reduction theorem" contradict the label "kernel-checked" on `thm:reduction`.
- Paper l.1043-1047 and l.1463: "V_4(k) is the free high-level extension of V_4(7)" and "the
  identity 4k^2+172k-571 is explained per instance for k >= 7". Proved (given the base) is only
  `V_4^q(k) subseteq W_k`, i.e. `dim <= 4k^2+172k-571` for k >= 8. Equality is a computation
  (k <= 12, dimension to k = 13), and the dimension count itself (N5 Theorem 7.2) is
  written-only. Say "bounded above by, and equal on the computed range".
- Paper `thm:bridge`: the "Equivalently V_D(F) = (<G> cap M_{<=D}) (+) V^q_D(F)" form is not
  in Lean; `chain_mem_closure_iff` is the iff only (N5 Cor. 3.3 says so). The derivation is
  three lines and correct, but the label covers more than the kernel does. Also `M_{<=D}` here
  means all polynomials of degree <= D, while elsewhere in the paper it is the squarefree
  weight-<=D monomials (e.g. `|M_{<=7}| = 63,004`); N5 writes `M^{ord}_{<=D}` for this reason.
- "S = V_4(7)" (paper l.1330, N7 Section 3): Lean defines `S` as the span of the certified
  echelon rows and proves `1 notin S`, `g in S`, `W_8(S)` closed. It does not prove
  `S = V_4(7)`; that identification is a computation (it is not needed for the bounds).

**F4. Unclear exposition (N5 note is stale and partly inaccurate).**
- N5 Section 9: "none of which is discharged in Lean", "not kernel-replayable without
  `native_decide` at present", "To register the modules, add ...", "Next concrete step" are
  all superseded by `N5Cert*` (modules are registered in `formal/ICFallDegrees.lean` l.78-100;
  audit lines exist). The note's own header promises that Section 9 names what is
  kernel-checked; it no longer does so correctly.
- N5 summary item 4 and Theorem 7.1 give `LFD >= 5` for a unit ideal; N7 Section 1 states the
  primary instance p = 40961 is a proper ideal, so that clause is vacuous for the primary
  instance. N5 does not say this.
- N5 cites `docs/IC_BOUNDS_PROGRAM.md`, which does not exist in the repository.
- Typo, N5 l.390, l.550, l.569: "829 x 3473". Depth 7 has `binom(14,<=4) = 1471` columns
  (paper l.1331, N6, N7 "829 rows, 642 free columns"; my recomputation). 3473 is
  `binom(14,<=5)`.
- N5 Section 3, "Remark (why D >= 2)" is wrong. For D <= 1, `<C> cap {deg <= 1} = 0` and the
  decomposition of Cor. 3.3 holds trivially (checked at D = 1). More generally the hypothesis
  `deg g <= D for g in G` of Theorem 3.1 is superfluous: in the transport step the chain
  equations used inside a derivation of budget `deg f <= D` already have degree <= D. So
  `D >= 2` in `thm:bridge` and `thm:restrict` is sufficient but not needed. Harmless.

**F5. Fixable gap (hypotheses of the dimension argument for proper ideals). N5 Cor. 7.3;
paper l.60-63, l.1454-1455.** Note first that N7 contains no dimension or evaluation-rank
argument; its lower bound is by an explicit witness. The evaluation argument is N5 Cor. 7.3,
labelled written-only (correctly: "the evaluation fact is not formalized").
- The argument is correct when `2^k | p-1`: then `R_k/<C_k>` is `K^{H^2}`, every ideal is the
  ideal of its zero set, `dim(I cap M_{<=D}) = |M_{<=D}| - rank(ev_Z on M_{<=D}) >=
  |M_{<=D}| - |Z|`, and a reduced ideal element of degree <= D outside `V_D` gives
  `lfd >= D+1` directly from the definition. I confirmed the equality exactly on 18 small
  instances (below). But `2^k | |K^x|` appears only inside a parenthesis; it is a hypothesis
  of the corollary and fails for the primary instance at k >= 14. State it. Without it Z must
  be counted over the splitting field.
- The corollary must also list its real hypotheses: (B2) and `g in S`, not (B1).
- `|Z| <= 2^{k+1}` ("by N4a", not in the files under review) needs a nonsingular curve: at
  `x_1 = t` the polynomial `S_3(t, x_2, t)` vanishes identically iff `t^3+at+b = 0` and
  `3t^2+a = 0`. State nonsingularity.
- The bound `|Z| <= 2^{k+1}` cannot give an all-depth statement: `binom(2k,<=4) - 2^{k+1}`
  falls below `4k^2+172k-571` for large k. Cor. 7.3 is a per-depth criterion needing the
  actual `|Z|`. The paper's wording "proper ideals when |Z| is below the dimension gap" is
  accurate; the all-depth result for the proper instance rests on the N7 witness only.

**F6. Unclear exposition (redundancy; not an error).** `thm:restrict` already gives
`1 notin V_4(k)` for all k >= 7 from `1 notin V_4(7)` alone, with no depth-8 certificate and
no `2 != 0`. So the exclusion clause of `thm:reduction` is subsumed; what the depth reduction
adds is the containment `V_4^q(k) subseteq W_k(S)` (dimension bound). The N7 witness bound
also follows from restriction: the witness uses only levels j < 7, is fixed by
`u_{i,j} -> 1` (j >= 7), so `v in V_4(k)` would force `v in V_4(7)`. That route would need a
depth-7 closedness check of S (14 x 115 memberships with wrap at 2^7) in place of the 3220
depth-8 memberships. This is my own argument, not machine-checked. The paper should say which
conclusion actually needs `thm:reduction`.

**F7. Unclear exposition. N7 Section 3 table.** "U_47 cap M_{<=4}: 1 basis vector outside S"
is basis-dependent (my basis has 326 and 361 such vectors for p = 11 and p = 40961). The
invariant, which I confirmed, is `dim((U_47 cap M_{<=4}) + S)/S = 1`. N7 concedes the basis
dependence only for E = 54.

**F8. Typo-level. Paper `thm:restrict`.** "for unit ideals the last fall degree is
nondecreasing in the depth" is not a Lean theorem as such; Lean has `lfd_gt_of_le` (unit at
the larger depth assumed). The claim follows since a unit ideal at depth k+1 restricts to one
at depth k and `lfd = min{d : 1 in V_d}` for unit ideals (`hasLFDBound_of_one`,
`one_of_hasLFDBound`). Fine, but it is a corollary, not one of the cited declarations.

## Steps verified, and how

By hand (line by line):
- N5 2.1 normal-form calculus; Lemma 2.2 (both bit cases, the wrap at j = k-1).
- N5 Theorem 3.1 both directions, including that the product rule is charged on derived rows
  (`deg a + deg rho(p) <= D`, variable-at-a-time induction with `deg < D`); Cor. 3.3; Thm 3.4.
- N5 Lemma 4.2 (three cases), carry lemmas C1-C5, Lemma 4.3 (needs k > k0, which is assumed).
- N5 Def. 5.1 is a subspace and equals the direct sum over patterns; (ii).
- N5 Theorem 6.1 (iii): box, weight budget (+1 per step against `< D`), high-variable
  permutation of components, and the low-variable identity
  `piece q w' = piece 0 (T embed0 piece_q w) + piece e_i (T embed0 piece_{q^-} w)` in both
  carry cases and outside the box; that `embed 0 (piece q w)` is an admissible input of the
  depth-(k0+1) hypothesis. The conclusion "for every k >= k0+1" is direct (no induction on k
  is needed or used). (i) is correct.
- N5 Theorem 7.1 assembly; Theorem 7.2 count (`829 + 230(k-7) + 2 binom(2k-14,2) =
  4k^2+172k-571`, 1061 at k = 8, 119 low rows of W_8, 3220 = 14 x 115 x 2).
- N7 Section 2 steps 1-5; LFD >= D+1 from one ideal element of degree <= D outside V_D follows
  from the definition in paper Section 2 and from `HasLFDBound`.
- `thm:restrict`: substitution fixes the generator for k >= 2, sends chain equations to chain
  equations or 0, does not raise degree; direction of the implication is right.

Lean, read but not compiled: `Closure.lean`, `QuotientLowerBound.lean`,
`ChainClosureBridge.lean` (statements and bridge proofs), `DepthReduction.lean` (`W`,
`Closed`, `piece_translate_low/high`, `W_closed`), `ChainModelMap.lean` (assembled theorems),
`ShiftWitness.lean` (all), `DepthRestrictionM2.lean` (all), instance files for p = 40961.
The statements match N5 Section 9 and N7 Sections 1, 4, 5. No `sorry`, `axiom` or
`native_decide` by grep.

By my own computation (numpy, own engine; scripts `ga.py`, `t1.py`-`t5.py` in the scratchpad;
no repo code imported):
- Bridge: direct ordinary-ring closure (multiplying derived rows by all monomials within
  budget) against the quotient closure, 36 cases, p in {3,5,7,13,17}, k in {2,3},
  D in {1,..,4}: `dim V_D = dim(<C> cap deg<=D) + dim V_D^q` and `1 in V_D iff 1 in V_D^q`
  in every case.
- Depth reduction: 18 random small instances (p in {11,13,17,97}, k0 in {2,3,4},
  D in {3,4,5}). Whenever `W_{k0+1}(V_D(k0))` was closed, `W_k` was closed and contained
  `V_D^q(k)` at every further depth tested (up to k0+3); no violation. Non-closed bases occur
  (so the hypothesis is not vacuous) and then containment fails, as expected.
- Certified instances, p = 40961 (21223, 9887, 4748) and p = 11 (2,5,7): `dim V_4(6) = 605`,
  profile (0,0,2,111,605); `dim V_4(7) = 829` in 1471 columns, profile (0,0,2,115,829),
  `1 notin`; `W_7(V_4(6))` not closed; `W_8` has dimension 1061 in 2517 columns, is closed,
  excludes 1; `V_4(8)` computed directly has dimension 1061 and lies in `W_8`; `W_9` has
  dimension 1301 and is closed. All as claimed.
- N7 witness: `dim(U_E cap M_{<=4})` = 270, 503, 534 for E = 30, 46, 47; all inside S for
  E <= 46; codimension 1 modulo S at E = 47, and those elements are outside `W_8` at depth 8.
  Both instances.
- Evaluation fact: 18 instances (p in {17,41,97,193}, k in {3,4,5}): `dim(A g) = 4^k - |Z|`
  and `dim(I cap M_{<=D}) = |M_{<=D}| - rank(ev_Z)` for D = 2..6, all exact.
- 40961 is prime, the curve is nonsingular, t is an abscissa; the p = 11 facts in F1.

## Unverified

- That the Lean development compiles and that the `decide +kernel` chunks pass (no Lean
  toolchain here). `N5CertSound.lean` (1072 lines), `N5Cert.lean`, `N5CertWitness.lean` and
  the generated data files were not audited; I checked only that `S` is the span of the data
  rows and that the top-level statements have the claimed shape.
- "N4a" (`|Z| <= 2^{k+1}`), beyond my sketch in F5.
- The 131 single-shift candidates, `U_6`, `U_14`, `U_48..U_62` rows of the N7 table, the
  cofactor term counts and timings.
- All closures at k >= 10, the D = 5 free-extension claims (paper l.339-349), and equality
  `V_4(k) = W_k` beyond k = 8.

## Verdict per theorem

- **`thm:bridge`**: the iff is proved as written (N5 Section 3; confirmed numerically). The
  direct-sum form is correct but written, not kernel-checked. `D >= 2` is not needed.
- **`thm:reduction`**: the general part (`W_closed`, `1 in W_k iff 1 in S`) is proved as
  written, for any subspace S and all k0, D, m. The "Consequently" part is proved given
  `2 != 0`, which the paper omits. The paper's statement is weaker than the note's (no
  containment) while the paper's later text uses the note's stronger form (F2). The two
  instance bounds `lfd >= 5, k >= 8` follow, subject to the scope caveat F1 and to the
  uncompiled kernel certificates; the base facts and a witness were reproduced independently.
- **`thm:restrict`**: proved as written (proof is in `U1-upper-bounds.md` and Lean, not in
  N5/N7). The monotonicity sentence is a correct corollary (F8).
- **N5 Cor. 7.3 / proper ideals by dimension**: correct under `2^k | p-1` and a nonsingular
  curve, which must be stated; written-only, labelled as such; not an all-depth result (F5).
- **N7 (`lfd_gt_of_lowIdeal`, `growth_schema`, the two instance theorems)**: argument
  correct as written; the statements match the Lean text; scope caveat F1 applies.
