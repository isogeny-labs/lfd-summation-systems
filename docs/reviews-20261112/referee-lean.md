# Static audit of `formal/ICFallDegrees` against `paper/curve-lfd-bounds.tex`

Repository: `/home/user/lfd-summation-systems` at commit `705279f`. Paths below are relative to it.

## Verdict

No escape hatch was found, every Lean declaration the paper names exists, and every kernel-checked theorem I read states what the paper's prose states (several are more general).
Nothing was compiled: `lake` and `lean` are not on PATH and there is no `formal/.lake`, so the axiom evidence is the committed log, not something I reproduced.
The findings are about interpretation and disclosure, not soundness.

## Findings, ranked

### Medium

**M1. Nothing was built here; the committed build log is not from-scratch evidence.**
- `experiments/ic-lean/build.log` contains only `Replayed ICFallDegrees.X` lines (19 modules) and `Build completed successfully`. "Replayed" is what Lake prints for an already-built module, so the log records a rebuild on a warm tree.
- The kernel-heavy modules (`N5CertKernel*`, `ExactK6Kernel*`) do not appear in it at all, which is consistent with a quiet warm build but shows nothing about them.
- What I could check statically: all 110 `source_sha256` entries in `experiments/ic-lean/verification.json` match the current files, and every file under `formal/ICFallDegrees/` is hashed.
- `axioms.log` and `verification.json` agree with each other (513 declarations, same names).

**M2. The `p = 11` instance of Theorem `thm:reduction` is not an elliptic-curve subgroup system at any depth, and the `p = 40961` instance is one only for `k ≤ 13`.**
- The Lean statements are faithful to the literal polynomial system: `N5CertInstanceP11.not_one_in_closure` and `N5CertWitnessInstanceP11.lfd_ge_five` quantify over all `k ≥ 8` with `(a,b,t) = (2,5,7)` over `ZMod 11`.
- For `p = 11`: `p − 1 = 10`, so `2^k ∤ p − 1` for every `k ≥ 2`, and there is no subgroup of order `2^k` in `F_11^*`.
- For `p = 11`: `t^3 + at + b = 10` is a non-residue mod 11 (my computation), so `t = 7` is not the abscissa of an `F_11`-rational point of `y^2 = x^3 + 2x + 5`.
- For `p = 40961`: `v_2(p − 1) = 13`, so "`lfd ≥ 5` at every depth `k ≥ 8`" is a subgroup statement only for `8 ≤ k ≤ 13`. The target is a rational abscissa (rhs is a residue).
- The paper calls `p = 11` "a small analogue" and gives the extension-field caveat for `p = 12289` and `p = 65537` under Theorem `thm:restrict`, but not for these two instances.
- For `p = 1153` (Theorem `thm:exactsix`): `64 | 1152`, the rhs is a residue and the discriminant is non-zero, so that instance is a genuine subgroup system at `k = 6`.

**M3. There is no curve in Lean; every theorem is about an explicitly written polynomial system.**
- No group law, no proof that `S3`/`S4` are summation polynomials, no hypothesis that `t` is an abscissa, no `2^k | p − 1`, no nonsingularity. `formal/README.md:131-133` says so ("The group law is not formalized here"), and the paper's introduction says the statements are "about specific polynomial encodings".
- `S3` is Semaev's polynomial verbatim. `M3CarryLower.lean:46`:
  `def S3 (a b u v w) := (u-v)^2*w^2-2*((u+v)*(u*v+a)+2*b)*w+(u*v-a)^2-4*b*(u+v)`
- The two-summand generator is written directly in normalized form (`PrimeM2General.normalized`, `PrimeM2General.lean:26`). I expanded `S3(u,w,t)` by hand with `u^2 → u_1`, `w^2 → w_1` and got the same top `u_1w_1 − 2(t^2+a)u_0w_0 − 2t(u_0w_1+u_1w_0)` and tail `t^2(u_1+w_1) + (−2at−4b)(u_0+w_0) + a^2−4bt`. Lean also proves `raw_normalization_identity` (raw − normalized is an explicit combination of two chain equations).
- `S4` is defined by the closed two-quadratic resultant formula (`M3CarryLower.lean:55`), `(AC'−A'C)^2 − (AB'−A'B)(BC'−B'C)`. It is not connected to Mathlib's `Polynomial.resultant`. The formula is the standard one; I could not machine-check it (no sympy). The `(3,3,3)` coefficient `−8(5t^3+3at+6b)` is kernel-proved (`coeff_s4_333`).
- For `m ≥ 4` no `S_{m+1}` exists in Lean. Corollary `cor:allm` rests on `ChainQuotient.reduced_inputs_lfd_bound`, which takes any set of reduced inputs; the paper says this.

### Low

**L1. The binary "normalized" constraint in Lean is only partly reduced.**
- `BinaryS3Bounds.normalized x y z xb yb zb b = xb*yb+xb*zb+yb*zb+x*y*z+C b` (`BinaryS3Bounds.lean`), where `xb` is the Frobenius-linearised square of the affine form `x`. Products are not reduced modulo `Z^2 − Z`, and the degree bound used is 3, not 2.
- This equals the paper's fully reduced coordinate polynomial (Section 2) only when the three affine forms have disjoint variable supports, which the theorem neither requires nor proves.
- The bound is unaffected: it holds for any input list containing the Boolean equations with degrees `≤ r+1`.
- `koblitz_descent_lfd_bound` is literally `binary_descent_lfd_bound` with `b := 1`. The coefficient `a` and Frobenius do not occur. The paper says the entry "is identical in content".
- The binary `raw` is `(xy)^2+(xz)^2+(yz)^2+xyz+b`, which is Semaev's binary `S3` in characteristic 2.

**L2. Hypotheses or forms that differ slightly from the prose.**
- Theorem `thm:carry`, general part: Lean also needs `(0 : P r K) ∉ G` (`CarryLower.not_hasLFDBound_of_carry`). The paper states `0 ∉ G` only in the finite-quotient paragraph.
- `GrowthSchema.growth_schema` requires `2 ≤ k0`; the paper says "any degree `D` and base depth `k_0`".
- Theorem `thm:restrict`: "for unit ideals the last fall degree is nondecreasing in the depth" is not a Lean declaration. Lean has `lfd_gt_of_le` for `D ≥ 2`; the monotonicity follows with a short extra argument (all inputs have degree 2).
- Theorem `thm:bridge`: Lean proves the iff (`chain_mem_closure_iff`). The direct-sum form "`V_D = (⟨G⟩ ∩ M_{≤D}) ⊕ V^q_D`" is not a separate declaration that I found.
- In the other direction, Lean is more general than the paper in several places: any field with `2 ≠ 0` instead of `F_p`, no `2^k | p − 1`, arbitrary (not necessarily finite) input sets.

**L3. One cited name is a module path, not the declaration's namespace.**
- The paper cites `QuotientLowerBound.closure_le_quotient_span`. The declaration is `ICFallDegrees.FiniteQuotient.closure_le_quotient_span` (file `QuotientLowerBound.lean:73`, audited at `axioms.log:74`).
- `eq_zero_of_eval_zero_at_prod_finset` is a Mathlib lemma, correctly described as such, and not in the project audit.

**L4. Stale audit README.**
- `experiments/ic-lean/README.md:6` says "339 audited declarations"; `verification.json`, `axioms.log`, `formal/README.md:45` and the paper (`\LEANDECLS` = 513) say 513.
- The same README (line 31) still lists the prime chain quotient instantiation as an open obligation; `verification.json` "limitations" says it is kernel-checked.

**L5. Trust base of the certificate theorems.**
- About 290 `decide +kernel` proofs (`N5CertKernelP11/P40961`, `N5CertWitnessInstance*`, `ExactK6KernelLower/UpperP1153`). These add no axiom but rely on the kernel's built-in big-number `Nat` arithmetic.
- The checker soundness theorems are ordinary Lean proofs, and the final statements mention only `PrimeM2General.inputs` with numeric constants, so the certificate data cannot weaken the statement.
- The verify script's forbidden-token scan (`scripts/verify_lean_fall_degrees.py:53`) omits `opaque`, `partial` and `csimp`; my grep found none of them.

## Claim-by-claim table

All line numbers are in `formal/ICFallDegrees/`. "Faithful" means I read the Lean statement and the definitions it uses.

| Paper claim | Lean declaration | Verdict |
|---|---|---|
| Section 2: `V_{F,D}`, `lfd(F) ≤ d` | `InClosure` (`Closure.lean:29`), `HasLFDBound` (`:118`), `IsLastFallDegree` (`:121`) | Faithful. Budget is `a.totalDegree + f.totalDegree ≤ D`, applied to derived rows; ordinary ring `MvPolynomial σ K`, not a quotient. |
| Finite-quotient theorem | `FiniteQuotient.finite_quotient_lfd_of_reduced_zero` | Faithful (`G ⊆ F`, `0 ∉ G`, reduced-ideal trivial, standard degree `≤ s`, input degree `≤ δ`, conclusion `max (s+1) δ`). |
| Boolean instance | `GridQuotient.boolean_lfd_bound` (`:95`) | Faithful. |
| Lemma `lem:leading` | `leading_chainPolynomial`, `standard_le_one`, `standard_degree_bound` (`ChainQuotient.lean`) | Faithful. |
| Theorem `thm:reduced` | `ChainQuotient.reduced_chain_ideal_trivial` (`:244`) | Faithful (`2 ≠ 0`, all `m k`). |
| Theorem `thm:chain` | `ChainQuotient.chain_lfd_bound`, `CurveBounds.prime_chain_lfd_bound` | Faithful; Lean does not need `m,k ≥ 1`. |
| Corollary `cor:m2` | `ChainQuotient.prime_m2_lfd_bound`, `prime_field_m2_lfd_bound`, `CurveBounds.prime_m2_lfd_bound`; `chainInputs_two_subset`; `PrimeM2General.degree_and_top` | Faithful (`5 ≤ p`, `2 ≤ k`, all `a b t`). Generator checked against Semaev `S3`, see M3. |
| Theorem `thm:three` | `PrimeM2Lower.no_one_at_three`, `lfd_ge_four_of_unit` | Faithful (`2 ≠ 0`, `2 ≤ k`, `3t^2+a ≠ 0`; unit ideal is an explicit hypothesis of the second). |
| FFD remark after `thm:three` | `PrimeM2General.exact_ordinary_ffd`, `PrimeFieldConclusions.ordinary_ffd` | Faithful: Koszul first fall `if 3*t^2+a=0 then 3 else 4`. |
| Corollary `cor:allm` | `ChainQuotient.reduced_inputs_lfd_bound` | Faithful as "chains plus any reduced inputs"; no `S_{m+1}` object for `m ≥ 4` (M3). |
| Table `tab:bounds`, binary and Koblitz | `CurveBounds.binary_descent_lfd_bound`, `koblitz_descent_lfd_bound` | Faithful with the caveats of L1 (`2 ≤ r`). |
| Table `tab:bounds`, Boolean sharpness | `BooleanSharp.exact_lfd` (`:128`) | Faithful: inputs are the field equations and `1 − ∏(1 − X_i)`, `IsLastFallDegree … (r+1)`, `2 ≤ r`. |
| Table `tab:twosided`, `= 4` family | `PrimeFieldConclusions.sharp_lfd`, `PrimeM2Sharp.exact_lfd_four`, `one_at_four` | Faithful (`b = t = 0`, `a ≠ 0`, `a^(2^k) ≠ 1`). |
| Theorem `thm:carry`, general part | `CarryLower.not_hasLFDBound_of_carry`, `chain_carry_lower` | Faithful with unstated `0 ∉ G` (L2). |
| Theorem `thm:carry`, three summands | `M3CarryLower.m3_lfd_ge_seven`, `prime_field_m3_lfd_ge_seven` (`:526`, `:550`) | Faithful: `∀ e≤6, ¬HasLFDBound (chainInputs K 3 k ∪ {chainS4 k a b t}) e`, `3 ≤ k`, `5t^3+3at+6b ≠ 0`. `chainS4_normal_form` proves the generator is reduced and congruent to `S4(u_{1,0},u_{2,0},u_{3,0},t)`. `S4` is the resultant closed form (M3). |
| Theorem `thm:bridge` | `ChainClosureBridge.chain_mem_closure_iff` (`:404`) | Faithful for the iff; direct-sum wording not a declaration (L2). |
| Theorem `thm:reduction`, general | `DepthReduction.W_closed` (`:517`), `single_zero_mem_W_iff` (`:391`), `ChainModelMap.n5_lfd_ge_five`, `prime_field_n5_lfd_ge_five` | Faithful. Base closedness, `1 ∉ S` and generator-in-base are explicit hypotheses. |
| Two explicit instances, `1 ∉ V_4(k)`, `k ≥ 8` | `N5CertInstanceP40961.not_one_in_closure`, `N5CertInstanceP11.not_one_in_closure` | Faithful. Constants match: `(21223, 9887, 4748)` over `ZMod 40961`, `(2, 5, 7)` over `ZMod 11`. See M2. |
| Same instances, `lfd ≥ 5` with no ideal hypothesis | `N5CertWitnessInstanceP40961.lfd_ge_five`, `N5CertWitnessInstanceP11.lfd_ge_five` | Faithful: `∀ e < 5, ¬HasLFDBound …`, no unit hypothesis. See M2. |
| Growth schema | `GrowthSchema.growth_schema` (`:142`) | Faithful with unstated `2 ≤ k0` (L2). |
| Theorem `thm:exactsix` | `ExactK6InstanceP1153.exact_lfd_five` | Faithful: `IsLastFallDegree (Set.range (inputs 6 hk (139 : ZMod 1153) 1094 185)) 5`; both halves proved. |
| Theorem `thm:restrict` | `DepthRestrictionM2.restrict_closure`, `not_one_of_le`, `lfd_gt_of_le` | Faithful (any field, `2 ≤ k`, `2 ≤ D`); the monotonicity sentence is a corollary, not a declaration (L2). |
| Consequence of `thm:restrict` for `p = 1153` | `ExactK6Deeper.not_one_in_closure_four_of_ge` | Faithful (`6 ≤ k`). |
| Selector bound tool (written theorem) | `LocalReductionBound.lfd_bound_of_local` | Exists and is a genuine local-to-global criterion; the paper correctly says the selector theorem itself is not formalized. |
| Not found | none | Every `\code{…}` Lean identifier in the paper resolves to a declaration or module. |

## Escape-hatch grep

Searched all `.lean` files under `formal/` for `sorry`, `admit`, `axiom`, `native_decide`, `unsafe`, `implemented_by`, `opaque`, `extern`, `ofReduceBool`, `reduceBool`, `csimp`, `partial`, `trustCompiler`, `skipKernel`, `macro`, `elab`, `syntax`, `run_cmd`, `#eval`.

| Hit | Matters? |
|---|---|
| `SparsePrimeCertificate.lean:10` "without native_decide" | No, a comment. |
| `QuadraticNoLow.lean:82` "admit no nonzero" | No, a docstring. |
| `N5CertWitness.lean:19,83`, `N5CertWitnessInstanceP11.lean:40`, `…P40961.lean:41`, `ChainQuotient.lean:22` "partial" | No, prose. |
| `Certificate.lean:7`, `SparseCertificate.lean:8`, `SparseDAG.lean:4`, `SparseLowerCertificate.lean:32` "external" | No, prose. |
| `SparseDemo.lean:7` `#eval` | No. The file is not imported by anything; the verify script runs it as a runtime diagnostic. |
| `set_option maxHeartbeats` (7, `M3CarryLower.lean`), `set_option Elab.async false` (kernel-certificate files) | No. |
| `decide +kernel` (about 290 uses) | Kernel evaluation, no axiom; see L5. |

No `sorry`, `admit`, `axiom`, `native_decide`, `unsafe`, `implemented_by`, `opaque` or `extern` as code anywhere.

## Import coverage

- `formal/ICFallDegrees.lean` imports 99 modules.
- Not imported by the root but reached transitively through `ExactK6InstanceP1153`: `ExactK6Cert`, `ExactK6DataP1153`, `ExactK6KernelLowerP1153`, `ExactK6KernelUpperP1153`, `ExactK6LowerSound`, `ExactK6UpperSound`.
- Imported by nothing: `Audit.lean` (513 `#print axioms` lines, run separately by the script) and `SparseDemo.lean`. Neither contains a theorem the paper cites.

## Audit records

- `experiments/ic-lean/axioms.log`: 513 entries; 506 `[propext, Classical.choice, Quot.sound]`, 5 `[propext, Quot.sound]`, 2 `[propext]`. No `sorryAx`, no `Lean.ofReduceBool`.
- Every paper-cited project declaration is in the log, including all of the table above, with only those three axioms.
- Toolchain `leanprover/lean4:v4.19.0`; `lake-manifest.json` pins Mathlib at `c44e0c8ee63ca166450922a373c7409c5d26b00b` (input rev `v4.19.0`), matching `verification.json`.

## What was and was not verified

- Verified by reading: every statement in the table and the definitions it depends on; the grep results; import coverage; name and axiom-log coverage; source hashes; the three instances' arithmetic side facts (Python).
- Checked by hand, not by machine: the expansion of `S3` into the normalized generator; the two-quadratic resultant formula.
- Not verified: that the project compiles, that `axioms.log` is the output of the current sources, any proof body beyond spot reading, and the running times the paper quotes.
