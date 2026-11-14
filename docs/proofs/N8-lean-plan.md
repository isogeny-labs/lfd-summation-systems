# N8 in Lean: the growth theorem, kernel-checked conditional on one analytic hypothesis

Status (2026-11-12). Evidence labels as in the repository: **kernel**, **proof**, **exact**,
**empirical**, **conjecture**. This note describes the Lean development
`formal/ICFallDegrees/Growth{Counting,Reduction,Spread,Main,Criterion,Existence}.lean` for the
written proof `docs/proofs/N8-growth.md` (T1, C1). The files are built one module at a time; they
are imported by `ICFallDegrees.lean` and listed in `Audit.lean` (35 lines; audit of 582 declarations, 2026-11-12).

## 1. What is kernel-checked, and under which hypothesis

**Result (kernel, conditional).** Let `K` be a field with `2 ≠ 0`, `k ≥ 2`, `a b t : K`. Suppose
a `CurveModel F₀ a b t p ξ` is given (Section 2) with `p = 0` or `p ≥ 2·2^k + 1`. Then for
`D ≥ 1` and `k ≥ D + Λ(D−1)`:

    ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D 1        (GrowthMain.growth_theorem)

and for a unit ideal `∀ e ≤ D, ¬ HasLFDBound … e` (`growth_lfd_gt`). The same holds with the
explicit threshold `k ≥ 14 D + 10` (`growth_theorem_linear`), so for `k ≥ 24` no
`e ≤ ⌊(k−10)/14⌋` is a bound (`growth_lfd_gt_linear`). With the kernel upper bound
`SelectorExtension.lfd_le_depth_add_two_of_two_ne`, every last fall degree `d` satisfies
`D < d ≤ k + 2` (`growth_lfd_window`).

All statements are about the audited definitions (`InClosure`, `HasLFDBound`,
`IsLastFallDegree`, `PrimeM2General.inputs`); the model spaces are linked to them through the
existing bridge `ChainClosureBridge.prime_m2_not_one`. Axioms of every theorem:
`propext, Classical.choice, Quot.sound`.

**The hypothesis is not proved in Lean.** The label of T1 and C1 becomes
"**kernel, conditional on `CurveModel`** (which is **proof**: N8-growth.md Sections 1 and 3)".

## 2. The analytic hypothesis, exactly

```lean
def mono (ξ : Fin 2 → F) (e : Fin 2 →₀ ℕ) : F := ξ 0 ^ e 0 * ξ 1 ^ e 1

def Circuit (F₀) (v : ι → F) (J : Finset ι) : Prop :=
  ¬ LinearIndepOn F₀ v (J : Set ι) ∧ ∀ J' ⊂ J, LinearIndepOn F₀ v (J' : Set ι)

def SpreadBound (F₀) (ξ : Fin 2 → F) (p : ℕ) : Prop :=
  ∀ J : Finset (Fin 2 →₀ ℕ), Circuit F₀ (mono ξ) J → ∀ e ∈ J, ∀ e' ∈ J, ∃ z₁ z₂ : ℤ,
    |2 * ((e 0 : ℤ) - e' 0) - p * z₁| + |2 * ((e 1 : ℤ) - e' 1) - p * z₂| ≤
      6 * ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2)

structure CurveModel (F₀) (a b t : K) (p : ℕ) (ξ : Fin 2 → F) : Prop where
  ne_zero : ∀ i, ξ i ≠ 0
  on_curve : S3at a b t ξ = 0          -- S_3(ξ₁, ξ₂, t) = 0 in F
  spread : SpreadBound F₀ ξ p          -- Theorem 3.3 of N8-growth.md
```

`K → F₀ → F` is a tower of fields (`IsScalarTower K F₀ F`). `SpreadBound` is Theorem 3.3 word
for word: "circuit" is a minimal `F₀`-dependent finite family of monomials `ξ^e`, indexed by
exponents (so two different exponents with equal monomials form a circuit of two terms, as in
the written proof), and `‖n‖_p ≤ B` is written `∃ z, |n − p z| ≤ B`.

What the written proof must supply, and Lean does not check:

1. the function field `F` of `E`, the functions `ξ₁ = x(P)`, `ξ₂ = x(R − P)` with `x(R) = t`, and
   the Semaev identity `S_3(ξ₁, ξ₂, t) = 0` (Section 1; gives `ne_zero`, `on_curve`);
2. Theorem 3.3 for `F₀` = constants (characteristic 0) or `F^p` (characteristic `p`), `p` the
   characteristic (Sections 3.1–3.2; this is where **nonsingular** and **t ≠ 0** are used).

Consequently "E nonsingular" and "t ≠ 0" do **not** appear as hypotheses of the Lean theorems:
they are the conditions under which the written proof produces a `CurveModel`. The parameter `p`
is not tied to `CharP K p` in Lean; only `p = 0 ∨ 2·2^k + 1 ≤ p` is used.

## 3. Map from the written proof to Lean

| written (N8-growth.md) | Lean | label |
|---|---|---|
| lift `f̃`, top-run part `f̃_hi`, FEI_k(w), S_D(k) | `GrowthReduction.{expand …, topPart, FEI, idealSpace}`; the lift is `ChainQuotient.expand K 2 k` on reduced chain polynomials | definitions |
| identity (1.1) | `GrowthReduction.lift_identity` (multiplied by `x_i^(N−2^j)`) | kernel |
| Prop 2.1 (FEI ⇒ closed) | `GrowthReduction.idealSpace_closed` | kernel |
| Prop 2.1 (equivalence) | `GrowthCriterion.closed_iff_fei` | kernel |
| Prop 2.2 | `GrowthReduction.not_one_of_fei`, `closure_normalForm_mem_idealSpace`, `lfd_gt_of_fei` | kernel |
| `g̃ = S_3` | `GrowthReduction.expand_normalized` | kernel |
| Theorem 3.1, Lemma 3.2, Theorem 3.3 | **hypothesis** `GrowthSpread.SpreadBound` | proof (not in Lean) |
| Lemma 3.4 | `GrowthSpread.split_of_no_crossing` (any family in any `F₀`-vector space) | kernel |
| Λ(w), Lemma 4.1 | `GrowthCounting.{Lam, Lam_spec, Lam_greatest, lambdaBound_Lam, Lam_le, lambdaBound_linear}`; values `GrowthMain.Lam_one … Lam_five` = 25, 36, 48, 59, 71 | kernel |
| Lemma 4.2 | `GrowthCounting.{window_cases, window_index}` | kernel |
| Lemma 4.3 | `GrowthCounting.{card_le_of_window, sigma_lt_of_lambda, dichotomy}`, `GrowthSpread.circuit_sigma_lt` | kernel (from `SpreadBound`) |
| Lemma 4.4 | `GrowthCounting.{wt_ge_of_close, wt_ge_of_top, straddle_contra}` | kernel |
| Theorem 5.1 | `GrowthSpread.fei_of_spread`; `GrowthMain.{growth_not_one, growth_theorem, growth_theorem_linear, growth_closure_vanishes, growth_theorem_small}` | kernel (from `CurveModel`) |
| Corollary 5.2 (a), (b) | `GrowthMain.{growth_lfd_gt, growth_lfd_gt_linear, growth_lfd_window, growth_lfd_gt_zmod}` | kernel (from `CurveModel`) |
| Corollary 5.2 (c), existence of unit instances | `GrowthExistence.{fermat_factor_dvd, unit_instance}` (unconditional); `growth_instance`, `growth_instance_window` | kernel; the last two from `CurveModel` |
| Prop 7.1 | `GrowthCriterion.{sharpSpace, Crit, not_one_of_crit, closure_normalForm_mem_sharpSpace, crit_of_fei}` | kernel |
| Corollary 5.3 (free extension of S_D), Prop 7.2, Sections 6–7 data | not formalized | proof / exact |

## 4. Differences from the written proof

1. **Lemma 4.1 is proved by another argument.** The written proof shows that
   `ρ(ℓ) = 2^(ℓ−1)/(7776·C(2ℓ, ≤w)²)` is increasing for `ℓ ≥ 4w` (the step that review finding F1
   corrected) and bounds `Λ(w) ≤ 13w + 23` with `C(n, ≤w) ≤ (en/w)^w` and real logarithms. Lean
   does not use the ratio argument. It proves `C(n, ≤w)·16^n ≤ 17^n·16^w` by Pascal induction
   (`Csum_chernoff`) and then `7776·C(2ℓ, ≤w)² ≤ 2^(ℓ−1)` for all `ℓ ≥ 13w + 24` by integer
   induction (`numeric_base`, `numeric_all`, `lambdaBound_linear`); the numeric steps are
   `17^52·16^2 ≤ 2^13·16^52` and `17^4 ≤ 2·16^4`. `Λ(w)` is then *defined* as the largest
   solution below `13w + 23` and proved to be the largest solution overall. The statement
   `Λ(w) ≥ 4w` of Lemma 4.1 is not proved (it is not used by Theorem 5.1). The bound holds for
   `w = 0` too.
2. **Window count.** `blocks` has five block indices per coordinate (the written lemma says at
   most six); the constant `36` is kept, so `Λ` is the written one.
3. **`D ≥ 1` instead of `D ≥ 2`.** The Lean theorems hold for `D = 1` (trivially interesting).
4. **Weight = degree.** The binary weight of an exponent pair is the total degree of the
   corresponding squarefree chain monomial (`ChainModelMap.digits_degree`); `S_D(k)` is a
   subspace of the reduced chain polynomials, so no separate group-algebra model is introduced.
5. No error was found in the written proof. The red-team finding F1 concerns a step that the
   Lean proof bypasses.

## 5. Not covered

* Theorem 3.1, Lemma 3.2, Theorem 3.3 and the function-field setting of Section 1 (the
  hypothesis `CurveModel`). Needs function fields of curves, places, differentials and
  Wronskians, absent from the pinned Mathlib.
* Corollary 5.3 (T4), Proposition 7.2, the statement `Λ(w) ≥ 4w`, the table of `Λ` beyond `w = 5`,
  the asymptotic slope, and all computations of Sections 6–7.
* The unequal-depth and coset remarks of Section 5.

Next step, if wanted: Corollary 5.3 reuses `split_of_no_crossing`, `circuit_sigma_lt` and the
digit lemmas with the cut at bit `k` of depth `k + 1`; it needs one more digit lemma and the
division by `x_i^N` (prime ideal), about 200 lines.

## 6. Build and integration

Each module was built alone through the lock
(`flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake build ICFallDegrees.<Module>'`);
all six build without errors or warnings from these files. `GrowthCounting` deliberately does
not import `Mathlib.Tactic` (with it the build exceeds the 7 GB address-space limit).
`GrowthMain.Lam_one … Lam_five` use `decide +kernel` on integers below `2^100` (no
`native_decide`).

Lines for `formal/ICFallDegrees.lean`:

```lean
import ICFallDegrees.GrowthCounting
import ICFallDegrees.GrowthReduction
import ICFallDegrees.GrowthSpread
import ICFallDegrees.GrowthMain
import ICFallDegrees.GrowthCriterion
import ICFallDegrees.GrowthExistence
```

Lines for `formal/ICFallDegrees/Audit.lean`:

```lean
#print axioms ICFallDegrees.GrowthCounting.lambdaBound_linear
#print axioms ICFallDegrees.GrowthCounting.lambdaBound_Lam
#print axioms ICFallDegrees.GrowthCounting.Lam_greatest
#print axioms ICFallDegrees.GrowthCounting.card_le_of_window
#print axioms ICFallDegrees.GrowthCounting.straddle_contra
#print axioms ICFallDegrees.GrowthReduction.lift_identity
#print axioms ICFallDegrees.GrowthReduction.idealSpace_closed
#print axioms ICFallDegrees.GrowthReduction.not_one_of_fei
#print axioms ICFallDegrees.GrowthReduction.closure_normalForm_mem_idealSpace
#print axioms ICFallDegrees.GrowthReduction.lfd_gt_of_fei
#print axioms ICFallDegrees.GrowthSpread.split_of_no_crossing
#print axioms ICFallDegrees.GrowthSpread.circuit_sigma_lt
#print axioms ICFallDegrees.GrowthSpread.fei_of_spread
#print axioms ICFallDegrees.GrowthMain.growth_not_one
#print axioms ICFallDegrees.GrowthMain.growth_theorem
#print axioms ICFallDegrees.GrowthMain.growth_theorem_linear
#print axioms ICFallDegrees.GrowthMain.growth_closure_vanishes
#print axioms ICFallDegrees.GrowthMain.growth_lfd_gt
#print axioms ICFallDegrees.GrowthMain.growth_lfd_gt_linear
#print axioms ICFallDegrees.GrowthMain.growth_lfd_window
#print axioms ICFallDegrees.GrowthMain.growth_lfd_gt_zmod
#print axioms ICFallDegrees.GrowthMain.Lam_one
#print axioms ICFallDegrees.GrowthMain.Lam_two
#print axioms ICFallDegrees.GrowthMain.Lam_three
#print axioms ICFallDegrees.GrowthMain.Lam_four
#print axioms ICFallDegrees.GrowthMain.Lam_five
#print axioms ICFallDegrees.GrowthMain.growth_theorem_small
#print axioms ICFallDegrees.GrowthCriterion.closed_iff_fei
#print axioms ICFallDegrees.GrowthCriterion.not_one_of_crit
#print axioms ICFallDegrees.GrowthCriterion.closure_normalForm_mem_sharpSpace
#print axioms ICFallDegrees.GrowthCriterion.crit_of_fei
#print axioms ICFallDegrees.GrowthExistence.fermat_factor_dvd
#print axioms ICFallDegrees.GrowthExistence.unit_instance
#print axioms ICFallDegrees.GrowthExistence.growth_instance
#print axioms ICFallDegrees.GrowthExistence.growth_instance_window
```

Each of these prints `[propext, Classical.choice, Quot.sound]` (checked from a scratch file with
`lake env lean`, under the lock, 2026-11-12; `Lam_two … Lam_four` are covered by
`growth_theorem_small`, which depends on them).
