import ICFallDegrees.GrowthSpread
import ICFallDegrees.UnitInstance
import ICFallDegrees.SelectorExtension

noncomputable section

/-!
# The growth theorem, conditional on the spread of circuits
(Theorem 5.1 / T1 and Corollary 5.2 / C1 of `docs/proofs/N8-growth.md`)

**Status: kernel, conditional on one named analytic hypothesis.**  Every theorem here is about
the audited objects `InClosure`, `HasLFDBound`, `IsLastFallDegree` and
`PrimeM2General.inputs k depth a b t`.  The hypothesis is the structure `CurveModel`: a field
`F ⊇ F₀ ⊇ K` with two nonzero elements `ξ₁, ξ₂` satisfying `S_3(ξ₁, ξ₂, t) = 0` and the
circuit spread bound `GrowthSpread.SpreadBound` (Theorem 3.3 of the written proof).  The written
proof supplies such a model for every nonsingular curve `y² = x³ + a x + b` and every `t ≠ 0`
(function field of the curve, `ξ₁ = x(P)`, `ξ₂ = x(R - P)`, `x(R) = t`, `p = char K`); that
construction and Theorem 3.3 are **not** formalized.  Consequently the hypotheses "nonsingular"
and "`t ≠ 0`" of the written statement do not appear below: they are what the written proof
uses to produce a `CurveModel`.

Main statements:

* `growth_not_one`: `1 ∉` the degree-`D` closure for `k ≥ D + lam`, any `lam ≥ 1` with
  `LambdaBound (D-1) lam`;
* `growth_theorem`: the same with `lam = Λ(D-1)` (`GrowthCounting.Lam`), i.e. `k ≥ k_D`;
* `growth_theorem_linear`: the same for `k ≥ 14 D + 10`;
* `growth_closure_vanishes`: every element of the closure has its lift in the ideal of the model;
* `growth_lfd_gt`, `growth_lfd_gt_linear`, `growth_lfd_window`: last-fall-degree consequences
  for unit ideals (with the selector upper bound `k + 2`);
* `growth_lfd_gt_zmod`: the prime-field instance with the unit ideal of `UnitInstance.unit_zmod`;
* `Lam_one … Lam_five`: `Λ(1), …, Λ(5) = 25, 36, 48, 59, 71` (so `k_2, …, k_6 = 27, 39, 52, 64, 77`).
-/

namespace ICFallDegrees.GrowthMain

open MvPolynomial GrowthCounting GrowthReduction GrowthSpread
open PrimeM2General (inputs)

variable {K F₀ F : Type*} [Field K] [Field F₀] [Field F] [Algebra K F₀] [Algebra F₀ F]
  [Algebra K F] [IsScalarTower K F₀ F]

/-- **The analytic input of the growth theorem (not proved in Lean).**

A model of the curve `S_3(x₁, x₂, t) = 0` in a field `F` over `K`, with an intermediate field
`F₀`: two nonzero elements `ξ₁, ξ₂ ∈ F` on the curve, such that the monomials `ξ₁^a ξ₂^b`
satisfy the circuit spread bound of Theorem 3.3 of `docs/proofs/N8-growth.md` over `F₀` with
the parameter `p` (the characteristic, in the written proof). -/
structure CurveModel (F₀ : Type*) {K F : Type*} [Field K] [Field F₀] [Field F] [Algebra F₀ F]
    [Algebra K F] (a b t : K) (p : ℕ) (ξ : Fin 2 → F) : Prop where
  /-- The coordinate functions are nonzero. -/
  ne_zero : ∀ i, ξ i ≠ 0
  /-- `S_3(ξ₁, ξ₂, t) = 0`. -/
  on_curve : S3at a b t ξ = 0
  /-- Theorem 3.3 of the written proof. -/
  spread : SpreadBound F₀ ξ p

section Model

variable {a b t : K} {p : ℕ} {ξ : Fin 2 → F}

omit [Field F₀] [Algebra K F₀] [Algebra F₀ F] [IsScalarTower K F₀ F] in
theorem X_not_mem_ker (hξ : ∀ i, ξ i ≠ 0) (i : Fin 2) :
    (X i : Q K) ∉ RingHom.ker (aeval ξ : Q K →ₐ[K] F) := by
  rw [RingHom.mem_ker, aeval_X]
  exact hξ i

omit [Field F₀] [Algebra K F₀] [Algebra F₀ F] [IsScalarTower K F₀ F] in
theorem S3poly_mem_ker (h : S3at a b t ξ = 0) :
    S3poly a b t ∈ RingHom.ker (aeval ξ : Q K →ₐ[K] F) := by
  rw [RingHom.mem_ker, aeval_S3poly]
  exact h

omit [Field F₀] [Algebra K F₀] [Algebra F₀ F] [IsScalarTower K F₀ F] in
instance ker_isPrime (ξ : Fin 2 → F) : (RingHom.ker (aeval ξ : Q K →ₐ[K] F)).IsPrime :=
  RingHom.ker_isPrime _

end Model

/-! ## The growth theorem -/

/-- **Growth theorem, general threshold.** Let `K` be a field with `2 ≠ 0`, `(F, F₀, ξ)` a
`CurveModel` for `(a, b, t)` with parameter `p`, `p = 0` or `p > 2^(k+1)`, `D ≥ 1`, and
`lam ≥ 1` such that `2^(ℓ-1) < 7776·C(2ℓ, ≤ D-1)^2` fails for all `ℓ > lam`.  If `k ≥ D + lam`,
then `1` is not in the degree-`D` closure of the two-summand chain presentation. -/
theorem growth_not_one (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D lam : ℕ) (hD : 1 ≤ D)
    (hlam : LambdaBound (D - 1) lam) (h1 : 1 ≤ lam) (hk : D + lam ≤ k) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : MvPolynomial (Fin (2 * k)) K) :=
  not_one_of_fei two_ne k depth a b t (RingHom.ker (aeval ξ : Q K →ₐ[K] F))
    (X_not_mem_ker model.ne_zero) (S3poly_mem_ker model.on_curve) D
    (fei_of_spread (F₀ := F₀) ξ p k (D - 1) lam hp model.spread hlam h1 (by omega))

/-- **Growth theorem (T1 of `docs/proofs/N8-growth.md`), conditional on `CurveModel`.**
For `D ≥ 1` and `k ≥ k_D = D + Λ(D-1)`: `1 ∉ V_D(k)`. -/
theorem growth_theorem (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : MvPolynomial (Fin (2 * k)) K) :=
  growth_not_one two_ne k depth a b t ξ p model hp D (Lam (D - 1)) hD (lambdaBound_Lam _)
    (Lam_pos _) hk

/-- The growth theorem with the explicit linear threshold `k ≥ 14 D + 10`. -/
theorem growth_theorem_linear (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ) (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : MvPolynomial (Fin (2 * k)) K) :=
  growth_not_one two_ne k depth a b t ξ p model hp D (13 * (D - 1) + 23) hD
    (lambdaBound_linear _) (by omega) (by omega)

/-- `k_D ≤ 14 D + 10`. -/
theorem kD_le (D : ℕ) (hD : 1 ≤ D) : D + Lam (D - 1) ≤ 14 * D + 10 := by
  have := Lam_le (D - 1)
  omega

/-- **T1, second sentence.** Under the hypotheses of the growth theorem the chain normal form of
every element of the degree-`D` closure is reduced, of degree `≤ D`, and its lift (chain
variables expanded as `u_{i,j} ↦ x_i^(2^j)`) vanishes at `(ξ₁, ξ₂)`. -/
theorem growth_closure_vanishes (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k)
    {f : MvPolynomial (Fin (2 * k)) K} (hf : InClosure (Set.range (inputs k depth a b t)) D f) :
    aeval ξ (ChainQuotient.expand K 2 k
      (ChainClosureBridge.normalForm (ChainQuotient.chainInputs K 2 k)
        (ChainQuotient.chainInputs_nonzero 2 k) f)) = 0 := by
  have h := closure_normalForm_mem_idealSpace two_ne k depth a b t
    (RingHom.ker (aeval ξ : Q K →ₐ[K] F)) (X_not_mem_ker model.ne_zero)
    (S3poly_mem_ker model.on_curve) D
    (fei_of_spread (F₀ := F₀) ξ p k (D - 1) (Lam (D - 1)) hp model.spread (lambdaBound_Lam _)
      (Lam_pos _) (by omega)) hf
  exact h.2.2

/-! ## Last fall degree (Corollary 5.2) -/

/-- **Corollary 5.2.** For a unit ideal, under the hypotheses of the growth theorem, no
`e ≤ D` is a last fall degree bound when `k ≥ k_D`. -/
theorem growth_lfd_gt (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k)
    (unit : (1 : MvPolynomial (Fin (2 * k)) K) ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k depth a b t)) e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  exact growth_theorem two_ne k depth a b t ξ p model hp D hD hk (derived.mono he)

/-- **Linear growth.** For a unit ideal with a `CurveModel` and `k ≥ 24`: no
`e ≤ ⌊(k - 10)/14⌋` is a last fall degree bound, i.e. `LFD(k) ≥ ⌊(k - 10)/14⌋ + 1`. -/
theorem growth_lfd_gt_linear (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (hk : 24 ≤ k)
    (unit : (1 : MvPolynomial (Fin (2 * k)) K) ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e ≤ (k - 10) / 14, ¬ HasLFDBound (Set.range (inputs k depth a b t)) e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  have hD : 1 ≤ (k - 10) / 14 := by omega
  have hk' : 14 * ((k - 10) / 14) + 10 ≤ k := by omega
  exact growth_theorem_linear two_ne k depth a b t ξ p model hp _ hD hk' (derived.mono he)

/-- **Two-sided window.** For a nonsingular curve, a unit ideal and a `CurveModel`, every last
fall degree `d` of the two-summand system at depth `k ≥ k_D` satisfies `D < d ≤ k + 2`
(lower bound: this file; upper bound: `SelectorExtension.lfd_le_depth_add_two_of_two_ne`). -/
theorem growth_lfd_window (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (nonsingular : 4 * a ^ 3 + 27 * b ^ 2 ≠ 0)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k)
    (unit : (1 : MvPolynomial (Fin (2 * k)) K) ∈ Ideal.span (Set.range (inputs k depth a b t)))
    (d : ℕ) (hd : IsLastFallDegree (Set.range (inputs k depth a b t)) d) :
    D < d ∧ d ≤ k + 2 := by
  constructor
  · by_contra hle
    exact growth_lfd_gt two_ne k depth a b t ξ p model hp D hD hk unit d (by omega) hd.1
  · by_contra hgt
    exact hd.2 (k + 2) (by omega)
      (SelectorExtension.lfd_le_depth_add_two_of_two_ne two_ne k depth a b t nonsingular)

/-- **Prime-field instance.** `p` prime with `2^(k+1) ∣ p - 1`, the curve `y² = x³ + a x`, a
nonsquare target abscissa `t` on the curve side (`Y² = t³ + a t`, `Y ≠ 0`): the ideal is the
unit ideal (`UnitInstance.unit_zmod`), and with a `CurveModel` for the characteristic `p` no
`e ≤ D` is a last fall degree bound when `k ≥ k_D`. -/
theorem growth_lfd_gt_zmod {F₀ F : Type*} [Field F₀] [Field F] (p : ℕ) [Fact p.Prime]
    [Algebra (ZMod p) F₀] [Algebra F₀ F] [Algebra (ZMod p) F] [IsScalarTower (ZMod p) F₀ F]
    (k : ℕ) (depth : 2 ≤ k) (hdvd : 2 ^ (k + 1) ∣ p - 1) {a t Y : ZMod p}
    (curve : Y ^ 2 = t ^ 3 + a * t) (hY : Y ≠ 0) (nonsquare : ¬ IsSquare t)
    (ξ : Fin 2 → F) (model : CurveModel F₀ a 0 t p ξ)
    (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k depth a 0 t)) e := by
  have prime : p.Prime := Fact.out
  have hpos : 0 < p - 1 := by have := prime.two_le; omega
  have hle : 2 ^ (k + 1) ≤ p - 1 := Nat.le_of_dvd hpos hdvd
  have h4 : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
  have hsucc : 2 ^ (k + 1) = 2 * 2 ^ k := pow_succ' 2 k
  have large : 5 ≤ p := by omega
  exact growth_lfd_gt (PrimeFieldConclusions.two_ne_zero p large) k depth a 0 t ξ p model
    (Or.inr (by omega)) D hD hk (UnitInstance.unit_zmod p k depth hdvd curve hY nonsquare)

/-! ## The first values of `Λ` -/

/-- A recursion for `C(n, ≤ w)` that the kernel evaluates quickly. -/
def csumFast (n : ℕ) : ℕ → ℕ
  | 0 => 1
  | w + 1 => csumFast n w + n.descFactorial (w + 1) / (w + 1).factorial

theorem csumFast_eq (n w : ℕ) : csumFast n w = Csum n w := by
  induction w with
  | zero => simp [csumFast, Csum]
  | succ w ih =>
    unfold Csum at ih ⊢
    rw [csumFast, ih, Finset.sum_range_succ _ (w + 1), Nat.choose_eq_descFactorial_div_factorial]

/-- The defining inequality of `Λ`, in kernel-evaluable form. -/
def lamTest (w ℓ : ℕ) : Bool := decide (2 ^ (ℓ - 1) < 7776 * csumFast (2 * ℓ) w ^ 2)

theorem Lam_eq_of (w n : ℕ) (hn : n ≤ 13 * w + 23) (hP : lamTest w n = true)
    (hnot : (List.range (13 * w + 23 - n)).all (fun i => !lamTest w (n + 1 + i)) = true) :
    Lam w = n := by
  unfold Lam
  rw [Nat.findGreatest_eq_iff]
  refine ⟨hn, fun _ => ?_, ?_⟩
  · simpa [lamTest, csumFast_eq] using hP
  · intro ℓ h1 h2
    rw [List.all_eq_true] at hnot
    have := hnot (ℓ - (n + 1)) (List.mem_range.mpr (by omega))
    rw [show n + 1 + (ℓ - (n + 1)) = ℓ by omega] at this
    simpa [lamTest, csumFast_eq] using this

theorem Lam_one : Lam 1 = 25 := Lam_eq_of 1 25 (by norm_num) (by decide +kernel) (by decide +kernel)
theorem Lam_two : Lam 2 = 36 := Lam_eq_of 2 36 (by norm_num) (by decide +kernel) (by decide +kernel)
theorem Lam_three : Lam 3 = 48 :=
  Lam_eq_of 3 48 (by norm_num) (by decide +kernel) (by decide +kernel)
theorem Lam_four : Lam 4 = 59 :=
  Lam_eq_of 4 59 (by norm_num) (by decide +kernel) (by decide +kernel)
theorem Lam_five : Lam 5 = 71 :=
  Lam_eq_of 5 71 (by norm_num) (by decide +kernel) (by decide +kernel)

/-- The growth theorem at its first thresholds: `1 ∉ V_D(k)` for
`(D, k ≥ k_D) = (2, 27), (3, 39), (4, 52), (5, 64), (6, 77)`. -/
theorem growth_theorem_small (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (ξ : Fin 2 → F) (p : ℕ) (model : CurveModel F₀ a b t p ξ)
    (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p) (D : ℕ)
    (hk : (D = 2 ∧ 27 ≤ k) ∨ (D = 3 ∧ 39 ≤ k) ∨ (D = 4 ∧ 52 ≤ k) ∨ (D = 5 ∧ 64 ≤ k) ∨
      (D = 6 ∧ 77 ≤ k)) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : MvPolynomial (Fin (2 * k)) K) := by
  rcases hk with ⟨rfl, h⟩ | ⟨rfl, h⟩ | ⟨rfl, h⟩ | ⟨rfl, h⟩ | ⟨rfl, h⟩
  · exact growth_theorem two_ne k depth a b t ξ p model hp 2 (by norm_num)
      (by rw [show 2 - 1 = 1 from rfl, Lam_one]; omega)
  · exact growth_theorem two_ne k depth a b t ξ p model hp 3 (by norm_num)
      (by rw [show 3 - 1 = 2 from rfl, Lam_two]; omega)
  · exact growth_theorem two_ne k depth a b t ξ p model hp 4 (by norm_num)
      (by rw [show 4 - 1 = 3 from rfl, Lam_three]; omega)
  · exact growth_theorem two_ne k depth a b t ξ p model hp 5 (by norm_num)
      (by rw [show 5 - 1 = 4 from rfl, Lam_four]; omega)
  · exact growth_theorem two_ne k depth a b t ξ p model hp 6 (by norm_num)
      (by rw [show 6 - 1 = 5 from rfl, Lam_five]; omega)

end ICFallDegrees.GrowthMain
