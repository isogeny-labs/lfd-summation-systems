import ICFallDegrees.GrowthMain
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.FieldTheory.Finite.Basic

noncomputable section

/-!
# Growth theorem: unit-ideal instances at every depth
(Corollary 5.2 (c) of `docs/proofs/N8-growth.md`)

For every depth `k ≥ 2`, every prime factor `p` of the Fermat number `2^(2^k) + 1` and all
non-squares `X, n ∈ F_p`, the two-summand chain presentation of the curve
`y² = x³ + (n - X²) x` with target `t = X` generates the unit ideal
(`unit_instance`, **kernel**, unconditional).  With a `GrowthMain.CurveModel` for that instance
(the analytic hypothesis, not proved in Lean) its last fall degree exceeds every `D` with
`k_D ≤ k` and lies in `(⌊(k-10)/14⌋, k + 2]` (`growth_instance`, `growth_instance_window`).
-/

namespace ICFallDegrees.GrowthExistence

open MvPolynomial GrowthCounting GrowthMain
open PrimeM2General (inputs)

/-- A prime factor `p` of `2^(2^k) + 1` satisfies `2^(k+1) ∣ p - 1`. -/
theorem fermat_factor_dvd (p : ℕ) [Fact p.Prime] (k : ℕ) (h : p ∣ 2 ^ 2 ^ k + 1) :
    2 ^ (k + 1) ∣ p - 1 := by
  have prime : p.Prime := Fact.out
  have hp2 : p ≠ 2 := by
    rintro rfl
    have hpos : 0 < 2 ^ k := Nat.two_pow_pos k
    have : 2 ∣ 2 ^ 2 ^ k := dvd_pow_self 2 (by omega)
    omega
  haveI : Fact (2 < p) := ⟨by have := prime.two_le; omega⟩
  have hz : ((2 ^ 2 ^ k + 1 : ℕ) : ZMod p) = 0 := (ZMod.natCast_zmod_eq_zero_iff_dvd _ p).mpr h
  have hpow : (2 : ZMod p) ^ 2 ^ k = -1 := by
    push_cast at hz
    exact eq_neg_of_add_eq_zero_left hz
  have hnot : ¬ (2 : ZMod p) ^ 2 ^ k = 1 := by
    rw [hpow]
    exact ZMod.neg_one_ne_one
  have hfin : (2 : ZMod p) ^ 2 ^ (k + 1) = 1 := by
    rw [pow_succ, pow_mul, hpow]
    norm_num
  have order := orderOf_eq_prime_pow hnot hfin
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by positivity)] at hfin
    exact zero_ne_one hfin
  rw [← order]
  exact ZMod.orderOf_dvd_card_sub_one two_ne

/-- **Existence of unit-ideal instances (Corollary 5.2 (c)), unconditional.** -/
theorem unit_instance (p : ℕ) [Fact p.Prime] (k : ℕ) (depth : 2 ≤ k) (h : p ∣ 2 ^ 2 ^ k + 1)
    {X n : ZMod p} (hX : ¬ IsSquare X) (hn : ¬ IsSquare n) :
    (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∈
      Ideal.span (Set.range (inputs k depth (n - X ^ 2) 0 X)) := by
  have hX0 : X ≠ 0 := by
    rintro rfl
    exact hX (IsSquare.zero)
  have hn0 : n ≠ 0 := by
    rintro rfl
    exact hn (IsSquare.zero)
  have hprod : IsSquare (X * n) := by
    rw [← quadraticChar_one_iff_isSquare (mul_ne_zero hX0 hn0), map_mul,
      quadraticChar_neg_one_iff_not_isSquare.mpr hX, quadraticChar_neg_one_iff_not_isSquare.mpr hn]
    norm_num
  obtain ⟨Y, hY⟩ := hprod
  have curve : Y ^ 2 = X ^ 3 + (n - X ^ 2) * X := by
    rw [pow_two, ← hY]
    ring
  have hY0 : Y ≠ 0 := by
    rintro rfl
    rw [mul_zero] at hY
    exact mul_ne_zero hX0 hn0 hY
  exact UnitInstance.unit_zmod p k depth (fermat_factor_dvd p k h) curve hY0 hX

/-- The curve of the instance is nonsingular: `4 a³ + 27·0² ≠ 0` with `a = n - X²`. -/
theorem instance_nonsingular (p : ℕ) [Fact p.Prime] (hp2 : p ≠ 2) {X n : ZMod p}
    (hn : ¬ IsSquare n) : 4 * (n - X ^ 2) ^ 3 + 27 * (0 : ZMod p) ^ 2 ≠ 0 := by
  have prime : p.Prime := Fact.out
  have ha : n - X ^ 2 ≠ 0 := by
    intro h0
    exact hn ⟨X, by rw [sub_eq_zero.mp h0, pow_two]⟩
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h2
    have dvd : p ∣ 2 := (ZMod.natCast_zmod_eq_zero_iff_dvd 2 p).mp (by exact_mod_cast h2)
    exact hp2 ((Nat.prime_dvd_prime_iff_eq prime Nat.prime_two).mp dvd)
  have four : (4 : ZMod p) = 2 * 2 := by norm_num
  rw [zero_pow (by norm_num), mul_zero, add_zero, four]
  exact mul_ne_zero (mul_ne_zero two_ne two_ne) (pow_ne_zero 3 ha)

variable {F₀ F : Type*} [Field F₀] [Field F]

/-- **Growth on the instances.** For the instances of `unit_instance`, with a `CurveModel` for
the characteristic `p`: no `e ≤ D` is a last fall degree bound when `k ≥ k_D`. -/
theorem growth_instance (p : ℕ) [Fact p.Prime] [Algebra (ZMod p) F₀] [Algebra F₀ F]
    [Algebra (ZMod p) F] [IsScalarTower (ZMod p) F₀ F] (k : ℕ) (depth : 2 ≤ k)
    (h : p ∣ 2 ^ 2 ^ k + 1) {X n : ZMod p} (hX : ¬ IsSquare X) (hn : ¬ IsSquare n)
    (ξ : Fin 2 → F) (model : CurveModel F₀ (n - X ^ 2) 0 X p ξ)
    (D : ℕ) (hD : 1 ≤ D) (hk : D + Lam (D - 1) ≤ k) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k depth (n - X ^ 2) 0 X)) e := by
  have hdvd := fermat_factor_dvd p k h
  have prime : p.Prime := Fact.out
  have hpos : 0 < p - 1 := by have := prime.two_le; omega
  have hle : 2 ^ (k + 1) ≤ p - 1 := Nat.le_of_dvd hpos hdvd
  have h4 : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
  have hsucc : 2 ^ (k + 1) = 2 * 2 ^ k := pow_succ' 2 k
  have large : 5 ≤ p := by omega
  exact growth_lfd_gt (PrimeFieldConclusions.two_ne_zero p large) k depth _ 0 X ξ p model
    (Or.inr (by omega)) D hD hk (unit_instance p k depth h hX hn)

/-- **Two-sided window on the instances** (`k ≥ 24`): every last fall degree `d` of the
instance satisfies `⌊(k - 10)/14⌋ < d ≤ k + 2`. -/
theorem growth_instance_window (p : ℕ) [Fact p.Prime] [Algebra (ZMod p) F₀] [Algebra F₀ F]
    [Algebra (ZMod p) F] [IsScalarTower (ZMod p) F₀ F] (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (h : p ∣ 2 ^ 2 ^ k + 1) {X n : ZMod p} (hX : ¬ IsSquare X) (hn : ¬ IsSquare n)
    (ξ : Fin 2 → F) (model : CurveModel F₀ (n - X ^ 2) 0 X p ξ)
    (d : ℕ) (hd : IsLastFallDegree (Set.range (inputs k depth (n - X ^ 2) 0 X)) d) :
    (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  have hdvd := fermat_factor_dvd p k h
  have prime : p.Prime := Fact.out
  have hpos : 0 < p - 1 := by have := prime.two_le; omega
  have hle : 2 ^ (k + 1) ≤ p - 1 := Nat.le_of_dvd hpos hdvd
  have h4 : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
  have hsucc : 2 ^ (k + 1) = 2 * 2 ^ k := pow_succ' 2 k
  have large : 5 ≤ p := by omega
  have two_ne := PrimeFieldConclusions.two_ne_zero p large
  constructor
  · by_contra hle'
    exact growth_lfd_gt_linear two_ne k depth _ 0 X ξ p model (Or.inr (by omega)) hk
      (unit_instance p k depth h hX hn) d (by omega) hd.1
  · by_contra hgt
    exact hd.2 (k + 2) (by omega)
      (SelectorExtension.lfd_le_depth_add_two_of_two_ne two_ne k depth _ 0 X
        (instance_nonsingular p (by omega) hn))

end ICFallDegrees.GrowthExistence
