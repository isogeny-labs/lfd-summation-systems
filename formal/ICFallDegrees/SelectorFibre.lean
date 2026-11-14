import ICFallDegrees.SelectorGrid
import ICFallDegrees.SelectorUnivariate
import Mathlib.RingTheory.MvPolynomial.Basic

noncomputable section

/-!
# Selectors and the fibre lemma of the selector bound

`E m k i α` is the character idempotent of the value `α` of summand `i`
(`docs/proofs/U1-upper-bounds.md`, Section 2.1): a polynomial of degree at most `k` in the
chain variables of summand `i` whose value at a chain point is `1` if `x i = α` and `0`
otherwise.

For the normalized two-summand presentation, the value of the generator at a chain point is
the summation polynomial `S_3(α, β, t) = c₂(α) β^2 + c₁(α) β + c₀(α)`. If the cubic is
nonsingular the three coefficients never vanish together (Lemma 2.2), and the fibre lemma
(Lemmas 2.3 and 2.4) derives, at degree `k + 2`, every row `E α * (w j - d₀ - d₁ * w 0)` and
`E α * (w 0 * w j - d₀ - d₁ * w 0)` whose bracket vanishes on the fibre of `α`.
-/

open Classical

namespace ICFallDegrees.SelectorFibre

open MvPolynomial FiniteQuotient ChainQuotient SelectorGrid
open PrimeM2General (inputs normalized active tail degree_and_top)
open GeneralQuadraticTop (top)

variable {K : Type*} [Field K]

/-! ## The selector -/

theorem geometric (γ : K) (k : ℕ) :
    (γ - 1)*∏ j ∈ Finset.range k, (1 + γ^2^j) = γ^2^k - 1 := by
  induction k with
  | zero => simp
  | succ k ih =>
      have sq : γ^2^(k+1) = (γ^2^k)^2 := by rw [pow_succ, pow_mul]
      rw [Finset.prod_range_succ, ← mul_assoc, ih, sq]
      ring

/-- The character idempotent of the value `α` of summand `i`. -/
def E (m k : ℕ) (i : Fin m) (α : K) : P (m*k) K :=
  C (((2 : K)^k)⁻¹)*∏ j : Fin k, (1 + C ((α⁻¹)^2^(j : ℕ))*X (var m k i j))

theorem E_totalDegree_le (m k : ℕ) (i : Fin m) (α : K) : (E m k i α).totalDegree ≤ k := by
  have factor : ∀ j : Fin k,
      (1 + C ((α⁻¹)^2^(j : ℕ))*X (var m k i j) : P (m*k) K).totalDegree ≤ 1 := by
    intro j
    refine (totalDegree_add _ _).trans (max_le (by simp) ?_)
    refine (totalDegree_mul _ _).trans ?_
    rw [totalDegree_C, zero_add, totalDegree_X]
  unfold E
  refine (totalDegree_mul _ _).trans ?_
  rw [totalDegree_C, zero_add]
  refine (totalDegree_finset_prod _ _).trans ?_
  calc ∑ j : Fin k, (1 + C ((α⁻¹)^2^(j : ℕ))*X (var m k i j) : P (m*k) K).totalDegree
      ≤ ∑ _j : Fin k, 1 := Finset.sum_le_sum (fun j _ => factor j)
    _ = k := by simp

theorem two_pow_ne_zero {k : ℕ} {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) : ((2 : K)^k) ≠ 0 := by
  haveI : NeZero (2^k : ℕ) := ⟨by positivity⟩
  have h : ((2^k : ℕ) : K) ≠ 0 := root.neZero'.ne
  push_cast at h
  exact h

theorem two_ne_zero_of_root {k : ℕ} (hk : 0 < k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) :
    (2 : K) ≠ 0 := by
  intro h
  apply two_pow_ne_zero root
  rw [h, zero_pow (by omega)]

theorem eval_E (m k : ℕ) {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) (i : Fin m) (α : K)
    (hα : α^2^k = 1) (x : Fin m → K) (hx : x i^2^k = 1) :
    eval (pt x) (E m k i α) = if x i = α then 1 else 0 := by
  have hN := two_pow_ne_zero root
  have αne : α ≠ 0 := by
    intro h
    rw [h, zero_pow (by positivity)] at hα
    exact zero_ne_one hα
  have factor : ∀ j : Fin k, eval (pt x) (1 + C ((α⁻¹)^2^(j : ℕ))*X (var m k i j)) =
      1 + (α⁻¹*x i)^2^(j : ℕ) := by
    intro j
    rw [map_add, map_one, map_mul, eval_C, eval_X, pt_var, mul_pow]
  unfold E
  rw [map_mul, eval_C, map_prod]
  simp only [factor]
  rw [Fin.prod_univ_eq_prod_range (fun j => 1 + (α⁻¹*x i)^2^j) k]
  have geo := geometric (α⁻¹*x i) k
  have γN : (α⁻¹*x i)^2^k = 1 := by rw [mul_pow, inv_pow, hα, hx, inv_one, one_mul]
  split_ifs with h
  · have γ1 : α⁻¹*x i = 1 := by rw [h]; exact inv_mul_cancel₀ αne
    rw [γ1]
    simp only [one_pow]
    rw [Finset.prod_const, Finset.card_range, one_add_one_eq_two]
    exact inv_mul_cancel₀ hN
  · have γne : α⁻¹*x i - 1 ≠ 0 := by
      intro h'
      apply h
      have one := sub_eq_zero.mp h'
      calc x i = α*(α⁻¹*x i) := by rw [← mul_assoc, mul_inv_cancel₀ αne, one_mul]
        _ = α := by rw [one, mul_one]
    rw [γN, sub_self] at geo
    rcases mul_eq_zero.mp geo with h' | h'
    · exact absurd h' γne
    · rw [h', mul_zero]

/-! ## The fibre polynomial of the normalized generator -/

/-- Coefficient of `β^2` in `S_3(α, β, t)`. -/
def c₂ (t α : K) : K := (α - t)^2
/-- Coefficient of `β` in `S_3(α, β, t)`. -/
def c₁ (a b t α : K) : K := -2*(t^2 + a)*α - 2*t*α^2 - 2*a*t - 4*b
/-- Constant coefficient of `S_3(α, β, t)` in `β`. -/
def c₀ (a b t α : K) : K := t^2*α^2 + (-2*a*t - 4*b)*α + a^2 - 4*b*t

/-- Lemma 2.2: for a nonsingular cubic the fibre polynomial is never identically zero.
Only `2 ≠ 0` is used of the characteristic. -/
theorem fibre_nonzero (two_ne : (2 : K) ≠ 0) (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0)
    (α : K) : ¬(c₂ t α = 0 ∧ c₁ a b t α = 0 ∧ c₀ a b t α = 0) := by
  rintro ⟨h2, h1, h0⟩
  have hα : α = t := by
    unfold c₂ at h2
    exact sub_eq_zero.mp ((pow_eq_zero_iff (by norm_num : (2 : ℕ) ≠ 0)).mp h2)
  rw [hα] at h1 h0
  unfold c₁ at h1
  unfold c₀ at h0
  have four : (4 : K) ≠ 0 := by
    have e : (4 : K) = 2*2 := by norm_num
    rw [e]
    exact mul_ne_zero two_ne two_ne
  have f0 : t^3 + a*t + b = 0 := by
    have h : (4 : K)*(t^3 + a*t + b) = 0 := by linear_combination -h1
    exact (mul_eq_zero.mp h).resolve_left four
  have g0 : 3*t^2 + a = 0 := by
    have h : (3*t^2 + a)^2 = 0 := by linear_combination h0 + 8*t*f0
    exact (pow_eq_zero_iff (by norm_num : (2 : ℕ) ≠ 0)).mp h
  apply nonsingular
  linear_combination (27*b - 27*a*t - 27*t^3)*f0 + (4*a^2 + 15*a*t^2 + 9*t^4)*g0

section Presentation

variable (k : ℕ) (depth : 2 ≤ k) (a b t : K)

theorem active_zero : active k depth 0 = var 2 k 0 ⟨0, by omega⟩ := by
  apply Fin.ext
  show (0 : ℕ) = 0 + k*0
  omega

theorem active_one : active k depth 1 = var 2 k 0 ⟨1, by omega⟩ := by
  apply Fin.ext
  show (1 : ℕ) = 1 + k*0
  omega

theorem active_two : active k depth 2 = var 2 k 1 ⟨0, by omega⟩ := by
  apply Fin.ext
  show k = 0 + k*1
  omega

theorem active_three : active k depth 3 = var 2 k 1 ⟨1, by omega⟩ := by
  apply Fin.ext
  show k+1 = 1 + k*1
  omega

/-- The normalized generator at a chain point is `S_3` of the two summand values. -/
theorem eval_normalized (x : Fin 2 → K) :
    eval (pt x) (normalized k depth a b t) =
      c₂ t (x 0)*(x 1)^2 + c₁ a b t (x 0)*x 1 + c₀ a b t (x 0) := by
  simp only [normalized, top, tail, active_zero, active_one, active_two, active_three, map_add,
    map_mul, eval_C, eval_X, pt_var, c₂, c₁, c₀, Fin.val_mk, pow_zero, pow_one]
  ring

/-- The chain of summand 1, with `w k = 1`. -/
def w (k : ℕ) (j : ℕ) : P (2*k) K := if h : j < k then X (var 2 k 1 ⟨j, h⟩) else 1

theorem w_lt {j : ℕ} (h : j < k) : (w k j : P (2*k) K) = X (var 2 k 1 ⟨j, h⟩) := dif_pos h

theorem w_top : (w k k : P (2*k) K) = 1 := dif_neg (lt_irrefl k)

theorem w_totalDegree_le (j : ℕ) : (w k j : P (2*k) K).totalDegree ≤ 1 := by
  unfold w
  split_ifs <;> simp

theorem eval_w (x : Fin 2 → K) (hx : x 1^2^k = 1) {j : ℕ} (hj : j ≤ k) :
    eval (pt x) (w k j) = x 1^2^j := by
  rcases Nat.lt_or_ge j k with h | h
  · rw [w_lt k h, eval_X, pt_var]
  · have e : j = k := le_antisymm hj h
    rw [e, w_top, map_one, hx]

/-- Rows that may still be multiplied by one chain variable at budget `k + 2`. -/
def lowSpace (k : ℕ) : Submodule K (P (2*k) K) := restrictTotalDegree (Fin (2*k)) K (k+1)

theorem mem_lowSpace (p : P (2*k) K) : p ∈ lowSpace k ↔ p.totalDegree ≤ k+1 := by
  rw [lowSpace, mem_restrictTotalDegree]

variable {k} {a b t}

/-- The fibre lemma (Lemmas 2.3 and 2.4 of the written proof) for the value `α` of the first
summand: rows of the block-2 fibre, multiplied by the selector, at degree `k + 2`. -/
theorem fibre {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) (nonsingular : 4*a^3 + 27*b^2 ≠ 0)
    (α : K) (hα : α^2^k = 1) :
    ∀ j ≤ k, ∀ d₀ d₁ : K,
      ((∀ β : K, β^2^k = 1 → c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α = 0 → d₀ + d₁*β = β^2^j) →
        InClosure (Set.range (inputs k depth a b t)) (k+2)
          (E 2 k 0 α*(w k j - C d₀ - C d₁*w k 0))) ∧
      ((∀ β : K, β^2^k = 1 → c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α = 0 →
          d₀ + d₁*β = β^(2^j+1)) →
        InClosure (Set.range (inputs k depth a b t)) (k+2)
          (E 2 k 0 α*(w k 0*w k j - C d₀ - C d₁*w k 0))) := by
  classical
  have two_ne := two_ne_zero_of_root (by omega : 0 < k) root
  obtain ⟨F, hF⟩ : ∃ F : Set (P (2*k) K), F = Set.range (inputs k depth a b t) := ⟨_, rfl⟩
  rw [← hF]
  obtain ⟨e, he⟩ : ∃ e : P (2*k) K, e = E 2 k 0 α := ⟨_, rfl⟩
  rw [← he]
  have edeg : e.totalDegree ≤ k := by rw [he]; exact E_totalDegree_le 2 k 0 α
  have eval_e : ∀ x : Fin 2 → K, (∀ i, x i^2^k = 1) →
      eval (pt x) e = if x 0 = α then 1 else 0 := by
    intro x hx
    rw [he]
    exact eval_E 2 k root 0 α hα x (hx 0)
  let Cl : Submodule K (P (2*k) K) := closureSpace F (k+2)
  let Low : Submodule K (P (2*k) K) := lowSpace k
  have included : chainInputs K 2 k ⊆ F := by
    rw [hF]; exact chainInputs_two_subset k depth a b t
  have gridMem : ∀ p : P (2*k) K, p.totalDegree ≤ k+2 →
      (∀ x : Fin 2 → K, (∀ i, x i^2^k = 1) → eval (pt x) p = 0) → p ∈ Cl :=
    fun p degree vanish => grid_mem 2 k root F included (k+2) p degree vanish
  have scale : ∀ (c : K) (y : P (2*k) K), y ∈ Cl → C c*y ∈ Cl := by
    intro c y hy
    rw [← smul_eq_C_mul]
    exact Cl.smul_mem c hy
  have wdeg : ∀ j, (w k j : P (2*k) K).totalDegree ≤ 1 := w_totalDegree_le k
  have lowE : e ∈ Low := (mem_lowSpace k e).mpr (by omega)
  have low : ∀ j, e*w k j ∈ Low := by
    intro j
    refine (mem_lowSpace k _).mpr ((totalDegree_mul _ _).trans ?_)
    have := wdeg j
    omega
  have mulH : ∀ j < k, ∀ y ∈ Cl, y ∈ Low → w k j*y ∈ Cl := by
    intro j hj y hy hlow
    rw [w_lt k hj]
    have degree := (mem_lowSpace k y).mp hlow
    exact InClosure.mul (X _) hy (by rw [totalDegree_X]; omega)
  have sqH : ∀ j < k, e*(w k j^2 - w k (j+1)) ∈ Cl := by
    intro j hj
    apply gridMem
    · refine (totalDegree_mul _ _).trans ?_
      have h1 : (w k j^2 - w k (j+1) : P (2*k) K).totalDegree ≤ 2 := by
        refine (totalDegree_sub_le _ _).trans (max_le ?_ ((wdeg _).trans (by norm_num)))
        refine (totalDegree_pow _ _).trans ?_
        have := wdeg j
        omega
      omega
    · intro x hx
      rw [map_mul, map_sub, map_pow, eval_w k x (hx 1) hj.le, eval_w k x (hx 1) hj,
        ← pow_mul, ← pow_succ, sub_self, mul_zero]
  have hk0 : 0 < k := by omega
  have hk1 : 1 < k := by omega
  -- the fibre row: the generator times the selector, with the first summand evaluated at α
  have fibre_row : e*(C (c₂ t α)*w k 1 + C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) ∈ Cl := by
    have generator : e*normalized k depth a b t ∈ Cl := by
      have member : inputs k depth a b t (Fin.last (2*k)) ∈ F := by
        rw [hF]; exact Set.mem_range_self _
      have degree := (degree_and_top k depth a b t (Fin.last (2*k))).1
      have last : inputs k depth a b t (Fin.last (2*k)) = normalized k depth a b t := by
        simp only [inputs, Fin.lastCases_last]
      rw [last] at member degree
      exact InClosure.mul e (.input member (by omega)) (by omega)
    have difference : e*normalized k depth a b t -
        e*(C (c₂ t α)*w k 1 + C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) ∈ Cl := by
      apply gridMem
      · refine (totalDegree_sub_le _ _).trans (max_le ?_ ?_)
        · refine (totalDegree_mul _ _).trans ?_
          have := (degree_and_top k depth a b t (Fin.last (2*k))).1
          simp only [inputs, Fin.lastCases_last] at this
          omega
        · refine (totalDegree_mul _ _).trans ?_
          have h1 : (C (c₂ t α)*w k 1 + C (c₁ a b t α)*w k 0 + C (c₀ a b t α) :
              P (2*k) K).totalDegree ≤ 1 := by
            refine (totalDegree_add _ _).trans (max_le ?_ (by simp))
            refine (totalDegree_add _ _).trans (max_le ?_ ?_)
            · refine (totalDegree_mul _ _).trans ?_
              rw [totalDegree_C, zero_add]; exact wdeg 1
            · refine (totalDegree_mul _ _).trans ?_
              rw [totalDegree_C, zero_add]; exact wdeg 0
          omega
      · intro x hx
        rw [map_sub, map_mul, map_mul, eval_e x hx, eval_normalized, map_add, map_add, map_mul,
          map_mul, eval_C, eval_C, eval_C, eval_w k x (hx 1) hk1.le, eval_w k x (hx 1) hk0.le]
        split_ifs with h
        · rw [h]
          simp only [pow_zero, pow_one]
          ring
        · ring
    have h := Cl.sub_mem generator difference
    rwa [sub_sub_cancel] at h
  -- a monic quadratic with the same roots among the roots of unity
  obtain ⟨s, n, quad, same⟩ : ∃ s n : K, e*(w k 1 - C s*w k 0 + C n) ∈ Cl ∧
      ∀ β : K, β^2^k = 1 → (β^2 - s*β + n = 0 ↔
        c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α = 0) := by
    have βne : ∀ β : K, β^2^k = 1 → β ≠ 0 := by
      intro β hβ h
      rw [h, zero_pow (by positivity)] at hβ
      exact zero_ne_one hβ
    by_cases h2 : c₂ t α = 0
    · by_cases h1 : c₁ a b t α = 0
      · have h0 : c₀ a b t α ≠ 0 := fun h0 => fibre_nonzero two_ne a b t nonsingular α ⟨h2, h1, h0⟩
        obtain ⟨i₀, hi₀⟩ : ∃ i₀ : K, i₀*c₀ a b t α = 1 := ⟨_, inv_mul_cancel₀ h0⟩
        have hC : (C i₀*C (c₀ a b t α) : P (2*k) K) = 1 := by rw [← C_mul, hi₀, C_1]
        rw [h2, h1] at fibre_row
        have he' : e ∈ Cl := by
          have h := scale i₀ _ fibre_row
          have ident : e = C i₀*(e*(C 0*w k 1 + C 0*w k 0 + C (c₀ a b t α))) := by
            simp only [map_zero, zero_mul, zero_add]
            linear_combination (-e)*hC
          rw [ident]; exact h
        refine ⟨0, 0, ?_, ?_⟩
        · have h := mulH 1 hk1 e he' lowE
          have ident : e*(w k 1 - C 0*w k 0 + C 0) = w k 1*e := by
            simp only [map_zero, zero_mul, sub_zero, add_zero]
            ring
          rw [ident]; exact h
        · intro β hβ
          rw [h2, h1]
          constructor
          · intro h
            have zero : β^2 = 0 := by linear_combination h
            exact absurd ((pow_eq_zero_iff (by norm_num : (2 : ℕ) ≠ 0)).mp zero) (βne β hβ)
          · intro h
            have zero : c₀ a b t α = 0 := by linear_combination h
            exact absurd zero h0
      · obtain ⟨i₁, hi₁⟩ : ∃ i₁ : K, i₁*c₁ a b t α = 1 := ⟨_, inv_mul_cancel₀ h1⟩
        have hC : (C i₁*C (c₁ a b t α) : P (2*k) K) = 1 := by rw [← C_mul, hi₁, C_1]
        have i₁ne : i₁ ≠ 0 := by
          intro h
          rw [h, zero_mul] at hi₁
          exact zero_ne_one hi₁
        rw [h2] at fibre_row
        have linear : e*(C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) ∈ Cl := by
          have ident : e*(C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) =
              e*(C 0*w k 1 + C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) := by
            simp only [map_zero, zero_mul, zero_add]
          rw [ident]; exact fibre_row
        have linear_low : e*(C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) ∈ Low := by
          have h := Low.add_mem (Low.smul_mem (c₁ a b t α) (low 0))
            (Low.smul_mem (c₀ a b t α) lowE)
          have ident : e*(C (c₁ a b t α)*w k 0 + C (c₀ a b t α)) =
              c₁ a b t α • (e*w k 0) + c₀ a b t α • e := by
            simp only [smul_eq_C_mul]
            ring
          rw [ident]; exact h
        refine ⟨-(i₁*c₀ a b t α), 0, ?_, ?_⟩
        · have h := Cl.sub_mem (scale i₁ _ (mulH 0 hk0 _ linear linear_low)) (sqH 0 hk0)
          have ident : e*(w k 1 - C (-(i₁*c₀ a b t α))*w k 0 + C 0) =
              C i₁*(w k 0*(e*(C (c₁ a b t α)*w k 0 + C (c₀ a b t α)))) -
                e*(w k 0^2 - w k (0+1)) := by
            simp only [map_neg, map_mul, map_zero, add_zero, zero_add]
            linear_combination (-(e*w k 0^2))*hC
          rw [ident]; exact h
        · intro β hβ
          rw [h2]
          have ident : β^2 - -(i₁*c₀ a b t α)*β + 0 =
              (i₁*β)*(0*β^2 + c₁ a b t α*β + c₀ a b t α) := by
            linear_combination (-β^2)*hi₁
          rw [ident]
          constructor
          · intro h
            rcases mul_eq_zero.mp h with h | h
            · exact absurd h (mul_ne_zero i₁ne (βne β hβ))
            · exact h
          · intro h
            rw [h, mul_zero]
    · obtain ⟨i₂, hi₂⟩ : ∃ i₂ : K, i₂*c₂ t α = 1 := ⟨_, inv_mul_cancel₀ h2⟩
      have hC : (C i₂*C (c₂ t α) : P (2*k) K) = 1 := by rw [← C_mul, hi₂, C_1]
      have i₂ne : i₂ ≠ 0 := by
        intro h
        rw [h, zero_mul] at hi₂
        exact zero_ne_one hi₂
      refine ⟨-(i₂*c₁ a b t α), i₂*c₀ a b t α, ?_, ?_⟩
      · have h := scale i₂ _ fibre_row
        have ident : e*(w k 1 - C (-(i₂*c₁ a b t α))*w k 0 + C (i₂*c₀ a b t α)) =
            C i₂*(e*(C (c₂ t α)*w k 1 + C (c₁ a b t α)*w k 0 + C (c₀ a b t α))) := by
          simp only [map_neg, map_mul]
          linear_combination (-(e*w k 1))*hC
        rw [ident]; exact h
      · intro β hβ
        have ident : β^2 - -(i₂*c₁ a b t α)*β + i₂*c₀ a b t α =
            i₂*(c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α) := by
          linear_combination (-β^2)*hi₂
        rw [ident]
        constructor
        · intro h
          rcases mul_eq_zero.mp h with h | h
          · exact absurd h i₂ne
          · exact h
        · intro h
          rw [h, mul_zero]
  have all := SelectorUnivariate.fibre_rows Cl Low e (w k) k s n hk0 root lowE low mulH sqH
    (by rw [algebraMap_eq]; exact quad) (w_top k)
  simp only [algebraMap_eq] at all
  intro j hj d₀ d₁
  obtain ⟨first, second⟩ := all j hj d₀ d₁
  constructor
  · intro values
    exact first (fun β hβ hQ => values β hβ ((same β hβ).mp hQ))
  · intro values
    exact second (fun β hβ hQ => values β hβ ((same β hβ).mp hQ))

end Presentation

end ICFallDegrees.SelectorFibre
