import ICFallDegrees.BinaryWeight
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Data.Nat.Log
import Mathlib.Data.Nat.Find
import Mathlib.Data.Finset.Prod
import Mathlib.Algebra.Order.Ring.Abs
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Growth theorem, counting part (Section 4 of `docs/proofs/N8-growth.md`)

Pure statements about natural numbers and binary digit sums (`BinaryWeight.wt`):

* `Csum n w = C(n, ≤ w)`, the count of numbers below `2^L` of weight `≤ w`
  (`card_wt_le`), and the count of pairs (`card_pairs_le`);
* `LambdaBound w lam`: the inequality `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2` fails for every `ℓ > lam`;
  `lambdaBound_linear`: it fails for every `ℓ > 13w + 23` (the bound `Λ(w) ≤ 13w + 23` of
  Lemma 4.1), proved here by the elementary estimate `C(n, ≤w)·16^n ≤ 17^n·16^w` and integer
  arithmetic, **not** by the ratio argument of the written proof (red-team finding F1 concerns
  that ratio argument; it is not used);
* `Lam w`, the largest `ℓ` satisfying the inequality (`Lam_spec`, `Lam_greatest`), which is the
  `Λ(w)` of the written proof, with `lambdaBound_Lam`;
* Lemma 4.2 (window) in the form needed: `window_index`;
* the count of Lemma 4.3: `card_le_of_window`, and its consequence `sigma_lt_of_lambda`;
* Lemma 4.4 (no straddling): `wt_ge_of_close`, `wt_ge_of_top`.

The distance `‖n‖_p` of the written proof is expressed by an existential multiple:
`‖n‖_p ≤ B` is `∃ z : ℤ, |n - p*z| ≤ B` (for `p = 0` this is `|n| ≤ B`).
-/

namespace ICFallDegrees.GrowthCounting

open BinaryWeight

/-! ## Partial binomial sums and the count of low-weight numbers -/

/-- `C(n, ≤ w) = Σ_{i ≤ w} binom(n, i)`. -/
def Csum (n w : ℕ) : ℕ := ∑ i ∈ Finset.range (w + 1), n.choose i

theorem Csum_zero_left (w : ℕ) : Csum 0 w = 1 := by
  unfold Csum
  rw [Finset.sum_range_succ']
  simp

theorem Csum_zero_right (n : ℕ) : Csum n 0 = 1 := by
  simp [Csum]

theorem Csum_succ_succ (n w : ℕ) : Csum (n + 1) (w + 1) = Csum n (w + 1) + Csum n w := by
  unfold Csum
  rw [Finset.sum_range_succ' _ (w + 1), Finset.sum_range_succ' _ (w + 1)]
  simp only [Nat.choose_zero_right, Nat.choose_succ_succ', Finset.sum_add_distrib]
  ring

theorem Csum_pos (n w : ℕ) : 1 ≤ Csum n w := by
  unfold Csum
  rw [Finset.sum_range_succ']
  simp

/-- The numbers below `2^L` of binary weight at most `w`. -/
def lowWeight (L w : ℕ) : Finset ℕ := (Finset.range (2 ^ L)).filter (fun x => wt x ≤ w)

theorem card_lowWeight_le : ∀ (L w : ℕ), (lowWeight L w).card ≤ Csum L w := by
  intro L
  induction L with
  | zero =>
    intro w
    rw [Csum_zero_left]
    unfold lowWeight
    calc _ ≤ (Finset.range (2 ^ 0)).card := Finset.card_filter_le _ _
      _ = 1 := by simp
  | succ L ih =>
    intro w
    have split : (lowWeight (L + 1) w).card =
        ((lowWeight (L + 1) w).filter (fun x => x % 2 = 0)).card +
          ((lowWeight (L + 1) w).filter (fun x => ¬ x % 2 = 0)).card :=
      (Finset.filter_card_add_filter_neg_card_eq_card _).symm
    have even : ((lowWeight (L + 1) w).filter (fun x => x % 2 = 0)).card ≤
        (lowWeight L w).card := by
      apply Finset.card_le_card_of_injOn (fun x => x / 2)
      · intro x hx
        simp only [lowWeight, Finset.coe_filter, Finset.mem_filter, Finset.mem_range,
          Set.mem_setOf_eq] at hx ⊢
        obtain ⟨⟨hlt, hwt⟩, _⟩ := hx
        rw [wt_eq] at hwt
        rw [pow_succ] at hlt
        exact ⟨by omega, by omega⟩
      · intro x hx y hy hxy
        simp only [Finset.coe_filter, Finset.mem_filter, Set.mem_setOf_eq] at hx hy
        have := hx.2
        have := hy.2
        simp only at hxy
        omega
    cases w with
    | zero =>
      have odd : ((lowWeight (L + 1) 0).filter (fun x => ¬ x % 2 = 0)).card = 0 := by
        rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
        intro x hx
        simp only [lowWeight, Finset.mem_filter, Finset.mem_range] at hx
        have hwt := hx.2
        rw [wt_eq] at hwt
        omega
      rw [split, odd, Csum_zero_right]
      have := ih 0
      rw [Csum_zero_right] at this
      omega
    | succ w =>
      have odd : ((lowWeight (L + 1) (w + 1)).filter (fun x => ¬ x % 2 = 0)).card ≤
          (lowWeight L w).card := by
        apply Finset.card_le_card_of_injOn (fun x => x / 2)
        · intro x hx
          simp only [lowWeight, Finset.coe_filter, Finset.mem_filter, Finset.mem_range,
            Set.mem_setOf_eq] at hx ⊢
          obtain ⟨⟨hlt, hwt⟩, hodd⟩ := hx
          rw [wt_eq] at hwt
          rw [pow_succ] at hlt
          exact ⟨by omega, by omega⟩
        · intro x hx y hy hxy
          simp only [Finset.coe_filter, Finset.mem_filter, Set.mem_setOf_eq] at hx hy
          have := hx.2
          have := hy.2
          simp only at hxy
          omega
      rw [split, Csum_succ_succ]
      have := ih (w + 1)
      have := ih w
      omega

/-- Pairs below `2^ℓ` of total weight at most `w`. -/
def lowPairs (ℓ w : ℕ) : Finset (ℕ × ℕ) :=
  (Finset.range (2 ^ ℓ) ×ˢ Finset.range (2 ^ ℓ)).filter (fun q => wt q.1 + wt q.2 ≤ w)

theorem card_lowPairs_le (ℓ w : ℕ) : (lowPairs ℓ w).card ≤ Csum (2 * ℓ) w := by
  refine le_trans ?_ (card_lowWeight_le (2 * ℓ) w)
  apply Finset.card_le_card_of_injOn (fun q => q.1 + 2 ^ ℓ * q.2)
  · intro q hq
    simp only [lowPairs, lowWeight, Finset.coe_filter, Finset.mem_filter, Finset.mem_product,
      Finset.mem_range, Set.mem_setOf_eq] at hq ⊢
    obtain ⟨⟨h1, h2⟩, hw⟩ := hq
    refine ⟨?_, ?_⟩
    · have hpow : 2 ^ (2 * ℓ) = 2 ^ ℓ * 2 ^ ℓ := by rw [two_mul, pow_add]
      have : 2 ^ ℓ * (q.2 + 1) ≤ 2 ^ ℓ * 2 ^ ℓ := Nat.mul_le_mul_left _ h2
      rw [hpow]
      nlinarith
    · rw [wt_add_mul_two_pow ℓ q.1 q.2 h1]
      exact hw
  · intro q hq q' hq' h
    simp only [lowPairs, Finset.coe_filter, Finset.mem_filter, Finset.mem_product,
      Finset.mem_range, Set.mem_setOf_eq] at hq hq'
    simp only at h
    have hpos : 0 < 2 ^ ℓ := Nat.two_pow_pos ℓ
    have m1 : (q.1 + 2 ^ ℓ * q.2) % 2 ^ ℓ = q.1 := by
      rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hq.1.1]
    have m2 : (q'.1 + 2 ^ ℓ * q'.2) % 2 ^ ℓ = q'.1 := by
      rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hq'.1.1]
    have d1 : (q.1 + 2 ^ ℓ * q.2) / 2 ^ ℓ = q.2 := by
      rw [Nat.add_mul_div_left _ _ hpos, Nat.div_eq_of_lt hq.1.1, zero_add]
    have d2 : (q'.1 + 2 ^ ℓ * q'.2) / 2 ^ ℓ = q'.2 := by
      rw [Nat.add_mul_div_left _ _ hpos, Nat.div_eq_of_lt hq'.1.1, zero_add]
    exact Prod.ext (by rw [← m1, ← m2, h]) (by rw [← d1, ← d2, h])

/-! ## The threshold `Λ(w)` (Lemma 4.1) -/

/-- `lam` is an upper bound for the set of `ℓ` with `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2`. -/
def LambdaBound (w lam : ℕ) : Prop :=
  ∀ ℓ, lam < ℓ → 7776 * Csum (2 * ℓ) w ^ 2 ≤ 2 ^ (ℓ - 1)

theorem LambdaBound.mono {w lam lam' : ℕ} (h : LambdaBound w lam) (hle : lam ≤ lam') :
    LambdaBound w lam' :=
  fun ℓ hℓ => h ℓ (lt_of_le_of_lt hle hℓ)

/-- An elementary tail estimate: `C(n, ≤w)·16^n ≤ 17^n·16^w`. -/
theorem Csum_chernoff : ∀ (n w : ℕ), Csum n w * 16 ^ n ≤ 17 ^ n * 16 ^ w := by
  intro n
  induction n with
  | zero =>
    intro w
    rw [Csum_zero_left]
    simpa using Nat.one_le_pow w 16 (by norm_num)
  | succ n ih =>
    intro w
    cases w with
    | zero =>
      rw [Csum_zero_right]
      simpa using Nat.pow_le_pow_left (by norm_num : 16 ≤ 17) (n + 1)
    | succ w =>
      rw [Csum_succ_succ]
      have h1 := ih (w + 1)
      have h2 := ih w
      calc (Csum n (w + 1) + Csum n w) * 16 ^ (n + 1)
          = 16 * (Csum n (w + 1) * 16 ^ n) + 16 * (Csum n w * 16 ^ n) := by ring
        _ ≤ 16 * (17 ^ n * 16 ^ (w + 1)) + 16 * (17 ^ n * 16 ^ w) := by gcongr
        _ = 17 ^ (n + 1) * 16 ^ (w + 1) := by ring

theorem mul_step (a x y r s A B R S : ℕ) (ih : a * x * y ≤ r * s) (step : A * B ≤ R * S) :
    a * (x * A) * (y * B) ≤ r * R * (s * S) := by
  calc a * (x * A) * (y * B) = (a * x * y) * (A * B) := by ring
    _ ≤ (r * s) * (R * S) := Nat.mul_le_mul ih step
    _ = r * R * (s * S) := by ring

theorem numeric_base : ∀ w : ℕ,
    15552 * 17 ^ (4 * (13 * w + 24)) * 16 ^ (2 * w) ≤ 2 ^ (13 * w + 24) * 16 ^ (4 * (13 * w + 24)) := by
  intro w
  induction w with
  | zero => norm_num
  | succ w ih =>
    have hA : (17 : ℕ) ^ (4 * (13 * (w + 1) + 24)) = 17 ^ (4 * (13 * w + 24)) * 17 ^ 52 := by
      rw [← pow_add]; congr 1
    have hB : (16 : ℕ) ^ (2 * (w + 1)) = 16 ^ (2 * w) * 16 ^ 2 := by
      rw [← pow_add]; congr 1
    have hC : (2 : ℕ) ^ (13 * (w + 1) + 24) = 2 ^ (13 * w + 24) * 2 ^ 13 := by
      rw [← pow_add]; congr 1
    have hD : (16 : ℕ) ^ (4 * (13 * (w + 1) + 24)) = 16 ^ (4 * (13 * w + 24)) * 16 ^ 52 := by
      rw [← pow_add]; congr 1
    have step : (17 : ℕ) ^ 52 * 16 ^ 2 ≤ 2 ^ 13 * 16 ^ 52 := by norm_num
    rw [hA, hB, hC, hD]
    exact mul_step _ _ _ _ _ _ _ _ _ ih step

theorem numeric_all (w j : ℕ) :
    15552 * 17 ^ (4 * (13 * w + 24 + j)) * 16 ^ (2 * w) ≤
      2 ^ (13 * w + 24 + j) * 16 ^ (4 * (13 * w + 24 + j)) := by
  induction j with
  | zero => simpa using numeric_base w
  | succ j ih =>
    have hA : (17 : ℕ) ^ (4 * (13 * w + 24 + (j + 1))) = 17 ^ (4 * (13 * w + 24 + j)) * 17 ^ 4 := by
      rw [← pow_add]; congr 1
    have hC : (2 : ℕ) ^ (13 * w + 24 + (j + 1)) = 2 ^ (13 * w + 24 + j) * 2 ^ 1 := by
      rw [← pow_add]; congr 1
    have hD : (16 : ℕ) ^ (4 * (13 * w + 24 + (j + 1))) = 16 ^ (4 * (13 * w + 24 + j)) * 16 ^ 4 := by
      rw [← pow_add]; congr 1
    have step : (17 : ℕ) ^ 4 * 1 ≤ 2 ^ 1 * 16 ^ 4 := by norm_num
    have := mul_step _ _ _ _ _ _ _ _ _ ih step
    rw [hA, hC, hD]
    simpa using this

/-- **Lemma 4.1, linear bound**: the inequality `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2` fails for every
`ℓ > 13w + 23`. -/
theorem lambdaBound_linear (w : ℕ) : LambdaBound w (13 * w + 23) := by
  intro ℓ hℓ
  obtain ⟨j, rfl⟩ : ∃ j, ℓ = 13 * w + 24 + j := ⟨ℓ - (13 * w + 24), by omega⟩
  set ℓ := 13 * w + 24 + j with hdef
  have ch := Csum_chernoff (2 * ℓ) w
  have sq : (Csum (2 * ℓ) w * 16 ^ (2 * ℓ)) ^ 2 ≤ (17 ^ (2 * ℓ) * 16 ^ w) ^ 2 :=
    Nat.pow_le_pow_left ch 2
  have num := numeric_all w j
  rw [← hdef] at num
  have pos : 0 < 16 ^ (4 * ℓ) := by positivity
  have hsplit : 2 ^ ℓ = 2 * 2 ^ (ℓ - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  have key : 2 * (7776 * Csum (2 * ℓ) w ^ 2) * 16 ^ (4 * ℓ) ≤ 2 * 2 ^ (ℓ - 1) * 16 ^ (4 * ℓ) := by
    calc 2 * (7776 * Csum (2 * ℓ) w ^ 2) * 16 ^ (4 * ℓ)
        = 15552 * (Csum (2 * ℓ) w * 16 ^ (2 * ℓ)) ^ 2 := by
          rw [show 4 * ℓ = 2 * ℓ + 2 * ℓ by ring, pow_add]; ring
      _ ≤ 15552 * (17 ^ (2 * ℓ) * 16 ^ w) ^ 2 := Nat.mul_le_mul_left _ sq
      _ = 15552 * 17 ^ (4 * ℓ) * 16 ^ (2 * w) := by
          rw [show 4 * ℓ = 2 * ℓ + 2 * ℓ by ring, show 2 * w = w + w by ring, pow_add, pow_add]
          ring
      _ ≤ 2 ^ ℓ * 16 ^ (4 * ℓ) := num
      _ = 2 * 2 ^ (ℓ - 1) * 16 ^ (4 * ℓ) := by rw [hsplit]
  have := Nat.le_of_mul_le_mul_right key pos
  omega

/-- `Λ(w)`: the largest `ℓ` with `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2` (it exists by
`lambdaBound_linear`; see `Lam_spec`, `Lam_greatest`). -/
def Lam (w : ℕ) : ℕ :=
  Nat.findGreatest (fun ℓ => 2 ^ (ℓ - 1) < 7776 * Csum (2 * ℓ) w ^ 2) (13 * w + 23)

theorem Lam_le (w : ℕ) : Lam w ≤ 13 * w + 23 := Nat.findGreatest_le _

theorem Lam_pos (w : ℕ) : 1 ≤ Lam w := by
  apply Nat.le_findGreatest (by omega)
  have := Csum_pos (2 * 1) w
  have : 1 ≤ Csum (2 * 1) w ^ 2 := Nat.one_le_pow _ _ this
  simp only [Nat.sub_self, pow_zero]
  omega

/-- `Λ(w)` satisfies the defining inequality. -/
theorem Lam_spec (w : ℕ) : 2 ^ (Lam w - 1) < 7776 * Csum (2 * Lam w) w ^ 2 := by
  have h := Nat.findGreatest_spec (P := fun ℓ => 2 ^ (ℓ - 1) < 7776 * Csum (2 * ℓ) w ^ 2)
    (m := 1) (n := 13 * w + 23) (by omega) (by
      have := Csum_pos (2 * 1) w
      have : 1 ≤ Csum (2 * 1) w ^ 2 := Nat.one_le_pow _ _ this
      simp only [Nat.sub_self, pow_zero]
      omega)
  exact h

/-- The inequality fails above `Λ(w)`. -/
theorem lambdaBound_Lam (w : ℕ) : LambdaBound w (Lam w) := by
  intro ℓ hℓ
  by_cases hle : ℓ ≤ 13 * w + 23
  · have := Nat.findGreatest_is_greatest (P := fun ℓ => 2 ^ (ℓ - 1) < 7776 * Csum (2 * ℓ) w ^ 2)
      hℓ hle
    exact not_lt.mp this
  · exact lambdaBound_linear w ℓ (by omega)

/-- `Λ(w)` is the largest solution of the inequality. -/
theorem Lam_greatest (w ℓ : ℕ) (h : 2 ^ (ℓ - 1) < 7776 * Csum (2 * ℓ) w ^ 2) : ℓ ≤ Lam w := by
  by_contra hlt
  exact absurd (lambdaBound_Lam w ℓ (by omega)) (not_le.mpr h)

/-! ## Digits: no straddling (Lemma 4.4) -/

theorem wt_two_pow_sub_one (n : ℕ) : wt (2 ^ n - 1) = n := by
  induction n with
  | zero => simp
  | succ n ih =>
    have hpos : 0 < 2 ^ n := Nat.two_pow_pos n
    have e : 2 ^ (n + 1) - 1 = 2 * (2 ^ n - 1) + 1 := by rw [pow_succ]; omega
    rw [e, wt_two_mul_add _ 1 (by norm_num), ih]
    omega

/-- A number in the top `2^s` of `[0, 2^k)` has its bits `s, …, k-1` equal to one. -/
theorem wt_ge_of_top (k s a : ℕ) (hs : s ≤ k) (ha : a < 2 ^ k) (htop : 2 ^ k ≤ a + 2 ^ s) :
    k - s ≤ wt a := by
  have hpos : 0 < 2 ^ s := Nat.two_pow_pos s
  have hsplit : 2 ^ k = 2 ^ s * 2 ^ (k - s) := by
    rw [← pow_add]
    congr 1
    omega
  have hdecomp : a = a % 2 ^ s + 2 ^ s * (a / 2 ^ s) := (Nat.mod_add_div a (2 ^ s)).symm
  have hr : a % 2 ^ s < 2 ^ s := Nat.mod_lt _ hpos
  have hq : a / 2 ^ s = 2 ^ (k - s) - 1 := by
    have hposq : 0 < 2 ^ (k - s) := Nat.two_pow_pos _
    have up : a / 2 ^ s < 2 ^ (k - s) := by
      rw [Nat.div_lt_iff_lt_mul hpos, mul_comm, ← hsplit]
      exact ha
    have low : 2 ^ (k - s) ≤ a / 2 ^ s + 1 := by
      by_contra hlt
      push_neg at hlt
      have : 2 ^ s * (a / 2 ^ s + 2) ≤ 2 ^ s * 2 ^ (k - s) := Nat.mul_le_mul_left _ (by omega)
      nlinarith
    omega
  rw [hdecomp, wt_add_mul_two_pow s _ _ hr, hq, wt_two_pow_sub_one]
  omega

/-- **Lemma 4.4, first case.** If `a' < N - 2^j ≤ a' + σ` with `σ < 2^(lam-1) ≤ 2^j`, then
`wt a' ≥ k - lam`. -/
theorem wt_ge_of_close (k lam j σ a' : ℕ) (h1 : 1 ≤ lam) (hlj : lam - 1 ≤ j) (hjk : j < k)
    (hσ : σ < 2 ^ (lam - 1)) (hlt : a' + 2 ^ j < 2 ^ k) (hge : 2 ^ k ≤ a' + 2 ^ j + σ) :
    k - lam ≤ wt a' := by
  have hmono : 2 ^ (lam - 1) ≤ 2 ^ j := Nat.pow_le_pow_right (by norm_num) hlj
  have top := wt_ge_of_top k (lam - 1) (a' + 2 ^ j) (by omega) hlt (by omega)
  have := wt_add_two_pow_le j a'
  omega

/-! ## The window (Lemma 4.2) -/

/-- Position of `a` relative to `a₀` when `2(a - a₀)` is within `2σ` of a multiple of `p`. -/
theorem window_cases (N p σ a a₀ : ℕ) (hp : p = 0 ∨ 2 * N + 1 ≤ p) (ha₀ : a₀ < N)
    (h : ∃ z : ℤ, |2 * ((a : ℤ) - a₀) - p * z| ≤ 2 * σ) :
    (a₀ ≤ a + σ ∧ a ≤ a₀ + σ) ∨ a < σ ∨ N < a + σ := by
  obtain ⟨z, hz⟩ := h
  rw [abs_le] at hz
  obtain ⟨hz1, hz2⟩ := hz
  rcases hp with rfl | hp
  · left
    simp only [Nat.cast_zero, zero_mul, sub_zero] at hz1 hz2
    omega
  · have hpz : (2 * (N : ℤ) + 1) ≤ p := by exact_mod_cast hp
    rcases lt_trichotomy z 0 with hneg | rfl | hpos
    · right; left
      have : (p : ℤ) * z ≤ -(p : ℤ) := by nlinarith
      omega
    · left
      simp only [mul_zero, sub_zero] at hz1 hz2
      omega
    · right; right
      have : (p : ℤ) ≤ (p : ℤ) * z := by nlinarith
      omega

/-- The block indices that `a` can have (Lemma 4.2: at most six aligned blocks; five here). -/
def blocks (N σ L a₀ : ℕ) : Finset ℕ :=
  {(a₀ - σ) / L, (a₀ - σ) / L + 1, 0, (N - σ) / L, (N - σ) / L + 1}

theorem card_blocks_le (N σ L a₀ : ℕ) : (blocks N σ L a₀).card ≤ 6 :=
  Finset.card_le_five.trans (by norm_num)

theorem div_mem_pair (L lo a : ℕ) (hL : 0 < L) (h1 : lo ≤ a) (h2 : a < lo + L) :
    a / L = lo / L ∨ a / L = lo / L + 1 := by
  have m1 : lo / L ≤ a / L := Nat.div_le_div_right h1
  have m2 : a / L ≤ (lo + L) / L := Nat.div_le_div_right h2.le
  rw [Nat.add_div_right _ hL] at m2
  omega

/-- **Lemma 4.2 (window).** -/
theorem window_index (N p σ ℓ a a₀ : ℕ) (hp : p = 0 ∨ 2 * N + 1 ≤ p) (hℓ : 2 * σ + 1 ≤ 2 ^ ℓ)
    (ha : a < N) (ha₀ : a₀ < N) (h : ∃ z : ℤ, |2 * ((a : ℤ) - a₀) - p * z| ≤ 2 * σ) :
    a / 2 ^ ℓ ∈ blocks N σ (2 ^ ℓ) a₀ := by
  have hL : 0 < 2 ^ ℓ := Nat.two_pow_pos ℓ
  simp only [blocks, Finset.mem_insert, Finset.mem_singleton]
  rcases window_cases N p σ a a₀ hp ha₀ h with ⟨h1, h2⟩ | h1 | h1
  · rcases div_mem_pair (2 ^ ℓ) (a₀ - σ) a hL (by omega) (by omega) with e | e
    · exact Or.inl e
    · exact Or.inr (Or.inl e)
  · right; right; left
    exact Nat.div_eq_of_lt (by omega)
  · rcases div_mem_pair (2 ^ ℓ) (N - σ) a hL (by omega) (by omega) with e | e
    · exact Or.inr (Or.inr (Or.inr (Or.inl e)))
    · exact Or.inr (Or.inr (Or.inr (Or.inr e)))

/-! ## The count of Lemma 4.3 -/

/-- A set of exponent pairs in `[0,N)²` of weight `≤ w`, all within the window of one pair,
has at most `36·C(2ℓ, ≤w)` elements. -/
theorem card_le_of_window (N p σ ℓ w : ℕ) (hp : p = 0 ∨ 2 * N + 1 ≤ p) (hℓ : 2 * σ + 1 ≤ 2 ^ ℓ)
    (J : Finset (ℕ × ℕ)) (a₀ b₀ : ℕ) (ha₀ : a₀ < N) (hb₀ : b₀ < N)
    (hJ : ∀ q ∈ J, q.1 < N ∧ q.2 < N ∧ wt q.1 + wt q.2 ≤ w ∧
      (∃ z : ℤ, |2 * ((q.1 : ℤ) - a₀) - p * z| ≤ 2 * σ) ∧
      (∃ z : ℤ, |2 * ((q.2 : ℤ) - b₀) - p * z| ≤ 2 * σ)) :
    J.card ≤ 36 * Csum (2 * ℓ) w := by
  have hL : 0 < 2 ^ ℓ := Nat.two_pow_pos ℓ
  let T : Finset ((ℕ × ℕ) × (ℕ × ℕ)) :=
    (blocks N σ (2 ^ ℓ) a₀ ×ˢ blocks N σ (2 ^ ℓ) b₀) ×ˢ lowPairs ℓ w
  have hcard : J.card ≤ T.card := by
    apply Finset.card_le_card_of_injOn
      (fun q => ((q.1 / 2 ^ ℓ, q.2 / 2 ^ ℓ), (q.1 % 2 ^ ℓ, q.2 % 2 ^ ℓ)))
    · intro q hq
      obtain ⟨h1, h2, hw, hz1, hz2⟩ := hJ q hq
      have d1 : q.1 = q.1 % 2 ^ ℓ + 2 ^ ℓ * (q.1 / 2 ^ ℓ) := (Nat.mod_add_div _ _).symm
      have d2 : q.2 = q.2 % 2 ^ ℓ + 2 ^ ℓ * (q.2 / 2 ^ ℓ) := (Nat.mod_add_div _ _).symm
      have w1 := wt_add_mul_two_pow ℓ (q.1 % 2 ^ ℓ) (q.1 / 2 ^ ℓ) (Nat.mod_lt _ hL)
      have w2 := wt_add_mul_two_pow ℓ (q.2 % 2 ^ ℓ) (q.2 / 2 ^ ℓ) (Nat.mod_lt _ hL)
      rw [← d1] at w1
      rw [← d2] at w2
      simp only [T, lowPairs, Finset.mem_product, Finset.mem_filter, Finset.mem_range]
      exact ⟨⟨window_index N p σ ℓ q.1 a₀ hp hℓ h1 ha₀ hz1,
        window_index N p σ ℓ q.2 b₀ hp hℓ h2 hb₀ hz2⟩,
        ⟨Nat.mod_lt _ hL, Nat.mod_lt _ hL⟩, by omega⟩
    · intro q _ q' _ h
      simp only [Prod.mk.injEq] at h
      obtain ⟨⟨e1, e2⟩, e3, e4⟩ := h
      have d1 := Nat.mod_add_div q.1 (2 ^ ℓ)
      have d2 := Nat.mod_add_div q.2 (2 ^ ℓ)
      have d1' := Nat.mod_add_div q'.1 (2 ^ ℓ)
      have d2' := Nat.mod_add_div q'.2 (2 ^ ℓ)
      rw [e1, e3] at d1
      rw [e2, e4] at d2
      exact Prod.ext (d1.symm.trans d1') (d2.symm.trans d2')
  have hT : T.card ≤ 36 * Csum (2 * ℓ) w := by
    simp only [T, Finset.card_product]
    have c1 := card_blocks_le N σ (2 ^ ℓ) a₀
    have c2 := card_blocks_le N σ (2 ^ ℓ) b₀
    have c3 := card_lowPairs_le ℓ w
    calc (blocks N σ (2 ^ ℓ) a₀).card * (blocks N σ (2 ^ ℓ) b₀).card * (lowPairs ℓ w).card
        ≤ 6 * 6 * Csum (2 * ℓ) w := Nat.mul_le_mul (Nat.mul_le_mul c1 c2) c3
      _ = 36 * Csum (2 * ℓ) w := by ring
  exact hcard.trans hT

/-- The arithmetic of Lemma 4.3: if a set of `m` elements fits `36·C(2ℓ, ≤w)` for every `ℓ`
with `2^ℓ ≥ 2σ + 1`, `σ = 3(m-1)(m-2)`, then `σ < 2^(lam-1)`. -/
theorem sigma_lt_of_lambda (w lam m : ℕ) (hlam : LambdaBound w lam) (h1 : 1 ≤ lam)
    (hm : ∀ ℓ, 2 * (3 * (m - 1) * (m - 2)) + 1 ≤ 2 ^ ℓ → m ≤ 36 * Csum (2 * ℓ) w) :
    3 * (m - 1) * (m - 2) < 2 ^ (lam - 1) := by
  by_contra hcon
  push_neg at hcon
  set σ := 3 * (m - 1) * (m - 2) with hσ
  have hsplit : 2 ^ lam = 2 * 2 ^ (lam - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  have hpos : 0 < 2 ^ (lam - 1) := Nat.two_pow_pos _
  set ℓ := Nat.clog 2 (2 * σ + 1) with hℓdef
  have hle : 2 * σ + 1 ≤ 2 ^ ℓ := Nat.le_pow_clog (by norm_num) _
  have hlt : lam < ℓ := by
    by_contra hnot
    push_neg at hnot
    have : 2 ^ ℓ ≤ 2 ^ lam := Nat.pow_le_pow_right (by norm_num) hnot
    omega
  have hpred : 2 ^ (ℓ - 1) < 2 * σ + 1 := by
    have := Nat.pow_pred_clog_lt_self (b := 2) (by norm_num) (x := 2 * σ + 1) (by omega)
    simpa [Nat.pred_eq_sub_one] using this
  have hcard := hm ℓ hle
  have hbound := hlam ℓ hlt
  have hm3 : 1 ≤ m := by
    by_contra h0
    have : m = 0 := by omega
    rw [this] at hσ
    simp at hσ
    omega
  have h6 : 2 * σ + 1 ≤ 6 * m ^ 2 := by
    obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
    rw [hσ]
    cases n with
    | zero => norm_num
    | succ n =>
      simp only [Nat.add_sub_cancel, show n + 1 + 1 - 2 = n by omega]
      nlinarith
  have hsq : m ^ 2 ≤ (36 * Csum (2 * ℓ) w) ^ 2 := Nat.pow_le_pow_left hcard 2
  have : 6 * m ^ 2 ≤ 7776 * Csum (2 * ℓ) w ^ 2 := by
    calc 6 * m ^ 2 ≤ 6 * (36 * Csum (2 * ℓ) w) ^ 2 := Nat.mul_le_mul_left _ hsq
      _ = 7776 * Csum (2 * ℓ) w ^ 2 := by ring
  omega

/-- The dichotomy of Lemma 4.3 for a pair `a' < a < N`: if `2(a - a')` is within `2σ` of a
multiple of `p` (`p = 0` or `p ≥ 2N + 1`, `2σ < N`), then `a - a' ≤ σ` or `a > N - σ`. -/
theorem dichotomy (N p σ a a' : ℕ) (hp : p = 0 ∨ 2 * N + 1 ≤ p) (haa : a' < a) (ha : a < N)
    (hσ : 2 * σ < N) (hz : ∃ z : ℤ, |2 * ((a : ℤ) - a') - p * z| ≤ 2 * σ) :
    a ≤ a' + σ ∨ N < a + σ := by
  obtain ⟨z, hz⟩ := hz
  rw [abs_le] at hz
  obtain ⟨hz1, hz2⟩ := hz
  rcases hp with rfl | hp
  · left
    simp only [Nat.cast_zero, zero_mul, sub_zero] at hz1 hz2
    omega
  · have hpz : (2 * (N : ℤ) + 1) ≤ p := by exact_mod_cast hp
    rcases lt_trichotomy z 0 with hneg | rfl | hpos
    · exfalso
      have : (p : ℤ) * z ≤ (p : ℤ) * (-1) :=
        mul_le_mul_of_nonneg_left (by omega) (Nat.cast_nonneg p)
      omega
    · left
      simp only [mul_zero, sub_zero] at hz1 hz2
      omega
    · right
      have : (p : ℤ) * 1 ≤ (p : ℤ) * z :=
        mul_le_mul_of_nonneg_left (by omega) (Nat.cast_nonneg p)
      omega

/-- **Lemma 4.4 as used in Theorem 5.1**: two exponents of weight `≤ w` on the two sides of a
top-run cut `N - 2^j` cannot be within `2σ` of each other modulo `p` (in the sense of
`2(a - a')`), when `σ < 2^(lam-1)` and `k ≥ w + lam + 1`. -/
theorem straddle_contra (k lam j σ p w a a' : ℕ) (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p)
    (h1 : 1 ≤ lam) (hk : w + lam + 1 ≤ k) (hj : j < k) (hσ : σ < 2 ^ (lam - 1))
    (ha : a < 2 ^ k) (htop : 2 ^ k - 2 ^ j ≤ a) (hnot : a' < 2 ^ k - 2 ^ j)
    (hwa : wt a ≤ w) (hwa' : wt a' ≤ w)
    (hz : ∃ z : ℤ, |2 * ((a : ℤ) - a') - p * z| ≤ 2 * σ) : False := by
  have hjk : 2 * 2 ^ j ≤ 2 ^ k := by
    rw [← pow_succ']
    exact Nat.pow_le_pow_right (by norm_num) hj
  have hpos : 0 < 2 ^ j := Nat.two_pow_pos j
  have wa := wt_ge_of_top k j a hj.le ha (by omega)
  have hlj : lam - 1 ≤ j := by omega
  have hmono : 2 ^ (lam - 1) ≤ 2 ^ j := Nat.pow_le_pow_right (by norm_num) hlj
  rcases dichotomy (2 ^ k) p σ a a' hp (by omega) ha (by omega) hz with hc | hf
  · have := wt_ge_of_close k lam j σ a' h1 hlj hj hσ (by omega) (by omega)
    omega
  · have := wt_ge_of_top k (lam - 1) a (by omega) ha (by omega)
    omega

end ICFallDegrees.GrowthCounting
