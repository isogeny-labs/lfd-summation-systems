import ICFallDegrees.ExactK6Cert
import ICFallDegrees.N5CertSound
import Mathlib.Algebra.Polynomial.Inductions
import Mathlib.Algebra.Polynomial.Eval.Defs

/-!
# Soundness of the depth-6 derivation checker (`1 ∈ V_5`)

`ExactK6Cert.rowCheck` certifies one row of a derivation in the quotient closure of degree 5 of
the two-summand chain presentation at depth 6.  Here the packed rows get their meaning and the
checks are shown to imply `InQuotClosure … 5 (vec (row n))` for every certified row.

* `slot B R c = R / B^c % B` is slot `c` of a packed natural; `slot_eval` is the uniqueness of
  base-`B` digits in the form "a polynomial over `ℕ` with coefficients `< B` is read back from its
  value at `B`".  All linear algebra on packed numbers goes through this one statement: the
  accumulator of `rowCheck` is the value at `B` of an explicit polynomial `accPoly`.
* `M c` is the squarefree chain monomial of slot `c = e₀ + 128 e₁`, with block exponents
  `(e₀ mod 64, e₁ mod 64)` (the wrap is built into `ChainTranslation.digits 6`), and
  `nf_M_mul : normalForm (M s * M c) = M (c + s)` when the low parts do not carry out of the
  128-window.
* `vec B R` is the vector of a stored row (slots with `e₀, e₁ < 64`).
* `scan_sound` turns the Bool scan into slot statements, `regroup` splits the `128 × 128` slots
  into four quadrants, and `row_sound` is the induction step: all earlier rows in the quotient
  closure with their degree bounds ⟹ the new row is.

Only `propext`, `Classical.choice`, `Quot.sound` are used.
-/

noncomputable section

namespace ICFallDegrees.ExactK6Upper

open ICFallDegrees.N5Cert ICFallDegrees.ExactK6Cert

/-! ## Slots of packed naturals -/

/-- Slot `c` of a natural packed in base `B`. -/
def slot (B R c : ℕ) : ℕ := R / B ^ c % B

theorem slot_zero (B R : ℕ) : slot B R 0 = R % B := by simp [slot]

theorem slot_div (B R m c : ℕ) : slot B (R / B ^ m) c = slot B R (c + m) := by
  unfold slot
  rw [Nat.div_div_eq_div_mul, ← pow_add, add_comm]

theorem slot_div_one (B R c : ℕ) : slot B (R / B) c = slot B R (c + 1) := by
  simpa using slot_div B R 1 c

theorem slot_mod {B R m c : ℕ} (h : c < m) : slot B (R % B ^ m) c = slot B R c := by
  unfold slot
  obtain ⟨k, rfl⟩ : ∃ k, m = c + (k + 1) := ⟨m - c - 1, by omega⟩
  rw [pow_add, Nat.mod_mul_right_div_self, pow_succ, Nat.mod_mod_of_dvd _ (Dvd.intro_left _ rfl)]

theorem slot_eq_zero {B R c : ℕ} (h : R < B ^ c) : slot B R c = 0 := by
  unfold slot
  rw [Nat.div_eq_of_lt h, Nat.zero_mod]

theorem slot_eq_zero_of_le {B R m c : ℕ} (hB : 0 < B) (h : R < B ^ m) (hc : m ≤ c) :
    slot B R c = 0 :=
  slot_eq_zero (lt_of_lt_of_le h (Nat.pow_le_pow_right hB hc))

/-- The slots below `n` recompose the residue modulo `B^n`. -/
theorem sum_slot (B R : ℕ) : ∀ n, ∑ c ∈ Finset.range n, slot B R c * B ^ c = R % B ^ n
  | 0 => by simp [Nat.mod_one]
  | n+1 => by
    rw [Finset.sum_range_succ, sum_slot B R n, Nat.mod_pow_succ]
    unfold slot
    ring

/-- Uniqueness of digits: a polynomial over `ℕ` with coefficients `< B` is read back from its
value at `B`. -/
theorem slot_eval (B : ℕ) : ∀ (c : ℕ) (f : Polynomial ℕ), (∀ i, f.coeff i < B) →
    slot B (f.eval B) c = f.coeff c := by
  intro c
  induction c with
  | zero =>
    intro f hf
    have h : f.divX.eval B * B + f.coeff 0 = f.eval B := by
      simpa using congrArg (Polynomial.eval B) (Polynomial.divX_mul_X_add f)
    rw [slot_zero, ← h, Nat.mul_add_mod_self_right, Nat.mod_eq_of_lt (hf 0)]
  | succ c ih =>
    intro f hf
    have hB : 0 < B := lt_of_le_of_lt (Nat.zero_le _) (hf 0)
    have h : f.divX.eval B * B + f.coeff 0 = f.eval B := by
      simpa using congrArg (Polynomial.eval B) (Polynomial.divX_mul_X_add f)
    have hdiv : f.eval B / B = f.divX.eval B := by
      rw [← h, add_comm, Nat.add_mul_div_right _ _ hB, Nat.div_eq_of_lt (hf 0), zero_add]
    rw [← slot_div_one, hdiv, ih f.divX (fun i => by rw [Polynomial.coeff_divX]; exact hf _),
      Polynomial.coeff_divX]

/-! ## Regrouping the `128 × 128` slots into four quadrants -/

section Regroup

variable {A : Type*} [AddCommMonoid A]

theorem sum_range_mul' (g : ℕ → A) (b : ℕ) : ∀ a : ℕ,
    ∑ c ∈ Finset.range (b * a), g c = ∑ j ∈ Finset.range a, ∑ i ∈ Finset.range b, g (i + b * j)
  | 0 => by simp
  | a+1 => by
    rw [Nat.mul_succ, Finset.sum_range_add, sum_range_mul' g b a, Finset.sum_range_succ]
    congr 1
    apply Finset.sum_congr rfl
    intro i _
    rw [add_comm]

theorem sum_range_half (f : ℕ → A) :
    ∑ i ∈ Finset.range 128, f i = ∑ i ∈ Finset.range 64, (f i + f (64 + i)) := by
  rw [show (128 : ℕ) = 64 + 64 from rfl, Finset.sum_range_add, Finset.sum_add_distrib]

/-- The four quadrants of the `128 × 128` window. -/
theorem regroup (g : ℕ → A) :
    ∑ c ∈ Finset.range 16384, g c = ∑ e1 ∈ Finset.range 64, ∑ e0 ∈ Finset.range 64,
      (g (e0 + 128 * e1) + g (64 + e0 + 128 * e1) + g (e0 + 128 * (64 + e1)) +
        g (64 + e0 + 128 * (64 + e1))) := by
  rw [show (16384 : ℕ) = 128 * 128 from rfl, sum_range_mul', sum_range_half]
  apply Finset.sum_congr rfl
  intro e1 _
  rw [sum_range_half, sum_range_half, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro e0 _
  abel

end Regroup

/-! ## Bits and chain exponents modulo `2^k` -/

open MvPolynomial FiniteQuotient ChainQuotient ChainTranslation ChainClosureBridge BinaryWeight
open ICFallDegrees.N5CertSound

theorem bit_mod {n j k : ℕ} (h : j < k) : bit (n % 2 ^ k) j = bit n j := by
  unfold bit
  obtain ⟨r, rfl⟩ : ∃ r, k = j + (r + 1) := ⟨k - j - 1, by omega⟩
  rw [pow_add, Nat.mod_mul_right_div_self, pow_succ, Nat.mod_mod_of_dvd _ (Dvd.intro_left _ rfl)]

/-- The chain exponent only depends on the exponent vector modulo `2^k`. -/
theorem digits_congr {m k : ℕ} {b b' : Fin m → ℕ} (h : ∀ i, b i % 2 ^ k = b' i % 2 ^ k) :
    digits k b = digits k b' := by
  ext v
  rw [digits_apply, digits_apply, ← bit_mod (depth v).isLt, h, bit_mod (depth v).isLt]

theorem digits_zero : digits (m := 2) 6 (0 : Fin 2 → ℕ) = 0 := by
  ext v
  rw [digits_apply]
  simp [bit]

variable {K : Type*} [Field K]

section Translation

variable (two_ne : (2 : K) ≠ 0)
include two_ne

/-- A power of a level-0 chain variable times a squarefree monomial. -/
theorem nf_X_pow_mul (i : Fin 2) : ∀ (n : ℕ) (b : Fin 2 → ℕ),
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6)
        (X (var 2 6 i ⟨0, by norm_num⟩) ^ n * monomial (digits 6 b) (1 : K)) =
      monomial (digits 6 (Function.update b i (b i + n))) 1
  | 0, b => by
    rw [pow_zero, one_mul, add_zero, Function.update_eq_self]
    exact normalForm_of_reduced _ _ (reduced_chain_ideal_trivial 2 6 two_ne) (digits_reduced _ _)
  | n+1, b => by
    rw [pow_succ', mul_assoc, normalForm_mul_normalForm _ _ (reduced_chain_ideal_trivial 2 6 two_ne),
      nf_X_pow_mul i n b]
    have hbox : ∀ i', (fun i' => Function.update b i (b i + n) i' % 2 ^ 6) i' < 2 ^ 6 :=
      fun i' => Nat.mod_lt _ (by norm_num)
    rw [digits_congr (b' := fun i' => Function.update b i (b i + n) i' % 2 ^ 6)
        (fun i' => (Nat.mod_mod _ _).symm),
      normalForm_X_mul_monomial two_ne i ⟨0, by norm_num⟩ _ hbox]
    congr 2
    apply digits_congr
    intro i'
    by_cases hi : i' = i
    · subst hi
      simp only [Function.update_self, pow_zero]
      omega
    · simp only [Function.update_of_ne hi, Nat.mod_mod]

/-- The monomial `u_{0,0}^{b₀} u_{1,0}^{b₁}`. -/
def pw (b : Fin 2 → ℕ) : P (2 * 6) K :=
  X (var 2 6 0 ⟨0, by norm_num⟩) ^ b 0 * X (var 2 6 1 ⟨0, by norm_num⟩) ^ b 1

theorem nf_pw (b : Fin 2 → ℕ) :
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6) (pw (K := K) b) =
      monomial (digits 6 b) 1 := by
  have h1 : pw (K := K) b = X (var 2 6 0 ⟨0, by norm_num⟩) ^ b 0 *
      (X (var 2 6 1 ⟨0, by norm_num⟩) ^ b 1 * monomial (digits 6 (0 : Fin 2 → ℕ)) (1 : K)) := by
    rw [digits_zero]
    show _ = _ * (_ * (1 : P (2 * 6) K))
    rw [mul_one]
    rfl
  rw [h1, normalForm_mul_normalForm _ _ (reduced_chain_ideal_trivial 2 6 two_ne),
    nf_X_pow_mul two_ne 1 (b 1) 0, nf_X_pow_mul two_ne 0 (b 0)]
  have hb : Function.update (Function.update (0 : Fin 2 → ℕ) 1 ((0 : Fin 2 → ℕ) 1 + b 1)) 0
      (Function.update (0 : Fin 2 → ℕ) 1 ((0 : Fin 2 → ℕ) 1 + b 1) 0 + b 0) = b := by
    funext i
    fin_cases i <;> simp
  rw [hb]

/-- The product of two squarefree monomials reduces to the monomial of the sum of the exponent
vectors (modulo `64`, inside `digits 6`). -/
theorem nf_digits_mul (b b' : Fin 2 → ℕ) :
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6)
        (monomial (digits 6 b) (1 : K) * monomial (digits 6 b') 1) =
      monomial (digits 6 (b + b')) 1 := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
  have e1 : pw (K := K) b' * pw b = pw (b + b') := by
    unfold pw
    simp only [Pi.add_apply, pow_add]
    ring
  rw [← nf_pw two_ne b, ← nf_pw two_ne b', ← normalForm_mul_normalForm _ _ zero,
    mul_comm (normalForm _ _ (pw b)) (pw b'), ← normalForm_mul_normalForm _ _ zero, e1,
    nf_pw two_ne (b + b')]

/-! ## The monomial of a slot -/

/-- The squarefree chain monomial of slot `c = e₀ + 128 e₁`: block exponents
`(e₀ mod 64, e₁ mod 64)`. -/
def M (c : ℕ) : P (2 * 6) K := monomial (digits 6 (col c)) 1

omit two_ne in
theorem col_add {c s : ℕ} (h : c % 128 + s % 128 < 128) : col (c + s) = col c + col s := by
  funext i
  fin_cases i
  · show (c + s) % 128 = c % 128 + s % 128
    omega
  · show (c + s) / 128 = c / 128 + s / 128
    omega

theorem nf_M_mul {s c : ℕ} (h : c % 128 + s % 128 < 128) :
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6) (M (K := K) s * M c) = M (c + s) := by
  unfold M
  rw [nf_digits_mul two_ne, col_add h, add_comm (col s)]

omit two_ne in
theorem M_congr {c c' : ℕ} (h0 : c % 128 % 64 = c' % 128 % 64)
    (h1 : c / 128 % 64 = c' / 128 % 64) : M (K := K) c = M c' := by
  unfold M
  have e : (2 : ℕ) ^ 6 = 64 := by norm_num
  rw [digits_congr (k := 6) (b := col c) (b' := col c')]
  intro i
  fin_cases i
  · show c % 128 % 2 ^ 6 = c' % 128 % 2 ^ 6
    rw [e]
    exact h0
  · show c / 128 % 2 ^ 6 = c' / 128 % 2 ^ 6
    rw [e]
    exact h1

omit two_ne in
theorem M_degree {c : ℕ} (h0 : c % 128 < 64) (h1 : c / 128 < 64) :
    (M (K := K) c).totalDegree ≤ popcount 6 (c % 128) + popcount 6 (c / 128) := by
  unfold M
  refine (totalDegree_monomial_le _ _).trans ?_
  have e : (2 : ℕ) ^ 6 = 64 := by norm_num
  have hd := ChainModelMap.digits_degree (k := 6) (col c) (by
    intro i
    fin_cases i
    · show c % 128 < 2 ^ 6
      omega
    · show c / 128 < 2 ^ 6
      omega)
  have hrfl : ((digits 6 (col c)).sum fun _ => id) = (digits 6 (col c)).degree := rfl
  rw [hrfl, hd, Fin.sum_univ_two, col_zero, col_one, popcount_eq_wt 6 _ (by omega),
    popcount_eq_wt 6 _ (by omega)]

omit two_ne in
theorem M_zero : M (K := K) 0 = 1 := by
  unfold M
  have : col 0 = 0 := (col_eq_zero_iff 0).mpr rfl
  rw [this, digits_zero]
  rfl

end Translation

/-! ## Vectors of packed rows -/

/-- The vector of a stored row: the slots with `e₀, e₁ < 64`. -/
def vec (B R : ℕ) : P (2 * 6) K :=
  ∑ e1 ∈ Finset.range 64, ∑ e0 ∈ Finset.range 64,
    (slot B R (e0 + 128 * e1) : K) • M (K := K) (e0 + 128 * e1)

/-- The digit polynomial of a stored row. -/
def tpoly (B R : ℕ) : Polynomial ℕ :=
  ∑ e1 ∈ Finset.range 64, ∑ e0 ∈ Finset.range 64,
    Polynomial.monomial (e0 + 128 * e1) (slot B R (e0 + 128 * e1))

variable (K) in
/-- The vector of a polynomial over `ℕ`: `X^c ↦ M c` on the window `c < 128 · 128`. -/
def Λ (f : Polynomial ℕ) : P (2 * 6) K :=
  ∑ c ∈ Finset.range 16384, (f.coeff c : K) • M (K := K) c

theorem Λ_apply (f : Polynomial ℕ) :
    Λ K f = ∑ c ∈ Finset.range 16384, (f.coeff c : K) • M (K := K) c := rfl

theorem Λ_zero : Λ K 0 = 0 := by
  rw [Λ_apply]
  simp

theorem Λ_add (f g : Polynomial ℕ) : Λ K (f + g) = Λ K f + Λ K g := by
  rw [Λ_apply, Λ_apply, Λ_apply, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro c _
  rw [Polynomial.coeff_add, Nat.cast_add, add_smul]

theorem Λ_sum {ι : Type*} (s : Finset ι) (g : ι → Polynomial ℕ) :
    Λ K (∑ i ∈ s, g i) = ∑ i ∈ s, Λ K (g i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [Λ_zero]
  | @insert a s has ih => rw [Finset.sum_insert has, Finset.sum_insert has, Λ_add, ih]

theorem Λ_monomial {c : ℕ} (hc : c < 16384) (a : ℕ) :
    Λ K (Polynomial.monomial c a) = (a : K) • M c := by
  rw [Λ_apply, Finset.sum_eq_single c]
  · rw [Polynomial.coeff_monomial, if_pos rfl]
  · intro b _ hb
    rw [Polynomial.coeff_monomial, if_neg (Ne.symm hb)]
    simp
  · intro h
    exact absurd (Finset.mem_range.mpr hc) h

theorem Λ_tpoly (B R : ℕ) : Λ K (tpoly B R) = vec B R := by
  unfold tpoly vec
  rw [Λ_sum]
  apply Finset.sum_congr rfl
  intro e1 h1
  rw [Λ_sum]
  apply Finset.sum_congr rfl
  intro e0 h0
  rw [Finset.mem_range] at h0 h1
  exact Λ_monomial (by omega) _

/-- What the scan establishes about a stored row. -/
structure RowFacts (B p d H : ℕ) : Prop where
  lt : H < B ^ 16384
  off1 : ∀ e1 < 64, ∀ e0 < 64, slot B H (64 + e0 + 128 * e1) = 0
  off2 : ∀ e1 < 64, ∀ e0 < 64, slot B H (e0 + 128 * (64 + e1)) = 0
  off3 : ∀ e1 < 64, ∀ e0 < 64, slot B H (64 + e0 + 128 * (64 + e1)) = 0
  ltp : ∀ e1 < 64, ∀ e0 < 64, slot B H (e0 + 128 * e1) < p
  wt : ∀ e1 < 64, ∀ e0 < 64, slot B H (e0 + 128 * e1) ≠ 0 → popcount 6 e0 + popcount 6 e1 ≤ d

/-- What the scan establishes about the accumulator. -/
structure AccFacts (B p acc H : ℕ) : Prop where
  lt : acc < B ^ 16384
  fold : ∀ e1 < 64, ∀ e0 < 64, (slot B acc (e0 + 128 * e1) + slot B acc (64 + e0 + 128 * e1) +
      slot B acc (e0 + 128 * (64 + e1)) + slot B acc (64 + e0 + 128 * (64 + e1)) +
      (p - 1) * slot B H (e0 + 128 * e1)) % p = 0

section RowFacts

variable {B p d H : ℕ}

/-- A sum over all slots of a stored row only sees the first quadrant. -/
theorem quad {A : Type*} [AddCommMonoid A] (h : RowFacts B p d H) (φ : ℕ → ℕ → A)
    (hφ : ∀ c, φ c 0 = 0) :
    ∑ c ∈ Finset.range 16384, φ c (slot B H c) =
      ∑ e1 ∈ Finset.range 64, ∑ e0 ∈ Finset.range 64, φ (e0 + 128 * e1) (slot B H (e0 + 128 * e1)) := by
  rw [regroup]
  apply Finset.sum_congr rfl
  intro e1 h1
  apply Finset.sum_congr rfl
  intro e0 h0
  rw [Finset.mem_range] at h0 h1
  rw [h.off1 e1 h1 e0 h0, h.off2 e1 h1 e0 h0, h.off3 e1 h1 e0 h0, hφ, hφ, hφ, add_zero, add_zero,
    add_zero]

theorem eval_tpoly (h : RowFacts B p d H) : (tpoly B H).eval B = H := by
  have hs := sum_slot B H 16384
  rw [Nat.mod_eq_of_lt h.lt, quad h (fun c x => x * B ^ c) (fun c => zero_mul _)] at hs
  unfold tpoly
  rw [Polynomial.eval_finset_sum]
  simp only [Polynomial.eval_finset_sum, Polynomial.eval_monomial]
  exact hs

theorem coeff_tpoly_lt (h : RowFacts B p d H) (hp : 0 < p) (c : ℕ) : (tpoly B H).coeff c < p := by
  unfold tpoly
  simp only [Polynomial.finset_sum_coeff, Polynomial.coeff_monomial]
  by_cases hc : c % 128 < 64 ∧ c / 128 < 64
  · rw [Finset.sum_eq_single (c / 128), Finset.sum_eq_single (c % 128)]
    · rw [if_pos (by omega)]
      exact h.ltp _ hc.2 _ hc.1
    · intro b _ hb
      rw [if_neg (by omega)]
    · intro hn
      exact absurd (Finset.mem_range.mpr hc.1) hn
    · intro b hb hne
      apply Finset.sum_eq_zero
      intro e0 he0
      rw [Finset.mem_range] at hb he0
      rw [if_neg (by omega)]
    · intro hn
      exact absurd (Finset.mem_range.mpr hc.2) hn
  · rw [Finset.sum_eq_zero]
    · exact hp
    · intro e1 he1
      apply Finset.sum_eq_zero
      intro e0 he0
      rw [Finset.mem_range] at he1 he0
      rw [if_neg (by omega)]

theorem vec_eq_sum {K : Type*} [Field K] (h : RowFacts B p d H) :
    vec (K := K) B H = ∑ c ∈ Finset.range 16384, (slot B H c : K) • M (K := K) c := by
  rw [quad h (fun c x => (x : K) • M (K := K) c) (fun c => by simp)]
  rfl

theorem vec_degree {K : Type*} [Field K] (h : RowFacts B p d H) :
    (vec (K := K) B H).totalDegree ≤ d := by
  unfold vec
  refine Finset.sum_induction _ (fun q : P (2 * 6) K => q.totalDegree ≤ d)
    (fun a b ha hb => (totalDegree_add a b).trans (max_le ha hb)) (by simp) ?_
  intro e1 h1
  refine Finset.sum_induction _ (fun q : P (2 * 6) K => q.totalDegree ≤ d)
    (fun a b ha hb => (totalDegree_add a b).trans (max_le ha hb)) (by simp) ?_
  intro e0 h0
  rw [Finset.mem_range] at h0 h1
  by_cases hz : slot B H (e0 + 128 * e1) = 0
  · rw [hz]
    simp
  · refine (totalDegree_smul_le _ _).trans ((M_degree (by omega) (by omega)).trans ?_)
    have := h.wt e1 h1 e0 h0 hz
    rw [show (e0 + 128 * e1) % 128 = e0 by omega, show (e0 + 128 * e1) / 128 = e1 by omega]
    exact this

end RowFacts

/-! ## Soundness of the scan -/

theorem scanSlots_sound (B p d w1 : ℕ) : ∀ (n e0 x1 x2 x3 x4 h : ℕ),
    scanSlots B p d w1 n e0 x1 x2 x3 x4 h = true →
    h / B ^ n = 0 ∧ ∀ i < n,
      (slot B x1 i + slot B x2 i + slot B x3 i + slot B x4 i + (p - 1) * slot B h i) % p = 0 ∧
      slot B h i < p ∧ (slot B h i ≠ 0 → popcount 6 (e0 + i) + w1 ≤ d)
  | 0, e0, x1, x2, x3, x4, h, hs => by
    simp only [scanSlots, beq_iff_eq] at hs
    subst hs
    exact ⟨by simp, fun i hi => absurd hi (Nat.not_lt_zero _)⟩
  | n+1, e0, x1, x2, x3, x4, h, hs => by
    simp only [scanSlots, Bool.and_eq_true, Bool.or_eq_true, beq_iff_eq, decide_eq_true_eq] at hs
    obtain ⟨⟨⟨h1, h2⟩, h3⟩, h4⟩ := hs
    obtain ⟨r0, r⟩ := scanSlots_sound B p d w1 n (e0 + 1) _ _ _ _ _ h4
    refine ⟨by rw [pow_succ', ← Nat.div_div_eq_div_mul]; exact r0, ?_⟩
    intro i hi
    cases i with
    | zero =>
      simp only [slot_zero, add_zero]
      exact ⟨h1, h2, fun hne => h3.resolve_left hne⟩
    | succ i =>
      have := r i (by omega)
      simp only [slot_div_one] at this
      rw [show e0 + 1 + i = e0 + (i + 1) by omega] at this
      exact this

/-- Block `j` (128 slots) of a packed natural. -/
def blk (B X j : ℕ) : ℕ := X / (B ^ 128) ^ j % B ^ 128

theorem scanBlocks_sound (B p d : ℕ) : ∀ (n e1 lo hi H : ℕ),
    scanBlocks B p d n e1 lo hi H = true →
    hi / (B ^ 128) ^ n = 0 ∧ H / (B ^ 128) ^ n = 0 ∧ ∀ j < n,
      blk B H j / B ^ 64 = 0 ∧ ∀ i < 64,
        (slot B (blk B lo j % B ^ 64) i + slot B (blk B lo j / B ^ 64) i +
          slot B (blk B hi j % B ^ 64) i + slot B (blk B hi j / B ^ 64) i +
          (p - 1) * slot B (blk B H j) i) % p = 0 ∧
        slot B (blk B H j) i < p ∧
        (slot B (blk B H j) i ≠ 0 → popcount 6 i + popcount 6 (e1 + j) ≤ d)
  | 0, e1, lo, hi, H, hs => by
    simp only [scanBlocks, Bool.and_eq_true, beq_iff_eq] at hs
    obtain ⟨rfl, rfl⟩ := hs
    exact ⟨by simp, by simp, fun j hj => absurd hj (Nat.not_lt_zero _)⟩
  | n+1, e1, lo, hi, H, hs => by
    simp only [scanBlocks, Bool.and_eq_true] at hs
    obtain ⟨hA, hB⟩ := hs
    obtain ⟨a0, a⟩ := scanSlots_sound _ _ _ _ _ _ _ _ _ _ _ hA
    obtain ⟨r1, r2, r⟩ := scanBlocks_sound B p d n (e1 + 1) _ _ _ hB
    refine ⟨by rw [pow_succ', ← Nat.div_div_eq_div_mul]; exact r1,
      by rw [pow_succ', ← Nat.div_div_eq_div_mul]; exact r2, ?_⟩
    intro j hj
    cases j with
    | zero =>
      have e : ∀ X, blk B X 0 = X % B ^ 128 := fun X => by simp [blk]
      simp only [e, add_zero]
      refine ⟨a0, fun i hi => ?_⟩
      have := a i hi
      rw [zero_add] at this
      exact this
    | succ j =>
      have e : ∀ X, blk B (X / B ^ 128) j = blk B X (j + 1) := fun X => by
        unfold blk
        rw [Nat.div_div_eq_div_mul, ← pow_succ']
      have := r j (by omega)
      simp only [e] at this
      rw [show e1 + 1 + j = e1 + (j + 1) by omega] at this
      exact this

theorem slot_blk {B X j i : ℕ} (hi : i < 128) : slot B (blk B X j) i = slot B X (i + 128 * j) := by
  unfold blk
  rw [slot_mod hi, ← pow_mul, slot_div]

theorem slot_blk_lo {B X j i : ℕ} (hi : i < 64) :
    slot B (blk B X j % B ^ 64) i = slot B X (i + 128 * j) := by
  rw [slot_mod hi, slot_blk (by omega)]

theorem slot_blk_hi {B X j i : ℕ} (hi : i < 64) :
    slot B (blk B X j / B ^ 64) i = slot B X (64 + i + 128 * j) := by
  rw [slot_div, slot_blk (by omega), add_comm i 64]

theorem lt_of_div_eq_zero' {a b : ℕ} (hb : 0 < b) (h : a / b = 0) : a < b := by
  by_contra hge
  have := Nat.div_pos (not_lt.mp hge) hb
  omega

set_option exponentiation.threshold 20000 in
theorem scan_sound {B p d acc H : ℕ} (hB : 0 < B)
    (hs : scanBlocks B p d 64 0 (acc % B ^ 8192) (acc / B ^ 8192) H = true) :
    RowFacts B p d H ∧ AccFacts B p acc H := by
  obtain ⟨r1, r2, r⟩ := scanBlocks_sound B p d 64 0 _ _ _ hs
  have e8 : (B ^ 128) ^ 64 = B ^ 8192 := by rw [← pow_mul]
  rw [e8] at r1 r2
  have hpos : ∀ n, 0 < B ^ n := fun n => pow_pos hB n
  have hH : H < B ^ 8192 := lt_of_div_eq_zero' (hpos _) r2
  have hhi : acc / B ^ 8192 < B ^ 8192 := lt_of_div_eq_zero' (hpos _) r1
  have hacc : acc < B ^ 16384 := by
    have := (Nat.div_lt_iff_lt_mul (hpos 8192)).mp hhi
    rw [← pow_add, show (8192 + 8192 : ℕ) = 16384 from rfl] at this
    exact this
  constructor
  · refine ⟨hH.trans_le (Nat.pow_le_pow_right hB (by norm_num)), ?_, ?_, ?_, ?_, ?_⟩
    · intro e1 h1 e0 h0
      rw [← slot_blk (by omega : 64 + e0 < 128)]
      exact slot_eq_zero_of_le hB (lt_of_div_eq_zero' (hpos _) (r e1 h1).1) (by omega)
    · intro e1 h1 e0 h0
      exact slot_eq_zero_of_le hB hH (by omega)
    · intro e1 h1 e0 h0
      exact slot_eq_zero_of_le hB hH (by omega)
    · intro e1 h1 e0 h0
      have := ((r e1 h1).2 e0 h0).2.1
      rwa [slot_blk (by omega)] at this
    · intro e1 h1 e0 h0
      have := ((r e1 h1).2 e0 h0).2.2
      rwa [slot_blk (by omega), zero_add] at this
  · refine ⟨hacc, ?_⟩
    intro e1 h1 e0 h0
    have := ((r e1 h1).2 e0 h0).1
    rw [slot_blk_lo h0, slot_blk_hi h0, slot_blk_lo h0, slot_blk_hi h0, slot_blk (by omega),
      slot_mod (by omega), slot_mod (by omega), slot_div, slot_div] at this
    rw [show e0 + 128 * e1 + 8192 = e0 + 128 * (64 + e1) by omega,
      show 64 + e0 + 128 * e1 + 8192 = 64 + e0 + 128 * (64 + e1) by omega] at this
    exact this

/-! ## The accumulator as a polynomial value -/

/-- The polynomial whose value at `B` is the accumulator. -/
def accPoly (u : UData) : List Term → Polynomial ℕ
  | [] => 0
  | t :: r => Polynomial.C t.1 * (Polynomial.X ^ t.2.1 * tpoly u.B (u.row t.2.2)) + accPoly u r

theorem eval_accPoly (u : UData) : ∀ terms : List Term,
    (∀ t ∈ terms, RowFacts u.B u.p (u.deg t.2.2) (u.row t.2.2)) →
    (accPoly u terms).eval u.B = ExactK6Cert.accum u terms
  | [], _ => by simp [accPoly, ExactK6Cert.accum]
  | t :: r, h => by
    rw [accPoly, ExactK6Cert.accum, Polynomial.eval_add, eval_accPoly u r (fun x hx => h x (List.mem_cons_of_mem _ hx)),
      Polynomial.eval_mul, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow,
      Polynomial.eval_X, eval_tpoly (h t (List.mem_cons_self)), Nat.shiftLeft_eq, UData.B, pow_mul,
      mul_comm ((2 ^ u.w) ^ t.2.1)]

theorem coeff_accPoly_le (u : UData) (hp : 0 < u.p) : ∀ terms : List Term,
    (∀ t ∈ terms, RowFacts u.B u.p (u.deg t.2.2) (u.row t.2.2)) →
    ∀ c, (accPoly u terms).coeff c ≤ coefSum terms * u.p
  | [], _, c => by simp [accPoly, coefSum]
  | t :: r, h, c => by
    have ih := coeff_accPoly_le u hp r (fun x hx => h x (List.mem_cons_of_mem _ hx)) c
    have ht := h t (List.mem_cons_self)
    rw [accPoly, coefSum, Polynomial.coeff_add, Polynomial.coeff_C_mul, Polynomial.coeff_X_pow_mul',
      add_mul]
    apply Nat.add_le_add _ ih
    apply Nat.mul_le_mul_left
    split_ifs
    · exact (coeff_tpoly_lt ht hp _).le
    · exact Nat.zero_le _

section Main

variable (two_ne : (2 : K) ≠ 0)
include two_ne

/-- One term: the translate of a stored row, as a reduced product. -/
theorem Λ_term (a s : ℕ) (hs0 : s % 128 < 64) (hs1 : s / 128 < 64) (B H : ℕ) :
    Λ K (Polynomial.C a * (Polynomial.X ^ s * tpoly B H)) =
      (a : K) • normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6) (M (K := K) s * vec B H) := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
  rw [← nfMap_apply _ _ zero]
  unfold tpoly vec
  simp only [Finset.mul_sum, Polynomial.X_pow_mul_monomial, Polynomial.C_mul_monomial, Λ_sum,
    map_sum, Finset.smul_sum]
  apply Finset.sum_congr rfl
  intro e1 h1
  apply Finset.sum_congr rfl
  intro e0 h0
  rw [Finset.mem_range] at h0 h1
  rw [Λ_monomial (by omega), mul_smul_comm, map_smul, nfMap_apply,
    nf_M_mul two_ne (by omega), Nat.cast_mul, mul_smul]

theorem Λ_accPoly (u : UData) : ∀ terms : List Term,
    (∀ t ∈ terms, t.2.1 % 128 < 64 ∧ t.2.1 / 128 < 64) →
    Λ K (accPoly u terms) = (terms.map fun t => (t.1 : K) •
      normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6)
        (M (K := K) t.2.1 * vec u.B (u.row t.2.2))).sum
  | [], _ => by simp [accPoly, Λ_zero]
  | t :: r, h => by
    rw [accPoly, Λ_add, Λ_accPoly u r (fun x hx => h x (List.mem_cons_of_mem _ hx)),
      Λ_term two_ne _ _ (h t (List.mem_cons_self)).1 (h t (List.mem_cons_self)).2, List.map_cons,
      List.sum_cons]

omit two_ne in
theorem cast_eq_of_mod {p a x : ℕ} (hp1 : 1 ≤ p) (hpK : (p : K) = 0)
    (h : (a + (p - 1) * x) % p = 0) : (a : K) = (x : K) := by
  obtain ⟨q, hq⟩ := Nat.dvd_of_mod_eq_zero h
  have := congrArg (Nat.cast : ℕ → K) hq
  rw [Nat.cast_add, Nat.cast_mul, Nat.cast_mul, Nat.cast_sub hp1, hpK, Nat.cast_one] at this
  linear_combination this

variable (F : Set (P (2 * 6) K))

/-- Row `i` is a well-formed row whose vector is in the quotient closure of degree 5. -/
def Good (u : UData) (i : ℕ) : Prop :=
  RowFacts u.B u.p (u.deg i) (u.row i) ∧
    InQuotClosure (chainInputs K 2 6) (chainInputs_nonzero 2 6) F 5 (vec (K := K) u.B (u.row i))

omit two_ne in
theorem B_pos (u : UData) : 0 < u.B := Nat.two_pow_pos _

set_option exponentiation.threshold 20000 in
/-- **The induction step.** A certified row is a good row if all earlier rows are. -/
theorem row_sound (u : UData) (hpK : (u.p : K) = 0) (hg : uGlobalOK u = true) (n : ℕ)
    (terms : List Term) (hc : rowCheck u n terms = true) (prev : ∀ i < n, Good F u i) :
    Good F u n := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
  simp only [uGlobalOK, Bool.and_eq_true, decide_eq_true_eq] at hg
  simp only [rowCheck, Bool.and_eq_true, List.all_eq_true, decide_eq_true_eq] at hc
  obtain ⟨⟨h1, h2⟩, h3⟩ := hc
  have hterm : ∀ t ∈ terms, t.2.2 < n ∧ t.2.1 % 128 < 64 ∧ t.2.1 / 128 < 64 ∧
      popcount 6 (t.2.1 % 128) + popcount 6 (t.2.1 / 128) + u.deg t.2.2 ≤ 5 := by
    intro t ht
    have := h1 t ht
    simp only [termOK, Bool.and_eq_true, decide_eq_true_eq] at this
    exact ⟨this.1.1.1, this.1.1.2, this.1.2, this.2⟩
  obtain ⟨rf, af⟩ := scan_sound (B_pos u) h3
  refine ⟨rf, ?_⟩
  have src : ∀ t ∈ terms, RowFacts u.B u.p (u.deg t.2.2) (u.row t.2.2) :=
    fun t ht => (prev _ (hterm t ht).1).1
  have hev := eval_accPoly u terms src
  have hco : ∀ c, (accPoly u terms).coeff c < u.B :=
    fun c => lt_of_le_of_lt (coeff_accPoly_le u (by omega) terms src c) h2
  have hslot : ∀ c, slot u.B (ExactK6Cert.accum u terms) c = (accPoly u terms).coeff c := fun c => by
    rw [← hev]
    exact slot_eval _ c _ hco
  have hvec : vec (K := K) u.B (u.row n) = Λ K (accPoly u terms) := by
    rw [Λ_apply]
    simp only [← hslot]
    rw [regroup]
    unfold vec
    apply Finset.sum_congr rfl
    intro e1 he1
    apply Finset.sum_congr rfl
    intro e0 he0
    rw [Finset.mem_range] at he0 he1
    rw [M_congr (K := K) (c := 64 + e0 + 128 * e1) (c' := e0 + 128 * e1) (by omega) (by omega),
      M_congr (K := K) (c := e0 + 128 * (64 + e1)) (c' := e0 + 128 * e1) (by omega) (by omega),
      M_congr (K := K) (c := 64 + e0 + 128 * (64 + e1)) (c' := e0 + 128 * e1) (by omega) (by omega),
      ← add_smul, ← add_smul, ← add_smul, ← Nat.cast_add, ← Nat.cast_add, ← Nat.cast_add]
    congr 1
    exact (cast_eq_of_mod (by omega) hpK (af.fold e1 he1 e0 he0)).symm
  rw [hvec, Λ_accPoly two_ne u terms (fun t ht => ⟨(hterm t ht).2.1, (hterm t ht).2.2.1⟩)]
  have key : ∀ L : List Term, (∀ t ∈ L, t ∈ terms) →
      InQuotClosure (chainInputs K 2 6) (chainInputs_nonzero 2 6) F 5
        (L.map fun t => (t.1 : K) • normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6)
          (M (K := K) t.2.1 * vec u.B (u.row t.2.2))).sum := by
    intro L
    induction L with
    | nil => intro _; exact .zero
    | cons t r ih =>
      intro hL
      rw [List.map_cons, List.sum_cons]
      refine .add (.smul _ ?_) (ih (fun x hx => hL x (List.mem_cons_of_mem _ hx)))
      obtain ⟨hlt, hs0, hs1, hbud⟩ := hterm t (hL t (List.mem_cons_self))
      obtain ⟨rfi, hin⟩ := prev _ hlt
      apply quot_mul_mem zero _ hin
      have d1 := M_degree (K := K) hs0 hs1
      have d2 := vec_degree (K := K) rfi
      omega
  exact key terms (fun t ht => ht)

/-! ## The generator row and the final row -/

/-- The polynomial of a sparse list. -/
def listPoly : List Entry → Polynomial ℕ
  | [] => 0
  | e :: r => Polynomial.monomial e.1 e.2 + listPoly r

omit two_ne in
theorem eval_listPoly (w : ℕ) : ∀ L : List Entry, (listPoly L).eval (2 ^ w) = packEntries w L
  | [] => by simp [listPoly, packEntries]
  | e :: r => by
    rw [listPoly, packEntries, Polynomial.eval_add, eval_listPoly w r, Polynomial.eval_monomial,
      Nat.shiftLeft_eq, pow_mul]

omit two_ne in
theorem coeff_listPoly_le : ∀ (L : List Entry) (c : ℕ), (listPoly L).coeff c ≤ entrySum L
  | [], c => by simp [listPoly, entrySum]
  | e :: r, c => by
    rw [listPoly, entrySum, Polynomial.coeff_add, Polynomial.coeff_monomial]
    apply Nat.add_le_add _ (coeff_listPoly_le r c)
    split_ifs
    · exact le_refl _
    · exact Nat.zero_le _

omit two_ne in
theorem Λ_listPoly : ∀ L : List Entry, (∀ e ∈ L, e.1 < 16384) →
    Λ K (listPoly L) = (L.map fun e => (e.2 : K) • M (K := K) e.1).sum
  | [], _ => by simp [listPoly, Λ_zero]
  | e :: r, h => by
    rw [listPoly, Λ_add, Λ_monomial (h e (List.mem_cons_self)),
      Λ_listPoly r (fun x hx => h x (List.mem_cons_of_mem _ hx)), List.map_cons, List.sum_cons]

set_option exponentiation.threshold 20000 in
omit two_ne in
/-- Row `0` is a well-formed row whose vector is the explicit generator sum. -/
theorem gen_sound (u : UData) (hgen : genOK u = true) :
    RowFacts u.B u.p (u.deg 0) (u.row 0) ∧
      vec (K := K) u.B (u.row 0) = (u.gen.map fun e => (e.2 : K) • M (K := K) e.1).sum := by
  simp only [genOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq, List.all_eq_true] at hgen
  obtain ⟨⟨⟨h1, h2⟩, h3⟩, h4⟩ := hgen
  obtain ⟨rf, _⟩ := scan_sound (B_pos u) h4
  refine ⟨rf, ?_⟩
  have hslot : ∀ c, slot u.B (u.row 0) c = (listPoly u.gen).coeff c := fun c => by
    rw [h1, ← eval_listPoly]
    exact slot_eval _ c _ (fun i => lt_of_le_of_lt (coeff_listPoly_le _ _) h2)
  rw [vec_eq_sum rf, ← Λ_listPoly u.gen h3, Λ_apply]
  simp only [hslot]

omit two_ne in
/-- A well-formed row equal to `1` has vector `1`. -/
theorem one_sound {B p d : ℕ} (hB : 2 ≤ B) (rf : RowFacts B p d 1) : vec (K := K) B 1 = 1 := by
  rw [vec_eq_sum rf, Finset.sum_eq_single 0]
  · rw [slot_zero, Nat.mod_eq_of_lt (by omega), Nat.cast_one, one_smul, M_zero]
  · intro c _ hc
    rw [slot_eq_zero (Nat.one_lt_pow hc (by omega)), Nat.cast_zero, zero_smul]
  · intro h
    exact absurd (Finset.mem_range.mpr (by norm_num)) h

/-- **Soundness of a derivation.** If the generator sum is in the quotient closure of degree 5,
every row `1 ≤ n ≤ N` passes `rowCheck`, and row `N` is `1`, then `1` is in the quotient
closure of degree 5. -/
theorem derivation_sound (u : UData) (hpK : (u.p : K) = 0) (hg : uGlobalOK u = true)
    (hgen : genOK u = true)
    (hin : InQuotClosure (chainInputs K 2 6) (chainInputs_nonzero 2 6) F 5
      (u.gen.map fun e => (e.2 : K) • M (K := K) e.1).sum)
    (N : ℕ) (terms : ℕ → List Term)
    (rows : ∀ n, 1 ≤ n → n ≤ N → rowCheck u n (terms n) = true) (fin : u.row N = 1) :
    InQuotClosure (chainInputs K 2 6) (chainInputs_nonzero 2 6) F 5 (1 : P (2 * 6) K) := by
  have all : ∀ n, n ≤ N → Good F u n := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro hn
      rcases Nat.eq_zero_or_pos n with rfl | hpos
      · obtain ⟨rf, hv⟩ := gen_sound (K := K) u hgen
        exact ⟨rf, by rw [hv]; exact hin⟩
      · exact row_sound two_ne F u hpK hg n (terms n) (rows n hpos hn)
          (fun i hi => ih i hi (by omega))
  obtain ⟨rf, hq⟩ := all N (le_refl _)
  rw [fin] at rf hq
  have hw : 2 ≤ u.B := by
    simp only [uGlobalOK, Bool.and_eq_true, decide_eq_true_eq] at hg
    have : 2 ^ 1 ≤ 2 ^ u.w := Nat.pow_le_pow_right (by norm_num) hg.2
    simpa [UData.B] using this
  rwa [one_sound hw rf] at hq

end Main

end ICFallDegrees.ExactK6Upper

end
