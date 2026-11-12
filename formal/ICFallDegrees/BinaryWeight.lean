import Mathlib.Data.Nat.Digits
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-!
# Binary weight and single-bit carries

`wt n` is the binary digit sum (popcount) of `n`.  The lemmas record what one
addition of a power of two does to the weight: it rises by at most one
(`wt_add_two_pow_le`), it does not rise when the added bit was set
(`wt_add_two_pow_le_of_odd`), a carry out of a `k`-bit window strictly lowers
the weight of the window remainder (`carry_wt`), and a wrap-around at the top of
a `(Z/2^n)`-block is harmless (`wt_mod_two_pow_add_le`).  `wt_add_mul_two_pow`
splits the weight of a number into its low-`k` and high parts, which is the
direct-sum decomposition by high pattern used in the depth reduction.
-/

namespace ICFallDegrees.BinaryWeight

/-- Binary digit sum. -/
def wt (n : ℕ) : ℕ := (Nat.digits 2 n).sum

@[simp] theorem wt_zero : wt 0 = 0 := by simp [wt]

/-- The defining recursion, valid for every `n`. -/
theorem wt_eq (n : ℕ) : wt n = n % 2 + wt (n / 2) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [wt]
  · unfold wt
    rw [Nat.digits_def' (by norm_num) hn, List.sum_cons]

theorem wt_one : wt 1 = 1 := by rw [wt_eq]; simp

theorem wt_two_mul_add (q b : ℕ) (hb : b < 2) : wt (2 * q + b) = b + wt q := by
  rw [wt_eq]
  have h1 : (2 * q + b) % 2 = b := by omega
  have h2 : (2 * q + b) / 2 = q := by omega
  rw [h1, h2]

/-- Adding one raises the weight by at most one. -/
theorem wt_succ_le (n : ℕ) : wt (n + 1) ≤ wt n + 1 := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    rw [wt_eq (n + 1), wt_eq n]
    rcases Nat.even_or_odd n with ⟨q, hq⟩ | ⟨q, hq⟩
    · have h1 : (n + 1) % 2 = 1 := by omega
      have h2 : (n + 1) / 2 = n / 2 := by omega
      have h3 : n % 2 = 0 := by omega
      rw [h1, h2, h3]
      omega
    · have h1 : (n + 1) % 2 = 0 := by omega
      have h2 : (n + 1) / 2 = n / 2 + 1 := by omega
      have h3 : n % 2 = 1 := by omega
      rw [h1, h2, h3]
      have := ih (n / 2) (by omega)
      omega

/-- Adding one bit raises the weight by at most one, whatever the carries. -/
theorem wt_add_two_pow_le (j a : ℕ) : wt (a + 2 ^ j) ≤ wt a + 1 := by
  induction j generalizing a with
  | zero => simpa using wt_succ_le a
  | succ j ih =>
    rw [wt_eq (a + 2 ^ (j + 1)), wt_eq a]
    have hp : 2 ^ (j + 1) = 2 ^ j * 2 := pow_succ 2 j
    have h1 : (a + 2 ^ (j + 1)) % 2 = a % 2 := by omega
    have h2 : (a + 2 ^ (j + 1)) / 2 = a / 2 + 2 ^ j := by omega
    rw [h1, h2]
    have := ih (a / 2)
    omega

/-- If the added bit is already set, the weight does not rise. -/
theorem wt_add_two_pow_le_of_odd (j a : ℕ) (h : a / 2 ^ j % 2 = 1) : wt (a + 2 ^ j) ≤ wt a := by
  induction j generalizing a with
  | zero =>
    simp only [pow_zero, Nat.div_one] at h ⊢
    rw [wt_eq (a + 1), wt_eq a]
    have h1 : (a + 1) % 2 = 0 := by omega
    have h2 : (a + 1) / 2 = a / 2 + 1 := by omega
    rw [h1, h2, h]
    have := wt_succ_le (a / 2)
    omega
  | succ j ih =>
    rw [wt_eq (a + 2 ^ (j + 1)), wt_eq a]
    have hp : 2 ^ (j + 1) = 2 ^ j * 2 := pow_succ 2 j
    have h1 : (a + 2 ^ (j + 1)) % 2 = a % 2 := by omega
    have h2 : (a + 2 ^ (j + 1)) / 2 = a / 2 + 2 ^ j := by omega
    have h3 : a / 2 / 2 ^ j % 2 = 1 := by
      rw [Nat.div_div_eq_div_mul, ← pow_succ']
      exact h
    rw [h1, h2]
    have := ih (a / 2) h3
    omega

/-- A set top bit above an `n`-bit remainder contributes exactly one. -/
theorem wt_two_pow_add (n r : ℕ) (hr : r < 2 ^ n) : wt (2 ^ n + r) = wt r + 1 := by
  induction n generalizing r with
  | zero =>
    have : r = 0 := by simpa using hr
    subst this
    simp [wt_one]
  | succ n ih =>
    rw [wt_eq (2 ^ (n + 1) + r), wt_eq r]
    have hp : 2 ^ (n + 1) = 2 ^ n * 2 := pow_succ 2 n
    have h1 : (2 ^ (n + 1) + r) % 2 = r % 2 := by omega
    have h2 : (2 ^ (n + 1) + r) / 2 = 2 ^ n + r / 2 := by omega
    rw [h1, h2, ih (r / 2) (by omega)]
    omega

/-- Weight splits across a low `n`-bit window and the part above it. -/
theorem wt_add_mul_two_pow (n l q : ℕ) (hl : l < 2 ^ n) : wt (l + 2 ^ n * q) = wt l + wt q := by
  induction n generalizing l with
  | zero =>
    have : l = 0 := by simpa using hl
    subst this
    simp
  | succ n ih =>
    rw [wt_eq (l + 2 ^ (n + 1) * q), wt_eq l]
    have hp : 2 ^ (n + 1) * q = 2 * (2 ^ n * q) := by rw [pow_succ]; ring
    have h1 : (l + 2 ^ (n + 1) * q) % 2 = l % 2 := by omega
    have h2 : (l + 2 ^ (n + 1) * q) / 2 = l / 2 + 2 ^ n * q := by omega
    rw [h1, h2, ih (l / 2) (by omega)]
    omega

/-- A carry out of an `n`-bit window happens only when the added bit was set. -/
theorem carry_odd (n j a : ℕ) (hj : j < n) (ha : a < 2 ^ n) (hc : 2 ^ n ≤ a + 2 ^ j) :
    a / 2 ^ j % 2 = 1 := by
  have hsplit : 2 ^ n = 2 ^ j * 2 ^ (n - j) := by
    rw [← pow_add]
    congr 1
    omega
  have hQ : 2 ^ (n - j) = 2 * 2 ^ (n - j - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  have hdiv : a / 2 ^ j = 2 ^ (n - j) - 1 := by
    apply Nat.div_eq_of_lt_le
    · rw [Nat.sub_one_mul, mul_comm, ← hsplit]
      omega
    · rw [Nat.sub_add_cancel (Nat.one_le_two_pow), mul_comm, ← hsplit]
      exact ha
  rw [hdiv, hQ]
  have : 0 < 2 ^ (n - j - 1) := Nat.two_pow_pos _
  omega

/-- A carry out of an `n`-bit window strictly lowers the weight of the window remainder. -/
theorem carry_wt (n j a : ℕ) (hj : j < n) (ha : a < 2 ^ n) (hc : 2 ^ n ≤ a + 2 ^ j) :
    wt (a + 2 ^ j - 2 ^ n) + 1 ≤ wt a := by
  have hr : a + 2 ^ j - 2 ^ n < 2 ^ n := by
    have : 2 ^ j < 2 ^ n := Nat.pow_lt_pow_right (by norm_num) hj
    omega
  have h1 := wt_two_pow_add n (a + 2 ^ j - 2 ^ n) hr
  have h2 : 2 ^ n + (a + 2 ^ j - 2 ^ n) = a + 2 ^ j := by omega
  rw [h2] at h1
  have h3 := wt_add_two_pow_le_of_odd j a (carry_odd n j a hj ha hc)
  omega

/-- Adding a bit inside a cyclic `n`-bit block raises the weight by at most one, even when the
addition wraps around at the top. -/
theorem wt_mod_two_pow_add_le (n j q : ℕ) (hj : j < n) (hq : q < 2 ^ n) :
    wt ((q + 2 ^ j) % 2 ^ n) ≤ wt q + 1 := by
  by_cases hlt : q + 2 ^ j < 2 ^ n
  · rw [Nat.mod_eq_of_lt hlt]
    exact wt_add_two_pow_le j q
  · have hjn : 2 ^ j < 2 ^ n := Nat.pow_lt_pow_right (by norm_num) hj
    have hr : q + 2 ^ j - 2 ^ n < 2 ^ n := by omega
    have hmod : (q + 2 ^ j) % 2 ^ n = q + 2 ^ j - 2 ^ n := by
      have hsum : q + 2 ^ j = (q + 2 ^ j - 2 ^ n) + 2 ^ n := by omega
      conv_lhs => rw [hsum]
      rw [Nat.add_mod_right, Nat.mod_eq_of_lt hr]
    rw [hmod]
    have h1 := wt_two_pow_add n (q + 2 ^ j - 2 ^ n) hr
    have h2 : 2 ^ n + (q + 2 ^ j - 2 ^ n) = q + 2 ^ j := by omega
    rw [h2] at h1
    have h3 := wt_add_two_pow_le j q
    omega

/-- Cyclic inverse of a one-bit step: adding back `d` after subtracting it modulo `N`. -/
theorem mod_add_sub_cancel (N x d : ℕ) (hx : x < N) (hd : d ≤ N) :
    ((x + N - d) % N + d) % N = x := by
  rw [Nat.mod_add_mod]
  have : x + N - d + d = x + N := by omega
  rw [this, Nat.add_mod_right, Nat.mod_eq_of_lt hx]

theorem mod_sub_add_cancel (N x d : ℕ) (hx : x < N) (hd : d ≤ N) :
    ((x + d) % N + N - d) % N = x := by
  have hpos : 0 < N := by omega
  have h1 : (x + d) % N + N - d = (x + d) % N + (N - d) := by omega
  rw [h1, Nat.mod_add_mod]
  have : x + d + (N - d) = x + N := by omega
  rw [this, Nat.add_mod_right, Nat.mod_eq_of_lt hx]

end ICFallDegrees.BinaryWeight
