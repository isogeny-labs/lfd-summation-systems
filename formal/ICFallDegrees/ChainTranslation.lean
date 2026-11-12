import ICFallDegrees.ChainClosureBridge
import ICFallDegrees.DepthReduction

noncomputable section

/-!
# The translation identity: chain normal form of a variable times a squarefree monomial

For an exponent vector `b ∈ [0,2^k)^m` let `digits k b` be the squarefree chain exponent
whose `(i,j)` entry is bit `j` of `b i`.  The chain variable `u_{i,j}` times the monomial of
`b` is congruent, modulo the chain ideal, to the monomial of `b` with `2^j` added to block `i`
modulo `2^k` (`translate_monomial_mem_span`): the proof is the carry induction, `u_{i,j}^2 ≡
u_{i,j+1}` pushing a set bit one level up, and `u_{i,k-1}^2 ≡ 1` discarding it at the top.
Hence the chain normal form of `u_{i,j} · m_b` is `m_{b + 2^j e_i mod 2^k}`
(`normalForm_X_mul_monomial`).  This is Lemma 2.2 of `docs/proofs/N5-reduction.md`.
-/

namespace ICFallDegrees.ChainTranslation

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight

variable {m k : ℕ} {K : Type*} [Field K]

/-! ## Bits -/

/-- Bit `j` of `n`, as a natural number. -/
def bit (n j : ℕ) : ℕ := n / 2 ^ j % 2

theorem bit_le_one (n j : ℕ) : bit n j ≤ 1 := by
  unfold bit
  have := Nat.mod_lt (n / 2 ^ j) (by norm_num : 0 < 2)
  omega

/-- Adding a higher bit does not change a lower bit. -/
theorem bit_add_two_pow_of_lt {j j' : ℕ} (h : j' < j) (n : ℕ) : bit (n + 2 ^ j) j' = bit n j' := by
  unfold bit
  have hsplit : 2 ^ j = 2 ^ (j - j' - 1) * 2 * 2 ^ j' := by
    rw [← pow_succ, ← pow_add]
    congr 1
    omega
  rw [hsplit, Nat.add_mul_div_right _ _ (Nat.two_pow_pos j'), Nat.add_mul_mod_self_right]

/-- Setting a clear bit. -/
theorem bit_add_two_pow_self {n j : ℕ} (h : bit n j = 0) : bit (n + 2 ^ j) j = 1 := by
  unfold bit at *
  rw [Nat.add_div_right _ (Nat.two_pow_pos j)]
  omega

/-- A clear bit `j` means the residue modulo `2^(j+1)` is below `2^j`. -/
theorem mod_lt_of_bit_zero {n j : ℕ} (h : bit n j = 0) : n % 2 ^ (j + 1) < 2 ^ j := by
  have hpos : 0 < 2 ^ j := Nat.two_pow_pos j
  have hsplit : 2 ^ (j + 1) = 2 * 2 ^ j := pow_succ' 2 j
  have hn : n = 2 ^ (j + 1) * (n / 2 ^ (j + 1)) + n % 2 ^ (j + 1) := (Nat.div_add_mod _ _).symm
  set A := n / 2 ^ (j + 1)
  set r := n % 2 ^ (j + 1) with hr
  have hr2 : r < 2 * 2 ^ j := by rw [← hsplit]; exact Nat.mod_lt _ (Nat.two_pow_pos _)
  have hdiv : n / 2 ^ j = 2 * A + r / 2 ^ j := by
    rw [hn, hsplit, show 2 * 2 ^ j * A = 2 ^ j * (2 * A) by ring, Nat.mul_add_div hpos]
  have hbit : bit n j = r / 2 ^ j % 2 := by
    unfold bit
    rw [hdiv, Nat.mul_add_mod]
  have hq : r / 2 ^ j < 2 := (Nat.div_lt_iff_lt_mul hpos).mpr hr2
  have hq0 : r / 2 ^ j = 0 := by
    rw [hbit, Nat.mod_eq_of_lt hq] at h
    exact h
  have := (Nat.div_lt_iff_lt_mul hpos).mp (by rw [hq0]; exact Nat.zero_lt_one)
  simpa using this

/-- With bit `j` clear, adding `2^j` does not disturb the quotient by `2^(j+1)`. -/
theorem div_succ_add_two_pow {n j : ℕ} (h : bit n j = 0) :
    (n + 2 ^ j) / 2 ^ (j + 1) = n / 2 ^ (j + 1) := by
  have hsplit : 2 ^ (j + 1) = 2 ^ j * 2 := pow_succ 2 j
  have hpos : 0 < 2 ^ (j + 1) := Nat.two_pow_pos _
  have hn : n = 2 ^ (j + 1) * (n / 2 ^ (j + 1)) + n % 2 ^ (j + 1) := (Nat.div_add_mod _ _).symm
  have hr := mod_lt_of_bit_zero h
  have hsmall : n % 2 ^ (j + 1) + 2 ^ j < 2 ^ (j + 1) := by omega
  conv_lhs => rw [hn]
  rw [add_assoc, Nat.mul_add_div hpos, Nat.div_eq_of_lt hsmall, add_zero]

/-- Adding a clear bit does not change the higher bits. -/
theorem bit_add_two_pow_of_gt {j j' : ℕ} (h : j < j') {n : ℕ} (h0 : bit n j = 0) :
    bit (n + 2 ^ j) j' = bit n j' := by
  unfold bit
  have hsplit : 2 ^ j' = 2 ^ (j + 1) * 2 ^ (j' - j - 1) := by
    rw [← pow_add]
    congr 1
    omega
  rw [hsplit, ← Nat.div_div_eq_div_mul, ← Nat.div_div_eq_div_mul, div_succ_add_two_pow h0]

/-- A set bit forces the number to be at least `2^j`. -/
theorem two_pow_le_of_bit_one {n j : ℕ} (h : bit n j = 1) : 2 ^ j ≤ n := by
  by_contra hlt
  push_neg at hlt
  unfold bit at h
  rw [Nat.div_eq_of_lt hlt] at h
  simp at h

/-- Clearing a set bit. -/
theorem bit_sub_two_pow {n j : ℕ} (h : bit n j = 1) : bit (n - 2 ^ j) j = 0 := by
  have hle := two_pow_le_of_bit_one h
  unfold bit at *
  have hdiv : (n - 2 ^ j) / 2 ^ j = n / 2 ^ j - 1 := by
    have := Nat.sub_mul_div n (2 ^ j) 1 (by simpa using hle)
    simpa using this
  rw [hdiv]
  have hpos : 0 < n / 2 ^ j := Nat.div_pos hle (Nat.two_pow_pos j)
  omega

/-- Adding a clear bit inside a `k`-bit window stays inside it. -/
theorem add_two_pow_lt_of_bit_zero {n j k : ℕ} (hj : j < k) (hn : n < 2 ^ k) (h : bit n j = 0) :
    n + 2 ^ j < 2 ^ k := by
  by_contra hge
  push_neg at hge
  have hsplit : 2 ^ k = 2 ^ (j + 1) * 2 ^ (k - j - 1) := by
    rw [← pow_add]
    congr 1
    omega
  have hzero : (n + 2 ^ j) / 2 ^ k = 0 := by
    rw [hsplit, ← Nat.div_div_eq_div_mul, div_succ_add_two_pow h, Nat.div_div_eq_div_mul,
      ← hsplit]
    exact Nat.div_eq_of_lt hn
  have hpos : 0 < (n + 2 ^ j) / 2 ^ k := Nat.div_pos hge (Nat.two_pow_pos k)
  omega

/-- Adding the top bit when it is set wraps around. -/
theorem add_top_bit_mod {n k : ℕ} (hk : 0 < k) (hn : n < 2 ^ k) (h : bit n (k - 1) = 1) :
    (n + 2 ^ (k - 1)) % 2 ^ k = n - 2 ^ (k - 1) := by
  have hle := two_pow_le_of_bit_one h
  have hsplit : 2 ^ k = 2 ^ (k - 1) * 2 := by
    rw [← pow_succ]
    congr 1
    omega
  have hsum : n + 2 ^ (k - 1) = (n - 2 ^ (k - 1)) + 2 ^ k := by omega
  rw [hsum, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]

/-! ## Squarefree chain exponents of an exponent vector -/

/-- The squarefree chain exponent of an exponent vector: entry `(i,j)` is bit `j` of `b i`. -/
def digits (k : ℕ) (b : Fin m → ℕ) : Fin (m * k) →₀ ℕ :=
  Finsupp.equivFunOnFinite.symm (fun v => bit (b (row v)) (depth v))

@[simp] theorem digits_apply (b : Fin m → ℕ) (v : Fin (m * k)) :
    digits k b v = bit (b (row v)) (depth v) := rfl

theorem digits_le_one (b : Fin m → ℕ) (v : Fin (m * k)) : digits k b v ≤ 1 := bit_le_one _ _

theorem digits_reduced (b : Fin m → ℕ) (c : K) :
    monomial (digits k b) c ∈ reducedSpace (chainInputs K m k) := by
  classical
  intro e he g hg
  obtain ⟨v, rfl⟩ := hg
  rw [leading_chainPolynomial, Finsupp.single_le_iff]
  have : e = digits k b := by
    by_cases hc : c = 0
    · simp [hc] at he
    · rw [support_monomial, if_neg hc, Finset.mem_singleton] at he
      exact he
  rw [this]
  have := digits_le_one (k := k) b v
  omega

theorem var_eq_iff {i : Fin m} {j : Fin k} {v : Fin (m * k)} :
    var m k i j = v ↔ row v = i ∧ depth v = j := by
  constructor
  · rintro rfl
    simp
  · rintro ⟨h1, h2⟩
    rw [← var_row_depth v, h1, h2]

/-- Setting a clear bit `j` in block `i` adds the single exponent `u_{i,j}`. -/
theorem digits_update_of_bit_zero {j : ℕ} (hj : j < k) (i : Fin m) (b : Fin m → ℕ)
    (h : bit (b i) j = 0) :
    digits k (Function.update b i (b i + 2 ^ j)) =
      digits k b + Finsupp.single (var m k i ⟨j, hj⟩) 1 := by
  classical
  ext v
  rw [Finsupp.add_apply, digits_apply, digits_apply, Finsupp.single_apply]
  by_cases hr : row v = i
  · rw [hr, Function.update_self]
    by_cases hd : (depth v : ℕ) = j
    · have hv : var m k i ⟨j, hj⟩ = v := var_eq_iff.mpr ⟨hr, Fin.ext hd⟩
      rw [if_pos hv, hd, bit_add_two_pow_self h, h]
    · have hv : ¬ var m k i ⟨j, hj⟩ = v := fun hv => hd (by
        have := (var_eq_iff.mp hv).2
        rw [this])
      rw [if_neg hv, add_zero]
      rcases Nat.lt_or_gt_of_ne hd with hlt | hgt
      · exact bit_add_two_pow_of_lt hlt _
      · exact bit_add_two_pow_of_gt hgt h
  · rw [Function.update_of_ne hr, if_neg (fun hv => hr (var_eq_iff.mp hv).1), add_zero]

/-- A set bit `j` in block `i` splits off the single exponent `u_{i,j}`. -/
theorem digits_of_bit_one {j : ℕ} (hj : j < k) (i : Fin m) (b : Fin m → ℕ) (h : bit (b i) j = 1) :
    digits k b =
      digits k (Function.update b i (b i - 2 ^ j)) + Finsupp.single (var m k i ⟨j, hj⟩) 1 := by
  have h0 : bit (Function.update b i (b i - 2 ^ j) i) j = 0 := by
    rw [Function.update_self]
    exact bit_sub_two_pow h
  have key := digits_update_of_bit_zero hj i _ h0
  rw [Function.update_self, Function.update_idem, Nat.sub_add_cancel (two_pow_le_of_bit_one h),
    Function.update_eq_self] at key
  exact key

/-! ## The carry congruence -/

theorem monomial_add_single_single (e : Fin (m * k) →₀ ℕ) (v : Fin (m * k)) :
    monomial (e + Finsupp.single v 1 + Finsupp.single v 1) (1 : K) = X v ^ 2 * monomial e 1 := by
  rw [X_pow_eq_monomial, monomial_mul, one_mul]
  congr 1
  rw [add_assoc, ← Finsupp.single_add, add_comm]

theorem X_mul_monomial_one (e : Fin (m * k) →₀ ℕ) (v : Fin (m * k)) :
    X v * monomial e (1 : K) = monomial (e + Finsupp.single v 1) 1 := by
  rw [X, monomial_mul, one_mul, add_comm]

/-- The carry congruence: `u_{i,j} · m_b ≡ m_{b + 2^j e_i mod 2^k}` modulo the chain ideal,
for `b` in the depth-`k` box. The induction is on `d = k - 1 - j`. -/
theorem translate_monomial_mem_span (i : Fin m) : ∀ (d j : ℕ) (hjd : j + d + 1 = k)
    (b : Fin m → ℕ) (_ : ∀ i', b i' < 2 ^ k),
    monomial (digits k b + Finsupp.single (var m k i ⟨j, by omega⟩) 1) (1 : K) -
      monomial (digits k (Function.update b i ((b i + 2 ^ j) % 2 ^ k))) 1 ∈
        Ideal.span (chainInputs K m k) := by
  intro d
  induction d with
  | zero =>
    intro j hjd b hb
    have hj : j < k := by omega
    by_cases h0 : bit (b i) j = 0
    · rw [Nat.mod_eq_of_lt (add_two_pow_lt_of_bit_zero hj (hb i) h0),
        digits_update_of_bit_zero hj i b h0, sub_self]
      exact (Ideal.span _).zero_mem
    · have h1 : bit (b i) j = 1 := by have := bit_le_one (b i) j; omega
      have hjk : j = k - 1 := by omega
      rw [digits_of_bit_one hj i b h1, monomial_add_single_single]
      have hwrap : (b i + 2 ^ j) % 2 ^ k = b i - 2 ^ j := by
        rw [hjk]
        exact add_top_bit_mod (by omega) (hb i) (hjk ▸ h1)
      rw [hwrap]
      have htail : ChainQuotient.chainTail K m k (var m k i ⟨j, hj⟩) = 1 := by
        unfold ChainQuotient.chainTail
        rw [dif_neg]
        simp only [depth_var]
        omega
      have hchain : X (var m k i ⟨j, hj⟩) ^ 2 = chainPolynomial K m k (var m k i ⟨j, hj⟩) + 1 := by
        unfold chainPolynomial
        rw [htail]
        ring
      rw [hchain, add_mul, one_mul, add_sub_cancel_right]
      exact (Ideal.span _).mul_mem_right _ (Ideal.subset_span ⟨_, rfl⟩)
  | succ d ih =>
    intro j hjd b hb
    have hj : j < k := by omega
    by_cases h0 : bit (b i) j = 0
    · rw [Nat.mod_eq_of_lt (add_two_pow_lt_of_bit_zero hj (hb i) h0),
        digits_update_of_bit_zero hj i b h0, sub_self]
      exact (Ideal.span _).zero_mem
    · have h1 : bit (b i) j = 1 := by have := bit_le_one (b i) j; omega
      have hle := two_pow_le_of_bit_one h1
      set b' := Function.update b i (b i - 2 ^ j) with hb'
      have hb'box : ∀ i', b' i' < 2 ^ k := by
        intro i'
        rcases eq_or_ne i' i with rfl | hne
        · rw [hb', Function.update_self]
          have := hb i'
          have := Nat.two_pow_pos j
          omega
        · rw [hb', Function.update_of_ne hne]
          exact hb i'
      have hj1 : j + 1 < k := by omega
      have htail : ChainQuotient.chainTail K m k (var m k i ⟨j, hj⟩) =
          X (var m k i ⟨j + 1, hj1⟩) := by
        unfold ChainQuotient.chainTail
        rw [dif_pos (by simp only [depth_var]; exact hj1)]
        simp only [row_var, depth_var]
      have hchain : X (var m k i ⟨j, hj⟩) ^ 2 =
          chainPolynomial K m k (var m k i ⟨j, hj⟩) + X (var m k i ⟨j + 1, hj1⟩) := by
        unfold chainPolynomial
        rw [htail]
        ring
      have hupd : Function.update b' i ((b' i + 2 ^ (j + 1)) % 2 ^ k) =
          Function.update b i ((b i + 2 ^ j) % 2 ^ k) := by
        rw [hb', Function.update_self, Function.update_idem]
        have e : b i - 2 ^ j + 2 ^ (j + 1) = b i + 2 ^ j := by rw [pow_succ]; omega
        rw [e]
      have step := ih (j + 1) (by omega) b' hb'box
      rw [hupd] at step
      have hdig : digits k b = digits k b' + Finsupp.single (var m k i ⟨j, hj⟩) 1 := by
        rw [hb']
        exact digits_of_bit_one hj i b h1
      rw [hdig, monomial_add_single_single, hchain, add_mul, X_mul_monomial_one]
      have identity : chainPolynomial K m k (var m k i ⟨j, hj⟩) * monomial (digits k b') 1 +
          monomial (digits k b' + Finsupp.single (var m k i ⟨j + 1, hj1⟩) 1) (1 : K) -
          monomial (digits k (Function.update b i ((b i + 2 ^ j) % 2 ^ k))) 1 =
          chainPolynomial K m k (var m k i ⟨j, hj⟩) * monomial (digits k b') 1 +
          (monomial (digits k b' + Finsupp.single (var m k i ⟨j + 1, hj1⟩) 1) (1 : K) -
            monomial (digits k (Function.update b i ((b i + 2 ^ j) % 2 ^ k))) 1) := by ring
      rw [identity]
      exact (Ideal.span _).add_mem
        ((Ideal.span _).mul_mem_right _ (Ideal.subset_span ⟨_, rfl⟩)) step

/-- The translation identity for the chain normal form (Lemma 2.2 of the written proof). -/
theorem normalForm_X_mul_monomial (two_ne : (2 : K) ≠ 0) (i : Fin m) (j : Fin k)
    (b : Fin m → ℕ) (hb : ∀ i', b i' < 2 ^ k) :
    ChainClosureBridge.normalForm (chainInputs K m k) (chainInputs_nonzero m k)
        (X (var m k i j) * monomial (digits k b) (1 : K)) =
      monomial (digits k (Function.update b i ((b i + 2 ^ (j : ℕ)) % 2 ^ k))) 1 := by
  apply ChainClosureBridge.normalForm_eq_of _ _ (reduced_chain_ideal_trivial m k two_ne)
    (digits_reduced _ _)
  rw [X_mul_monomial_one]
  have := translate_monomial_mem_span (K := K) i (k - 1 - j) j (by have := j.isLt; omega) b hb
  simpa using this

end ICFallDegrees.ChainTranslation
