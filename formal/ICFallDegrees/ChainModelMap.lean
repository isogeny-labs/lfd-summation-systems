import ICFallDegrees.ChainTranslation
import ICFallDegrees.DepthReductionBridge

noncomputable section

/-!
# The model map: exponent model → reduced chain polynomials

`encode k0 (l,q) i = l i + 2^k0 q i` recovers the exponent vector of a two-level exponent, and
`Φ k0 k = lmapDomain (digits k ∘ encode k0)` sends the basis vector `δ_(l,q)` of the exponent
model to the squarefree chain monomial of that exponent vector.  On the depth-`k` box, `Φ`
is injective, preserves the weight filtration (`Φ_degree`), and intertwines the model
translation with reduced multiplication by the chain variable (`Φ_step`), the last by the
carry congruence `ChainTranslation.normalForm_X_mul_monomial` and the model soundness
`encode_step` (Lemma 4.2 of `docs/proofs/N5-reduction.md`).  Hence `Φ` is a `ModelMap`
(`modelMap`), and N5 assembles without any interface hypothesis
(`n5_lfd_ge_five`, `prime_field_n5_lfd_ge_five`).
-/

namespace ICFallDegrees.ChainModelMap

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight DepthReduction ChainTranslation
open DepthReductionBridge

variable {m k : ℕ} {K : Type*} [Field K]

/-! ## Sums of bits -/

theorem bit_succ (n j : ℕ) : bit n (j + 1) = bit (n / 2) j := by
  unfold bit
  rw [Nat.div_div_eq_div_mul, ← pow_succ']

theorem bit_zero' (n : ℕ) : bit n 0 = n % 2 := by simp [bit]

/-- A number below `2^k` is the weighted sum of its `k` bits. -/
theorem sum_bits (k : ℕ) : ∀ n, n < 2 ^ k → ∑ j : Fin k, 2 ^ (j : ℕ) * bit n j = n := by
  induction k with
  | zero =>
    intro n hn
    have : n = 0 := by simpa using hn
    subst this
    simp
  | succ k ih =>
    intro n hn
    rw [Fin.sum_univ_succ]
    simp only [Fin.val_zero, pow_zero, one_mul, Fin.val_succ, bit_zero']
    have h1 : ∀ j : Fin k, 2 ^ ((j : ℕ) + 1) * bit n ((j : ℕ) + 1) =
        2 * (2 ^ (j : ℕ) * bit (n / 2) j) := by
      intro j
      rw [bit_succ, pow_succ]
      ring
    rw [Finset.sum_congr rfl (fun j _ => h1 j), ← Finset.mul_sum,
      ih (n / 2) (by rw [pow_succ] at hn; omega)]
    omega

/-- The weight is the number of set bits in a `k`-bit window. -/
theorem wt_eq_sum_bits (k : ℕ) : ∀ n, n < 2 ^ k → wt n = ∑ j : Fin k, bit n j := by
  induction k with
  | zero =>
    intro n hn
    have : n = 0 := by simpa using hn
    subst this
    simp
  | succ k ih =>
    intro n hn
    rw [Fin.sum_univ_succ, wt_eq n]
    simp only [Fin.val_zero, Fin.val_succ, bit_zero']
    rw [Finset.sum_congr rfl (fun (j : Fin k) _ => bit_succ n j),
      ih (n / 2) (by rw [pow_succ] at hn; omega)]

/-- Reindexing a sum over chain variables as a double sum over blocks and levels. -/
theorem sum_var (g : Fin m → ℕ → ℕ) :
    ∑ v : Fin (m * k), g (row v) (depth v) = ∑ i : Fin m, ∑ j : Fin k, g i j := by
  rw [← Equiv.sum_comp finProdFinEquiv, Fintype.sum_prod_type]
  simp only [row_var, depth_var]

theorem bitExp_digits (b : Fin m → ℕ) (hb : ∀ i, b i < 2 ^ k) : bitExp (digits k b) = b := by
  funext i
  rw [bitExp_apply]
  simp only [digits_apply, row_var, depth_var]
  exact sum_bits k (b i) (hb i)

theorem digits_injOn {b b' : Fin m → ℕ} (hb : ∀ i, b i < 2 ^ k) (hb' : ∀ i, b' i < 2 ^ k)
    (h : digits k b = digits k b') : b = b' := by
  rw [← bitExp_digits b hb, ← bitExp_digits b' hb', h]

theorem digits_degree (b : Fin m → ℕ) (hb : ∀ i, b i < 2 ^ k) :
    (digits k b).degree = ∑ i : Fin m, wt (b i) := by
  have sum : (digits k b).degree = ∑ v : Fin (m * k), digits k b v :=
    Finsupp.sum_fintype _ (fun _ n => n) (fun _ => rfl)
  rw [sum]
  simp only [digits_apply]
  rw [sum_var (fun i j => bit (b i) j)]
  apply Finset.sum_congr rfl
  intro i _
  exact (wt_eq_sum_bits k (b i) (hb i)).symm

theorem digits_zero_iff (b : Fin m → ℕ) (hb : ∀ i, b i < 2 ^ k) : digits k b = 0 ↔ b = 0 := by
  constructor
  · intro h
    apply digits_injOn hb (fun i => Nat.two_pow_pos k)
    rw [h]
    ext v
    simp [bit]
  · rintro rfl
    ext v
    simp [bit]

/-! ## Encoding the two-level exponents -/

/-- The exponent vector of a two-level exponent. -/
def encode (k0 : ℕ) (a : Ex m) : Fin m → ℕ := fun i => a.1 i + 2 ^ k0 * a.2 i

theorem two_pow_split (k0 k : ℕ) (hk : k0 ≤ k) : 2 ^ k0 * 2 ^ (k - k0) = 2 ^ k := by
  rw [← pow_add]
  congr 1
  omega

theorem encode_lt {k0 : ℕ} (hk : k0 ≤ k) {a : Ex m} (ha : InBox k0 k a) (i : Fin m) :
    encode k0 a i < 2 ^ k := by
  have hPN := two_pow_split k0 k hk
  have h1 : 2 ^ k0 * (a.2 i + 1) ≤ 2 ^ k0 * 2 ^ (k - k0) := Nat.mul_le_mul_left _ (ha.2 i)
  have h2 : 2 ^ k0 * (a.2 i + 1) = 2 ^ k0 * a.2 i + 2 ^ k0 := by ring
  have hl := ha.1 i
  unfold encode
  omega

theorem encode_injOn {k0 : ℕ} {a a' : Ex m} (ha : InBox k0 k a) (ha' : InBox k0 k a')
    (h : encode k0 a = encode k0 a') : a = a' := by
  have hP : 0 < 2 ^ k0 := Nat.two_pow_pos k0
  have key : ∀ i, a.1 i = a'.1 i ∧ a.2 i = a'.2 i := by
    intro i
    have hi := congrFun h i
    unfold encode at hi
    have e1 : (a.1 i + 2 ^ k0 * a.2 i) % 2 ^ k0 = a.1 i := by
      rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt (ha.1 i)]
    have e2 : (a'.1 i + 2 ^ k0 * a'.2 i) % 2 ^ k0 = a'.1 i := by
      rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt (ha'.1 i)]
    have d1 : (a.1 i + 2 ^ k0 * a.2 i) / 2 ^ k0 = a.2 i := by
      rw [mul_comm, Nat.add_mul_div_right _ _ hP, Nat.div_eq_of_lt (ha.1 i), zero_add]
    have d2 : (a'.1 i + 2 ^ k0 * a'.2 i) / 2 ^ k0 = a'.2 i := by
      rw [mul_comm, Nat.add_mul_div_right _ _ hP, Nat.div_eq_of_lt (ha'.1 i), zero_add]
    constructor
    · rw [← e1, ← e2, hi]
    · rw [← d1, ← d2, hi]
  exact Prod.ext (funext fun i => (key i).1) (funext fun i => (key i).2)

theorem encode_zero (k0 : ℕ) : encode k0 ((0 : Fin m → ℕ), (0 : Fin m → ℕ)) = 0 := by
  funext i
  simp [encode]

/-- Model soundness (Lemma 4.2): the model step is addition of `2^j` in block `i` modulo `2^k`. -/
theorem encode_step {k0 : ℕ} {i : Fin m} {j : ℕ} (hk : k0 < k) {a : Ex m}
    (ha : InBox k0 k a) :
    encode k0 (step k0 k i j a) =
      Function.update (encode k0 a) i ((encode k0 a i + 2 ^ j) % 2 ^ k) := by
  obtain ⟨l, q⟩ := a
  have hl : l i < 2 ^ k0 := ha.1 i
  have hq : q i < 2 ^ (k - k0) := ha.2 i
  have hPN := two_pow_split k0 k hk.le
  have hP : 0 < 2 ^ k0 := Nat.two_pow_pos k0
  have hN : 0 < 2 ^ (k - k0) := Nat.two_pow_pos _
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos j
  funext i'
  rcases eq_or_ne i' i with rfl | hne
  · rw [Function.update_self]
    unfold step
    split_ifs with hjl
    · have hjk : 2 ^ j < 2 ^ k0 := Nat.pow_lt_pow_right (by norm_num) hjl
      unfold stepLow
      split_ifs with hc
      · have hc' : l i' + 2 ^ j < 2 ^ k0 := hc
        simp only [encode, Function.update_self]
        have h1 : 2 ^ k0 * (q i' + 1) ≤ 2 ^ k0 * 2 ^ (k - k0) := Nat.mul_le_mul_left _ hq
        have h2 : 2 ^ k0 * (q i' + 1) = 2 ^ k0 * q i' + 2 ^ k0 := by ring
        rw [Nat.mod_eq_of_lt (by omega)]
        ring
      · have hc' : ¬ l i' + 2 ^ j < 2 ^ k0 := hc
        simp only [encode, Function.update_self]
        by_cases hq1 : q i' + 1 < 2 ^ (k - k0)
        · rw [Nat.mod_eq_of_lt hq1]
          have h1 : 2 ^ k0 * (q i' + 2) ≤ 2 ^ k0 * 2 ^ (k - k0) := Nat.mul_le_mul_left _ (by omega)
          have h2 : 2 ^ k0 * (q i' + 2) = 2 ^ k0 * q i' + 2 ^ k0 + 2 ^ k0 := by ring
          have h3 : 2 ^ k0 * (q i' + 1) = 2 ^ k0 * q i' + 2 ^ k0 := by ring
          rw [Nat.mod_eq_of_lt (by omega)]
          omega
        · have hqN : q i' + 1 = 2 ^ (k - k0) := by omega
          rw [hqN, Nat.mod_self, mul_zero, add_zero]
          have h3 : 2 ^ k0 * (q i' + 1) = 2 ^ k0 * q i' + 2 ^ k0 := by ring
          have hsum : l i' + 2 ^ k0 * q i' + 2 ^ j = (l i' + 2 ^ j - 2 ^ k0) + 2 ^ k := by
            rw [← hPN, ← hqN]
            omega
          rw [hsum, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
    · unfold stepHigh
      simp only [encode, Function.update_self]
      have h2j : 2 ^ j = 2 ^ k0 * 2 ^ (j - k0) := by
        rw [← pow_add]
        congr 1
        omega
      have key : ∀ x, (l i' + 2 ^ k0 * x) % 2 ^ k = l i' + 2 ^ k0 * (x % 2 ^ (k - k0)) := by
        intro x
        have hx : x = 2 ^ (k - k0) * (x / 2 ^ (k - k0)) + x % 2 ^ (k - k0) :=
          (Nat.div_add_mod _ _).symm
        have hsplit : l i' + 2 ^ k0 * x =
            (l i' + 2 ^ k0 * (x % 2 ^ (k - k0))) + 2 ^ k * (x / 2 ^ (k - k0)) := by
          rw [← hPN]
          conv_lhs => rw [hx]
          ring
        rw [hsplit, Nat.add_mul_mod_self_left]
        have hm : x % 2 ^ (k - k0) < 2 ^ (k - k0) := Nat.mod_lt _ hN
        have h1 : 2 ^ k0 * (x % 2 ^ (k - k0) + 1) ≤ 2 ^ k0 * 2 ^ (k - k0) :=
          Nat.mul_le_mul_left _ hm
        have h2 : 2 ^ k0 * (x % 2 ^ (k - k0) + 1) = 2 ^ k0 * (x % 2 ^ (k - k0)) + 2 ^ k0 := by
          ring
        rw [Nat.mod_eq_of_lt (by omega)]
      rw [h2j, show l i' + 2 ^ k0 * q i' + 2 ^ k0 * 2 ^ (j - k0) =
        l i' + 2 ^ k0 * (q i' + 2 ^ (j - k0)) by ring, key]
  · rw [Function.update_of_ne hne]
    unfold step stepLow stepHigh
    split_ifs <;> simp [encode, Function.update_of_ne hne]

/-! ## The model map -/

/-- The model map: `δ_(l,q) ↦` the squarefree chain monomial of `l + 2^k0 q`. -/
def Φ (k0 k : ℕ) : Al m K →ₗ[K] FiniteQuotient.P (m * k) K :=
  Finsupp.lmapDomain K K (fun a : Ex m => digits k (encode k0 a))

variable {k0 : ℕ}

theorem Φ_single (a : Ex m) (c : K) :
    Φ (K := K) k0 k (Finsupp.single a c) = monomial (digits k (encode k0 a)) c := by
  show Finsupp.mapDomain _ (Finsupp.single a c) = _
  rw [Finsupp.mapDomain_single, single_eq_monomial]

theorem exponent_injOn (hk : k0 ≤ k) :
    Set.InjOn (fun a : Ex m => digits k (encode k0 a)) {a | InBox k0 k a} := by
  intro a ha a' ha' h
  exact encode_injOn ha ha' (digits_injOn (encode_lt hk ha) (encode_lt hk ha') h)

theorem Φ_coeff (hk : k0 ≤ k) {w : Al m K} (hw : Boxed k0 k w) {a : Ex m} (ha : InBox k0 k a) :
    coeff (digits k (encode k0 a)) (Φ k0 k w) = w a := by
  change Finsupp.mapDomain (fun a : Ex m => digits k (encode k0 a)) w (digits k (encode k0 a)) = w a
  exact Finsupp.mapDomain_apply' {a : Ex m | InBox k0 k a} w (fun a ha => hw a ha)
    (exponent_injOn hk) ha

theorem inBox_zero (k0 k : ℕ) : InBox k0 k ((0 : Fin m → ℕ), (0 : Fin m → ℕ)) :=
  ⟨lowBox_zero _, lowBox_zero _⟩

/-- `Φ` reflects `1` on boxed elements. -/
theorem Φ_one (hk : k0 ≤ k) {w : Al m K} (hw : Boxed k0 k w)
    (h : Φ k0 k w = (1 : FiniteQuotient.P (m * k) K)) :
    w = Finsupp.single ((0 : Fin m → ℕ), (0 : Fin m → ℕ)) 1 := by
  classical
  ext a
  rw [Finsupp.single_apply]
  by_cases ha : InBox k0 k a
  · rw [← Φ_coeff hk hw ha, h, coeff_one]
    have hiff : (0 : Fin (m * k) →₀ ℕ) = digits k (encode k0 a) ↔ (0, 0) = a := by
      constructor
      · intro h0
        have hz : digits k (encode k0 a) = 0 := h0.symm
        rw [digits_zero_iff _ (encode_lt hk ha)] at hz
        rw [← encode_zero k0] at hz
        exact encode_injOn (inBox_zero k0 k) ha hz.symm
      · rintro rfl
        rw [encode_zero]
        ext v
        simp [bit]
    by_cases h0 : (0 : Fin (m * k) →₀ ℕ) = digits k (encode k0 a)
    · rw [if_pos h0, if_pos (hiff.mp h0)]
    · rw [if_neg h0, if_neg (fun h' => h0 (hiff.mpr h'))]
  · have hz : w a = 0 := Finsupp.not_mem_support_iff.mp (fun h' => ha (hw a h'))
    rw [hz, if_neg]
    rintro rfl
    exact ha (inBox_zero k0 k)

theorem pwt_eq_degree (hk : k0 ≤ k) {a : Ex m} (ha : InBox k0 k a) :
    pwt a = (digits k (encode k0 a)).degree := by
  rw [digits_degree _ (encode_lt hk ha)]
  unfold pwt lowWt encode
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i _
  exact (wt_add_mul_two_pow k0 _ _ (ha.1 i)).symm

/-- `Φ` reflects the degree filtration on boxed elements. -/
theorem Φ_degree (hk : k0 ≤ k) {w : Al m K} (hw : Boxed k0 k w) (D : ℕ)
    (h : (Φ k0 k w).totalDegree < D) : WtLT D w := by
  intro a ha
  have hmem : digits k (encode k0 a) ∈ (Φ k0 k w).support := by
    rw [mem_support_iff, Φ_coeff hk hw (hw a ha)]
    exact Finsupp.mem_support_iff.mp ha
  have hle := le_totalDegree hmem
  have hdeg : (digits k (encode k0 a)).degree = (digits k (encode k0 a)).sum (fun _ e => e) :=
    rfl
  rw [pwt_eq_degree hk (hw a ha), hdeg]
  omega

theorem translate_single (i : Fin m) (j : ℕ) (a : Ex m) (c : K) :
    translate (K := K) k0 k i j (Finsupp.single a c) = Finsupp.single (step k0 k i j a) c := by
  show Finsupp.mapDomain _ (Finsupp.single a c) = _
  rw [Finsupp.mapDomain_single]

theorem linearMap_ext_support {M : Type*} [AddCommMonoid M] [Module K M] (w : Al m K)
    (L₁ L₂ : Al m K →ₗ[K] M)
    (h : ∀ a ∈ w.support, L₁ (Finsupp.single a (w a)) = L₂ (Finsupp.single a (w a))) :
    L₁ w = L₂ w := by
  rw [← Finsupp.sum_single w]
  unfold Finsupp.sum
  rw [map_sum, map_sum]
  exact Finset.sum_congr rfl h

/-- `Φ` intertwines the model translation with reduced multiplication by the chain variable. -/
theorem Φ_step (two_ne : (2 : K) ≠ 0) (hk : k0 < k) {w : Al m K} (hw : Boxed k0 k w) (i : Fin m)
    (j : Fin k) :
    ChainClosureBridge.normalForm (chainInputs K m k) (chainInputs_nonzero m k)
        (X (var m k i j) * Φ k0 k w) =
      Φ k0 k (translate k0 k i (j : ℕ) w) := by
  classical
  have zero := reduced_chain_ideal_trivial m k two_ne
  let L₁ : Al m K →ₗ[K] FiniteQuotient.P (m * k) K :=
    (ChainClosureBridge.nfMap (chainInputs K m k) (chainInputs_nonzero m k) zero).comp
      ((LinearMap.mulLeft K (X (var m k i j))).comp (Φ k0 k))
  let L₂ : Al m K →ₗ[K] FiniteQuotient.P (m * k) K := (Φ k0 k).comp (translate k0 k i (j : ℕ))
  change L₁ w = L₂ w
  apply linearMap_ext_support w L₁ L₂
  intro a ha
  have hbox := hw a ha
  simp only [L₁, L₂, LinearMap.comp_apply, LinearMap.mulLeft_apply,
    ChainClosureBridge.nfMap_apply]
  have hm : ∀ e : Fin (m * k) →₀ ℕ, monomial e (w a) = C (w a) * monomial e (1 : K) := by
    intro e
    rw [C_mul_monomial, mul_one]
  rw [translate_single, Φ_single, Φ_single, hm, hm, mul_left_comm, ← smul_eq_C_mul,
    ChainClosureBridge.normalForm_smul _ _ zero,
    normalForm_X_mul_monomial two_ne i j _ (encode_lt hk.le hbox), ← encode_step hk hbox,
    smul_eq_C_mul]

/-! ## N5 assembled with the concrete model map -/

/-- `Φ 7 k` is a model map for the two-summand presentation at every depth `k ≥ 8`. -/
theorem modelMap (two_ne : (2 : K) ≠ 0) (hk : k0 < k) (D : ℕ) :
    ModelMap k0 k D (Φ (m := 2) (K := K) k0 k) where
  one := fun _ hw h => Φ_one hk.le hw h
  degree := fun _ hw h => Φ_degree hk.le hw D h
  step := fun _ hw i j => Φ_step two_ne hk hw i j

/-! ## The normalized generator is reduced -/

theorem reduced_of_le_one (p : FiniteQuotient.P (m * k) K)
    (h : ∀ e ∈ p.support, ∀ v, e v ≤ 1) : p ∈ reducedSpace (chainInputs K m k) := by
  intro e he g hg
  obtain ⟨v, rfl⟩ := hg
  rw [leading_chainPolynomial, Finsupp.single_le_iff]
  have := h e he v
  omega

theorem C_reduced (c : K) : (C c : FiniteQuotient.P (m * k) K) ∈ reducedSpace (chainInputs K m k) := by
  classical
  apply reduced_of_le_one
  intro e he v
  rw [C_apply, support_monomial] at he
  split_ifs at he with hc
  · simp at he
  · rw [Finset.mem_singleton] at he
    simp [he]

theorem X_reduced (i : Fin (m * k)) :
    (X i : FiniteQuotient.P (m * k) K) ∈ reducedSpace (chainInputs K m k) := by
  classical
  apply reduced_of_le_one
  intro e he v
  rw [support_X, Finset.mem_singleton] at he
  subst he
  rw [Finsupp.single_apply]
  split_ifs <;> omega

theorem X_mul_X_reduced {i j : Fin (m * k)} (hij : i ≠ j) :
    (X i * X j : FiniteQuotient.P (m * k) K) ∈ reducedSpace (chainInputs K m k) := by
  classical
  apply reduced_of_le_one
  intro e he v
  rw [X, X, monomial_mul, one_mul, support_monomial, if_neg one_ne_zero, Finset.mem_singleton] at he
  subst he
  rw [Finsupp.add_apply, Finsupp.single_apply, Finsupp.single_apply]
  split_ifs <;> omega

theorem C_mul_reduced (c : K) {p : FiniteQuotient.P (m * k) K}
    (hp : p ∈ reducedSpace (chainInputs K m k)) :
    C c * p ∈ reducedSpace (chainInputs K m k) := by
  rw [← smul_eq_C_mul]
  exact (reducedSpace _).smul_mem c hp

open PrimeM2General GeneralQuadraticTop in
/-- The normalized two-summand generator is a reduced (squarefree) polynomial. -/
theorem normalized_reduced (k : ℕ) (hk2 : 2 ≤ k) (a b t : K) :
    normalized k hk2 a b t ∈ reducedSpace (chainInputs K 2 k) := by
  have inj := active_injective k hk2
  have ne : ∀ i j : Fin 4, i ≠ j → active k hk2 i ≠ active k hk2 j := fun i j h => inj.ne h
  unfold normalized top tail
  refine (reducedSpace _).add_mem ((reducedSpace _).add_mem ((reducedSpace _).add_mem ?_ ?_) ?_)
    ((reducedSpace _).add_mem ((reducedSpace _).add_mem ?_ ?_) ?_)
  · exact C_mul_reduced _ (X_mul_X_reduced (ne 0 2 (by decide)))
  · exact C_mul_reduced _ ((reducedSpace _).add_mem (X_mul_X_reduced (ne 0 3 (by decide)))
      (X_mul_X_reduced (ne 1 2 (by decide))))
  · exact X_mul_X_reduced (ne 1 3 (by decide))
  · exact C_mul_reduced _ ((reducedSpace _).add_mem (X_reduced _) (X_reduced _))
  · exact C_mul_reduced _ ((reducedSpace _).add_mem (X_reduced _) (X_reduced _))
  · exact C_reduced _

open PrimeM2General in
theorem normalForm_normalized (two_ne : (2 : K) ≠ 0) (k : ℕ) (hk2 : 2 ≤ k) (a b t : K) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (normalized k hk2 a b t) = normalized k hk2 a b t :=
  ChainClosureBridge.normalForm_of_reduced _ _ (reduced_chain_ideal_trivial 2 k two_ne)
    (normalized_reduced k hk2 a b t)

/-! ## N5 assembled with the concrete model map -/

open PrimeM2General in
/-- The generator lies in the base: a low element `v ∈ S` in the box of weight `≤ D` whose
monomial image is the normalized generator. -/
def GeneratorInBase' (k0 D : ℕ) (S : Submodule K (Lo 2 K)) (k : ℕ) (hk2 : 2 ≤ k) (a b t : K) :
    Prop :=
  (normalized k hk2 a b t).totalDegree ≤ D →
    ∃ v ∈ S, (∀ l ∈ v.support, LowBox k0 l) ∧ (∀ l ∈ v.support, lowWt l ≤ D) ∧
      Φ k0 k (embed 0 v) = normalized k hk2 a b t

open PrimeM2General in
theorem generatorInBase_of (two_ne : (2 : K) ≠ 0) {k0 D : ℕ} {S : Submodule K (Lo 2 K)} {k : ℕ}
    {hk2 : 2 ≤ k} {a b t : K} (gen : GeneratorInBase' k0 D S k hk2 a b t) :
    GeneratorInBase k0 D S k hk2 a b t (Φ k0 k) := by
  intro hd
  obtain ⟨v, hv, h1, h2, h3⟩ := gen hd
  exact ⟨v, hv, h1, h2, h3.trans (normalForm_normalized two_ne k hk2 a b t).symm⟩

open PrimeM2General in
/-- **N5, assembled with the concrete model map.** Base at depth `k0` (`1 ∉ S`, the free
extension closed at depth `k0+1`, generator in the base): `1` is not in the degree-`D` closure of
the two-summand presentation at any depth `k ≥ k0+1`. -/
theorem not_one_of_base (k0 D : ℕ) (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (gen : GeneratorInBase' k0 D S k hk2 a b t) :
    ¬ InClosure (Set.range (inputs k hk2 a b t)) D (1 : FiniteQuotient.P (2*k) K) :=
  DepthReductionBridge.not_one_of_base k0 D S base one k hk hk2 two_ne a b t (Φ k0 k)
    (modelMap two_ne (by omega) D) (generatorInBase_of two_ne gen)

open PrimeM2General in
theorem lfd_gt_of_base (k0 D : ℕ) (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (gen : GeneratorInBase' k0 D S k hk2 a b t)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e :=
  DepthReductionBridge.lfd_gt_of_base k0 D S base one k hk hk2 two_ne a b t (Φ k0 k)
    (modelMap two_ne (by omega) D) (generatorInBase_of two_ne gen) unit

open PrimeM2General in
/-- N5 with the depth-7 base and degree 4: a unit-ideal two-summand instance whose base
excludes `1`, contains the generator, and whose depth-8 free extension is closed has
`LFD ≥ 5` at every depth `k ≥ 8`. -/
theorem n5_lfd_ge_five (S : Submodule K (Lo 2 K))
    (base : Closed 7 8 4 (W 7 4 S 8))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (gen : GeneratorInBase' 7 4 S k hk2 a b t)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  intro e he
  exact lfd_gt_of_base 7 4 S base one k hk hk2 two_ne a b t gen unit e (by omega)

open PrimeM2General in
/-- Prime-field entry point. -/
theorem prime_field_n5_lfd_ge_five (p : ℕ) [Fact p.Prime] (large : 5 ≤ p)
    (S : Submodule (ZMod p) (Lo 2 (ZMod p)))
    (base : Closed 7 8 4 (W 7 4 S 8))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : ZMod p) ∉ S)
    (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) (a b t : ZMod p)
    (gen : GeneratorInBase' 7 4 S k hk2 a b t)
    (unit : (1 : FiniteQuotient.P (2*k) (ZMod p)) ∈
      Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e :=
  n5_lfd_ge_five S base one k hk hk2 (PrimeFieldConclusions.two_ne_zero p large) a b t gen unit

end ICFallDegrees.ChainModelMap
