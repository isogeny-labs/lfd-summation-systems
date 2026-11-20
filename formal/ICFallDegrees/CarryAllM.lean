import ICFallDegrees.M3CarryLower

noncomputable section

/-!
# The carry lower bound for every number of summands

Generalizes `M3CarryLower` from three summands to every `m ≥ 3`, for an
arbitrary generator `S` in the summand coordinates `x_1, …, x_m` (the
summation polynomial `S_{m+1}(x_1, …, x_m, t)` is the intended instance).

* `lift m k S` is the chain normal form of `S(u_{1,0}, …, u_{m,0})`: the monomial
  `x^e` goes to the squarefree chain monomial whose digit set in block `i` is the
  set of binary digits `j < k` of `e i` (`lift_reduced`, `lift_congr`; the
  congruence holds for every exponent, because `u_{i,0}^(2^k) ≡ 1`).
* **Onset bound** (`lfd_ge_onset`). Let `m ≥ 3`, `k ≥ m`, and let `S` have degree at
  most `2^(m-1)` in each variable with a nonzero coefficient at
  `x_1^(2^(m-1)-1) ⋯ x_m^(2^(m-1)-1)`. Then the chain normal form `g` has total
  degree exactly `m(m-1)`, the carried row `h = NF(u_{1,0} g)` has degree at most
  `m(m-1)` (an exponent `e ≤ 2^(m-1)+1 < 2^m - 1` has at most `m-1` binary
  digits), and `h` is not a multiple of `g` (the lift is injective below `2^k`).
  Hence `LFD(chains ∪ {g}) ≥ m(m-1)+1` (`CarryLower.chain_carry_lower`).
* **Full-degree bound** (`lfd_ge_full`). For every `m, k ≥ 1` and every `S`: if the
  chain normal form has the full degree `m k` and `NF(u_{i,0} g)` is not a multiple
  of `g` for one block `i`, then `LFD ≥ m k + 1` (with `ChainQuotient.reduced_inputs_lfd_bound`
  this is exact).
* `m3_onset`: the three-summand instance `S = S_4(x1,x2,x3,t)` of the onset bound,
  from `M3CarryLower.s4_support` and `M3CarryLower.coeff_s4_333`.

Reuse (lean-instance-gate): the lower bound itself is `CarryLower.chain_carry_lower` used as
stated; reducedness and uniqueness come from `ChainQuotient.leading_chainPolynomial` and
`ChainQuotient.reduced_chain_ideal_trivial` (inside `chain_carry_lower`). New here: the digit
weights (`weight_*`, Mathlib has no digit-sum bound of this form for `Nat.testBit`), and the
`m`-block versions of the `M3CarryLower` lift lemmas (`maskExp`, `lift`, `lift_congr`), whose
`m = 3` originals are left unchanged; `lift_congr` now holds for every exponent (the
`M3CarryLower` version needs exponents `≤ 5`).
-/

namespace ICFallDegrees.CarryAllM

open MvPolynomial FiniteQuotient ChainQuotient

variable {K : Type*} [Field K]

/-- Polynomials in the `m` summand coordinates. -/
abbrev Z (m : ℕ) (K : Type*) [Field K] := MvPolynomial (Fin m) K

/-! ## Binary digit weights -/

/-- The number of binary digits of `x` among the lowest `k`. -/
def weight (k x : ℕ) : ℕ := ∑ j : Fin k, if x.testBit j then 1 else 0

theorem weight_eq_range (k x : ℕ) :
    weight k x=∑ j ∈ Finset.range k, if x.testBit j then 1 else 0 :=
  Fin.sum_univ_eq_sum_range (fun j => if x.testBit j then 1 else 0) k

/-- Digits at or above `n` of a number below `2^n` vanish. -/
theorem weight_eq_low {n k x : ℕ} (hx : x<2^n) (hk : n≤k) :
    weight k x=∑ j ∈ Finset.range n, if x.testBit j then 1 else 0 := by
  rw [weight_eq_range,← Finset.sum_range_add_sum_Ico _ hk]
  have high : ∑ j ∈ Finset.Ico n k, (if x.testBit j then 1 else 0)=0 := by
    apply Finset.sum_eq_zero
    intro j hj
    have hj' : n≤j := (Finset.mem_Ico.mp hj).1
    have large : x<2^j := lt_of_lt_of_le hx (Nat.pow_le_pow_right (by norm_num) hj')
    rw [Nat.testBit_lt_two_pow large]
    rfl
  rw [high,add_zero]

/-- `2^n - 1` has exactly `n` binary digits. -/
theorem weight_two_pow_sub_one {n k : ℕ} (hk : n≤k) : weight k (2^n-1)=n := by
  have hx : 2^n-1<2^n := Nat.sub_lt (Nat.two_pow_pos n) (by norm_num)
  rw [weight_eq_low hx hk]
  have ones : ∀ j ∈ Finset.range n, (if (2^n-1).testBit j then 1 else 0)=1 := by
    intro j hj
    rw [Nat.testBit_two_pow_sub_one]
    simp [Finset.mem_range.mp hj]
  rw [Finset.sum_congr rfl ones]
  simp

/-- A number below `2^n` other than `2^n - 1` has at most `n - 1` binary digits. -/
theorem weight_le_pred {n k x : ℕ} (hx : x<2^n) (hne : x≠2^n-1) (hk : n≤k) :
    weight k x≤n-1 := by
  rw [weight_eq_low hx hk]
  -- Some digit below `n` is zero.
  have gap : ∃ j<n, x.testBit j=false := by
    by_contra full
    push_neg at full
    apply hne
    apply Nat.eq_of_testBit_eq
    intro j
    rw [Nat.testBit_two_pow_sub_one]
    by_cases hj : j<n
    · have := full j hj
      cases h : x.testBit j
      · exact absurd h this
      · simp [hj]
    · have large : x<2^j := lt_of_lt_of_le hx (Nat.pow_le_pow_right (by norm_num) (by omega))
      rw [Nat.testBit_lt_two_pow large]
      simp [hj]
  obtain ⟨j0,hj0,zero⟩ := gap
  have mem : j0∈Finset.range n := Finset.mem_range.mpr hj0
  rw [← Finset.add_sum_erase _ _ mem,zero]
  simp only [Bool.false_eq_true,↓reduceIte,zero_add]
  calc ∑ j ∈ (Finset.range n).erase j0, (if x.testBit j then 1 else 0)
        ≤∑ _j ∈ (Finset.range n).erase j0, 1 := by
          apply Finset.sum_le_sum
          intro j _
          split <;> omega
    _ = n-1 := by
          rw [Finset.sum_const,Finset.card_erase_of_mem mem,Finset.card_range,smul_eq_mul,
            mul_one]

theorem weight_le (k x : ℕ) : weight k x≤k := by
  unfold weight
  calc ∑ j : Fin k, (if x.testBit j then 1 else 0)≤∑ _j : Fin k, 1 := by
        apply Finset.sum_le_sum
        intro j _
        split <;> omega
    _ = k := by simp

/-! ## The chain exponent of a monomial -/

/-- The chain exponent of `x^e`: bit `j` of `e i` at the chain variable `u_{i,j}`. -/
def maskExp (m k : ℕ) (e : Fin m →₀ ℕ) : Fin (m*k) →₀ ℕ :=
  Finsupp.equivFunOnFinite.symm fun v => if (e (row v)).testBit (depth v) then 1 else 0

theorem maskExp_apply (m k : ℕ) (e : Fin m →₀ ℕ) (v : Fin (m*k)) :
    maskExp m k e v=if (e (row v)).testBit (depth v) then 1 else 0 := rfl

theorem maskExp_var (m k : ℕ) (e : Fin m →₀ ℕ) (i : Fin m) (j : Fin k) :
    maskExp m k e (var m k i j)=if (e i).testBit j then 1 else 0 := by
  rw [maskExp_apply,row_var,depth_var]

theorem maskExp_equiv (m k : ℕ) (e : Fin m →₀ ℕ) (i : Fin m) (j : Fin k) :
    maskExp m k e (finProdFinEquiv (i,j))=if (e i).testBit j then 1 else 0 :=
  maskExp_var m k e i j

theorem maskExp_le_one (m k : ℕ) (e : Fin m →₀ ℕ) (v : Fin (m*k)) : maskExp m k e v≤1 := by
  rw [maskExp_apply]; split <;> omega

theorem maskExp_injective {m k : ℕ} {e e' : Fin m →₀ ℕ} (he : ∀ i, e i<2^k)
    (he' : ∀ i, e' i<2^k) (h : maskExp m k e=maskExp m k e') : e=e' := by
  ext i
  apply Nat.eq_of_testBit_eq
  intro j
  by_cases hj : j<k
  · have := congrArg (fun s : Fin (m*k) →₀ ℕ => s (var m k i ⟨j,hj⟩)) h
    simp only [maskExp_var] at this
    cases h1 : (e i).testBit j <;> cases h2 : (e' i).testBit j <;> simp_all
  · have large : 2^k≤2^j := Nat.pow_le_pow_right (by norm_num) (by omega)
    rw [Nat.testBit_lt_two_pow ((he i).trans_le large),
      Nat.testBit_lt_two_pow ((he' i).trans_le large)]

/-- The degree of a chain exponent is the sum of the digit weights of its blocks. -/
theorem maskExp_sum (m k : ℕ) (e : Fin m →₀ ℕ) :
    ((maskExp m k e).sum fun _ n => n)=∑ i : Fin m, weight k (e i) := by
  rw [Finsupp.sum_fintype _ _ (fun _ => rfl)]
  rw [← Fintype.sum_equiv finProdFinEquiv (fun ij => maskExp m k e (finProdFinEquiv ij)) _
    (fun _ => rfl)]
  rw [Fintype.sum_prod_type]
  simp only [maskExp_equiv]
  rfl

/-! ## The lift -/

/-- The chain normal form of a polynomial in the summand coordinates. -/
def lift (m k : ℕ) (p : Z m K) : P (m*k) K :=
  ∑ e∈p.support, monomial (maskExp m k e) (coeff e p)

theorem mem_support_lift {m k : ℕ} {p : Z m K} {s : Fin (m*k) →₀ ℕ}
    (hs : s∈(lift m k p).support) : ∃ e∈p.support, s=maskExp m k e := by
  classical
  obtain ⟨e,he,hs⟩ := Finset.mem_biUnion.mp (support_sum hs)
  exact ⟨e,he,Finset.mem_singleton.mp (support_monomial_subset hs)⟩

theorem lift_reduced (m k : ℕ) (p : Z m K) : lift m k p∈reducedSpace (chainInputs K m k) := by
  intro s hs g hg
  obtain ⟨e,-,rfl⟩ := mem_support_lift hs
  obtain ⟨v,rfl⟩ := hg
  rw [leading_chainPolynomial,Finsupp.single_le_iff]
  have := maskExp_le_one m k e v
  omega

theorem coeff_lift {m k : ℕ} {p : Z m K} {e : Fin m →₀ ℕ}
    (hp : ∀ s∈p.support, ∀ i, s i<2^k) (he : ∀ i, e i<2^k) :
    coeff (maskExp m k e) (lift m k p)=coeff e p := by
  classical
  unfold lift
  rw [coeff_sum]
  simp only [coeff_monomial]
  by_cases hmem : e∈p.support
  · rw [Finset.sum_eq_single e]
    · simp
    · intro s hs hne
      rw [if_neg]
      intro heq
      exact hne (maskExp_injective (hp s hs) he heq)
    · intro h
      exact absurd hmem h
  · rw [Finset.sum_eq_zero,eq_comm]
    · exact not_mem_support_iff.mp hmem
    · intro s hs
      rw [if_neg]
      intro heq
      exact hmem (maskExp_injective (hp s hs) he heq ▸ hs)

/-- Degree bound for a lift: the sum of per-block digit-weight bounds. -/
theorem lift_totalDegree_le (m k w : ℕ) (p : Z m K)
    (hp : ∀ s∈p.support, ∀ i, weight k (s i)≤w) : (lift m k p).totalDegree≤m*w := by
  apply (totalDegree_finset_sum _ _).trans
  apply Finset.sup_le
  intro e he
  apply (totalDegree_monomial_le _ _).trans
  have key : ((maskExp m k e).sum fun _ => id)=∑ i : Fin m, weight k (e i) :=
    maskExp_sum m k e
  rw [key]
  calc ∑ i : Fin m, weight k (e i)≤∑ _i : Fin m, w :=
        Finset.sum_le_sum (fun i _ => hp e he i)
    _ = m*w := by simp

/-! ## The lift is congruent to the substitution `x_i ↦ u_{i,0}` -/

/-- The chain substitution `x_i ↦ u_{i,0}`. -/
def toChain (m k : ℕ) (hk : 0<k) : Z m K →ₐ[K] P (m*k) K :=
  rename fun i => var m k i ⟨0,hk⟩

/-- The class of a chain variable in the quotient by the chains (`1` past the end). -/
def chainClass (m k : ℕ) (i : Fin m) (j : ℕ) :
    P (m*k) K ⧸ Ideal.span (chainInputs K m k) :=
  if h : j<k then Ideal.Quotient.mk _ (X (var m k i ⟨j,h⟩)) else 1

theorem chainClass_square (m k : ℕ) (i : Fin m) (j : ℕ) (hj : j+1<k) :
    (chainClass (K := K) m k i j)^2=chainClass m k i (j+1) := by
  have hj' : j<k := by omega
  have mem : chainPolynomial K m k (var m k i ⟨j,hj'⟩)∈Ideal.span (chainInputs K m k) :=
    Ideal.subset_span ⟨_,rfl⟩
  have tail : chainTail K m k (var m k i ⟨j,hj'⟩)=X (var m k i ⟨j+1,hj⟩) := by
    unfold chainTail
    split_ifs with hlt
    · congr 2
      · exact row_var m k i _
      · apply Fin.ext
        simp only [depth_var]
    · exfalso
      simp only [depth_var] at hlt
      omega
  unfold chainClass
  rw [dif_pos hj',dif_pos hj,← map_pow,Ideal.Quotient.eq]
  simpa [chainPolynomial,tail] using mem

theorem chainClass_last (m k : ℕ) (i : Fin m) (hk : 0<k) :
    (chainClass (K := K) m k i (k-1))^2=1 := by
  have hj : k-1<k := by omega
  have mem : chainPolynomial K m k (var m k i ⟨k-1,hj⟩)∈Ideal.span (chainInputs K m k) :=
    Ideal.subset_span ⟨_,rfl⟩
  have tail : chainTail K m k (var m k i ⟨k-1,hj⟩)=1 := by
    unfold chainTail
    split_ifs with hlt
    · exfalso
      simp only [depth_var] at hlt
      omega
    · rfl
  unfold chainClass
  rw [dif_pos hj,← map_pow,← map_one (Ideal.Quotient.mk (Ideal.span (chainInputs K m k))),
    Ideal.Quotient.eq]
  simpa [chainPolynomial,tail] using mem

/-- Iterated squaring along a chain. -/
theorem chainClass_pow (m k : ℕ) (i : Fin m) :
    ∀ j, j<k → (chainClass (K := K) m k i 0)^(2^j)=chainClass m k i j := by
  intro j
  induction j with
  | zero => intro _; simp
  | succ j ih =>
    intro hj
    rw [pow_succ (2 : ℕ) j,pow_mul,ih (by omega),chainClass_square m k i j hj]

theorem chainClass_period (m k : ℕ) (i : Fin m) (hk : 0<k) :
    (chainClass (K := K) m k i 0)^(2^k)=1 := by
  obtain ⟨n,rfl⟩ : ∃ n, k=n+1 := ⟨k-1,by omega⟩
  rw [pow_succ (2 : ℕ) n,pow_mul,chainClass_pow m (n+1) i n (by omega)]
  have := chainClass_last (K := K) m (n+1) i (by omega)
  simpa using this

/-- Binary expansion in a ring: if `y (j+1) = y j ^ 2` below `k`, then the product of the
`y j` over the binary digits `j < k` of `x < 2^k` is `y 0 ^ x`. -/
theorem prod_digits {R : Type*} [CommRing R] :
    ∀ (k : ℕ) (y : ℕ → R), (∀ j, j+1<k → y j^2=y (j+1)) → ∀ x, x<2^k →
      ∏ j : Fin k, y j^(if x.testBit j then 1 else 0)=y 0^x := by
  intro k
  induction k with
  | zero =>
    intro y _ x hx
    have : x=0 := by simpa using hx
    subst this
    simp
  | succ k ih =>
    intro y hy x hx
    rw [Fin.prod_univ_succ]
    have shift : ∏ j : Fin k, y ((j : ℕ)+1)^(if (x/2).testBit j then 1 else 0)=y 1^(x/2) :=
      ih (fun j => y (j+1)) (fun j hj => hy (j+1) (by omega)) (x/2)
        (by rw [pow_succ] at hx; omega)
    simp only [Fin.val_succ,Fin.val_zero]
    have digits : ∀ j : Fin k, (if x.testBit (j+1) then 1 else 0)=
        (if (x/2).testBit j then 1 else 0) := by
      intro j
      rw [Nat.testBit_div_two]
    simp only [digits]
    rw [shift]
    have bit0 : (if x.testBit 0 then 1 else 0)=x%2 := by
      rw [Nat.testBit_zero]
      rcases Nat.mod_two_eq_zero_or_one x with h | h <;> simp [h]
    rw [bit0]
    by_cases half : x/2=0
    · rw [half,pow_zero,mul_one]
      congr 1
      omega
    · have hpos : 0<k := by
        rcases Nat.eq_zero_or_pos k with h0 | h0
        · subst h0
          have : x<2 := by simpa using hx
          omega
        · exact h0
      have h1 : y 1=y 0^2 := (hy 0 (by omega)).symm
      rw [h1,← pow_mul,← pow_add]
      congr 1
      omega

/-- In the quotient by the chains, block `i` of a chain monomial is `u_{i,0}^x`. -/
theorem block (m k : ℕ) (hk : 0<k) (i : Fin m) (x : ℕ) :
    ∏ j : Fin k, Ideal.Quotient.mk (Ideal.span (chainInputs K m k))
        (X (var m k i j))^(if x.testBit j then 1 else 0)=
      Ideal.Quotient.mk (Ideal.span (chainInputs K m k)) (X (var m k i ⟨0,hk⟩))^x := by
  have classes : ∀ j : Fin k, Ideal.Quotient.mk (Ideal.span (chainInputs K m k))
      (X (var m k i j))=chainClass m k i j := by
    intro j
    unfold chainClass
    rw [dif_pos j.isLt]
  have first : Ideal.Quotient.mk (Ideal.span (chainInputs K m k)) (X (var m k i ⟨0,hk⟩))=
      chainClass m k i 0 := by
    unfold chainClass
    rw [dif_pos hk]
  simp only [classes,first]
  -- Reduce the exponent modulo `2^k`.
  have digits : ∀ j : Fin k, (if x.testBit j then 1 else 0)=
      (if (x%2^k).testBit j then 1 else 0) := by
    intro j
    rw [Nat.testBit_mod_two_pow]
    simp [j.isLt]
  simp only [digits]
  rw [prod_digits k (chainClass m k i) (fun j hj => chainClass_square m k i j hj) (x%2^k)
    (Nat.mod_lt _ (Nat.two_pow_pos k))]
  conv_rhs => rw [← Nat.mod_add_div x (2^k),pow_add,pow_mul,chainClass_period m k i hk,
    one_pow,mul_one]

theorem monomial_maskExp (m k : ℕ) (e : Fin m →₀ ℕ) (c : K) :
    monomial (maskExp m k e) c=C c*∏ i : Fin m, ∏ j : Fin k,
      X (var m k i j)^(if (e i).testBit j then 1 else 0) := by
  rw [monomial_eq,Finsupp.prod_fintype _ _ (fun _ => pow_zero _)]
  congr 1
  rw [← Fintype.prod_equiv finProdFinEquiv (fun ij => X (finProdFinEquiv ij)^maskExp m k e
    (finProdFinEquiv ij)) _ (fun _ => rfl)]
  rw [Fintype.prod_prod_type]
  simp only [maskExp_var,maskExp_equiv]

theorem monomial_congr (m k : ℕ) (hk : 0<k) (e : Fin m →₀ ℕ) (c : K) :
    Ideal.Quotient.mk (Ideal.span (chainInputs K m k)) (monomial (maskExp m k e) c)=
      Ideal.Quotient.mk (Ideal.span (chainInputs K m k)) (toChain m k hk (monomial e c)) := by
  rw [monomial_maskExp,toChain,rename_monomial,monomial_eq,
    Finsupp.prod_mapDomain_index (fun _ => pow_zero _) (fun _ _ _ => pow_add _ _ _),
    Finsupp.prod_fintype _ _ (fun _ => pow_zero _)]
  simp only [map_mul,map_prod,map_pow]
  congr 1
  apply Finset.prod_congr rfl
  intro i _
  exact block m k hk i (e i)

/-- The lift is congruent to the chain substitution modulo the chains, for every polynomial. -/
theorem lift_congr (m k : ℕ) (hk : 0<k) (p : Z m K) :
    lift m k p-toChain m k hk p∈Ideal.span (chainInputs K m k) := by
  rw [← Ideal.Quotient.eq]
  conv_rhs => rw [p.as_sum]
  unfold lift
  rw [map_sum,map_sum,map_sum]
  apply Finset.sum_congr rfl
  intro e _
  exact monomial_congr m k hk e _

/-- The carried row `NF(u_{i,0} * NF(S))` is the lift of `x_i * S`. -/
theorem carry_congr (m k : ℕ) (hk : 0<k) (p : Z m K) (i : Fin m) :
    lift m k (X i*p)-X (var m k i ⟨0,hk⟩)*lift m k p∈Ideal.span (chainInputs K m k) := by
  have hcong := lift_congr (K := K) m k hk (X i*p)
  have gcong := lift_congr (K := K) m k hk p
  have identity : lift m k (X i*p)-X (var m k i ⟨0,hk⟩)*lift m k p=
      (lift m k (X i*p)-toChain m k hk (X i*p))-
        X (var m k i ⟨0,hk⟩)*(lift m k p-toChain m k hk p) := by
    simp only [map_mul,toChain,rename_X]
    ring
  rw [identity]
  exact Ideal.sub_mem _ hcong (Ideal.mul_mem_left _ _ gcong)

/-! ## The full-degree bound -/

/-- **Full-degree bound.** For every `m, k ≥ 1`, field with `2 ≠ 0` and generator `S`: if
the chain normal form `g = lift S` has the full degree `m k` and, for one block `i`,
`NF(u_{i,0} g)` is not a scalar multiple of `g`, then no `e ≤ m k` is an LFD bound of
`chains ∪ {g}`. -/
theorem lfd_ge_full (two_ne : (2 : K)≠0) (m k : ℕ) (hk : 0<k) (S : Z m K) (i : Fin m)
    (full : (lift m k S).totalDegree=m*k)
    (not_multiple : ∀ c : K, lift m k (X i*S)≠c • lift m k S) :
    ∀ e≤m*k, ¬HasLFDBound (chainInputs K m k∪{lift m k S}) e := by
  apply CarryLower.chain_carry_lower m k two_ne (lift m k S) (lift m k (X i*S))
    (var m k i ⟨0,hk⟩) (m*k) (lift_reduced _ _ _) full (lift_reduced _ _ _) ?_
    (carry_congr m k hk S i) not_multiple
  apply (lift_totalDegree_le m k k _ (fun s _ i => weight_le k (s i))).trans
  exact le_refl _

/-- **Exact value `m k + 1`.** Under the hypotheses of `lfd_ge_full`, the last fall degree of
`chains ∪ {g}` is exactly `m k + 1`: the upper bound is `ChainQuotient.reduced_inputs_lfd_bound`. -/
theorem lfd_eq_full (two_ne : (2 : K)≠0) (m k : ℕ) (hk : 0<k) (S : Z m K) (i : Fin m)
    (full : (lift m k S).totalDegree=m*k)
    (not_multiple : ∀ c : K, lift m k (X i*S)≠c • lift m k S) :
    IsLastFallDegree (chainInputs K m k∪{lift m k S}) (m*k+1) := by
  refine ⟨?_, ?_⟩
  · apply reduced_inputs_lfd_bound m k two_ne
    rintro f hf
    rw [Set.mem_singleton_iff.mp hf]
    exact lift_reduced m k S
  · intro e he
    exact lfd_ge_full two_ne m k hk S i full not_multiple e (by omega)

/-! ## The onset bound -/

/-- The exponent `(2^(m-1)-1, …, 2^(m-1)-1)`. -/
def topExp (m : ℕ) : Fin m →₀ ℕ := Finsupp.equivFunOnFinite.symm fun _ => 2^(m-1)-1

theorem topExp_apply (m : ℕ) (i : Fin m) : topExp m i=2^(m-1)-1 := rfl

theorem two_pow_pred_add_one_lt {m : ℕ} (hm : 3≤m) : 2^(m-1)+1<2^m-1 := by
  have h : 2^m=2*2^(m-1) := by
    rw [← pow_succ']
    congr 1
    omega
  have four : 4≤2^(m-1) := by
    calc 4=2^2 := by norm_num
      _ ≤ 2^(m-1) := Nat.pow_le_pow_right (by norm_num) (by omega)
  omega

theorem two_pow_le {m k : ℕ} (hk : m≤k) : 2^m≤2^k := Nat.pow_le_pow_right (by norm_num) hk

/-- An exponent at most `2^(m-1)+1` has at most `m-1` binary digits (`m ≥ 3`, `k ≥ m`). -/
theorem weight_small {m k x : ℕ} (hm : 3≤m) (hk : m≤k) (hx : x≤2^(m-1)+1) :
    weight k x≤m-1 := by
  have lt := two_pow_pred_add_one_lt hm
  exact weight_le_pred (by omega) (by omega) hk

theorem X_mul_support {m : ℕ} {S : Z m K} (i : Fin m) {b : ℕ}
    (bound : ∀ s∈S.support, ∀ l, s l≤b) {s : Fin m →₀ ℕ} (hs : s∈(X i*S).support) :
    ∀ l, s l≤b+1 := by
  rw [support_X_mul] at hs
  obtain ⟨s',hs',rfl⟩ := Finset.mem_map.mp hs
  intro l
  have := bound s' hs' l
  simp only [addLeftEmbedding_apply,Finsupp.coe_add,Pi.add_apply,Finsupp.single_apply]
  split_ifs <;> omega

/-- **Onset bound.** Let `m ≥ 3`, `k ≥ m`, `2 ≠ 0` in `K`, and let `S` have degree at most
`2^(m-1)` in each of the `m` summand coordinates, with a nonzero coefficient at
`x_1^(2^(m-1)-1) ⋯ x_m^(2^(m-1)-1)`. Then the direct chain system `chains ∪ {NF(S)}` has
`LFD ≥ m(m-1)+1`: no `e ≤ m(m-1)` is an LFD bound. -/
theorem lfd_ge_onset (two_ne : (2 : K)≠0) (m k : ℕ) (hm : 3≤m) (hk : m≤k) (S : Z m K)
    (bound : ∀ s∈S.support, ∀ i, s i≤2^(m-1))
    (top : coeff (topExp m) S≠0) :
    ∀ e≤m*(m-1), ¬HasLFDBound (chainInputs K m k∪{lift m k S}) e := by
  have hk0 : 0<k := by omega
  have i0 : Fin m := ⟨0,by omega⟩
  have below : 2^(m-1)+1<2^k :=
    lt_of_lt_of_le (lt_of_lt_of_le (two_pow_pred_add_one_lt hm) (Nat.sub_le _ _)) (two_pow_le hk)
  have small : ∀ s∈S.support, ∀ i, s i<2^k := fun s hs i => by have := bound s hs i; omega
  have small' : ∀ s∈(X i0*S).support, ∀ i, s i<2^k :=
    fun s hs i => by have := X_mul_support i0 bound hs i; omega
  -- The generator has degree exactly `m(m-1)`.
  have g_degree : (lift m k S).totalDegree=m*(m-1) := by
    apply le_antisymm
    · exact lift_totalDegree_le m k (m-1) S
        (fun s hs i => weight_small hm hk ((bound s hs i).trans (Nat.le_succ _)))
    · have topSmall : ∀ i, topExp m i<2^k := by
        intro i
        rw [topExp_apply]
        exact lt_of_le_of_lt (Nat.sub_le _ _) (lt_trans (Nat.lt_succ_self _) below)
      have coefficient : coeff (maskExp m k (topExp m)) (lift m k S)≠0 := by
        rw [coeff_lift small topSmall]
        exact top
      have := le_totalDegree (mem_support_iff.mpr coefficient)
      rw [maskExp_sum] at this
      simp only [topExp_apply,weight_two_pow_sub_one (show m-1≤k by omega),
        Finset.sum_const,Finset.card_univ,Fintype.card_fin,smul_eq_mul] at this
      exact this
  -- The carried row has degree at most `m(m-1)`.
  have h_degree : (lift m k (X i0*S)).totalDegree≤m*(m-1) :=
    lift_totalDegree_le m k (m-1) _
      (fun s hs i => weight_small hm hk (X_mul_support i0 bound hs i))
  -- The carried row is not a multiple of the generator.
  have not_multiple : ∀ c : K, lift m k (X i0*S)≠c • lift m k S := by
    intro c equal
    have polynomial : X i0*S=C c*S := by
      ext e
      by_cases he : ∀ i, e i<2^k
      · have := congrArg (coeff (maskExp m k e)) equal
        rw [coeff_smul,coeff_lift small' he,coeff_lift small he] at this
        rw [this,coeff_C_mul,smul_eq_mul]
      · push_neg at he
        obtain ⟨i,hi⟩ := he
        have z1 : coeff e (X i0*S)=0 := by
          rw [← not_mem_support_iff]
          intro mem
          exact absurd (small' e mem i) (not_lt.mpr hi)
        have z2 : coeff e S=0 := by
          rw [← not_mem_support_iff]
          intro mem
          exact absurd (small e mem i) (not_lt.mpr hi)
        rw [z1,coeff_C_mul,z2,mul_zero]
    have factor : (X i0-C c)*S=0 := by rw [sub_mul,polynomial,sub_self]
    have linear : (X i0-C c : Z m K)≠0 := by
      intro zero
      have := congrArg (coeff (Finsupp.single i0 1)) zero
      rw [coeff_sub,coeff_X,coeff_C,if_neg] at this
      · simp at this
      · intro h
        have := congrArg (fun s : Fin m →₀ ℕ => s i0) h
        simp at this
    have nonzero : S≠0 := by
      intro zero
      rw [zero,coeff_zero] at top
      exact top rfl
    exact nonzero ((mul_eq_zero.mp factor).resolve_left linear)
  exact CarryLower.chain_carry_lower m k two_ne (lift m k S) (lift m k (X i0*S))
    (var m k i0 ⟨0,hk0⟩) (m*(m-1)) (lift_reduced _ _ _) g_degree (lift_reduced _ _ _) h_degree
    (carry_congr m k hk0 S i0) not_multiple

/-! ## Three summands as an instance -/

theorem topExp_three : topExp 3=M3CarryLower.expo 3 3 3 := by
  ext i
  fin_cases i <;> simp [topExp_apply,M3CarryLower.expo,Finsupp.single_apply]

/-- The three-summand instance of the onset bound: for `k ≥ 3` and `5t^3+3at+6b ≠ 0`, the
direct system with `g = NF(S_4(x1,x2,x3,t))` has `LFD ≥ 7`. -/
theorem m3_onset (two_ne : (2 : K)≠0) (k : ℕ) (hk : 3≤k) (a b t : K)
    (nondegenerate : 5*t^3+3*a*t+6*b≠0) :
    ∀ e≤6, ¬HasLFDBound (chainInputs K 3 k∪{lift 3 k (M3CarryLower.s4 a b t)}) e := by
  apply lfd_ge_onset two_ne 3 k le_rfl hk
  · intro s hs i
    exact (M3CarryLower.s4_support hs i).trans (by norm_num)
  · rw [topExp_three,M3CarryLower.coeff_s4_333]
    have eight : (-8 : K)≠0 := by
      have : (-8 : K)=-(2^3) := by norm_num
      rw [this]
      exact neg_ne_zero.mpr (pow_ne_zero _ two_ne)
    exact mul_ne_zero eight nondegenerate

end ICFallDegrees.CarryAllM
