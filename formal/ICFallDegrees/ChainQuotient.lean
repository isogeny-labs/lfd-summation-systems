import ICFallDegrees.GridQuotient
import ICFallDegrees.PrimeFieldConclusions
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.FieldTheory.Separable
import Mathlib.Algebra.GeomSum
import Mathlib.Algebra.MvPolynomial.Monad

noncomputable section

/-!
# Subgroup-chain presentations: a Groebner subsystem and the `m*k+1` bound

The two-power subgroup presentation of a prime-field factor base writes each
summand coordinate `x_i` through a chain `u_{i,0}, …, u_{i,k-1}` with
`u_{i,j}^2 = u_{i,j+1}` and `u_{i,k-1}^2 = 1`, so that `x_i = u_{i,0}` ranges
over the subgroup of order `2^k`.  Here the chain equations for `m` summands
of depth `k` are shown to be a reduced-ideal-trivial subsystem of every input
list containing them, in the sense of `FiniteQuotient`: a polynomial of the
chain ideal whose monomials are all squarefree is zero.  The proof expands
each chain variable as a two-power of a single summand variable,
`u_{i,j} ↦ z_i^(2^j)`; squarefree chain monomials become the distinct
monomials of partial degree below `2^k`, and the chain ideal maps into the
ideal of `z_i^(2^k) - 1`, which vanishes on the grid of `2^k`-th roots of
unity in an algebraic closure.  Mathlib's grid vanishing theorem then forces
the image to vanish.  No Groebner basis, normal form, closure rank, or
solution count is assumed.

Consequences: for every `m`, `k`, every field of characteristic other than
two, and every input list containing the chain equations, the last fall
degree is at most `max (m*k+1) delta`, where `delta` bounds the input degrees.
For the normalized two-summand presentation of `PrimeM2General` this gives
`LFD ≤ 2*k+1` for every depth `k ≥ 2`, prime `p ≥ 5`, curve, and target.
The bound concerns the specified chain encoding; it says nothing about
relation yield or cryptographic cost.
-/

namespace ICFallDegrees.ChainQuotient

open MvPolynomial MonomialOrder FiniteQuotient
open scoped MonomialOrder

section Variables

/-- The summand (row) index of a chain variable. -/
def row {m k : ℕ} (v : Fin (m*k)) : Fin m := (finProdFinEquiv.symm v).1

/-- The squaring depth of a chain variable. -/
def depth {m k : ℕ} (v : Fin (m*k)) : Fin k := (finProdFinEquiv.symm v).2

/-- The chain variable of summand `i` at depth `j`. -/
abbrev var (m k : ℕ) (i : Fin m) (j : Fin k) : Fin (m*k) := finProdFinEquiv (i,j)

@[simp] theorem row_var (m k : ℕ) (i : Fin m) (j : Fin k) : row (var m k i j) = i := by
  simp [row, var]

@[simp] theorem depth_var (m k : ℕ) (i : Fin m) (j : Fin k) : depth (var m k i j) = j := by
  simp [depth, var]

theorem var_row_depth {m k : ℕ} (v : Fin (m*k)) : var m k (row v) (depth v) = v :=
  finProdFinEquiv.apply_symm_apply v

theorem row_val {m k : ℕ} (v : Fin (m*k)) : (row v : ℕ) = v / k := Fin.coe_divNat v

theorem depth_val {m k : ℕ} (v : Fin (m*k)) : (depth v : ℕ) = v % k := Fin.coe_modNat v

theorem var_val (m k : ℕ) (i : Fin m) (j : Fin k) : (var m k i j : ℕ) = j + k * i := rfl

end Variables

section Chain

variable (K : Type*) [Field K]

/-- The square of a chain variable equals the next chain variable, or one at
the end of a chain. -/
def chainTail (m k : ℕ) (v : Fin (m*k)) : P (m*k) K :=
  if h : (depth v).val+1 < k then X (var m k (row v) ⟨(depth v).val+1,h⟩) else 1

def chainPolynomial (m k : ℕ) (v : Fin (m*k)) : P (m*k) K :=
  X v^2-chainTail K m k v

/-- All `m*k` chain equations of `m` summands at depth `k`. -/
def chainInputs (m k : ℕ) : Set (P (m*k) K) := Set.range (chainPolynomial K m k)

/-- Expansion of each chain variable as a two-power of its summand variable. -/
def expand (m k : ℕ) : P (m*k) K →ₐ[K] MvPolynomial (Fin m) K :=
  bind₁ (fun v => X (row v)^2^(depth v : ℕ))

end Chain

variable {K : Type*} [Field K]

theorem chainTail_totalDegree_le (m k : ℕ) (v : Fin (m*k)) :
    (chainTail K m k v).totalDegree≤1 := by
  unfold chainTail
  split_ifs <;> simp

theorem lt_single_two_of_degree_le_one {r : ℕ} (c : Fin r →₀ ℕ) (hc : c.degree≤1) (v : Fin r) :
    c≺[degLex] Finsupp.single v 2 := by
  rw [degLex_lt_iff, Finsupp.DegLex.lt_iff]
  left
  simp only [ofDegLex_toDegLex, Finsupp.degree_single]
  omega

theorem leading_chainPolynomial (m k : ℕ) (v : Fin (m*k)) :
    degLex.degree (chainPolynomial K m k v)=Finsupp.single v 2 := by
  classical
  unfold chainPolynomial
  have square : degLex.degree (X v^2 : P (m*k) K)=Finsupp.single v 2 := by
    rw [X_pow_eq_monomial, degree_monomial]
    simp
  rw [degree_sub_of_lt, square]
  rw [square, degree_lt_iff]
  · intro c hc
    apply lt_single_two_of_degree_le_one
    have bound := (le_totalDegree hc).trans (chainTail_totalDegree_le (K := K) m k v)
    simpa [Finsupp.degree, Finsupp.sum] using bound
  · exact lt_single_two_of_degree_le_one 0 (by simp) v

theorem chainPolynomial_nonzero (m k : ℕ) (v : Fin (m*k)) : chainPolynomial K m k v≠0 := by
  intro zero
  have lead := leading_chainPolynomial (K := K) m k v
  rw [zero, degree_zero] at lead
  exact (Finsupp.single_ne_zero.mpr (by norm_num : (2 : ℕ)≠0)) lead.symm

theorem chainInputs_nonzero (m k : ℕ) : (0 : P (m*k) K)∉chainInputs K m k := by
  rintro ⟨v,hv⟩
  exact chainPolynomial_nonzero m k v hv

/-- Standard exponents of the chain subsystem are squarefree. -/
theorem standard_le_one (m k : ℕ) (e : Fin (m*k) →₀ ℕ)
    (standard : Standard (chainInputs K m k) e) (v : Fin (m*k)) : e v≤1 := by
  have h := standard (chainPolynomial K m k v) ⟨v,rfl⟩
  rw [leading_chainPolynomial, Finsupp.single_le_iff] at h
  omega

theorem standard_degree_bound (m k : ℕ) (e : Fin (m*k) →₀ ℕ)
    (standard : Standard (chainInputs K m k) e) : e.degree≤m*k := by
  have sum : e.degree=∑ i : Fin (m*k), e i :=
    Finsupp.sum_fintype e (fun _ n => n) (fun _ => rfl)
  rw [sum]
  calc ∑ i : Fin (m*k), e i≤∑ _i : Fin (m*k), 1 :=
        Finset.sum_le_sum (fun i _ => standard_le_one m k e standard i)
    _ = m*k := by simp

theorem expand_X (m k : ℕ) (v : Fin (m*k)) :
    expand K m k (X v)=X (row v)^2^(depth v : ℕ) := bind₁_X_right _ _

/-- The chain equations expand to zero or to `z_i^(2^k) - 1`. -/
theorem expand_chainPolynomial (m k : ℕ) (v : Fin (m*k)) :
    expand K m k (chainPolynomial K m k v)=
      if (depth v).val+1<k then 0 else X (row v)^2^k-1 := by
  unfold chainPolynomial chainTail
  split_ifs with h
  · simp only [map_sub, map_pow, expand_X, row_var, depth_var, sub_eq_zero]
    rw [← pow_mul, ← pow_succ]
  · have top : (depth v).val+1=k := by have := (depth v).isLt; omega
    simp only [map_sub, map_pow, expand_X, map_one]
    rw [← pow_mul, ← pow_succ, top]

/-- The exponent of the expanded squarefree monomial: binary digits per summand. -/
def bitExp {m k : ℕ} (e : Fin (m*k) →₀ ℕ) : Fin m →₀ ℕ :=
  ∑ v ∈ e.support, Finsupp.single (row v) (2^(depth v : ℕ)*e v)

theorem expand_monomial (m k : ℕ) (e : Fin (m*k) →₀ ℕ) (c : K) :
    expand K m k (monomial e c)=monomial (bitExp e) c := by
  classical
  unfold expand bitExp
  rw [bind₁_monomial]
  have product : ∏ v ∈ e.support, (X (row v)^2^(depth v : ℕ))^e v=
      monomial (∑ v ∈ e.support, Finsupp.single (row v) (2^(depth v : ℕ)*e v)) (1 : K) := by
    rw [monomial_sum_one]
    apply Finset.prod_congr rfl
    intro v _
    rw [← pow_mul, X_pow_eq_monomial]
  rw [product, C_mul_monomial, mul_one]

theorem bitExp_apply {m k : ℕ} (e : Fin (m*k) →₀ ℕ) (i : Fin m) :
    bitExp e i=∑ j : Fin k, 2^(j : ℕ)*e (var m k i j) := by
  classical
  unfold bitExp
  rw [Finsupp.finset_sum_apply]
  simp only [Finsupp.single_apply]
  rw [Finset.sum_subset (Finset.subset_univ e.support)]
  · rw [← Equiv.sum_comp finProdFinEquiv, Fintype.sum_prod_type]
    simp only [row_var, depth_var]
    rw [Finset.sum_comm]
    simp [Finset.sum_ite_eq']
  · intro v _ hv
    simp [Finsupp.not_mem_support_iff.mp hv]

theorem geom_lt (k : ℕ) : ∑ j : Fin k, 2^(j : ℕ)<2^k := by
  rw [Fin.sum_univ_eq_sum_range (fun j => 2^j) k]
  exact Nat.geomSum_lt (le_refl 2) (fun j hj => Finset.mem_range.mp hj)

theorem bitExp_lt {m k : ℕ} (e : Fin (m*k) →₀ ℕ) (he : ∀ v, e v≤1) (i : Fin m) :
    bitExp e i<2^k := by
  rw [bitExp_apply]
  calc ∑ j : Fin k, 2^(j : ℕ)*e (var m k i j)≤∑ j : Fin k, 2^(j : ℕ) := by
        apply Finset.sum_le_sum
        intro j _
        calc 2^(j : ℕ)*e (var m k i j)≤2^(j : ℕ)*1 := Nat.mul_le_mul_left _ (he _)
          _ = 2^(j : ℕ) := mul_one _
    _ < 2^k := geom_lt k

/-- Binary digits are determined by their weighted sum. -/
theorem binary_sum_injective : ∀ (k : ℕ) (b b' : Fin k → ℕ), (∀ j, b j≤1) → (∀ j, b' j≤1) →
    ∑ j : Fin k, 2^(j : ℕ)*b j=∑ j : Fin k, 2^(j : ℕ)*b' j → ∀ j, b j=b' j := by
  intro k
  induction k with
  | zero => intro b b' _ _ _ j; exact j.elim0
  | succ k ih =>
    intro b b' hb hb' hsum
    rw [Fin.sum_univ_succ, Fin.sum_univ_succ] at hsum
    have shift : ∀ c : Fin (k+1) → ℕ,
        ∑ j : Fin k, 2^((j.succ : Fin (k+1)) : ℕ)*c j.succ=2*∑ j : Fin k, 2^(j : ℕ)*c j.succ := by
      intro c
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro j _
      rw [Fin.val_succ, pow_succ]
      ring
    rw [shift b, shift b'] at hsum
    simp only [Fin.val_zero, pow_zero, one_mul] at hsum
    have h0 := hb 0
    have h0' := hb' 0
    have head : b 0=b' 0 := by omega
    have tails : ∑ j : Fin k, 2^(j : ℕ)*b j.succ=∑ j : Fin k, 2^(j : ℕ)*b' j.succ := by omega
    have rest := ih (fun j => b j.succ) (fun j => b' j.succ) (fun j => hb _) (fun j => hb' _) tails
    intro j
    exact Fin.cases head (fun j => rest j) j

theorem bitExp_injective {m k : ℕ} (e e' : Fin (m*k) →₀ ℕ)
    (he : ∀ v, e v≤1) (he' : ∀ v, e' v≤1) (h : bitExp e=bitExp e') : e=e' := by
  ext v
  have hi := congrArg (fun f : Fin m →₀ ℕ => f (row v)) h
  simp only [bitExp_apply] at hi
  have digits := binary_sum_injective k (fun j => e (var m k (row v) j))
    (fun j => e' (var m k (row v) j)) (fun j => he _) (fun j => he' _) hi (depth v)
  simpa [var_row_depth] using digits

/-- The chain subsystem is reduced-ideal trivial: a polynomial of the chain
ideal with only squarefree monomials is zero. -/
theorem reduced_chain_ideal_trivial (m k : ℕ) (two_ne : (2 : K)≠0) :
    ReducedIdealTrivial (chainInputs K m k) := by
  classical
  intro p hp reduced
  have multilinear : ∀ e ∈ p.support, ∀ v, e v≤1 :=
    fun e he v => standard_le_one m k e (reduced e he) v
  let L := AlgebraicClosure K
  let Φ : P (m*k) K →+* MvPolynomial (Fin m) L :=
    (MvPolynomial.map (algebraMap K L)).comp (expand K m k).toRingHom
  have expansion : Φ p=∑ e ∈ p.support, monomial (bitExp e) (algebraMap K L (coeff e p)) := by
    conv_lhs => rw [p.as_sum]
    rw [map_sum]
    apply Finset.sum_congr rfl
    intro e _
    simp only [Φ, RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe,
      expand_monomial, map_monomial]
  have hN : ((2^k : ℕ) : K)≠0 := by
    rw [Nat.cast_pow]
    exact pow_ne_zero _ two_ne
  let q : Polynomial K := Polynomial.X^2^k-Polynomial.C 1
  have sep : q.Separable := Polynomial.separable_X_pow_sub_C 1 hN one_ne_zero
  have spl : Polynomial.Splits (algebraMap K L) q := IsAlgClosed.splits_codomain q
  let S : Finset L := (q.rootSet L).toFinset
  have cardS : S.card=2^k := by
    rw [Set.toFinset_card, Polynomial.card_rootSet_eq_natDegree sep spl]
    exact Polynomial.natDegree_X_pow_sub_C
  have rootS : ∀ x ∈ S, x^2^k=1 := by
    intro x hx
    rw [Set.mem_toFinset, Polynomial.mem_rootSet] at hx
    have := hx.2
    simpa [q, sub_eq_zero] using this
  have vanish : ∀ x : Fin m → L, (∀ i, x i ∈ S) → eval x (Φ p)=0 := by
    intro x hx
    have le : Ideal.span (chainInputs K m k)≤RingHom.ker ((eval x).comp Φ) := by
      apply Ideal.span_le.mpr
      rintro g ⟨v,rfl⟩
      rw [SetLike.mem_coe, RingHom.mem_ker, RingHom.comp_apply]
      simp only [Φ, RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe,
        expand_chainPolynomial]
      split_ifs
      · simp
      · simp [eval_map, rootS _ (hx _)]
    exact le hp
  have degree : ∀ i, (Φ p).degreeOf i<S.card := by
    intro i
    rw [cardS, degreeOf_lt_iff (by positivity)]
    intro d hd
    rw [expansion] at hd
    obtain ⟨e,he,hd⟩ := Finset.mem_biUnion.mp (support_sum hd)
    have single := support_monomial_subset hd
    rw [Finset.mem_singleton] at single
    rw [single]
    exact bitExp_lt e (multilinear e he) i
  have zero : Φ p=0 := eq_zero_of_eval_zero_at_prod_finset (Φ p) (fun _ => S) degree vanish
  by_contra nonzero
  have nonempty : p.support.Nonempty :=
    Finset.nonempty_iff_ne_empty.mpr (fun h => nonzero (support_eq_empty.mp h))
  obtain ⟨e₀,he₀⟩ := nonempty
  have coefficient := congrArg (coeff (bitExp e₀)) zero
  rw [expansion, coeff_sum, coeff_zero] at coefficient
  simp only [coeff_monomial] at coefficient
  rw [Finset.sum_eq_single e₀] at coefficient
  · rw [if_pos rfl] at coefficient
    exact (mem_support_iff.mp he₀) ((map_eq_zero_iff _ (algebraMap K L).injective).mp coefficient)
  · intro e he hne
    rw [if_neg]
    intro heq
    exact hne (bitExp_injective e e₀ (multilinear e he) (multilinear e₀ he₀) heq)
  · intro h
    exact absurd he₀ h

/-- The chain bound: every input list containing the `m*k` chain equations
has last fall degree at most `max (m*k+1) delta`, where `delta` bounds all
input degrees. Nothing about the other inputs is assumed. -/
theorem chain_lfd_bound (m k : ℕ) (two_ne : (2 : K)≠0) (F : Set (P (m*k) K))
    (included : chainInputs K m k⊆F) (delta : ℕ)
    (input_degree : ∀ f ∈ F, f.totalDegree≤delta) :
    HasLFDBound F (max (m*k+1) delta) :=
  finite_quotient_lfd_of_reduced_zero F (chainInputs K m k) included (chainInputs_nonzero m k)
    (reduced_chain_ideal_trivial m k two_ne) (m*k) delta (standard_degree_bound m k) input_degree

/-- With normalized inputs of degree at most `m*k+1`, the bound is `m*k+1`. -/
theorem normalized_chain_lfd_bound (m k : ℕ) (two_ne : (2 : K)≠0) (F : Set (P (m*k) K))
    (included : chainInputs K m k⊆F)
    (input_degree : ∀ f ∈ F, f.totalDegree≤m*k+1) : HasLFDBound F (m*k+1) := by
  simpa using chain_lfd_bound m k two_ne F included (m*k+1) input_degree

/-- Every polynomial of the reduced space of the chains (squarefree in the `m*k` chain
variables) has total degree at most `m*k`. -/
theorem reduced_totalDegree_le (m k : ℕ) (f : P (m*k) K)
    (hf : f ∈ reducedSpace (chainInputs K m k)) : f.totalDegree≤m*k := by
  apply Finset.sup_le
  intro e he
  exact standard_degree_bound m k e (hf e he)

theorem chainPolynomial_totalDegree_le (m k : ℕ) (v : Fin (m*k)) :
    (chainPolynomial K m k v).totalDegree≤2 := by
  unfold chainPolynomial
  apply (totalDegree_sub_le _ _).trans
  apply max_le
  · simp
  · exact (chainTail_totalDegree_le (K := K) m k v).trans (by norm_num)

/-- The unconditional dimension bound for every number of summands: the chain equations together
with any inputs in the reduced space of the chains (for example the normal form of a summation
polynomial in any number of summands) have `LFD ≤ m*k+1`. -/
theorem reduced_inputs_lfd_bound (m k : ℕ) (two_ne : (2 : K)≠0) (R : Set (P (m*k) K))
    (reduced : R⊆reducedSpace (chainInputs K m k)) :
    HasLFDBound (chainInputs K m k∪R) (m*k+1) := by
  apply normalized_chain_lfd_bound m k two_ne _ Set.subset_union_left
  rintro f (⟨v,rfl⟩ | hf)
  · have hv := v.isLt
    exact (chainPolynomial_totalDegree_le (K := K) m k v).trans (by omega)
  · exact (reduced_totalDegree_le m k f (reduced hf)).trans (Nat.le_succ _)

/-! ## The two-summand prime presentation at every depth -/

theorem chainInputs_two_subset (k : ℕ) (depth : 2≤k) (a b t : K) :
    chainInputs K 2 k⊆Set.range (PrimeM2General.inputs k depth a b t) := by
  rintro g ⟨v,rfl⟩
  refine ⟨Fin.castSucc v,?_⟩
  simp only [PrimeM2General.inputs, Fin.lastCases_castSucc]
  unfold chainPolynomial chainTail PrimeM2General.chainTail
  have hv := v.isLt
  have decomposition := Nat.mod_add_div (v : ℕ) k
  have remainder := Nat.mod_lt (v : ℕ) (by omega : 0<k)
  have quotient : (v : ℕ)/k<2 := (Nat.div_lt_iff_lt_mul (by omega)).mpr hv
  rw [depth_val]
  have cases : (v : ℕ)/k=0 ∨ (v : ℕ)/k=1 := by
    rcases (Nat.lt_succ_iff.mp quotient).lt_or_eq with q | q
    · exact Or.inl (Nat.lt_one_iff.mp q)
    · exact Or.inr q
  split_ifs with h₁ h₂ h₂
  · exfalso
    rcases cases with q | q <;> rw [q] at decomposition <;>
      simp only [Nat.mul_zero, Nat.mul_one, Nat.add_zero] at decomposition <;> omega
  · rw [sub_eq_add_neg]
  · rw [sub_eq_add_neg]
    congr 2
    congr 1
    apply Fin.ext
    rw [var_val, row_val]
    simp only [Fin.val_mk]
    linarith [Nat.mod_add_div (v : ℕ) k]
  · exfalso
    rcases cases with q | q <;> rw [q] at decomposition <;>
      simp only [Nat.mul_zero, Nat.mul_one, Nat.add_zero] at decomposition <;> omega

/-- For every depth `k ≥ 2`, the normalized two-summand chain presentation
has last fall degree at most `2*k+1`, over every field of characteristic
other than two and for every curve and target parameter. -/
theorem prime_m2_lfd_bound (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a b t : K) :
    HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (2*k+1) := by
  have degrees : ∀ f ∈ Set.range (PrimeM2General.inputs k depth a b t), f.totalDegree≤2 := by
    rintro f ⟨i,rfl⟩
    exact (PrimeM2General.degree_and_top k depth a b t i).1.le
  have bound := chain_lfd_bound 2 k two_ne _ (chainInputs_two_subset k depth a b t) 2 degrees
  have collapse : max (2*k+1) 2=2*k+1 := by omega
  rwa [collapse] at bound

/-- Prime-field entry point: `p ≥ 5`, every depth, curve and target. -/
theorem prime_field_m2_lfd_bound (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (a b t : ZMod p) :
    HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (2*k+1) :=
  prime_m2_lfd_bound (PrimeFieldConclusions.two_ne_zero p large) k depth a b t

end ICFallDegrees.ChainQuotient
