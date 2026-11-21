import ICFallDegrees.AffineTransport
import ICFallDegrees.GrowthRationalUnit
import ICFallDegrees.GrowthRationalPrimeNative
import ICFallDegrees.ChainClosureBridge
import ICFallDegrees.BinaryWeight

/-!
Degree-controlled specialization for growth lower bounds.

These transport lemmas do not assert that a Semaev specialization satisfies
their input hypotheses, or that its source ideal is the unit ideal.
The preprocessing budget is explicit and may exceed the source budget.
-/
noncomputable section
namespace ICFallDegrees.SummandGrowthTransfer
open MvPolynomial
variable {σ τ K : Type*} [CommRing K]

theorem map_closure_budget
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D E : ℕ) (hDE : D ≤ E)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G E (φ f))
    {f} (hf : InClosure F D f) : InClosure G E (φ f) := by
  induction hf with
  | input member degree => exact inputs _ member degree
  | zero => simpa using (InClosure.zero : InClosure G E 0)
  | add _ _ ih ih' => simpa only [map_add] using ih.add ih'
  | smul c _ ih => simpa only [map_smul] using ih.smul c
  | mul a _ budget ih =>
      simpa only [map_mul] using ih.mul (φ a)
        (((Nat.add_le_add (AffineTransport.degree_le φ linear a)
          (AffineTransport.degree_le φ linear _)).trans budget).trans hDE)

/-- A fixed preprocessing cost does not become a cost at every derivation step. -/
theorem not_one_of_specialization
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D C : ℕ)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G (max D C) (φ f))
    (excluded : ¬ InClosure G (max D C) 1) : ¬ InClosure F D 1 := by
  intro h
  apply excluded
  simpa using map_closure_budget φ linear F G D (max D C)
    (Nat.le_max_left D C) inputs h

/-- Source unit-ideal status is a separate hypothesis, not implied by specialization. -/
theorem no_lfd_bound_of_specialization
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D C : ℕ)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G (max D C) (φ f))
    (excluded : ¬ InClosure G (max D C) 1)
    (unit : (1 : MvPolynomial σ K) ∈ Ideal.span F) : ¬ HasLFDBound F D := by
  intro bound
  exact not_one_of_specialization φ linear F G D C inputs excluded
    (one_of_hasLFDBound unit bound)

/-- Explicit prime-field transfer at arbitrary depth, conditional on the actual
specialized inputs. This makes no higher-summand Semaev identity assumption implicit. -/
theorem prime_growth_transfer {σ : Type*}
    (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k) (a : ZMod p)
    (F : Set (MvPolynomial σ (ZMod p)))
    (φ : MvPolynomial σ (ZMod p) →ₐ[ZMod p]
      MvPolynomial (Fin (2*k)) (ZMod p))
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (D C : ℕ) (positive : 1 ≤ max D C)
    (hk : 14 * max D C + 10 ≤ k) (hp : 2 * 2^k + 1 ≤ p)
    (ha : a + 3 ≠ 0)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D →
      InClosure (Set.range (PrimeM2General.inputs k depth a (-1-a) 1))
        (max D C) (φ f)) : ¬ InClosure F D 1 := by
  apply not_one_of_specialization φ linear F
    (Set.range (PrimeM2General.inputs k depth a (-1-a) 1)) D C inputs
  exact GrowthRationalPrimeNative.not_one_of_large_prime p k a depth
    (max D C) positive hk hp ha

/-- The determinant formula for the formal-degree-(2,2) Sylvester resultant. -/
def quadraticSylvester (a b c d e f : K) : K :=
  (a*f-c*d)^2-(a*e-b*d)*(b*f-c*e)

theorem quadraticSylvester_constant (a b c f : K) :
    quadraticSylvester a b c 0 0 f = a^2*f^2 := by
  unfold quadraticSylvester
  ring

/-- Specialize the universal quadratic Sylvester formula, retaining degree two.
This is a formula identity; it does not redefine resultant at the reduced degree. -/
theorem order_two_pair_specialization (a x y : K) :
    quadraticSylvester ((x-y)^2)
      (-2*((x+y)*(x*y+a)+2*(-1-a)))
      ((x*y-a)^2-4*(-1-a)*(x+y))
      ((1-1)^2) (-2*((1+1)*(1*1+a)+2*(-1-a)))
      ((1*1-a)^2-4*(-1-a)*(1+1)) =
      (a+3)^4*(x-y)^4 := by
  unfold quadraticSylvester
  ring

/-- Scalar recurrence produced by repeatedly removing two order-two arguments.
The connection to universal Semaev polynomials is a separate written argument. -/
def pairIterate (c s : K) : ℕ → K
  | 0 => s
  | r+1 => c^(4^(r+1)) * pairIterate c s r ^ 4

theorem pairIterate_closed (c s : K) (r : ℕ) :
    pairIterate c s r = c^(r*4^r) * s^(4^r) := by
  induction r with
  | zero => simp [pairIterate]
  | succ r ih =>
      rw [pairIterate, ih, mul_pow, ← pow_mul, ← pow_mul, ← mul_assoc, ← pow_add]
      have h1 : 4^(r+1) + r*4^r*4 = (r+1)*4^(r+1) := by
        rw [pow_succ]
        ring
      have h2 : 4^r*4 = 4^(r+1) := by rw [pow_succ]
      rw [h1, h2]

section Counting
variable {A : Type*} [AddCommMonoid A] [DecidableEq A]

/-- The points in B need not be rational over the field of the target set. -/
def sumTargets (B : Finset A) (m : ℕ) : Finset A := by
  classical
  exact Finset.univ.image (fun xs : Fin m → {x // x ∈ B} => ∑ i, (xs i : A))

theorem card_sumTargets_le (B : Finset A) (m : ℕ) :
    (sumTargets B m).card ≤ B.card ^ m := by
  classical
  unfold sumTargets
  calc
    _ ≤ (Finset.univ : Finset (Fin m → {x // x ∈ B})).card := Finset.card_image_le
    _ = B.card ^ m := by simp

theorem exists_target_outside_sums (B Q Z : Finset A) (m : ℕ)
    (enough : B.card ^ m + Z.card < Q.card) :
    ∃ R ∈ Q, R ∉ Z ∧ ∀ xs : Fin m → {x // x ∈ B}, (∑ i, (xs i : A)) ≠ R := by
  classical
  have not_subset : ¬ Q ⊆ sumTargets B m ∪ Z := by
    intro h
    have h1 := Finset.card_le_card h
    have h2 := Finset.card_union_le (sumTargets B m) Z
    have h3 := card_sumTargets_le B m
    omega
  obtain ⟨R, hR, outside⟩ := Finset.not_subset.mp not_subset
  have hsum : R ∉ sumTargets B m := fun h => outside (Finset.mem_union_left Z h)
  have hZ : R ∉ Z := fun h => outside (Finset.mem_union_right _ h)
  refine ⟨R, hR, hZ, ?_⟩
  intro xs eq
  apply hsum
  unfold sumTargets
  exact Finset.mem_image.mpr ⟨xs, Finset.mem_univ _, eq⟩
end Counting

/-- Repeated squaring pays the maximum degree budget, not the sum of all
budgets. Concrete chain-normal-form corrections must still be supplied. -/
theorem normalized_square_tower (F : Set (MvPolynomial σ K))
    (h : ℕ → MvPolynomial σ K)
    (degree : ∀ j, (h j).totalDegree ≤ 2*(j+1))
    (base : InClosure F 2 (h 0))
    (correction : ∀ j, InClosure F (4*(j+1)) (h j ^ 2 - h (j+1)))
    (n : ℕ) : InClosure F (max 2 (4*n)) (h n) := by
  induction n with
  | zero => simpa using base
  | succ n ih =>
      have prior : InClosure F (4*(n+1)) (h n) := ih.mono (by omega)
      have square : InClosure F (4*(n+1)) (h n ^ 2) := by
        simpa only [pow_two] using prior.mul (h n) (by have := degree n; omega)
      have next := square.sub (correction n)
      have identity : h n ^ 2 - (h n ^ 2 - h (n+1)) = h (n+1) := by ring
      rw [identity] at next
      simpa only [max_eq_right (by omega : 2 ≤ 4*(n+1))] using next

section ActualNormalForms
open FiniteQuotient ChainClosureBridge
variable {r : ℕ} {L : Type*} [Field L]

theorem nf_square (G : Set (P r L)) (nz : (0 : P r L) ∉ G)
    (z : ReducedIdealTrivial G) (p : P r L) :
    normalForm G nz ((normalForm G nz p)^2) = normalForm G nz (p^2) := by
  have h1 := normalForm_mul_normalForm G nz z p p
  have h2 := normalForm_mul_normalForm G nz z (normalForm G nz p) p
  calc
    normalForm G nz ((normalForm G nz p)^2) =
        normalForm G nz (normalForm G nz p * p) := by
      simpa only [pow_two] using h2.symm
    _ = normalForm G nz (p * normalForm G nz p) := by rw [mul_comm]
    _ = normalForm G nz (p^2) := by simpa only [pow_two] using h1.symm

/-- The normalization corrections are constructed by graded division, not
supplied as closure assumptions. Only the support-degree estimate remains. -/
theorem actual_normal_power_tower (G F : Set (P r L))
    (nz : (0 : P r L) ∉ G) (z : ReducedIdealTrivial G) (included : G ⊆ F)
    (g : P r L) (base : InClosure F 2 (normalForm G nz g))
    (degree : ∀ j, (normalForm G nz (g^(2^j))).totalDegree ≤ 2*(j+1))
    (n : ℕ) : InClosure F (max 2 (4*n)) (normalForm G nz (g^(2^n))) := by
  induction n with
  | zero => simpa using base
  | succ n ih =>
    let h := normalForm G nz (g^(2^n))
    have prior : InClosure F (4*(n+1)) h := ih.mono (by omega)
    have hd : h.totalDegree ≤ 2*(n+1) := degree n
    have budget : (h^2).totalDegree ≤ 4*(n+1) := by
      have := totalDegree_mul h h
      rw [← pow_two] at this
      omega
    have square : InClosure F (4*(n+1)) (h^2) := by
      simpa only [pow_two] using prior.mul h (by omega)
    have correction : InClosure F (4*(n+1)) (h^2-normalForm G nz (h^2)) :=
      InClosure.transport budget
        (fun f hf hdeg => .input (included hf) (hdeg.trans budget))
        (normalForm_reduction G nz (h^2))
    have next := square.sub correction
    have id : h^2-(h^2-normalForm G nz (h^2)) = normalForm G nz (h^2) := by ring
    rw [id] at next
    dsimp only [h] at next
    rw [nf_square G nz z, ← pow_mul, ← pow_succ] at next
    simpa only [max_eq_right (by omega : 2 ≤ 4*(n+1))] using next
end ActualNormalForms

section DigitBounds
open BinaryWeight

theorem binary_weight_lt_two_pow (j n : ℕ) (hn : n < 2^j) : wt n ≤ j := by
  induction j generalizing n with
  | zero =>
    have : n = 0 := by simpa using hn
    subst n
    simp
  | succ j ih =>
    rw [wt_eq]
    have hd : n / 2 < 2^j := by rw [pow_succ] at hn; omega
    have := ih (n/2) hd
    have := Nat.mod_lt n (by omega : 0 < 2)
    omega

theorem binary_weight_mod_le (k n : ℕ) : wt (n % 2^k) ≤ wt n := by
  have h := wt_add_mul_two_pow k (n % 2^k) (n / 2^k)
    (Nat.mod_lt n (Nat.two_pow_pos k))
  rw [Nat.mod_add_div] at h
  omega

/-- The support estimate needed at squaring level j, including cyclic wrap. -/
theorem binary_weight_power_support (j k n : ℕ) (hn : n ≤ 2^(j+1)) :
    wt (n % 2^k) ≤ j+1 := by
  apply (binary_weight_mod_le k n).trans
  rcases lt_or_eq_of_le hn with h | rfl
  · exact binary_weight_lt_two_pow (j+1) n h
  · have h := wt_two_pow_add (j+1) 0 (Nat.two_pow_pos (j+1))
    simp only [Nat.add_zero, wt_zero, zero_add] at h
    omega
end DigitBounds

section PrimeInputPowers
open FiniteQuotient ChainQuotient ChainClosureBridge PrimeM2General
variable {L : Type*} [Field L]

/-- The real normalized two-summand input supplies the base and the chain
corrections. The sole remaining estimate is the support degree of its powers. -/
theorem prime_input_power_tower (two_ne : (2 : L) ≠ 0)
    (k : ℕ) (hk : 2 ≤ k) (a b t : L)
    (degree : ∀ j, (normalForm (chainInputs L 2 k) (chainInputs_nonzero 2 k)
      ((normalized k hk a b t)^(2^j))).totalDegree ≤ 2*(j+1))
    (n : ℕ) : InClosure (Set.range (inputs k hk a b t)) (max 2 (4*n))
      (normalForm (chainInputs L 2 k) (chainInputs_nonzero 2 k)
        ((normalized k hk a b t)^(2^n))) := by
  apply actual_normal_power_tower (chainInputs L 2 k)
    (Set.range (inputs k hk a b t)) (chainInputs_nonzero 2 k)
    (reduced_chain_ideal_trivial 2 k two_ne) (chainInputs_two_subset k hk a b t)
    (normalized k hk a b t) _ degree n
  rw [ChainModelMap.normalForm_normalized two_ne k hk a b t]
  have hd := (degree_and_top k hk a b t (Fin.last (2*k))).1.le
  have hd' : (normalized k hk a b t).totalDegree ≤ 2 := by
    simpa only [inputs, Fin.lastCases_last] using hd
  exact .input ⟨Fin.last (2*k), by simp only [inputs, Fin.lastCases_last]⟩ hd'
end PrimeInputPowers

section BinaryLift
open FiniteQuotient ChainTranslation ChainModelMap
open scoped BigOperators
variable {m : ℕ} {L : Type*} [Field L]

/-- Explicit termwise binary representative of a raw polynomial.
Its equality with graded-division normal form is not built into this definition. -/
def binaryLift (k : ℕ) (p : MvPolynomial (Fin m) L) : P (m*k) L :=
  ∑ e ∈ p.support,
    monomial (digits k (fun i => e i % 2^k)) (coeff e p)

theorem binaryLift_degree (k j : ℕ) (p : MvPolynomial (Fin m) L)
    (box : ∀ e ∈ p.support, ∀ i, e i ≤ 2^(j+1)) :
    (binaryLift k p).totalDegree ≤ m*(j+1) := by
  classical
  unfold binaryLift
  apply totalDegree_finsetSum_le
  intro e he
  apply (totalDegree_monomial_le _ _).trans
  change (digits k (fun i => e i % 2^k)).degree ≤ m*(j+1)
  rw [digits_degree _ (fun i => Nat.mod_lt _ (Nat.two_pow_pos k))]
  calc
    (∑ i : Fin m, BinaryWeight.wt (e i % 2^k)) ≤ ∑ _i : Fin m, (j+1) :=
      Finset.sum_le_sum (fun i _ => binary_weight_power_support j k (e i) (box e he i))
    _ = m*(j+1) := by simp

theorem binaryLift_reduced (k : ℕ) (p : MvPolynomial (Fin m) L) :
    binaryLift k p ∈ reducedSpace (ChainQuotient.chainInputs L m k) := by
  classical
  unfold binaryLift
  apply (reducedSpace _).sum_mem
  intro e _
  exact digits_reduced _ _

theorem s3_coordinate_degree (a b t : L) (i : Fin 2) :
    (GrowthReduction.S3poly a b t).degreeOf i ≤ 2 := by
  classical
  have hx : (X (0 : Fin 2) : MvPolynomial (Fin 2) L).degreeOf i ≤ 1 := by
    rw [degreeOf_X]; split_ifs <;> omega
  have hy : (X (1 : Fin 2) : MvPolynomial (Fin 2) L).degreeOf i ≤ 1 := by
    rw [degreeOf_X]; split_ifs <;> omega
  have hxy : ((X 0 : MvPolynomial (Fin 2) L)*X 1).degreeOf i ≤ 1 := by
    apply (degreeOf_mul_le i _ _).trans
    fin_cases i <;> norm_num [degreeOf_X]
  have hsum : ((X 0 : MvPolynomial (Fin 2) L)+X 1).degreeOf i ≤ 1 :=
    (degreeOf_add_le i _ _).trans (max_le hx hy)
  have hsub : ((X 0 : MvPolynomial (Fin 2) L)-X 1).degreeOf i ≤ 1 :=
    (degreeOf_sub_le i _ _).trans (max_le hx hy)
  have hpa : ((X 0*X 1 : MvPolynomial (Fin 2) L)+C a).degreeOf i ≤ 1 := by
    exact (degreeOf_add_le i _ _).trans (max_le hxy (by simp))
  have hma : ((X 0*X 1 : MvPolynomial (Fin 2) L)-C a).degreeOf i ≤ 1 := by
    exact (degreeOf_sub_le i _ _).trans (max_le hxy (by simp))
  have hm : (((X 0+X 1)*(X 0*X 1+C a) : MvPolynomial (Fin 2) L)+C (2*b)).degreeOf i ≤ 2 := by
    apply (degreeOf_add_le i _ _).trans
    apply max_le
    · exact (degreeOf_mul_le i _ _).trans (by omega)
    · simp only [degreeOf_C]
      omega
  have hfirst : (((X 0-X 1)^2 : MvPolynomial (Fin 2) L)*C (t^2)).degreeOf i ≤ 2 := by
    apply (degreeOf_mul_C_le _ i _).trans
    exact (degreeOf_pow_le i _ 2).trans (by omega)
  have hsecond : ((2 : MvPolynomial (Fin 2) L)*((X 0+X 1)*(X 0*X 1+C a)+C (2*b))*C t).degreeOf i ≤ 2 := by
    apply (degreeOf_mul_C_le _ i _).trans
    change (C (2:L)*_).degreeOf i ≤ 2
    exact (degreeOf_C_mul_le _ i _).trans hm
  have hthird : (((X 0*X 1-C a)^2 : MvPolynomial (Fin 2) L)).degreeOf i ≤ 2 :=
    (degreeOf_pow_le i _ 2).trans (by omega)
  have hlast : ((C (4*b)*(X 0+X 1) : MvPolynomial (Fin 2) L)).degreeOf i ≤ 2 :=
    (degreeOf_C_mul_le _ i _).trans (hsum.trans (by omega))
  unfold GrowthReduction.S3poly
  exact (degreeOf_sub_le i _ _).trans (max_le
    ((degreeOf_add_le i _ _).trans (max_le
      ((degreeOf_sub_le i _ _).trans (max_le hfirst hsecond)) hthird)) hlast)

theorem binaryLift_s3_power_degree (k j : ℕ) (a b t : L) :
    (binaryLift k ((GrowthReduction.S3poly a b t)^(2^j))).totalDegree ≤ 2*(j+1) := by
  apply binaryLift_degree
  intro e he i
  have h := (degreeOf_le_iff.mp (degreeOf_pow_le i (GrowthReduction.S3poly a b t) (2^j))) e he
  have h2 := s3_coordinate_degree a b t i
  have h3 : 2^j * (GrowthReduction.S3poly a b t).degreeOf i ≤ 2^j*2 :=
    Nat.mul_le_mul_left _ h2
  rw [← pow_succ] at h3
  exact h.trans h3
end BinaryLift

#print axioms map_closure_budget
#print axioms not_one_of_specialization
#print axioms no_lfd_bound_of_specialization
#print axioms prime_growth_transfer
#print axioms order_two_pair_specialization
#print axioms pairIterate_closed
#print axioms card_sumTargets_le
#print axioms exists_target_outside_sums
#print axioms normalized_square_tower
#print axioms nf_square
#print axioms actual_normal_power_tower
#print axioms binary_weight_power_support
#print axioms prime_input_power_tower
#print axioms binaryLift_degree
#print axioms binaryLift_reduced
#print axioms binaryLift_s3_power_degree
end ICFallDegrees.SummandGrowthTransfer
