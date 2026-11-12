import ICFallDegrees.PrimeM2General
import ICFallDegrees.SquareChainRefutation
import ICFallDegrees.QuadraticNoLow

/-! A growing-depth exact LFD family from actual normalized S3 equations
for y²=x³+a*x and target (0,0). The theorem applies to prime-field subgroup
instances whenever their usual geometric/subgroup hypotheses hold. -/
noncomputable section
open Classical
namespace ICFallDegrees.PrimeM2Sharp
open MvPolynomial HomogeneousFirstFall GeneralQuadraticTop GeneralQuadraticClassification
open PrimeM2General
variable {K : Type*} [Field K]
abbrev P (k : ℕ) (K : Type*) [Field K] := PrimeM2General.P k K

def index (k : ℕ) (b : Bool) (j : Fin k) : Fin (2*k) :=
  ⟨if b then k+j.val else j.val,by cases b <;> simp only [Bool.false_eq_true,↓reduceIte]; all_goals omega⟩

def value (k : ℕ) (b : Bool) (j : ℕ) : P k K :=
  if small : j<k then X (index k b ⟨j,small⟩) else 1

def product (k : ℕ) (j : ℕ) : P k K := value k false j*value k true j

theorem value_degree (k : ℕ) (b : Bool) (j : ℕ) : (value k b j : P k K).totalDegree≤1 := by
  unfold value
  split_ifs <;> simp

theorem product_degree (k : ℕ) (j : ℕ) : (product k j : P k K).totalDegree≤2 := by
  exact (totalDegree_mul _ _).trans (Nat.add_le_add (value_degree _ _ _) (value_degree _ _ _))

theorem chain_identity (k : ℕ) (depth : 2≤k) (a : K) (b : Bool) (j : ℕ) (small : j<k) :
    value k b j^2-value k b (j+1)=inputs k depth a 0 0 (index k b ⟨j,small⟩).castSucc := by
  simp only [inputs,Fin.lastCases_castSucc]
  by_cases next : j+1<k
  · have e1 : j+1≠k := by omega
    have e2 : j+1≠2*k := by omega
    have e3 : k+(j+1)≠k := by omega
    have e4 : k+(j+1)≠2*k := by omega
    cases b <;> simp [value,index,inputs,chainTail,small,next,e1,e2,e3,e4,sub_eq_add_neg,Nat.add_assoc]
  · have last : j+1=k := by omega
    cases b <;> simp [value,index,inputs,chainTail,small,next,last,sub_eq_add_neg,Nat.add_assoc,two_mul]

theorem chain_at_four (k : ℕ) (depth : 2≤k) (a : K) (b : Bool) (j : ℕ) (small : j<k) :
    InClosure (Set.range (inputs k depth a 0 0)) 4 (value k b j^2-value k b (j+1)) := by
  rw [chain_identity k depth a b j small]
  exact .input ⟨_,rfl⟩ (by rw [(degree_and_top k depth a 0 0 _).1]; decide)

theorem product_defect (k : ℕ) (depth : 2≤k) (a : K) (j : ℕ) (bounded : j≤k) :
    InClosure (Set.range (inputs k depth a 0 0)) 4 (product k j^2-product k (j+1)) := by
  by_cases small : j<k
  · have first := chain_at_four k depth a false j small
    have second := chain_at_four k depth a true j small
    have df : (value k false j^2-value k false (j+1) : P k K).totalDegree≤2 := by
      exact (totalDegree_sub_le _ _).trans (max_le
        ((totalDegree_pow _ _).trans (Nat.mul_le_mul_left 2 (value_degree _ _ _)))
        ((value_degree _ _ _).trans (by decide)))
    have ds : (value k true j^2-value k true (j+1) : P k K).totalDegree≤2 := by
      exact (totalDegree_sub_le _ _).trans (max_le
        ((totalDegree_pow _ _).trans (Nat.mul_le_mul_left 2 (value_degree _ _ _)))
        ((value_degree _ _ _).trans (by decide)))
    have first_mult := first.mul (value k true j^2) (by
      have h := (totalDegree_pow (value k true j : P k K) 2).trans
        (Nat.mul_le_mul_left 2 (value_degree k true j))
      omega)
    have second_mult := second.mul (value k false (j+1)) (by
      have := value_degree (K:=K) k false (j+1)
      omega)
    have identity : (product k j : P k K)^2-product k (j+1)=
        value k true j^2*(value k false j^2-value k false (j+1))+
          value k false (j+1)*(value k true j^2-value k true (j+1)) := by unfold product; ring
    rw [identity]
    exact first_mult.add second_mult
  · have equal : j=k := by omega
    subst j
    simpa [product,value] using (InClosure.zero : InClosure (Set.range (inputs k depth a 0 0)) 4 0)

theorem initial_identity (k : ℕ) (depth : 2≤k) (a : K) :
    SquareChainRefutation.row (product k) a 0=normalized k depth a 0 0 := by
  have h0 : 0<k := by omega
  have h1 : 1<k := by omega
  simp [SquareChainRefutation.row,product,value,index,active,normalized,top,tail,h0,h1]
  ring

theorem one_at_four (k : ℕ) (depth : 2≤k) (a : K) (outside : a^(2^k)≠1) :
    InClosure (Set.range (inputs k depth a 0 0)) 4 1 := by
  apply SquareChainRefutation.one_at_four _ (product k) (product_degree k) a k
  · exact product_defect k depth a
  · rw [initial_identity k depth a]
    have input : normalized k depth a 0 0=inputs k depth a 0 0 (Fin.last (2*k)) := by simp [inputs]
    rw [input]
    exact .input ⟨_,rfl⟩ (by rw [(degree_and_top k depth a 0 0 _).1]; decide)
  · simp [product,value]
  · simp [product,value]
  · exact outside

theorem no_one_at_three (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a : K) (nonzero : a≠0) :
    ¬InClosure (Set.range (inputs k depth a 0 0)) 3 1 := by
  let tails : Fin (2*k+1) → P k K := Fin.lastCases (tail k depth a 0 0) (chainTail k)
  have presentation : inputs k depth a 0 0=(fun i => q (active k depth) (-2*a) 0 i+tails i) := by
    funext i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · simp [inputs,normalized,tails,q]
    · simp [inputs,q,tails]
  rw [presentation]
  apply QuadraticNoLow.no_one_at_three
  · exact homogeneous_q _ _ _
  · intro i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · have bound := tail_small k depth a 0 0
      simpa only [tails,Fin.lastCases_last] using (show (tail k depth a 0 0).totalDegree≤1 by omega)
    · have bound := chain_tail_small (K:=K) k i
      simpa only [tails,Fin.lastCases_castSucc] using (show (chainTail k i : P k K).totalDegree≤1 by omega)
  · intro h relation
    exact low_syzygy_zero _ (pivot (active k depth))
      (coefficient_pivots _ (active_injective k depth) (-2*a) 0) 2 2 le_rfl relation
  · intro h relation
    apply degree_three_syzygy_zero _ (active_injective k depth) (-2*a) 0 _ relation
    simpa only [zero_pow (by decide : 2≠0),sub_zero] using mul_ne_zero (neg_ne_zero.mpr two_ne) nonzero

/-- Matching bounds for every chain depth k≥2, on a narrow order-two-target
curve family. This is not a sharpness theorem for the general mk+1 bound. -/
theorem exact_lfd_four (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k)
    (a : K) (nonzero : a≠0) (outside : a^(2^k)≠1) :
    IsLastFallDegree (Set.range (inputs k depth a 0 0)) 4 := by
  apply exactLFD_of_first_one (by decide) (one_at_four k depth a outside)
  exact no_one_at_three two_ne k depth a nonzero

end ICFallDegrees.PrimeM2Sharp
