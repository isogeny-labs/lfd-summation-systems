import ICFallDegrees.CyclicInverseRewrite
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

/-! The all-length even alternating obstruction has an explicit polynomial
bridge. Every multiplier is charged before normalization. No trace, rank,
normalization oracle, or equality of V-spaces is assumed. -/
noncomputable section
namespace ICFallDegrees.CyclicProductBridge
open MvPolynomial
open scoped BigOperators
variable {σ K : Type*} [CommRing K]

private theorem product_degree (a : ℕ → (MvPolynomial σ K)) (ha : ∀ i, (a i).totalDegree≤1) (n : ℕ) :
    (∏ i∈Finset.range n, a i).totalDegree≤n := by
  apply (totalDegree_finset_prod _ _).trans
  calc
    (∑ i∈Finset.range n, (a i).totalDegree)≤∑ _i∈Finset.range n, 1 :=
      Finset.sum_le_sum (fun i _ => ha i)
    _=n := by simp

private theorem path_degree (a b c : ℕ → (MvPolynomial σ K))
    (ha : ∀ i, (a i).totalDegree≤1) (hb : ∀ i, (b i).totalDegree≤1)
    (hc : ∀ i, (c i).totalDegree≤1) (n : ℕ) :
    (c 0*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*c n).totalDegree≤n+1 := by
  apply (totalDegree_sub_le _ _).trans
  apply max_le
  · exact (totalDegree_mul _ _).trans (by have h := product_degree a ha n; have j := hc 0; omega)
  · exact (totalDegree_mul _ _).trans (by have h := product_degree b hb n; have j := hc n; omega)

/-- A chain of quadratic carries lifts to arbitrary length at its ordinary
product degree, retaining the degree-three setup cost. -/
theorem conveyor (F : Set (MvPolynomial σ K)) (a b c : ℕ → (MvPolynomial σ K))
    (ha : ∀ i, (a i).totalDegree≤1) (hb : ∀ i, (b i).totalDegree≤1)
    (hc : ∀ i, (c i).totalDegree≤1)
    (rows : ∀ i, InClosure F 3 (c i*a i-b i*c (i+1))) (n : ℕ) :
    InClosure F (max 3 (n+1))
      (c 0*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*c n) := by
  induction n with
  | zero => simp only [Finset.range_zero,Finset.prod_empty,mul_one,one_mul,sub_self]; exact .zero
  | succ n ih =>
      have prior := ih.mono (by omega : max 3 (n+1)≤max 3 (n+1+1))
      have degree := path_degree a b c ha hb hc n
      have next := prior.mul (a n) (by have h := ha n; omega)
      have quadratic : (c n*a n-b n*c (n+1)).totalDegree≤2 := by
        apply (totalDegree_sub_le _ _).trans
        exact max_le ((totalDegree_mul _ _).trans (by have h := hc n; have j := ha n; omega))
          ((totalDegree_mul _ _).trans (by have h := hb n; have j := hc (n+1); omega))
      have correction := CyclicInverseRewrite.multiply_quadratic F _
        (∏ i∈Finset.range n, b i) (max 3 (n+1+1)) (Nat.le_max_left _ _) (rows n) quadratic
        (by have h := product_degree b hb n; omega)
      have identity :
          c 0*(∏ i∈Finset.range (n+1), a i)-(∏ i∈Finset.range (n+1), b i)*c (n+1)=
          a n*(c 0*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*c n)+
          (∏ i∈Finset.range n, b i)*(c n*a n-b n*c (n+1)) := by
        rw [Finset.prod_range_succ,Finset.prod_range_succ]
        ring
      rw [identity]
      exact next.add correction

theorem endpoint_bridge (F : Set (MvPolynomial σ K)) (a b c : ℕ → (MvPolynomial σ K)) (u w last : (MvPolynomial σ K))
    (ha : ∀ i, (a i).totalDegree≤1) (hb : ∀ i, (b i).totalDegree≤1)
    (hc : ∀ i, (c i).totalDegree≤1)
    (hu : u.totalDegree≤1) (hw : w.totalDegree≤1) (hl : last.totalDegree≤1)
    (first : InClosure F 3 (u-c 0*w))
    (rows : ∀ i, InClosure F 3 (c i*a i-b i*c (i+1)))
    (n : ℕ) (final : InClosure F 3 (c n*w-last)) :
    InClosure F (max 3 (n+2))
      (u*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*last) := by
  have first_degree : (u-c 0*w).totalDegree≤2 := by
    exact (totalDegree_sub_le _ _).trans (max_le (by omega)
      ((totalDegree_mul _ _).trans (by have h := hc 0; omega)))
  have final_degree : (c n*w-last).totalDegree≤2 := by
    exact (totalDegree_sub_le _ _).trans (max_le
      ((totalDegree_mul _ _).trans (by have h := hc n; omega)) (by omega))
  have initial := CyclicInverseRewrite.multiply_quadratic F _
    (∏ i∈Finset.range n, a i) (max 3 (n+2)) (Nat.le_max_left _ _) first first_degree
    (by have h := product_degree a ha n; omega)
  have ending := CyclicInverseRewrite.multiply_quadratic F _
    (∏ i∈Finset.range n, b i) (max 3 (n+2)) (Nat.le_max_left _ _) final final_degree
    (by have h := product_degree b hb n; omega)
  have middle := (conveyor F a b c ha hb hc rows n).mono
    (by omega : max 3 (n+1)≤max 3 (n+2))
  have lifted := middle.mul w (by have h := path_degree a b c ha hb hc n; omega)
  have identity : u*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*last=
      (∏ i∈Finset.range n, a i)*(u-c 0*w)+
      w*(c 0*(∏ i∈Finset.range n, a i)-(∏ i∈Finset.range n, b i)*c n)+
      (∏ i∈Finset.range n, b i)*(c n*w-last) := by ring
  rw [identity]
  exact (initial.add lifted).add ending

/-- For a cycle of length 2(n+1), the two alternating degree-(n+1)
monomials differ by a derivation at max(3,n+2). The premise lists only
square-chain and inverse input derivations and the closing cyclic identity. -/
theorem even_alternating_bridge (F : Set (MvPolynomial σ K)) (u v : ℕ → (MvPolynomial σ K))
    (hu : ∀ i, (u i).totalDegree≤1) (hv : ∀ i, (v i).totalDegree≤1)
    (su : ∀ i, InClosure F 2 ((u i)^2-u (i+1)))
    (sv : ∀ i, InClosure F 2 ((v i)^2-v (i+1)))
    (inverse : ∀ i, InClosure F 2 (u i*v i-1))
    (n : ℕ) (closed : v (2*(n+1))=v 0) :
    InClosure F (max 3 (n+2))
      ((∏ i∈Finset.range (n+1), u (2*i))-(∏ i∈Finset.range (n+1), v (2*i+1))) := by
  let a : ℕ → (MvPolynomial σ K) := fun i => u (2*i+2)
  let b : ℕ → (MvPolynomial σ K) := fun i => v (2*i+1)
  let c : ℕ → (MvPolynomial σ K) := fun i => u (2*i+1)
  have first : InClosure F 3 (u 0-c 0*v 0) := by
    simpa [c,neg_sub] using (CyclicInverseRewrite.carry F (u 0) (v 0) (u 1) (hu 0) (hv 0) (hu 1)
      (su 0) (inverse 0)).neg
  have rows : ∀ i, InClosure F 3 (c i*a i-b i*c (i+1)) := by
    intro i
    have carry := CyclicInverseRewrite.carry F (u (2*i+1)) (v (2*i+1)) (u (2*i+1+1))
      (hu _) (hv _) (hu _) (su _) (inverse _)
    have next := CyclicInverseRewrite.adjacent F (u (2*i+1)) (v (2*i+1))
      (u (2*i+1+1)) (u (2*i+1+1+1)) (hu _) (hv _) (hu _) (hu _) carry (su _)
    simpa [a,b,c,Nat.mul_add,Nat.add_assoc] using next
  have final : InClosure F 3 (c n*v 0-v (2*n+1)) := by
    have reversed : InClosure F 2 (v (2*n+1)*u (2*n+1)-1) := by simpa [mul_comm] using inverse (2*n+1)
    have carry := CyclicInverseRewrite.carry F (v (2*n+1)) (u (2*n+1)) (v (2*n+1+1))
      (hv _) (hu _) (hv _) (sv _) reversed
    have index : 2*n+1+1=2*(n+1) := by omega
    rw [index,closed] at carry
    simpa [c,mul_comm] using carry
  have bridge := endpoint_bridge F a b c (u 0) (v 0) (v (2*n+1))
    (fun i => hu _) (fun i => hv _) (fun i => hu _) (hu _) (hv _) (hv _) first rows n final
  have left : u 0*(∏ i∈Finset.range n, a i)=∏ i∈Finset.range (n+1), u (2*i) := by
    rw [Finset.prod_range_succ']
    simp only [a,Nat.mul_add,Nat.mul_one,Nat.mul_zero]
    ring
  have right : (∏ i∈Finset.range n, b i)*v (2*n+1)=∏ i∈Finset.range (n+1), v (2*i+1) := by
    rw [Finset.prod_range_succ]
  rw [left,right] at bridge
  exact bridge

end ICFallDegrees.CyclicProductBridge
