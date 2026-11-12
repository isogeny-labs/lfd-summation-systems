import ICFallDegrees.BinaryDescent
import ICFallDegrees.GridQuotient
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Algebra.Field.ZMod

/-! Target-independent bounds for arbitrary collections of binary S3
constraints, including split encodings with any number of summands. All
Boolean equations are explicit; basis extraction does not increase degree. -/
noncomputable section
namespace ICFallDegrees.BinaryS3Bounds
open MvPolynomial BinaryDescent GridQuotient
open scoped BigOperators
variable {r n s : ℕ} {K : Type*} [Field K] [Algebra (ZMod 2) K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

theorem component_degree (basis : Basis (Fin n) (ZMod 2) K) (i : Fin n) (f : P r K) :
    (component basis i f).totalDegree≤f.totalDegree := by
  classical
  apply Finset.sup_le
  intro e he
  apply le_totalDegree
  by_contra absent
  have zero : coeff e f=0 := not_mem_support_iff.mp absent
  have extracted : coeff e (component basis i f)=0 := by
    change basis.coord i (coeff e f)=0
    rw [zero,map_zero]
  exact (mem_support_iff.mp he) extracted

def affine (c : K) (a : Fin r → K) : P r K := C c+∑ i, C (a i)*X i

theorem affine_degree (c : K) (a : Fin r → K) : (affine c a).totalDegree≤1 := by
  apply (totalDegree_add _ _).trans
  apply max_le (by rw [totalDegree_C]; omega)
  apply totalDegree_finsetSum_le
  intro i _
  exact ((isHomogeneous_X K i).C_mul (a i)).totalDegree_le

/-- In characteristic two, squaring an affine coordinate form reduces
by an explicit combination of the Boolean input equations. -/
theorem affine_square_identity [CharP K 2] (c : K) (a : Fin r → K) :
    (affine c a)^2-affine (c^2) (fun i => (a i)^2)=
      ∑ i, C ((a i)^2)*((X i : P r K)^2-X i) := by
  classical
  unfold affine
  rw [add_pow_char,sum_pow_char]
  simp only [mul_pow,←map_pow]
  rw [add_sub_add_left_eq_sub,←Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

def raw (x y z : P r K) (b : K) : P r K :=
  (x*y)^2+(x*z)^2+(y*z)^2+x*y*z+C b

def normalized (x y z xb yb zb : P r K) (b : K) : P r K :=
  xb*yb+xb*zb+yb*zb+x*y*z+C b

private theorem product_degree {x y : P r K} {a b : ℕ}
    (hx : x.totalDegree≤a) (hy : y.totalDegree≤b) : (x*y).totalDegree≤a+b :=
  (totalDegree_mul _ _).trans (Nat.add_le_add hx hy)

private theorem add_degree {x y : P r K} {a : ℕ}
    (hx : x.totalDegree≤a) (hy : y.totalDegree≤a) : (x+y).totalDegree≤a :=
  (totalDegree_add _ _).trans (max_le hx hy)

theorem raw_degree (x y z : P r K) (b : K)
    (hx : x.totalDegree≤1) (hy : y.totalDegree≤1) (hz : z.totalDegree≤1) :
    (raw x y z b).totalDegree≤4 := by
  have square {u : P r K} (hu : u.totalDegree≤2) : (u^2).totalDegree≤4 := by
    have h := totalDegree_pow u 2
    omega
  unfold raw
  exact add_degree (add_degree (add_degree
    (add_degree (square (product_degree hx hy)) (square (product_degree hx hz)))
    (square (product_degree hy hz))) ((product_degree (product_degree hx hy) hz).trans (by decide)))
    (by rw [totalDegree_C]; omega)

theorem normalized_degree (x y z xb yb zb : P r K) (b : K)
    (hx : x.totalDegree≤1) (hy : y.totalDegree≤1) (hz : z.totalDegree≤1)
    (hxb : xb.totalDegree≤1) (hyb : yb.totalDegree≤1) (hzb : zb.totalDegree≤1) :
    (normalized x y z xb yb zb b).totalDegree≤3 := by
  unfold normalized
  exact add_degree (add_degree (add_degree
    (add_degree ((product_degree hxb hyb).trans (by decide))
      ((product_degree hxb hzb).trans (by decide)))
    ((product_degree hyb hzb).trans (by decide))) (product_degree (product_degree hx hy) hz))
    (by rw [totalDegree_C]; omega)

/-- Displaying the square defects keeps raw preprocessing separate from
an uncharged Boolean reduction inside a closure derivation. -/
theorem normalization_identity (x y z xb yb zb : P r K) (b : K) :
    raw x y z b-normalized x y z xb yb zb b=
      (y^2+z^2)*(x^2-xb)+(xb+z^2)*(y^2-yb)+(xb+yb)*(z^2-zb) := by
  unfold raw normalized
  ring

def rawConstraint (constants : Fin s → Fin 3 → K)
    (coefficients : Fin s → Fin 3 → Fin r → K) (b : K) (j : Fin s) : P r K :=
  raw (affine (constants j 0) (coefficients j 0))
    (affine (constants j 1) (coefficients j 1))
    (affine (constants j 2) (coefficients j 2)) b

def normalizedConstraint (constants : Fin s → Fin 3 → K)
    (coefficients : Fin s → Fin 3 → Fin r → K) (b : K) (j : Fin s) : P r K :=
  normalized (affine (constants j 0) (coefficients j 0))
    (affine (constants j 1) (coefficients j 1))
    (affine (constants j 2) (coefficients j 2))
    (affine ((constants j 0)^2) (fun i => (coefficients j 0 i)^2))
    (affine ((constants j 1)^2) (fun i => (coefficients j 1 i)^2))
    (affine ((constants j 2)^2) (fun i => (coefficients j 2 i)^2)) b

def descendedInputs (basis : Basis (Fin n) (ZMod 2) K)
    (constraints : Fin s → P r K) : Set (P r (ZMod 2)) :=
  Set.range (fun i : Fin r => X i^2-X i) ∪
    {p | ∃ j i, p=component basis i (constraints j)}

private theorem boolean_degree (i : Fin r) : ((X i : P r (ZMod 2))^2-X i).totalDegree≤2 := by
  exact (totalDegree_sub_le _ _).trans (by simp)

/-- Every target and every S3 incidence graph are allowed. In a left-chain
m-summand encoding one takes s=m-1 and r=sum(base dimensions)+(m-2)*n. -/
theorem raw_split_lfd_bound (basis : Basis (Fin n) (ZMod 2) K)
    (constants : Fin s → Fin 3 → K) (coefficients : Fin s → Fin 3 → Fin r → K) (b : K) :
    HasLFDBound (descendedInputs basis (rawConstraint constants coefficients b)) (max (r+1) 4) := by
  apply boolean_lfd_bound _ (fun i => Or.inl ⟨i,rfl⟩) 4
  rintro f (⟨i,rfl⟩ | ⟨j,i,rfl⟩)
  · exact (boolean_degree i).trans (by decide)
  · exact (component_degree basis i _).trans (raw_degree _ _ _ b
      (affine_degree _ _) (affine_degree _ _) (affine_degree _ _))

theorem normalized_split_lfd_bound (basis : Basis (Fin n) (ZMod 2) K)
    (constants : Fin s → Fin 3 → K) (coefficients : Fin s → Fin 3 → Fin r → K) (b : K) :
    HasLFDBound (descendedInputs basis (normalizedConstraint constants coefficients b)) (max (r+1) 3) := by
  apply boolean_lfd_bound _ (fun i => Or.inl ⟨i,rfl⟩) 3
  rintro f (⟨i,rfl⟩ | ⟨j,i,rfl⟩)
  · exact (boolean_degree i).trans (by decide)
  · exact (component_degree basis i _).trans (normalized_degree _ _ _ _ _ _ b
      (affine_degree _ _) (affine_degree _ _) (affine_degree _ _)
      (affine_degree _ _) (affine_degree _ _) (affine_degree _ _))

end ICFallDegrees.BinaryS3Bounds
