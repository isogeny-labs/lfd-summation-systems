import Mathlib.RingTheory.MvPolynomial.Basic
import Mathlib.Algebra.MvPolynomial.CommRing
import Mathlib.Algebra.CharP.Two
import Mathlib.LinearAlgebra.Basis.Defs
import Mathlib.Tactic.Ring

/-! Coefficient descent through an actual basis, and explicit Boolean
preprocessing identities for the characteristic-two S3 polynomial. -/
noncomputable section
namespace ICFallDegrees.BinaryDescent
open MvPolynomial
open scoped BigOperators

section Coefficients
variable {σ k K : Type*} {n : ℕ} [CommRing k] [CommRing K] [Algebra k K]

def component (basis : Basis (Fin n) k K) (i : Fin n) :
    MvPolynomial σ K →ₗ[k] MvPolynomial σ k where
  toFun p := Finsupp.mapRange (basis.coord i) (basis.coord i).map_zero p
  map_add' p q := by
    ext e
    change (basis.coord i) (coeff e (p+q))=
      (basis.coord i) (coeff e p)+(basis.coord i) (coeff e q)
    rw [coeff_add,map_add]
  map_smul' c p := by
    ext e
    change (basis.coord i) (coeff e (c • p))=c • ((basis.coord i) (coeff e p))
    rw [coeff_smul,map_smul]

theorem component_scaled_map (basis : Basis (Fin n) k K) (i : Fin n)
    (a : K) (p : MvPolynomial σ k) :
    component basis i (C a*map (algebraMap k K) p)=C (basis.coord i a)*p := by
  ext e
  change (basis.coord i) (coeff e (C a*map (algebraMap k K) p))=
    coeff e (C (basis.coord i a)*p)
  rw [coeff_C_mul,coeff_map,coeff_C_mul]
  calc (basis.coord i) (a*algebraMap k K (coeff e p))=
      (basis.coord i) ((coeff e p) • a) := by rw [Algebra.smul_def,mul_comm]
       _=(coeff e p) • (basis.coord i a) := (basis.coord i).map_smul _ _
       _=(basis.coord i a)*(coeff e p) := by rw [smul_eq_mul,mul_comm]

theorem component_of_expansion (basis : Basis (Fin n) k K)
    (coordinates : Fin n → MvPolynomial σ k) (p : MvPolynomial σ K)
    (expansion : p=∑ j, C (basis j)*map (algebraMap k K) (coordinates j)) (i : Fin n) :
    component basis i p=coordinates i := by
  classical
  rw [expansion,map_sum]
  simp [component_scaled_map,Basis.coord_apply,Finsupp.single_apply,apply_ite]

end Coefficients

section Normalization
variable {K : Type*} [CommRing K] [CharP K 2]
abbrev P (K : Type*) [CommRing K] := MvPolynomial (Fin 4) K

def x (a : K) : P K := X 0+C a*X 1
def y (a : K) : P K := X 2+C a*X 3
def xbar (a : K) : P K := X 0+(C a)^2*X 1
def ybar (a : K) : P K := X 2+(C a)^2*X 3
def boolean (i : Fin 4) : P K := X i^2-X i

def s3 (a t b : K) : P K :=
  (x a*y a)^2+(x a*C t)^2+(y a*C t)^2+x a*y a*C t+C b

def normalized (a t b : K) : P K :=
  xbar a*ybar a+(C t)^2*(xbar a+ybar a)+C t*x a*y a+C b

/-- Literal raw-to-normalized correction. Degree-four products are displayed;
this identity is not a free degree reduction inside a certificate. -/
theorem normalization_identity (a t b : K) : s3 a t b-normalized a t b=
    ((y a)^2+(C t)^2)*(boolean 0+(C a)^2*boolean 1)+
      (xbar a+(C t)^2)*(boolean 2+(C a)^2*boolean 3) := by
  have identity : s3 a t b-normalized a t b-
      (((y a)^2+(C t)^2)*(boolean 0+(C a)^2*boolean 1)+
       (xbar a+(C t)^2)*(boolean 2+(C a)^2*boolean 3))=
      2*(C a*X 0*X 1*((y a)^2+(C t)^2)+
        C a*X 2*X 3*(xbar a+(C t)^2)) := by
    simp only [s3,normalized,x,y,xbar,ybar,boolean]
    ring
  rw [CharTwo.two_eq_zero (R:=P K),zero_mul] at identity
  exact sub_eq_zero.mp identity

/-- The usual squared-sum S3 and its separated-squares expression agree in
characteristic two. -/
theorem s3_squared_sum (a t b : K) :
    (x a*y a+x a*C t+y a*C t)^2+x a*y a*C t+C b=s3 a t b := by
  have identity : (x a*y a+x a*C t+y a*C t)^2+x a*y a*C t+C b-s3 a t b=
      2*(x a*y a*x a*C t+x a*y a*y a*C t+x a*C t*y a*C t) := by
    unfold s3
    ring
  rw [CharTwo.two_eq_zero (R:=P K),zero_mul] at identity
  exact sub_eq_zero.mp identity

end Normalization
end ICFallDegrees.BinaryDescent
