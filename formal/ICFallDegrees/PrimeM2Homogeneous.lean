import ICFallDegrees.QuadraticHomogeneousFamily

/-! Exact ordinary Koszul first fall for the derivative-zero branch of the
normalized two-summand S3 system, with three power-chain variables per
summand. Coefficients and the coefficient ring are arbitrary. Geometric
curve admissibility may be imposed separately; the algebraic theorem is
stronger and does not assume it. -/
noncomputable section
namespace ICFallDegrees.PrimeM2Homogeneous
open MvPolynomial HomogeneousFirstFall
variable {K : Type*} [CommRing K] [Nontrivial K]
abbrev P (K : Type*) [CommRing K] := MvPolynomial (Fin 6) K

def normalized (a b t : K) : P K :=
  X 1*X 4+C (-2*t)*(X 1*X 3+X 0*X 4)+C (-2*(t^2+a))*X 0*X 3+
    C (t^2)*(X 1+X 4)+C (-2*a*t-4*b)*(X 0+X 3)+C (a^2-4*b*t)

def tail (a b t : K) : P K :=
  C (t^2)*(X 1+X 4)+C (-2*a*t-4*b)*(X 0+X 3)+C (a^2-4*b*t)

def inputs (a b t : K) : Fin 7 → P K := ![X 0^2-X 1,X 1^2-X 2,X 2^2-1,
  X 3^2-X 4,X 4^2-X 5,X 5^2-1,normalized a b t]

def rawS3 (a b t : K) : P K :=
  (X 0-X 3)^2*C (t^2)-2*((X 0+X 3)*(X 0*X 3+C a)+C (2*b))*C t+
    (X 0*X 3-C a)^2-C (4*b)*(X 0+X 3)

theorem raw_normalization_identity (a b t : K) : rawS3 a b t-normalized a b t=
    (X 3-C t)^2*(X 0^2-X 1)+
      (X 1-C (2*t)*X 0+C (t^2))*(X 3^2-X 4) := by
  simp only [rawS3,normalized,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring

theorem tail_small (a b t : K) : (tail a b t).totalDegree<2 := by
  have first := ((isHomogeneous_X K (1 : Fin 6)).add (isHomogeneous_X K 4)).C_mul (t^2)
  have second := ((isHomogeneous_X K (0 : Fin 6)).add (isHomogeneous_X K 3)).C_mul (-2*a*t-4*b)
  have bound : (tail a b t).totalDegree≤1 := by
    unfold tail
    apply (totalDegree_add _ _).trans
    apply max_le
    · exact (totalDegree_add _ _).trans (max_le first.totalDegree_le second.totalDegree_le)
    · rw [totalDegree_C]
      decide
  omega

theorem derivative_zero_family_identity (a b t : K) (condition : 3*t^2+a=0) :
    QuadraticHomogeneousFamily.system (-2*t) (tail a b t)=inputs a b t := by
  have aeq : a= -(3*t^2) := by
    calc a=(3*t^2+a)-3*t^2 := by ring
         _= -(3*t^2) := by rw [condition]; simp
  have coefficient : -2*(t^2+a)=(-2*t)^2 := by rw [aeq]; ring
  funext i
  fin_cases i <;>
    simp [QuadraticHomogeneousFamily.system,QuadraticHomogeneousFamily.q,
      QuadraticHomogeneousFamily.lower,QuadraticHomogeneousFamily.top,
      inputs,normalized,tail,coefficient,map_pow,sub_eq_add_neg] <;> ring

/-- Fixed m=2,k=3, arbitrary coefficient ring and coefficients: under the
explicit derivative-zero condition, the actual normalized S3 presentation
has exact ordinary homogeneous-syzygy first fall three modulo Koszul. -/
theorem derivative_zero_ffd_three (a b t : K) (condition : 3*t^2+a=0) :
    IsKoszulFirstFall (inputs a b t) 3 := by
  rw [←derivative_zero_family_identity a b t condition]
  exact QuadraticHomogeneousFamily.exact_first_fall_three (-2*t) (tail a b t) (tail_small a b t)

end ICFallDegrees.PrimeM2Homogeneous
