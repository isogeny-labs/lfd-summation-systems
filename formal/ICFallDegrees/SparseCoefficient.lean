import ICFallDegrees.SparsePolynomial

namespace ICFallDegrees.SparsePolynomial

open MvPolynomial

variable {r : ℕ} {K : Type*} [CommRing K] [DecidableEq K]

def coefficient (e : Exponent r) : Sparse r K → K
  | [] => 0
  | (f,c)::p => (if sameExponent f e then c else 0)+coefficient e p

def sum : List (Sparse r K) → Sparse r K
  | [] => []
  | p::ps => add p (sum ps)

def variableData (i : Fin r) : Sparse r K :=
  [(fun j => if j=i then 1 else 0,1)]

theorem toExponent_injective : Function.Injective (toExponent : Exponent r → Fin r →₀ ℕ) := by
  intro e f equal
  funext i
  exact congrArg (fun v : Fin r →₀ ℕ => v i) equal

theorem meaning_coefficient (e : Exponent r) (p : Sparse r K) :
    coeff (toExponent e) (meaning p)=coefficient e p := by
  classical
  induction p with
  | nil => simp [meaning,coefficient]
  | cons term p ih =>
      obtain ⟨f,c⟩ := term
      simp [meaning,coefficient,coeff_monomial,toExponent_injective.eq_iff,ih]

theorem meaning_sum (ps : List (Sparse r K)) : meaning (sum ps)=(ps.map meaning).sum := by
  induction ps with
  | nil => simp [sum,meaning]
  | cons p ps ih => simp [sum,meaning_add,ih]

theorem meaning_variable (i : Fin r) : meaning (variableData i : Sparse r K)=X i := by
  have exponent : toExponent (fun j : Fin r => if j=i then 1 else 0)=Finsupp.single i 1 := by
    ext j
    simp [toExponent,Finsupp.single_apply,eq_comm]
  simp [variableData,meaning,exponent,X]

end ICFallDegrees.SparsePolynomial
