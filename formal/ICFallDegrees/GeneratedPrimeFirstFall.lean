import ICFallDegrees.GeneratedPrimeLower
import ICFallDegrees.Prime257System
import ICFallDegrees.LowerCertificateSupport
import ICFallDegrees.SparseCertificate
import ICFallDegrees.FirstFall
import Mathlib.Tactic.FinCases

namespace ICFallDegrees.GeneratedPrimeFirstFall
open SparsePolynomial GeneratedPrimeLower LowerCertificateSupport

def basisData : List (Sparse 6 K) := [
  [((exponent 2 0 0 0 0 0),1), ((exponent 0 1 0 0 0 0),256)],
  [((exponent 1 0 0 1 0 0),1), ((exponent 1 0 0 0 1 0),52), ((exponent 0 1 0 1 0 0),52), ((exponent 0 1 0 0 1 0),134), ((exponent 1 0 0 0 0 0),120), ((exponent 0 1 0 0 0 0),193), ((exponent 0 0 0 1 0 0),120), ((exponent 0 0 0 0 1 0),193), ((exponent 0 0 0 0 0 0),195)],
  [((exponent 0 2 0 0 0 0),1), ((exponent 0 0 1 0 0 0),256)],
  [((exponent 0 0 2 0 0 0),1), ((exponent 0 0 0 0 0 0),256)],
  [((exponent 0 0 0 2 0 0),1), ((exponent 0 0 0 0 1 0),256)],
  [((exponent 0 0 0 0 2 0),1), ((exponent 0 0 0 0 0 1),256)],
  [((exponent 0 0 0 0 0 2),1), ((exponent 0 0 0 0 0 0),256)]
]
def pivotData : List (Exponent 6) := [
  (exponent 2 0 0 0 0 0),
  (exponent 1 0 0 1 0 0),
  (exponent 0 2 0 0 0 0),
  (exponent 0 0 2 0 0 0),
  (exponent 0 0 0 2 0 0),
  (exponent 0 0 0 0 2 0),
  (exponent 0 0 0 0 0 2)
]
def basis (i : Fin 7) : Sparse 6 K := basisData.get i
def pivot (i : Fin 7) : Exponent 6 := pivotData.get i
def linA : Sparse 6 K := [(exponent 0 1 0 0 0 0,1),(exponent 1 0 0 0 0 0,84)]
def linB : Sparse 6 K := [(exponent 0 0 0 0 1 0,1),(exponent 0 0 0 1 0 0,173)]
def certificate : SparseCertificate.Proof 6 K :=
  .add (.mul linA (.input 6))
    (.scale 256 (.mul linB (.add (.input 1) (.scale 140 (.input 0)))))
def witnessData : Sparse 6 K := (SparseCertificate.check inputs 3 certificate).getD []
def separatingExponent : Exponent 6 := exponent 0 0 1 0 1 0

set_option maxRecDepth 20000
set_option maxHeartbeats 100000000
theorem checked_space : SparseLowerCertificate.check inputs basis pivot 2=true := by rfl
theorem checked_fall : SparseCertificate.check inputs 3 certificate=some witnessData := by rfl
theorem witness_degree : degree witnessData=2 := by rfl
theorem witness_coefficient : coefficient separatingExponent witnessData=1 := by rfl
theorem basis_coefficient_zero (i : Fin 7) : coefficient separatingExponent (basis i)=0 := by
  fin_cases i <;> rfl
theorem pivot_degree (i : Fin 7) : SparseLowerCertificate.weight (pivot i)=2 := by
  fin_cases i <;> rfl

theorem witness_not_at_two : ¬InClosure {p | p∈inputs.map meaning} 2 (meaning witnessData) := by
  intro derived
  have member := checked_closure_containment inputs basis pivot 2 checked_space derived
  have zero := rowSpace_coefficient_zero (fun i => meaning (basis i))
    (toExponent separatingExponent) (fun i => by rw [meaning_coefficient]; exact basis_coefficient_zero i) member
  rw [meaning_coefficient,witness_coefficient] at zero
  exact one_ne_zero zero

theorem affine_first_strict_fall_three : IsFirstStrictFall {p | p∈inputs.map meaning} 3 := by
  have upper := SparseCertificate.check_sound inputs 3 certificate checked_fall
  have degree : (meaning witnessData).totalDegree<3 := by
    have bound := meaning_degree_le witnessData
    rw [witness_degree] at bound
    omega
  refine ⟨⟨meaning witnessData,upper,degree,witness_not_at_two⟩,?_⟩
  intro d hd
  rintro ⟨p,hp,dp,absent⟩
  have member := checked_closure_containment inputs basis pivot 2 checked_space (hp.mono (by omega))
  have pivots : ∀ i : Fin 7, 2≤(toExponent (pivot i)).degree := by
    intro i
    rw [SparseLowerCertificate.weight_toExponent,pivot_degree]
  have zero := rowSpace_no_low (fun i => meaning (basis i)) (fun i => toExponent (pivot i))
    (checked_rref inputs basis pivot 2 checked_space) 2 pivots member (by omega)
  subst p
  exact absent InClosure.zero

theorem actual_curve_first_strict_fall_three : IsFirstStrictFall Prime257System.normalizedInputs 3 := by
  rw [←Prime257System.input_set_semantics]
  exact affine_first_strict_fall_three

end ICFallDegrees.GeneratedPrimeFirstFall
