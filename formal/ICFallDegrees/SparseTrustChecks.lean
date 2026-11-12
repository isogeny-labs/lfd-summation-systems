import ICFallDegrees.SparseDAG
import ICFallDegrees.SparseLowerCertificate
import Mathlib.Data.ZMod.Basic

namespace ICFallDegrees.SparseTrustChecks
open SparsePolynomial
abbrev K := ZMod 5
def x : Sparse 1 K := [(fun _ => 1,1)]
def one : Sparse 1 K := SparseCertificate.oneData
def basis (_ : Fin 1) : Sparse 1 K := x
def pivot (_ : Fin 1) : Exponent 1 := fun _ => 1

/-- Small acceptance and mutation checks for the two executable trust paths. -/
theorem lower_accepted : SparseLowerCertificate.check [x] basis pivot 1=true := by rfl
theorem zero_pivot_rejected :
    SparseLowerCertificate.check [x] basis (fun _ _ => 0) 1=false := by rfl
theorem wrong_coefficient_rejected :
    SparseLowerCertificate.check [x] (fun _ => scale 2 x) pivot 1=false := by rfl
theorem missing_input_span_rejected :
    SparseLowerCertificate.check [one] basis pivot 1=false := by rfl
theorem dag_accepted : SparseDAG.checkClaimedUnit [one] 0 [(.input 0,one)]=true := by rfl
theorem tampered_row_rejected : SparseDAG.checkClaimedUnit [one] 0 [(.input 0,x)]=false := by rfl
theorem forward_reference_rejected :
    SparseDAG.checkClaimedUnit [one] 0 [(.linear [(0,1)],one)]=false := by rfl
theorem excessive_budget_rejected :
    SparseDAG.checkEntry [x] 1 [x] (.mul x 0,x)=false := by rfl
end ICFallDegrees.SparseTrustChecks
