import ICFallDegrees.SparseCertificate
import ICFallDegrees.PrimeCounterexample

/-!
# Kernel-replayed sparse prime-field certificates

Executable sparse data encode the actual F_5 example. The checker derives 1
at normalized budget 3 and raw budget 4; lower bounds are the independently
proved containing-space/evaluation arguments. Acceptance is proved by
kernel reduction, without native_decide, external rank assumptions, or an
unverified arithmetic implementation override.
-/

namespace ICFallDegrees.SparsePrimeCertificate

open MvPolynomial SparsePolynomial SparseCertificate

abbrev K := PrimeCounterexample.K
abbrev P := QuadraticFamily.P K

def exponent (a b : ℕ) : Exponent 2 := fun i => if i=0 then a else b
def term (a b : ℕ) (c : K) : Sparse 2 K := [(exponent a b,c)]
def cxData : Sparse 2 K := [(exponent 2 0,1),(exponent 0 0,4)]
def cyData : Sparse 2 K := [(exponent 0 2,1),(exponent 0 0,4)]
def hData : Sparse 2 K := term 1 1 1
def rawData : Sparse 2 K := [(exponent 2 2,1),(exponent 1 1,1),(exponent 0 0,4)]
def normalizedInputs : List (Sparse 2 K) := [cxData,cyData,hData]
def rawInputs : List (Sparse 2 K) := [cxData,cyData,rawData]

def xCertificate (hc : Proof 2 K) : Proof 2 K :=
  .add (.mul (term 0 1 1) hc) (.scale 4 (.mul (term 1 0 1) (.input 1)))
def unitCertificate (hc : Proof 2 K) : Proof 2 K :=
  .add (.mul (term 1 0 1) (xCertificate hc)) (.scale 4 (.input 0))
def normalizedUnit : Proof 2 K := unitCertificate (.input 2)
def normalizedFromRaw : Proof 2 K :=
  .add (.add (.input 2) (.scale 4 (.mul (term 0 2 1) (.input 0)))) (.scale 4 (.input 1))
def rawUnit : Proof 2 K := unitCertificate normalizedFromRaw

theorem normalized_accepted : check normalizedInputs 3 normalizedUnit=some oneData := by
  have computed : check normalizedInputs 3 normalizedUnit=some [(exponent 0 0,1)] := by rfl
  have zero : exponent 0 0=(fun _ => 0) := by funext i; simp [exponent]
  rw [zero] at computed
  exact computed
theorem raw_accepted : check rawInputs 4 rawUnit=some oneData := by
  have computed : check rawInputs 4 rawUnit=some [(exponent 0 0,1)] := by rfl
  have zero : exponent 0 0=(fun _ => 0) := by funext i; simp [exponent]
  rw [zero] at computed
  exact computed
theorem normalized_lower_budget_rejected : check normalizedInputs 2 normalizedUnit=none := by decide
theorem raw_lower_budget_rejected : check rawInputs 3 rawUnit=none := by decide
theorem invalid_input_rejected : check normalizedInputs 3 (.input 99)=none := by decide
theorem cancellation_lowers_computed_degree :
    degree (normalize ([(exponent 4 4,1),(exponent 4 4,4),(exponent 0 0,1)] : Sparse 2 K))=0 := by
  decide
theorem raw_degree_is_charged : degree (normalize rawData)=4 := by decide

theorem toExponent_exponent (a b : ℕ) :
    toExponent (exponent a b)=Finsupp.single (0 : Fin 2) a+Finsupp.single 1 b := by
  ext i
  have hi : i=0 ∨ i=1 := by omega
  rcases hi with rfl | rfl <;> simp [toExponent,exponent]

theorem meaning_cx : meaning cxData=(QuadraticFamily.cx : P) := by
  have four : (4 : K)= -1 := by decide
  simp [cxData,meaning,toExponent_exponent,QuadraticFamily.cx,QuadraticFamily.x,
    X_pow_eq_monomial,four,sub_eq_add_neg]

theorem meaning_cy : meaning cyData=(QuadraticFamily.cy : P) := by
  have four : (4 : K)= -1 := by decide
  simp [cyData,meaning,toExponent_exponent,QuadraticFamily.cy,QuadraticFamily.y,
    X_pow_eq_monomial,four,sub_eq_add_neg]

theorem meaning_h : meaning hData=QuadraticFamily.h (0 : K) := by
  simp [hData,term,meaning,toExponent_exponent,QuadraticFamily.h,
    QuadraticFamily.x,QuadraticFamily.y,X,monomial_mul]

theorem meaning_raw : meaning rawData=QuadraticFamily.raw (0 : K) := by
  have four : (4 : K)= -1 := by decide
  rw [QuadraticFamily.raw_form]
  simp [rawData,meaning,toExponent_exponent,QuadraticFamily.h,
    QuadraticFamily.x,QuadraticFamily.y,X_pow_eq_monomial,X,monomial_mul,
    monomial_pow,four,sub_eq_add_neg,add_assoc]

theorem meaning_raw_is_curveS3 : meaning rawData=PrimeCounterexample.curveS3 :=
  meaning_raw.trans PrimeCounterexample.raw_is_curveS3

theorem normalized_input_semantics :
    {p : P | p∈normalizedInputs.map meaning}=QuadraticFamily.inputs (0 : K) := by
  ext p
  simp [normalizedInputs,meaning_cx,meaning_cy,meaning_h,QuadraticFamily.inputs]

theorem raw_input_semantics :
    {p : P | p∈rawInputs.map meaning}=QuadraticFamily.rawInputs (0 : K) := by
  ext p
  simp [rawInputs,meaning_cx,meaning_cy,meaning_raw,QuadraticFamily.rawInputs]

theorem normalized_one_from_checker : InClosure (QuadraticFamily.inputs (0 : K)) 3 1 := by
  rw [←normalized_input_semantics]
  simpa only [meaning_oneData] using check_sound normalizedInputs 3 normalizedUnit normalized_accepted

theorem raw_one_from_checker : InClosure (QuadraticFamily.rawInputs (0 : K)) 4 1 := by
  rw [←raw_input_semantics]
  simpa only [meaning_oneData] using check_sound rawInputs 4 rawUnit raw_accepted

theorem normalized_exact_lfd_from_checker : IsLastFallDegree (QuadraticFamily.inputs (0 : K)) 3 := by
  apply exactLFD_of_first_one (by omega) normalized_one_from_checker
  exact QuadraticFamily.no_one_at_two 0

theorem raw_exact_lfd_from_checker : IsLastFallDegree (QuadraticFamily.rawInputs (0 : K)) 4 := by
  apply exactLFD_of_first_one (by omega) raw_one_from_checker
  exact QuadraticFamily.raw_no_one_at_three 0

end ICFallDegrees.SparsePrimeCertificate
