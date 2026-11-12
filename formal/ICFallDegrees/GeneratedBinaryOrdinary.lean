import ICFallDegrees.SparseLowerCertificate
import ICFallDegrees.SparseDAG
import Mathlib.Algebra.Field.ZMod
import Mathlib.Tactic.FinCases

/-! Generated candidate data. Each identity is verified by kernel reduction. -/
namespace ICFallDegrees.GeneratedBinaryOrdinary
open SparsePolynomial SparseLowerCertificate SparseDAG
abbrev K := ZMod 2
def exponent (a b c d : ℕ) : Exponent 4 := fun i => [a,b,c,d].get i
def inputs : List (Sparse 4 K) := [
  [((exponent 2 0 0 0),1), ((exponent 1 0 0 0),1)],
  [((exponent 0 2 0 0),1), ((exponent 0 1 0 0),1)],
  [((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)],
  [((exponent 0 0 0 2),1), ((exponent 0 0 0 1),1)],
  [((exponent 1 0 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 1),1)],
  [((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)],
  [((exponent 1 0 1 0),1), ((exponent 1 0 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 0 1 0),1)]
]
def basisData : List (Sparse 4 K) := [
  [((exponent 2 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)],
  [((exponent 1 1 0 0),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 0),1)],
  [((exponent 1 0 1 0),1)],
  [((exponent 1 0 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 1),1)],
  [((exponent 0 2 0 0),1), ((exponent 0 1 0 0),1)],
  [((exponent 0 1 1 0),1), ((exponent 0 0 1 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 0),1)],
  [((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1)],
  [((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)],
  [((exponent 0 0 0 2),1), ((exponent 0 0 0 1),1)],
  [((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)]
]
def pivotData : List (Exponent 4) := [
  (exponent 2 0 0 0),
  (exponent 1 1 0 0),
  (exponent 1 0 1 0),
  (exponent 1 0 0 1),
  (exponent 0 2 0 0),
  (exponent 0 1 1 0),
  (exponent 0 1 0 1),
  (exponent 0 0 2 0),
  (exponent 0 0 0 2),
  (exponent 1 0 0 0)
]
def basis (i : Fin 10) : Sparse 4 K := basisData.get i
def pivot (i : Fin 10) : Exponent 4 := pivotData.get i
set_option maxRecDepth 20000
set_option maxHeartbeats 100000000
theorem lower_checked : SparseLowerCertificate.check inputs basis pivot 2=true := by rfl
theorem one_not_at_two : ¬InClosure {p | p∈inputs.map meaning} 2 (1 : MvPolynomial (Fin 4) K) :=
  check_excludes_one inputs basis pivot 2 lower_checked
def claim0 : Sparse 4 K := [((exponent 2 0 0 0),1), ((exponent 1 0 0 0),1)]
def claim1 : Sparse 4 K := [((exponent 0 2 0 0),1), ((exponent 0 1 0 0),1)]
def claim2 : Sparse 4 K := [((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)]
def claim3 : Sparse 4 K := [((exponent 0 0 0 2),1), ((exponent 0 0 0 1),1)]
def claim4 : Sparse 4 K := [((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 1 0 0 0),1), ((exponent 0 1 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 1),1)]
def claim5 : Sparse 4 K := [((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)]
def claim6 : Sparse 4 K := [((exponent 1 0 1 0),1), ((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1)]
def claim7 : Sparse 4 K := [((exponent 1 0 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1)]
def claim8 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 1 0 1 0),1)]
def claim9 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1)]
def claim10 : Sparse 4 K := [((exponent 0 2 1 0),1), ((exponent 0 1 1 0),1)]
def claim11 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 1 0 1 0),1)]
def claim12 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1)]
def claim13 : Sparse 4 K := [((exponent 1 0 0 2),1), ((exponent 1 0 0 1),1)]
def claim14 : Sparse 4 K := [((exponent 1 0 0 2),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)]
def claim15 : Sparse 4 K := [((exponent 1 1 0 1),1), ((exponent 0 2 1 0),1), ((exponent 1 1 0 0),1), ((exponent 0 2 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1)]
def claim16 : Sparse 4 K := [((exponent 1 1 0 1),1), ((exponent 1 1 0 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1)]
def claim17 : Sparse 4 K := [((exponent 1 0 0 2),1), ((exponent 0 1 1 1),1), ((exponent 1 0 0 1),1), ((exponent 0 1 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 2),1)]
def claim18 : Sparse 4 K := [((exponent 0 1 1 1),1), ((exponent 0 1 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 1),1)]
def claim19 : Sparse 4 K := [((exponent 1 0 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 1),1)]
def claim20 : Sparse 4 K := [((exponent 0 1 1 0),1), ((exponent 0 0 1 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 0),1)]
def claim21 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 1 1 0 1),1), ((exponent 1 1 0 0),1), ((exponent 1 0 0 1),1)]
def claim22 : Sparse 4 K := [((exponent 0 0 1 1),1)]
def claim23 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 0 1 1 1),1), ((exponent 0 1 1 0),1), ((exponent 0 0 1 1),1)]
def claim24 : Sparse 4 K := [((exponent 0 0 0 0),1)]
def program : List (Entry 4 K) := [
  (.input 0,claim0),
  (.input 1,claim1),
  (.input 2,claim2),
  (.input 3,claim3),
  (.input 4,claim4),
  (.input 5,claim5),
  (.input 6,claim6),
  (.linear [(4,1), (6,1)],claim7),
  (.mul (variableData 2) 0,claim8),
  (.linear [(7,1), (8,1)],claim9),
  (.mul (variableData 2) 1,claim10),
  (.mul (variableData 0) 2,claim11),
  (.linear [(7,1), (11,1)],claim12),
  (.mul (variableData 0) 3,claim13),
  (.linear [(4,1), (5,1), (13,1)],claim14),
  (.mul (variableData 1) 4,claim15),
  (.linear [(1,1), (10,1), (15,1)],claim16),
  (.mul (variableData 3) 4,claim17),
  (.linear [(3,1), (4,1), (5,1), (14,1), (17,1)],claim18),
  (.mul (variableData 3) 5,claim19),
  (.linear [(4,1), (5,1), (19,1)],claim20),
  (.mul (variableData 0) 7,claim21),
  (.linear [(4,1), (5,1), (9,1), (16,1), (20,1), (21,1)],claim22),
  (.mul (variableData 2) 7,claim23),
  (.linear [(12,1), (18,1), (20,1), (22,1), (23,1)],claim24)
]
def prefix0 : List (Sparse 4 K) := []
def prefix1 : List (Sparse 4 K) := prefix0++[claim0]
def prefix2 : List (Sparse 4 K) := prefix1++[claim1]
def prefix3 : List (Sparse 4 K) := prefix2++[claim2]
def prefix4 : List (Sparse 4 K) := prefix3++[claim3]
def prefix5 : List (Sparse 4 K) := prefix4++[claim4]
def prefix6 : List (Sparse 4 K) := prefix5++[claim5]
def prefix7 : List (Sparse 4 K) := prefix6++[claim6]
def prefix8 : List (Sparse 4 K) := prefix7++[claim7]
def prefix9 : List (Sparse 4 K) := prefix8++[claim8]
def prefix10 : List (Sparse 4 K) := prefix9++[claim9]
def prefix11 : List (Sparse 4 K) := prefix10++[claim10]
def prefix12 : List (Sparse 4 K) := prefix11++[claim11]
def prefix13 : List (Sparse 4 K) := prefix12++[claim12]
def prefix14 : List (Sparse 4 K) := prefix13++[claim13]
def prefix15 : List (Sparse 4 K) := prefix14++[claim14]
def prefix16 : List (Sparse 4 K) := prefix15++[claim15]
def prefix17 : List (Sparse 4 K) := prefix16++[claim16]
def prefix18 : List (Sparse 4 K) := prefix17++[claim17]
def prefix19 : List (Sparse 4 K) := prefix18++[claim18]
def prefix20 : List (Sparse 4 K) := prefix19++[claim19]
def prefix21 : List (Sparse 4 K) := prefix20++[claim20]
def prefix22 : List (Sparse 4 K) := prefix21++[claim21]
def prefix23 : List (Sparse 4 K) := prefix22++[claim22]
def prefix24 : List (Sparse 4 K) := prefix23++[claim23]
def prefix25 : List (Sparse 4 K) := prefix24++[claim24]
theorem step_checked0 : checkEntry inputs 3 prefix0 (.input 0,claim0)=true := by rfl
theorem step_checked1 : checkEntry inputs 3 prefix1 (.input 1,claim1)=true := by rfl
theorem step_checked2 : checkEntry inputs 3 prefix2 (.input 2,claim2)=true := by rfl
theorem step_checked3 : checkEntry inputs 3 prefix3 (.input 3,claim3)=true := by rfl
theorem step_checked4 : checkEntry inputs 3 prefix4 (.input 4,claim4)=true := by rfl
theorem step_checked5 : checkEntry inputs 3 prefix5 (.input 5,claim5)=true := by rfl
theorem step_checked6 : checkEntry inputs 3 prefix6 (.input 6,claim6)=true := by rfl
theorem step_checked7 : checkEntry inputs 3 prefix7 (.linear [(4,1), (6,1)],claim7)=true := by rfl
theorem step_checked8 : checkEntry inputs 3 prefix8 (.mul (variableData 2) 0,claim8)=true := by rfl
theorem step_checked9 : checkEntry inputs 3 prefix9 (.linear [(7,1), (8,1)],claim9)=true := by rfl
theorem step_checked10 : checkEntry inputs 3 prefix10 (.mul (variableData 2) 1,claim10)=true := by rfl
theorem step_checked11 : checkEntry inputs 3 prefix11 (.mul (variableData 0) 2,claim11)=true := by rfl
theorem step_checked12 : checkEntry inputs 3 prefix12 (.linear [(7,1), (11,1)],claim12)=true := by rfl
theorem step_checked13 : checkEntry inputs 3 prefix13 (.mul (variableData 0) 3,claim13)=true := by rfl
theorem step_checked14 : checkEntry inputs 3 prefix14 (.linear [(4,1), (5,1), (13,1)],claim14)=true := by rfl
theorem step_checked15 : checkEntry inputs 3 prefix15 (.mul (variableData 1) 4,claim15)=true := by rfl
theorem step_checked16 : checkEntry inputs 3 prefix16 (.linear [(1,1), (10,1), (15,1)],claim16)=true := by rfl
theorem step_checked17 : checkEntry inputs 3 prefix17 (.mul (variableData 3) 4,claim17)=true := by rfl
theorem step_checked18 : checkEntry inputs 3 prefix18 (.linear [(3,1), (4,1), (5,1), (14,1), (17,1)],claim18)=true := by rfl
theorem step_checked19 : checkEntry inputs 3 prefix19 (.mul (variableData 3) 5,claim19)=true := by rfl
theorem step_checked20 : checkEntry inputs 3 prefix20 (.linear [(4,1), (5,1), (19,1)],claim20)=true := by rfl
theorem step_checked21 : checkEntry inputs 3 prefix21 (.mul (variableData 0) 7,claim21)=true := by rfl
theorem step_checked22 : checkEntry inputs 3 prefix22 (.linear [(4,1), (5,1), (9,1), (16,1), (20,1), (21,1)],claim22)=true := by rfl
theorem step_checked23 : checkEntry inputs 3 prefix23 (.mul (variableData 2) 7,claim23)=true := by rfl
theorem step_checked24 : checkEntry inputs 3 prefix24 (.linear [(12,1), (18,1), (20,1), (22,1), (23,1)],claim24)=true := by rfl
theorem all_steps (i : Fin 25) : checkEntry inputs 3 ((program.take i.val).map Prod.snd) (program.get i)=true := by
  fin_cases i
  · exact step_checked0
  · exact step_checked1
  · exact step_checked2
  · exact step_checked3
  · exact step_checked4
  · exact step_checked5
  · exact step_checked6
  · exact step_checked7
  · exact step_checked8
  · exact step_checked9
  · exact step_checked10
  · exact step_checked11
  · exact step_checked12
  · exact step_checked13
  · exact step_checked14
  · exact step_checked15
  · exact step_checked16
  · exact step_checked17
  · exact step_checked18
  · exact step_checked19
  · exact step_checked20
  · exact step_checked21
  · exact step_checked22
  · exact step_checked23
  · exact step_checked24
theorem upper_checked : checkClaimedUnit inputs 3 program=true := by
  have all := checkEntriesFrom_of_each inputs 3 program [] all_steps
  change (checkEntriesFrom inputs 3 [] program && (add claim24 (scale (-1) SparseCertificate.oneData)).isEmpty)=true
  rw [all]
  rfl
theorem one_at_three : InClosure {p | p∈inputs.map meaning} 3 (1 : MvPolynomial (Fin 4) K) :=
  checkClaimedUnit_sound inputs 3 program upper_checked
theorem exact_lfd_three : IsLastFallDegree {p | p∈inputs.map meaning} 3 := by
  apply exactLFD_of_first_one (by omega) one_at_three
  exact one_not_at_two
end ICFallDegrees.GeneratedBinaryOrdinary
