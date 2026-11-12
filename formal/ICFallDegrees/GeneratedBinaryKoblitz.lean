import ICFallDegrees.SparseLowerCertificate
import ICFallDegrees.SparseDAG
import Mathlib.Algebra.Field.ZMod
import Mathlib.Tactic.FinCases

/-! Generated candidate data. Each identity is verified by kernel reduction. -/
namespace ICFallDegrees.GeneratedBinaryKoblitz
open SparsePolynomial SparseLowerCertificate SparseDAG
abbrev K := ZMod 2
def exponent (a b c d : ℕ) : Exponent 4 := fun i => [a,b,c,d].get i
def inputs : List (Sparse 4 K) := [
  [((exponent 2 0 0 0),1), ((exponent 1 0 0 0),1)],
  [((exponent 0 2 0 0),1), ((exponent 0 1 0 0),1)],
  [((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)],
  [((exponent 0 0 0 2),1), ((exponent 0 0 0 1),1)],
  [((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)],
  [((exponent 1 0 1 0),1), ((exponent 1 0 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 1),1)],
  [((exponent 1 0 1 0),1), ((exponent 1 0 0 0),1), ((exponent 0 1 0 1),1), ((exponent 0 0 1 0),1)],
  [((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1)]
]
def basisData : List (Sparse 4 K) := [
  [((exponent 2 0 0 0),1), ((exponent 1 0 0 0),1)],
  [((exponent 1 1 0 0),1), ((exponent 0 0 1 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)],
  [((exponent 1 0 1 0),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1)],
  [((exponent 1 0 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)],
  [((exponent 0 2 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)],
  [((exponent 0 1 1 0),1), ((exponent 0 0 1 1),1), ((exponent 0 0 1 0),1)],
  [((exponent 0 1 0 1),1)],
  [((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)],
  [((exponent 0 0 0 2),1), ((exponent 0 0 0 1),1)],
  [((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)]
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
  (exponent 0 1 0 0)
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
def claim4 : Sparse 4 K := [((exponent 0 1 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)]
def claim5 : Sparse 4 K := [((exponent 1 0 1 0),1), ((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 1 0 0 0),1), ((exponent 0 1 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 1),1)]
def claim6 : Sparse 4 K := [((exponent 1 0 1 0),1), ((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)]
def claim7 : Sparse 4 K := [((exponent 1 0 1 0),1), ((exponent 0 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1)]
def claim8 : Sparse 4 K := [((exponent 1 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 0 0 0),1)]
def claim9 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 1 0 1 0),1)]
def claim10 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 0 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1)]
def claim11 : Sparse 4 K := [((exponent 2 0 0 1),1), ((exponent 1 0 0 1),1)]
def claim12 : Sparse 4 K := [((exponent 2 0 0 1),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 0 0 0),1)]
def claim13 : Sparse 4 K := [((exponent 0 2 1 0),1), ((exponent 0 1 1 0),1)]
def claim14 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 1 0 1 0),1)]
def claim15 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 0 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 1 0),1)]
def claim16 : Sparse 4 K := [((exponent 0 1 2 0),1), ((exponent 0 1 1 0),1)]
def claim17 : Sparse 4 K := [((exponent 1 0 0 2),1), ((exponent 1 0 0 1),1)]
def claim18 : Sparse 4 K := [((exponent 1 0 0 2),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 0 0 0 0),1)]
def claim19 : Sparse 4 K := [((exponent 1 1 0 0),1), ((exponent 1 0 0 1),1), ((exponent 1 0 0 0),1)]
def claim20 : Sparse 4 K := [((exponent 1 1 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 0 0),1)]
def claim21 : Sparse 4 K := [((exponent 0 2 0 0),1), ((exponent 0 1 0 1),1), ((exponent 0 1 0 0),1)]
def claim22 : Sparse 4 K := [((exponent 0 1 0 1),1)]
def claim23 : Sparse 4 K := [((exponent 0 1 1 0),1), ((exponent 0 0 1 1),1), ((exponent 0 0 1 0),1)]
def claim24 : Sparse 4 K := [((exponent 2 0 1 0),1), ((exponent 2 0 0 1),1), ((exponent 1 1 1 0),1), ((exponent 2 0 0 0),1), ((exponent 1 0 1 0),1), ((exponent 1 0 0 0),1)]
def claim25 : Sparse 4 K := [((exponent 1 1 1 0),1), ((exponent 0 0 1 1),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 0),1)]
def claim26 : Sparse 4 K := [((exponent 1 1 1 0),1), ((exponent 1 1 0 1),1), ((exponent 0 2 1 0),1), ((exponent 1 1 0 0),1), ((exponent 0 1 1 0),1), ((exponent 0 1 0 0),1)]
def claim27 : Sparse 4 K := [((exponent 1 1 0 1),1), ((exponent 1 0 0 0),1), ((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)]
def claim28 : Sparse 4 K := [((exponent 1 0 2 0),1), ((exponent 1 0 1 1),1), ((exponent 0 1 2 0),1), ((exponent 1 0 1 0),1), ((exponent 0 0 2 0),1), ((exponent 0 0 1 0),1)]
def claim29 : Sparse 4 K := [((exponent 1 0 1 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 1 0),1)]
def claim30 : Sparse 4 K := [((exponent 1 0 1 1),1), ((exponent 1 0 0 2),1), ((exponent 0 1 1 1),1), ((exponent 1 0 0 1),1), ((exponent 0 0 1 1),1), ((exponent 0 0 0 1),1)]
def claim31 : Sparse 4 K := [((exponent 0 1 1 1),1), ((exponent 0 0 1 0),1), ((exponent 0 0 0 1),1)]
def claim32 : Sparse 4 K := [((exponent 2 0 0 1),1), ((exponent 1 1 1 0),1), ((exponent 1 1 0 1),1), ((exponent 1 0 0 0),1)]
def claim33 : Sparse 4 K := [((exponent 0 0 0 1),1), ((exponent 0 0 0 0),1)]
def claim34 : Sparse 4 K := [((exponent 1 0 1 1),1), ((exponent 0 1 2 0),1), ((exponent 0 1 1 1),1), ((exponent 0 0 1 0),1)]
def claim35 : Sparse 4 K := [((exponent 0 0 0 0),1)]
def program : List (Entry 4 K) := [
  (.input 0,claim0),
  (.input 1,claim1),
  (.input 2,claim2),
  (.input 3,claim3),
  (.input 4,claim4),
  (.input 5,claim5),
  (.linear [(4,1), (5,1)],claim6),
  (.input 6,claim7),
  (.linear [(6,1), (7,1)],claim8),
  (.mul (variableData 2) 0,claim9),
  (.linear [(6,1), (8,1), (9,1)],claim10),
  (.mul (variableData 3) 0,claim11),
  (.linear [(8,1), (11,1)],claim12),
  (.mul (variableData 2) 1,claim13),
  (.mul (variableData 0) 2,claim14),
  (.linear [(6,1), (8,1), (14,1)],claim15),
  (.mul (variableData 1) 2,claim16),
  (.mul (variableData 0) 3,claim17),
  (.linear [(8,1), (17,1)],claim18),
  (.mul (variableData 0) 4,claim19),
  (.linear [(8,1), (19,1)],claim20),
  (.mul (variableData 1) 4,claim21),
  (.linear [(1,1), (21,1)],claim22),
  (.mul (variableData 2) 4,claim23),
  (.mul (variableData 0) 6,claim24),
  (.linear [(0,1), (6,1), (8,1), (10,1), (12,1), (22,1), (23,1), (24,1)],claim25),
  (.mul (variableData 1) 6,claim26),
  (.linear [(4,1), (13,1), (20,1), (22,1), (23,1), (25,1), (26,1)],claim27),
  (.mul (variableData 2) 6,claim28),
  (.linear [(2,1), (6,1), (8,1), (15,1), (16,1), (23,1), (28,1)],claim29),
  (.mul (variableData 3) 6,claim30),
  (.linear [(8,1), (18,1), (29,1), (30,1)],claim31),
  (.mul (variableData 0) 8,claim32),
  (.linear [(12,1), (22,1), (23,1), (25,1), (27,1), (32,1)],claim33),
  (.mul (variableData 2) 8,claim34),
  (.linear [(16,1), (23,1), (29,1), (31,1), (33,1), (34,1)],claim35)
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
def prefix26 : List (Sparse 4 K) := prefix25++[claim25]
def prefix27 : List (Sparse 4 K) := prefix26++[claim26]
def prefix28 : List (Sparse 4 K) := prefix27++[claim27]
def prefix29 : List (Sparse 4 K) := prefix28++[claim28]
def prefix30 : List (Sparse 4 K) := prefix29++[claim29]
def prefix31 : List (Sparse 4 K) := prefix30++[claim30]
def prefix32 : List (Sparse 4 K) := prefix31++[claim31]
def prefix33 : List (Sparse 4 K) := prefix32++[claim32]
def prefix34 : List (Sparse 4 K) := prefix33++[claim33]
def prefix35 : List (Sparse 4 K) := prefix34++[claim34]
def prefix36 : List (Sparse 4 K) := prefix35++[claim35]
theorem step_checked0 : checkEntry inputs 3 prefix0 (.input 0,claim0)=true := by rfl
theorem step_checked1 : checkEntry inputs 3 prefix1 (.input 1,claim1)=true := by rfl
theorem step_checked2 : checkEntry inputs 3 prefix2 (.input 2,claim2)=true := by rfl
theorem step_checked3 : checkEntry inputs 3 prefix3 (.input 3,claim3)=true := by rfl
theorem step_checked4 : checkEntry inputs 3 prefix4 (.input 4,claim4)=true := by rfl
theorem step_checked5 : checkEntry inputs 3 prefix5 (.input 5,claim5)=true := by rfl
theorem step_checked6 : checkEntry inputs 3 prefix6 (.linear [(4,1), (5,1)],claim6)=true := by rfl
theorem step_checked7 : checkEntry inputs 3 prefix7 (.input 6,claim7)=true := by rfl
theorem step_checked8 : checkEntry inputs 3 prefix8 (.linear [(6,1), (7,1)],claim8)=true := by rfl
theorem step_checked9 : checkEntry inputs 3 prefix9 (.mul (variableData 2) 0,claim9)=true := by rfl
theorem step_checked10 : checkEntry inputs 3 prefix10 (.linear [(6,1), (8,1), (9,1)],claim10)=true := by rfl
theorem step_checked11 : checkEntry inputs 3 prefix11 (.mul (variableData 3) 0,claim11)=true := by rfl
theorem step_checked12 : checkEntry inputs 3 prefix12 (.linear [(8,1), (11,1)],claim12)=true := by rfl
theorem step_checked13 : checkEntry inputs 3 prefix13 (.mul (variableData 2) 1,claim13)=true := by rfl
theorem step_checked14 : checkEntry inputs 3 prefix14 (.mul (variableData 0) 2,claim14)=true := by rfl
theorem step_checked15 : checkEntry inputs 3 prefix15 (.linear [(6,1), (8,1), (14,1)],claim15)=true := by rfl
theorem step_checked16 : checkEntry inputs 3 prefix16 (.mul (variableData 1) 2,claim16)=true := by rfl
theorem step_checked17 : checkEntry inputs 3 prefix17 (.mul (variableData 0) 3,claim17)=true := by rfl
theorem step_checked18 : checkEntry inputs 3 prefix18 (.linear [(8,1), (17,1)],claim18)=true := by rfl
theorem step_checked19 : checkEntry inputs 3 prefix19 (.mul (variableData 0) 4,claim19)=true := by rfl
theorem step_checked20 : checkEntry inputs 3 prefix20 (.linear [(8,1), (19,1)],claim20)=true := by rfl
theorem step_checked21 : checkEntry inputs 3 prefix21 (.mul (variableData 1) 4,claim21)=true := by rfl
theorem step_checked22 : checkEntry inputs 3 prefix22 (.linear [(1,1), (21,1)],claim22)=true := by rfl
theorem step_checked23 : checkEntry inputs 3 prefix23 (.mul (variableData 2) 4,claim23)=true := by rfl
theorem step_checked24 : checkEntry inputs 3 prefix24 (.mul (variableData 0) 6,claim24)=true := by rfl
theorem step_checked25 : checkEntry inputs 3 prefix25 (.linear [(0,1), (6,1), (8,1), (10,1), (12,1), (22,1), (23,1), (24,1)],claim25)=true := by rfl
theorem step_checked26 : checkEntry inputs 3 prefix26 (.mul (variableData 1) 6,claim26)=true := by rfl
theorem step_checked27 : checkEntry inputs 3 prefix27 (.linear [(4,1), (13,1), (20,1), (22,1), (23,1), (25,1), (26,1)],claim27)=true := by rfl
theorem step_checked28 : checkEntry inputs 3 prefix28 (.mul (variableData 2) 6,claim28)=true := by rfl
theorem step_checked29 : checkEntry inputs 3 prefix29 (.linear [(2,1), (6,1), (8,1), (15,1), (16,1), (23,1), (28,1)],claim29)=true := by rfl
theorem step_checked30 : checkEntry inputs 3 prefix30 (.mul (variableData 3) 6,claim30)=true := by rfl
theorem step_checked31 : checkEntry inputs 3 prefix31 (.linear [(8,1), (18,1), (29,1), (30,1)],claim31)=true := by rfl
theorem step_checked32 : checkEntry inputs 3 prefix32 (.mul (variableData 0) 8,claim32)=true := by rfl
theorem step_checked33 : checkEntry inputs 3 prefix33 (.linear [(12,1), (22,1), (23,1), (25,1), (27,1), (32,1)],claim33)=true := by rfl
theorem step_checked34 : checkEntry inputs 3 prefix34 (.mul (variableData 2) 8,claim34)=true := by rfl
theorem step_checked35 : checkEntry inputs 3 prefix35 (.linear [(16,1), (23,1), (29,1), (31,1), (33,1), (34,1)],claim35)=true := by rfl
theorem all_steps (i : Fin 36) : checkEntry inputs 3 ((program.take i.val).map Prod.snd) (program.get i)=true := by
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
  · exact step_checked25
  · exact step_checked26
  · exact step_checked27
  · exact step_checked28
  · exact step_checked29
  · exact step_checked30
  · exact step_checked31
  · exact step_checked32
  · exact step_checked33
  · exact step_checked34
  · exact step_checked35
theorem upper_checked : checkClaimedUnit inputs 3 program=true := by
  have all := checkEntriesFrom_of_each inputs 3 program [] all_steps
  change (checkEntriesFrom inputs 3 [] program && (add claim35 (scale (-1) SparseCertificate.oneData)).isEmpty)=true
  rw [all]
  rfl
theorem one_at_three : InClosure {p | p∈inputs.map meaning} 3 (1 : MvPolynomial (Fin 4) K) :=
  checkClaimedUnit_sound inputs 3 program upper_checked
theorem exact_lfd_three : IsLastFallDegree {p | p∈inputs.map meaning} 3 := by
  apply exactLFD_of_first_one (by omega) one_at_three
  exact one_not_at_two
end ICFallDegrees.GeneratedBinaryKoblitz
