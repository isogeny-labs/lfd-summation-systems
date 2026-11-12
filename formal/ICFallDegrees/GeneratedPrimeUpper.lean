import ICFallDegrees.GeneratedPrimeLower
import ICFallDegrees.SparseDAG
import Mathlib.Tactic.FinCases

/-! Generated untrusted derivation data; every row identity is checked. -/
namespace ICFallDegrees.GeneratedPrimeUpper
open SparsePolynomial SparseDAG GeneratedPrimeLower
set_option maxRecDepth 20000
set_option maxHeartbeats 100000000
def claim0 : Sparse 6 K := [((exponent 2 0 0 0 0 0),1), ((exponent 0 1 0 0 0 0),256)]
def claim1 : Sparse 6 K := [((exponent 0 2 0 0 0 0),1), ((exponent 0 0 1 0 0 0),256)]
def claim2 : Sparse 6 K := [((exponent 0 0 2 0 0 0),1), ((exponent 0 0 0 0 0 0),256)]
def claim3 : Sparse 6 K := [((exponent 0 0 0 2 0 0),1), ((exponent 0 0 0 0 1 0),256)]
def claim4 : Sparse 6 K := [((exponent 0 0 0 0 2 0),1), ((exponent 0 0 0 0 0 1),256)]
def claim5 : Sparse 6 K := [((exponent 0 0 0 0 0 2),1), ((exponent 0 0 0 0 0 0),256)]
def claim6 : Sparse 6 K := [((exponent 1 0 0 1 0 0),117), ((exponent 1 0 0 0 1 0),173), ((exponent 0 1 0 1 0 0),173), ((exponent 0 1 0 0 1 0),1), ((exponent 1 0 0 0 0 0),162), ((exponent 0 1 0 0 0 0),222), ((exponent 0 0 0 1 0 0),162), ((exponent 0 0 0 0 1 0),222), ((exponent 0 0 0 0 0 0),199)]
def claim7 : Sparse 6 K := [((exponent 1 0 0 1 0 0),1), ((exponent 1 0 0 0 1 0),52), ((exponent 0 1 0 1 0 0),52), ((exponent 0 1 0 0 1 0),134), ((exponent 1 0 0 0 0 0),120), ((exponent 0 1 0 0 0 0),193), ((exponent 0 0 0 1 0 0),120), ((exponent 0 0 0 0 1 0),193), ((exponent 0 0 0 0 0 0),195)]
def claim8 : Sparse 6 K := [((exponent 2 0 1 0 0 0),1), ((exponent 0 1 1 0 0 0),256)]
def claim9 : Sparse 6 K := [((exponent 2 0 0 1 0 0),1), ((exponent 0 1 0 1 0 0),256)]
def claim10 : Sparse 6 K := [((exponent 2 0 0 0 1 0),1), ((exponent 0 1 0 0 1 0),256)]
def claim11 : Sparse 6 K := [((exponent 2 0 0 0 0 1),1), ((exponent 0 1 0 0 0 1),256)]
def claim12 : Sparse 6 K := [((exponent 1 2 0 0 0 0),1), ((exponent 1 0 1 0 0 0),256)]
def claim13 : Sparse 6 K := [((exponent 0 2 0 1 0 0),1), ((exponent 0 0 1 1 0 0),256)]
def claim14 : Sparse 6 K := [((exponent 0 2 0 0 1 0),1), ((exponent 0 0 1 0 1 0),256)]
def claim15 : Sparse 6 K := [((exponent 1 0 2 0 0 0),1), ((exponent 1 0 0 0 0 0),256)]
def claim16 : Sparse 6 K := [((exponent 0 1 2 0 0 0),1), ((exponent 0 1 0 0 0 0),256)]
def claim17 : Sparse 6 K := [((exponent 0 0 2 1 0 0),1), ((exponent 0 0 0 1 0 0),256)]
def claim18 : Sparse 6 K := [((exponent 0 0 2 0 1 0),1), ((exponent 0 0 0 0 1 0),256)]
def claim19 : Sparse 6 K := [((exponent 0 0 2 0 0 1),1), ((exponent 0 0 0 0 0 1),256)]
def claim20 : Sparse 6 K := [((exponent 1 0 0 2 0 0),1), ((exponent 1 0 0 0 1 0),256)]
def claim21 : Sparse 6 K := [((exponent 0 1 0 2 0 0),1), ((exponent 0 1 0 0 1 0),256)]
def claim22 : Sparse 6 K := [((exponent 0 0 1 2 0 0),1), ((exponent 0 0 1 0 1 0),256)]
def claim23 : Sparse 6 K := [((exponent 0 0 0 2 1 0),1), ((exponent 0 0 0 0 2 0),256)]
def claim24 : Sparse 6 K := [((exponent 0 0 0 2 1 0),1), ((exponent 0 0 0 0 0 1),256)]
def claim25 : Sparse 6 K := [((exponent 0 0 0 2 0 1),1), ((exponent 0 0 0 0 1 1),256)]
def claim26 : Sparse 6 K := [((exponent 1 0 0 0 2 0),1), ((exponent 1 0 0 0 0 1),256)]
def claim27 : Sparse 6 K := [((exponent 0 1 0 0 2 0),1), ((exponent 0 1 0 0 0 1),256)]
def claim28 : Sparse 6 K := [((exponent 0 0 1 0 2 0),1), ((exponent 0 0 1 0 0 1),256)]
def claim29 : Sparse 6 K := [((exponent 0 0 0 1 2 0),1), ((exponent 0 0 0 1 0 1),256)]
def claim30 : Sparse 6 K := [((exponent 1 0 0 0 0 2),1), ((exponent 1 0 0 0 0 0),256)]
def claim31 : Sparse 6 K := [((exponent 0 1 0 0 0 2),1), ((exponent 0 1 0 0 0 0),256)]
def claim32 : Sparse 6 K := [((exponent 0 0 1 0 0 2),1), ((exponent 0 0 1 0 0 0),256)]
def claim33 : Sparse 6 K := [((exponent 0 0 0 1 0 2),1), ((exponent 0 0 0 1 0 0),256)]
def claim34 : Sparse 6 K := [((exponent 0 0 0 0 1 2),1), ((exponent 0 0 0 0 1 0),256)]
def claim35 : Sparse 6 K := [((exponent 2 0 0 1 0 0),1), ((exponent 2 0 0 0 1 0),52), ((exponent 1 1 0 1 0 0),52), ((exponent 1 1 0 0 1 0),134), ((exponent 2 0 0 0 0 0),120), ((exponent 1 1 0 0 0 0),193), ((exponent 1 0 0 1 0 0),120), ((exponent 1 0 0 0 1 0),193), ((exponent 1 0 0 0 0 0),195)]
def claim36 : Sparse 6 K := [((exponent 1 1 0 1 0 0),1), ((exponent 1 1 0 0 1 0),52), ((exponent 1 1 0 0 0 0),236), ((exponent 1 0 0 0 1 0),116), ((exponent 0 1 0 1 0 0),53), ((exponent 0 1 0 0 1 0),186), ((exponent 1 0 0 0 0 0),226), ((exponent 0 1 0 0 0 0),150), ((exponent 0 0 0 1 0 0),158), ((exponent 0 0 0 0 1 0),207), ((exponent 0 0 0 0 0 0),64)]
def claim37 : Sparse 6 K := [((exponent 1 1 0 1 0 0),1), ((exponent 1 1 0 0 1 0),52), ((exponent 0 2 0 1 0 0),52), ((exponent 0 2 0 0 1 0),134), ((exponent 1 1 0 0 0 0),120), ((exponent 0 2 0 0 0 0),193), ((exponent 0 1 0 1 0 0),120), ((exponent 0 1 0 0 1 0),193), ((exponent 0 1 0 0 0 0),195)]
def claim38 : Sparse 6 K := [((exponent 1 1 0 0 0 0),1), ((exponent 1 0 0 0 1 0),1), ((exponent 0 1 0 1 0 0),139), ((exponent 0 1 0 0 1 0),237), ((exponent 0 0 1 1 0 0),35), ((exponent 0 0 1 0 1 0),21), ((exponent 1 0 0 0 0 0),95), ((exponent 0 1 0 0 0 0),55), ((exponent 0 0 1 0 0 0),36), ((exponent 0 0 0 1 0 0),121), ((exponent 0 0 0 0 1 0),4), ((exponent 0 0 0 0 0 0),36)]
def claim39 : Sparse 6 K := [((exponent 1 0 1 1 0 0),1), ((exponent 1 0 1 0 1 0),52), ((exponent 0 1 1 1 0 0),52), ((exponent 0 1 1 0 1 0),134), ((exponent 1 0 1 0 0 0),120), ((exponent 0 1 1 0 0 0),193), ((exponent 0 0 1 1 0 0),120), ((exponent 0 0 1 0 1 0),193), ((exponent 0 0 1 0 0 0),195)]
def claim40 : Sparse 6 K := [((exponent 1 0 0 2 0 0),1), ((exponent 1 0 0 1 1 0),52), ((exponent 0 1 0 2 0 0),52), ((exponent 0 1 0 1 1 0),134), ((exponent 1 0 0 1 0 0),120), ((exponent 0 1 0 1 0 0),193), ((exponent 0 0 0 2 0 0),120), ((exponent 0 0 0 1 1 0),193), ((exponent 0 0 0 1 0 0),195)]
def claim41 : Sparse 6 K := [((exponent 1 0 0 1 1 0),1), ((exponent 0 1 0 1 1 0),52), ((exponent 1 0 0 0 1 0),53), ((exponent 0 1 0 1 0 0),116), ((exponent 0 1 0 0 1 0),186), ((exponent 0 0 0 1 1 0),236), ((exponent 1 0 0 0 0 0),158), ((exponent 0 1 0 0 0 0),207), ((exponent 0 0 0 1 0 0),226), ((exponent 0 0 0 0 1 0),150), ((exponent 0 0 0 0 0 0),64)]
def claim42 : Sparse 6 K := [((exponent 1 0 0 1 1 0),1), ((exponent 1 0 0 0 2 0),52), ((exponent 0 1 0 1 1 0),52), ((exponent 0 1 0 0 2 0),134), ((exponent 1 0 0 0 1 0),120), ((exponent 0 1 0 0 1 0),193), ((exponent 0 0 0 1 1 0),120), ((exponent 0 0 0 0 2 0),193), ((exponent 0 0 0 0 1 0),195)]
def claim43 : Sparse 6 K := [((exponent 1 0 0 0 1 0),1), ((exponent 1 0 0 0 0 1),89), ((exponent 0 1 0 1 0 0),98), ((exponent 0 1 0 0 1 0),96), ((exponent 0 1 0 0 0 1),2), ((exponent 0 0 0 1 1 0),98), ((exponent 1 0 0 0 0 0),36), ((exponent 0 1 0 0 0 0),135), ((exponent 0 0 0 1 0 0),58), ((exponent 0 0 0 0 1 0),250), ((exponent 0 0 0 0 0 1),187), ((exponent 0 0 0 0 0 0),187)]
def claim44 : Sparse 6 K := [((exponent 1 0 0 1 0 1),1), ((exponent 1 0 0 0 1 1),52), ((exponent 0 1 0 1 0 1),52), ((exponent 0 1 0 0 1 1),134), ((exponent 1 0 0 0 0 1),120), ((exponent 0 1 0 0 0 1),193), ((exponent 0 0 0 1 0 1),120), ((exponent 0 0 0 0 1 1),193), ((exponent 0 0 0 0 0 1),195)]
def claim45 : Sparse 6 K := [((exponent 1 2 0 1 0 0),1), ((exponent 1 0 1 1 0 0),256)]
def claim46 : Sparse 6 K := [((exponent 1 2 0 1 0 0),1), ((exponent 1 0 1 0 1 0),52), ((exponent 0 1 1 1 0 0),52), ((exponent 0 1 1 0 1 0),134), ((exponent 1 0 1 0 0 0),120), ((exponent 0 1 1 0 0 0),193), ((exponent 0 0 1 1 0 0),120), ((exponent 0 0 1 0 1 0),193), ((exponent 0 0 1 0 0 0),195)]
def claim47 : Sparse 6 K := [((exponent 1 2 0 0 1 0),1), ((exponent 1 0 1 0 1 0),256)]
def claim48 : Sparse 6 K := [((exponent 1 0 2 0 0 1),1), ((exponent 1 0 0 0 0 1),256)]
def claim49 : Sparse 6 K := [((exponent 0 1 2 1 0 0),1), ((exponent 0 1 0 1 0 0),256)]
def claim50 : Sparse 6 K := [((exponent 0 1 2 0 1 0),1), ((exponent 0 1 0 0 1 0),256)]
def claim51 : Sparse 6 K := [((exponent 0 1 2 0 0 1),1), ((exponent 0 1 0 0 0 1),256)]
def claim52 : Sparse 6 K := [((exponent 0 0 2 1 1 0),1), ((exponent 0 0 0 1 1 0),256)]
def claim53 : Sparse 6 K := [((exponent 0 0 2 1 0 1),1), ((exponent 0 0 0 1 0 1),256)]
def claim54 : Sparse 6 K := [((exponent 0 0 2 0 1 1),1), ((exponent 0 0 0 0 1 1),256)]
def claim55 : Sparse 6 K := [((exponent 1 1 0 2 0 0),1), ((exponent 1 1 0 0 1 0),256)]
def claim56 : Sparse 6 K := [((exponent 1 0 1 2 0 0),1), ((exponent 1 0 1 0 1 0),256)]
def claim57 : Sparse 6 K := [((exponent 0 1 1 2 0 0),1), ((exponent 0 1 1 0 1 0),256)]
def claim58 : Sparse 6 K := [((exponent 0 1 0 2 1 0),1), ((exponent 0 1 0 0 2 0),256)]
def claim59 : Sparse 6 K := [((exponent 0 1 0 2 1 0),1), ((exponent 0 1 0 0 0 1),256)]
def claim60 : Sparse 6 K := [((exponent 0 0 1 2 1 0),1), ((exponent 0 0 1 0 2 0),256)]
def claim61 : Sparse 6 K := [((exponent 0 0 1 2 1 0),1), ((exponent 0 0 1 0 0 1),256)]
def claim62 : Sparse 6 K := [((exponent 0 0 1 2 0 1),1), ((exponent 0 0 1 0 1 1),256)]
def claim63 : Sparse 6 K := [((exponent 0 0 0 2 1 1),1), ((exponent 0 0 0 0 0 2),256)]
def claim64 : Sparse 6 K := [((exponent 0 0 0 2 1 1),1), ((exponent 0 0 0 0 0 0),256)]
def claim65 : Sparse 6 K := [((exponent 1 1 0 0 2 0),1), ((exponent 1 1 0 0 0 1),256)]
def claim66 : Sparse 6 K := [((exponent 1 0 1 0 2 0),1), ((exponent 1 0 1 0 0 1),256)]
def claim67 : Sparse 6 K := [((exponent 1 0 0 1 2 0),1), ((exponent 1 0 0 1 0 1),256)]
def claim68 : Sparse 6 K := [((exponent 1 0 0 1 2 0),1), ((exponent 1 0 0 0 1 1),52), ((exponent 0 1 0 1 0 1),52), ((exponent 0 1 0 0 1 1),134), ((exponent 1 0 0 0 0 1),120), ((exponent 0 1 0 0 0 1),193), ((exponent 0 0 0 1 0 1),120), ((exponent 0 0 0 0 1 1),193), ((exponent 0 0 0 0 0 1),195)]
def claim69 : Sparse 6 K := [((exponent 0 1 1 0 2 0),1), ((exponent 0 1 1 0 0 1),256)]
def claim70 : Sparse 6 K := [((exponent 0 1 0 1 2 0),1), ((exponent 0 1 0 1 0 1),256)]
def claim71 : Sparse 6 K := [((exponent 0 0 1 1 2 0),1), ((exponent 0 0 1 1 0 1),256)]
def claim72 : Sparse 6 K := [((exponent 1 1 0 0 0 2),1), ((exponent 1 1 0 0 0 0),256)]
def claim73 : Sparse 6 K := [((exponent 1 1 0 0 0 2),1), ((exponent 1 0 0 0 0 1),168), ((exponent 0 1 0 1 0 0),41), ((exponent 0 1 0 0 1 0),141), ((exponent 0 1 0 0 0 1),255), ((exponent 0 0 1 1 0 0),35), ((exponent 0 0 1 0 1 0),21), ((exponent 0 0 0 1 1 0),159), ((exponent 1 0 0 0 0 0),59), ((exponent 0 1 0 0 0 0),177), ((exponent 0 0 1 0 0 0),36), ((exponent 0 0 0 1 0 0),63), ((exponent 0 0 0 0 1 0),11), ((exponent 0 0 0 0 0 1),70), ((exponent 0 0 0 0 0 0),106)]
def claim74 : Sparse 6 K := [((exponent 1 0 1 0 0 2),1), ((exponent 1 0 1 0 0 0),256)]
def claim75 : Sparse 6 K := [((exponent 0 1 1 0 0 2),1), ((exponent 0 1 1 0 0 0),256)]
def claim76 : Sparse 6 K := [((exponent 0 1 0 1 0 2),1), ((exponent 0 1 0 1 0 0),256)]
def claim77 : Sparse 6 K := [((exponent 0 1 0 0 1 2),1), ((exponent 0 1 0 0 1 0),256)]
def claim78 : Sparse 6 K := [((exponent 0 0 1 1 0 2),1), ((exponent 0 0 1 1 0 0),256)]
def claim79 : Sparse 6 K := [((exponent 0 0 1 0 1 2),1), ((exponent 0 0 1 0 1 0),256)]
def claim80 : Sparse 6 K := [((exponent 0 0 0 1 1 2),1), ((exponent 0 0 0 1 1 0),256)]
def claim81 : Sparse 6 K := [((exponent 1 2 0 1 0 0),1), ((exponent 1 2 0 0 1 0),52), ((exponent 1 2 0 0 0 0),236), ((exponent 1 1 0 0 1 0),116), ((exponent 0 2 0 1 0 0),53), ((exponent 0 2 0 0 1 0),186), ((exponent 1 1 0 0 0 0),226), ((exponent 0 2 0 0 0 0),150), ((exponent 0 1 0 1 0 0),158), ((exponent 0 1 0 0 1 0),207), ((exponent 0 1 0 0 0 0),64)]
def claim82 : Sparse 6 K := [((exponent 1 1 0 0 1 0),1), ((exponent 0 1 1 1 0 0),35), ((exponent 0 1 1 0 1 0),21), ((exponent 1 0 1 0 0 0),1), ((exponent 1 0 0 0 0 1),231), ((exponent 0 1 1 0 0 0),36), ((exponent 0 1 0 1 0 0),81), ((exponent 0 1 0 0 1 0),230), ((exponent 0 1 0 0 0 1),190), ((exponent 0 0 1 1 0 0),155), ((exponent 0 0 1 0 1 0),41), ((exponent 0 0 0 1 1 0),58), ((exponent 1 0 0 0 0 0),49), ((exponent 0 1 0 0 0 0),183), ((exponent 0 0 1 0 0 0),233), ((exponent 0 0 0 1 0 0),183), ((exponent 0 0 0 0 1 0),240), ((exponent 0 0 0 0 0 1),32), ((exponent 0 0 0 0 0 0),210)]
def claim83 : Sparse 6 K := [((exponent 1 1 1 1 0 0),1), ((exponent 1 1 1 0 1 0),52), ((exponent 1 1 1 0 0 0),236), ((exponent 1 0 1 0 1 0),116), ((exponent 0 1 1 1 0 0),53), ((exponent 0 1 1 0 1 0),186), ((exponent 1 0 1 0 0 0),226), ((exponent 0 1 1 0 0 0),150), ((exponent 0 0 1 1 0 0),158), ((exponent 0 0 1 0 1 0),207), ((exponent 0 0 1 0 0 0),64)]
def claim84 : Sparse 6 K := [((exponent 1 1 0 2 0 0),1), ((exponent 1 1 0 1 1 0),52), ((exponent 1 1 0 1 0 0),236), ((exponent 1 0 0 1 1 0),116), ((exponent 0 1 0 2 0 0),53), ((exponent 0 1 0 1 1 0),186), ((exponent 1 0 0 1 0 0),226), ((exponent 0 1 0 1 0 0),150), ((exponent 0 0 0 2 0 0),158), ((exponent 0 0 0 1 1 0),207), ((exponent 0 0 0 1 0 0),64)]
def claim85 : Sparse 6 K := [((exponent 1 1 0 1 1 0),1), ((exponent 0 1 1 1 0 0),149), ((exponent 0 1 1 0 1 0),38), ((exponent 0 1 0 1 1 0),194), ((exponent 1 0 1 0 0 0),63), ((exponent 1 0 0 0 0 1),40), ((exponent 0 1 1 0 0 0),212), ((exponent 0 1 0 1 0 0),32), ((exponent 0 1 0 0 1 0),98), ((exponent 0 1 0 0 0 1),24), ((exponent 0 0 1 1 0 0),24), ((exponent 0 0 1 0 1 0),28), ((exponent 0 0 0 1 1 0),184), ((exponent 1 0 0 0 0 0),102), ((exponent 0 1 0 0 0 0),196), ((exponent 0 0 1 0 0 0),19), ((exponent 0 0 0 1 0 0),76), ((exponent 0 0 0 0 1 0),35), ((exponent 0 0 0 0 0 1),188), ((exponent 0 0 0 0 0 0),196)]
def claim86 : Sparse 6 K := [((exponent 1 1 0 1 1 0),1), ((exponent 1 1 0 0 2 0),52), ((exponent 1 1 0 0 1 0),236), ((exponent 1 0 0 0 2 0),116), ((exponent 0 1 0 1 1 0),53), ((exponent 0 1 0 0 2 0),186), ((exponent 1 0 0 0 1 0),226), ((exponent 0 1 0 0 1 0),150), ((exponent 0 0 0 1 1 0),158), ((exponent 0 0 0 0 2 0),207), ((exponent 0 0 0 0 1 0),64)]
def claim87 : Sparse 6 K := [((exponent 1 1 0 0 0 1),1), ((exponent 0 1 1 1 0 0),120), ((exponent 0 1 1 0 1 0),72), ((exponent 0 1 0 1 1 0),22), ((exponent 1 0 1 0 0 0),187), ((exponent 1 0 0 0 0 1),217), ((exponent 0 1 1 0 0 0),50), ((exponent 0 1 0 1 0 0),135), ((exponent 0 1 0 0 1 0),161), ((exponent 0 1 0 0 0 1),170), ((exponent 0 0 1 1 0 0),245), ((exponent 0 0 1 0 1 0),189), ((exponent 0 0 0 1 1 0),111), ((exponent 1 0 0 0 0 0),64), ((exponent 0 1 0 0 0 0),32), ((exponent 0 0 1 0 0 0),242), ((exponent 0 0 0 1 0 0),23), ((exponent 0 0 0 0 1 0),34), ((exponent 0 0 0 0 0 1),105), ((exponent 0 0 0 0 0 0),237)]
def claim88 : Sparse 6 K := [((exponent 1 1 0 1 0 1),1), ((exponent 1 1 0 0 1 1),52), ((exponent 1 1 0 0 0 1),236), ((exponent 1 0 0 0 1 1),116), ((exponent 0 1 0 1 0 1),53), ((exponent 0 1 0 0 1 1),186), ((exponent 1 0 0 0 0 1),226), ((exponent 0 1 0 0 0 1),150), ((exponent 0 0 0 1 0 1),158), ((exponent 0 0 0 0 1 1),207), ((exponent 0 0 0 0 0 1),64)]
def claim89 : Sparse 6 K := [((exponent 1 1 0 1 0 1),1), ((exponent 1 1 0 0 1 1),52), ((exponent 1 0 0 0 1 1),116), ((exponent 0 1 1 1 0 0),207), ((exponent 0 1 1 0 1 0),227), ((exponent 0 1 0 1 1 0),205), ((exponent 0 1 0 1 0 1),53), ((exponent 0 1 0 0 1 1),186), ((exponent 1 0 1 0 0 0),72), ((exponent 1 0 0 0 0 1),157), ((exponent 0 1 1 0 0 0),22), ((exponent 0 1 0 1 0 0),8), ((exponent 0 1 0 0 1 0),40), ((exponent 0 1 0 0 0 1),122), ((exponent 0 0 1 1 0 0),5), ((exponent 0 0 1 0 1 0),114), ((exponent 0 0 0 1 1 0),18), ((exponent 0 0 0 1 0 1),158), ((exponent 0 0 0 0 1 1),207), ((exponent 1 0 0 0 0 0),59), ((exponent 0 1 0 0 0 0),158), ((exponent 0 0 1 0 0 0),199), ((exponent 0 0 0 1 0 0),226), ((exponent 0 0 0 0 1 0),200), ((exponent 0 0 0 0 0 1),213), ((exponent 0 0 0 0 0 0),94)]
def claim90 : Sparse 6 K := [((exponent 1 1 1 0 0 0),1), ((exponent 1 0 1 0 1 0),1), ((exponent 0 1 1 1 0 0),139), ((exponent 0 1 1 0 1 0),237), ((exponent 0 0 2 1 0 0),35), ((exponent 0 0 2 0 1 0),21), ((exponent 1 0 1 0 0 0),95), ((exponent 0 1 1 0 0 0),55), ((exponent 0 0 2 0 0 0),36), ((exponent 0 0 1 1 0 0),121), ((exponent 0 0 1 0 1 0),4), ((exponent 0 0 1 0 0 0),36)]
def claim91 : Sparse 6 K := [((exponent 1 1 1 0 0 0),1), ((exponent 1 0 1 0 1 0),1), ((exponent 0 1 1 1 0 0),139), ((exponent 0 1 1 0 1 0),237), ((exponent 1 0 1 0 0 0),95), ((exponent 0 1 1 0 0 0),55), ((exponent 0 0 1 1 0 0),121), ((exponent 0 0 1 0 1 0),4), ((exponent 0 0 1 0 0 0),36), ((exponent 0 0 0 1 0 0),35), ((exponent 0 0 0 0 1 0),21), ((exponent 0 0 0 0 0 0),36)]
def claim92 : Sparse 6 K := [((exponent 1 1 0 1 0 0),1), ((exponent 1 0 0 1 1 0),1), ((exponent 0 1 0 2 0 0),139), ((exponent 0 1 0 1 1 0),237), ((exponent 0 0 1 2 0 0),35), ((exponent 0 0 1 1 1 0),21), ((exponent 1 0 0 1 0 0),95), ((exponent 0 1 0 1 0 0),55), ((exponent 0 0 1 1 0 0),36), ((exponent 0 0 0 2 0 0),121), ((exponent 0 0 0 1 1 0),4), ((exponent 0 0 0 1 0 0),36)]
def claim93 : Sparse 6 K := [((exponent 0 1 1 1 0 0),1), ((exponent 0 1 1 0 1 0),52), ((exponent 0 1 0 1 1 0),70), ((exponent 0 0 1 1 1 0),1), ((exponent 1 0 1 0 0 0),235), ((exponent 1 0 0 0 0 1),138), ((exponent 0 1 1 0 0 0),236), ((exponent 0 1 0 1 0 0),19), ((exponent 0 1 0 0 1 0),208), ((exponent 0 1 0 0 0 1),237), ((exponent 0 0 1 1 0 0),118), ((exponent 0 0 1 0 1 0),21), ((exponent 0 0 0 1 1 0),245), ((exponent 1 0 0 0 0 0),11), ((exponent 0 1 0 0 0 0),227), ((exponent 0 0 1 0 0 0),235), ((exponent 0 0 0 1 0 0),100), ((exponent 0 0 0 0 1 0),3), ((exponent 0 0 0 0 0 1),186), ((exponent 0 0 0 0 0 0),241)]
def claim94 : Sparse 6 K := [((exponent 1 1 0 0 1 0),1), ((exponent 1 0 0 0 2 0),1), ((exponent 0 1 0 1 1 0),139), ((exponent 0 1 0 0 2 0),237), ((exponent 0 0 1 1 1 0),35), ((exponent 0 0 1 0 2 0),21), ((exponent 1 0 0 0 1 0),95), ((exponent 0 1 0 0 1 0),55), ((exponent 0 0 1 0 1 0),36), ((exponent 0 0 0 1 1 0),121), ((exponent 0 0 0 0 2 0),4), ((exponent 0 0 0 0 1 0),36)]
def claim95 : Sparse 6 K := [((exponent 0 1 0 1 1 0),1), ((exponent 0 0 1 1 1 0),166), ((exponent 0 1 0 1 0 0),190), ((exponent 0 1 0 0 1 0),178), ((exponent 0 1 0 0 0 1),145), ((exponent 0 0 1 1 0 0),101), ((exponent 0 0 1 0 1 0),79), ((exponent 0 0 1 0 0 1),204), ((exponent 0 0 0 1 1 0),154), ((exponent 0 1 0 0 0 0),234), ((exponent 0 0 1 0 0 0),96), ((exponent 0 0 0 1 0 0),101), ((exponent 0 0 0 0 1 0),138), ((exponent 0 0 0 0 0 1),96), ((exponent 0 0 0 0 0 0),66)]
def claim96 : Sparse 6 K := [((exponent 1 1 0 0 0 1),1), ((exponent 1 0 0 0 1 1),1), ((exponent 0 1 0 1 0 1),139), ((exponent 0 1 0 0 1 1),237), ((exponent 0 0 1 1 0 1),35), ((exponent 0 0 1 0 1 1),21), ((exponent 1 0 0 0 0 1),95), ((exponent 0 1 0 0 0 1),55), ((exponent 0 0 1 0 0 1),36), ((exponent 0 0 0 1 0 1),121), ((exponent 0 0 0 0 1 1),4), ((exponent 0 0 0 0 0 1),36)]
def claim97 : Sparse 6 K := [((exponent 1 0 0 0 1 1),1), ((exponent 0 1 0 1 0 1),139), ((exponent 0 1 0 0 1 1),237), ((exponent 0 0 1 1 1 0),256), ((exponent 0 0 1 1 0 1),35), ((exponent 0 0 1 0 1 1),21), ((exponent 1 0 0 0 0 1),247), ((exponent 0 1 0 1 0 0),127), ((exponent 0 1 0 0 1 0),214), ((exponent 0 1 0 0 0 1),84), ((exponent 0 0 1 1 0 0),160), ((exponent 0 0 1 0 1 0),188), ((exponent 0 0 1 0 0 1),231), ((exponent 0 0 0 1 1 0),176), ((exponent 0 0 0 1 0 1),121), ((exponent 0 0 0 0 1 1),4), ((exponent 1 0 0 0 0 0),228), ((exponent 0 1 0 0 0 0),167), ((exponent 0 0 1 0 0 0),67), ((exponent 0 0 0 1 0 0),21), ((exponent 0 0 0 0 1 0),148), ((exponent 0 0 0 0 0 1),14), ((exponent 0 0 0 0 0 0),15)]
def claim98 : Sparse 6 K := [((exponent 1 0 1 2 0 0),1), ((exponent 1 0 1 1 1 0),52), ((exponent 0 1 1 2 0 0),52), ((exponent 0 1 1 1 1 0),134), ((exponent 1 0 1 1 0 0),120), ((exponent 0 1 1 1 0 0),193), ((exponent 0 0 1 2 0 0),120), ((exponent 0 0 1 1 1 0),193), ((exponent 0 0 1 1 0 0),195)]
def claim99 : Sparse 6 K := [((exponent 1 0 1 1 1 0),1), ((exponent 0 1 1 1 1 0),52), ((exponent 1 0 1 0 1 0),53), ((exponent 0 1 1 0 1 0),65), ((exponent 0 0 1 1 1 0),75), ((exponent 1 0 1 0 0 0),140), ((exponent 1 0 0 0 0 1),183), ((exponent 0 1 1 0 0 0),73), ((exponent 0 1 0 1 0 0),138), ((exponent 0 1 0 0 1 0),22), ((exponent 0 1 0 0 0 1),90), ((exponent 0 0 1 1 0 0),192), ((exponent 0 0 1 0 1 0),35), ((exponent 0 0 1 0 0 1),115), ((exponent 0 0 0 1 1 0),25), ((exponent 1 0 0 0 0 0),9), ((exponent 0 1 0 0 0 0),218), ((exponent 0 0 1 0 0 0),85), ((exponent 0 0 0 1 0 0),255), ((exponent 0 0 0 0 1 0),206), ((exponent 0 0 0 0 0 1),51), ((exponent 0 0 0 0 0 0),132)]
def claim100 : Sparse 6 K := [((exponent 1 0 1 1 1 0),1), ((exponent 1 0 1 0 2 0),52), ((exponent 0 1 1 1 1 0),52), ((exponent 0 1 1 0 2 0),134), ((exponent 1 0 1 0 1 0),120), ((exponent 0 1 1 0 1 0),193), ((exponent 0 0 1 1 1 0),120), ((exponent 0 0 1 0 2 0),193), ((exponent 0 0 1 0 1 0),195)]
def claim101 : Sparse 6 K := [((exponent 1 0 1 0 1 0),1), ((exponent 1 0 1 0 0 1),89), ((exponent 0 1 1 0 1 0),140), ((exponent 0 1 1 0 0 1),2), ((exponent 0 0 1 1 1 0),250), ((exponent 1 0 1 0 0 0),136), ((exponent 1 0 0 0 0 1),97), ((exponent 0 1 1 0 0 0),137), ((exponent 0 1 0 1 0 0),90), ((exponent 0 1 0 0 1 0),249), ((exponent 0 1 0 0 0 1),14), ((exponent 0 0 1 1 0 0),47), ((exponent 0 0 1 0 1 0),175), ((exponent 0 0 1 0 0 1),5), ((exponent 0 0 0 1 1 0),61), ((exponent 1 0 0 0 0 0),207), ((exponent 0 1 0 0 0 0),131), ((exponent 0 0 1 0 0 0),156), ((exponent 0 0 0 1 0 0),211), ((exponent 0 0 0 0 1 0),112), ((exponent 0 0 0 0 0 1),145), ((exponent 0 0 0 0 0 0),209)]
def claim102 : Sparse 6 K := [((exponent 1 0 1 1 0 1),1), ((exponent 1 0 1 0 1 1),52), ((exponent 0 1 1 1 0 1),52), ((exponent 0 1 1 0 1 1),134), ((exponent 1 0 1 0 0 1),120), ((exponent 0 1 1 0 0 1),193), ((exponent 0 0 1 1 0 1),120), ((exponent 0 0 1 0 1 1),193), ((exponent 0 0 1 0 0 1),195)]
def claim103 : Sparse 6 K := [((exponent 1 0 0 1 2 0),1), ((exponent 0 1 0 1 2 0),52), ((exponent 1 0 0 0 2 0),53), ((exponent 0 1 0 1 1 0),116), ((exponent 0 1 0 0 2 0),186), ((exponent 0 0 0 1 2 0),236), ((exponent 1 0 0 0 1 0),158), ((exponent 0 1 0 0 1 0),207), ((exponent 0 0 0 1 1 0),226), ((exponent 0 0 0 0 2 0),150), ((exponent 0 0 0 0 1 0),64)]
def claim104 : Sparse 6 K := [((exponent 0 1 0 1 0 1),1), ((exponent 0 1 0 0 1 1),140), ((exponent 0 0 1 1 1 0),7), ((exponent 0 0 1 1 0 1),89), ((exponent 0 0 1 0 1 1),2), ((exponent 0 1 0 1 0 0),126), ((exponent 0 1 0 0 1 0),14), ((exponent 0 1 0 0 0 1),171), ((exponent 0 0 1 1 0 0),183), ((exponent 0 0 1 0 1 0),244), ((exponent 0 0 1 0 0 1),182), ((exponent 0 0 0 1 1 0),38), ((exponent 0 0 0 1 0 1),136), ((exponent 0 0 0 0 1 1),137), ((exponent 0 1 0 0 0 0),154), ((exponent 0 0 1 0 0 0),50), ((exponent 0 0 0 1 0 0),249), ((exponent 0 0 0 0 1 0),82), ((exponent 0 0 0 0 0 1),27), ((exponent 0 0 0 0 0 0),195)]
def claim105 : Sparse 6 K := [((exponent 1 0 0 0 1 1),1), ((exponent 1 0 0 0 0 2),89), ((exponent 0 1 0 1 0 1),98), ((exponent 0 1 0 0 1 1),96), ((exponent 0 1 0 0 0 2),2), ((exponent 0 0 0 1 1 1),98), ((exponent 1 0 0 0 0 1),36), ((exponent 0 1 0 0 0 1),135), ((exponent 0 0 0 1 0 1),58), ((exponent 0 0 0 0 1 1),250), ((exponent 0 0 0 0 0 2),187), ((exponent 0 0 0 0 0 1),187)]
def claim106 : Sparse 6 K := [((exponent 0 1 0 0 1 1),1), ((exponent 0 0 1 1 1 0),177), ((exponent 0 0 1 1 0 1),224), ((exponent 0 0 1 0 1 1),83), ((exponent 0 0 0 1 1 1),87), ((exponent 1 0 0 0 0 1),130), ((exponent 0 1 0 1 0 0),128), ((exponent 0 1 0 0 1 0),157), ((exponent 0 1 0 0 0 1),180), ((exponent 0 0 1 1 0 0),2), ((exponent 0 0 1 0 1 0),186), ((exponent 0 0 1 0 0 1),233), ((exponent 0 0 0 1 1 0),73), ((exponent 0 0 0 1 0 1),82), ((exponent 0 0 0 0 1 1),99), ((exponent 1 0 0 0 0 0),110), ((exponent 0 1 0 0 0 0),248), ((exponent 0 0 1 0 0 0),6), ((exponent 0 0 0 1 0 0),254), ((exponent 0 0 0 0 1 0),21), ((exponent 0 0 0 0 0 1),187), ((exponent 0 0 0 0 0 0),230)]
def claim107 : Sparse 6 K := [((exponent 1 1 1 0 1 0),1), ((exponent 0 1 2 1 0 0),35), ((exponent 0 1 2 0 1 0),21), ((exponent 1 0 2 0 0 0),1), ((exponent 1 0 1 0 0 1),231), ((exponent 0 1 2 0 0 0),36), ((exponent 0 1 1 1 0 0),81), ((exponent 0 1 1 0 1 0),230), ((exponent 0 1 1 0 0 1),190), ((exponent 0 0 2 1 0 0),155), ((exponent 0 0 2 0 1 0),41), ((exponent 0 0 1 1 1 0),58), ((exponent 1 0 1 0 0 0),49), ((exponent 0 1 1 0 0 0),183), ((exponent 0 0 2 0 0 0),233), ((exponent 0 0 1 1 0 0),183), ((exponent 0 0 1 0 1 0),240), ((exponent 0 0 1 0 0 1),32), ((exponent 0 0 1 0 0 0),210)]
def claim108 : Sparse 6 K := [((exponent 1 1 1 0 1 0),1), ((exponent 1 0 1 0 0 1),231), ((exponent 0 1 1 0 1 0),130), ((exponent 0 1 1 0 0 1),190), ((exponent 0 0 1 1 1 0),63), ((exponent 1 0 1 0 0 0),32), ((exponent 1 0 0 0 0 1),130), ((exponent 0 1 1 0 0 0),85), ((exponent 0 1 0 1 0 0),251), ((exponent 0 1 0 0 1 0),156), ((exponent 0 1 0 0 0 1),85), ((exponent 0 0 1 1 0 0),208), ((exponent 0 0 1 0 1 0),60), ((exponent 0 0 1 0 0 1),212), ((exponent 0 0 0 1 1 0),95), ((exponent 1 0 0 0 0 0),138), ((exponent 0 1 0 0 0 0),42), ((exponent 0 0 1 0 0 0),187), ((exponent 0 0 0 1 0 0),96), ((exponent 0 0 0 0 1 0),207), ((exponent 0 0 0 0 0 1),91), ((exponent 0 0 0 0 0 0),15)]
def claim109 : Sparse 6 K := [((exponent 1 1 0 1 1 0),1), ((exponent 0 1 1 2 0 0),35), ((exponent 0 1 1 1 1 0),21), ((exponent 1 0 1 1 0 0),1), ((exponent 1 0 0 1 0 1),231), ((exponent 0 1 1 1 0 0),36), ((exponent 0 1 0 2 0 0),81), ((exponent 0 1 0 1 1 0),230), ((exponent 0 1 0 1 0 1),190), ((exponent 0 0 1 2 0 0),155), ((exponent 0 0 1 1 1 0),41), ((exponent 0 0 0 2 1 0),58), ((exponent 1 0 0 1 0 0),49), ((exponent 0 1 0 1 0 0),183), ((exponent 0 0 1 1 0 0),233), ((exponent 0 0 0 2 0 0),183), ((exponent 0 0 0 1 1 0),240), ((exponent 0 0 0 1 0 1),32), ((exponent 0 0 0 1 0 0),210)]
def claim110 : Sparse 6 K := [((exponent 0 1 1 1 1 0),1), ((exponent 1 0 1 0 0 1),98), ((exponent 0 1 1 0 1 0),198), ((exponent 0 1 1 0 0 1),213), ((exponent 0 0 1 1 1 0),107), ((exponent 1 0 1 0 0 0),94), ((exponent 1 0 0 0 0 1),186), ((exponent 0 1 1 0 0 0),106), ((exponent 0 1 0 1 0 0),73), ((exponent 0 1 0 0 1 0),209), ((exponent 0 1 0 0 0 1),148), ((exponent 0 0 1 1 0 0),50), ((exponent 0 0 1 0 1 0),1), ((exponent 0 0 1 0 0 1),163), ((exponent 0 0 0 1 1 0),245), ((exponent 1 0 0 0 0 0),63), ((exponent 0 1 0 0 0 0),37), ((exponent 0 0 1 0 0 0),2), ((exponent 0 0 0 1 0 0),85), ((exponent 0 0 0 0 1 0),215), ((exponent 0 0 0 0 0 1),72), ((exponent 0 0 0 0 0 0),74)]
def claim111 : Sparse 6 K := [((exponent 1 1 0 0 2 0),1), ((exponent 0 1 1 1 1 0),35), ((exponent 0 1 1 0 2 0),21), ((exponent 1 0 1 0 1 0),1), ((exponent 1 0 0 0 1 1),231), ((exponent 0 1 1 0 1 0),36), ((exponent 0 1 0 1 1 0),81), ((exponent 0 1 0 0 2 0),230), ((exponent 0 1 0 0 1 1),190), ((exponent 0 0 1 1 1 0),155), ((exponent 0 0 1 0 2 0),41), ((exponent 0 0 0 1 2 0),58), ((exponent 1 0 0 0 1 0),49), ((exponent 0 1 0 0 1 0),183), ((exponent 0 0 1 0 1 0),233), ((exponent 0 0 0 1 1 0),183), ((exponent 0 0 0 0 2 0),240), ((exponent 0 0 0 0 1 1),32), ((exponent 0 0 0 0 1 0),210)]
def claim112 : Sparse 6 K := [((exponent 1 0 1 0 0 1),1), ((exponent 0 1 1 0 1 0),207), ((exponent 0 1 1 0 0 1),36), ((exponent 0 0 1 1 1 0),108), ((exponent 1 0 1 0 0 0),77), ((exponent 1 0 0 0 0 1),68), ((exponent 0 1 1 0 0 0),153), ((exponent 0 1 0 1 0 0),181), ((exponent 0 1 0 0 1 0),174), ((exponent 0 1 0 0 0 1),158), ((exponent 0 0 1 1 0 0),80), ((exponent 0 0 1 0 1 0),61), ((exponent 0 0 1 0 0 1),60), ((exponent 0 0 0 1 1 0),218), ((exponent 1 0 0 0 0 0),217), ((exponent 0 1 0 0 0 0),16), ((exponent 0 0 1 0 0 0),188), ((exponent 0 0 0 1 0 0),151), ((exponent 0 0 0 0 1 0),255), ((exponent 0 0 0 0 0 1),14), ((exponent 0 0 0 0 0 0),215)]
def claim113 : Sparse 6 K := [((exponent 1 1 0 0 1 1),1), ((exponent 0 1 1 1 0 1),35), ((exponent 0 1 1 0 1 1),21), ((exponent 1 0 1 0 0 1),1), ((exponent 1 0 0 0 0 2),231), ((exponent 0 1 1 0 0 1),36), ((exponent 0 1 0 1 0 1),81), ((exponent 0 1 0 0 1 1),230), ((exponent 0 1 0 0 0 2),190), ((exponent 0 0 1 1 0 1),155), ((exponent 0 0 1 0 1 1),41), ((exponent 0 0 0 1 1 1),58), ((exponent 1 0 0 0 0 1),49), ((exponent 0 1 0 0 0 1),183), ((exponent 0 0 1 0 0 1),233), ((exponent 0 0 0 1 0 1),183), ((exponent 0 0 0 0 1 1),240), ((exponent 0 0 0 0 0 2),32), ((exponent 0 0 0 0 0 1),210)]
def claim114 : Sparse 6 K := [((exponent 1 1 0 0 1 1),1), ((exponent 0 1 1 1 0 1),35), ((exponent 0 1 1 0 1 1),21), ((exponent 0 1 1 0 1 0),50), ((exponent 0 0 1 1 1 0),2), ((exponent 0 0 1 1 0 1),251), ((exponent 0 0 1 0 1 1),150), ((exponent 0 0 0 1 1 1),51), ((exponent 1 0 1 0 0 0),180), ((exponent 1 0 0 0 0 1),198), ((exponent 0 1 1 0 0 0),104), ((exponent 0 1 0 1 0 0),249), ((exponent 0 1 0 0 1 0),245), ((exponent 0 1 0 0 0 1),135), ((exponent 0 0 1 1 0 0),121), ((exponent 0 0 1 0 1 0),144), ((exponent 0 0 1 0 0 1),206), ((exponent 0 0 0 1 1 0),240), ((exponent 0 0 0 1 0 1),173), ((exponent 0 0 0 0 1 1),124), ((exponent 1 0 0 0 0 0),79), ((exponent 0 1 0 0 0 0),19), ((exponent 0 0 1 0 0 0),228), ((exponent 0 0 0 1 0 0),63), ((exponent 0 0 0 0 1 0),253), ((exponent 0 0 0 0 0 1),47), ((exponent 0 0 0 0 0 0),162)]
def claim115 : Sparse 6 K := [((exponent 1 1 1 0 0 1),1), ((exponent 0 1 2 1 0 0),120), ((exponent 0 1 2 0 1 0),72), ((exponent 0 1 1 1 1 0),22), ((exponent 1 0 2 0 0 0),187), ((exponent 1 0 1 0 0 1),217), ((exponent 0 1 2 0 0 0),50), ((exponent 0 1 1 1 0 0),135), ((exponent 0 1 1 0 1 0),161), ((exponent 0 1 1 0 0 1),170), ((exponent 0 0 2 1 0 0),245), ((exponent 0 0 2 0 1 0),189), ((exponent 0 0 1 1 1 0),111), ((exponent 1 0 1 0 0 0),64), ((exponent 0 1 1 0 0 0),32), ((exponent 0 0 2 0 0 0),242), ((exponent 0 0 1 1 0 0),23), ((exponent 0 0 1 0 1 0),34), ((exponent 0 0 1 0 0 1),105), ((exponent 0 0 1 0 0 0),237)]
def claim116 : Sparse 6 K := [((exponent 1 1 1 0 0 1),1), ((exponent 0 1 1 0 1 0),32), ((exponent 0 1 1 0 0 1),10), ((exponent 0 0 1 1 1 0),121), ((exponent 1 0 1 0 0 0),181), ((exponent 1 0 0 0 0 1),162), ((exponent 0 1 1 0 0 0),110), ((exponent 0 1 0 1 0 0),56), ((exponent 0 1 0 0 1 0),13), ((exponent 0 1 0 0 0 1),159), ((exponent 0 0 1 1 0 0),56), ((exponent 0 0 1 0 1 0),28), ((exponent 0 0 1 0 0 1),79), ((exponent 0 0 0 1 1 0),188), ((exponent 1 0 0 0 0 0),197), ((exponent 0 1 0 0 0 0),201), ((exponent 0 0 1 0 0 0),175), ((exponent 0 0 0 1 0 0),56), ((exponent 0 0 0 0 1 0),253), ((exponent 0 0 0 0 0 1),185), ((exponent 0 0 0 0 0 0),252)]
def claim117 : Sparse 6 K := [((exponent 1 1 0 1 0 1),1), ((exponent 0 1 1 2 0 0),120), ((exponent 0 1 1 1 1 0),72), ((exponent 0 1 0 2 1 0),22), ((exponent 1 0 1 1 0 0),187), ((exponent 1 0 0 1 0 1),217), ((exponent 0 1 1 1 0 0),50), ((exponent 0 1 0 2 0 0),135), ((exponent 0 1 0 1 1 0),161), ((exponent 0 1 0 1 0 1),170), ((exponent 0 0 1 2 0 0),245), ((exponent 0 0 1 1 1 0),189), ((exponent 0 0 0 2 1 0),111), ((exponent 1 0 0 1 0 0),64), ((exponent 0 1 0 1 0 0),32), ((exponent 0 0 1 1 0 0),242), ((exponent 0 0 0 2 0 0),23), ((exponent 0 0 0 1 1 0),34), ((exponent 0 0 0 1 0 1),105), ((exponent 0 0 0 1 0 0),237)]
def claim118 : Sparse 6 K := [((exponent 0 1 1 1 0 1),1), ((exponent 0 1 1 0 1 1),52), ((exponent 0 1 1 0 1 0),185), ((exponent 0 0 1 1 1 0),160), ((exponent 0 0 1 1 0 1),60), ((exponent 0 0 1 0 1 1),152), ((exponent 0 0 0 1 1 1),9), ((exponent 1 0 1 0 0 0),152), ((exponent 1 0 0 0 0 1),146), ((exponent 0 1 1 0 0 0),25), ((exponent 0 1 0 1 0 0),34), ((exponent 0 1 0 0 1 0),240), ((exponent 0 1 0 0 0 1),182), ((exponent 0 0 1 1 0 0),37), ((exponent 0 0 1 0 1 0),81), ((exponent 0 0 1 0 0 1),55), ((exponent 0 0 0 1 1 0),39), ((exponent 0 0 0 1 0 1),99), ((exponent 0 0 0 0 1 1),120), ((exponent 1 0 0 0 0 0),24), ((exponent 0 1 0 0 0 0),43), ((exponent 0 0 1 0 0 0),159), ((exponent 0 0 0 1 0 0),121), ((exponent 0 0 0 0 1 0),58), ((exponent 0 0 0 0 0 1),112), ((exponent 0 0 0 0 0 0),209)]
def claim119 : Sparse 6 K := [((exponent 1 1 0 0 0 2),1), ((exponent 0 1 1 1 0 1),120), ((exponent 0 1 1 0 1 1),72), ((exponent 0 1 0 1 1 1),22), ((exponent 1 0 1 0 0 1),187), ((exponent 1 0 0 0 0 2),217), ((exponent 0 1 1 0 0 1),50), ((exponent 0 1 0 1 0 1),135), ((exponent 0 1 0 0 1 1),161), ((exponent 0 1 0 0 0 2),170), ((exponent 0 0 1 1 0 1),245), ((exponent 0 0 1 0 1 1),189), ((exponent 0 0 0 1 1 1),111), ((exponent 1 0 0 0 0 1),64), ((exponent 0 1 0 0 0 1),32), ((exponent 0 0 1 0 0 1),242), ((exponent 0 0 0 1 0 1),23), ((exponent 0 0 0 0 1 1),34), ((exponent 0 0 0 0 0 2),105), ((exponent 0 0 0 0 0 1),237)]
def claim120 : Sparse 6 K := [((exponent 0 1 0 1 1 1),1), ((exponent 0 0 1 1 1 0),57), ((exponent 0 0 1 1 0 1),152), ((exponent 0 0 1 0 1 1),194), ((exponent 0 0 0 1 1 1),161), ((exponent 1 0 0 0 0 1),103), ((exponent 0 1 0 1 0 0),200), ((exponent 0 1 0 0 1 0),105), ((exponent 0 1 0 0 0 1),163), ((exponent 0 0 1 1 0 0),71), ((exponent 0 0 1 0 1 0),89), ((exponent 0 0 1 0 0 1),24), ((exponent 0 0 0 1 1 0),39), ((exponent 0 0 0 1 0 1),215), ((exponent 0 0 0 0 1 1),218), ((exponent 1 0 0 0 0 0),153), ((exponent 0 1 0 0 0 0),130), ((exponent 0 0 1 0 0 0),94), ((exponent 0 0 0 1 0 0),115), ((exponent 0 0 0 0 1 0),239), ((exponent 0 0 0 0 0 1),210), ((exponent 0 0 0 0 0 0),80)]
def claim121 : Sparse 6 K := [((exponent 1 1 1 1 0 0),1), ((exponent 1 0 1 1 1 0),1), ((exponent 0 1 1 2 0 0),139), ((exponent 0 1 1 1 1 0),237), ((exponent 1 0 1 1 0 0),95), ((exponent 0 1 1 1 0 0),55), ((exponent 0 0 1 2 0 0),121), ((exponent 0 0 1 1 1 0),4), ((exponent 0 0 1 1 0 0),36), ((exponent 0 0 0 2 0 0),35), ((exponent 0 0 0 1 1 0),21), ((exponent 0 0 0 1 0 0),36)]
def claim122 : Sparse 6 K := [((exponent 0 1 1 0 1 0),1), ((exponent 0 1 1 0 0 1),89), ((exponent 0 0 1 1 1 0),38), ((exponent 1 0 1 0 0 0),103), ((exponent 1 0 0 0 0 1),161), ((exponent 0 1 1 0 0 0),215), ((exponent 0 1 0 1 0 0),139), ((exponent 0 1 0 0 1 0),222), ((exponent 0 1 0 0 0 1),222), ((exponent 0 0 1 1 0 0),234), ((exponent 0 0 1 0 1 0),176), ((exponent 0 0 1 0 0 1),132), ((exponent 0 0 0 1 1 0),190), ((exponent 1 0 0 0 0 0),74), ((exponent 0 1 0 0 0 0),93), ((exponent 0 0 1 0 0 0),202), ((exponent 0 0 0 1 0 0),206), ((exponent 0 0 0 0 1 0),52), ((exponent 0 0 0 0 0 1),129), ((exponent 0 0 0 0 0 0),158)]
def claim123 : Sparse 6 K := [((exponent 1 1 1 0 1 0),1), ((exponent 1 0 1 0 2 0),1), ((exponent 0 1 1 1 1 0),139), ((exponent 0 1 1 0 2 0),237), ((exponent 1 0 1 0 1 0),95), ((exponent 0 1 1 0 1 0),55), ((exponent 0 0 1 1 1 0),121), ((exponent 0 0 1 0 2 0),4), ((exponent 0 0 1 0 1 0),36), ((exponent 0 0 0 1 1 0),35), ((exponent 0 0 0 0 2 0),21), ((exponent 0 0 0 0 1 0),36)]
def claim124 : Sparse 6 K := [((exponent 0 1 1 0 0 1),1), ((exponent 0 0 1 1 1 0),75), ((exponent 1 0 1 0 0 0),160), ((exponent 1 0 0 0 0 1),160), ((exponent 0 1 1 0 0 0),193), ((exponent 0 1 0 1 0 0),163), ((exponent 0 1 0 0 1 0),25), ((exponent 0 1 0 0 0 1),213), ((exponent 0 0 1 1 0 0),62), ((exponent 0 0 1 0 1 0),155), ((exponent 0 0 1 0 0 1),231), ((exponent 0 0 0 1 1 0),98), ((exponent 1 0 0 0 0 0),123), ((exponent 0 1 0 0 0 0),104), ((exponent 0 0 1 0 0 0),253), ((exponent 0 0 0 1 0 0),71), ((exponent 0 0 0 0 1 0),242), ((exponent 0 0 0 0 0 1),150), ((exponent 0 0 0 0 0 0),215)]
def claim125 : Sparse 6 K := [((exponent 1 1 1 0 0 1),1), ((exponent 1 0 1 0 1 1),1), ((exponent 0 1 1 1 0 1),139), ((exponent 0 1 1 0 1 1),237), ((exponent 1 0 1 0 0 1),95), ((exponent 0 1 1 0 0 1),55), ((exponent 0 0 1 1 0 1),121), ((exponent 0 0 1 0 1 1),4), ((exponent 0 0 1 0 0 1),36), ((exponent 0 0 0 1 0 1),35), ((exponent 0 0 0 0 1 1),21), ((exponent 0 0 0 0 0 1),36)]
def claim126 : Sparse 6 K := [((exponent 1 0 1 0 1 1),1), ((exponent 0 1 1 0 1 1),205), ((exponent 0 0 1 1 1 0),131), ((exponent 0 0 1 1 0 1),5), ((exponent 0 0 1 0 1 1),207), ((exponent 0 0 0 1 1 1),34), ((exponent 1 0 1 0 0 0),100), ((exponent 1 0 0 0 0 1),169), ((exponent 0 1 1 0 0 0),7), ((exponent 0 1 0 1 0 0),221), ((exponent 0 1 0 0 1 0),65), ((exponent 0 1 0 0 0 1),239), ((exponent 0 0 1 1 0 0),140), ((exponent 0 0 1 0 1 0),113), ((exponent 0 0 1 0 0 1),159), ((exponent 0 0 0 1 1 0),215), ((exponent 0 0 0 1 0 1),152), ((exponent 0 0 0 0 1 1),46), ((exponent 1 0 0 0 0 0),252), ((exponent 0 1 0 0 0 0),36), ((exponent 0 0 1 0 0 0),30), ((exponent 0 0 0 1 0 0),112), ((exponent 0 0 0 0 1 0),216), ((exponent 0 0 0 0 0 1),172), ((exponent 0 0 0 0 0 0),171)]
def claim127 : Sparse 6 K := [((exponent 0 1 1 1 1 0),1), ((exponent 0 1 1 0 2 0),52), ((exponent 0 1 0 1 2 0),70), ((exponent 0 0 1 1 2 0),1), ((exponent 1 0 1 0 1 0),235), ((exponent 1 0 0 0 1 1),138), ((exponent 0 1 1 0 1 0),236), ((exponent 0 1 0 1 1 0),19), ((exponent 0 1 0 0 2 0),208), ((exponent 0 1 0 0 1 1),237), ((exponent 0 0 1 1 1 0),118), ((exponent 0 0 1 0 2 0),21), ((exponent 0 0 0 1 2 0),245), ((exponent 1 0 0 0 1 0),11), ((exponent 0 1 0 0 1 0),227), ((exponent 0 0 1 0 1 0),235), ((exponent 0 0 0 1 1 0),100), ((exponent 0 0 0 0 2 0),3), ((exponent 0 0 0 0 1 1),186), ((exponent 0 0 0 0 1 0),241)]
def claim128 : Sparse 6 K := [((exponent 0 0 1 1 1 0),1), ((exponent 0 0 1 1 0 1),166), ((exponent 0 0 1 0 1 1),147), ((exponent 0 0 0 1 1 1),188), ((exponent 1 0 0 0 0 1),83), ((exponent 0 1 0 1 0 0),159), ((exponent 0 1 0 0 1 0),253), ((exponent 0 1 0 0 0 1),47), ((exponent 0 0 1 1 0 0),13), ((exponent 0 0 1 0 1 0),128), ((exponent 0 0 1 0 0 1),97), ((exponent 0 0 0 1 1 0),146), ((exponent 0 0 0 1 0 1),209), ((exponent 0 0 0 0 1 1),55), ((exponent 1 0 0 0 0 0),90), ((exponent 0 1 0 0 0 0),13), ((exponent 0 0 1 0 0 0),54), ((exponent 0 0 0 1 0 0),244), ((exponent 0 0 0 0 1 0),189), ((exponent 0 0 0 0 0 1),115), ((exponent 0 0 0 0 0 0),146)]
def claim129 : Sparse 6 K := [((exponent 0 1 1 1 0 1),1), ((exponent 0 1 1 0 1 1),52), ((exponent 0 1 0 1 1 1),70), ((exponent 0 0 1 1 1 1),1), ((exponent 1 0 1 0 0 1),235), ((exponent 1 0 0 0 0 2),138), ((exponent 0 1 1 0 0 1),236), ((exponent 0 1 0 1 0 1),19), ((exponent 0 1 0 0 1 1),208), ((exponent 0 1 0 0 0 2),237), ((exponent 0 0 1 1 0 1),118), ((exponent 0 0 1 0 1 1),21), ((exponent 0 0 0 1 1 1),245), ((exponent 1 0 0 0 0 1),11), ((exponent 0 1 0 0 0 1),227), ((exponent 0 0 1 0 0 1),235), ((exponent 0 0 0 1 0 1),100), ((exponent 0 0 0 0 1 1),3), ((exponent 0 0 0 0 0 2),186), ((exponent 0 0 0 0 0 1),241)]
def claim130 : Sparse 6 K := [((exponent 0 0 1 1 1 1),1), ((exponent 0 0 1 1 0 1),186), ((exponent 0 0 1 0 1 1),163), ((exponent 0 0 0 1 1 1),154), ((exponent 1 0 0 0 0 1),184), ((exponent 0 1 0 1 0 0),79), ((exponent 0 1 0 0 1 0),125), ((exponent 0 1 0 0 0 1),99), ((exponent 0 0 1 1 0 0),224), ((exponent 0 0 1 0 1 0),98), ((exponent 0 0 1 0 0 1),172), ((exponent 0 0 0 1 1 0),29), ((exponent 0 0 0 1 0 1),175), ((exponent 0 0 0 0 1 1),119), ((exponent 1 0 0 0 0 0),212), ((exponent 0 1 0 0 0 0),107), ((exponent 0 0 1 0 0 0),157), ((exponent 0 0 0 1 0 0),38), ((exponent 0 0 0 0 1 0),249), ((exponent 0 0 0 0 0 1),218), ((exponent 0 0 0 0 0 0),15)]
def claim131 : Sparse 6 K := [((exponent 0 1 0 1 1 1),1), ((exponent 0 0 1 1 1 1),166), ((exponent 0 1 0 1 0 1),190), ((exponent 0 1 0 0 1 1),178), ((exponent 0 1 0 0 0 2),145), ((exponent 0 0 1 1 0 1),101), ((exponent 0 0 1 0 1 1),79), ((exponent 0 0 1 0 0 2),204), ((exponent 0 0 0 1 1 1),154), ((exponent 0 1 0 0 0 1),234), ((exponent 0 0 1 0 0 1),96), ((exponent 0 0 0 1 0 1),101), ((exponent 0 0 0 0 1 1),138), ((exponent 0 0 0 0 0 2),96), ((exponent 0 0 0 0 0 1),66)]
def claim132 : Sparse 6 K := [((exponent 0 0 1 1 0 1),1), ((exponent 0 0 1 0 1 1),241), ((exponent 0 0 0 1 1 1),18), ((exponent 1 0 0 0 0 1),182), ((exponent 0 1 0 1 0 0),230), ((exponent 0 1 0 0 1 0),110), ((exponent 0 1 0 0 0 1),80), ((exponent 0 0 1 1 0 0),243), ((exponent 0 0 1 0 1 0),219), ((exponent 0 0 1 0 0 1),150), ((exponent 0 0 0 1 1 0),39), ((exponent 0 0 0 1 0 1),63), ((exponent 0 0 0 0 1 1),231), ((exponent 1 0 0 0 0 0),247), ((exponent 0 1 0 0 0 0),141), ((exponent 0 0 1 0 0 0),86), ((exponent 0 0 0 1 0 0),195), ((exponent 0 0 0 0 1 0),217), ((exponent 0 0 0 0 0 1),139), ((exponent 0 0 0 0 0 0),62)]
def claim133 : Sparse 6 K := [((exponent 1 0 1 1 1 0),1), ((exponent 1 0 1 1 0 1),89), ((exponent 0 1 1 1 1 0),140), ((exponent 0 1 1 1 0 1),2), ((exponent 0 0 1 2 1 0),250), ((exponent 1 0 1 1 0 0),136), ((exponent 1 0 0 1 0 1),97), ((exponent 0 1 1 1 0 0),137), ((exponent 0 1 0 2 0 0),90), ((exponent 0 1 0 1 1 0),249), ((exponent 0 1 0 1 0 1),14), ((exponent 0 0 1 2 0 0),47), ((exponent 0 0 1 1 1 0),175), ((exponent 0 0 1 1 0 1),5), ((exponent 0 0 0 2 1 0),61), ((exponent 1 0 0 1 0 0),207), ((exponent 0 1 0 1 0 0),131), ((exponent 0 0 1 1 0 0),156), ((exponent 0 0 0 2 0 0),211), ((exponent 0 0 0 1 1 0),112), ((exponent 0 0 0 1 0 1),145), ((exponent 0 0 0 1 0 0),209)]
def claim134 : Sparse 6 K := [((exponent 0 1 1 0 1 1),1), ((exponent 0 0 1 0 1 1),217), ((exponent 0 0 0 1 1 1),167), ((exponent 1 0 1 0 0 0),156), ((exponent 1 0 0 0 0 1),184), ((exponent 0 1 1 0 0 0),42), ((exponent 0 1 0 1 0 0),153), ((exponent 0 1 0 0 1 0),96), ((exponent 0 1 0 0 0 1),140), ((exponent 0 0 1 1 0 0),200), ((exponent 0 0 1 0 1 0),33), ((exponent 0 0 1 0 0 1),220), ((exponent 0 0 0 1 1 0),208), ((exponent 0 0 0 1 0 1),11), ((exponent 0 0 0 0 1 1),148), ((exponent 1 0 0 0 0 0),209), ((exponent 0 1 0 0 0 0),18), ((exponent 0 0 1 0 0 0),166), ((exponent 0 0 0 1 0 0),101), ((exponent 0 0 0 0 1 0),230), ((exponent 0 0 0 0 0 1),188), ((exponent 0 0 0 0 0 0),156)]
def claim135 : Sparse 6 K := [((exponent 1 0 1 0 1 1),1), ((exponent 1 0 1 0 0 2),89), ((exponent 0 1 1 0 1 1),140), ((exponent 0 1 1 0 0 2),2), ((exponent 0 0 1 1 1 1),250), ((exponent 1 0 1 0 0 1),136), ((exponent 1 0 0 0 0 2),97), ((exponent 0 1 1 0 0 1),137), ((exponent 0 1 0 1 0 1),90), ((exponent 0 1 0 0 1 1),249), ((exponent 0 1 0 0 0 2),14), ((exponent 0 0 1 1 0 1),47), ((exponent 0 0 1 0 1 1),175), ((exponent 0 0 1 0 0 2),5), ((exponent 0 0 0 1 1 1),61), ((exponent 1 0 0 0 0 1),207), ((exponent 0 1 0 0 0 1),131), ((exponent 0 0 1 0 0 1),156), ((exponent 0 0 0 1 0 1),211), ((exponent 0 0 0 0 1 1),112), ((exponent 0 0 0 0 0 2),145), ((exponent 0 0 0 0 0 1),209)]
def claim136 : Sparse 6 K := [((exponent 0 0 1 0 1 1),1), ((exponent 0 0 0 1 1 1),31), ((exponent 1 0 1 0 0 0),202), ((exponent 1 0 0 0 0 1),231), ((exponent 0 1 1 0 0 0),126), ((exponent 0 1 0 1 0 0),140), ((exponent 0 1 0 0 1 0),154), ((exponent 0 1 0 0 0 1),190), ((exponent 0 0 1 1 0 0),18), ((exponent 0 0 1 0 1 0),13), ((exponent 0 0 1 0 0 1),71), ((exponent 0 0 0 1 1 0),69), ((exponent 0 0 0 1 0 1),43), ((exponent 0 0 0 0 1 1),238), ((exponent 1 0 0 0 0 0),140), ((exponent 0 1 0 0 0 0),24), ((exponent 0 0 1 0 0 0),101), ((exponent 0 0 0 1 0 0),176), ((exponent 0 0 0 0 1 0),115), ((exponent 0 0 0 0 0 1),70), ((exponent 0 0 0 0 0 0),116)]
def claim137 : Sparse 6 K := [((exponent 0 1 0 1 0 2),1), ((exponent 0 1 0 0 1 2),140), ((exponent 0 0 1 1 1 1),7), ((exponent 0 0 1 1 0 2),89), ((exponent 0 0 1 0 1 2),2), ((exponent 0 1 0 1 0 1),126), ((exponent 0 1 0 0 1 1),14), ((exponent 0 1 0 0 0 2),171), ((exponent 0 0 1 1 0 1),183), ((exponent 0 0 1 0 1 1),244), ((exponent 0 0 1 0 0 2),182), ((exponent 0 0 0 1 1 1),38), ((exponent 0 0 0 1 0 2),136), ((exponent 0 0 0 0 1 2),137), ((exponent 0 1 0 0 0 1),154), ((exponent 0 0 1 0 0 1),50), ((exponent 0 0 0 1 0 1),249), ((exponent 0 0 0 0 1 1),82), ((exponent 0 0 0 0 0 2),27), ((exponent 0 0 0 0 0 1),195)]
def claim138 : Sparse 6 K := [((exponent 0 0 0 1 1 1),1), ((exponent 1 0 0 0 0 1),161), ((exponent 0 1 0 1 0 0),177), ((exponent 0 1 0 0 1 0),42), ((exponent 0 1 0 0 0 1),209), ((exponent 0 0 1 1 0 0),108), ((exponent 0 0 1 0 1 0),233), ((exponent 0 0 1 0 0 1),103), ((exponent 0 0 0 1 1 0),59), ((exponent 0 0 0 1 0 1),210), ((exponent 0 0 0 0 1 1),105), ((exponent 1 0 0 0 0 0),156), ((exponent 0 1 0 0 0 0),242), ((exponent 0 0 1 0 0 0),33), ((exponent 0 0 0 1 0 0),87), ((exponent 0 0 0 0 1 0),226), ((exponent 0 0 0 0 0 1),239), ((exponent 0 0 0 0 0 0),39)]
def claim139 : Sparse 6 K := [((exponent 1 0 2 0 0 1),1), ((exponent 0 1 2 0 1 0),207), ((exponent 0 1 2 0 0 1),36), ((exponent 0 0 2 1 1 0),108), ((exponent 1 0 2 0 0 0),77), ((exponent 1 0 1 0 0 1),68), ((exponent 0 1 2 0 0 0),153), ((exponent 0 1 1 1 0 0),181), ((exponent 0 1 1 0 1 0),174), ((exponent 0 1 1 0 0 1),158), ((exponent 0 0 2 1 0 0),80), ((exponent 0 0 2 0 1 0),61), ((exponent 0 0 2 0 0 1),60), ((exponent 0 0 1 1 1 0),218), ((exponent 1 0 1 0 0 0),217), ((exponent 0 1 1 0 0 0),16), ((exponent 0 0 2 0 0 0),188), ((exponent 0 0 1 1 0 0),151), ((exponent 0 0 1 0 1 0),255), ((exponent 0 0 1 0 0 1),14), ((exponent 0 0 1 0 0 0),215)]
def claim140 : Sparse 6 K := [((exponent 1 0 1 0 0 0),1), ((exponent 1 0 0 0 0 1),237), ((exponent 0 1 1 0 0 0),114), ((exponent 0 1 0 1 0 0),177), ((exponent 0 1 0 0 1 0),92), ((exponent 0 1 0 0 0 1),81), ((exponent 0 0 1 1 0 0),107), ((exponent 0 0 1 0 1 0),205), ((exponent 0 0 1 0 0 1),92), ((exponent 0 0 0 1 1 0),19), ((exponent 0 0 0 1 0 1),249), ((exponent 0 0 0 0 1 1),23), ((exponent 1 0 0 0 0 0),33), ((exponent 0 1 0 0 0 0),127), ((exponent 0 0 1 0 0 0),7), ((exponent 0 0 0 1 0 0),254), ((exponent 0 0 0 0 1 0),3), ((exponent 0 0 0 0 0 1),102), ((exponent 0 0 0 0 0 0),82)]
def claim141 : Sparse 6 K := [((exponent 1 0 1 1 0 1),1), ((exponent 0 1 1 1 1 0),207), ((exponent 0 1 1 1 0 1),36), ((exponent 0 0 1 2 1 0),108), ((exponent 1 0 1 1 0 0),77), ((exponent 1 0 0 1 0 1),68), ((exponent 0 1 1 1 0 0),153), ((exponent 0 1 0 2 0 0),181), ((exponent 0 1 0 1 1 0),174), ((exponent 0 1 0 1 0 1),158), ((exponent 0 0 1 2 0 0),80), ((exponent 0 0 1 1 1 0),61), ((exponent 0 0 1 1 0 1),60), ((exponent 0 0 0 2 1 0),218), ((exponent 1 0 0 1 0 0),217), ((exponent 0 1 0 1 0 0),16), ((exponent 0 0 1 1 0 0),188), ((exponent 0 0 0 2 0 0),151), ((exponent 0 0 0 1 1 0),255), ((exponent 0 0 0 1 0 1),14), ((exponent 0 0 0 1 0 0),215)]
def claim142 : Sparse 6 K := [((exponent 1 0 0 0 0 1),1), ((exponent 0 1 1 0 0 0),4), ((exponent 0 1 0 1 0 0),46), ((exponent 0 1 0 0 1 0),84), ((exponent 0 1 0 0 0 1),144), ((exponent 0 0 1 1 0 0),233), ((exponent 0 0 1 0 1 0),150), ((exponent 0 0 1 0 0 1),244), ((exponent 0 0 0 1 1 0),89), ((exponent 0 0 0 1 0 1),181), ((exponent 0 0 0 0 1 1),235), ((exponent 1 0 0 0 0 0),250), ((exponent 0 1 0 0 0 0),156), ((exponent 0 0 1 0 0 0),81), ((exponent 0 0 0 1 0 0),233), ((exponent 0 0 0 0 1 0),234), ((exponent 0 0 0 0 0 1),253), ((exponent 0 0 0 0 0 0),107)]
def claim143 : Sparse 6 K := [((exponent 1 0 1 0 0 2),1), ((exponent 0 1 1 0 1 1),207), ((exponent 0 1 1 0 0 2),36), ((exponent 0 0 1 1 1 1),108), ((exponent 1 0 1 0 0 1),77), ((exponent 1 0 0 0 0 2),68), ((exponent 0 1 1 0 0 1),153), ((exponent 0 1 0 1 0 1),181), ((exponent 0 1 0 0 1 1),174), ((exponent 0 1 0 0 0 2),158), ((exponent 0 0 1 1 0 1),80), ((exponent 0 0 1 0 1 1),61), ((exponent 0 0 1 0 0 2),60), ((exponent 0 0 0 1 1 1),218), ((exponent 1 0 0 0 0 1),217), ((exponent 0 1 0 0 0 1),16), ((exponent 0 0 1 0 0 1),188), ((exponent 0 0 0 1 0 1),151), ((exponent 0 0 0 0 1 1),255), ((exponent 0 0 0 0 0 2),14), ((exponent 0 0 0 0 0 1),215)]
def claim144 : Sparse 6 K := [((exponent 0 1 1 0 0 0),1), ((exponent 0 1 0 1 0 0),162), ((exponent 0 1 0 0 1 0),32), ((exponent 0 1 0 0 0 1),122), ((exponent 0 0 1 1 0 0),163), ((exponent 0 0 1 0 1 0),198), ((exponent 0 0 1 0 0 1),65), ((exponent 0 0 0 1 1 0),52), ((exponent 0 0 0 1 0 1),137), ((exponent 0 0 0 0 1 1),138), ((exponent 1 0 0 0 0 0),150), ((exponent 0 1 0 0 0 0),232), ((exponent 0 0 1 0 0 0),36), ((exponent 0 0 0 1 0 0),161), ((exponent 0 0 0 0 1 0),99), ((exponent 0 0 0 0 0 1),217), ((exponent 0 0 0 0 0 0),67)]
def claim145 : Sparse 6 K := [((exponent 1 1 1 0 1 0),1), ((exponent 1 1 1 0 0 1),89), ((exponent 1 0 1 1 1 0),38), ((exponent 2 0 1 0 0 0),103), ((exponent 2 0 0 0 0 1),161), ((exponent 1 1 1 0 0 0),215), ((exponent 1 1 0 1 0 0),139), ((exponent 1 1 0 0 1 0),222), ((exponent 1 1 0 0 0 1),222), ((exponent 1 0 1 1 0 0),234), ((exponent 1 0 1 0 1 0),176), ((exponent 1 0 1 0 0 1),132), ((exponent 1 0 0 1 1 0),190), ((exponent 2 0 0 0 0 0),74), ((exponent 1 1 0 0 0 0),93), ((exponent 1 0 1 0 0 0),202), ((exponent 1 0 0 1 0 0),206), ((exponent 1 0 0 0 1 0),52), ((exponent 1 0 0 0 0 1),129), ((exponent 1 0 0 0 0 0),158)]
def claim146 : Sparse 6 K := [((exponent 0 1 0 1 0 0),1), ((exponent 0 1 0 0 1 0),177), ((exponent 0 1 0 0 0 1),104), ((exponent 0 0 1 1 0 0),147), ((exponent 0 0 1 0 1 0),81), ((exponent 0 0 1 0 0 1),26), ((exponent 0 0 0 1 1 0),115), ((exponent 0 0 0 1 0 1),117), ((exponent 0 0 0 0 1 1),192), ((exponent 1 0 0 0 0 0),29), ((exponent 0 1 0 0 0 0),146), ((exponent 0 0 1 0 0 0),162), ((exponent 0 0 0 1 0 0),91), ((exponent 0 0 0 0 1 0),189), ((exponent 0 0 0 0 0 1),47), ((exponent 0 0 0 0 0 0),144)]
def claim147 : Sparse 6 K := [((exponent 0 1 1 0 2 0),1), ((exponent 0 1 1 0 1 1),89), ((exponent 0 0 1 1 2 0),38), ((exponent 1 0 1 0 1 0),103), ((exponent 1 0 0 0 1 1),161), ((exponent 0 1 1 0 1 0),215), ((exponent 0 1 0 1 1 0),139), ((exponent 0 1 0 0 2 0),222), ((exponent 0 1 0 0 1 1),222), ((exponent 0 0 1 1 1 0),234), ((exponent 0 0 1 0 2 0),176), ((exponent 0 0 1 0 1 1),132), ((exponent 0 0 0 1 2 0),190), ((exponent 1 0 0 0 1 0),74), ((exponent 0 1 0 0 1 0),93), ((exponent 0 0 1 0 1 0),202), ((exponent 0 0 0 1 1 0),206), ((exponent 0 0 0 0 2 0),52), ((exponent 0 0 0 0 1 1),129), ((exponent 0 0 0 0 1 0),158)]
def claim148 : Sparse 6 K := [((exponent 0 1 0 0 1 0),1), ((exponent 0 1 0 0 0 1),210), ((exponent 0 0 1 1 0 0),54), ((exponent 0 0 1 0 1 0),208), ((exponent 0 0 1 0 0 1),125), ((exponent 0 0 0 1 1 0),36), ((exponent 0 0 0 1 0 1),204), ((exponent 0 0 0 0 1 1),227), ((exponent 1 0 0 0 0 0),249), ((exponent 0 1 0 0 0 0),130), ((exponent 0 0 1 0 0 0),256), ((exponent 0 0 0 1 0 0),173), ((exponent 0 0 0 0 1 0),81), ((exponent 0 0 0 0 0 1),70), ((exponent 0 0 0 0 0 0),130)]
def claim149 : Sparse 6 K := [((exponent 0 1 1 0 1 1),1), ((exponent 0 1 1 0 0 2),89), ((exponent 0 0 1 1 1 1),38), ((exponent 1 0 1 0 0 1),103), ((exponent 1 0 0 0 0 2),161), ((exponent 0 1 1 0 0 1),215), ((exponent 0 1 0 1 0 1),139), ((exponent 0 1 0 0 1 1),222), ((exponent 0 1 0 0 0 2),222), ((exponent 0 0 1 1 0 1),234), ((exponent 0 0 1 0 1 1),176), ((exponent 0 0 1 0 0 2),132), ((exponent 0 0 0 1 1 1),190), ((exponent 1 0 0 0 0 1),74), ((exponent 0 1 0 0 0 1),93), ((exponent 0 0 1 0 0 1),202), ((exponent 0 0 0 1 0 1),206), ((exponent 0 0 0 0 1 1),52), ((exponent 0 0 0 0 0 2),129), ((exponent 0 0 0 0 0 1),158)]
def claim150 : Sparse 6 K := [((exponent 0 1 0 0 0 1),1), ((exponent 0 0 1 1 0 0),117), ((exponent 0 0 1 0 1 0),163), ((exponent 0 0 1 0 0 1),30), ((exponent 0 0 0 1 1 0),21), ((exponent 0 0 0 1 0 1),4), ((exponent 0 0 0 0 1 1),220), ((exponent 1 0 0 0 0 0),4), ((exponent 0 1 0 0 0 0),108), ((exponent 0 0 1 0 0 0),51), ((exponent 0 0 0 1 0 0),188), ((exponent 0 0 0 0 1 0),60), ((exponent 0 0 0 0 0 1),142), ((exponent 0 0 0 0 0 0),93)]
def claim151 : Sparse 6 K := [((exponent 0 1 1 0 0 2),1), ((exponent 0 0 1 1 1 1),75), ((exponent 1 0 1 0 0 1),160), ((exponent 1 0 0 0 0 2),160), ((exponent 0 1 1 0 0 1),193), ((exponent 0 1 0 1 0 1),163), ((exponent 0 1 0 0 1 1),25), ((exponent 0 1 0 0 0 2),213), ((exponent 0 0 1 1 0 1),62), ((exponent 0 0 1 0 1 1),155), ((exponent 0 0 1 0 0 2),231), ((exponent 0 0 0 1 1 1),98), ((exponent 1 0 0 0 0 1),123), ((exponent 0 1 0 0 0 1),104), ((exponent 0 0 1 0 0 1),253), ((exponent 0 0 0 1 0 1),71), ((exponent 0 0 0 0 1 1),242), ((exponent 0 0 0 0 0 2),150), ((exponent 0 0 0 0 0 1),215)]
def claim152 : Sparse 6 K := [((exponent 0 0 1 1 0 0),1), ((exponent 0 0 1 0 1 0),184), ((exponent 0 0 1 0 0 1),33), ((exponent 0 0 0 1 1 0),81), ((exponent 0 0 0 1 0 1),99), ((exponent 0 0 0 0 1 1),22), ((exponent 1 0 0 0 0 0),84), ((exponent 0 1 0 0 0 0),94), ((exponent 0 0 1 0 0 0),34), ((exponent 0 0 0 1 0 0),1), ((exponent 0 0 0 0 1 0),161), ((exponent 0 0 0 0 0 1),184), ((exponent 0 0 0 0 0 0),172)]
def claim153 : Sparse 6 K := [((exponent 0 0 1 1 1 1),1), ((exponent 0 0 1 1 0 2),166), ((exponent 0 0 1 0 1 2),147), ((exponent 0 0 0 1 1 2),188), ((exponent 1 0 0 0 0 2),83), ((exponent 0 1 0 1 0 1),159), ((exponent 0 1 0 0 1 1),253), ((exponent 0 1 0 0 0 2),47), ((exponent 0 0 1 1 0 1),13), ((exponent 0 0 1 0 1 1),128), ((exponent 0 0 1 0 0 2),97), ((exponent 0 0 0 1 1 1),146), ((exponent 0 0 0 1 0 2),209), ((exponent 0 0 0 0 1 2),55), ((exponent 1 0 0 0 0 1),90), ((exponent 0 1 0 0 0 1),13), ((exponent 0 0 1 0 0 1),54), ((exponent 0 0 0 1 0 1),244), ((exponent 0 0 0 0 1 1),189), ((exponent 0 0 0 0 0 2),115), ((exponent 0 0 0 0 0 1),146)]
def claim154 : Sparse 6 K := [((exponent 0 0 1 0 1 0),1), ((exponent 0 0 1 0 0 1),203), ((exponent 0 0 0 1 1 0),58), ((exponent 0 0 0 1 0 1),80), ((exponent 0 0 0 0 1 1),175), ((exponent 1 0 0 0 0 0),221), ((exponent 0 1 0 0 0 0),128), ((exponent 0 0 1 0 0 0),116), ((exponent 0 0 0 1 0 0),251), ((exponent 0 0 0 0 1 0),207), ((exponent 0 0 0 0 0 1),96), ((exponent 0 0 0 0 0 0),138)]
def claim155 : Sparse 6 K := [((exponent 0 0 2 1 0 1),1), ((exponent 0 0 2 0 1 1),241), ((exponent 0 0 1 1 1 1),18), ((exponent 1 0 1 0 0 1),182), ((exponent 0 1 1 1 0 0),230), ((exponent 0 1 1 0 1 0),110), ((exponent 0 1 1 0 0 1),80), ((exponent 0 0 2 1 0 0),243), ((exponent 0 0 2 0 1 0),219), ((exponent 0 0 2 0 0 1),150), ((exponent 0 0 1 1 1 0),39), ((exponent 0 0 1 1 0 1),63), ((exponent 0 0 1 0 1 1),231), ((exponent 1 0 1 0 0 0),247), ((exponent 0 1 1 0 0 0),141), ((exponent 0 0 2 0 0 0),86), ((exponent 0 0 1 1 0 0),195), ((exponent 0 0 1 0 1 0),217), ((exponent 0 0 1 0 0 1),139), ((exponent 0 0 1 0 0 0),62)]
def claim156 : Sparse 6 K := [((exponent 0 0 1 0 0 1),1), ((exponent 0 0 0 1 1 0),49), ((exponent 0 0 0 1 0 1),107), ((exponent 0 0 0 0 1 1),114), ((exponent 1 0 0 0 0 0),138), ((exponent 0 1 0 0 0 0),172), ((exponent 0 0 1 0 0 0),5), ((exponent 0 0 0 1 0 0),245), ((exponent 0 0 0 0 1 0),121), ((exponent 0 0 0 0 0 1),44), ((exponent 0 0 0 0 0 0),128)]
def claim157 : Sparse 6 K := [((exponent 0 0 1 2 0 1),1), ((exponent 0 0 1 1 1 1),241), ((exponent 0 0 0 2 1 1),18), ((exponent 1 0 0 1 0 1),182), ((exponent 0 1 0 2 0 0),230), ((exponent 0 1 0 1 1 0),110), ((exponent 0 1 0 1 0 1),80), ((exponent 0 0 1 2 0 0),243), ((exponent 0 0 1 1 1 0),219), ((exponent 0 0 1 1 0 1),150), ((exponent 0 0 0 2 1 0),39), ((exponent 0 0 0 2 0 1),63), ((exponent 0 0 0 1 1 1),231), ((exponent 1 0 0 1 0 0),247), ((exponent 0 1 0 1 0 0),141), ((exponent 0 0 1 1 0 0),86), ((exponent 0 0 0 2 0 0),195), ((exponent 0 0 0 1 1 0),217), ((exponent 0 0 0 1 0 1),139), ((exponent 0 0 0 1 0 0),62)]
def claim158 : Sparse 6 K := [((exponent 0 0 0 1 1 0),1), ((exponent 0 0 0 1 0 1),34), ((exponent 0 0 0 0 1 1),177), ((exponent 1 0 0 0 0 0),219), ((exponent 0 1 0 0 0 0),248), ((exponent 0 0 1 0 0 0),236), ((exponent 0 0 0 0 1 0),83), ((exponent 0 0 0 0 0 1),163), ((exponent 0 0 0 0 0 0),228)]
def claim159 : Sparse 6 K := [((exponent 0 0 2 0 1 1),1), ((exponent 0 0 1 1 1 1),31), ((exponent 1 0 2 0 0 0),202), ((exponent 1 0 1 0 0 1),231), ((exponent 0 1 2 0 0 0),126), ((exponent 0 1 1 1 0 0),140), ((exponent 0 1 1 0 1 0),154), ((exponent 0 1 1 0 0 1),190), ((exponent 0 0 2 1 0 0),18), ((exponent 0 0 2 0 1 0),13), ((exponent 0 0 2 0 0 1),71), ((exponent 0 0 1 1 1 0),69), ((exponent 0 0 1 1 0 1),43), ((exponent 0 0 1 0 1 1),238), ((exponent 1 0 1 0 0 0),140), ((exponent 0 1 1 0 0 0),24), ((exponent 0 0 2 0 0 0),101), ((exponent 0 0 1 1 0 0),176), ((exponent 0 0 1 0 1 0),115), ((exponent 0 0 1 0 0 1),70), ((exponent 0 0 1 0 0 0),116)]
def claim160 : Sparse 6 K := [((exponent 0 0 0 1 0 1),1), ((exponent 0 0 0 0 1 1),255), ((exponent 1 0 0 0 0 0),40), ((exponent 0 1 0 0 0 0),222), ((exponent 0 0 1 0 0 0),104), ((exponent 0 0 0 1 0 0),4), ((exponent 0 0 0 0 1 0),121), ((exponent 0 0 0 0 0 1),124), ((exponent 0 0 0 0 0 0),150)]
def claim161 : Sparse 6 K := [((exponent 0 0 1 1 1 1),1), ((exponent 0 0 0 2 1 1),31), ((exponent 1 0 1 1 0 0),202), ((exponent 1 0 0 1 0 1),231), ((exponent 0 1 1 1 0 0),126), ((exponent 0 1 0 2 0 0),140), ((exponent 0 1 0 1 1 0),154), ((exponent 0 1 0 1 0 1),190), ((exponent 0 0 1 2 0 0),18), ((exponent 0 0 1 1 1 0),13), ((exponent 0 0 1 1 0 1),71), ((exponent 0 0 0 2 1 0),69), ((exponent 0 0 0 2 0 1),43), ((exponent 0 0 0 1 1 1),238), ((exponent 1 0 0 1 0 0),140), ((exponent 0 1 0 1 0 0),24), ((exponent 0 0 1 1 0 0),101), ((exponent 0 0 0 2 0 0),176), ((exponent 0 0 0 1 1 0),115), ((exponent 0 0 0 1 0 1),70), ((exponent 0 0 0 1 0 0),116)]
def claim162 : Sparse 6 K := [((exponent 0 0 0 0 1 1),1), ((exponent 1 0 0 0 0 0),143), ((exponent 0 1 0 0 0 0),119), ((exponent 0 0 1 0 0 0),86), ((exponent 0 0 0 1 0 0),93), ((exponent 0 0 0 0 1 0),19), ((exponent 0 0 0 0 0 1),170), ((exponent 0 0 0 0 0 0),74)]
def claim163 : Sparse 6 K := [((exponent 1 0 1 1 0 0),1), ((exponent 1 0 0 1 0 1),237), ((exponent 0 1 1 1 0 0),114), ((exponent 0 1 0 2 0 0),177), ((exponent 0 1 0 1 1 0),92), ((exponent 0 1 0 1 0 1),81), ((exponent 0 0 1 2 0 0),107), ((exponent 0 0 1 1 1 0),205), ((exponent 0 0 1 1 0 1),92), ((exponent 0 0 0 2 1 0),19), ((exponent 0 0 0 2 0 1),249), ((exponent 0 0 0 1 1 1),23), ((exponent 1 0 0 1 0 0),33), ((exponent 0 1 0 1 0 0),127), ((exponent 0 0 1 1 0 0),7), ((exponent 0 0 0 2 0 0),254), ((exponent 0 0 0 1 1 0),3), ((exponent 0 0 0 1 0 1),102), ((exponent 0 0 0 1 0 0),82)]
def claim164 : Sparse 6 K := [((exponent 1 0 0 0 0 0),1), ((exponent 0 1 0 0 0 0),203), ((exponent 0 0 1 0 0 0),18), ((exponent 0 0 0 1 0 0),41), ((exponent 0 0 0 0 1 0),109), ((exponent 0 0 0 0 0 1),82), ((exponent 0 0 0 0 0 0),242)]
def claim165 : Sparse 6 K := [((exponent 1 0 1 0 0 1),1), ((exponent 1 0 0 0 0 2),237), ((exponent 0 1 1 0 0 1),114), ((exponent 0 1 0 1 0 1),177), ((exponent 0 1 0 0 1 1),92), ((exponent 0 1 0 0 0 2),81), ((exponent 0 0 1 1 0 1),107), ((exponent 0 0 1 0 1 1),205), ((exponent 0 0 1 0 0 2),92), ((exponent 0 0 0 1 1 1),19), ((exponent 0 0 0 1 0 2),249), ((exponent 0 0 0 0 1 2),23), ((exponent 1 0 0 0 0 1),33), ((exponent 0 1 0 0 0 1),127), ((exponent 0 0 1 0 0 1),7), ((exponent 0 0 0 1 0 1),254), ((exponent 0 0 0 0 1 1),3), ((exponent 0 0 0 0 0 2),102), ((exponent 0 0 0 0 0 1),82)]
def claim166 : Sparse 6 K := [((exponent 0 1 0 0 0 0),1), ((exponent 0 0 1 0 0 0),204), ((exponent 0 0 0 1 0 0),223), ((exponent 0 0 0 0 1 0),211), ((exponent 0 0 0 0 0 1),124), ((exponent 0 0 0 0 0 0),40)]
def claim167 : Sparse 6 K := [((exponent 1 0 0 0 0 2),1), ((exponent 0 1 1 0 0 1),4), ((exponent 0 1 0 1 0 1),46), ((exponent 0 1 0 0 1 1),84), ((exponent 0 1 0 0 0 2),144), ((exponent 0 0 1 1 0 1),233), ((exponent 0 0 1 0 1 1),150), ((exponent 0 0 1 0 0 2),244), ((exponent 0 0 0 1 1 1),89), ((exponent 0 0 0 1 0 2),181), ((exponent 0 0 0 0 1 2),235), ((exponent 1 0 0 0 0 1),250), ((exponent 0 1 0 0 0 1),156), ((exponent 0 0 1 0 0 1),81), ((exponent 0 0 0 1 0 1),233), ((exponent 0 0 0 0 1 1),234), ((exponent 0 0 0 0 0 2),253), ((exponent 0 0 0 0 0 1),107)]
def claim168 : Sparse 6 K := [((exponent 0 0 1 0 0 0),1), ((exponent 0 0 0 1 0 0),235), ((exponent 0 0 0 0 1 0),71), ((exponent 0 0 0 0 0 1),171), ((exponent 0 0 0 0 0 0),80)]
def claim169 : Sparse 6 K := [((exponent 0 1 0 1 0 1),1), ((exponent 0 1 0 0 1 1),177), ((exponent 0 1 0 0 0 2),104), ((exponent 0 0 1 1 0 1),147), ((exponent 0 0 1 0 1 1),81), ((exponent 0 0 1 0 0 2),26), ((exponent 0 0 0 1 1 1),115), ((exponent 0 0 0 1 0 2),117), ((exponent 0 0 0 0 1 2),192), ((exponent 1 0 0 0 0 1),29), ((exponent 0 1 0 0 0 1),146), ((exponent 0 0 1 0 0 1),162), ((exponent 0 0 0 1 0 1),91), ((exponent 0 0 0 0 1 1),189), ((exponent 0 0 0 0 0 2),47), ((exponent 0 0 0 0 0 1),144)]
def claim170 : Sparse 6 K := [((exponent 0 0 0 1 0 0),1), ((exponent 0 0 0 0 1 0),102), ((exponent 0 0 0 0 0 1),230), ((exponent 0 0 0 0 0 0),223)]
def claim171 : Sparse 6 K := [((exponent 0 1 0 0 1 1),1), ((exponent 0 1 0 0 0 2),210), ((exponent 0 0 1 1 0 1),54), ((exponent 0 0 1 0 1 1),208), ((exponent 0 0 1 0 0 2),125), ((exponent 0 0 0 1 1 1),36), ((exponent 0 0 0 1 0 2),204), ((exponent 0 0 0 0 1 2),227), ((exponent 1 0 0 0 0 1),249), ((exponent 0 1 0 0 0 1),130), ((exponent 0 0 1 0 0 1),256), ((exponent 0 0 0 1 0 1),173), ((exponent 0 0 0 0 1 1),81), ((exponent 0 0 0 0 0 2),70), ((exponent 0 0 0 0 0 1),130)]
def claim172 : Sparse 6 K := [((exponent 0 0 0 0 1 0),1), ((exponent 0 0 0 0 0 1),107), ((exponent 0 0 0 0 0 0),98)]
def claim173 : Sparse 6 K := [((exponent 0 0 1 1 0 1),1), ((exponent 0 0 0 2 1 0),49), ((exponent 0 0 0 2 0 1),107), ((exponent 0 0 0 1 1 1),114), ((exponent 1 0 0 1 0 0),138), ((exponent 0 1 0 1 0 0),172), ((exponent 0 0 1 1 0 0),5), ((exponent 0 0 0 2 0 0),245), ((exponent 0 0 0 1 1 0),121), ((exponent 0 0 0 1 0 1),44), ((exponent 0 0 0 1 0 0),128)]
def claim174 : Sparse 6 K := [((exponent 0 0 0 0 0 1),1), ((exponent 0 0 0 0 0 0),45)]
def claim175 : Sparse 6 K := [((exponent 1 0 0 0 0 1),1), ((exponent 0 1 0 0 0 1),203), ((exponent 0 0 1 0 0 1),18), ((exponent 0 0 0 1 0 1),41), ((exponent 0 0 0 0 1 1),109), ((exponent 0 0 0 0 0 2),82), ((exponent 0 0 0 0 0 1),242)]
def claim176 : Sparse 6 K := [((exponent 0 0 0 0 0 0),1)]
def program : List (Entry 6 K) := [
  (.input 0,claim0),
  (.input 1,claim1),
  (.input 2,claim2),
  (.input 3,claim3),
  (.input 4,claim4),
  (.input 5,claim5),
  (.input 6,claim6),
  (.linear [(6,134)],claim7),
  (.mul (variableData 2) 0,claim8),
  (.mul (variableData 3) 0,claim9),
  (.mul (variableData 4) 0,claim10),
  (.mul (variableData 5) 0,claim11),
  (.mul (variableData 0) 1,claim12),
  (.mul (variableData 3) 1,claim13),
  (.mul (variableData 4) 1,claim14),
  (.mul (variableData 0) 2,claim15),
  (.mul (variableData 1) 2,claim16),
  (.mul (variableData 3) 2,claim17),
  (.mul (variableData 4) 2,claim18),
  (.mul (variableData 5) 2,claim19),
  (.mul (variableData 0) 3,claim20),
  (.mul (variableData 1) 3,claim21),
  (.mul (variableData 2) 3,claim22),
  (.mul (variableData 4) 3,claim23),
  (.linear [(4,1), (23,1)],claim24),
  (.mul (variableData 5) 3,claim25),
  (.mul (variableData 0) 4,claim26),
  (.mul (variableData 1) 4,claim27),
  (.mul (variableData 2) 4,claim28),
  (.mul (variableData 3) 4,claim29),
  (.mul (variableData 0) 5,claim30),
  (.mul (variableData 1) 5,claim31),
  (.mul (variableData 2) 5,claim32),
  (.mul (variableData 3) 5,claim33),
  (.mul (variableData 4) 5,claim34),
  (.mul (variableData 0) 7,claim35),
  (.linear [(0,57), (7,57), (9,84), (10,256), (35,173)],claim36),
  (.mul (variableData 1) 7,claim37),
  (.linear [(1,221), (13,222), (14,236), (36,113), (37,144)],claim38),
  (.mul (variableData 2) 7,claim39),
  (.mul (variableData 3) 7,claim40),
  (.linear [(3,57), (7,57), (20,84), (21,256), (40,173)],claim41),
  (.mul (variableData 4) 7,claim42),
  (.linear [(4,70), (26,168), (27,255), (41,23), (42,234)],claim43),
  (.mul (variableData 5) 7,claim44),
  (.mul (variableData 3) 12,claim45),
  (.linear [(39,1), (45,1)],claim46),
  (.mul (variableData 4) 12,claim47),
  (.mul (variableData 5) 15,claim48),
  (.mul (variableData 3) 16,claim49),
  (.mul (variableData 4) 16,claim50),
  (.mul (variableData 5) 16,claim51),
  (.mul (variableData 4) 17,claim52),
  (.mul (variableData 5) 17,claim53),
  (.mul (variableData 5) 18,claim54),
  (.mul (variableData 1) 20,claim55),
  (.mul (variableData 2) 20,claim56),
  (.mul (variableData 2) 21,claim57),
  (.mul (variableData 4) 21,claim58),
  (.linear [(27,1), (58,1)],claim59),
  (.mul (variableData 4) 22,claim60),
  (.linear [(28,1), (60,1)],claim61),
  (.mul (variableData 5) 22,claim62),
  (.mul (variableData 5) 24,claim63),
  (.linear [(5,1), (63,1)],claim64),
  (.mul (variableData 1) 26,claim65),
  (.mul (variableData 2) 26,claim66),
  (.mul (variableData 3) 26,claim67),
  (.linear [(44,1), (67,1)],claim68),
  (.mul (variableData 2) 27,claim69),
  (.mul (variableData 3) 27,claim70),
  (.mul (variableData 3) 28,claim71),
  (.mul (variableData 1) 30,claim72),
  (.linear [(38,1), (43,256), (72,1)],claim73),
  (.mul (variableData 2) 30,claim74),
  (.mul (variableData 2) 31,claim75),
  (.mul (variableData 3) 31,claim76),
  (.mul (variableData 4) 31,claim77),
  (.mul (variableData 3) 32,claim78),
  (.mul (variableData 4) 32,claim79),
  (.mul (variableData 4) 33,claim80),
  (.mul (variableData 1) 36,claim81),
  (.linear [(1,12), (12,60), (13,179), (14,56), (38,162), (43,95), (46,144), (47,35), (81,113)],claim82),
  (.mul (variableData 2) 36,claim83),
  (.mul (variableData 3) 36,claim84),
  (.linear [(3,165), (7,223), (21,83), (36,35), (38,221), (41,235), (43,195), (55,84), (82,63), (84,173)],claim85),
  (.mul (variableData 4) 36,claim86),
  (.linear [(4,169), (26,235), (27,204), (43,223), (65,256), (82,35), (85,84), (86,173)],claim87),
  (.mul (variableData 5) 36,claim88),
  (.linear [(87,21), (88,1)],claim89),
  (.mul (variableData 2) 38,claim90),
  (.linear [(2,221), (17,222), (18,236), (90,1)],claim91),
  (.mul (variableData 3) 38,claim92),
  (.linear [(3,239), (7,228), (21,128), (22,84), (36,208), (38,256), (41,208), (43,24), (82,235), (92,49)],claim93),
  (.mul (variableData 4) 38,claim94),
  (.linear [(4,108), (26,27), (27,231), (28,53), (43,252), (82,27), (93,83), (94,230)],claim95),
  (.mul (variableData 5) 38,claim96),
  (.linear [(87,256), (93,120), (95,103), (96,1)],claim97),
  (.mul (variableData 3) 39,claim98),
  (.linear [(22,57), (39,57), (56,84), (57,256), (93,141), (95,153), (98,173)],claim99),
  (.mul (variableData 4) 39,claim100),
  (.linear [(28,70), (66,168), (69,255), (99,23), (100,234)],claim101),
  (.mul (variableData 5) 39,claim102),
  (.mul (variableData 4) 41,claim103),
  (.linear [(4,172), (26,167), (27,203), (29,89), (43,236), (68,8), (70,159), (95,157), (97,98), (103,249)],claim104),
  (.mul (variableData 5) 43,claim105),
  (.linear [(5,209), (30,39), (31,229), (97,243), (104,60), (105,14)],claim106),
  (.mul (variableData 2) 82,claim107),
  (.linear [(2,24), (15,256), (16,221), (17,102), (18,216), (49,222), (50,236), (93,176), (95,16), (107,1)],claim108),
  (.mul (variableData 3) 82,claim109),
  (.linear [(3,28), (7,169), (21,143), (22,115), (24,242), (39,208), (43,207), (44,246), (57,84), (85,208), (93,118), (95,256), (97,58), (101,235), (104,162), (109,49)],claim110),
  (.mul (variableData 4) 82,claim111),
  (.linear [(4,36), (27,163), (28,19), (29,240), (43,123), (65,13), (69,16), (87,13), (93,239), (95,228), (97,176), (101,13), (104,208), (110,198), (111,244)],claim112),
  (.mul (variableData 5) 82,claim113),
  (.linear [(5,225), (30,26), (31,67), (104,176), (106,59), (112,256), (113,1)],claim114),
  (.mul (variableData 2) 87,claim115),
  (.linear [(2,15), (15,70), (16,207), (17,12), (18,68), (49,137), (50,185), (93,122), (95,198), (110,235), (112,140), (115,1)],claim116),
  (.mul (variableData 3) 87,claim117),
  (.linear [(3,158), (7,205), (21,67), (22,74), (24,215), (39,89), (43,134), (44,161), (57,31), (59,207), (89,208), (93,238), (95,145), (97,139), (101,255), (104,241), (106,13), (110,70), (114,235), (117,49)],claim118),
  (.mul (variableData 5) 87,claim119),
  (.linear [(5,77), (30,142), (31,39), (73,35), (104,99), (106,256), (112,120), (118,88), (119,222)],claim120),
  (.mul (variableData 3) 91,claim121),
  (.linear [(3,112), (22,233), (39,47), (57,85), (83,106), (91,170), (93,95), (95,32), (99,106), (101,32), (108,142), (110,78), (112,139), (121,151)],claim122),
  (.mul (variableData 4) 91,claim123),
  (.linear [(4,248), (28,35), (66,73), (69,82), (101,253), (108,73), (110,124), (112,198), (122,222), (123,184)],claim124),
  (.mul (variableData 5) 91,claim125),
  (.linear [(112,162), (116,256), (118,118), (122,180), (124,205), (125,1)],claim126),
  (.mul (variableData 4) 93,claim127),
  (.linear [(4,237), (27,241), (28,117), (29,80), (43,98), (69,253), (70,133), (71,79), (95,216), (97,108), (101,61), (104,27), (106,141), (110,79), (112,193), (127,178)],claim128),
  (.mul (variableData 5) 93,claim129),
  (.linear [(5,71), (30,119), (31,20), (104,238), (106,139), (112,22), (118,256), (120,187), (128,177), (129,1)],claim130),
  (.mul (variableData 5) 95,claim131),
  (.linear [(5,140), (31,8), (32,169), (104,170), (106,125), (120,55), (128,21), (130,135), (131,202)],claim132),
  (.mul (variableData 3) 101,claim133),
  (.linear [(3,195), (7,22), (21,166), (22,41), (24,4), (39,228), (43,141), (44,19), (61,147), (93,173), (95,137), (97,40), (99,236), (101,51), (102,187), (104,97), (106,94), (110,208), (112,182), (122,207), (124,155), (126,42), (128,128), (132,209), (133,21)],claim134),
  (.mul (variableData 5) 101,claim135),
  (.linear [(5,52), (30,111), (31,122), (32,117), (74,78), (75,201), (104,50), (106,163), (112,47), (122,37), (124,174), (126,229), (128,167), (130,196), (132,82), (134,21), (135,28)],claim136),
  (.mul (variableData 5) 104,claim137),
  (.linear [(5,137), (31,11), (32,162), (33,138), (34,105), (76,224), (77,6), (78,147), (79,191), (104,211), (106,67), (128,28), (130,26), (132,34), (137,33)],claim138),
  (.mul (variableData 2) 112,claim139),
  (.linear [(2,80), (15,164), (16,169), (17,127), (18,190), (19,31), (48,159), (50,17), (51,70), (52,210), (93,252), (95,93), (112,18), (122,42), (124,176), (128,176), (132,82), (136,112), (139,98)],claim140),
  (.mul (variableData 3) 112,claim141),
  (.linear [(3,212), (7,80), (21,152), (22,97), (24,78), (39,103), (43,209), (44,121), (61,41), (93,249), (95,212), (97,133), (101,41), (102,255), (104,91), (106,177), (110,100), (112,155), (118,32), (122,169), (124,17), (126,104), (128,57), (132,66), (134,157), (136,3), (138,129), (140,37), (141,2)],claim142),
  (.mul (variableData 5) 112,claim143),
  (.linear [(5,47), (30,8), (31,200), (32,128), (74,242), (75,231), (104,112), (106,214), (112,130), (122,75), (124,228), (128,79), (130,179), (132,115), (134,236), (136,162), (138,214), (140,104), (142,27), (143,15)],claim144),
  (.mul (variableData 0) 122,claim145),
  (.linear [(0,124), (7,241), (8,249), (11,242), (36,94), (38,53), (39,149), (41,110), (43,212), (82,110), (87,115), (91,138), (93,138), (95,80), (99,22), (101,99), (108,217), (110,141), (112,92), (116,38), (122,51), (124,235), (128,127), (132,249), (136,221), (140,243), (142,246), (144,241), (145,40)],claim146),
  (.mul (variableData 4) 122,claim147),
  (.linear [(4,131), (27,154), (28,48), (29,192), (43,137), (69,47), (71,244), (95,108), (97,114), (101,215), (104,88), (106,137), (112,140), (122,112), (124,29), (128,187), (132,194), (134,71), (136,16), (138,45), (140,171), (142,105), (144,163), (146,191), (147,210)],claim148),
  (.mul (variableData 5) 122,claim149),
  (.linear [(5,232), (30,174), (31,208), (32,82), (75,176), (104,246), (106,206), (112,247), (122,14), (124,186), (128,71), (130,156), (132,250), (134,207), (136,140), (138,251), (140,231), (142,132), (144,66), (146,255), (148,168), (149,50)],claim150),
  (.mul (variableData 5) 124,claim151),
  (.linear [(5,217), (30,43), (31,46), (32,144), (75,34), (104,145), (106,82), (112,43), (122,94), (124,246), (128,210), (130,237), (132,90), (136,17), (138,56), (140,239), (142,97), (144,76), (146,232), (148,102), (150,134), (151,223)],claim152),
  (.mul (variableData 5) 128,claim153),
  (.linear [(5,178), (30,124), (31,256), (32,244), (33,176), (34,141), (78,248), (79,232), (80,253), (104,188), (106,80), (128,201), (130,82), (132,36), (136,207), (138,216), (140,77), (142,193), (144,91), (146,191), (148,222), (150,124), (152,111), (153,175)],claim154),
  (.mul (variableData 2) 132,claim155),
  (.linear [(2,36), (17,30), (18,8), (19,9), (53,108), (54,71), (93,168), (95,62), (112,124), (122,92), (124,100), (128,204), (130,145), (132,197), (136,177), (138,189), (140,200), (142,245), (144,203), (148,95), (150,42), (152,242), (154,225), (155,149)],claim156),
  (.mul (variableData 3) 132,claim157),
  (.linear [(3,42), (7,156), (21,10), (22,167), (24,214), (25,148), (43,112), (44,142), (62,190), (64,79), (95,83), (97,69), (104,24), (106,66), (128,117), (130,44), (132,63), (136,178), (138,154), (140,24), (142,182), (144,65), (146,122), (148,187), (150,143), (152,253), (154,173), (156,11), (157,67)],claim158),
  (.mul (variableData 2) 136,claim159),
  (.linear [(2,12), (15,24), (16,188), (17,137), (18,256), (19,212), (54,79), (93,9), (95,141), (112,2), (122,233), (124,112), (128,32), (130,136), (132,31), (136,136), (138,134), (140,52), (142,38), (144,143), (146,40), (148,216), (150,43), (152,212), (154,205), (156,48), (158,237), (159,178)],claim160),
  (.mul (variableData 3) 136,claim161),
  (.linear [(3,234), (7,81), (21,81), (22,205), (24,229), (25,247), (39,216), (43,157), (44,18), (64,196), (93,226), (95,126), (97,92), (101,76), (104,62), (112,175), (122,76), (124,148), (128,121), (130,197), (132,178), (136,238), (138,180), (140,45), (142,91), (144,199), (146,136), (148,170), (150,14), (152,207), (154,123), (156,11), (158,216), (160,5), (161,60)],claim162),
  (.mul (variableData 3) 140,claim163),
  (.linear [(3,48), (7,243), (21,252), (22,87), (24,210), (25,128), (39,241), (43,214), (44,63), (93,36), (95,120), (97,65), (101,61), (104,14), (106,150), (112,225), (122,155), (124,85), (128,138), (132,179), (136,89), (138,146), (140,223), (142,202), (144,75), (146,16), (148,28), (150,155), (152,244), (154,96), (156,246), (158,167), (160,116), (162,164), (163,16)],claim164),
  (.mul (variableData 5) 140,claim165),
  (.linear [(5,243), (30,164), (31,4), (32,68), (33,117), (34,17), (104,142), (106,234), (112,146), (122,104), (124,76), (128,16), (132,83), (136,205), (138,86), (140,34), (142,218), (144,6), (146,67), (148,193), (150,104), (152,252), (154,103), (156,181), (158,114), (160,255), (162,46), (164,252), (165,111)],claim166),
  (.mul (variableData 5) 142,claim167),
  (.linear [(5,175), (30,149), (31,125), (32,119), (33,241), (34,63), (104,172), (106,1), (124,82), (128,179), (132,8), (136,107), (138,219), (140,218), (142,130), (144,61), (146,74), (148,105), (150,170), (152,128), (154,16), (156,159), (158,167), (160,210), (162,20), (164,112), (166,217), (167,108)],claim168),
  (.mul (variableData 5) 146,claim169),
  (.linear [(5,158), (31,27), (32,71), (33,191), (34,208), (104,151), (106,190), (128,8), (132,79), (136,102), (138,144), (140,213), (142,22), (144,43), (146,54), (148,97), (150,176), (152,244), (154,102), (156,225), (158,182), (160,188), (162,95), (164,212), (166,46), (168,95), (169,106)],claim170),
  (.mul (variableData 5) 148,claim171),
  (.linear [(5,223), (31,155), (32,233), (33,224), (34,88), (106,117), (128,108), (132,218), (136,181), (138,174), (140,189), (142,239), (144,181), (146,87), (148,199), (150,24), (152,9), (154,4), (156,174), (158,192), (160,70), (162,184), (164,69), (166,173), (168,168), (170,6), (171,140)],claim172),
  (.mul (variableData 3) 156,claim173),
  (.linear [(3,93), (7,87), (24,70), (25,6), (43,102), (132,185), (136,133), (138,16), (140,119), (142,92), (144,148), (146,145), (148,21), (150,34), (152,254), (154,94), (156,199), (158,86), (160,120), (162,126), (164,245), (166,39), (168,5), (170,21), (172,231), (173,72)],claim174),
  (.mul (variableData 5) 164,claim175),
  (.linear [(5,191), (142,18), (144,185), (146,42), (148,40), (150,162), (152,43), (154,211), (156,149), (158,31), (160,206), (162,65), (164,89), (166,27), (168,225), (170,27), (172,102), (174,3), (175,239)],claim176)
]
def prefix0 : List (Sparse 6 K) := []
def prefix1 : List (Sparse 6 K) := prefix0++[claim0]
def prefix2 : List (Sparse 6 K) := prefix1++[claim1]
def prefix3 : List (Sparse 6 K) := prefix2++[claim2]
def prefix4 : List (Sparse 6 K) := prefix3++[claim3]
def prefix5 : List (Sparse 6 K) := prefix4++[claim4]
def prefix6 : List (Sparse 6 K) := prefix5++[claim5]
def prefix7 : List (Sparse 6 K) := prefix6++[claim6]
def prefix8 : List (Sparse 6 K) := prefix7++[claim7]
def prefix9 : List (Sparse 6 K) := prefix8++[claim8]
def prefix10 : List (Sparse 6 K) := prefix9++[claim9]
def prefix11 : List (Sparse 6 K) := prefix10++[claim10]
def prefix12 : List (Sparse 6 K) := prefix11++[claim11]
def prefix13 : List (Sparse 6 K) := prefix12++[claim12]
def prefix14 : List (Sparse 6 K) := prefix13++[claim13]
def prefix15 : List (Sparse 6 K) := prefix14++[claim14]
def prefix16 : List (Sparse 6 K) := prefix15++[claim15]
def prefix17 : List (Sparse 6 K) := prefix16++[claim16]
def prefix18 : List (Sparse 6 K) := prefix17++[claim17]
def prefix19 : List (Sparse 6 K) := prefix18++[claim18]
def prefix20 : List (Sparse 6 K) := prefix19++[claim19]
def prefix21 : List (Sparse 6 K) := prefix20++[claim20]
def prefix22 : List (Sparse 6 K) := prefix21++[claim21]
def prefix23 : List (Sparse 6 K) := prefix22++[claim22]
def prefix24 : List (Sparse 6 K) := prefix23++[claim23]
def prefix25 : List (Sparse 6 K) := prefix24++[claim24]
def prefix26 : List (Sparse 6 K) := prefix25++[claim25]
def prefix27 : List (Sparse 6 K) := prefix26++[claim26]
def prefix28 : List (Sparse 6 K) := prefix27++[claim27]
def prefix29 : List (Sparse 6 K) := prefix28++[claim28]
def prefix30 : List (Sparse 6 K) := prefix29++[claim29]
def prefix31 : List (Sparse 6 K) := prefix30++[claim30]
def prefix32 : List (Sparse 6 K) := prefix31++[claim31]
def prefix33 : List (Sparse 6 K) := prefix32++[claim32]
def prefix34 : List (Sparse 6 K) := prefix33++[claim33]
def prefix35 : List (Sparse 6 K) := prefix34++[claim34]
def prefix36 : List (Sparse 6 K) := prefix35++[claim35]
def prefix37 : List (Sparse 6 K) := prefix36++[claim36]
def prefix38 : List (Sparse 6 K) := prefix37++[claim37]
def prefix39 : List (Sparse 6 K) := prefix38++[claim38]
def prefix40 : List (Sparse 6 K) := prefix39++[claim39]
def prefix41 : List (Sparse 6 K) := prefix40++[claim40]
def prefix42 : List (Sparse 6 K) := prefix41++[claim41]
def prefix43 : List (Sparse 6 K) := prefix42++[claim42]
def prefix44 : List (Sparse 6 K) := prefix43++[claim43]
def prefix45 : List (Sparse 6 K) := prefix44++[claim44]
def prefix46 : List (Sparse 6 K) := prefix45++[claim45]
def prefix47 : List (Sparse 6 K) := prefix46++[claim46]
def prefix48 : List (Sparse 6 K) := prefix47++[claim47]
def prefix49 : List (Sparse 6 K) := prefix48++[claim48]
def prefix50 : List (Sparse 6 K) := prefix49++[claim49]
def prefix51 : List (Sparse 6 K) := prefix50++[claim50]
def prefix52 : List (Sparse 6 K) := prefix51++[claim51]
def prefix53 : List (Sparse 6 K) := prefix52++[claim52]
def prefix54 : List (Sparse 6 K) := prefix53++[claim53]
def prefix55 : List (Sparse 6 K) := prefix54++[claim54]
def prefix56 : List (Sparse 6 K) := prefix55++[claim55]
def prefix57 : List (Sparse 6 K) := prefix56++[claim56]
def prefix58 : List (Sparse 6 K) := prefix57++[claim57]
def prefix59 : List (Sparse 6 K) := prefix58++[claim58]
def prefix60 : List (Sparse 6 K) := prefix59++[claim59]
def prefix61 : List (Sparse 6 K) := prefix60++[claim60]
def prefix62 : List (Sparse 6 K) := prefix61++[claim61]
def prefix63 : List (Sparse 6 K) := prefix62++[claim62]
def prefix64 : List (Sparse 6 K) := prefix63++[claim63]
def prefix65 : List (Sparse 6 K) := prefix64++[claim64]
def prefix66 : List (Sparse 6 K) := prefix65++[claim65]
def prefix67 : List (Sparse 6 K) := prefix66++[claim66]
def prefix68 : List (Sparse 6 K) := prefix67++[claim67]
def prefix69 : List (Sparse 6 K) := prefix68++[claim68]
def prefix70 : List (Sparse 6 K) := prefix69++[claim69]
def prefix71 : List (Sparse 6 K) := prefix70++[claim70]
def prefix72 : List (Sparse 6 K) := prefix71++[claim71]
def prefix73 : List (Sparse 6 K) := prefix72++[claim72]
def prefix74 : List (Sparse 6 K) := prefix73++[claim73]
def prefix75 : List (Sparse 6 K) := prefix74++[claim74]
def prefix76 : List (Sparse 6 K) := prefix75++[claim75]
def prefix77 : List (Sparse 6 K) := prefix76++[claim76]
def prefix78 : List (Sparse 6 K) := prefix77++[claim77]
def prefix79 : List (Sparse 6 K) := prefix78++[claim78]
def prefix80 : List (Sparse 6 K) := prefix79++[claim79]
def prefix81 : List (Sparse 6 K) := prefix80++[claim80]
def prefix82 : List (Sparse 6 K) := prefix81++[claim81]
def prefix83 : List (Sparse 6 K) := prefix82++[claim82]
def prefix84 : List (Sparse 6 K) := prefix83++[claim83]
def prefix85 : List (Sparse 6 K) := prefix84++[claim84]
def prefix86 : List (Sparse 6 K) := prefix85++[claim85]
def prefix87 : List (Sparse 6 K) := prefix86++[claim86]
def prefix88 : List (Sparse 6 K) := prefix87++[claim87]
def prefix89 : List (Sparse 6 K) := prefix88++[claim88]
def prefix90 : List (Sparse 6 K) := prefix89++[claim89]
def prefix91 : List (Sparse 6 K) := prefix90++[claim90]
def prefix92 : List (Sparse 6 K) := prefix91++[claim91]
def prefix93 : List (Sparse 6 K) := prefix92++[claim92]
def prefix94 : List (Sparse 6 K) := prefix93++[claim93]
def prefix95 : List (Sparse 6 K) := prefix94++[claim94]
def prefix96 : List (Sparse 6 K) := prefix95++[claim95]
def prefix97 : List (Sparse 6 K) := prefix96++[claim96]
def prefix98 : List (Sparse 6 K) := prefix97++[claim97]
def prefix99 : List (Sparse 6 K) := prefix98++[claim98]
def prefix100 : List (Sparse 6 K) := prefix99++[claim99]
def prefix101 : List (Sparse 6 K) := prefix100++[claim100]
def prefix102 : List (Sparse 6 K) := prefix101++[claim101]
def prefix103 : List (Sparse 6 K) := prefix102++[claim102]
def prefix104 : List (Sparse 6 K) := prefix103++[claim103]
def prefix105 : List (Sparse 6 K) := prefix104++[claim104]
def prefix106 : List (Sparse 6 K) := prefix105++[claim105]
def prefix107 : List (Sparse 6 K) := prefix106++[claim106]
def prefix108 : List (Sparse 6 K) := prefix107++[claim107]
def prefix109 : List (Sparse 6 K) := prefix108++[claim108]
def prefix110 : List (Sparse 6 K) := prefix109++[claim109]
def prefix111 : List (Sparse 6 K) := prefix110++[claim110]
def prefix112 : List (Sparse 6 K) := prefix111++[claim111]
def prefix113 : List (Sparse 6 K) := prefix112++[claim112]
def prefix114 : List (Sparse 6 K) := prefix113++[claim113]
def prefix115 : List (Sparse 6 K) := prefix114++[claim114]
def prefix116 : List (Sparse 6 K) := prefix115++[claim115]
def prefix117 : List (Sparse 6 K) := prefix116++[claim116]
def prefix118 : List (Sparse 6 K) := prefix117++[claim117]
def prefix119 : List (Sparse 6 K) := prefix118++[claim118]
def prefix120 : List (Sparse 6 K) := prefix119++[claim119]
def prefix121 : List (Sparse 6 K) := prefix120++[claim120]
def prefix122 : List (Sparse 6 K) := prefix121++[claim121]
def prefix123 : List (Sparse 6 K) := prefix122++[claim122]
def prefix124 : List (Sparse 6 K) := prefix123++[claim123]
def prefix125 : List (Sparse 6 K) := prefix124++[claim124]
def prefix126 : List (Sparse 6 K) := prefix125++[claim125]
def prefix127 : List (Sparse 6 K) := prefix126++[claim126]
def prefix128 : List (Sparse 6 K) := prefix127++[claim127]
def prefix129 : List (Sparse 6 K) := prefix128++[claim128]
def prefix130 : List (Sparse 6 K) := prefix129++[claim129]
def prefix131 : List (Sparse 6 K) := prefix130++[claim130]
def prefix132 : List (Sparse 6 K) := prefix131++[claim131]
def prefix133 : List (Sparse 6 K) := prefix132++[claim132]
def prefix134 : List (Sparse 6 K) := prefix133++[claim133]
def prefix135 : List (Sparse 6 K) := prefix134++[claim134]
def prefix136 : List (Sparse 6 K) := prefix135++[claim135]
def prefix137 : List (Sparse 6 K) := prefix136++[claim136]
def prefix138 : List (Sparse 6 K) := prefix137++[claim137]
def prefix139 : List (Sparse 6 K) := prefix138++[claim138]
def prefix140 : List (Sparse 6 K) := prefix139++[claim139]
def prefix141 : List (Sparse 6 K) := prefix140++[claim140]
def prefix142 : List (Sparse 6 K) := prefix141++[claim141]
def prefix143 : List (Sparse 6 K) := prefix142++[claim142]
def prefix144 : List (Sparse 6 K) := prefix143++[claim143]
def prefix145 : List (Sparse 6 K) := prefix144++[claim144]
def prefix146 : List (Sparse 6 K) := prefix145++[claim145]
def prefix147 : List (Sparse 6 K) := prefix146++[claim146]
def prefix148 : List (Sparse 6 K) := prefix147++[claim147]
def prefix149 : List (Sparse 6 K) := prefix148++[claim148]
def prefix150 : List (Sparse 6 K) := prefix149++[claim149]
def prefix151 : List (Sparse 6 K) := prefix150++[claim150]
def prefix152 : List (Sparse 6 K) := prefix151++[claim151]
def prefix153 : List (Sparse 6 K) := prefix152++[claim152]
def prefix154 : List (Sparse 6 K) := prefix153++[claim153]
def prefix155 : List (Sparse 6 K) := prefix154++[claim154]
def prefix156 : List (Sparse 6 K) := prefix155++[claim155]
def prefix157 : List (Sparse 6 K) := prefix156++[claim156]
def prefix158 : List (Sparse 6 K) := prefix157++[claim157]
def prefix159 : List (Sparse 6 K) := prefix158++[claim158]
def prefix160 : List (Sparse 6 K) := prefix159++[claim159]
def prefix161 : List (Sparse 6 K) := prefix160++[claim160]
def prefix162 : List (Sparse 6 K) := prefix161++[claim161]
def prefix163 : List (Sparse 6 K) := prefix162++[claim162]
def prefix164 : List (Sparse 6 K) := prefix163++[claim163]
def prefix165 : List (Sparse 6 K) := prefix164++[claim164]
def prefix166 : List (Sparse 6 K) := prefix165++[claim165]
def prefix167 : List (Sparse 6 K) := prefix166++[claim166]
def prefix168 : List (Sparse 6 K) := prefix167++[claim167]
def prefix169 : List (Sparse 6 K) := prefix168++[claim168]
def prefix170 : List (Sparse 6 K) := prefix169++[claim169]
def prefix171 : List (Sparse 6 K) := prefix170++[claim170]
def prefix172 : List (Sparse 6 K) := prefix171++[claim171]
def prefix173 : List (Sparse 6 K) := prefix172++[claim172]
def prefix174 : List (Sparse 6 K) := prefix173++[claim173]
def prefix175 : List (Sparse 6 K) := prefix174++[claim174]
def prefix176 : List (Sparse 6 K) := prefix175++[claim175]
def prefix177 : List (Sparse 6 K) := prefix176++[claim176]
theorem step_checked0 : checkEntry inputs 4 prefix0 (.input 0,claim0)=true := by rfl
theorem step_checked1 : checkEntry inputs 4 prefix1 (.input 1,claim1)=true := by rfl
theorem step_checked2 : checkEntry inputs 4 prefix2 (.input 2,claim2)=true := by rfl
theorem step_checked3 : checkEntry inputs 4 prefix3 (.input 3,claim3)=true := by rfl
theorem step_checked4 : checkEntry inputs 4 prefix4 (.input 4,claim4)=true := by rfl
theorem step_checked5 : checkEntry inputs 4 prefix5 (.input 5,claim5)=true := by rfl
theorem step_checked6 : checkEntry inputs 4 prefix6 (.input 6,claim6)=true := by rfl
theorem step_checked7 : checkEntry inputs 4 prefix7 (.linear [(6,134)],claim7)=true := by rfl
theorem step_checked8 : checkEntry inputs 4 prefix8 (.mul (variableData 2) 0,claim8)=true := by rfl
theorem step_checked9 : checkEntry inputs 4 prefix9 (.mul (variableData 3) 0,claim9)=true := by rfl
theorem step_checked10 : checkEntry inputs 4 prefix10 (.mul (variableData 4) 0,claim10)=true := by rfl
theorem step_checked11 : checkEntry inputs 4 prefix11 (.mul (variableData 5) 0,claim11)=true := by rfl
theorem step_checked12 : checkEntry inputs 4 prefix12 (.mul (variableData 0) 1,claim12)=true := by rfl
theorem step_checked13 : checkEntry inputs 4 prefix13 (.mul (variableData 3) 1,claim13)=true := by rfl
theorem step_checked14 : checkEntry inputs 4 prefix14 (.mul (variableData 4) 1,claim14)=true := by rfl
theorem step_checked15 : checkEntry inputs 4 prefix15 (.mul (variableData 0) 2,claim15)=true := by rfl
theorem step_checked16 : checkEntry inputs 4 prefix16 (.mul (variableData 1) 2,claim16)=true := by rfl
theorem step_checked17 : checkEntry inputs 4 prefix17 (.mul (variableData 3) 2,claim17)=true := by rfl
theorem step_checked18 : checkEntry inputs 4 prefix18 (.mul (variableData 4) 2,claim18)=true := by rfl
theorem step_checked19 : checkEntry inputs 4 prefix19 (.mul (variableData 5) 2,claim19)=true := by rfl
theorem step_checked20 : checkEntry inputs 4 prefix20 (.mul (variableData 0) 3,claim20)=true := by rfl
theorem step_checked21 : checkEntry inputs 4 prefix21 (.mul (variableData 1) 3,claim21)=true := by rfl
theorem step_checked22 : checkEntry inputs 4 prefix22 (.mul (variableData 2) 3,claim22)=true := by rfl
theorem step_checked23 : checkEntry inputs 4 prefix23 (.mul (variableData 4) 3,claim23)=true := by rfl
theorem step_checked24 : checkEntry inputs 4 prefix24 (.linear [(4,1), (23,1)],claim24)=true := by rfl
theorem step_checked25 : checkEntry inputs 4 prefix25 (.mul (variableData 5) 3,claim25)=true := by rfl
theorem step_checked26 : checkEntry inputs 4 prefix26 (.mul (variableData 0) 4,claim26)=true := by rfl
theorem step_checked27 : checkEntry inputs 4 prefix27 (.mul (variableData 1) 4,claim27)=true := by rfl
theorem step_checked28 : checkEntry inputs 4 prefix28 (.mul (variableData 2) 4,claim28)=true := by rfl
theorem step_checked29 : checkEntry inputs 4 prefix29 (.mul (variableData 3) 4,claim29)=true := by rfl
theorem step_checked30 : checkEntry inputs 4 prefix30 (.mul (variableData 0) 5,claim30)=true := by rfl
theorem step_checked31 : checkEntry inputs 4 prefix31 (.mul (variableData 1) 5,claim31)=true := by rfl
theorem step_checked32 : checkEntry inputs 4 prefix32 (.mul (variableData 2) 5,claim32)=true := by rfl
theorem step_checked33 : checkEntry inputs 4 prefix33 (.mul (variableData 3) 5,claim33)=true := by rfl
theorem step_checked34 : checkEntry inputs 4 prefix34 (.mul (variableData 4) 5,claim34)=true := by rfl
theorem step_checked35 : checkEntry inputs 4 prefix35 (.mul (variableData 0) 7,claim35)=true := by rfl
theorem step_checked36 : checkEntry inputs 4 prefix36 (.linear [(0,57), (7,57), (9,84), (10,256), (35,173)],claim36)=true := by rfl
theorem step_checked37 : checkEntry inputs 4 prefix37 (.mul (variableData 1) 7,claim37)=true := by rfl
theorem step_checked38 : checkEntry inputs 4 prefix38 (.linear [(1,221), (13,222), (14,236), (36,113), (37,144)],claim38)=true := by rfl
theorem step_checked39 : checkEntry inputs 4 prefix39 (.mul (variableData 2) 7,claim39)=true := by rfl
theorem step_checked40 : checkEntry inputs 4 prefix40 (.mul (variableData 3) 7,claim40)=true := by rfl
theorem step_checked41 : checkEntry inputs 4 prefix41 (.linear [(3,57), (7,57), (20,84), (21,256), (40,173)],claim41)=true := by rfl
theorem step_checked42 : checkEntry inputs 4 prefix42 (.mul (variableData 4) 7,claim42)=true := by rfl
theorem step_checked43 : checkEntry inputs 4 prefix43 (.linear [(4,70), (26,168), (27,255), (41,23), (42,234)],claim43)=true := by rfl
theorem step_checked44 : checkEntry inputs 4 prefix44 (.mul (variableData 5) 7,claim44)=true := by rfl
theorem step_checked45 : checkEntry inputs 4 prefix45 (.mul (variableData 3) 12,claim45)=true := by rfl
theorem step_checked46 : checkEntry inputs 4 prefix46 (.linear [(39,1), (45,1)],claim46)=true := by rfl
theorem step_checked47 : checkEntry inputs 4 prefix47 (.mul (variableData 4) 12,claim47)=true := by rfl
theorem step_checked48 : checkEntry inputs 4 prefix48 (.mul (variableData 5) 15,claim48)=true := by rfl
theorem step_checked49 : checkEntry inputs 4 prefix49 (.mul (variableData 3) 16,claim49)=true := by rfl
theorem step_checked50 : checkEntry inputs 4 prefix50 (.mul (variableData 4) 16,claim50)=true := by rfl
theorem step_checked51 : checkEntry inputs 4 prefix51 (.mul (variableData 5) 16,claim51)=true := by rfl
theorem step_checked52 : checkEntry inputs 4 prefix52 (.mul (variableData 4) 17,claim52)=true := by rfl
theorem step_checked53 : checkEntry inputs 4 prefix53 (.mul (variableData 5) 17,claim53)=true := by rfl
theorem step_checked54 : checkEntry inputs 4 prefix54 (.mul (variableData 5) 18,claim54)=true := by rfl
theorem step_checked55 : checkEntry inputs 4 prefix55 (.mul (variableData 1) 20,claim55)=true := by rfl
theorem step_checked56 : checkEntry inputs 4 prefix56 (.mul (variableData 2) 20,claim56)=true := by rfl
theorem step_checked57 : checkEntry inputs 4 prefix57 (.mul (variableData 2) 21,claim57)=true := by rfl
theorem step_checked58 : checkEntry inputs 4 prefix58 (.mul (variableData 4) 21,claim58)=true := by rfl
theorem step_checked59 : checkEntry inputs 4 prefix59 (.linear [(27,1), (58,1)],claim59)=true := by rfl
theorem step_checked60 : checkEntry inputs 4 prefix60 (.mul (variableData 4) 22,claim60)=true := by rfl
theorem step_checked61 : checkEntry inputs 4 prefix61 (.linear [(28,1), (60,1)],claim61)=true := by rfl
theorem step_checked62 : checkEntry inputs 4 prefix62 (.mul (variableData 5) 22,claim62)=true := by rfl
theorem step_checked63 : checkEntry inputs 4 prefix63 (.mul (variableData 5) 24,claim63)=true := by rfl
theorem step_checked64 : checkEntry inputs 4 prefix64 (.linear [(5,1), (63,1)],claim64)=true := by rfl
theorem step_checked65 : checkEntry inputs 4 prefix65 (.mul (variableData 1) 26,claim65)=true := by rfl
theorem step_checked66 : checkEntry inputs 4 prefix66 (.mul (variableData 2) 26,claim66)=true := by rfl
theorem step_checked67 : checkEntry inputs 4 prefix67 (.mul (variableData 3) 26,claim67)=true := by rfl
theorem step_checked68 : checkEntry inputs 4 prefix68 (.linear [(44,1), (67,1)],claim68)=true := by rfl
theorem step_checked69 : checkEntry inputs 4 prefix69 (.mul (variableData 2) 27,claim69)=true := by rfl
theorem step_checked70 : checkEntry inputs 4 prefix70 (.mul (variableData 3) 27,claim70)=true := by rfl
theorem step_checked71 : checkEntry inputs 4 prefix71 (.mul (variableData 3) 28,claim71)=true := by rfl
theorem step_checked72 : checkEntry inputs 4 prefix72 (.mul (variableData 1) 30,claim72)=true := by rfl
theorem step_checked73 : checkEntry inputs 4 prefix73 (.linear [(38,1), (43,256), (72,1)],claim73)=true := by rfl
theorem step_checked74 : checkEntry inputs 4 prefix74 (.mul (variableData 2) 30,claim74)=true := by rfl
theorem step_checked75 : checkEntry inputs 4 prefix75 (.mul (variableData 2) 31,claim75)=true := by rfl
theorem step_checked76 : checkEntry inputs 4 prefix76 (.mul (variableData 3) 31,claim76)=true := by rfl
theorem step_checked77 : checkEntry inputs 4 prefix77 (.mul (variableData 4) 31,claim77)=true := by rfl
theorem step_checked78 : checkEntry inputs 4 prefix78 (.mul (variableData 3) 32,claim78)=true := by rfl
theorem step_checked79 : checkEntry inputs 4 prefix79 (.mul (variableData 4) 32,claim79)=true := by rfl
theorem step_checked80 : checkEntry inputs 4 prefix80 (.mul (variableData 4) 33,claim80)=true := by rfl
theorem step_checked81 : checkEntry inputs 4 prefix81 (.mul (variableData 1) 36,claim81)=true := by rfl
theorem step_checked82 : checkEntry inputs 4 prefix82 (.linear [(1,12), (12,60), (13,179), (14,56), (38,162), (43,95), (46,144), (47,35), (81,113)],claim82)=true := by rfl
theorem step_checked83 : checkEntry inputs 4 prefix83 (.mul (variableData 2) 36,claim83)=true := by rfl
theorem step_checked84 : checkEntry inputs 4 prefix84 (.mul (variableData 3) 36,claim84)=true := by rfl
theorem step_checked85 : checkEntry inputs 4 prefix85 (.linear [(3,165), (7,223), (21,83), (36,35), (38,221), (41,235), (43,195), (55,84), (82,63), (84,173)],claim85)=true := by rfl
theorem step_checked86 : checkEntry inputs 4 prefix86 (.mul (variableData 4) 36,claim86)=true := by rfl
theorem step_checked87 : checkEntry inputs 4 prefix87 (.linear [(4,169), (26,235), (27,204), (43,223), (65,256), (82,35), (85,84), (86,173)],claim87)=true := by rfl
theorem step_checked88 : checkEntry inputs 4 prefix88 (.mul (variableData 5) 36,claim88)=true := by rfl
theorem step_checked89 : checkEntry inputs 4 prefix89 (.linear [(87,21), (88,1)],claim89)=true := by rfl
theorem step_checked90 : checkEntry inputs 4 prefix90 (.mul (variableData 2) 38,claim90)=true := by rfl
theorem step_checked91 : checkEntry inputs 4 prefix91 (.linear [(2,221), (17,222), (18,236), (90,1)],claim91)=true := by rfl
theorem step_checked92 : checkEntry inputs 4 prefix92 (.mul (variableData 3) 38,claim92)=true := by rfl
theorem step_checked93 : checkEntry inputs 4 prefix93 (.linear [(3,239), (7,228), (21,128), (22,84), (36,208), (38,256), (41,208), (43,24), (82,235), (92,49)],claim93)=true := by rfl
theorem step_checked94 : checkEntry inputs 4 prefix94 (.mul (variableData 4) 38,claim94)=true := by rfl
theorem step_checked95 : checkEntry inputs 4 prefix95 (.linear [(4,108), (26,27), (27,231), (28,53), (43,252), (82,27), (93,83), (94,230)],claim95)=true := by rfl
theorem step_checked96 : checkEntry inputs 4 prefix96 (.mul (variableData 5) 38,claim96)=true := by rfl
theorem step_checked97 : checkEntry inputs 4 prefix97 (.linear [(87,256), (93,120), (95,103), (96,1)],claim97)=true := by rfl
theorem step_checked98 : checkEntry inputs 4 prefix98 (.mul (variableData 3) 39,claim98)=true := by rfl
theorem step_checked99 : checkEntry inputs 4 prefix99 (.linear [(22,57), (39,57), (56,84), (57,256), (93,141), (95,153), (98,173)],claim99)=true := by rfl
theorem step_checked100 : checkEntry inputs 4 prefix100 (.mul (variableData 4) 39,claim100)=true := by rfl
theorem step_checked101 : checkEntry inputs 4 prefix101 (.linear [(28,70), (66,168), (69,255), (99,23), (100,234)],claim101)=true := by rfl
theorem step_checked102 : checkEntry inputs 4 prefix102 (.mul (variableData 5) 39,claim102)=true := by rfl
theorem step_checked103 : checkEntry inputs 4 prefix103 (.mul (variableData 4) 41,claim103)=true := by rfl
theorem step_checked104 : checkEntry inputs 4 prefix104 (.linear [(4,172), (26,167), (27,203), (29,89), (43,236), (68,8), (70,159), (95,157), (97,98), (103,249)],claim104)=true := by rfl
theorem step_checked105 : checkEntry inputs 4 prefix105 (.mul (variableData 5) 43,claim105)=true := by rfl
theorem step_checked106 : checkEntry inputs 4 prefix106 (.linear [(5,209), (30,39), (31,229), (97,243), (104,60), (105,14)],claim106)=true := by rfl
theorem step_checked107 : checkEntry inputs 4 prefix107 (.mul (variableData 2) 82,claim107)=true := by rfl
theorem step_checked108 : checkEntry inputs 4 prefix108 (.linear [(2,24), (15,256), (16,221), (17,102), (18,216), (49,222), (50,236), (93,176), (95,16), (107,1)],claim108)=true := by rfl
theorem step_checked109 : checkEntry inputs 4 prefix109 (.mul (variableData 3) 82,claim109)=true := by rfl
theorem step_checked110 : checkEntry inputs 4 prefix110 (.linear [(3,28), (7,169), (21,143), (22,115), (24,242), (39,208), (43,207), (44,246), (57,84), (85,208), (93,118), (95,256), (97,58), (101,235), (104,162), (109,49)],claim110)=true := by rfl
theorem step_checked111 : checkEntry inputs 4 prefix111 (.mul (variableData 4) 82,claim111)=true := by rfl
theorem step_checked112 : checkEntry inputs 4 prefix112 (.linear [(4,36), (27,163), (28,19), (29,240), (43,123), (65,13), (69,16), (87,13), (93,239), (95,228), (97,176), (101,13), (104,208), (110,198), (111,244)],claim112)=true := by rfl
theorem step_checked113 : checkEntry inputs 4 prefix113 (.mul (variableData 5) 82,claim113)=true := by rfl
theorem step_checked114 : checkEntry inputs 4 prefix114 (.linear [(5,225), (30,26), (31,67), (104,176), (106,59), (112,256), (113,1)],claim114)=true := by rfl
theorem step_checked115 : checkEntry inputs 4 prefix115 (.mul (variableData 2) 87,claim115)=true := by rfl
theorem step_checked116 : checkEntry inputs 4 prefix116 (.linear [(2,15), (15,70), (16,207), (17,12), (18,68), (49,137), (50,185), (93,122), (95,198), (110,235), (112,140), (115,1)],claim116)=true := by rfl
theorem step_checked117 : checkEntry inputs 4 prefix117 (.mul (variableData 3) 87,claim117)=true := by rfl
theorem step_checked118 : checkEntry inputs 4 prefix118 (.linear [(3,158), (7,205), (21,67), (22,74), (24,215), (39,89), (43,134), (44,161), (57,31), (59,207), (89,208), (93,238), (95,145), (97,139), (101,255), (104,241), (106,13), (110,70), (114,235), (117,49)],claim118)=true := by rfl
theorem step_checked119 : checkEntry inputs 4 prefix119 (.mul (variableData 5) 87,claim119)=true := by rfl
theorem step_checked120 : checkEntry inputs 4 prefix120 (.linear [(5,77), (30,142), (31,39), (73,35), (104,99), (106,256), (112,120), (118,88), (119,222)],claim120)=true := by rfl
theorem step_checked121 : checkEntry inputs 4 prefix121 (.mul (variableData 3) 91,claim121)=true := by rfl
theorem step_checked122 : checkEntry inputs 4 prefix122 (.linear [(3,112), (22,233), (39,47), (57,85), (83,106), (91,170), (93,95), (95,32), (99,106), (101,32), (108,142), (110,78), (112,139), (121,151)],claim122)=true := by rfl
theorem step_checked123 : checkEntry inputs 4 prefix123 (.mul (variableData 4) 91,claim123)=true := by rfl
theorem step_checked124 : checkEntry inputs 4 prefix124 (.linear [(4,248), (28,35), (66,73), (69,82), (101,253), (108,73), (110,124), (112,198), (122,222), (123,184)],claim124)=true := by rfl
theorem step_checked125 : checkEntry inputs 4 prefix125 (.mul (variableData 5) 91,claim125)=true := by rfl
theorem step_checked126 : checkEntry inputs 4 prefix126 (.linear [(112,162), (116,256), (118,118), (122,180), (124,205), (125,1)],claim126)=true := by rfl
theorem step_checked127 : checkEntry inputs 4 prefix127 (.mul (variableData 4) 93,claim127)=true := by rfl
theorem step_checked128 : checkEntry inputs 4 prefix128 (.linear [(4,237), (27,241), (28,117), (29,80), (43,98), (69,253), (70,133), (71,79), (95,216), (97,108), (101,61), (104,27), (106,141), (110,79), (112,193), (127,178)],claim128)=true := by rfl
theorem step_checked129 : checkEntry inputs 4 prefix129 (.mul (variableData 5) 93,claim129)=true := by rfl
theorem step_checked130 : checkEntry inputs 4 prefix130 (.linear [(5,71), (30,119), (31,20), (104,238), (106,139), (112,22), (118,256), (120,187), (128,177), (129,1)],claim130)=true := by rfl
theorem step_checked131 : checkEntry inputs 4 prefix131 (.mul (variableData 5) 95,claim131)=true := by rfl
theorem step_checked132 : checkEntry inputs 4 prefix132 (.linear [(5,140), (31,8), (32,169), (104,170), (106,125), (120,55), (128,21), (130,135), (131,202)],claim132)=true := by rfl
theorem step_checked133 : checkEntry inputs 4 prefix133 (.mul (variableData 3) 101,claim133)=true := by rfl
theorem step_checked134 : checkEntry inputs 4 prefix134 (.linear [(3,195), (7,22), (21,166), (22,41), (24,4), (39,228), (43,141), (44,19), (61,147), (93,173), (95,137), (97,40), (99,236), (101,51), (102,187), (104,97), (106,94), (110,208), (112,182), (122,207), (124,155), (126,42), (128,128), (132,209), (133,21)],claim134)=true := by rfl
theorem step_checked135 : checkEntry inputs 4 prefix135 (.mul (variableData 5) 101,claim135)=true := by rfl
theorem step_checked136 : checkEntry inputs 4 prefix136 (.linear [(5,52), (30,111), (31,122), (32,117), (74,78), (75,201), (104,50), (106,163), (112,47), (122,37), (124,174), (126,229), (128,167), (130,196), (132,82), (134,21), (135,28)],claim136)=true := by rfl
theorem step_checked137 : checkEntry inputs 4 prefix137 (.mul (variableData 5) 104,claim137)=true := by rfl
theorem step_checked138 : checkEntry inputs 4 prefix138 (.linear [(5,137), (31,11), (32,162), (33,138), (34,105), (76,224), (77,6), (78,147), (79,191), (104,211), (106,67), (128,28), (130,26), (132,34), (137,33)],claim138)=true := by rfl
theorem step_checked139 : checkEntry inputs 4 prefix139 (.mul (variableData 2) 112,claim139)=true := by rfl
theorem step_checked140 : checkEntry inputs 4 prefix140 (.linear [(2,80), (15,164), (16,169), (17,127), (18,190), (19,31), (48,159), (50,17), (51,70), (52,210), (93,252), (95,93), (112,18), (122,42), (124,176), (128,176), (132,82), (136,112), (139,98)],claim140)=true := by rfl
theorem step_checked141 : checkEntry inputs 4 prefix141 (.mul (variableData 3) 112,claim141)=true := by rfl
theorem step_checked142 : checkEntry inputs 4 prefix142 (.linear [(3,212), (7,80), (21,152), (22,97), (24,78), (39,103), (43,209), (44,121), (61,41), (93,249), (95,212), (97,133), (101,41), (102,255), (104,91), (106,177), (110,100), (112,155), (118,32), (122,169), (124,17), (126,104), (128,57), (132,66), (134,157), (136,3), (138,129), (140,37), (141,2)],claim142)=true := by rfl
theorem step_checked143 : checkEntry inputs 4 prefix143 (.mul (variableData 5) 112,claim143)=true := by rfl
theorem step_checked144 : checkEntry inputs 4 prefix144 (.linear [(5,47), (30,8), (31,200), (32,128), (74,242), (75,231), (104,112), (106,214), (112,130), (122,75), (124,228), (128,79), (130,179), (132,115), (134,236), (136,162), (138,214), (140,104), (142,27), (143,15)],claim144)=true := by rfl
theorem step_checked145 : checkEntry inputs 4 prefix145 (.mul (variableData 0) 122,claim145)=true := by rfl
theorem step_checked146 : checkEntry inputs 4 prefix146 (.linear [(0,124), (7,241), (8,249), (11,242), (36,94), (38,53), (39,149), (41,110), (43,212), (82,110), (87,115), (91,138), (93,138), (95,80), (99,22), (101,99), (108,217), (110,141), (112,92), (116,38), (122,51), (124,235), (128,127), (132,249), (136,221), (140,243), (142,246), (144,241), (145,40)],claim146)=true := by rfl
theorem step_checked147 : checkEntry inputs 4 prefix147 (.mul (variableData 4) 122,claim147)=true := by rfl
theorem step_checked148 : checkEntry inputs 4 prefix148 (.linear [(4,131), (27,154), (28,48), (29,192), (43,137), (69,47), (71,244), (95,108), (97,114), (101,215), (104,88), (106,137), (112,140), (122,112), (124,29), (128,187), (132,194), (134,71), (136,16), (138,45), (140,171), (142,105), (144,163), (146,191), (147,210)],claim148)=true := by rfl
theorem step_checked149 : checkEntry inputs 4 prefix149 (.mul (variableData 5) 122,claim149)=true := by rfl
theorem step_checked150 : checkEntry inputs 4 prefix150 (.linear [(5,232), (30,174), (31,208), (32,82), (75,176), (104,246), (106,206), (112,247), (122,14), (124,186), (128,71), (130,156), (132,250), (134,207), (136,140), (138,251), (140,231), (142,132), (144,66), (146,255), (148,168), (149,50)],claim150)=true := by rfl
theorem step_checked151 : checkEntry inputs 4 prefix151 (.mul (variableData 5) 124,claim151)=true := by rfl
theorem step_checked152 : checkEntry inputs 4 prefix152 (.linear [(5,217), (30,43), (31,46), (32,144), (75,34), (104,145), (106,82), (112,43), (122,94), (124,246), (128,210), (130,237), (132,90), (136,17), (138,56), (140,239), (142,97), (144,76), (146,232), (148,102), (150,134), (151,223)],claim152)=true := by rfl
theorem step_checked153 : checkEntry inputs 4 prefix153 (.mul (variableData 5) 128,claim153)=true := by rfl
theorem step_checked154 : checkEntry inputs 4 prefix154 (.linear [(5,178), (30,124), (31,256), (32,244), (33,176), (34,141), (78,248), (79,232), (80,253), (104,188), (106,80), (128,201), (130,82), (132,36), (136,207), (138,216), (140,77), (142,193), (144,91), (146,191), (148,222), (150,124), (152,111), (153,175)],claim154)=true := by rfl
theorem step_checked155 : checkEntry inputs 4 prefix155 (.mul (variableData 2) 132,claim155)=true := by rfl
theorem step_checked156 : checkEntry inputs 4 prefix156 (.linear [(2,36), (17,30), (18,8), (19,9), (53,108), (54,71), (93,168), (95,62), (112,124), (122,92), (124,100), (128,204), (130,145), (132,197), (136,177), (138,189), (140,200), (142,245), (144,203), (148,95), (150,42), (152,242), (154,225), (155,149)],claim156)=true := by rfl
theorem step_checked157 : checkEntry inputs 4 prefix157 (.mul (variableData 3) 132,claim157)=true := by rfl
theorem step_checked158 : checkEntry inputs 4 prefix158 (.linear [(3,42), (7,156), (21,10), (22,167), (24,214), (25,148), (43,112), (44,142), (62,190), (64,79), (95,83), (97,69), (104,24), (106,66), (128,117), (130,44), (132,63), (136,178), (138,154), (140,24), (142,182), (144,65), (146,122), (148,187), (150,143), (152,253), (154,173), (156,11), (157,67)],claim158)=true := by rfl
theorem step_checked159 : checkEntry inputs 4 prefix159 (.mul (variableData 2) 136,claim159)=true := by rfl
theorem step_checked160 : checkEntry inputs 4 prefix160 (.linear [(2,12), (15,24), (16,188), (17,137), (18,256), (19,212), (54,79), (93,9), (95,141), (112,2), (122,233), (124,112), (128,32), (130,136), (132,31), (136,136), (138,134), (140,52), (142,38), (144,143), (146,40), (148,216), (150,43), (152,212), (154,205), (156,48), (158,237), (159,178)],claim160)=true := by rfl
theorem step_checked161 : checkEntry inputs 4 prefix161 (.mul (variableData 3) 136,claim161)=true := by rfl
theorem step_checked162 : checkEntry inputs 4 prefix162 (.linear [(3,234), (7,81), (21,81), (22,205), (24,229), (25,247), (39,216), (43,157), (44,18), (64,196), (93,226), (95,126), (97,92), (101,76), (104,62), (112,175), (122,76), (124,148), (128,121), (130,197), (132,178), (136,238), (138,180), (140,45), (142,91), (144,199), (146,136), (148,170), (150,14), (152,207), (154,123), (156,11), (158,216), (160,5), (161,60)],claim162)=true := by rfl
theorem step_checked163 : checkEntry inputs 4 prefix163 (.mul (variableData 3) 140,claim163)=true := by rfl
theorem step_checked164 : checkEntry inputs 4 prefix164 (.linear [(3,48), (7,243), (21,252), (22,87), (24,210), (25,128), (39,241), (43,214), (44,63), (93,36), (95,120), (97,65), (101,61), (104,14), (106,150), (112,225), (122,155), (124,85), (128,138), (132,179), (136,89), (138,146), (140,223), (142,202), (144,75), (146,16), (148,28), (150,155), (152,244), (154,96), (156,246), (158,167), (160,116), (162,164), (163,16)],claim164)=true := by rfl
theorem step_checked165 : checkEntry inputs 4 prefix165 (.mul (variableData 5) 140,claim165)=true := by rfl
theorem step_checked166 : checkEntry inputs 4 prefix166 (.linear [(5,243), (30,164), (31,4), (32,68), (33,117), (34,17), (104,142), (106,234), (112,146), (122,104), (124,76), (128,16), (132,83), (136,205), (138,86), (140,34), (142,218), (144,6), (146,67), (148,193), (150,104), (152,252), (154,103), (156,181), (158,114), (160,255), (162,46), (164,252), (165,111)],claim166)=true := by rfl
theorem step_checked167 : checkEntry inputs 4 prefix167 (.mul (variableData 5) 142,claim167)=true := by rfl
theorem step_checked168 : checkEntry inputs 4 prefix168 (.linear [(5,175), (30,149), (31,125), (32,119), (33,241), (34,63), (104,172), (106,1), (124,82), (128,179), (132,8), (136,107), (138,219), (140,218), (142,130), (144,61), (146,74), (148,105), (150,170), (152,128), (154,16), (156,159), (158,167), (160,210), (162,20), (164,112), (166,217), (167,108)],claim168)=true := by rfl
theorem step_checked169 : checkEntry inputs 4 prefix169 (.mul (variableData 5) 146,claim169)=true := by rfl
theorem step_checked170 : checkEntry inputs 4 prefix170 (.linear [(5,158), (31,27), (32,71), (33,191), (34,208), (104,151), (106,190), (128,8), (132,79), (136,102), (138,144), (140,213), (142,22), (144,43), (146,54), (148,97), (150,176), (152,244), (154,102), (156,225), (158,182), (160,188), (162,95), (164,212), (166,46), (168,95), (169,106)],claim170)=true := by rfl
theorem step_checked171 : checkEntry inputs 4 prefix171 (.mul (variableData 5) 148,claim171)=true := by rfl
theorem step_checked172 : checkEntry inputs 4 prefix172 (.linear [(5,223), (31,155), (32,233), (33,224), (34,88), (106,117), (128,108), (132,218), (136,181), (138,174), (140,189), (142,239), (144,181), (146,87), (148,199), (150,24), (152,9), (154,4), (156,174), (158,192), (160,70), (162,184), (164,69), (166,173), (168,168), (170,6), (171,140)],claim172)=true := by rfl
theorem step_checked173 : checkEntry inputs 4 prefix173 (.mul (variableData 3) 156,claim173)=true := by rfl
theorem step_checked174 : checkEntry inputs 4 prefix174 (.linear [(3,93), (7,87), (24,70), (25,6), (43,102), (132,185), (136,133), (138,16), (140,119), (142,92), (144,148), (146,145), (148,21), (150,34), (152,254), (154,94), (156,199), (158,86), (160,120), (162,126), (164,245), (166,39), (168,5), (170,21), (172,231), (173,72)],claim174)=true := by rfl
theorem step_checked175 : checkEntry inputs 4 prefix175 (.mul (variableData 5) 164,claim175)=true := by rfl
theorem step_checked176 : checkEntry inputs 4 prefix176 (.linear [(5,191), (142,18), (144,185), (146,42), (148,40), (150,162), (152,43), (154,211), (156,149), (158,31), (160,206), (162,65), (164,89), (166,27), (168,225), (170,27), (172,102), (174,3), (175,239)],claim176)=true := by rfl
theorem all_steps (i : Fin 177) : checkEntry inputs 4 ((program.take i.val).map Prod.snd) (program.get i)=true := by
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
  · exact step_checked36
  · exact step_checked37
  · exact step_checked38
  · exact step_checked39
  · exact step_checked40
  · exact step_checked41
  · exact step_checked42
  · exact step_checked43
  · exact step_checked44
  · exact step_checked45
  · exact step_checked46
  · exact step_checked47
  · exact step_checked48
  · exact step_checked49
  · exact step_checked50
  · exact step_checked51
  · exact step_checked52
  · exact step_checked53
  · exact step_checked54
  · exact step_checked55
  · exact step_checked56
  · exact step_checked57
  · exact step_checked58
  · exact step_checked59
  · exact step_checked60
  · exact step_checked61
  · exact step_checked62
  · exact step_checked63
  · exact step_checked64
  · exact step_checked65
  · exact step_checked66
  · exact step_checked67
  · exact step_checked68
  · exact step_checked69
  · exact step_checked70
  · exact step_checked71
  · exact step_checked72
  · exact step_checked73
  · exact step_checked74
  · exact step_checked75
  · exact step_checked76
  · exact step_checked77
  · exact step_checked78
  · exact step_checked79
  · exact step_checked80
  · exact step_checked81
  · exact step_checked82
  · exact step_checked83
  · exact step_checked84
  · exact step_checked85
  · exact step_checked86
  · exact step_checked87
  · exact step_checked88
  · exact step_checked89
  · exact step_checked90
  · exact step_checked91
  · exact step_checked92
  · exact step_checked93
  · exact step_checked94
  · exact step_checked95
  · exact step_checked96
  · exact step_checked97
  · exact step_checked98
  · exact step_checked99
  · exact step_checked100
  · exact step_checked101
  · exact step_checked102
  · exact step_checked103
  · exact step_checked104
  · exact step_checked105
  · exact step_checked106
  · exact step_checked107
  · exact step_checked108
  · exact step_checked109
  · exact step_checked110
  · exact step_checked111
  · exact step_checked112
  · exact step_checked113
  · exact step_checked114
  · exact step_checked115
  · exact step_checked116
  · exact step_checked117
  · exact step_checked118
  · exact step_checked119
  · exact step_checked120
  · exact step_checked121
  · exact step_checked122
  · exact step_checked123
  · exact step_checked124
  · exact step_checked125
  · exact step_checked126
  · exact step_checked127
  · exact step_checked128
  · exact step_checked129
  · exact step_checked130
  · exact step_checked131
  · exact step_checked132
  · exact step_checked133
  · exact step_checked134
  · exact step_checked135
  · exact step_checked136
  · exact step_checked137
  · exact step_checked138
  · exact step_checked139
  · exact step_checked140
  · exact step_checked141
  · exact step_checked142
  · exact step_checked143
  · exact step_checked144
  · exact step_checked145
  · exact step_checked146
  · exact step_checked147
  · exact step_checked148
  · exact step_checked149
  · exact step_checked150
  · exact step_checked151
  · exact step_checked152
  · exact step_checked153
  · exact step_checked154
  · exact step_checked155
  · exact step_checked156
  · exact step_checked157
  · exact step_checked158
  · exact step_checked159
  · exact step_checked160
  · exact step_checked161
  · exact step_checked162
  · exact step_checked163
  · exact step_checked164
  · exact step_checked165
  · exact step_checked166
  · exact step_checked167
  · exact step_checked168
  · exact step_checked169
  · exact step_checked170
  · exact step_checked171
  · exact step_checked172
  · exact step_checked173
  · exact step_checked174
  · exact step_checked175
  · exact step_checked176
theorem checked : checkClaimedUnit inputs 4 program=true := by
  have all := checkEntriesFrom_of_each inputs 4 program [] all_steps
  change (checkEntriesFrom inputs 4 [] program && (add claim176 (scale (-1) SparseCertificate.oneData)).isEmpty)=true
  rw [all]
  rfl
theorem one_at_four : InClosure {p | p∈inputs.map meaning} 4 (1 : MvPolynomial (Fin 6) K) :=
  checkClaimedUnit_sound inputs 4 program checked
theorem exact_lfd_four : IsLastFallDegree {p | p∈inputs.map meaning} 4 := by
  apply exactLFD_of_first_one (by omega) one_at_four
  exact GeneratedPrimeLower.excludes_one
end ICFallDegrees.GeneratedPrimeUpper
