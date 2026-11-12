/-!
# Kernel-checkable base certificate for N5 (checker functions)

Computable checker for the three per-instance hypotheses of
`ChainModelMap.n5_lfd_ge_five` / `ChainModelMap.not_one_of_base` with `k0 = 7`, `D = 4`:
(B1) `1 ∉ S`, (B2) `W_8(S)` closed, (G) the generator lies in `S`, where `S` is the row space
of a reduced echelon basis of the depth-7 degree-4 closure `V_4(7)` of one instance.

Everything here is plain `Nat`/`List` code with structural recursion only (no well-founded
recursion, no `Array`, no `String`), so that the Lean kernel can evaluate it by `decide +kernel`.
The soundness proofs are in `N5CertSound`; the instance data in `N5CertData`.

## Encoding

A monomial of the depth-7 two-summand algebra is the *code* `m = l₀ + 128 l₁` (`l₀, l₁ < 128`),
read as the low exponent vector `![m % 128, m / 128]` (`N5CertSound.col`).  The reduced echelon
basis has `nRows` rows; its pivot columns and free columns are codes.

* `pivEnc`  : 16-bit slots, slot `idx` is the pivot code of row `idx`;
* `freeEnc` : 16-bit slots, slot `f` is the code of free column `f` (`f < nFree`);
* `tbl`     : 16-bit slots indexed by code `m < 16384`: `idx+1` if `m` is the pivot of row `idx`,
  `1024 + f` if `m` is free column `f`, `0` otherwise (a hint; every use is re-verified);
* `rows`    : a binary tree whose entry `idx` is row `idx` *packed*: free column `f` carries
  the coefficient in bits `[w f, w (f+1))`, `B = 2^w`;
* `low`     : the rows whose pivot has weight `≤ 3`, each with its free-column entries as a
  sparse list `(code, coefficient)` (consistency with the packed row is checked);
* `gen`     : the generator as a sparse list.

Membership of a sparse vector `u` in the row space is the echelon test
`u = ∑_{pivot entries (m,x) of u} x · row_m`: on the free columns this is checked by
accumulating `∑ x · packedRow + ∑_{free entries} (p − x) · B^f` in one natural number (no slot
overflows, since all summands are `< p²` and `|u| ≤ 512`, and `512 p² < B`) and testing every
slot for divisibility by `p`.
-/

namespace ICFallDegrees.N5Cert

/-- Binary tree of packed rows, indexed by `Tree.get`. -/
inductive Tree where
  | leaf : Tree
  | node : Tree → Nat → Tree → Tree

/-- Heap-style indexing: `0` is the root, index `i > 0` goes left if odd, right if even. -/
def Tree.get : Tree → Nat → Nat
  | .leaf, _ => 0
  | .node l v r, i => if i = 0 then v else if i % 2 = 1 then l.get (i / 2) else r.get (i / 2 - 1)

abbrev Entry := Nat × Nat

structure Data where
  p : Nat
  w : Nat
  nRows : Nat
  nFree : Nat
  pivEnc : Nat
  freeEnc : Nat
  tbl : Nat
  rows : Tree
  low : List (Nat × List Entry)
  gen : List Entry

/-- 16-bit slot `i` of a packed natural. -/
def get16 (enc i : Nat) : Nat := (enc >>> (16 * i)) &&& 65535

def Data.piv (d : Data) (i : Nat) : Nat := get16 d.pivEnc i
def Data.freeCode (d : Data) (f : Nat) : Nat := get16 d.freeEnc f
def Data.tblGet (d : Data) (m : Nat) : Nat := get16 d.tbl m
def Data.B (d : Data) : Nat := 2 ^ d.w
def Data.row (d : Data) (i : Nat) : Nat := d.rows.get i

/-- Popcount with fuel (structural). -/
def popcount : Nat → Nat → Nat
  | 0, _ => 0
  | fuel+1, n => n % 2 + popcount fuel (n / 2)

/-- `∀ i < n, P i`, as a Bool loop. -/
def allBelow (P : Nat → Bool) : Nat → Bool
  | 0 => true
  | n+1 => P n && allBelow P n

/-- `∀ i, lo ≤ i < lo + n → P i`. -/
def rangeAll (P : Nat → Bool) (lo : Nat) : Nat → Bool
  | 0 => true
  | n+1 => P (lo + n) && rangeAll P lo n

/-- The `n` slots of a packed natural (least significant first). -/
def unpack (B : Nat) : Nat → Nat → List Nat
  | 0, _ => []
  | n+1, R => (R % B) :: unpack B n (R / B)

/-- Every one of the `n` slots is divisible by `p`, and nothing is left above them. -/
def slotsZero (B p : Nat) : Nat → Nat → Bool
  | 0, R => R == 0
  | n+1, R => (R % B % p == 0) && slotsZero B p n (R / B)

/-- Every one of the `n` slots is `< bound`, and nothing is left above them. -/
def slotsLt (B bound : Nat) : Nat → Nat → Bool
  | 0, R => R == 0
  | n+1, R => (R % B < bound) && slotsLt B bound n (R / B)

/-- A sparse entry is admissible: `0 < x < p` and its code is a verified pivot or free column. -/
def entryOK (d : Data) (e : Entry) : Bool :=
  0 < e.2 && e.2 < d.p && e.1 < 16384 &&
    (match d.tblGet e.1 with
      | 0 => false
      | v+1 =>
        if v < 1024 then v < d.nRows && d.piv v == e.1
        else v + 1 - 1024 < d.nFree && d.freeCode (v + 1 - 1024) == e.1)

/-- A free-column entry: `0 < x < p`, code verified as a free column. -/
def freeEntryOK (d : Data) (e : Entry) : Bool :=
  0 < e.2 && e.2 < d.p && e.1 < 16384 &&
    (match d.tblGet e.1 with
      | 0 => false
      | v+1 => 1024 ≤ v && v + 1 - 1024 < d.nFree && d.freeCode (v + 1 - 1024) == e.1)

/-- The echelon accumulator: `∑_{pivot entries} x · row + ∑_{free entries} (p − x) · B^f`. -/
def accum (d : Data) : List Entry → Nat
  | [] => 0
  | e :: rest =>
    (match d.tblGet e.1 with
      | 0 => 0
      | v+1 => if v < 1024 then e.2 * d.row v else (d.p - e.2) <<< (d.w * (v + 1 - 1024))) +
      accum d rest

/-- Pack free-column entries: `∑ c · B^f`. -/
def packFree (d : Data) : List Entry → Nat
  | [] => 0
  | e :: rest =>
    (match d.tblGet e.1 with
      | 0 => 0
      | v+1 => e.2 <<< (d.w * (v + 1 - 1024))) + packFree d rest

/-- Membership of a sparse vector in the row space. -/
def member (d : Data) (u : List Entry) : Bool :=
  u.length ≤ 512 && u.all (entryOK d) && slotsZero d.B d.p d.nFree (accum d u)

/-- One translation `u_{i,j}` (`j < 7`) applied to a monomial code at depth 8: the no-carry
part (`carry = false`) or the carried part with the level-7 bit removed (`carry = true`). -/
def shiftEntry (i j : Nat) (carry : Bool) (e : Entry) : Option Entry :=
  let l0 := e.1 % 128
  let l1 := e.1 / 128
  if i = 0 then
    let s := l0 + 2 ^ j
    if s < 128 then (if carry then none else some (s + 128 * l1, e.2))
    else (if carry then some (s - 128 + 128 * l1, e.2) else none)
  else
    let s := l1 + 2 ^ j
    if s < 128 then (if carry then none else some (l0 + 128 * s, e.2))
    else (if carry then some (l0 + 128 * (s - 128), e.2) else none)

def shiftList (i j : Nat) (carry : Bool) (L : List Entry) : List Entry :=
  L.filterMap (shiftEntry i j carry)

/-- The full sparse row of a low row: pivot with coefficient `1`, then the free entries. -/
def lowFull (d : Data) (r : Nat × List Entry) : List Entry := (d.piv r.1, 1) :: r.2

/-- All `2 · 7 · 2` shifted vectors of one low row lie in the row space. -/
def shiftsOK (d : Data) (r : Nat × List Entry) : Bool :=
  allBelow (fun j => member d (shiftList 0 j false (lowFull d r)) &&
    member d (shiftList 0 j true (lowFull d r)) &&
    member d (shiftList 1 j false (lowFull d r)) &&
    member d (shiftList 1 j true (lowFull d r))) 7

/-- A low row is consistent with its packed row and passes all shift checks. -/
def lowRowOK (d : Data) (r : Nat × List Entry) : Bool :=
  r.1 < d.nRows && r.2.length ≤ 256 && r.2.all (freeEntryOK d) &&
    packFree d r.2 == d.row r.1 && shiftsOK d r

/-- Row `idx` is well formed: its pivot code is registered in the table, and its packed form
has all `nFree` slots `< p` and nothing above them. -/
def rowOK (d : Data) (idx : Nat) : Bool :=
  d.tblGet (d.piv idx) == idx + 1 && d.piv idx < 16384 && slotsLt d.B d.p d.nFree (d.row idx)

/-- Free column `f` is registered as free (hence is not a pivot) and is a valid code. -/
def freeOK (d : Data) (f : Nat) : Bool :=
  1024 ≤ d.tblGet (d.freeCode f) && d.freeCode f < 16384

/-- The low list is exactly the list of rows whose pivot has weight `≤ 3`, in order. -/
def lowIdxOK (d : Data) : Bool :=
  (List.range d.nRows).filter (fun i => popcount 16 (d.piv i) ≤ 3) == d.low.map Prod.fst

/-- Global sanity: the slot width is large enough for the accumulator bound `512 p² < B`,
the constant monomial (code 0) is a free column, and `p ≥ 3`. -/
def globalOK (d : Data) : Bool :=
  512 * d.p * d.p < d.B && 1024 ≤ d.tblGet 0 && 3 ≤ d.p

end ICFallDegrees.N5Cert
