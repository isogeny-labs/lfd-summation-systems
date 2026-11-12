import ICFallDegrees.N5Cert

/-!
# Kernel-checkable certificates for an exact last fall degree at depth 6 (checker functions)

Two computable checkers for one two-summand prime chain instance at depth `k = 6`
(`PrimeM2General.inputs 6 _ a b t` over `ZMod p`), both plain `Nat`/`List` code with structural
recursion only, so that the Lean kernel evaluates them by `decide +kernel`.  The soundness proofs
are in `ExactK6LowerSound` and `ExactK6UpperSound`.

## Lower half: `1 ∉ V_4`

The containing space is the row space `S` of a reduced echelon basis of the depth-6 degree-4
quotient closure `V_4(6)`, in the data format of `N5Cert.Data` (packed rows on the free columns,
membership by the echelon accumulator `N5Cert.member`).  A squarefree chain monomial with block
exponents `(e₀, e₁)`, `e₀, e₁ < 64`, has the *complement code* `m = 4095 - (e₀ + 64 e₁)`; the
complement makes every column of weight `≤ 4` a code of popcount `≥ 8`, so that the depth-7
specific low-row list of `N5Cert.Data` is empty and `N5CertSound.Valid` holds with `low = []`.
The rows of pivot weight `≤ 3` (popcount `≥ 9`) are given separately as sparse lists; for each
of them the `12` exact translations of depth 6 (`shiftCode`, with the wrap `2^6 ≡ 0`) are tested
for membership.

## Upper half: `1 ∈ V_5`

A derivation of `1` in the quotient closure of degree 5.  A vector of the depth-6 group algebra
is one natural number with `128 × 128` slots of width `w` (`B = 2^w`): slot `e₀ + 128 e₁` carries
the coefficient of the monomial with block exponents `(e₀ mod 64, e₁ mod 64)`.  The stored rows
occupy the slots with `e₀, e₁ < 64` only, so the translate of a row by `(s₀, s₁)`, `s₀, s₁ < 64`,
is the single shift `row <<< (w * (s₀ + 128 s₁))` with no wrap; the wrap is done by the final
scan, which folds the four quadrants.  Row `n` is certified by a list of terms
`(c, s, i)`, `i < n`, `wt(s) + deg(row i) ≤ 5`, with

  `∑ c · (row i <<< (w s))  ≡  row n   (mod p, after folding)`.

The accumulation costs three big-number operations per term whatever the density of the rows; the
scan reads every slot once.
-/

namespace ICFallDegrees.ExactK6Cert

open ICFallDegrees.N5Cert

/-! ## Lower half -/

/-- The translation `u_{i,j}` at depth 6 on complement codes `m = 4095 - (e₀ + 64 e₁)`:
`e_i ↦ (e_i + 2^j) mod 64`. -/
def shiftCode (i j m : Nat) : Nat :=
  if i = 0 then 4095 - (((4095 - m) % 64 + 2 ^ j) % 64 + 64 * ((4095 - m) / 64))
  else 4095 - ((4095 - m) % 64 + 64 * (((4095 - m) / 64 + 2 ^ j) % 64))

def shiftList (i j : Nat) (L : List Entry) : List Entry := L.map fun e => (shiftCode i j e.1, e.2)

/-- All `2 · 6` translates of one low row lie in the row space. -/
def shiftsOK (d : Data) (r : Nat × List Entry) : Bool :=
  allBelow (fun j => member d (shiftList 0 j (lowFull d r)) &&
    member d (shiftList 1 j (lowFull d r))) 6

/-- A low row is consistent with its packed row and passes all shift checks. -/
def lowRowOK (d : Data) (r : Nat × List Entry) : Bool :=
  r.1 < d.nRows && r.2.length ≤ 256 && r.2.all (freeEntryOK d) &&
    packFree d r.2 == d.row r.1 && shiftsOK d r

/-- The low list is exactly the list of rows whose pivot has weight `≤ 3` (popcount `≥ 9` of the
complement code), in order. -/
def lowIdxOK (d : Data) (low : List (Nat × List Entry)) : Bool :=
  (List.range d.nRows).filter (fun i => 9 ≤ popcount 16 (d.piv i)) == low.map Prod.fst

/-- All pivot and free codes are codes of the depth-6 box. -/
def boxOK (d : Data) : Bool :=
  allBelow (fun i => d.piv i < 4096) d.nRows && allBelow (fun f => d.freeCode f < 4096) d.nFree

/-- The constant monomial (complement code 4095) is registered as a free column. -/
def constOK (d : Data) : Bool := 1024 ≤ d.tblGet 4095

/-! ## Upper half -/

/-- A term `(c, s, i)`: coefficient, translation `s = s₀ + 128 s₁`, source row. -/
abbrev Term := Nat × Nat × Nat

structure UData where
  p : Nat
  w : Nat
  /-- the rows, packed with `128 × 128` slots of width `w` -/
  rows : Tree
  /-- 4-bit slots: the degree bound of row `n` -/
  degEnc : Nat
  /-- the generator as a sparse list `(slot, coefficient)` -/
  gen : List Entry

def UData.B (u : UData) : Nat := 2 ^ u.w
def UData.row (u : UData) (n : Nat) : Nat := u.rows.get n
def UData.deg (u : UData) (n : Nat) : Nat := (u.degEnc >>> (4 * n)) &&& 15

/-- `∑ c · (row i <<< (w s))`. -/
def accum (u : UData) : List Term → Nat
  | [] => 0
  | t :: r => t.1 * (u.row t.2.2 <<< (u.w * t.2.1)) + accum u r

def coefSum : List Term → Nat
  | [] => 0
  | t :: r => t.1 + coefSum r

/-- A term of row `n` is admissible: earlier source row, translation in the box, and the degree
budget `wt(s) + deg(row i) ≤ 5`. -/
def termOK (u : UData) (n : Nat) (t : Term) : Bool :=
  t.2.2 < n && t.2.1 % 128 < 64 && t.2.1 / 128 < 64 &&
    popcount 6 (t.2.1 % 128) + popcount 6 (t.2.1 / 128) + u.deg t.2.2 ≤ 5

/-- Scan of the 64 slots `e₀ = 0, …, 63` of one block `e₁` (`w1 = popcount e₁`): the four
quadrant slots `x1 … x4` of the accumulator fold to the slot of the row block `h` modulo `p`,
the row slot is `< p`, a nonzero row slot has weight `≤ d`, and nothing of `h` is left. -/
def scanSlots (B p d w1 : Nat) : Nat → Nat → Nat → Nat → Nat → Nat → Nat → Bool
  | 0, _, _, _, _, _, h => h == 0
  | n+1, e0, x1, x2, x3, x4, h =>
    ((x1 % B + x2 % B + x3 % B + x4 % B + (p - 1) * (h % B)) % p == 0) && decide (h % B < p) &&
      (h % B == 0 || decide (popcount 6 e0 + w1 ≤ d)) &&
      scanSlots B p d w1 n (e0 + 1) (x1 / B) (x2 / B) (x3 / B) (x4 / B) (h / B)

/-- Scan of the blocks `e₁ = 0, …, 63`: `lo`, `hi` are the low and high halves of the
accumulator, `H` the row. -/
def scanBlocks (B p d : Nat) : Nat → Nat → Nat → Nat → Nat → Bool
  | 0, _, _, hi, H => hi == 0 && H == 0
  | n+1, e1, lo, hi, H =>
    scanSlots B p d (popcount 6 e1) 64 0 (lo % B ^ 128 % B ^ 64) (lo % B ^ 128 / B ^ 64)
        (hi % B ^ 128 % B ^ 64) (hi % B ^ 128 / B ^ 64) (H % B ^ 128) &&
      scanBlocks B p d n (e1 + 1) (lo / B ^ 128) (hi / B ^ 128) (H / B ^ 128)

def rowScan (u : UData) (n acc : Nat) : Bool :=
  scanBlocks u.B u.p (u.deg n) 64 0 (acc % u.B ^ 8192) (acc / u.B ^ 8192) (u.row n)

/-- Row `n` is the folded reduction of the combination given by `terms`. -/
def rowCheck (u : UData) (n : Nat) (terms : List Term) : Bool :=
  terms.all (termOK u n) && decide (coefSum terms * u.p < u.B) && rowScan u n (accum u terms)

/-- `∑ x · B^c` over a sparse list. -/
def packEntries (w : Nat) : List Entry → Nat
  | [] => 0
  | e :: r => (e.2 <<< (w * e.1)) + packEntries w r

def entrySum : List Entry → Nat
  | [] => 0
  | e :: r => e.2 + entrySum r

/-- Row `0` is the packed generator list, and it is a well-formed row. -/
def genOK (u : UData) : Bool :=
  u.row 0 == packEntries u.w u.gen && decide (entrySum u.gen < u.B) &&
    u.gen.all (fun e => decide (e.1 < 16384)) && rowScan u 0 (u.row 0)

def uGlobalOK (u : UData) : Bool := decide (2 ≤ u.p) && decide (1 ≤ u.w)

end ICFallDegrees.ExactK6Cert
