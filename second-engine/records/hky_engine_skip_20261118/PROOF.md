# The lead-cover criterion: which products the hky3 closure may skip, and why the space stays V_D

Date 2026-11-18. Engine: `src/cryptanalysis/hky_fp_tower_fast2.rs`, option `Fast2Options::skip`
(`--hky3-skip on`; default **off**), with `skip_order` (κ, below) and `skip_witness` (the witness
spaces of §2.3). Every number elsewhere in this directory is `[pcc:unverified]`; this file is a proof,
not a measurement. No literature was available in this sandbox (research-kb, the paper collection
and ecdb are not reachable here); the argument below is self-contained. Its shape is that of the
F5 "rewritten" idea (x_b·(x_a·s) = x_a·(x_b·s), keep one ordering), but the statement and the proof
are written for the degree-truncated linear closure with degree falls, where F5's signature order
does not apply directly; it is **not** a Gröbner-basis criterion and claims nothing beyond this engine.

## 1. Setting and notation (exactly what the engine computes)

**Ring.** R = F_p[y_0..y_{n−1}, z_1..z_f]/(y_j² − a_j y_j y_{nx(j)} − b_j y_{nx(j)} − c_j y_j − d_j), the
tower rules of `HkyRing` (`nx(j)` the next variable of the chain, absent at a chain's top, where
a_j = b_j = 0). Under the lex order with y_0 > y_1 > … (every `nx(j) > j`) the leading monomial of each
relation is y_j², and these are pairwise coprime, so the relations are a Gröbner basis (Buchberger's
first criterion): the monomials square-free in the y's form a basis of R, and `mul_tower_var` (the
reference's normal form of y_j·m) is the multiplication of the commutative ring R. For the paper
instances (`--instance`: Kummer tower, scale 1, top constant 1) a_j = c_j = 0, b_j = 1 below the top,
d = 1 at the top: R = F_p[(Z/2^k)²], the group algebra, and y_{i,j} = X_i^{2^j}.

For a variable v let T_v : R → R be multiplication by v. Two facts are used:

- **(C) commutation**: T_a T_b = T_b T_a (R is commutative). Unit test `variable_products_commute`
  checks it exhaustively on every monomial for t = 3, two blocks, Kummer and general rules.
- **(G) degree**: deg T_v x ≤ deg x + 1, deg = number of variables of a monomial (`key_degree`;
  free variables count their exponent). From `mul_tower_var`: if y_j ∤ m the result is y_j m
  (degree +1); otherwise it is m/y_j · (a y_j y_nx + b y_nx + c y_j + d), whose terms are m·y_nx
  (recursively ≤ deg m + 1), (m/y_j)·y_nx (≤ deg m), m, m/y_j. Same unit test.

**Columns and leads.** M_{≤d} = span of the basis monomials of degree ≤ d. The engine's columns are
sorted by key, and the key's highest bits are the degree, so the columns of higher degree come
first. A vector's *lead* is its first non-zero column; hence lead(x) has degree deg x, and for
x ∈ M_{≤d} with non-zero degree-d part top_d(x), lead(x) = lead(top_d(x)).

**Symmetry (split mode).** σ is an involution of the tower variables proved by `check_symmetry` to
map every rule to the rule of the image (so σ extends to a ring automorphism, degree preserving,
σ T_v σ = T_{σv}) and every generator to ± itself. M_{≤c} = M⁺ ⊕ M⁻ (σ-eigenspaces, both
M_{≤d}-graded), the engine's parts are these eigenspaces in orbit coordinates, and π_± = (1 ± σ)/2.
Products are formed only for the representatives v with σ(v) ≥ v (`Space::applies`); `multiply`
returns 2π₊(T_v x) and 2π₋(T_v x) (p odd, so the factor 2 changes neither spans nor leads). Without
the symmetry there is one part, σ = id, π = id, every variable is a representative. Write ε for a
part.

**The closure.** V = V_{F,c} is the least subspace of M_{≤c} containing the generators (of degree ≤ c)
and closed under x ↦ T_v x for every variable v and every x ∈ V with deg x < c. (G) puts T_v x in
M_{≤c}. V is σ-invariant (σV is the least such space for σF = ±F, i.e. V).

**What the engine maintains** (unchanged by this round except where it skips): an echelon W =
W₊ ⊕ W₋ (rows with distinct leads per part, monic, each row with a fixed lead; interreduction and
the second pass replace a row by a *version*: the same lead plus a combination of rows with later
leads). Generator rows are inserted first. When a row r of degree < c is inserted, the items
(r, v) for every representative v are queued; serving an item forms both projections of T_v r′
(r′ the version at that time), reduces them against W and inserts what survives. The run ends
when the queue is empty, or when 1 becomes a pivot (then V = M_{≤c}: 1 ∈ V and products of 1 by
distinct chain variables reach every monomial; unchanged by this round). Let W be the final span,
U = W ∩ M_{<c}, U_{≤d} = W ∩ M_{≤d}, U_{ε,≤d} = W_ε ∩ M_{≤d}. Since every version of every row lies
in W and its later-lead part lies in M_{≤deg r},

> U_{≤d} = U_{+,≤d} ⊕ U_{−,≤d}, and the final rows of part ε of degree ≤ d are a basis of U_{ε,≤d}. (E)

## 2. The criterion

### 2.1 Order and covers

Fix an injective order κ: representatives → {0, …, m−1} (`skip_order` 1, the default: κ decreasing
in the variable index; 0: increasing). For a part ε and a column ℓ of degree d (1 ≤ d ≤ c − 1) and a
value κ₀, say **ℓ is covered with κ₀ in part ε** when the engine holds a vector w with

> w = Σ_i λ_i · π_ε T_{a_i} s_i,  s_i ∈ W (versions of rows) of degree d − 1, κ(a_i) ≤ κ₀, lead(w) = ℓ. (W)

(d − 1 ≤ c − 2 is then automatic.) The engine records the least such κ₀ seen so far in
`cover[ε][ℓ]` (never increased). Two sources of such w, both implemented:

1. *Formed products.* When an item (s, a) with deg s = d − 1 ≤ c − 2 is formed and a projection
   π_ε T_a s′ has its lead ℓ in degree d, that projection is a w (one term).
2. *Witness spaces* (`skip_witness`, §2.3): when a row s of degree d − 1 ≤ c − 2 is inserted, the
   degree-d parts h_{a,ε} = top_d(π_ε T_a s) for every representative a are inserted with label κ(a)
   into a labelled semi-echelon (`Witness`); a row placed at column ℓ with label κ₀ is
   Σ λ_i h_{a_i,ε}(s_i) with κ(a_i) ≤ κ₀ and lead ℓ. Its preimage w = Σ λ_i π_ε T_{a_i} s_i has
   top_d(w) equal to that row (top_d is linear, and top_d(π_ε T_a s) only depends on top_{d−1}(s) by (G)),
   so lead(w) = ℓ: (W) holds.

### 2.2 The rule

> **Skip.** When an item (r, b) is served, r a row of part ε with lead ℓ, it is dropped (not formed,
> counted in `skipped_items`) iff ℓ is covered in part ε with some κ₀ < κ(b).

Nothing else changes: every other item is formed and reduced as before; generator rows are never
skipped (they are not items); rows of degree c have no items. Faults 6 (also skip when κ₀ = κ(b)),
7 (covers from products whose lead has any degree, and full products in the witness spaces) and 8
(skip whenever ℓ is covered at all) are negative controls outside the rule.

### 2.3 The witness structure keeps (W)

Inserting v with label k (sequentially) reduces v at its lead only: against the row at that column if its label is
≤ k (the result stays a combination of vectors of label ≤ k); if the row there has a label > k, v
(normalised) takes the column with label k and the old row minus v continues with the old label
(still a combination of vectors of label ≤ the old label, with its lead removed). By induction on
the operations, every row ever placed at column ℓ with label k is a combination of inserted vectors
of label ≤ k whose lead is ℓ; that is all (W) needs. (The structure is moreover a basis of each
E_κ = span of the vectors of label ≤ κ, so a column gets the least κ whose E_κ has it as a lead; that
is only why it covers as much as it does.)

The engine inserts each batch's vectors by parallel rounds (`Witness::insert_batch`), with the same
three operations, so the same invariant: (a) `prereduce(v, k)` subtracts from v, at its current lead,
the row there while that row's label is ≤ k (v stays a combination of vectors of label ≤ k); it stops
at a free column or at a row of label > k; (b) per leading column one survivor is placed (the
smallest label): at a free column as a new row with its own label and lead, or in place of a row of
larger label ol, which is replaced and goes on as `old − new` (both monic at that column, so the lead
is removed) with label ol — a combination of vectors of label ≤ ol; (c) the others return to (a).
Two economies only make fewer covers and are therefore sound by the same argument: vectors are
built only for variables with κ(a) < `skip_witness_max` (when set), and inserting them in increasing
label order. (A first, sequential version of this structure gave covers identical to the batched one
on the t = 10 c = 6 proxy: same `skipped_items`, same muladds — README §3.)

## 3. Theorem

**Theorem.** With the skip on, if the run ends with an empty queue and 1 not a pivot, then W = V
(hence the same `dims` = profile and `contains_one` = false as without the skip). If 1 becomes a pivot,
V = M_{≤c}. In either case the output equals that of the run without the skip.

*Proof.*

**(i) W ⊆ V.** Unchanged by skipping: every inserted vector is a generator row or a reduced
projection π_ε T_v x with x ∈ W ∩ M_{<c}; by induction on insertions W ⊆ V (V is closed, σ-invariant,
so π_ε T_v x = (T_v x + ε T_{σv} σx)/2 ∈ V). In particular if 1 becomes a pivot, 1 ∈ V and V = M_{≤c}.
Conversely, if 1 ∈ V then (by (ii) below) 1 ∈ W, so 1 is a pivot of W₊ (the only element of W
supported on the last column is a multiple of 1): the run cannot end with an empty queue without
reporting 1. So the two runs report the same `contains_one`, and when it is true both report M_{≤c}.

**(ii) V ⊆ W when the queue empties.** W contains the generators. It suffices that
T_v U ⊆ W for every variable v; then W is a subspace of M_{≤c} with the defining closure property,
so V ⊆ W by minimality. We show it for the representatives in increasing κ, and get the others from
σ.

*Induction hypothesis* H(κ₁): T_a U ⊆ W and T_{σa} U ⊆ W for every representative a with κ(a) < κ₁.
H(0) is empty. Fix b with κ(b) = κ₁ and assume H(κ₁). Let

> K = { x ∈ U : T_b x ∈ W },

a subspace of U. We show U_{≤d} ⊆ K for d = −1, 0, …, c − 1 by induction on d (U_{≤−1} = 0).

Assume U_{≤d−1} ⊆ K. Fix a part ε. For every final row r of part ε with deg r ≤ d we name an element
w_r ∈ K ∩ U_{ε,≤d} with lead(w_r) = lead(r):

1. **deg r < d.** w_r = r ∈ U_{≤d−1} ⊆ K.
2. **deg r = d and the item (r, b) was formed** (from the version r′ then current). Both projections of
   T_b r′ were reduced into W, so T_b r′ = π₊T_b r′ + π₋T_b r′ ∈ W: r′ ∈ K. r′ ∈ U_{ε,≤d} and
   lead(r′) = lead(r). Take w_r = r′.
3. **deg r = d and (r, b) was skipped.** Then lead(r) = ℓ was covered with some κ₀ < κ(b): there is
   w = Σ λ_i π_ε T_{a_i} s_i as in (W), κ(a_i) ≤ κ₀ < κ(b), deg s_i = d − 1, lead(w) = ℓ. We check
   w ∈ K ∩ U_{ε,≤d}:
   - Each s_i ∈ W ∩ M_{≤d−1} = U_{≤d−1} ⊆ K (d − 1 ≤ c − 2), so T_b s_i ∈ W; by (G) T_b s_i ∈ M_{≤d} and
     d ≤ c − 1, so T_b s_i ∈ U. s_i lies in one part ε_i, so σ s_i = ε_i s_i and T_b σ s_i = ε_i T_b s_i ∈ U.
   - w ∈ U: T_{a_i} s_i ∈ W by H(κ₁) (s_i ∈ U, κ(a_i) < κ₁), and π_ε maps W into W (W = W₊ ⊕ W₋), and
     by (G) w ∈ M_{≤d} ⊆ M_{<c}; it lies in part ε.
   - T_b w ∈ W: π_ε T_a s = ½(T_a s + ε σ T_a s) = ½(T_a s + ε T_{σa} σ s), so by (C)
     T_b π_ε T_a s = ½(T_a (T_b s) + ε T_{σa}(T_b σ s)), and both terms are in W by H(κ₁) applied to
     T_b s, T_b σs ∈ U (the hypothesis covers a and σa). Summing over i: T_b w ∈ W.
   Take w_r = w.

The w_r have pairwise distinct leads, so they are linearly independent; there are as many as final
rows of part ε of degree ≤ d, i.e. dim U_{ε,≤d} by (E). They lie in the subspace K ∩ U_{ε,≤d}, so
U_{ε,≤d} ⊆ K. With both parts and (E), U_{≤d} ⊆ K. This ends the induction on d: U = U_{≤c−1} ⊆ K,
i.e. T_b U ⊆ W.

Then T_{σb} U = σ T_b σ U = σ T_b U ⊆ σ W = W (σU = U and σW = W: W = W₊ ⊕ W₋ with W_± in the
σ-eigenspaces, and σ preserves degree). So H(κ₁ + 1) holds. After the last representative,
T_v U ⊆ W for all v (every variable is a representative or the image of one). ∎

### Remarks on the proof (what is and is not used)

- **Truncation at degree c.** Used exactly in three places: the witnesses come from rows of degree
  ≤ c − 2 (so T_b s ∈ U, which (G) guarantees), K lives in U = W ∩ M_{<c}, and the covered row r has
  degree d ≤ c − 1. Products of rows of degree c − 1 *can* be skipped (their lead may be covered by
  products of degree-(c − 2) rows); products that leave M_{≤c} are never needed.
- **Degree falls / mutants.** No assumption on how a row arose: a row r of degree d obtained by a
  fall (a product of degree d + 1 whose top cancelled) is skipped only if its *lead* is covered by
  products of degree-(d − 1) rows; whether r itself came from such a product is irrelevant. The
  proof never follows a row's history, only the final echelon (E) and the cover certificates (W).
- **Normal form by the chain relations.** Enters only through (C) and (G), both properties of the
  multiplication the engine uses, and through σ being a ring automorphism (`check_symmetry`).
- **Versions, interreduction, epochs, sweep, batches, threads.** Not used: case 2 takes whatever
  version was multiplied; witnesses may use any version (all are in W); the schedule may change
  arbitrarily (it does: counters differ), the final W is V.
- **Timing.** A cover is used only after it is recorded, and covers only decrease, so every skip is
  backed by a certificate that exists at the end of the run.
- **Why κ₀ < κ(b) strictly.** Case 3 needs closure under T_{a_i} *before* b. With κ₀ = κ(b), i.e.
  a_i = b, T_b w needs T_b(T_b s) ∈ W, which is the closure being proved (fault 6). With no order at
  all (fault 8) the argument is circular. With a witness whose lead is of degree < deg s + 1
  (fault 7), s ∉ U_{≤d−1} at the stage d of the lead, so s ∈ K is not yet known.
- **Exactness.** Nothing is sampled; no probabilistic step.

## 4. What the gate checks against this proof

The proof is the reason to trust the skip; the gate (README §4) checks the implementation:
equality of dims/profile/contains_one with the definitional reference on every gate cell, the
independent closedness check (`--hky3-verify-closed on`: every variable times every final row of
degree < c and every generator reduces to zero, computed without the skip, the epochs or the sweep)
on every closure without 1, and the negative controls (faults 6–8) to show the checks can see a
skip the proof does not cover.
