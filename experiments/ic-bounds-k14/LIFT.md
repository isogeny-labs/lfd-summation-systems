# Degree lift through the quotient: V_D -> V_{D+1} without a closure in M_<=D+1

Labels: **proof** (written out below), **exact** (certified computation, command named), **empirical**,
**estimate**, **unverified**.  Nothing here has been run on a real GPU; see section 9.

| file | role |
|---|---|
| `csd/ic_bounds_lift.py` | the lift (imports `ic_bounds_big` from the same directory) |
| `scripts/ic_bounds_lift_run.py` | driver: `lift`, `plan`, `selftest` (no repository imports) |
| `scripts/ic_bounds_lift_validate.py` | local validation harness (oracle: the direct closure of `csd/ic_bounds_big.py`) |
| `experiments/ic-bounds-k14/validation/lift-*.json` | validation records |

## 1. Statement

Notation as in `README.md`: F = F_p (p odd), monomials (a1, a2) of weight popcount(a1) + popcount(a2), the 2k
shifts u = u_{i,j}; a shift maps a monomial of weight w to a monomial of weight <= w + 1.  L = span(M_<=D).
V_n is the least subspace of span(M_<=n) that contains g and is closed under f -> u f for f in V_n supported on
weight < n.  sigma swaps the summands, g is sigma-invariant, V_n = V_n^+ (+) V_n^- (parts q = +, -), and the
engine stores each part in orbit coordinates as a reduced echelon basis for the weight-descending column order.

**Definitions.**  A weight-D monomial is *regular* if its sigma-orbit is a pivot column of both parts of V_D (a
self-paired monomial: of the + part), otherwise *free-type*.  A weight-(D+1) monomial is *special* if one of its
D+1 weight-D divisors is free-type, otherwise *regular*.  Both notions are sigma-stable.  For a part q:
Reg_q, Sp_q = the regular / special weight-(D+1) columns of the part, Fr_q = the free columns of V_D^q,
T_q = F^{Sp_q} (+) F^{Fr_q}.  red_q : A^q -> F^{Reg_q} (+) T_q is the identity on the weight-(D+1) coordinates
and the residual against the basis of V_D^q on the others; it is onto with kernel V_D^q.

**Theorem (proof, section 2).**  Let V_D be the closed degree-D closure.  For each part q let K_q be the
subspace of T_q produced by the algorithm of section 3 when it stops ("fixed point"), and let F_q be the set
of elements of K_q with zero Sp-part (the *falls*), read as vectors on the free columns of V_D^q.  Then

1. V_{D+1}^q n L = V_D^q (+) F_q   (F_q as elements supported on the free columns);
2. dim V_{D+1}^q = dim V_D^q + |Reg_q| + dim K_q, and dim V_{D+1} is the sum over the two parts;
3. 1 is in V_{D+1} iff 1 is in V_D or the class of 1 (a free column of the + part) lies in F_+;
4. the reduced echelon basis of V_{D+1}^q n L for the degree-D column order has pivot set
   pivots(V_D^q) u pivots(F_q); it is obtained by adding the echelonised falls to the basis of V_D^q.

If the algorithm stops without 1 in F_+ (and 1 is not in V_D), that is a complete proof that 1 is not in V_{D+1}:
every generating vector listed in section 3 is processed, none is sampled.

## 2. Proof

**Lemma 1 (fixed point).**  Put W = V_{D+1} n L, W_0 = V_D, X_i = W_i + sum_u u W_i, W_{i+1} = X_i n L.  Then
W_i c W_{i+1} c W and X_i c V_{D+1} for all i, and if W_{i+1} = W_i then W_i = W and X_i = V_{D+1}.

*Proof.*  V_{D+1} n L contains g and is closed under the shifts of its elements of weight < D (they lie in
V_{D+1} and have weight <= D), so it contains the least such space V_D: W_0 c W.  If W_i c W then every element
of W_i has weight <= D < D + 1, so u W_i c V_{D+1}, X_i c V_{D+1} and W_{i+1} = X_i n L c W; W_i c X_i n L is
clear.  If W_{i+1} = W_i: X_i contains g, and an element of X_i supported on weight < D + 1 lies in
X_i n L = W_i, whose shifts lie in X_i.  So X_i is closed, hence contains V_{D+1}.  QED

**Lemma 2 (which rows).**  If V_D is closed then V_D + sum_u u V_D = V_D + span{u b : b a basis row of V_D with
pivot weight D, u a shift}.  *Proof.*  A basis row with pivot weight < D is supported on weight < D (column
order), so its shifts lie in V_D.  QED  (`--include-low` shifts the other rows too; the result must not change.)

**Lemma 3 (split).**  W_i and X_i are sigma-stable, X_i^q = W_i^q + sum_{j, p} pi_q u_{1,j} W_i^p with
pi_q = (1 + s_q sigma)/2, and W_{i+1}^q = X_i^q n L^q.  *Proof.*  sigma u_{1,j} sigma = u_{2,j}, so by induction
all spaces are sigma-stable; for w in W_i^p, pi_q u_{2,j} w = pi_q sigma u_{1,j} sigma w = +- pi_q u_{1,j} w.  QED
In orbit coordinates the engine's "lift, shift, project" map of part p into part q is 2 pi_q u_{1,j} on A^p.

**Lemma 4 (shape of the generators).**  Let [m] be a regular weight-D orbit with representative m (a1 <= a2),
b^+ and b^- its basis rows.  Put l_m = (lift b^+ + lift b^-)/2 and l_{sigma m} = (lift b^+ - lift b^-)/2
(self-paired: l_m = lift b^+).  Then l_x = e_x + phi(x) with phi(x) supported on monomials whose orbit is a free
column of at least one part, i.e. on free-type monomials of weight D and on monomials of weight < D.  For a
shift u_{1,j}:
G(x, j, q) := 2 pi_q u_{1,j} l_x has, in the coordinates of part q, at most one nonzero regular top
coordinate: the column of u_{1,j} x, with coefficient 1, -1 or 2 (2: self-paired top, + part), if u_{1,j} x has
weight D + 1, is regular and has a column in part q (a self-paired top has no column in the - part).  Every other
top coordinate is special, because it is a shift of a free-type monomial of weight D.  The images of the two
basis rows are y^+ = G(m) + G(sigma m), y^- = G(m) - G(sigma m).  A basis row of pivot weight D whose orbit is
free-type (pivot in exactly one part) and every element supported on the free columns of a part have special
top coordinates only.

**Lemma 5 (elimination).**  Let E = F^Reg (+) T and let S be a finite subset of E whose elements have at most
one nonzero Reg coordinate.  For r in Reg let n_1, ..., n_s (s = s_r) be the elements with nonzero r-coordinate,
scaled to r-coordinate 1.  Let K be spanned by the differences n_{i+1} - n_i (all r, 1 <= i < s_r) and by the
elements of S without Reg coordinate.  Then K c T, span S = (direct sum of F n_1^(r) over the r with s_r >= 1) (+) K,
and span S n T = K.  *Proof.*  Triangular: the n_1^(r) have distinct leading regular coordinates and K has none.  QED

**Proof of the theorem.**  By Lemma 1 and Lemma 3 it suffices to follow X_i^q modulo V_D^q.  Write
W_i^q = V_D^q (+) Phi_i^q with Phi_i^q supported on the free columns (the residuals of W_i^q; Phi_0 = 0).  By
Lemmas 2, 3 and 4, red_q(X_i^q) is spanned by S_0 u S_i', where S_0 consists of red_q G(x, j, q) for all regular
weight-D monomials x and all j, and of red_q y for the rows y of free-type pivot orbits, and S_i' of the red_q of
the images of a basis of Phi_i^p, p = +, -.  All of these have at most one regular coordinate, and S_i' has none.
For a regular column r of part q with representative M = (A1, A2), the elements of S_0 with nonzero r-coordinate
are exactly the G(x, j, q) with u_{1,j} x in {M, sigma M}: one for each of the D + 1 divisor variables of M when
A1 != A2 (x = M / u_{1,j} for j in A1, x = sigma M / u_{1,j} for j in A2), and popcount(A1) = (D+1)/2 when
A1 = A2; in particular s_r >= 1.  Lemma 5 gives red_q(X_i^q) = (reductors) (+) K_q^(i) with
K_q^(i) = K_q^(0) + span S_i', and W_{i+1}^q / V_D^q = red_q(X_i^q n L^q) = red_q(X_i^q) n F^{Fr_q} = K_q^(i) n F^{Fr_q}
(an element with zero top coordinates has no regular coordinate, so lies in K_q^(i), with zero Sp-part).  So
Phi_{i+1}^q = F_q^(i+1) := {elements of K_q^(i) with zero Sp-part}.  The algorithm adds to K exactly the
generators of K^(0), then the images of the falls until every fall row has been shifted (the argument that a row
processed once at its state then is enough is the engine's closure invariant: a fall row is later changed only by
multiples of fall rows with later pivots, which are processed in their turn).  At that point
K_q n F^{Fr_q} = F_q and K_q contains the images of F: W_{i+1} = W_i.  Lemma 1 gives statements 1 and 3
(1 is the last column, of weight 0, of the + part), Lemma 5 gives 2 (red_q is injective modulo V_D^q), and 4
holds because a fall row is supported on free columns after its pivot: adding the echelonised falls to the
reduced basis of V_D^q changes no existing leading term.  QED

**Corrections to the outline that started this work.**
(a) No image has to be reduced against V_D individually: red_q of the image of a basis row e_c + B[c] is
red_q(C e_c) + B[c] @ R, where R holds red_q(C e_f) for the free columns f, and red_q of a unit vector is a
unit vector or minus one basis row (a gather).  The whole cost is the products B[rows] @ R (section 4).
(b) Only rows of pivot weight D are shifted (Lemma 2); the low rows are covered by the closedness of V_D.
(c) Split: a self-paired weight-(D+1) monomial has a column only in the + part (coefficient 2, (D+1)/2
generators); in the - part its generators have no regular coordinate and go directly into K_-.  The count is
dim V_{D+1} = dim V_D + sum_q (|Reg_q| + dim K_q) with Reg_q counted in orbit columns of part q.
(d) T is larger than the 80,000 quoted in the outline: at the target |Sp_q| <= 30,910 and dim T_q <= 56,046 per part
(bound from 1,405 free weight-6 columns x 22 multiples; the true value is printed at the start of the run).
(e) What was called "about 12.6 million reductions" becomes 10.56 million rows of matrix products (section 6).

## 3. Algorithm and split bookkeeping

Index bookkeeping uses the engine's `Layout(k, D+1, split=True)`: the low columns of each part are exactly the
degree-D columns in the same order (asserted at start: representatives and weights equal), offset by the number
`s0` of weight-(D+1) columns.  All tables are host numpy; only matrices are on the device.

* **Unit images** (`Lift.slots`).  For a degree-D column c of part p with orbit representative h and swap
  h' = sigma h, C = project_q . u_{1,j} . lift_p sends e_c to
  proj_q[u h] e_{col_q(u h)} + s_p proj_q[u h'] e_{col_q(u h')} (second term absent for a self-paired column),
  with the engine's tables `lift_idx`, `proj_coef` (1, s_q, 2 on the + diagonal, 0 on the - diagonal).  This is
  the same map as the engine's `shift_map`, written per source column.
* **red_q of a unit vector** (`Lift.eterm`): a top column is regular (kept apart as (row, regular id,
  coefficient)) or special (unit vector of T); a low free column is a unit vector of T; a low pivot column is
  minus its basis row (gathered from the V_D^q tiles, on the Fr columns).  The result is reduced modulo K_q.
* **Followers.**  R_{p,j,q} = (red_q C e_f for f free in part p), a (#Fr_p x dim T_q) matrix, has no regular
  coordinate (asserted).  It is stored float32 on the physical columns of K_q and kept reduced modulo K_q:
  K_q is in reduced echelon form, so after K_q gains rows with pivots P the reduced follower is
  M - M[:, P] @ rows(P) (`Follower.sync`, lazy, forced before K_q compacts its columns).  The image of a basis row
  e_c + B[c] modulo V_D^q and K_q is then eterm(c) + B[c] @ follower.  Followers live in a
  least-recently-used cache under `--cache-gb`; an evicted one is rebuilt from gathers.
* **Units.**  For a regular weight-D orbit o (rows r+ in part 0 and r- in part 1; only r+ if self-paired) and a
  shift j: Y+ = B_0[r+] @ follower(0, j), Y- = B_1[r-] @ follower(1, j),
  g1 = (Y+ + Y-)/2 + eterm(first term of the part-0 column), g2 = (Y+ - Y-)/2 + eterm(second term): by Lemma 4
  each has a single explicit top term (first term: image of the representative m; second: image of sigma m,
  whose part-0 coefficient equals the coefficient of G(sigma m)).  Self-paired o: g1 = Y+ + eterm(first term).
* **J-groups.**  J(o, j) = bits(a1 | a2) | {j} (a set of bit positions).  All D + 1 generators of a regular
  column M have J = bits(A1 | A2).  Units are sorted by J and a batch is a union of whole J-groups, so the
  generators of a regular column never straddle two batches.  In a batch the regular rows are scaled to top
  coefficient 1, sorted by column, and consecutive rows of the same column are subtracted (Lemma 5).  **Checked
  for every regular column: it is met in exactly one batch with exactly D + 1 generators (popcount(A1) for a
  self-paired column); after the pass every regular column must have been met** (assertions, disabled only by
  the negative-control hook).
* **Direct rows**: basis rows of pivot weight D whose orbit is free-type (and the low rows with
  `--include-low`): image = eterm(both terms) + B[c] @ follower; no regular term (asserted).
* **K_q** is an engine `Part` on the T columns ordered Sp first, then Fr in the degree-D order; candidates
  (differences and direct vectors, already reduced modulo K_q) are added with `Part.add_residuals`.  The rows
  of K_q with pivot in the Fr block are the falls, already in reduced echelon form for the degree-D order.
* **Fall phase.**  Round r takes every fall row of K_+ and K_- not yet shifted into target q (both targets),
  reads it as a vector on Fr (`fall_block`), marks it, and adds fall @ follower(q0, j) for all j to K_q.
  Rounds repeat until a round starts with nothing pending: `fixed_point: true`.
* **Finish.**  The fall rows are added to the V_D parts (`Part.add_residuals`); `Part.sha()` and
  `Part.profile()` of the result are the canonical hash and low profile of V_{D+1} n M_<=D;
  `contains_one` is "the weight-0 column is a pivot of the + part".  dim V_{D+1} = sum_q rank V_D^q + |Reg_q|
  + rank K_q; profile = (low profile, dim V_{D+1}).
* **`--stop-on-one`** (off by default; `--to-fixed-point` cancels it) stops when the T-column of 1 becomes a
  pivot of K_+.  The record then has `fixed_point: false`, `bounds_only: true`, `contains_one: true` (exact:
  the space found is inside V_{D+1}) and dimension, profile, hashes of the subspace found so far only.

Exactness: same rules as the engine (float32 storage of entries in [0, p), float64 products with inner dimension
at most `max_inner(p)` between reductions, `modp_`).  The additional elementwise products are
(Y+ +- Y-) x 2^-1 (|value| < 2 p^2 <= 2^49) and row scalings by a residue (< p^2).

## 4. Why this layout (table of reductors or pairs of divisor variables?)

The outline offered a table of reductor vectors (590,000 regular orbits x 25,000 floats per part) or passes
organised by pairs of divisor variables.  Neither is used:

* with followers, a generator is a row of a matrix product and costs the same wherever it is computed, so the
  only question is which generators must be in memory together;
* grouping by J puts all generators of a regular column in one batch with **each product row computed once**
  (a pair-of-variables pass computes each row about twice, 2D/(D+1); a table computes each once but must hold
  one reduced vector per regular column and be updated whenever K changes);
* memory: a table has #Reg_q x (columns of K_q) entries, the complete set of followers for one target has
  2k x #Fr x (columns of K_q): 592,020 against 2 x 14 x 25,136 = 703,808 rows at the target, the same size.
  Followers need no update when K changes other than the lazy `sync`, and any of them can be evicted.
* Because everything is reduced modulo the current K_q, the work per row is #Fr_p x (physical columns of K_q):
  it falls as K_q grows and is zero once K_q = T_q (then only the coverage check of the regular columns runs).

## 5. Validation (local; `ulimit -v 1500000`, one BLAS thread; oracle = direct closure of the engine)

For every case the harness computes V_{D+1} directly with `csd/ic_bounds_big.py` (closed, full sweep clean), then
V_D with the same engine and the lift from it, and requires the same dim V_{D+1}, the same profile
(dim of V_{D+1} n M_<=j, j <= D+1), the same `contains_one`, the same dimension per part, and the same canonical
hash of V_{D+1} n M_<=D per part.  Modes: `plain` (defaults, follower cache 0.25 GB), `low` (`--include-low`,
falls shifted only at the end), `tiny` (K tiles of 7 rows, batches of 40 units, 5 direct rows, column block 33,
follower cache 4 MB, falls shifted every 2 batches, V_D tiles of 11 rows).  Devices: `cpu` and `fake-gpu` (the
engine's stand-in that raises on any mixing of host and device arrays).

**The hash from the direct closure (exact, `lift-claim-lowrows.json`).**  Claim: the basis rows of the direct
closure with pivot weight <= D, restricted to the degree-D columns, are the reduced echelon basis of
V_{D+1} n M_<=D.  Proof: an element of V supported on the low columns is the combination of the rows of its
pivots, all low; the low rows are reduced among themselves and vanish on the top columns.  Check: for 3 cases
(6 parts) the intersection was recomputed independently in dense integer arithmetic from a scrambled basis
(left kernel of the top block, then a plain echelon form): dimension and hash equal in all 6.

| case | solutions | dim V_D | direct: dim V_{D+1}, 1 in | lift: dim, 1 in | profile equal | hashes equal (direct = lift, + / - part) | dim T, dim K (per part) | rounds | runs (device:mode) all equal | lift s / direct s (cpu, plain) |
|---|---|---|---|---|---|---|---|---|---|---|
| smoke: unit-12289-k3-D2to3 (p = 12289) | 0 | 1 | 7, no | 7, no | yes (0, 0, 1, 7) | yes 7a002c9e/7a002c9e, 419ce84f/419ce84f | [22, 19], [3, 3] | 0 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.0 / 0.0 |
| smoke: unit-12289-k4-D2to3 (p = 12289) | 0 | 1 | 9, no | 9, no | yes (0, 0, 1, 9) | yes 997762d6/997762d6, f23d672b/f23d672b | [48, 44], [4, 4] | 0 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.0 / 0.0 |
| smoke: unit-12289-k3-D3to4 (p = 12289) | 0 | 7 | 57, yes | 57, yes | yes (1, 7, 22, 42, 57) | yes 632d419b/632d419b, 6f4a8242/6f4a8242 | [27, 22], [27, 22] | 3 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.0 / 0.0 |
| smoke: unit-12289-k4-D3to4 (p = 12289) | 0 | 9 | 163, yes | 163, yes | yes (1, 9, 37, 93, 163) | yes 335830dc/335830dc, 18ce1749/18ce1749 | [81, 72], [81, 72] | 5 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.0 / 0.0 |
| smoke: unit-12289-k5-D3to4 (p = 12289) | 0 | 11 | 386, yes | 386, yes | yes (1, 11, 56, 176, 386) | yes 51994086/51994086, 9a7f9304/9a7f9304 | [194, 180], [194, 180] | 6 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.1 / 0.1 |
| smoke: unit-12289-k6-D3to4 (p = 12289) | 0 | 13 | 605, no | 605, no | yes (0, 0, 2, 111, 605) | yes e9293547/e9293547, dea12bce/dea12bce | [400, 380], [304, 287] | 6 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.2 / 0.1 |
| smoke: unit-12289-k6-D4to5 (p = 12289) | 0 | 605 | 1586, yes | 1586, yes | yes (1, 13, 79, 299, 794, 1586) | yes c13c0eb4/c13c0eb4, 3583762b/3583762b | [100, 97], [100, 97] | 2 | cpu:plain, cpu:tiny, cpu:low, fake-gpu:plain, fake-gpu:tiny, fake-gpu:low: YES | 0.1 / 0.4 |
| enter: unit-65537-k3-D3to4 (p = 65537) | 0 | 7 | 57, yes | 57, yes | yes (1, 7, 22, 42, 57) | yes 632d419b/632d419b, 6f4a8242/6f4a8242 | [27, 22], [27, 22] | 3 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.0 / 0.0 |
| enter: unit-65537-k4-D3to4 (p = 65537) | 0 | 9 | 163, yes | 163, yes | yes (1, 9, 37, 93, 163) | yes 335830dc/335830dc, 18ce1749/18ce1749 | [81, 72], [81, 72] | 5 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.0 / 0.0 |
| enter: unit-65537-k5-D3to4 (p = 65537) | 0 | 11 | 386, yes | 386, yes | yes (1, 11, 56, 176, 386) | yes 51994086/51994086, 9a7f9304/9a7f9304 | [194, 180], [194, 180] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.1 / 0.1 |
| enter: unit-65537-k6-D4to5 (p = 65537) | 0 | 605 | 1586, yes | 1586, yes | yes (1, 13, 79, 299, 794, 1586) | yes c13c0eb4/c13c0eb4, 3583762b/3583762b | [100, 97], [100, 97] | 2 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.1 / 0.3 |
| enter: unit-65537-k7-D4to5 (p = 65537) | 0 | 829 | 3473, yes | 3473, yes | yes (1, 15, 106, 470, 1471, 3473) | yes 46eee628/46eee628, 32a2dade/32a2dade | [787, 779], [787, 779] | 3 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 1.4 / 1.2 |
| enter: unit-65537-k8-D4to5 (p = 65537) | 0 | 1061 | 6885, yes | 6885, yes | yes (1, 17, 137, 697, 2517, 6885) | yes f5c360b6/f5c360b6, 5a1d2e32/5a1d2e32 | [2233, 2219], [2233, 2219] | 4 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 9.1 / 5.3 |
| enter: unit-65537-k9-D4to5 (p = 65537) | 0 | 1301 | 12616, yes | 12616, yes | yes (1, 19, 172, 988, 4048, 12616) | yes 6baf9ebf/6baf9ebf, 4d99d146/4d99d146 | [4835, 4814], [4835, 4814] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 102.4 / 22.3 |
| stay: unit-65537-k2-D2to3 (p = 65537) | 0 | 1 | 5, no | 5, no | yes (0, 0, 1, 5) | yes a3d73cb3/a3d73cb3, a1e03200/a1e03200 | [8, 6], [2, 2] | 0 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.0 / 0.0 |
| stay: unit-65537-k3-D2to3 (p = 65537) | 0 | 1 | 7, no | 7, no | yes (0, 0, 1, 7) | yes 3f468fa4/3f468fa4, 419ce84f/419ce84f | [22, 19], [3, 3] | 0 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.0 / 0.0 |
| stay: unit-65537-k4-D2to3 (p = 65537) | 0 | 1 | 9, no | 9, no | yes (0, 0, 1, 9) | yes a1f129f6/a1f129f6, f23d672b/f23d672b | [48, 44], [4, 4] | 0 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.0 / 0.0 |
| stay: unit-65537-k6-D3to4 (p = 65537) | 0 | 13 | 605, no | 605, no | yes (0, 0, 2, 111, 605) | yes 28d687bf/28d687bf, 4a9d70bc/4a9d70bc | [400, 380], [304, 287] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.2 / 0.1 |
| stay: unit-65537-k7-D3to4 (p = 65537) | 0 | 15 | 829, no | 829, no | yes (0, 0, 2, 115, 829) | yes 9d596790/9d596790, 8186e3b7/8186e3b7 | [741, 714], [416, 397] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.5 / 0.2 |
| stay: unit-65537-k8-D3to4 (p = 65537) | 0 | 17 | 1061, no | 1061, no | yes (0, 0, 2, 119, 1061) | yes b2d75689/b2d75689, b812a031/b812a031 | [1267, 1232], [532, 511] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 1.6 / 0.4 |
| proper: proper-65537-X2-k6-D3to4 (p = 65537) | 2 | 13 | 605, no | 605, no | yes (0, 0, 2, 111, 605) | yes d08c0872/d08c0872, 1c06cf3d/1c06cf3d | [400, 380], [304, 287] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.2 / 0.1 |
| proper: proper-65537-X2-k6-D4to5 (p = 65537) | 2 | 605 | 1584, no | 1584, no | yes (0, 11, 77, 297, 792, 1584) | yes a67c5760/a67c5760, 5b5b4a7a/5b5b4a7a | [100, 97], [99, 96] | 2 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.1 / 0.4 |
| proper: proper-65537-X2-k7-D3to4 (p = 65537) | 2 | 15 | 829, no | 829, no | yes (0, 0, 2, 115, 829) | yes 52a5f310/52a5f310, a3a9deae/a3a9deae | [741, 714], [416, 397] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.5 / 0.2 |
| proper: proper-65537-X2-k7-D4to5 (p = 65537) | 2 | 829 | 3471, no | 3471, no | yes (0, 13, 104, 468, 1469, 3471) | yes 4aa54063/4aa54063, 5e4462fc/5e4462fc | [787, 779], [786, 778] | 3 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 1.3 / 1.2 |
| proper: proper-65537-X2-k8-D3to4 (p = 65537) | 4 | 17 | 1061, no | 1061, no | yes (0, 0, 2, 119, 1061) | yes 7fd9018e/7fd9018e, 736a8393/736a8393 | [1267, 1232], [532, 511] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 1.5 / 0.4 |
| proper: proper-65537-X2-k8-D4to5 (p = 65537) | 4 | 1061 | 6881, no | 6881, no | yes (0, 13, 133, 693, 2513, 6881) | yes 381f9fd9/381f9fd9, 59e48349/59e48349 | [2233, 2219], [2231, 2217] | 4 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 8.8 / 4.8 |
| proper: proper-65537-X1424-k9-D4to5 (p = 65537) | 2 | 1301 | 12614, no | 12614, no | yes (0, 17, 170, 986, 4046, 12614) | yes 6a37595d/6a37595d, 58cf8074/58cf8074 | [4835, 4814], [4834, 4813] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 102.1 / 22.9 |
| proper: proper-15361-k6-D4to5 (p = 15361) | 2 | 605 | 1584, no | 1584, no | yes (0, 11, 77, 297, 792, 1584) | yes 559de8b3/559de8b3, 69490952/69490952 | [100, 97], [99, 96] | 2 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.1 / 0.4 |
| proper: proper-15361-k7-D4to5 (p = 15361) | 4 | 829 | 3469, no | 3469, no | yes (0, 11, 102, 466, 1467, 3469) | yes 20d00e19/20d00e19, 90ee6687/90ee6687 | [787, 779], [785, 777] | 3 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 1.5 / 1.3 |
| proper: proper-15361-k7-D3to4 (p = 15361) | 4 | 15 | 829, no | 829, no | yes (0, 0, 2, 115, 829) | yes 92f7a73b/92f7a73b, bc0422fe/bc0422fe | [741, 714], [416, 397] | 6 | cpu:plain, cpu:low, fake-gpu:plain, fake-gpu:low: YES | 0.6 / 0.2 |
| tiny: unit-65537-k4-D3to4 (p = 65537) | 0 | 9 | 163, yes | 163, yes | yes (1, 9, 37, 93, 163) | yes 335830dc/335830dc, 18ce1749/18ce1749 | [81, 72], [81, 72] | 4 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: unit-65537-k5-D3to4 (p = 65537) | 0 | 11 | 386, yes | 386, yes | yes (1, 11, 56, 176, 386) | yes 51994086/51994086, 9a7f9304/9a7f9304 | [194, 180], [194, 180] | 5 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: unit-65537-k6-D3to4 (p = 65537) | 0 | 13 | 605, no | 605, no | yes (0, 0, 2, 111, 605) | yes 28d687bf/28d687bf, 4a9d70bc/4a9d70bc | [400, 380], [304, 287] | 6 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: unit-65537-k6-D4to5 (p = 65537) | 0 | 605 | 1586, yes | 1586, yes | yes (1, 13, 79, 299, 794, 1586) | yes c13c0eb4/c13c0eb4, 3583762b/3583762b | [100, 97], [100, 97] | 2 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: unit-65537-k7-D3to4 (p = 65537) | 0 | 15 | 829, no | 829, no | yes (0, 0, 2, 115, 829) | yes 9d596790/9d596790, 8186e3b7/8186e3b7 | [741, 714], [416, 397] | 6 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: unit-65537-k7-D4to5 (p = 65537) | 0 | 829 | 3473, yes | 3473, yes | yes (1, 15, 106, 470, 1471, 3473) | yes 46eee628/46eee628, 32a2dade/32a2dade | [787, 779], [787, 779] | 4 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: proper-65537-X2-k6-D4to5 (p = 65537) | 2 | 605 | 1584, no | 1584, no | yes (0, 11, 77, 297, 792, 1584) | yes a67c5760/a67c5760, 5b5b4a7a/5b5b4a7a | [100, 97], [99, 96] | 2 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: proper-65537-X2-k7-D4to5 (p = 65537) | 2 | 829 | 3471, no | 3471, no | yes (0, 13, 104, 468, 1469, 3471) | yes 4aa54063/4aa54063, 5e4462fc/5e4462fc | [787, 779], [786, 778] | 4 | cpu:tiny, fake-gpu:tiny: YES | - |
| tiny: proper-15361-k7-D4to5 (p = 15361) | 4 | 829 | 3469, no | 3469, no | yes (0, 11, 102, 466, 1467, 3469) | yes 20d00e19/20d00e19, 90ee6687/90ee6687 | [787, 779], [785, 777] | 4 | cpu:tiny, fake-gpu:tiny: YES | - |
| k10d4: unit-65537-k10-D4to5 (p = 65537) | 0 | 1549 | 18569, no | 18569, no | yes (0, 0, 2, 138, 3163, 18569) | yes 50831915/50831915, c443b473/c443b473 | [9110, 9081], [7543, 7517] | 6 | cpu:plain, fake-gpu:plain: YES | 480.8 / 95.1 |
| k10d4-proper: proper-65537-X1424-k10-D4to5 (p = 65537) | 14 | 1549 | 18569, no | 18569, no | yes (0, 0, 2, 138, 3163, 18569) | yes ab41edaa/ab41edaa, 94905423/94905423 | [9110, 9081], [7543, 7517] | 6 | cpu:plain: YES | 458.6 / 94.2 |
| bigp: p16777213-a1-X5-k5-D3to4 (p = 16777213) | None | 11 | 386, yes | 386, yes | yes (1, 11, 56, 176, 386) | yes 51994086/51994086, 9a7f9304/9a7f9304 | [194, 180], [194, 180] | 6 | cpu:plain, cpu:tiny, fake-gpu:plain, fake-gpu:tiny: YES | 0.1 / 0.1 |
| bigp: p16777213-a1-X5-k6-D3to4 (p = 16777213) | None | 13 | 605, no | 605, no | yes (0, 0, 2, 111, 605) | yes 1b30156e/1b30156e, 36a26de1/36a26de1 | [400, 380], [304, 287] | 6 | cpu:plain, cpu:tiny, fake-gpu:plain, fake-gpu:tiny: YES | 0.5 / 0.2 |
| bigp: p16777213-a1-X5-k6-D4to5 (p = 16777213) | None | 605 | 1586, yes | 1586, yes | yes (1, 13, 79, 299, 794, 1586) | yes c13c0eb4/c13c0eb4, 3583762b/3583762b | [100, 97], [100, 97] | 2 | cpu:plain, cpu:tiny, fake-gpu:plain, fake-gpu:tiny: YES | 0.2 / 0.7 |
| k10d5: unit-65537-k10-D5to6 (p = 65537) | 0 | 18569 | 60460, yes | 60460, yes | yes (1, 21, 211, 1351, 6196, 21700, 60460) | yes a2f9245b/a2f9245b, 11874086/11874086 | [1965, 1954], [1965, 1954] | 2 | cpu:plain, fake-gpu:plain: YES | 23.6 / 0.0 |
| k11d4: unit-65537-k11-D4to5 (p = 65537) | 0 | 1805 | 25033, no | 25033, no | yes (0, 0, 2, 142, 3441, 25033) | yes 7e993959/7e993959, f9befbc4/f9befbc4 | [15711, 15673], [10502, 10472] | 7 | cpu:plain: YES | 1827.5 / 257.1 |

The row `k10d5` (k = 10, D = 5 -> 6) is compared with the on-card result of `RESULTS.md` (V_6(10) is the full
space, 60,460 columns, closed, 1 in V_6), not with a local direct closure, which does not fit the light limits:
the expected hash is that of the full space (`full_space_sha`, all columns pivots).  V_5(10) (rank 18,569 in
21,700 columns) was computed locally by the engine.  The rows `k10d4*` and `k11d4` ran under the heavy lock (peak 1.3 GB and 2.2 GB resident).
Columns: "runs" lists every (device:mode) run of the case, all of which must agree with the direct closure;
"rounds" and the times are those of the first run.

Instances: `unit-65537` = y^2 = x^3 + x over F_65537, X = 4271 (the unit instance of `instances/`, no
solutions); `unit-12289` = same curve over F_12289, X = 9362; `proper-65537-X2` and `-X1424` = y^2 = x^3 + 4021 x
+ 9917 over F_65537 with target abscissa 2 resp. 1424 (the curve of `proper-p65537-s12`; the number of solutions
at that depth is in the table); `proper-15361` = the same curve over F_15361 (X = 4491 at k = 6, 7776 at k = 7).
The generator is `chain_generator(p, a, b, X, k)` of the engine.

**Negative controls (exact, `lift-negative.json`; 8 cases x 3 controls = 24 runs, test hooks, exit code 4,
records marked `fixed_point: false`).**  (i) `--test-no-fall-shifts` (the fall rounds after the main phase are
skipped; a round interleaved in the main phase, one per 8 unit batches, still runs, so some falls are shifted once): the answer changed in 8 of 8 cases; `contains_one` flipped from true to false in the three cases
where 1 enters at D + 1 (unit k = 6, 7, 8, D = 4 -> 5: dimension 1523 / 2844 / 4617 instead of 1586 / 3473 /
6885), and the dimension fell in the others (e.g. unit k = 8, D = 3 -> 4: 137 instead of 1061).  So the
iteration of the falls is necessary and the comparison detects its absence.  (ii) `--test-drop-shift j` (all
rows of one shift index dropped, j = 0 and j = k - 1): the dimension changed in 2 of 16 runs (unit k = 7 and 8,
D = 3 -> 4, j = k - 1: 828 instead of 829, 1060 instead of 1061); in the other 14 the dropped rows were
redundant and the result was unchanged, as it must be then.  With the hook off, a dropped generator of a regular
column is caught by the coverage assertion.

**Engine checkpoint, device budget, resume (exact, `lift-ckpt-k8-D4.json`; unit instance, k = 8, D = 4 -> 5).**
The V_4 checkpoint was written by `ic_bounds_big_run.py closure --ckpt` and read by `ic_bounds_lift_run.py lift
--ckpt`: (a) cpu, (b) fake-gpu with the basis on the device, (c) fake-gpu with `--basis-on host --dev-gb 0.08
--work-gb 0.05` (the stand-in raises above 80 MB of live device arrays: peak 75 MB, 243 follower evictions),
(d) two SIGTERMs (exit 75, 75) and `--resume` with a lift checkpoint after every item: all four equal the direct
closure (dimension, profile, contains_one, both hashes).

**`selftest --device fake-gpu` (exact, `lift-selftest-fake-gpu.json`)**: 7 built-in cases (3 with 1 entering,
4 without; p = 12289, 15361, 65537), each on cpu, on the device and on the device with tiny settings, against
embedded values taken from the direct closure (`lift-selftest-expected.json`): passed.

**Added 2026-11-11 (`lift-k10d5-proper.json`, set `k10d5-proper`, cpu:plain, no memory limit, 4 BLAS threads).**
proper-65537-X1424, k = 10, D = 5 -> 6, 14 solutions: direct dim V_6 = 60,446, 1 not in; lift 60,446, 1 not in; profile
(0, 7, 197, 1337, 6182, 21686, 60446) equal; hashes equal (12711885/12711885, dd697484/dd697484); dim T [1965, 1954],
dim K [1958, 1947]; 2 rounds; lift 36.7 s against 760.6 s direct. Same lift and engine hashes as below; the harness
hash changed because this set was added. This case is not counted in the totals that follow.

Totals: 170 lift runs over 37 distinct (instance, k, D) cases compared with the direct closure (smoke 42,
enter 28, stay 24, proper 40, tiny 18, k10d4 2, k10d5 2, k10d4-proper 1, k11d4 1, bigp 12), 0 mismatches.  Every
record carries the sha256 of `ic_bounds_lift.py` (93e8d19b6863... in all of them), of the engine and of the harness
(two harness hashes: the set `bigp` was added to the harness after the other sets had run; nothing else changed).

## 6. Cost model and measured counts

Counts per run are in the record: `counts.images` (rows of B @ follower), `counts.eterm_rows` (unit-vector
images reduced), `reductions` (their sum), `madds` (`gemm`, `k_reduce_update_rref`, `follower_sync`).

**Measured (exact counts, one CPU thread, `lift-k10d5.json`, `lift-k10d4.json`).**

| case | regime | rows of B @ follower | unit images | multiply-adds: products / follower sync / K echelon / total | lift seconds | for comparison |
|---|---|---|---|---|---|---|
| k = 10, D = 5 -> 6 (largest local case in the regime of the target: 37,972 regular against 788 special top columns, dim T = 1965 + 1954) | K = T after 14 % of the schedule; the remaining batches only run the coverage check | 45,876 (of 317,345 scheduled) + 30,920 fall rows | 141,045 | 6.62e10 / 1.20e11 / 1.90e10 / 2.05e11 | 23.6 | closure of V_5(10): 5.40e11 multiply-adds, 97 s (same thread); direct closure of V_6(10) on the H200: 9.34e12 (RESULTS.md) |
| k = 10, D = 4 -> 5 (tiny V_4 = 1549 inside 6196 columns: 1960 regular against 13,544 special top columns, dim T = 9110 + 9081) | K ends at 7543 + 7517 | 60,720 + 32,280 | 1,924,996 | 7.47e11 / 3.61e11 / 4.91e12 / 6.01e12 | 480.8 | direct closure of V_5(10): 5.40e11, 95.1 s |

In the regime "almost everything special" the lift is the direct closure in another basis and is slower (11
times more multiply-adds at k = 10, D = 4 -> 5); it is exact there, which is all that was asked.  In the regime
of the target (few special columns) it used 2.05e11 multiply-adds where the direct closure used 9.34e12.

**Extrapolation to k = 14, D = 6 -> 7 (ESTIMATE, not a measurement).**  Sizes from section 7: #Fr = 25,128,
dim T <= 56,046 per part, 56 followers, 10.56 million product rows.

* Scenario A, K_q becomes all of T_q early (what happened at k = 10, D = 5 -> 6, and what the conjecture
  "V_7(14) is the full space" predicts).  Scaling each measured term with its own size law:
  follower sync 0.51 x 56 x #Fr x T^2 / 2 = 1.1e15; products until K is full (19.6 x dim T rows, at 0.28 T
  columns on average, as measured) 0.9e15; echelon of K 1.26 x 2 T^3 = 0.4e15; adding the 50,256 falls to the
  V_6 basis 0.14e15.  **About 3e15 multiply-adds (2e15 to 4e15)**, 0.6 of the V_6(14) closure run
  (4.69e15 in 8077 s on the H200): **about 1 to 1.5 hours at that run's overall rate** (0.58e12 per second), less
  if the large products run nearer the DGEMM rate.  The laws are fitted to one small case; the constant can be
  off by a factor of a few.
* Scenario B, K_q stops at codimension c in T_q: add 5.28e6 rows x 25,128 x c per part, i.e. 1.3e15 per part
  for c = 10,000 and 7.4e15 per part for c = 56,046 (K empty, the worst case: 1.5e16 in total, 7 hours at
  0.58e12 per second).
* Reading the 45 GB checkpoint and verifying its file hashes (twice) comes on top (minutes).

## 7. Memory at the target (k = 14, D = 6 -> 7, p = 65537), from `plan` (**estimate**; sizes are upper bounds)

`python3 ic_bounds_lift_run.py plan --k 14 --D 6 --free 25136,25120 --free-top 703,702 --work-gb 16 --dev-gb 130 --basis-on host`
(the split 703 + 702 of the 1,405 free weight-6 columns is a guess; only the sum enters):

| quantity | value |
|---|---|
| rank of V_6 per part / free columns | 224,688 + 224,234 / 25,136 + 25,120 |
| weight-7 orbit columns per part | 592,020 |
| special weight-7 columns per part (bound: 1,405 x 22) | <= 30,910 |
| dim T per part | <= 56,046 |
| regular weight-6 orbits, units (orbit, shift) per target | 187,849; 2,629,886 |
| rows of B @ follower in total (both targets) | 10,558,884 |
| one follower (25,136 x 56,046 float32) | <= 5.25 GB |
| all 56 followers / those of one J-group for one target (14) | <= 294 GB / <= 73.4 GB, shrinking as K grows |
| K per part, worst case over its rank | 2.9 GB |
| V_6 basis (float32) | 21.0 + 21.0 GB |
| batch buffers (`--work-gb 16`: 2,274 units per batch) | 16 GB |
| device, `--dev-gb 130 --basis-on host`: cache 109.6 + work 16 + K 2.9 + 1.5 | 130 GB (by construction) |
| device during the final step (one V_6 part + work) | 37 GB |
| host: V_6 basis + layouts and index tables + one tile in transit | about 48 GB |

The follower cache is filled as needed and evicts the least recently used follower when `--cache-gb` (derived
from `--dev-gb`) would be exceeded; an evicted follower is rebuilt from gathers and one reduction.  While K is
small the cache cannot hold all followers of both targets (294 GB); the schedule alternates the two targets,
so there will be evictions in the first batches (each rebuild: #Fr x rank K x columns of K multiply-adds).
With `--basis-on device` (basis 42 GB on the card) the cache would be 67.6 GB, less than one J-group needs
at the start: use `--basis-on host` on the 141 GB card.

## 8. Remote command lines

Remote layout: `/workspace/code` holds `ic_bounds_big.py`, `ic_bounds_lift.py`, `ic_bounds_big_run.py`,
`ic_bounds_lift_run.py`; the V_6 checkpoint is the directory given to the closure's `--ckpt` (the lift reads
its newest complete `ck-*` slot, verifying every file hash, and refuses a checkpoint whose manifest does not
say `closure_complete`).

```
cd /workspace/code
python3 ic_bounds_lift_run.py selftest --device gpu --out /workspace/out/lift-selftest-gpu.json      # must exit 0
python3 ic_bounds_lift_run.py plan --k 14 --D 6 --free 25136,25120 --free-top 703,702 --work-gb 16 --dev-gb 130 --basis-on host

# small real-device check before the target: k = 10, D = 5 -> 6 (V_5 computed first; expect dim 60460, contains_one true,
# dim_K [1965, 1954], low hashes a2f9245b18b2e968... / 118740867a8b80f4... as in validation/lift-k10d5.json)
python3 ic_bounds_lift_run.py lift --instance /workspace/in/unit-p65537-a1-s10-k10.json --D 5 --device gpu --work-gb 4 --out /workspace/out/lift-k10-D5to6.json

# the target: V_6(14) -> V_7(14)
python3 ic_bounds_lift_run.py lift --ckpt /workspace/ckpt/k14-D6 --instance /workspace/in/unit-p65537-a1-s10-k14.json \
  --device gpu --basis-on host --dev-gb 130 --work-gb 16 \
  --lift-ckpt /workspace/ckpt/k14-D6-lift --ckpt-minutes 20 --resume \
  --status /workspace/out/k14-D6to7.status.json --out /workspace/out/k14-D6to7.json
```

Exit codes: 0 fixed point (`fixed_point: true`: `contains_one`, `dim_VD1`, `profile`, `low_basis_sha256` are
statements about V_7); 2 stopped on 1 (`--stop-on-one` only); 75 SIGTERM (run the same command again: `--resume`
continues from the lift checkpoint).  Add `--stop-on-one` to stop as soon as 1 is found (the record then gives
lower bounds only).  Add `--include-low` to make the result independent of the closedness of the V_6 basis
(the 73,587 low rows: about 20 % more rows).  The log is JSON lines (`stage.pos` of `of` schedule items, `rank_K`, `ncols_K`, `cache_mb`,
`eta_main_s`); `--status` holds the latest line.  The checkpoint must come from the closure run whose record
says `closed: true` (the checkpoint written at the end of the closure phase precedes the sweep; the lift
records the checkpoint's file hashes as `vd.basis_id` and cannot know the sweep's verdict).

## 9. Not verified on a real device

Nothing in `ic_bounds_lift.py` has run under cupy.  Everything was exercised through `--device fake-gpu`
(numpy semantics, host/device mixing forbidden, live device bytes counted and, in the checkpoint test, limited).
The engine calls it reuses (`Part.accumulate`, `add_residuals`, `_rows`, `compact`, `modp_`, tile spilling) ran
on an H200 according to `README.md` section 12.  New cupy uses, all **unverified**:

* indexing: `X[rows_dev] = block` and `X[rows_dev]` (row gather/scatter with one index array), `Yt[ri_dev, pos_dev]
  += vals` (two index arrays, augmented assignment, unique pairs), `Yt[ri_dev, f0:] += G` (index array and slice,
  augmented), `M[a:b] = E.astype(float32)`, `M[:, idx_dev]` (compaction of a follower: a copy of up to 5 GB),
  `M[:, idx_dev] = 0.0`, `A[arange_dev, piv_dev] = 1.0`, `C[:ns] = Cd`;
* `out[:, a:b]` passed as a view and updated in place by `+=`, `-=` and `modp_` (the engine does the same);
* `G *= vec_dev[:, None]`, `g *= python_float`; `float32 -> float64` casts of column blocks of a follower;
* `bool(asnumpy(any(X != 0)))` on a 0-d result; products with an empty operand (n x 0 matrices) when K_q = T_q;
* the memory behaviour: followers are reallocated at every compaction of K (when 10 % of its columns are dead)
  and after evictions, so the cupy pool can fragment (the k = 14 closure was run with the pool disabled);
  `--dev-gb` budgets the cache, the K worst case, `--work-gb` and 1.5 GB of slack, but the peak of the batch
  buffers (2,274 units: about 5 copies of 4,548 x 56,046 float64 plus the product blocks) is a formula, not a
  measurement at that size; if the first log lines show `dev_used_mb` near the card's size, restart with a
  smaller `--work-gb` (the lift checkpoint is bound to the batch settings: use a new `--lift-ckpt` directory);
* speed: host-to-device transfer of the batch rows with `--basis-on host` (about 0.5 GB per batch from pageable
  memory), the DGEMM rate on 25,000 x 15,000 blocks, the cost of follower rebuilds while K is small.

`selftest --device gpu` (7 cases, cpu against device against embedded values, plus the engine's arithmetic test)
and the k = 10, D = 5 -> 6 run of section 8 exercise all of these calls except the spill of the basis
(`--basis-on host`), which applies only to a basis read from a checkpoint.  To exercise it on the card first:

```
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k10.json --D 5 --device gpu --work-gb 8 --ckpt /workspace/ckpt/k10-D5 --out /workspace/out/k10-D5.json
python3 ic_bounds_lift_run.py lift --ckpt /workspace/ckpt/k10-D5 --device gpu --basis-on host --dev-gb 6 --work-gb 2 --out /workspace/out/lift-k10-D5to6-host.json
```

## 9a. Runs on a real card (added 2026-11-11)

Section 9 was written before any run on a device. Since then, on H200 cards (records in `runs/`): `selftest --device gpu`
passed twice; k = 10, D = 5 -> 6 on the unit instance (60,460, 1 in); k = 14, D = 6 -> 7 on the unit instance over
p = 65537 (1,683,218 = the whole space, 1 in; 751 s); and k = 14, D = 6 -> 7 on the random curve over p = 114689 with
2,370 solutions (1,680,848 = the ideal bound exactly, 1 not in; 4,679 s). The last is the soundness control at the
target size: the lifted space lies in the ideal, so any spurious element would have pushed the dimension above the
bound. See `RESULTS.md`.

## 10. Limits

* The theorem is conditional on its input: V_D must be the closed degree-D closure (the engine's `closed: true`).
  A checkpoint does not carry that verdict; `--include-low` removes the dependence on closedness (the input must
  then only be a subspace of V_D containing g; Lemma 1 holds from any such start), at about 20 % more rows.
* `fixed_point: true` is an exact statement only as far as the arithmetic is exact on the device that ran it
  (same caveat as the engine; `guard_hits` is in the record and must be 0).
* The lift does not produce a basis of V_{D+1}, only its dimension, profile, the basis of V_{D+1} n M_<=D and
  membership of 1.  There is no independent certificate of the result (no sweep in M_<=D+1): its correctness
  rests on the proof of section 2, the assertions (every regular column met once with all its generators; no
  regular term in followers and direct rows) and the comparison with the direct closure on the 37 cases above.
* Largest validated case in the regime of the target: k = 10, D = 5 -> 6 (60,460 columns), and there the
  oracle is the on-card record, not a local closure.  Largest case validated against a local direct closure:
  k = 11, D = 4 -> 5 (35,443 columns, dim V_5 = 25,033, 1 not in V_5; lift 1828 s against 257 s direct, under the
  heavy lock, peak 2.2 GB) and k = 10, D = 4 -> 5 (21,700 columns), unit and proper instance.
* Primes: 12289, 15361, 65537 and 16777213 (`bigp`: y^2 = x^3 + x, X = 5, `max_inner` = 4).  For these the inner dimension never exceeds `max_inner(p)`, so the split
  of `addmul` is exercised only by `selftest` (against Python integers, chunks of 4 and 7) and by the p = 16777213 runs.
* The lift is slower than the direct closure when most weight-(D+1) columns are special (tiny V_D).
* Resume was tested at k = 8 only.  The fall rounds interleaved in the main phase (`--fall-every`) change
  the order of work, not the result (tested with 2, 8 and never-before-the-end).
* `plan` gives upper bounds from the free-column counts; the true number of special columns is printed in the
  `start` line of the run.  `eta_main_s` is a linear extrapolation over schedule items and overestimates once
  K has grown.
* The time and multiply-add predictions of section 6 are extrapolations from one small case.
