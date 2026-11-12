# Fall spaces of the two-summand prime chain closure (bounds program, item 5)

Structure of `V_4` and of the fall space `F_4 = V_4 ∩ M_<=3` of the two-summand prime
chain presentation `{NF(S_3)} ∪ chains` in the group algebra `A_k = F_p[(Z/2^k)^2]`, for
generic instances (random curve, random target of large order), k = 6..12.

Headline (exact for the instances and depths computed; the reduction lemma is a written
proof): **for k >= 7 the degree-4 closure is the free high-level extension of the depth-7
closure,**

    V_4(k) = W_k := (+)_{h}  tau_h ( V_4(7) ∩ M_{<= 4 - wt h} ),         h on levels >= 7, wt h <= 2,

so that `dim V_4(k) = 829 + 230 (k-7) + 2 C(2k-14, 2) = 4k^2 + 172k - 571` and
`dim F_4(k) = 115 + 4 (k-7) = 4k + 87`, and `1 ∉ W_k` for every k because `1 ∉ V_4(7)`.
`W_k` is the explicit, k-uniform closed containing space asked for: it is defined from the
single finite object `V_4(7)` (an 829 x 3473 matrix over F_p), it contains every admissible
generator row, it is closed under every admissible variable product for all k >= 8 as soon
as `W_8` is (reduction lemma, proved below), and it excludes 1. The base cannot be taken at
depth 6: `V_4(6)` carries 42 dimensions produced by wrap-around at level 5.

Evidence labels: **kernel** (Lean-checked), **proof** (written argument in full here),
**exact** (certified computation, certificate named), **empirical** (pattern in data),
**conjecture**.

## 1. Setting

Modulo the chains a squarefree chain monomial is an exponent pair `(a1, a2) ∈ (Z/2^k)^2`, its
degree is the Hamming weight `wt(a1) + wt(a2)`, and the chain variable `u_{i,j}` acts as the
translation `tau_{2^j e_i}` (addition of `2^j` to block i modulo `2^k`; a carry lowers the
weight). The normalized generator

    g = NF(S_3(x1, x2, t)) = x1^2 x2^2 - 2t (x1^2 x2 + x1 x2^2) + t^2 (x1^2 + x2^2) - 2(t^2 + a) x1 x2
        - (2at + 4b)(x1 + x2) + (a^2 - 4bt)

is supported on `{0,1,2}^2` (levels 0 and 1 only) and has weight 2. `V_D` is the least
subspace of `M_<=D` containing every `tau_mu g` with `wt(mu) <= D-2` (the generator rows,
`2k^2 + k + 1` of them at D = 4) and closed under `tau_{2^j e_i}` applied to its elements of
weight < D (`formal/ICFallDegrees/Closure.lean`; working modulo the chains equals the ordinary
closure for D >= 2). `Gen_k` denotes the span of the generator rows at depth k. Levels are the
bit positions j of the exponents; the *level set* of a monomial is the set of levels set in
either block.

Engine: `csd/ic_bounds_fall.py` (once-per-pivot shift closure as in
`csd/ic_chain_shift_closure.py`, with the birth round of every pivot, verified closed by
re-reducing all 2k shifts of all rows of weight < D; reduced row echelon form with columns
sorted by weight, so that the rows whose pivot has weight <= d span exactly `V_D ∩ M_<=d`).
The free extension `W_k` is built by `free_extension`, and `check_containing_space` verifies,
without using `V_4(k)`, that a given spanning set contains all generator rows, is closed under
every admissible product and does not contain 1 (the exclusion certificate of
`ic_fall_certificates` in the group algebra).

## 2. Instances

| tag | p | curve `y^2 = x^3 + a x + b` | target (x, y) | target order | solutions at k = 6..10 |
|---|---|---|---|---|---|
| `p40961-s7` (primary, exported) | 40961 = 2^13·5 + 1 | a = 21223, b = 9887 | (4748, 524) | 40817 (= curve order, prime) | 2, 2, 2, 6, 26 |
| `p1073153-s1` (unit ideal, multipliers) | 1073153 = 2^13·131 + 1 | a = 281783, b = 132345 | (247294, 569236) | 536971 | 0 at k = 6, 7 |
| `12289:3` (cross-check) | 12289 | random (seed 3) | random | see `analysis.json` | |
| `12289:j1728` (cross-check) | 12289 | a = 1, b = 0 | random target, seed 10 | | depth-10 probe family |

Instances are drawn by `csd.ic_bounds_fall.random_instance(p, 6, seed)`: random nonzero
a, b with nonsingular curve, random rational target with y ≠ 0 and point order > 2^8. The
generator is depth-independent as exponent pairs; the chain-mask encoding and the solution set
are recomputed at every depth (`encode_generator`, `prime_points` with the depth-k subgroup
generator). Hashes: `p40961-s7/instance.json` (`instance_hash`), per-basis sha256 in the
`k*-D4.json` files.

## 3. Dimensions (exact)

Primary instance, degree 4 (`k*-D4.json`; each closure verified closed, 9 rounds of
shift-and-fall at every depth):

| k | \|M_<=4\| | dim V_4 | 4k^2+172k-571 | dim V_4 ∩ M_<=3 | 4k+87 | dim V_4 ∩ M_<=2 | 1 ∈ V_4 |
|---|---|---|---|---|---|---|---|
| 6 | 794 | 605 | 605 | 111 | 111 | 2 | no |
| 7 | 1471 | 829 | 829 | 115 | 115 | 2 | no |
| 8 | 2517 | 1061 | 1061 | 119 | 119 | 2 | no |
| 9 | 4048 | 1301 | 1301 | 123 | 123 | 2 | no |
| 10 | 6196 | 1549 | 1549 | 127 | 127 | 2 | no |
| 11 | 9109 | 1805 | 1805 | 131 | 131 | 2 | no |
| 12 | 12951 | 2069 | 2069 | 135 | 135 | 2 | no |

Degree 5 (`k*-D5.json`; `V_5 ∩ M_<=4` exported as `k*-D5-low.npz`):

| k | \|M_<=5\| - \|Z\| | dim V_5 | profile dim V_5 ∩ M_<=d, d = 0..5 | certificate |
|---|---|---|---|---|
| 6 | 1584 | 1584 | 0, 11, 77, 297, 792, 1584 | V_5 = I ∩ M_<=5 (rank = target, M_<5 onto Z): LFD <= 5, V_5 ∩ M_<=4 = I ∩ M_<=4 |
| 7 | 3471 | 3471 | 0, 13, 104, 468, 1469, 3471 | same |
| 8 | 6883 | 6883 | 0, 15, 135, 695, 2515, 6883 | same (13 rounds, 1274 s) |

At k <= 9 the degree-5 closure is the whole degree-5 part of the ideal, so `V_5 ∩ M_<=4` is
`I ∩ M_<=4` (every weight-<=4 polynomial vanishing on Z) and carries no closure structure of its
own; the first depth where `V_5` separates from the ideal is k = 10 (`experiments/ic-chain-depth10`).

## 4. The k-uniform description (main result)

### 4.1 Statement

Fix an instance and let `S := V_4(7) ⊂ A_7`, read as polynomials with exponents in `[0, 128)^2`
inside `A_k` for k >= 7. Let `High_k` be the exponent pairs supported on levels >= 7 (both
blocks) of weight <= 4 and set

    W_k := sum_{h ∈ High_k} tau_h ( S ∩ M_{<= 4 - wt h} )  ⊂ M_<=4 ⊂ A_k.

The sum is direct (distinct h give disjoint monomial supports), so

    dim W_k = s_4 + 2(k-7) s_3 + C(2k-14, 2) s_2 + C(2k-14, 3) s_1 + C(2k-14, 4) s_0,
    s_d := dim (S ∩ M_<=d);   here (s_0, .., s_4) = (0, 0, 2, 115, 829),

    dim W_k = 829 + 230 (k-7) + 2 C(2k-14, 2) = 4k^2 + 172k - 571.

**Theorem (reduction; proof).** For every k >= 7:
(i) `W_k ⊇ Gen_k`; (ii) `1 ∈ W_k ⟺ 1 ∈ S`; (iii) if `W_8` is closed under admissible
variable products in `A_8`, then `W_k` is closed in `A_k` for every k >= 8.
Hence, under the hypothesis of (iii), `V_4(k) ⊆ W_k` for all k >= 8, so `1 ∉ V_4(k)` and
`dim V_4(k) <= 4k^2 + 172k - 571` for every k >= 8.

**Verified (exact)** on the primary instance: `V_4(k) = W_k` as row spaces for k = 7..12,
and `check_containing_space(W_k)` passes (contains the generator rows, closed, 1 ∉ W_k) for
k = 7..12 independently of `V_4(k)`; in particular the hypothesis of (iii) holds. Also
`V_4(k) ∩ M_{levels <= 6} = S` for k = 8..12 (the low part is stable), and the same holds
on the cross-check instances (Section 8). The hypothesis of (iii) is a finite computation at
depth 8 for the instance (1061 rows, 2517 columns, 16·119 shift reductions); the theorem
removes k from the problem.

### 4.2 Proof of the reduction theorem

Write every exponent pair as `low + high` with `low ∈ [0,128)^2` and `high ≡ 0 mod 128`
blockwise; `wt(low + high) = wt(low) + wt(high)`. An element of `A_k` decomposes uniquely by
its high pattern, and `W_k` is the set of elements whose h-component lies in
`tau_h(S ∩ M_{<= 4 - wt h})` for every h. Since components of distinct h have disjoint
supports, an element of `W_k` has weight <= 3 iff each of its components does.

(i) For `wt(mu) <= 2` write `mu = mu_low + mu_high`. Then `tau_mu g = tau_{mu_high}(tau_{mu_low} g)`;
`tau_{mu_low} g` is a generator row of `A_7` (no wrap: `mu_low + s <= 96 + 2 < 128` for
`wt(mu_low) <= 2`, `s ∈ {0,1,2}^2`), so it lies in `S` with weight <= 2 + wt(mu_low) =
4 - wt(mu_high). Hence `tau_mu g ∈ W_k`.

(ii) The h = 0 component of 1 is 1, and `W_k ∩ (h = 0 component) = S`.

(iii) By linearity it suffices to show `tau_u(tau_h v) ∈ W_k` for `h ∈ High_k`,
`v ∈ S ∩ M_{<= 3 - wt h}` and a variable `u = 2^j e_i`.
*Case j >= 7.* `tau_u tau_h v = tau_{h+u} v` and `h + u` is again supported on levels >= 7
(sum of two such), with `wt(h+u) <= wt(h) + 1`: adding one bit either sets a clear bit, or
carries through a run of set bits (clearing them, setting the next clear one), or wraps out
at the top level (weight decreases). So `wt(h+u) + wt(v) <= 4` and `tau_{h+u} v ∈ W_k`.
*Case j <= 6.* `tau_u tau_h v = tau_h(tau_u v)`. The exponents of v are < 128, so adding
`2^j` (j <= 6) gives values < 256: `tau_u v = v_0 + tau_{128 e_i} v_1`, where `v_0` collects the
monomials with no carry out of level 6 and `v_1` the carried ones with the bit 128 removed;
both have exponents < 128, and the decomposition is the same in `A_8` and in `A_k` (the carry
into level 7 is exact in both). Weights: `wt(v_0) <= wt(v) + 1`; a carried monomial m of v
had bits j..6 of block i set and after the addition has them cleared and bit 7 set, so
`wt(v_1) <= wt(m) - (7-j) <= wt(v) - 1`. Now `tau_u v = tau_u(tau_0 v) ∈ W_8` by the
hypothesis (`tau_0 v ∈ W_8` has weight <= 3), and by uniqueness of the high decomposition in
`A_8`, `v_0 ∈ S` and `v_1 ∈ S ∩ M_<=3`. Therefore `tau_h v_0 ∈ tau_h(S ∩ M_{<= 4 - wt h})`
because `wt(v_0) <= wt(v) + 1 <= 4 - wt(h)`, and `tau_{h + 128 e_i} v_1 ∈ W_k` because
`h + 128 e_i ∈ High_k` has weight <= wt(h) + 1 and `wt(v_1) <= wt(v) - 1 <= 2 - wt(h)`, so
`wt(h + 128 e_i) + wt(v_1) <= 3`. Both summands are in `W_k`, hence so is `tau_u tau_h v`. ∎

Consequence: `V_4(k)` is the least closed subspace containing `Gen_k`, so `V_4(k) ⊆ W_k`
(given (iii)'s hypothesis), and `W_k ⊆ V_4(k)` holds whenever `S ⊆ V_4(k)` (then
`tau_h v ∈ V_4(k)` by shifting v along the bits of h, each intermediate weight being
<= 3). The computation shows both inclusions for k <= 12.

### 4.3 Why the base is depth 7 and not depth 6

`Free_k^{(6)}(V_4(6))` (levels >= 6 free) has the right dimension at every k but is **not**
`V_4(k)`: for k = 7, 8, 9 the two spaces each have exactly 42 dimensions the other lacks
(`free_extension_from_k0_6_fails` in `analysis.json`), and the free space is not closed.
The stable low part `V_4(k) ∩ M_{levels <= 5}` has dimension 563 = 605 - 42 for every
k = 7, 8, 9 and is contained in `V_4(6)`: the depth-6 closure contains 42 weight-4 elements
produced by shifts that wrap around at level 5 (`u_{i,5}` applied to fall elements with bit
5 set; the fall space has a 33-dimensional piece on the level set {0,..,5}), and
`V_4(k)` for k >= 7 contains 42 elements coupling levels 5 and 6 through honest carries that
are not of the free form `tau_{64 e_i}(S_6 ∩ M_<=3)`. The first depth whose closure is
wrap-free and whose top level is touched by the fall space only through the free translates
`tau_{64 e_i} g`, `tau_{64 e_i} w` is k = 7; from there on everything above is free.
(`analysis.json`: the fall space of `V_4(6)` restricted to the top-level (bit 5) columns has rank
37, that of `V_4(7)` restricted to the bit-6 columns has rank 4, exactly the four free translates.)

## 5. The fall space F_4 = V_4 ∩ M_<=3 (dimension 4k + 87)

### 5.1 Level-set decomposition (exact, direct sum)

For a level set L let `M_L` be the span of the monomials whose level set is contained in L.
`F_4` is the direct sum of its pieces `F_4 ∩ M_L` (verified: the level-restricted subspaces
span `F_4`, `decomposes_into_level_pieces = true`), and the Moebius counts of the exact
level sets are, for every k = 6, 7, 8 computed:

| level set | count | what it is |
|---|---|---|
| {0,1} | 3 | g, tau_{e1} g, tau_{e2} g (V_3) |
| {0,1,2} | 13 | tau_{2e_i} g, tau_{4e_i} g (V_3, 4); w; tau_{e_i} w, tau_{2e_i} w (4); tau_{(3,0)} g, tau_{(0,3)} g (weight-3 generator rows); 2 further combinations |
| {0,1,j}, j >= 3 | 2 each | tau_{2^j e_i} g (V_3) |
| {0,1,2,3} | 23 | core |
| {0,1,2,j}, j >= 4 | 2 each | tau_{2^j e_i} w |
| {0,2,3,4} | 1 | core (multiplier supported on the 8-dimensional cube {0..15}^2 minus (15,15)) |
| {1,2,3,4} | 5 | core |
| {0,1,2,3,4} | 23 | core |
| {0,1,2,3,4,5} | 33 | core (touches level 5; source of the depth-6 wrap-around) |

Sum: `3 + 13 + 2(k-3) + 23 + 2(k-4) + 1 + 5 + 23 + 33 = 4k + 87`. Equivalently (verified as
row spaces for k = 8..12, `fall_space_structure`):

    F_4(k) = F_4(7) (+) (+)_{j = 7}^{k-1} (+)_{i = 1,2} span( tau_{2^j e_i} g, tau_{2^j e_i} w ),

so the k-dependence of the fall space is carried entirely by the free translates of the two
weight-2 elements g and w; the 101-dimensional core lives on levels {0,..,5} and is the
same at every depth >= 6.

### 5.2 The weight-2 fall element w (proof + exact)

`V_4 ∩ M_<=2` is 2-dimensional at every depth computed: it is spanned by g and by

    w = Q_w · g,   Q_w = x1^2 x2^2 + 2t (x1^2 x2 + x1 x2^2) + 2(a + 5t^2) x1 x2
                        - t (11at + 8b + 17t^3)/(a + 3t^2) · (x1^2 + x2^2)
                        + 2 (5a^2 t + 2ab + 15 a t^3 - 2b t^2 + 16 t^5)/(a + 3t^2) · (x1 + x2),

a combination of the nine translates `tau_h g`, `h ∈ {0,1,2}^2` (all generator rows of weight
<= 4), whose product with g has no monomial of weight >= 3. Derivation: `Q g` for Q supported
on `{0,1,2}^2` lives on `{0,..,4}^2`, whose weight->=3 exponent pairs are exactly
(3,1), (3,2), (3,4), (1,3), (2,3), (4,3), (3,3); killing these seven coefficients is a 7x9
linear system over Q(a,b,t) whose kernel is 2-dimensional, spanned by `q = delta_0` (giving g)
and the vector above (`symbolic_w` in the analysis, solved with sympy; `w_formula_on_instance`
evaluates it mod p on the instance: weight 2, in V_4, and g, w span V_4 ∩ M_<=2). The
denominator `a + 3t^2` is the quantity of the kernel-checked Theorem "no constant at degree
three" (`3t^2 + a ≠ 0`): where it vanishes the weight-2 fall degenerates, which is one
component of the exceptional locus. w has 17 monomials, exponents in `{0,1,2,4}^2 ∪ {(3,0),(0,3)}`
minus (2,2) after normalisation, level set {0,1,2}, and it is born in round 0 (a combination of
generator rows).

### 5.3 Iteration depth and birth rounds (exact)

The closure needs 9 rounds of shift-and-fall at every depth 6..12 (rank history at k = 6:
79, 123, 212, 272, 363, 422, 491, 557, 598, 605). The birth-round histograms of the pivots
(`birth_histogram`, `low_birth_histogram` in `k*-D4.json`) are identical for the weight-<=3
pivots at every depth except for the two free families, which are born in rounds 0 (g
translates) and 1 (w translates). The fall space is therefore genuinely iterated: 9 nested
layers of falls, each layer the weight-<=3 part of the shifts of the previous one.

### 5.4 Multiplier filtration (unit ideal, exact)

For a unit ideal g is invertible in `A_k` and every v ∈ V_4 has a unique multiplier
`q_v = v g^{-1}` (computed by the 2-adic Fourier transform over F_p, `multipliers`). With
`B_w` the span of the Dirac masses of weight <= w and `Q = V_4 g^{-1}`:

| k | dim(Q ∩ B_w), w = 0, 1, 2, ... | max weight needed | kappa_4(w) = dim{q ∈ B_w : wt(qg) <= 4} |
|---|---|---|---|
| 6 | 1, 13, 79, 123, 212, 280, 377, 465, 546, 588, 605, 605, 605 | 10 | 1, 13, 79, 123, 212, 280, 377, 465, 583, 717, 781, 793, 794 |
| 7 | 1, 15, 106, 160, 282, 368, 497, 611, 724, 796, 829, 829, 829, 829, 829 | 10 | 1, 15, 106, 160, 282, 368, 497 (w <= 6, matrix capped) |

For the fall part at k = 6, `dim(Q_fall ∩ B_w)` = 1, 13, 18, 34, 43, 59, 72, 88, 103, 111 (w = 0..9)
and `kappa_3(w)` = 1, 13, 18, 34, 43, 59, 72, 88, 103, 222, ... : the fall space is "everything of
weight <= 3 reachable with translates of weight <= 8", and from weight 9 on most weight-<=3
combinations of translates are not derivable. So `Q ∩ B_w = kappa_4(w)` for w <= 7 at k = 6: every weight-<=4 combination of translates of
g by elements of weight <= 7 lies in V_4; from w = 8 on there are weight-<=4 combinations of
translates that are *not* in V_4 (583 > 546), and V_4 needs translates of weight up to 10 =
2k - 2 (not only the weight-<=2 translates of the generator rows). This is the quantitative
form of "V_4 is far from everything reachable": the closure never recovers the deep
translates.

## 6. What a Lean proof of LFD >= 5 (k >= 7, generic) now needs

1. The combinatorial reduction theorem of Section 4.2, which is instance-free: a statement
   about carry-free translations of a weight-filtered subspace of `A_7` into `A_k`. It needs
   only the arithmetic of binary carries and the direct-sum decomposition by high pattern.
2. A base certificate per instance: `V_4(7)` (829 rows, 3473 columns modulo chains), the fact
   `1 ∉ V_4(7)`, and the closedness of `W_8` in `A_8` (1061 rows, 2517 columns). Both are
   exclusion certificates in the sense of `check_closed_exclusion_certificate`, finite and
   replayable; they are larger than the depth-6 objects the kernel could not replay without
   `native_decide`, so the formal route is the reduction theorem plus a trusted (or
   `native_decide`) base check.
3. Genericity of the base: the base data depend on (p, a, b, t). The symbolic form of w shows
   the condition `a + 3t^2 ≠ 0`; the complete exceptional locus of the depth-7 closure (which
   must contain the order-two targets and the (0,0)-target j = 1728 family) is item 2 of the
   bounds program, not settled here. Across the instances computed (Section 8) the pivot set
   of `V_4(7)` in the weight-sorted RREF is identical on all three instances (p = 40961 random curve,
   p = 12289 random curve, p = 12289 j = 1728 curve with the depth-10 probe target), which is the
   combinatorial content a genericity theorem would have to certify.

What is *not* found: a closed form for the 829-dimensional base `V_4(7)` itself (equivalently
for the 101-dimensional core of the fall space on levels {0,..,5}) as explicit polynomial
families in (a, b, t). What is known about it: it is generated from g in 9 rounds, its
weight-<=2 part is span(g, w) with w explicit, its fall space splits into the level pieces of
Section 5.1, its multipliers reach translation weight 10, and the 42-dimensional coupling
between levels 5 and 6 is what forces the base to depth 7.

## 7. Prediction for degree 5 (conjecture, heavy run)

The same mechanism predicts stabilisation of `V_5(k)` once a depth `k_0` is reached whose
closure is wrap-free and whose fall space touches the top level only freely; the free
formula then reads `dim V_5(k) = s_5 + 2(k-k_0) s_4 + C(2(k-k_0),2) s_3 + C(2(k-k_0),3) s_2 + ...`
with `s_d = dim V_5(k_0) ∩ M_<=d`. Because `1 ∈ V_5(k)` for k <= 9 and `1 ∉ V_5(10)`, the base
is at least k_0 = 10 (the first depth where the degree-5 closure separates from the ideal is
the analogue of k = 6 for degree 4, and the first wrap-free base would then be k_0 = 11).
Testing this needs `V_5(11)` (35,443 columns) and the low profile of `V_5(10)`; see
`heavy_runs_needed` in the task result.

## 8. Cross-instance checks (exact)

| instance | dims V_4, k = 6..9 | pivot set of V_4(7) equals primary | V_4(k) = W_k and certificate, k = 8, 9 |
|---|---|---|---|
| `12289:3`: p = 12289, a = 3899, b = 9710, target (9895, 11610), order 6151 | 605, 829, 1061, 1301 | yes (symmetric difference 0) | yes, yes |
| `12289:j1728`: p = 12289, a = 1, b = 0, target (9362, 1126) — the depth-10 probe instance of `experiments/ic-chain-depth10/probe.json` | 605, 829, 1061, 1301 | yes | yes, yes |
| `p1073153-s1` (unit ideal) | 605, 829 (k = 6, 7) | not compared | multiplier filtration above |

The depth-10 probe instance (unit ideal, `dim V_4(10) = 1549`, `LFD >= 6` from the degree-5
closure) therefore has the same degree-4 structure: `V_4(k) = W_k` for its base `V_4(7)`.

## 9. Reproduction

    ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python3 scripts/ic_bounds_fall_export.py --prime 40961 --seed 7 --k 6 7 8 9 10 --k5 6 7 8
    python3 scripts/ic_bounds_fall_analyze.py --prime 40961 --seed 7 --kmax 10          # analysis.json
    python3 scripts/ic_bounds_fall_analyze.py --prime 40961 --seed 7 --kmax 12 --skip-moebius --skip-w \
        --unit-prime 0 --extra-instances '' --analysis-name analysis-k12.json           # k = 11, 12 (recomputed, not stored)

Files: `p40961-s7/instance.json`, `k{6..10}-D4.npz` (RREF basis uint16, pivots, pivot
weights, birth rounds, masks), `k{6..10}-D4.json`, `k{6,7,8}-D5-low.npz` (`V_5 ∩ M_<=4`),
`k{6,7,8}-D5.json`, `analysis.json`, `analysis-k12.json`. Running times on this machine: the
degree-4 closures take 0.6 s (k = 6) to 21 s (k = 12); the degree-5 closure at k = 8
takes 1274 s (about 1.8 GB resident).
