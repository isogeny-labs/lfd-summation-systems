# N6: genericity of the depth-7 base. What is true, what is provable, what is proved

> **Note (2026-11-12).** Links to `docs/IC_*.md` point to notes of the original development that are not included in this repository; only `docs/IC_FALL_COUNTEREXAMPLES.md` is. The statements they describe are in the Lean sources and in `paper/`.

Date: 2026-11-08. Branch `claude/lfd-proofs-20261108`. Companion data and commands:
`experiments/ic-bounds-genericity/` (README there). This document settles the five questions of
the N6 task (`docs/IC_BOUNDS_PROGRAM.md`, Section 4) for the two-summand prime chain closure at
degree 4, and states the theorem that can honestly be claimed.

Evidence labels: **kernel** (Lean-checked), **proof** (written out here), **exact** (certified
computation; command in the companion README), **empirical**, **conjecture**. A closure stall is
never a bound; every dimension quoted comes from a closure that was re-verified closed.

## 0. Summary

1. **The recommended theorem (proof, Section 6).** Fix an odd prime `p` and a point
   `P = (a, b, t) ∈ F_p^3` at which the finite depth-7/8 certificate `Cert(P)` passes
   (Definition 1.4: a subspace `S ⊆ M_{<=4}` of the depth-7 algebra containing the depth-7
   generator rows, excluding `1`, whose free extension to depth 8 is closed). Then for **every**
   `μ ∈ F_p^*` and **every** depth `k >= 8`, the system with parameters `(μ^2 a, μ^3 b, μ t)` has
   `1 ∉ V_4(k)` and `dim V_4(k) <= Σ_h dim(S ∩ M_{<= 4 - wt h})` (`= 4k^2 + 172k - 571` for the
   generic profile `(0, 0, 2, 115, 829)`). Hence `LFD >= 5` at every depth `k >= 8` on the whole
   `G_m`-orbit of `P`, for a unit ideal unconditionally and for a proper ideal whenever
   `|Z| < |M_{<=4}| - dim W_k`; for the generic profile the right-hand side is
   `C(2k, <=4) - (4k^2 + 172k - 571) >= 1456` for `k >= 8` and increasing, so this holds for every
   actual instance with few solutions, and the universal bound `|Z| <= 2^{k+1}` (N4a) guarantees it
   for `8 <= k <= 12` only (`2^{k+1}` outgrows the quartic from `k = 13`). This is N5 strengthened from the `H`-orbit (`2^k` points) to the full
   `G_m`-orbit (`p - 1` points) by a twisted-algebra argument (Theorem 3.3). It is unconditional
   given a per-point certificate; the certificate is exact on three instances (item 5), and its
   conclusion was cross-checked at depth 7 on the 288 coset points of Section 7.
2. **The genericity statement itself remains a conjecture** (Conjecture 6.2): "`Cert` passes outside
   a proper Zariski-closed `G_m`-invariant set `E`". It is *true* in the following precise sense
   (Theorem 2.4): if the depth-7 closure over the function field `K = F(A, B)` (t = 1) excludes `1`
   and its depth-8 free extension is closed over `K`, then there is a nonzero polynomial
   `D ∈ F[A, B]`, the common denominator of the generic reduced echelon basis, such that `Cert`
   passes at every point with `D ≠ 0`. What is missing is the *proof of the two generic facts*,
   and the only known route to them is the generic (parametric) computation itself: no
   semicontinuity argument transfers them from the computed points (Section 2.2 shows that both
   strict inequalities `dim V_4(P) > d_gen` and `dim V_4(P) < d_gen` occur at exceptional points,
   so neither closure dimension nor "`1 ∉ V_4`" is monotone under specialisation).
3. **Size of the parametric certificate (exact at k = 6, Section 4).** Along a random line in the
   `(A, B)`-plane the denominator `D` of the generic depth-6 degree-4 basis has degree `>= 480`
   (lcm of 360 reconstructed entries, 21 of which have denominator degree `>= 400`), numerators
   reach degree 538, and the basis has 68,803 nonzero entries; at `k = 7` the basis has 93,011
   nonzero entries and, along the same line, the **same** denominator degree 480 (numerators up to
   533): the degree of the exceptional divisor does not grow from depth 6 to depth 7, consistent
   with the complexity living in the level-`<= 5` core of item 5. A bivariate certificate therefore
   has of the order of `10^10` coefficients at `k = 6` and `k = 7` alike. It is out of reach of a kernel check,
   and of this machine. The fall space alone is much smaller (denominators of degree `<= 115` at
   `k = 6`), but the fall space certifies nothing without a generic rank statement that again needs
   the full generic computation (Section 4.3).
4. **Weighted homogeneity (Section 3).** N3(c) is proved: `(a, b, t) ~ (μ^2 a, μ^3 b, μ t)` for
   `μ ∈ H` have equal closures up to a weight-preserving diagonal map (same dimensions, profiles
   and pivot sets). For `μ ∉ H` the closure of `μ·P` is the closure of `P` in the twisted algebra
   `F[x]/(x_i^{2^k} - c)`, `c = μ^{-2^k}` (the coset `μ^{-1} H` as factor base). **Full invariance
   of the exact closure dimension is false**: two counterexamples from existing exact data
   (Section 3.2). **Full invariance of the lower bound is proved** for `k >= 8` (Theorem 3.3).
   Exact: the depth-7 closure is the *same subspace* for every twist `c` (9 values including
   `c = 0`) on three instances, while at depth 6 the subspace moves with `c` but its dimension,
   profile and pivot set do not; all 192 cosets at `k = 6` and all 96 at `k = 7` (p = 12289)
   give identical dimensions and pivots.
5. **C2 and depth (Section 5).** C2 is in the exceptional set of the depth-7 base (its closure is
   larger, 857 > 829, so `D` must vanish there by Theorem 2.4), yet it has `LFD = 5` at `k = 7` on
   the one target computed. So `E` is a strict over-approximation of the true failure set, the
   true failure sets are not monotone in `k`, and the theorem must be stated with `E` attached to
   the base depth, which Theorem 6.1 does.
6. **Plan (Section 8).** (a) A deterministic proof of the generic facts at `k = 7` needs the
   parametric basis: about `(d+1)^2/2` closures with `d ≈ 2·(deg D + 300)`, i.e. `>= 10^6`
   closures and `>= 10^10` coefficients; not feasible here, barely feasible on a cluster, never
   kernel-checkable. (b) A per-prime exhaustive theorem at `p = 193` (37,249 plane points, about
   5 s each, 52 CPU-hours) gives the complete list `E_193` and an exact, complete statement for
   that prime at every depth `k >= 8`. (c) A probabilistic genericity check (one random point over
   a 61-bit prime) gives confidence, not a theorem. (d) Lean can check N5, Theorem 2.4 and
   Theorem 3.3 as abstract statements with `Cert(P)` as a decidable hypothesis; the hypothesis
   itself needs `native_decide` (829 × 1471 echelon modulo `p` plus 1904 product reductions).

## 1. Setting

### 1.1 The algebra and the closure

Fix a field `F` (for the summation system `F = F_p`, `p` odd; the closure is defined over any
field and its dimensions are invariant under field extension, so the construction makes sense
for every odd `p` even when `2^k ∤ p - 1`: the system is then a summation system over
`F_{p^e}` with `2^k | p^e - 1`, as remarked after N5). For a unit `c ∈ F^*` or `c = 0` let

    A_k^{(c)} := F[x_1, x_2] / (x_1^{2^k} - c, x_2^{2^k} - c),      A_k := A_k^{(1)}.

Its monomial basis is `x^e`, `e = (e_1, e_2) ∈ [0, 2^k)^2`; the weight is
`wt(e) = popcount(e_1) + popcount(e_2)`; `M_{<=D}` is the span of the monomials of weight `<= D`.
The chain variable `u_{i,j}` acts as the operator

    τ^{(c)}_{2^j e_i} x^e = x^{e + 2^j e_i}            if e_i + 2^j < 2^k,
                          = c · x^{e + 2^j e_i - 2^k e_i} otherwise (wrap).

Weights rise by at most one under any `τ`. The generator is
`g(a, b, t) = S_3(x_1, x_2, t)` for `y^2 = x^3 + a x + b`,

    g = x_1^2 x_2^2 - 2t (x_1^2 x_2 + x_1 x_2^2) + t^2 (x_1^2 + x_2^2) - 2(t^2 + a) x_1 x_2
        - (2at + 4b)(x_1 + x_2) + (a^2 - 4bt),

with coefficients in `Z[a, b, t]`; its exponents are `<= 2`, so it is the same element of every
`A_k^{(c)}`, `k >= 2`. The generator rows at degree `D` are `Gen_k = {τ_μ g : wt μ <= D - 2}`.

**Lemma 1.1 (generator rows never wrap for k >= 4).** For `wt μ <= 2` and `k >= 4` every
monomial of `τ_μ g` has exponents `<= 2^{k-1} + 2^{k-2} + 2 < 2^k`, so `τ^{(c)}_μ g = τ_μ g` for
every `c`. *Proof.* `2^{k-1} + 2^{k-2} + 2 <= 2^k - 1` iff `3 <= 2^{k-2}`. ∎

**Definition 1.2 (closure).** `V_D^{(c)}(P)` is the least subspace `W ⊆ M_{<=D}` with
`Gen ⊆ W` and `τ^{(c)}_v(W ∩ M_{<=D-1}) ⊆ W` for all `2k` variables `v`. It is the union of the
chain `V^{(0)} = span Gen`, `V^{(r+1)} = V^{(r)} + Σ_v τ_v(V^{(r)} ∩ M_{<=D-1})`. Write
`V_D(P) = V_D^{(1)}(P)`, `low_d = dim(V_D ∩ M_{<=d})`, and "profile" for `(low_0, …, low_D)`.
Modulo the chains this is the degree-`D` closure of the repository (equal to the ordinary closure
for `D >= 2`; repository premise, cross-checked in items 2–4).

**Lemma 1.3 (comparison).** If `W ⊆ M_{<=D}` contains `Gen` and is closed
(`τ_v(W ∩ M_{<=D-1}) ⊆ W` for all `v`), then `V_D ⊆ W`. *Proof.* Induction on the chain. ∎

Column order: monomials sorted by weight, heaviest first, the constant last. In a reduced row
echelon basis with respect to this order, a row whose pivot has weight `<= d` lies in `M_{<=d}`,
so the rows with pivot weight `<= d` span exactly `W ∩ M_{<=d}`, and `1 ∈ W` iff the constant
column is a pivot (a vector lies in the row space iff it equals the combination of the rows given
by its pivot-column entries; `e_const` has zero pivot entries unless the constant is a pivot,
and a row with the last column as pivot is `e_const` itself).

### 1.2 The certificate

**Definition 1.4 (depth-7/8 certificate Cert(P)).** A subspace `S ⊆ M_{<=4} ⊆ A_7` such that
(i) `Gen_7(P) ⊆ S`; (ii) `1 ∉ S`; (iii) the free extension
`W_8(S) := ⊕_{h ∈ {0,128}^2} τ_h (S ∩ M_{<= 4 - wt h})` is closed in `A_8` (c = 1). The natural
choice is `S = V_4(7)(P)`, and then (i) is automatic. `Cert(P)` passes if such an `S` exists;
item 5 verified it with `S = V_4(7)` on three instances (829 rows, 1471 columns; `W_8` has 1061
rows and 2517 columns), independently re-verified in its `VERIFY.md`.

**N5 (reduction theorem; proof in `experiments/ic-bounds-fall/README.md` §4.2, checked by its
verifier).** With `W_k(S) := ⊕_h τ_h(S ∩ M_{<= 4 - wt h})`, `h` supported on levels `>= 7`,
`wt h <= 2`: (i) `W_k ⊇ Gen_k` for `k >= 7` [uses only (i) of Cert and Lemma 1.1]; (ii)
`1 ∈ W_k ⟺ 1 ∈ S`; (iii) if `W_8` is closed in `A_8`, `W_k` is closed in `A_k` for every
`k >= 8`. Consequently `V_4(k)(P) ⊆ W_k(S)`, `1 ∉ V_4(k)(P)` and
`dim V_4(k)(P) <= Σ_h dim(S ∩ M_{<= 4 - wt h})` for all `k >= 8`.

## 2. Semicontinuity (question 1)

### 2.1 What is semicontinuous and what is not

Let `R` be an integral domain of parameters (`Z[a, b, t]`, `F[a, b, t]`, or `F[A, B]` for the
slice `t = 1`), `K = Frac R`, and `V_D^K` the closure over `K` with weight-sorted RREF basis `B`
(unique). For a ring homomorphism `ρ : R → L` into a field ("the point `P`"), `V_D(ρ)` is the
closure over `L` of the specialised generator.

The closure is iterative and each round intersects with `M_{<=D-1}`. For a *fixed* polynomial
matrix the rank is lower semicontinuous (it can only drop at a point) and the dimension of the
kernel of a fixed map is upper semicontinuous (it can only jump up). The low part
`V^{(r)} ∩ M_{<=D-1}` is the kernel of the projection to the heavy columns *restricted to a space
that itself may shrink*, so its dimension is `rank(V^{(r)}) - rank(proj V^{(r)})` with both terms
dropping at special points: **it is not monotone in either direction**, and neither is the
closure dimension.

**Proposition 2.1 (both strict inequalities occur; exact, item 2 data at k = 6, p = 193).**
Generic `dim V_4 = 605`, `low_3 = 111`, `low_2 = 1` at degree 3.
* Larger closure through more falls: order-two targets (`t^3 + at + b = 0`): `607`, `low_4 = 113`
  (the extra fall is `q = x_1 x_2 - t(x_1 + x_2) - a - 2t^2` with `g = q^2`, N3(b)).
* Smaller closure despite more falls: the tangent line `3t^2 + a = 0`: `dim(V_3 ∩ M_{<=2}) = 3 > 1`
  (N3(a), proof) and yet `dim V_4 = 603 < 605`.
* Smaller closure: C1, `602`; special points `600`, `599`.

**Toy mechanism.** Two rows `r_1 = (1, s)`, `r_2 = (1, 0)` with the first column heavy and the
second light: generically `rank 2` with a one-dimensional light part `(0, 1)`; at `s = 0` the
rank is `1` and the light part is `0`. Conversely `r_1 = (s, 1)`, `r_2 = (0, 0)`: the light
part is `0` generically and `(0, 1)` at `s = 0`. Both effects occur in the closure, and the
second feeds the next round, which is why exceptional points can have larger closures.

### 2.2 Constructibility; "1 ∉ V_4" is neither open nor closed

**Proposition 2.2 (proof).** The functions `P ↦ dim V_D(P)`, `P ↦ profile`, `P ↦ pivot set` and
the predicates `1 ∉ V_D(P)`, `Cert(P)` are constructible in `P` (finite Boolean combinations of
polynomial equations). *Proof.* The closure algorithm with a fixed pivoting rule performs
finitely many eliminations; each branch is determined by the vanishing or non-vanishing of
finitely many polynomials in the parameters (the candidate pivots), and on each branch the result
is a fixed matrix of rational functions. The parameter space is therefore a finite disjoint union
of locally closed strata on each of which all of the above are constant. ∎

**Proposition 2.3 ("1 ∉ V_4" is not open; kernel + exact).** On the line `b = t = 0`
(`j = 1728`, target `(0, 0)`) the kernel-checked K4 gives `1 ∈ V_4` at every depth when
`a^{2^k} ≠ 1`, while for `a ∈ H` item 2 finds `1 ∉ V_4` (`dim V_4 = 730`, `|Z| = 64` at
`k = 6`, re-run here: `experiments/ic-bounds-genericity`). Inside the line, `{1 ∉ V_4}` is the
finite set `a ∈ H`, closed and not open; its complement in the line is open in the line but not
closed in the ambient space. So `{1 ∉ V_4}` is constructible, contains a dense open set if the
generic statement holds, and is neither open nor closed. The same example shows that
`dim V_4` is not `G_m`-invariant (Section 3.2). "`W_8(V_4(7)(P))` closed" is constructible by
Proposition 2.2; no openness claim is made for it.

### 2.3 The one statement that does hold: specialisation of a parametric basis

**Theorem 2.4 (specialisation; proof).** Let `B̃ ∈ R^{d × N}` (columns = `M_{<=D}` of `A_{k_0}`
in weight order), `Δ ∈ R \ {0}`, and pivot columns `π_1 < … < π_d` with
`B̃[i, π_j] = δ_{ij} Δ` and `B̃[i, c] = 0` for `c < π_i`. Assume the identities in `R^N`:

* (I1) for every generator row `γ` (entries in `R`): `Δ γ = Σ_j γ[π_j] B̃_j`;
* (I2) for every `i` with `wt π_i <= D - 1` and every variable `v`:
  `Δ τ_v(B̃_i) = Σ_j (τ_v B̃_i)[π_j] B̃_j`;
* (I3) the constant column is not among the `π_j`.

Then for every `ρ : R → L` with `ρ(Δ) ≠ 0`, the rows `ρ(B̃_i)/ρ(Δ)` form a weight-sorted RREF
basis with pivots `π` of a subspace `W_ρ ⊆ M_{<=D} ⊗ L`; `W_ρ ⊇ Gen(ρ)`; `W_ρ` is closed;
`1 ∉ W_ρ`. Hence `V_D(ρ) ⊆ W_ρ`, `1 ∉ V_D(ρ)`, `dim V_D(ρ) <= d`, and
`dim(V_D(ρ) ∩ M_{<=e}) <= #{j : wt π_j <= e}` for every `e`.

*Proof.* The specialised rows have entry `ρ(Δ) ≠ 0` at their own pivot and `0` at the other
pivots and to the left of their pivot, so after division by `ρ(Δ)` they are in RREF with pivots
`π` and are independent. Specialising (I1) and dividing by `ρ(Δ)` expresses every generator row
as a combination of the rows. By the column order, `W_ρ ∩ M_{<=D-1}` is spanned by the rows with
`wt π_i <= D - 1`, and specialising (I2) shows that every variable shift of such a row is a
combination of the rows; closedness follows by linearity. By (I3) the constant column is not a
pivot, so `1 ∉ W_ρ`. Lemma 1.3 gives `V_D(ρ) ⊆ W_ρ`, and the low parts are contained in the low
parts. ∎

**Where the data come from.** If `V_D^K` has RREF basis `B` with common denominator `Δ`
(`B̃ = Δ B`), then (I1) and (I2) hold over `K` because `V_D^K` contains `Gen` and is closed, and
clearing denominators gives the `R`-identities; (I3) holds iff `1 ∉ V_D^K`. So Theorem 2.4
answers question 1 in the required form: **the exceptional set of the parametric certificate is
`{Δ = 0}`, the zero set of the common denominator of the generic RREF basis**, equivalently the
product (up to multiplicity) of the pivots met by any elimination that produces it, and
closedness and exclusion hold at every point where `Δ ≠ 0`. Two remarks:

* The iterative nature of the closure plays no role in the certificate: only the final basis and
  the identities matter. The branch structure of Proposition 2.2 is replaced by a single
  polynomial `Δ`.
* `Δ(P) ≠ 0` does **not** imply `V_D(P) = W_ρ`: the specialised generator rows may be dependent
  (toy example above) and then `V_D(P) ⊊ W_ρ`. For the lower bound only the containment is needed.

**Corollary 2.5 (direction of the inequalities).** If `dim V_D(P) > dim V_D^K`, or
`dim(V_D(P) ∩ M_{<=e}) > dim(V_D^K ∩ M_{<=e})` for some `e`, then `Δ(P) = 0`. If
`dim V_D(P) < dim V_D^K`, nothing follows: such a point may or may not lie in `{Δ = 0}`. (The
tangent line does: the weight-2 rows of the generic basis have denominator `A + 3` along the
sampled line, exact in Section 4.2 and as predicted by the closed form of `w = Q_w g` with
denominator `a + 3t^2`, so `(A + 3) | Δ`. Whether C1 does is open; Section 8 says what would
decide it.)

### 2.4 Why nothing stronger is available

The lattice limit `Λ̄_P` of the generic closure along a curve through `P` (reduction of
`V^K ∩ R_{(P)}^N` for a discrete valuation) has the generic dimension and contains no information
about `V_D(P)` in general: Proposition 2.1 gives points with `V_D(P) ⊄ Λ̄_P` (order two,
`607 > 605`) and the toy example gives `V_D(P) ⊊ Λ̄_P`. In particular **the generic facts
(`1 ∉ V_4^K(7)`, `W_8^K` closed) cannot be inferred from their truth at computed points**: rank
is lower semicontinuous, so a point certifies lower bounds on generic ranks, but the generic facts
are upper bounds on the closure (containment in a space excluding `1`). This is the structural
reason why N6 needs a generic computation and not more instances.

## 3. Weighted homogeneity (question 3)

### 3.1 H-invariance and the twist correspondence

For `μ ∈ F^*` let `φ_μ : A_k → A_k` be the linear bijection `x^e ↦ μ^{-(e_1 + e_2)} x^e`
(`e_i ∈ [0, 2^k)` as integers). It preserves every `M_{<=d}`, is diagonal (so it preserves pivot
columns of row spaces) and fixes the line of the constant.

**Proposition 3.1 (proof).** (a) `g(μ^2 a, μ^3 b, μ t) = μ^4 φ_μ(g(a, b, t))` (weighted
homogeneity of `S_3` of weight 4 with `x_i, t` of weight 1, `a` of weight 2, `b` of weight 3;
identity in `Z[a, b, t, x_1, x_2, μ]`, re-checked by sympy in item 2's VERIFY).
(b) `φ_μ^{-1} τ_u φ_μ = μ^{|u|} τ_u^{(c)}` with `c = μ^{-2^k}`, `|u| = 2^j`: without wrap the
exponent sum grows by `|u|`, with wrap by `|u| - 2^k`.
(c) Hence for `k >= 4` and every `D`:

    V_D(μ^2 a, μ^3 b, μ t) = φ_μ ( V_D^{(c)}(a, b, t) ),        c = μ^{-2^k}.

*Proof of (c).* By (a), (b) and Lemma 1.1, `Gen(μ·P) = μ^4 μ^{|ν|} φ_μ(τ^{(c)}_ν g(P)) =
(scalars) · φ_μ(Gen(P))`, and a subspace `W` is `τ^{(c)}`-closed iff `φ_μ(W)` is `τ`-closed
(conjugation by `φ_μ` sends each `τ_u` to a nonzero multiple of `τ_u^{(c)}`). The least closed
space containing `φ_μ(Gen(P))` is therefore `φ_μ` of the least `τ^{(c)}`-closed space containing
`Gen(P)`. ∎

**Corollary 3.2 (N3(c); proof).** For `μ ∈ H` (`μ^{2^k} = 1`), `c = 1` and
`V_D(μ·P) = φ_μ(V_D(P))`: equal dimensions, profiles and pivot sets at every degree and depth.
This is the H-orbit invariance of item 2 (its proof is the same map). The general case says:
**scaling by `μ ∉ H` is the same system with the coset `μ^{-1} H` as factor base**, i.e. the
chain with `x_{i,k} = c` instead of `1`.

### 3.2 Full invariance of the exact closure dimension is false

**Proposition 3.3 (exact + kernel).** `dim V_4(μ·P) = dim V_4(P)` fails in general.
* `P = (a, 0, 0)` with `a ∈ H`, `μ` with `μ^2 ∉ H` (`c ≠ ±1`): `dim V_4(P) = 730`, `1 ∉ V_4`,
  but `dim V_4(μ·P) = |M_{<=4}| = 794`, `1 ∈ V_4` (K4, kernel; both values re-run here at
  `p = 193`, `k = 6`: `a = 13`, `μ^2 a = 52`). In the twisted picture: the solution set
  `x_1 x_2 = a`, `x_i^{2^k} = c` is empty when `c^2 ≠ a^{2^k} = 1`.
* A C2 point with `t ≠ 0`: `(A, B) = (20, 0)` at `p = 193`, `k = 6`. On the sheet `t = 1`
  (`a = 20, b = 0, t = 1`): `|Z| = 10`, `dim V_4 = 784`; on the sheet `t = 96` (`a = 5, b = 0`):
  `|Z| = 32`, `dim V_4 = 762` (item 2 data, `plane/p193-k6-B0.jsonl` and
  `data/p193-k6-a5-b0.jsonl`). Both are saturated (`V_4 = I_4`, `dim = 794 - |Z|`), so the
  *class* is invariant but the *dimension* is not: it moves with `|Z|`, which depends on the coset.

So the statement "the class depends only on `(A, B)`" (item 2, empirical) cannot be upgraded to
invariance of the closure. What is invariant is weaker and what is provable is the bound:

### 3.3 Full invariance of the lower bound (proof)

**Theorem 3.3.** Suppose `Cert(P)` passes with the subspace `S` (Definition 1.4). Then for every
`c ∈ F` (including `c = 0`) and every `k >= 8`, `W_k(S)` is closed in `A_k^{(c)}` and contains
`Gen_k`; hence `V_4^{(c)}(k)(P) ⊆ W_k(S)` and, for every `μ ∈ F^*`,

    1 ∉ V_4(k)(μ^2 a, μ^3 b, μ t),     dim V_4(k)(μ·P) <= dim W_k(S),
    dim(V_4(k)(μ·P) ∩ M_{<=e}) <= dim(W_k(S) ∩ M_{<=e}).

*Proof.* Run through the proof of N5(iii) with `τ^{(c)}`. Elements of `W_k` decompose by high
pattern `h` (levels `>= 7`); it suffices to treat `τ^{(c)}_u τ_h v`, `v ∈ S ∩ M_{<= 3 - wt h}`.
*Case `u = 2^j e_i`, `j <= 6`:* the operators are multiplications in the commutative algebra
`A_k^{(c)}`, so `τ^{(c)}_u τ_h v = τ^{(c)}_h (τ^{(c)}_u v)`. Since `v` has exponents `< 128`,
`τ^{(c)}_u v = τ_u v = v_0 + τ_{128 e_i} v_1` is computed without any wrap (exponents `< 256`), and
the N5 argument (closedness of `W_8` in `A_8`, where this shift never wraps either, and uniqueness
of the high decomposition) gives `v_0 ∈ S ∩ M_{<= wt v + 1}`, `v_1 ∈ S ∩ M_{<= wt v - 1}`. Then
`τ_h v_0` involves no wrap (exponents of `v_0` are `< 128`, `h` is `0 mod 128` blockwise), and
`τ^{(c)}_h τ_{128 e_i} v_1 = c^ε τ_{h'} v_1` with `h' = (h + 128 e_i) mod 2^k` supported on levels
`>= 7` of weight `<= wt h + 1` (a carry run clears bits and sets one; a wrap clears a bit), so both
summands lie in `W_k`. *Case `j >= 7`:*
`τ^{(c)}_u τ_h v = c^ε τ_{(h + u) mod 2^k} v` with `ε ∈ {0, 1}`, and `(h + u) mod 2^k` is again
supported on levels `>= 7` with weight `<= wt h + 1` (a carry run clears bits and sets one; a wrap
clears a bit), so the result lies in `τ_{h'}(S ∩ M_{<= 4 - wt h'})` ⊆ `W_k`. `Gen_k ⊆ W_k` is
N5(i), which uses only `Gen_7 ⊆ S` and Lemma 1.1 (the generator rows are the same in every
`A_k^{(c)}`). Lemma 1.3 gives `V_4^{(c)}(k)(P) ⊆ W_k(S)`; `1 ∉ W_k(S)` by N5(ii). Finally
Proposition 3.1(c): `V_4(k)(μ·P) = φ_μ(V_4^{(c)}(k)(P)) ⊆ φ_μ(W_k(S))` with `c = μ^{-2^k}`,
and `φ_μ` preserves weights and the constant line. ∎

**Remark (why k >= 8 and not 7).** At depth 7 the wrap at level 6 does occur inside `A_7`, so
closedness of `S` in `A_7` does not transfer to `A_7^{(c)}` by this argument. Section 7 shows
exactly that `V_4(7)` is nevertheless the same subspace for every `c` on the three instances.

**Remark (the exceptional set can be taken `G_m`-invariant).** The conclusion of Theorem 3.3 at
`μ·P` follows from `Cert(P)`; so the set of points where the *conclusion* holds for all `k >= 8`
contains `G_m · {Cert passes}`, and the exceptional set of the recommended theorem is contained in
the `G_m`-invariant set `∩_μ μ · {Cert fails}`. In the normalized plane this is a statement about
`(A, B) = (a/t^2, b/t^3)` for `t ≠ 0`; the plane `{t = 0}` is a closed set treated separately
(its own two-parameter family, with the K4 line inside it).

### 3.4 Exact evidence on the closure itself (Section 7 for commands)

* **Depth 7 is wrap-insensitive (exact, 3 instances).** For `P` = `p40961-s7`, `p12289-r3`,
  `p12289-j1728` and `c ∈ {1, 0, 2, 3, -1, four random}`: `V_4^{(c)}(7)(P)` is the **same
  subspace** of `M_{<=4}` as `V_4(7)(P)` (row spaces equal; 829 rows, profile `(0,0,2,115,829)`,
  closed, `1 ∉`). Consequence (proof from this exact datum): `V_4(7)` is closed under both
  `τ^{(0)}` and `τ^{(1)}`, hence under `τ^{(c)} = τ^{(0)} + c(τ^{(1)} - τ^{(0)})` for every `c`, so
  `V_4^{(c)}(7) ⊆ V_4(7)` for *every* `c`, and `V_4(7)(μ·P) ⊆ φ_μ(V_4(7)(P))` for every
  `μ ∈ F_p^*`; equality holds for the tested `c`.
* **Depth 6 is not (exact).** Same instances and `c`: dimension `605`, profile `(0,0,2,111,605)`
  and pivot set identical for all `c`, but the subspace differs from `V_4(6)` for every `c ≠ 1`
  (neither inclusion holds). These are the 42 wrap-around dimensions of item 5 §4.3, now seen as
  a `c`-dependent family of equal dimension. The same holds at the order-two point
  `(A, B) = (22, 170)`, `p = 193` (607/835), and the tangent point `(190, 23)` (603/828).
* **Full coset sweep (exact).** `p = 12289`, instance `p12289-r3`: all 192 cosets of `H` at
  `k = 6` and all 96 at `k = 7` give the same dimension, profile and pivot set as `μ = 1`
  (`coset_test.json`). Item 2's 1268/1268 test covered only the 3 cosets that exist at `p = 193`.
* **Twisted depth-8 certificate (exact, consistent with Theorem 3.3).** On `p40961-s7`, `W_8`
  built from `V_4(7)` is closed under `τ^{(c)}` for `c ∈ {1, 0, 2, 777}` and equals
  `V_4^{(c)}(8)` for each (`extension_twist_test.json`).

## 4. The certificate a Lean lower bound needs, and its size (question 2)

### 4.1 The identities

By Theorem 2.4 and N5, a parametric certificate for all `k >= 8` on `{Δ ≠ 0}` consists of
`(B̃, Δ, π)` at depth 7 and three families of identities in `R^N`:

* (I1) `2·7^2 + 7 + 1 = 106` identities, one per generator row, each over `N_7 = 1471` columns
  (the generator rows have 9 nonzero entries, so each identity has at most 9 terms on the right);
* (I2') closedness of the free extension at depth 8: for every row `τ_h B̃_i` of `W_8` with
  `wt h + wt π_i <= 3` (`4·8 + 87 = 119` rows for the generic profile) and each of the 16
  variables: `Δ τ_v τ_h B̃_i = Σ (τ_v τ_h B̃_i)[τ_{h'} π_j] τ_{h'} B̃_j`, i.e. `1904` identities over
  `N_8 = 2517` columns (the sum has at most two blocks `h'`, and as many terms as `τ_v τ_h B̃_i`
  has pivot-column entries);
* (I3) the constant is not a pivot of `B̃` (a statement about `π`, free).

Closedness of `S` itself in `A_7` is not needed (N5 never uses it), so no depth-7 product
identities are required for `k >= 8`. For a theorem that also covers `k = 7` one adds the
`14 × 115` depth-7 identities (with the wraps at level 6), or uses the exact wrap-insensitivity of
Section 3.4 per instance. Theorem 2.4's proof applies verbatim to the reduced spanning set
`{τ_h B̃_i}` of `W_8`: its pivot columns `τ_h π_i` each carry a single nonzero entry (`Δ`) among
the rows, so at a point with `Δ ≠ 0` the specialised rows are independent, their low part is
spanned by the rows with `wt h + wt π_i <= 3`, and (I2') specialises to closedness.

### 4.2 Sizes (exact at k = 6; k = 7 as computed)

Nonzero entries of the RREF basis at a random point (68,803 at two points over `p = 40961` and
`p = 33554393`, 68,801 at a third over `p = 12289`: the generic support is the set of entries that
are nonzero rational functions, and at a point a few numerators vanish; the counts below are
representative, not invariants):

| k | shape | nonzero entries | density | nonzero entries in the 4k+87 fall rows |
|---|---|---|---|---|
| 6 | 605 × 794 | 68,803 | 14.3 % | 8,841 (111 rows × 299 light columns) |
| 7 | 829 × 1471 | 93,011 | 7.6 % | 8,893 (115 rows × 470 light columns) |

Degrees along a random line `(A, B) = (A_0 + s A_1, B_0 + s B_1)`, `t = 1`, over
`p = 33554393`, with every entry reconstructed as a reduced rational function of `s` by Newton
interpolation and the extended Euclidean algorithm (`line_degree_k6.json`; 1200 sample points, the
pivot set is the generic one at all 1200, each reconstruction verified with a slack of at least
170 points):

| k | sampled entries | reconstructed | max deg numerator | max deg denominator | deg lcm of denominators | `D` at the order-two point of the line | `D` at the tangent point |
|---|---|---|---|---|---|---|---|
| 6 | 360 (60 rows × 6) | 360 | 538 | 480 | **480** | 0 | 0 |
| 7 | 360 (60 rows × 6) | 360 | 533 | 480 | **480** | 0 | 0 |

Per pivot weight at `k = 6`: the two weight-2 rows have denominators of degree `<= 1` (`A + 3`
for `w`, as the closed form predicts); the 28 sampled weight-3 (fall) rows have denominator
degrees with quartiles 7 / 18 / 79 and maximum 115; the 30 sampled weight-4 rows have quartiles
30 / 82 / 361 and maximum 480, with 21 of 180 sampled entries at degree `>= 400`. The vanishing
of the reconstructed `D` at the intersection of the line with `1 + A + B = 0` and with `A = -3`
is exactly what Corollary 2.5 predicts (larger closure, resp. the denominator of `w`), and is an
independent check of the reconstruction.

For a generic line the degree of `D|_line` is the total degree of `D(A, B)`, so (exact lower
bound from the sample, equality likely) **`deg D >= 480` at `k = 6`**. A bivariate polynomial of
degree 480 has up to 116k coefficients and the numerators (degree up to 538) up to 145k. The
parametric basis at `k = 6` therefore holds of the order of `68,803 × 1.4·10^5 ≈ 10^10`
coefficients (tens of gigabytes modulo a word-sized prime, more over `Z`), and the identities
(I1), (I2') are products of such polynomials. At `k = 7` the counts are larger
(93,011 entries, the same degrees: `k = 7` used 1400 sample points on the same line, with a slack of
at least 370 points per reconstruction; per pivot weight the denominator quartiles are 4 / 18 / 78,
max 111 for the fall rows and 15 / 73 / 99, max 480 for the weight-4 rows). The equality of the
degrees at `k = 6` and `k = 7` (480 and 480; numerators 538 and 533) says that the exceptional
divisor of the base does not grow with the depth, as the free-extension structure predicts; whether
the two lcm polynomials coincide was not checked. This settles question 2 negatively as far as a Lean-checkable or
even machine-storable explicit `E` is concerned.

### 4.3 Can the fall space alone carry the certificate?

The fall rows are far smaller (denominator degree `<= 115` against 480 at `k = 6`), and `V_4` is
spanned by `Gen ∪ F_4 ∪ {τ_v f : f ∈ F_4}` (every element of the closure is produced from a fall
by one shift or is a generator row). So one may ask for a certificate consisting of a parametric
fall space `F̃` only. At a *point* `P` this works: with `W_P := span(Gen(P) ∪ F(P) ∪ τ_v F(P))`,
`W_P` is closed as soon as `W_P ∩ M_{<=3} = F(P)`, which is a rank condition
(`rank` of the weight-4 projection of the spanning rows `= dim W_P - dim F(P)`), cheap to check
at `P`. For the *generic* statement one needs `dim W^K - rank_K(proj_4) = dim F^K`, an **upper**
bound on a generic rank (`dim W^K <= rank proj_4 + dim F^K`), and upper bounds on generic ranks
are not certified by points (Section 2.4); they need the explicit dependencies among the
`1716 + 115` spanning rows, i.e. the full generic basis again. So the fall space reduces the
*per-point* certificate but not the *generic* one. (It does suggest the heavy-run order in
Section 8: the parametric fall space at `k = 7` is a feasible object, about `10^8` coefficients,
and would make the per-point check a rank computation.)

## 5. C2 and depth dependence (question 4)

Facts (item 2 and its VERIFY): at `p = 193`, `k = 6`, C2 points are saturated
(`V_4 = I_4`, `LFD = 4`); at `p = 257`, `k = 7`, C2 points have `dim V_4 = 857 > 829` and the one
target computed to degree 5 has `V_5 = I_5`, `LFD = 5`; at `k = 5` no C1/C2 deficiency exists.

Consequences for the form of the theorem:

1. **C2 ⊂ E for the depth-7 base**, by Corollary 2.5: at C2 the depth-7 closure is larger than
   the generic one, so `Δ_7 = 0` there (and `Cert` with `S = V_4(7)` may or may not pass; whether
   the free extension of the 857-dimensional space is closed is a finite check per point, not
   done). Yet `LFD = 5` holds at the computed C2 target. So **`E` is a sufficient-condition set,
   not the failure set**, and the theorem must be phrased "outside `E`", never "exactly off `E`".
2. **The true failure sets move with `k`** (`{LFD <= 4}` contains C2 at `k = 6` and, on the
   evidence, not at `k = 7`), while the theorem's `E` is attached to the base depth and is the
   same for all `k >= 8` (N5 removes `k`). The theorem is therefore necessarily of the form
   "there is a depth-independent `E_7` such that for all `k >= 8` …", with `E_7` an
   over-approximation at every depth. A sharper `E_k` per depth would need a base at depth `k`.
3. **Per prime or over `Q`?** If the generic facts hold over `Q(A, B)` with `Δ ∈ Z[A, B]`
   primitive, Theorem 2.4 applies to every prime at once with `E_p = {Δ ≡ 0 mod p}`, whose
   `F_p`-points vary with `p` exactly as C1 and C2 are observed to (a fixed polynomial of degree
   `>= 19` acquires different rational points modulo different primes; the verifier's remark that
   190 points carry no evidence for a *curve* is correct, and nothing here decides it). For
   primes `p <= deg Δ`, `Δ mod p` may vanish on all of `F_p^2`, and a separate computation over
   `F_p(A, B)` is needed: the statement is then genuinely per prime. With `deg Δ >= 480` at
   `k = 6` this concerns every prime used so far (`193, 257, 449, 577`), and Schwartz–Zippel
   gives no nontrivial bound on `|E_p|` at those primes. The complete-plane count at `p = 193`
   (about 760 deviating points of 37,249) is consistent with a divisor having four absolutely
   irreducible components over `F_193` plus components without rational points.
4. Nothing in the depth dependence contradicts Theorem 6.1: the theorem is about `k >= 8` and a
   fixed base; C2 is simply inside `E`.

## 6. The theorem

**Theorem 6.1 (proof: N5 + Theorem 3.3; the hypothesis is exact on the instances named).**
Let `p` be an odd prime, `k >= 8`, and `P = (a, b, t) ∈ F_p^3` such that `Cert(P)` passes with a
subspace `S` of profile `(s_0, …, s_4)`. Then for every `μ ∈ F_p^*` the two-summand chain system
with parameters `(μ^2 a, μ^3 b, μ t)` at depth `k` satisfies

    1 ∉ V_4(k),    dim V_4(k) <= s_4 + 2(k-7) s_3 + C(2k-14, 2) s_2 + C(2k-14, 3) s_1 + C(2k-14, 4) s_0.

In particular, if the ideal is the unit ideal then `LFD >= 5`; if the ideal is proper then
`LFD >= 5` whenever `|M_{<=4}| - |Z|` exceeds the right-hand side (then `V_4(k) ≠ I_4`, so the
closure has not reached the ideal at degree 4; N4a). For the generic profile `(0, 0, 2, 115, 829)`
the right-hand side is `4k^2 + 172k - 571` and the condition reads
`|Z| < C(2k, <=4) - 4k^2 - 172k + 571`, a quantity equal to 1456 at `k = 8`, 2747 at `k = 9` and
increasing; the universal bound `|Z| <= 2^{k+1}` (N4a) guarantees it for `8 <= k <= 12`
(`C(2k,<=4) - 2^{k+1} = 2005, 3024, 4148, 5013, 4759` against `1061, 1301, 1549, 1805, 2069`),
and fails to do so from `k = 13` on, where the instance's actual `|Z|` must be used (it is small for
every instance computed: 2, 2, 2, 6, 26 at `k = 6..10` on the primary instance). The same conclusion
holds over every extension `F_{p^e}` (closure dimensions are extension-invariant).

The hypothesis is exact for the three instances of item 5 (`p = 40961` random curve;
`p = 12289` random curve; `p = 12289`, `j = 1728` probe target) and, by Corollary 3.2 and the
twisted depth-8 check, for all their `G_m`-orbits; Section 7 adds the 288 coset points, where the
depth-7 closure was recomputed (the depth-8 certificate on those points follows from Theorem 3.3).

**Conjecture 6.2 (genericity; what N6 asked for).** For every odd prime `p` (and over `Q`) there
is a nonzero polynomial `Δ_7 ∈ F_p[A, B]` (resp. `Z[A, B]`) such that `Cert(A, B, 1)` passes with
`S = V_4(7)` of profile `(0, 0, 2, 115, 829)` whenever `Δ_7(A, B) ≠ 0`. Equivalently (Theorem 2.4):
`1 ∉ V_4^{K}(7)` and `W_8(V_4^K(7))` is closed over `K = F_p(A, B)`, with the generic profile.
`E := {t = 0} ∪ {Δ_7(a/t^2, b/t^3) = 0}` then contains `{1 + A + B = 0}`, `{A = -3}`, the cusp,
the special point `(32, 103)`, and the points of C2 (all have larger depth-6 or depth-7
closures); whether it contains C1 is open. **Status: conjecture**, supported by exact
computations at every point tried (dimension 829, identical pivot set and identical support of the
RREF basis across primes 193, 12289, 40961, 33554393 and across all cosets), and *not* provable
from such points (Section 2.4).

**What is proved beyond N5.** (a) Theorem 2.4: the exact shape of the open condition and the
identities; (b) Proposition 3.1 and Corollary 3.2: N3(c) and its extension to the coset/twist
correspondence; (c) Theorem 3.3 and Theorem 6.1: the all-depth bound on the full `G_m`-orbit of
every certified instance, a factor `(p-1)/2^k` more parameter points per certificate than N5
(320 at `p = 40961`, `k = 7`); (d) Propositions 2.1–2.3 and 3.3: the negative results that fix
the form of any genericity theorem (no semicontinuity, not open, not `G_m`-invariant in
dimension).

## 7. Computations (exact; all commands in `experiments/ic-bounds-genericity/README.md`)

All runs: `ulimit -v 3000000`, one BLAS thread, under 20 minutes each, peak RSS below 0.5 GB,
outputs 0.2 MB. Engine: `experiments/ic-bounds-genericity/genericity_engine.py` (twisted shift
on top of `csd/ic_bounds_fall.py`'s grid and `csd/elim.fp_echelon`; it reproduces the recorded
rank histories `[79, 123, …, 605]` and `[106, 160, …, 829]` of item 5 entry by entry).

| file | content | result |
|---|---|---|
| `twist_test.json` | `V_4^{(c)}` for 9 values of `c`, 6 instances, `k = 6, 7` | depth 7: same subspace for all `c`; depth 6: same dim/profile/pivots, different subspace for `c ≠ 1` |
| `coset_test.json` | all cosets of `H` at `p = 12289`: 192 at `k = 6`, 96 at `k = 7` | all identical to `μ = 1` (dim, profile, pivots) |
| `extension_twist_test.json` | twisted closedness of `W_8`, cross-prime pivots, density | `W_8` closed for `c ∈ {1,0,2,777}` and `= V_4^{(c)}(8)`; pivot sets equal across 4 primes and random points; densities above |
| `line_degree_k6.json` | 1200 closures along a line, 360 entries reconstructed | `deg D >= 480`, numerators `<= 538`, `D = 0` at the order-two and tangent points |
| `line_degree_k7.json` | 1400 closures along the same line at `k = 7`, 360 entries reconstructed | `deg D >= 480` again, numerators `<= 533`, `D = 0` at the order-two and tangent points |

Negative controls built into the runs: the depth-6 twist shows that the engine does detect
`c`-dependence where it exists; the reconstructed `D` vanishes where theory says it must.

## 8. A feasible plan (question 5), with costs

Costs are for this machine (one closure: 0.65 s at `k = 6`, 1.45 s at `k = 7`, 3.4 s for the
depth-8 free-extension check; 15 GB RAM, under 1 GB free disk).

1. **Prove the generic facts at `k = 7` by a parametric computation (the only route to
   Conjecture 6.2).** Evaluation/interpolation over `F_p[A, B]` (two variables; the slice `t = 1`
   suffices by Theorem 3.3 for `k >= 8`): the numerators and denominators have degrees of the
   order of `d ≈ 500` (480/538 at `k = 6`, 480/533 at `k = 7`), so a dense
   bivariate reconstruction needs about `(d_num + d_den + 2)^2 / 2 >= 5·10^5` closures (`>= 8`
   CPU-days at `k = 7`), storage of `9.3·10^4 × 5·10^5` values (`>= 180` GB) and a reconstruction
   of `10^10` coefficients. Not feasible here; feasible on a cluster with a compiled engine;
   the resulting object is not kernel-checkable. Symbolic (fraction-free) elimination over
   `F_p[A, B]` on an 829 × 1471 matrix with entries of degree in the hundreds is slower still.
   **Recommendation: do not pursue as a Lean target.**
2. **Per-prime exhaustive theorem at `p = 193` (feasible heavy run, about 52 CPU-hours, 4 cores:
   13 h).** For each of the 37,249 points `(A, B)` of the plane: `V_4(7)`, `1 ∉`, free extension,
   depth-8 closedness (Definition 1.4 with `S = V_4(7)`). Output: the exact list `E_193` of
   failing points and the profiles. Combined with Theorem 6.1 this gives the complete exact
   statement "for `p = 193`, every `(a, b, t)` with `t ≠ 0`, `(a/t^2, b/t^3) ∉ E_193`, every
   `k >= 8`: `LFD >= 5`" (unit and proper, by the dimension bound), and it decides whether C1,
   the tangent line and the special points pass the per-point certificate (Corollary 2.5 leaves
   that open). It is also the strongest *exact* statement about `E` obtainable at a small prime,
   since `deg Δ >= 480 > 193` makes Schwartz–Zippel vacuous there. Store only the failing points
   and the class values (a few hundred kB). The same at `p = 257` costs 92 CPU-hours.
3. **Probabilistic genericity check (cheap, not a proof).** At one random point over a 61-bit
   prime (engine change: the current kernel is limited to `p < 2^31`; use two 31-bit primes and
   CRT, or a 64-bit kernel): `dim V_4(7) = 829`, generic pivot set, certificate passes. By
   Schwartz–Zippel this is the generic behaviour with error at most `deg/p ≈ 10^3/2^61` *if* the
   degree of the relevant minors is bounded as estimated. Worth doing as a sanity check of
   Conjecture 6.2 over large primes; label the result empirical.
4. **The parametric fall space `F_4(7)` (optional, feasible on a larger machine).** About
   `(155 + 115 + 2)^2 / 2 ≈ 3.7·10^4` closures at `k = 7` (15 CPU-hours), `8,893` entries ×
   `≈ 10^4` coefficients `≈ 10^8` coefficients. It turns the per-point certificate into a rank
   check (Section 4.3) and would exhibit `E ∩ {fall space degenerates}` as the zero set of a
   degree-`<= 115` polynomial, the natural candidate for "the" exceptional curve. It does not
   prove Conjecture 6.2.
5. **What Lean could check without `native_decide`.** (i) N5 (instance-free carry combinatorics);
   (ii) Theorem 3.3 (the twisted version of the same combinatorics plus the diagonal map
   `φ_μ`); (iii) Theorem 2.4 as an abstract statement about a matrix over a commutative ring with
   a unit-pivot structure and the identities (I1)–(I3); (iv) Theorem 6.1 in the form
   "`∀ p k a b t μ, 8 <= k → Cert p a b t → LFD-lower-bound (μ^2 a) (μ^3 b) (μ t) k`" with `Cert`
   a decidable predicate. The kernel cannot evaluate `Cert` at an instance (829 × 1471 modulo
   `p` plus 1904 reductions over 2517 columns), so the instance-level fact enters either as
   `native_decide` or as a hypothesis. There is no Lean route to Conjecture 6.2: its certificate
   is the object of item 1.

## 9. The strongest honest statement

* **Proved (this document):** Theorem 6.1. Every point passing the finite depth-7/8 certificate
  gives `LFD >= 5` at every depth `k >= 8` on its whole `G_m`-orbit; the exceptional set of any
  genericity theorem is the zero set of the common denominator of the generic basis (Theorem 2.4)
  and can be taken `G_m`-invariant.
* **Exact:** the hypothesis of Theorem 6.1 at the three item-5 instances and the 288 coset
  points; wrap-insensitivity of the depth-7 base; `deg D >= 480` at `k = 6`.
* **Conjecture:** Conjecture 6.2, the generic statement. It cannot be obtained from more
  instances; it needs the parametric computation of Section 8, item 1, whose size (`>= 10^10`
  coefficients) puts it beyond this machine and beyond any kernel check. The explicit `E` is out
  of reach; a complete per-prime `E_193` is within reach (Section 8, item 2).
