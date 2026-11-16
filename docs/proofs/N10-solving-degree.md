# N10: from the last fall degree to the solving degree

**Status (2026-11-14): written proofs.** Nothing in this file is kernel-checked. The bridge is a published
theorem (Caminata–Gorla); the rest is elementary and built on it. Labels as in the other notes: **kernel**,
**proof** (written out here), **exact** (certified computation), **empirical**, **conjecture**.

## 0. Conventions

- `A = K[u_1, …, u_r]`, `F ⊂ A` finite, `I = (F)`, and `V_{F,D}` the closure of the paper (`sec:lfd`). This is
  Caminata–Gorla's (CG) Definition 1.5 verbatim (multipliers `g ∈ R_{≤ D − deg h}`). `LFD(F)` is CG's `d_F`.
- For a coset factor base `cH` with `c ∉ K`, the generator has coefficients in `K(c)`. Then `K` below means
  `K(c)` and `A = K(c)[u]`. Every statement is over that field.
- **Solving degree** `sd_σ(F)`, CG Definition 1.1: the least `d` such that the row space of the degree-`d`
  Macaulay matrix contains a `σ`-Gröbner basis of `I`. Here "the degree-`d` Macaulay matrix" means the matrix
  after the iteration CG describe in their Section 1. That is Lazard's algorithm: after each reduction, every
  row of degree below `d` is multiplied by every monomial that keeps the degree at most `d`, until nothing new
  appears. It is not the step degree of an F4 run.
- `max.GB.deg_σ(F)` is the largest degree of an element of the reduced `σ`-Gröbner basis of `I`. It is 0 when
  `I = A`, because then the basis is `{1}`.
- `σ` is any degree-compatible term order unless stated otherwise.
- `G ⊆ F` is the **cube subsystem**, in one of two forms:
  - for prime chains, the chain equations `u_v^2 − u_{v'}` and `u_v^2 − 1` (char `K ≠ 2`);
  - for binary descents, the Boolean equations `Z_v^2 − Z_v`.

  In both cases `r` is the number of variables.

## 1. The bridge (published)

**Theorem 1 (CG, arXiv:2112.05579, Thm 3.1).** For every degree-compatible `σ`:
`sd_σ(F) = max{ LFD(F), max.GB.deg_σ(F) }`.

The proof uses CG's equality `rowsp(M_d) = V_{F,d}` (their eq. (1), proved in their reference [15]). That is
exactly why the solving degree here is the mutant-iterated one.

## 2. Gröbner-basis degree on a cube

**Lemma 2 (staircase; proof).** Let `G ⊆ F` be a cube subsystem and `σ` degree-compatible. Then:

1. Every minimal generator of the initial ideal `in_σ(I)` is either some `u_v^2` or a squarefree monomial.
2. If `I ≠ A`, a squarefree minimal generator of degree `d` satisfies `2^d − 1 ≤ dim_K A/I`.
3. Hence `max.GB.deg_σ(F) ≤ max(2, ⌊log_2(dim_K A/I + 1)⌋)`, and also `max.GB.deg_σ(F) ≤ max(2, r)`.

*Proof.*
- **Leading terms of `G`.** Each element of `G` is `u_v^2` minus a polynomial of degree at most 1, so its
  `σ`-leading monomial is `u_v^2`. These leading monomials are pairwise coprime, so `G` is a Gröbner basis
  (Buchberger's first criterion), and `in(G) = (u_1^2, …, u_r^2)`.
- **Part 1.** Since `I ⊇ (G)`, also `in(I) ⊇ in(G)`. Let `μ` be a minimal generator of `in(I)`. If some
  `u_v^2` divides `μ`, then `u_v^2 ∈ in(I)`, and minimality forces `μ = u_v^2`. Otherwise `μ` is squarefree.
- **Part 2.** If `μ` is squarefree of degree `d`, its proper divisors are the `2^d − 1` squarefree monomials
  on proper subsets of its support. None of them lies in `in(I)`: one that did would contradict the minimality
  of `μ`, because `in(I)` is a monomial ideal. So all of them are standard. The standard monomials are a
  `K`-basis of `A/I` (Macaulay), which gives `2^d − 1 ≤ dim_K A/I`.
- **Part 3.** A squarefree monomial has degree at most `r`. For degree-compatible `σ`, an element of the
  reduced basis has the degree of its leading monomial, and those leading monomials are exactly the minimal
  generators of `in(I)`. ∎

**Lemma 3 (dimension = number of points; proof).** If `G ⊆ F` is a cube subsystem, then `A/I` is reduced, and
`dim_K A/I = |Z|`, where `Z` is the zero set of `F` over an algebraic closure `L` of `K`.

*Proof.*
- **`A/(G)` is reduced.** By Lemma 2, `dim_K A/(G) = 2^r` (squarefree standard monomials).
  - *Chains:* `u_{i,0}` determines the other chain variables of row `i`. The row's zeros correspond to the
    roots of `X^{2^k} − 1`, which is separable since `2 ≠ 0` in `K`. So `G` has exactly `∏ 2^k = 2^r` zeros
    over `L`.
  - *Boolean:* the zeros are the `2^r` points of `{0,1}^r`.

  In both cases, the number of zeros equals the dimension, so `A/(G) ⊗ L ≅ L^{2^r}` is reduced, and hence so
  is `A/(G)`.
- **Passing to `I`.** `(A/I) ⊗_K L` is a quotient of `(A/(G)) ⊗_K L ≅ L^{2^r}`. A quotient of a finite
  product of copies of `L` is the product of the copies indexed by the points it keeps. Here those are the
  points of `Z`, so `(A/I) ⊗_K L ≅ L^{|Z|}` and `dim_K A/I = |Z|`. The argument is made after base change on
  purpose: over an imperfect `K`, a reduced algebra can have fewer `L`-points than its dimension.

This is the statement "I is the vanishing ideal of `Z`" of `sec:groupalgebra`, extended to any `F ⊇ G`. ∎

**Corollary 4 (proof).** If `G ⊆ F` and `σ` is degree-compatible, then
`LFD(F) ≤ sd_σ(F) ≤ max{ LFD(F), 2, ⌊log_2(|Z| + 1)⌋ }` and `sd_σ(F) ≤ max{ LFD(F), 2, r }`.

Each of the following implies `sd_σ(F) = LFD(F)`, for every degree-compatible `σ`:

- (a) `Z = ∅`;
- (b) `LFD(F) ≥ max(2, ⌊log_2(|Z|+1)⌋)`;
- (c) `LFD(F) ≥ max(2, r)`.

*Proof.* Combine Theorem 1 with Lemmas 2 and 3. For (a): if `Z = ∅`, then `I = A` (Nullstellensatz), and the
reduced basis is `{1}`. ∎

**Remark (index calculus; not a theorem).** `Z` counts ordered tuples, and `S_{m+1}` is symmetric. So a
random target has `E|Z| ≈ |B|^m / #E(F_p)`, where `|B|` is the factor-base size. That is `m!` times the
number of unordered relations.

- **Systems with no relation.** These are the usual case when `|B|^m` is well below `#E(F_p)`. Case (a)
  applies, and `sd = LFD`.
- **Systems that do give relations.** One relation with distinct `x_i` already gives `|Z| ≥ m!`. Condition
  (b) then needs `LFD ≥ ⌊log_2(|Z|+1)⌋`. For a single relation that threshold is:
  - 2 for `m = 2, 3`;
  - 4 for `m = 4`;
  - 6 for `m = 5`.

  Every exact LFD in the paper with `m ≥ 3` is far above these thresholds.
- **What this does not say.** It is a per-instance condition. It does not say that the solving degree equals
  the LFD for every system an index-calculus run meets.

## 3. Two summands

**Lemma 5 (solution count; proof).** Let `p > 3`, let `E: y^2 = x^3 + ax + b` be nonsingular over `F_p`, let
`t ∈ F_p` be arbitrary, and let `H ⊂ L^*` be the group of `2^k`-th roots of unity with `c ∈ L^*` (coset
`cH`). The two-summand chain system has `|Z| ≤ 2^{k+1}`.

*Proof.*
- **Setup.** The zeros of the chain equations correspond bijectively to pairs `(x_1, x_2) ∈ (cH)^2`. Since the
  normalized generator agrees with `S_3(x_1, x_2, t)` on these points, `Z` is in bijection with
  `{(x_1, x_2) ∈ (cH)^2 : S_3(x_1, x_2, t) = 0}`.
- **The polynomial in `X`.** Write
  `S_3(x_1, x_2, x_3) = (x_1 − x_2)^2 x_3^2 − 2((x_1 + x_2)(x_1 x_2 + a) + 2b) x_3 + (x_1 x_2 − a)^2 − 4b(x_1 + x_2)`.
  As a polynomial in `X = x_2`, the coefficient of `X^2` in `S_3(x_1, X, t)` is `(x_1 − t)^2`.
- **Case `x_1 ≠ t`.** The polynomial has degree 2, so it has at most two roots.
- **Case `x_1 = t`.** Then `S_3(t, X, t) = −4(t^3 + at + b)·X + (t^4 − 2at^2 − 8bt + a^2)`.
  - If `t^3 + at + b ≠ 0`, the polynomial is linear and has one root.
  - If `t^3 + at + b = 0`, substituting `b` makes the constant `9t^4 + 6at^2 + a^2 = (3t^2 + a)^2`. This is
    nonzero, because `t` would otherwise be a double root of `x^3 + ax + b`, contradicting nonsingularity. So
    there is no root.
- **Total.** Summing over the at most `2^k` values of `x_1` gives `|Z| ≤ 2·2^k`. ∎

**Theorem 6 (two summands; proof on top of kernel).** Assume char `≠ 2, 3`, `k ≥ 2`, `E` nonsingular, and any
target abscissa `t`. Then for every degree-compatible `σ`, `sd_σ(F) ≤ max{ LFD(F), k+1 } ≤ k + 2`, and
`sd_σ(F) ≥ LFD(F)`.

*Proof.*
- **Upper bound.** Lemma 5 gives `|Z| ≤ 2^{k+1}`. By Lemma 2, `2^d − 1 ≤ 2^{k+1}` forces `d ≤ k+1`. The
  bound `LFD ≤ k+2` is the kernel-checked theorem `SelectorExtension.lfd_le_depth_add_two_of_two_ne` for
  subgroup bases. Theorem 1 then gives the upper bound.
- **Coset bases.** The kernel theorem takes no coset parameter. A direct computation gives
  `S_3(cx_1, cx_2, x_3; a, b) = c^4 · S_3(x_1, x_2, x_3/c; a/c^2, b/c^3)`. Multiplying an input by the unit
  `c^4` changes neither the ideal nor any `V_{F,D}`. The rescaled curve has discriminant `(4a^3 + 27b^2)/c^6 ≠ 0`.
  So a coset system over `K` is a subgroup system over `K(c)` for a nonsingular curve, and the kernel theorem
  applies there.
- **Lower bound.** `sd_σ(F) ≥ LFD(F)` is part of Theorem 1. ∎

**Lemma 5′ (binary two summands; proof).** Let `K = F_{2^n}`, `E: y^2 + xy = x^3 + ax^2 + b` with `b ≠ 0`, and
let `W_1, W_2` be `F_2`-subspaces of dimensions `k ≤ ℓ`. Then the Boolean two-summand system has
`|Z| ≤ 2^{k+1}`.

*Proof.*
- **The polynomial in `X`.** In characteristic 2,
  `S_3(x_1, x_2, x_3) = x_1^2 x_2^2 + x_1^2 x_3^2 + x_2^2 x_3^2 + x_1 x_2 x_3 + b`.
  For `x_1 ∈ W_1`, the polynomial `S_3(x_1, X, t)` has `X^2`-coefficient `(x_1 + t)^2`.
- **Case `x_1 ≠ t`:** at most two roots.
- **Case `x_1 = t ≠ 0`:** `t^2 X + t^4 + b` is linear, with one root.
- **Case `x_1 = t = 0`:** the polynomial is the constant `b ≠ 0`, with no root.
- **Total.** Boolean coordinates determine `(x_1, x_2)` uniquely through the bases, so `|Z| ≤ 2·|W_1|`. ∎

**Corollary 6′ (proof on top of a written bound).** The Boolean system lives over `F_2`. For binary
two-summand systems with `min(k,ℓ) ≥ 1`, `sd_σ(F) ≤ max{ LFD(F), min(k,ℓ)+1 }`. With the written selector
bound `LFD ≤ min(k,ℓ)+2` of the companion manuscript, this gives `sd_σ(F) ≤ min(k,ℓ)+2`.

## 4. Several summands, direct presentation

**Theorem 7 (proof on top of kernel).** For the direct `m`-summand prime chain system, `r = mk`, and
`δ = deg NF(S_{m+1})`, every degree-compatible `σ` satisfies `sd_σ(F) ≤ max{ LFD(F), mk } ≤ max(mk+1, δ)`.

In particular, wherever `LFD(F) = mk + 1`, the solving degree is `mk + 1` for every degree-compatible order.

*Proof.* Corollary 4 gives `sd_σ(F) ≤ max(LFD, 2, mk)`, which is `max(LFD, mk)` since `mk ≥ 2`. The
normalized generator is squarefree, so `δ ≤ mk`. The kernel-checked finite-quotient theorem (`s = mk`) then
gives `LFD ≤ mk + 1`. ∎ More generally, `sd = LFD` whenever `LFD ≥ mk`.

The same argument holds for binary systems with `r` Boolean coordinates:
`sd_σ(F) ≤ max{ LFD(F), r } ≤ max(r+1, δ)`.

## 5. What this does and does not give

- **What it gives.** Every exact LFD of the paper that meets (a), (b) or (c) is also the exact solving degree,
  for every degree-compatible order. Every proved LFD upper bound becomes a solving-degree upper bound, up to
  the Gröbner-basis terms above.
- **What it does not give.**
  - It does not bound the step degree of F4/F5 implementations. Those follow other selection strategies.
    Magma's "degree of regularity" in PKM is such a step degree. Despite the shared name, it is not the
    degree of regularity of CG Definition 1.4.
  - It does not bound the cost of index calculus. That also needs the number of relation systems and the
    linear algebra.
- **Where the bound is weak.** At large depth with small `p`, `|Z|` is far above `2^{LFD}`, and Corollary 4
  leaves `max.GB.deg` open between its trivial lower bound and `⌊log_2(|Z|+1)⌋`. There, the solving degree is
  a computation (`experiments/ic-solving-degree/`), not a consequence of this note.
- **Novelty.** Not claimed for Lemmas 2 and 3. The argument is elementary; see the literature check in
  `experiments/ic-solving-degree/README.md`.
