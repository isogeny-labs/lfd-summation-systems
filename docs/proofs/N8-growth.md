# N8: the last fall degree of the two-summand prime chain systems is unbounded in the depth

Status (2026-11-12, revised after review). Evidence labels: **kernel** (Lean-checked), **proof** (written
out in full here), **exact** (certified computation, command and output named), **empirical**, **conjecture**.
Scripts: `scripts/ic_bounds_growth_*.py`; outputs: `experiments/ic-bounds-growth/`.

Reviewed on 2026-11-12 by an independent red-team session: `docs/reviews-N8-growth-redteam-20261112.md`
(no fatal error; findings F1–F14). This revision answers them; the list is in Section 10. The review applies
to the earlier blob `6c98619e…`. A second review, pinned to commit `09dd9760` and blob `0e58d23c` of this file
(`docs/reviews-N8-growth-redteam2-20261112.md`), found no fatal error; its wording findings R2-1 to R2-13 are applied
in §11. Labels: **T1, C1 (a), (b), T4, Proposition 7.2: proof; their reduction and counting steps (Props 2.1, 2.2,
7.1, Lemmas 3.4, 4.1–4.4, Theorem 5.1, Corollary 5.2) are kernel-checked conditional on the hypothesis
`GrowthMain.CurveModel` (Theorem 3.3 for this curve, a written proof, not proved in Lean); T3 and Proposition 7.1:
kernel; C1 (c), existence of unit instances: kernel, unconditional (`GrowthExistence.unit_instance`).**

## 0. Results

| # | statement | label |
|---|---|---|
| T1 | **Growth theorem.** Let K be a field of characteristic 0 or p > 2^(k+1), E: y² = x³ + ax + b nonsingular over K, t ∈ K, **t ≠ 0**. Let D ≥ 2 and k ≥ k_D := D + Λ(D−1) (Λ in Section 4; k_3 = 39, k_4 = 52, k_5 = 64, k_6 = 77; **k_D ≤ 14D + 10** for all D). Then 1 ∉ V_D(k). More precisely every element of V_D(k), read as a polynomial with exponents in [0, 2^k)², has its lift in the prime ideal 𝔭 of Section 1 (the kernel of x_i ↦ ξ_i into the function field of the curve; the form the Lean statement uses). | **proof** (Sections 2–5) |
| C1 | Let K have characteristic 0 or p > 2^(k+1), E nonsingular, t ≠ 0. (a) For D ≥ 2 and k ≥ k_D, 1 is not in the ordinary degree-D closure of the chain presentation (bridge: needs only 2 ≠ 0 in K). (b) If moreover the presentation generates the unit ideal, no e ≤ D is a last fall degree bound; with the kernel upper bound, **⌊(k − 10)/14⌋ + 1 ≤ LFD(k) ≤ k + 2 for k ≥ 38** (and 1 + max{D : k_D ≤ k} ≤ LFD(k) ≤ k + 2 for k ≥ 27). (c) Such unit-ideal instances exist at every depth k ≥ 2 (over F_p, p a prime factor of 2^(2^k) + 1, curves y² = x³ + ax; kernel lemma `UnitInstance.unit_zmod`), so the last fall degree of the two-summand chain systems is unbounded, with lower slope ≥ 1/14 (→ 1/13 asymptotically with the exact Λ). The condition 2^k ∣ p − 1 of the repository is the natural setting but is not used. | **proof** (Corollary 5.2), given T1; the bridge, the upper bound and the unit lemma are **kernel** |
| T2 | Vanishing sums on a curve (function-field abc): Theorem 3.1. Not new: in characteristic 0 it gives the same pairwise bound as Brownawell–Masser (1986), Theorem B; the characteristic-p form with independence over F^p goes back to Voloch (1985). A self-contained proof is written out in Section 3 because neither statement was retrieved in the original (Section 3.1). | **proof** (Section 3) |
| T3 | Reduction: S_D(k) is closed under the degree-D rule of A_k if and only if the splitting property FEI_k(D−1) holds (every instance), and then V_D(k) ⊂ S_D(k), 1 ∉ V_D(k). | **proof** (Propositions 2.1, 2.2); also recomputed both ways on 33 cases (Section 6) |
| T4 | The spaces S_D(k) are free extensions in the depth: S_D(k+1) = ⊕_q x^(2^k q)·(S_D(k) ∩ M_≤D−wt q) for k ≥ D + Λ(D) (char 0 or p > 2^(k+2)). | **proof** (Corollary 5.3) |
| E1 | S_D(k) is closed (equivalently FEI_k(D−1)) first at k = 4, 8, 12 for D = 3, 4, 5, and at every larger depth computed (D = 3, 4: up to k = 9; D = 5: k = 12 only): on the data FEI_k(w) ⇔ k ≥ 4w − 4. The case D = 5, k = 12 is a certificate of 1 ∉ V_5(12) from 12,951 monomials on the curve, without a closure. | **exact** (Section 6.2); the pattern is a **conjecture** |
| E2 | Low-part criterion Crit_D(k) (Proposition 7.1, **proof** that it implies 1 ∉ V_D(k)): in the computed range (D = 3, k ≤ 7; D = 4, k ≤ 8) it holds exactly for k ≥ 4D − 10; at k = 4D − 10 the part of V_D(k) of weight < D **equals** S_{D−1}(k). That Crit_D(k) holds iff k ≥ 4D − 10 in general is a pattern in the data. | **exact** in the computed range (Section 6.2); the general statement **empirical** |
| E3 | Every element of V_4(7) is derived without any wrap-around (it is the part below 2^7 of a closure in the polynomial ring). | **exact** (Section 7.2) |
| E4 | dim S_3(k) = 4k + 98 for k = 7, 8, 9 (instance p = 12289); dim S_4(k) = 3163, 3569, 3855 for k = 10, 11, 12 (instance p = 40961). Of these only dim S_4(10) = 3163 is itself an entry of a certified profile (p_4 of V_5(10)); 4k + 98 is the entry p_3 of the certified V_5(k). Extrapolating S_4 by the free-extension recursion gives 4451 at k = 14, the entry p_4 of V_6(14) (instances p = 65537, 163841). | **exact** dimensions, on different instances; the identifications **empirical** |
| R1 | Refuted: "V_D(k) equals the no-wrap ideal space S_D(k) above the base depth". S_D is strictly larger (Section 6). | **exact** |
| O1 | Not proved: the sharp statements 1 ∉ V_D(4D − 10) (reduced to Crit_D(4D − 10), a statement about the ideal of the curve) and 1 ∈ V_D(4D − 11); any upper bound with slope below 1; the case t = 0. | open (Section 7) |

The hypothesis t ≠ 0 is used once (Lemma 3.2): it makes the two coordinate functions have their
poles at different points with the other coordinate a unit there. The relevant evidence that the
case is different is the run at t = 0 with b ≠ 0 (`ideal-D3-p12289-t0.json`, Section 6.3): three
ideal elements of weight 2 instead of two, S_3(k) closed from k = 5 instead of 4. The repository
record `D4-p12289-t0-X0` (`experiments/ic-bounds-freeext/table.md`), for which no free-extension
base was found, is the doubly special case t = 0, b = 0: there R = (0, 0) is a point of order two,
S_3(x₁, x₂, 0) = (x₁x₂ − a)² (checked with `s3_coeffs` for a = 1) and ξ₁ξ₂ = a, so the two-term
elements x₁^(2^j) x₂^(2^j) − a^(2^j) lie in 𝔭 with weight 2 and the method must exclude it (the m = 2
case of Theorem 3.3 fails). That record is evidence about (t, b) = (0, 0), not about t = 0 alone.

The bound k_D is far from the measured 4D − 10. The proof shows **why** the degree must grow
(a low-weight polynomial vanishing on the curve cannot connect exponents that are far apart, and
wrap-around at the top connects exponents that are 2^k apart), not the constant 4.

---

## 1. Setting

K is a field, k ≥ 2, N = 2^k, A_k = K[x₁, x₂]/(x₁^N − 1, x₂^N − 1) with the monomial basis
x₁^(a₁) x₂^(a₂), 0 ≤ a_i < N. The **weight** of a monomial is wt(a₁) + wt(a₂) (binary digit
sums), the weight of an element is the largest weight in its support, M_≤D is the span of the
monomials of weight ≤ D. The variable u_{i,j} acts as multiplication by x_i^(2^j), 0 ≤ j < k.

g ∈ A_k is the image of

    S_3(x₁, x₂, t) = (x₁−x₂)²t² − 2((x₁+x₂)(x₁x₂+a) + 2b)t + (x₁x₂−a)² − 4b(x₁+x₂),

supported on exponents {0,1,2}², of weight 2. V_D(k) ⊂ M_≤D is the least subspace containing g
and closed under f ↦ u_{i,j}·f for f ∈ V_D(k) of weight ≤ D − 1 (the **degree-D rule**).
This is the model of `docs/proofs/U1-upper-bounds.md` §1; the kernel bridge
(`ChainClosureBridge.chain_mem_closure_iff`) identifies "1 ∈ V_D(k)" with "1 lies in the ordinary
degree-D closure of the chain presentation".

**Lift.** For f ∈ A_k let f̃ ∈ K[x₁, x₂] be the polynomial with the same coefficients on the
exponents in [0, N)². The lift is linear, not multiplicative: for f ∈ A_k and v = x_i^(2^j),

    (v·f)~ = v·f̃ − (x_i^N − 1)·x_i^(2^j − N)·f̃_hi,        f̃_hi = Σ_{a_i ≥ N − 2^j} c_a x^a,          (1.1)

where f̃ = Σ c_a x^a: the monomials with a_i + 2^j ≥ N (those whose bits j, …, k−1 of a_i are all
one, the **top run**) wrap, the others do not.

**The curve.** Let k̄ be an algebraic closure of K, R ∈ E(k̄) a point with x(R) = t, and
F = k̄(E). Put

    ξ₁ = x(P),   ξ₂ = x(R − P)   ∈ F     (functions of P ∈ E),

and let 𝔭 be the kernel of K[x₁, x₂] → F, x_i ↦ ξ_i. It is a prime ideal, it contains S_3 (Semaev:
S_3(x(P₁), x(P₂), x(P₃)) = 0 when P₁ + P₂ + P₃ = O; take P₃ = −R), and it contains neither 1 nor
x_i nor x_i^N − 1 (ξ_i is not constant). No irreducibility of S_3 is needed.
(Numerical check of the Semaev identity for this formula: `ic_bounds_growth_lemmas.py semaev`.)

**Definition.** S_D(k) = { f ∈ M_≤D(A_k) : f̃ ∈ 𝔭 }, the **no-wrap ideal space**. It is a
subspace, g ∈ S_D(k) (k ≥ 2, so g̃ = S_3) and 1 ∉ S_D(k).

---

## 2. Reduction to a splitting property of the ideal 𝔭

**Definition (FEI_k(w)).** For every f̃ ∈ 𝔭 with exponents in [0, N)² and weight ≤ w, every
i ∈ {1, 2} and every j < k: the top-run part f̃_hi = Σ_{a_i ≥ N − 2^j} c_a x^a lies in 𝔭.

**Proposition 2.1.** S_D(k) is closed under the degree-D rule of A_k if and only if FEI_k(D − 1)
holds.

*Proof.* Let f ∈ S_D(k) have weight ≤ D − 1 and v = x_i^(2^j). Then v f ∈ M_≤D, v f̃ ∈ 𝔭, and by
(1.1) (v f)~ ∈ 𝔭 iff (x_i^N − 1)·x_i^(2^j − N)·f̃_hi ∈ 𝔭. As 𝔭 is prime and contains neither
x_i^N − 1 nor x_i, this holds iff f̃_hi ∈ 𝔭 (note x_i^(N − 2^j) divides f̃_hi). The elements f̃ with
f ∈ S_D(k), wt f ≤ D − 1 are exactly the elements of 𝔭 with exponents in [0,N)² and weight ≤ D − 1. ∎

**Proposition 2.2.** If FEI_k(D − 1) holds then V_D(k) ⊂ S_D(k); in particular 1 ∉ V_D(k).

*Proof.* S_D(k) contains g and is closed (Proposition 2.1); V_D(k) is the least such space;
1 ∉ 𝔭. ∎

So an exclusion certificate for 1 can be an **ideal-theoretic** statement about low-weight
polynomials vanishing on one fixed curve, with no closure computation. Sections 3 and 4 prove
FEI_k(w) for k > w + Λ(w).

---

## 3. Vanishing sums on a curve

### 3.1 A Wronskian bound with p-th power coefficients

Let k̄ be algebraically closed of characteristic p ≥ 0, X a smooth projective curve of genus γ
over k̄, F = k̄(X). Let F₀ = k̄ if p = 0 and F₀ = F^p if p > 0. For v₁, …, v_m ∈ F^* put

    H(v₁, …, v_m) = − Σ_P min_i ord_P(v_i)        (sum over the points P of X).

By the product formula H(λv) = H(v) for λ ∈ F^*, and for any i, j

    H(v) = Σ_P max_l (−ord_P(v_l/v_j)) ≥ Σ_P max(0, −ord_P(v_i/v_j)) = ½ Σ_P |ord_P(v_i/v_j)|.     (3.1)

**Theorem 3.1.** Let 𝒮 be a finite set of points of X and v₁, …, v_m ∈ F^* (m ≥ 2) such that

1. v₁ + ⋯ + v_m = 0;
2. v₁, …, v_{m−1} are linearly independent over F₀;
3. for every P ∉ 𝒮 and every i: ord_P(v_i) = 0 if p = 0, and ord_P(v_i) ≡ 0 mod p if p > 0.

Then H(v₁, …, v_m) ≤ ((m−1)(m−2)/2)·(2γ − 2 + |𝒮|).

**Attribution (this theorem is not new; revised after review, F2).**

* *Characteristic 0.* Brownawell–Masser, Vanishing sums in function fields, Math. Proc. Cambridge
  Philos. Soc. 100 (1986), Theorem B, in the form quoted in arXiv:2201.10964 (Theorem 5): F/ℂ of
  genus g, u₁, …, u_k S-units with 1 + u₁ + ⋯ + u_k = 0 and **no vanishing proper subsum**; then
  max_i H(u_i) ≤ binom(k, 2)·(|S| + max(0, 2g − 2)). The original paper was not read here; the
  statement is second-hand. Its hypothesis (no vanishing proper subsum) and its conclusion
  (heights of the u_i, with max(0, 2g − 2)) differ in form from Theorem 3.1 (independence,
  projective height, 2γ − 2). For the use made below they agree: a circuit over the constants with
  m terms v_e has no vanishing proper subsum; dividing by one term v_{e′} gives
  1 + Σ_{e≠e′} v_e/v_{e′} = 0 with k = m − 1, and for g = 1, |S| ≤ 6 the quoted statement is
  H(v_e/v_{e′}) ≤ binom(m−1, 2)·6 = 3(m−1)(m−2), exactly the pairwise bound used in Theorem 3.3.
* *Characteristic p.* The Wronskian method with linear independence over F^p and ordinary
  derivatives is attributed in the literature to Voloch (Diagonal equations over function fields,
  Bol. Soc. Brasil. Mat. 16 (1985)). **The exact statement and constants of Voloch's paper were not
  retrieved**, neither by the author nor by the reviewer. Hypothesis 3 (orders divisible by p
  outside 𝒮, instead of units) is a harmless variant needed for coefficients in F^p.

No novelty is claimed for Theorem 3.1. The proof below is kept so that the file is
self-contained; it is the usual Wronskian argument.

*Proof.* Fix z ∈ F with F/k̄(z) separable and let δ = d/dz, a derivation of F with kernel F₀
(for p = 0 the kernel is the algebraic closure of k̄ in F, which is k̄ = F₀; for p > 0, ker δ is a
subfield with F^p ⊂ ker δ ⊊ F, since δ(z) = 1 ≠ 0, and [F : F^p] = p is prime, so there is no
intermediate field and ker δ = F^p = F₀ exactly. This is what makes hypothesis 2, independence over
F₀, the hypothesis of the Wronskian criterion (a)).

*(a) Wronskian criterion.* For w₁, …, w_r ∈ F put W_δ(w) = det(δ^i w_j)_{0≤i<r, 1≤j≤r}. If
w₁, …, w_r are linearly independent over ker δ then W_δ(w) ≠ 0. Induction on r; r = 1 is clear.
Suppose W_δ(w₁, …, w_r) = 0 and, by induction, W_δ(w₂, …, w_r) ≠ 0. The columns of the matrix are
dependent over F with a nonzero coefficient at the first column (the last r − 1 columns are
independent because their top (r−1) × (r−1) minor is W_δ(w₂, …, w_r)): there are c₂, …, c_r ∈ F with
δ^i w₁ + Σ_{j≥2} c_j δ^i w_j = 0 for 0 ≤ i ≤ r − 1. Applying δ to the equation with index i ≤ r − 2 and
subtracting the equation with index i + 1 gives Σ_{j≥2} δ(c_j)·δ^i w_j = 0 for 0 ≤ i ≤ r − 2, hence
δ(c_j) = 0 as W_δ(w₂, …, w_r) ≠ 0. So w₁ + Σ c_j w_j = 0 is a dependence over ker δ. (No
assumption on the characteristic or on r is used.)

*(b) Change of derivation.* Let P be a point, t = t_P a uniformizer and δ_t = d/dt (t is a
separating element). Then δ = c·δ_t with c = δ(t) = (dz/dt)^(−1), δ^i = c^i δ_t^i + (an F-combination of
δ_t^l, l < i), and by row operations

    W_δ(w₁, …, w_r) = c^(r(r−1)/2)·W_{δ_t}(w₁, …, w_r).

*(c) Local bounds.* In the completion k̄((t)) the derivation δ_t is the formal derivative, so
ord_P(δ_t^l w) ≥ ord_P(w) − l. Expanding the determinant,

    ord_P W_{δ_t}(w₁, …, w_r) ≥ Σ_j ord_P(w_j) − r(r−1)/2          for every P.                 (3.2)

If p > 0 and ord_P(w_j) = p·e_j for all j, write w_j = t^(p e_j)·y_j with y_j a unit of k̄[[t]]; then
δ_t^l(w_j) = t^(p e_j)·δ_t^l(y_j) because δ_t(t^p) = 0, so ord_P(δ_t^l w_j) ≥ ord_P(w_j), and

    ord_P W_{δ_t}(w₁, …, w_r) ≥ Σ_j ord_P(w_j)                         for such P.                 (3.3)

For p = 0 and all w_j units at P, (3.3) holds for the same reason (e_j = 0).

*(d) Conclusion.* Put r = m − 1, n₀ = (m−1)(m−2)/2 and W = W_δ(v₁, …, v_{m−1}), nonzero by (a) and
hypothesis 2. By hypothesis 1 and multilinearity, W = ± W_δ(v_j : j ≠ i) for every i. For each point
P choose i = i_P with ord_P(v_i) = μ_P := min_l ord_P(v_l), and let A_P = Σ_{l=1}^m ord_P(v_l).
By (b), (3.2) for P ∈ 𝒮 and (3.3) for P ∉ 𝒮 (hypothesis 3),

    ord_P(W) ≥ −n₀·ord_P(dz) + (A_P − μ_P) − n₀·[P ∈ 𝒮],

where ord_P(dz) = ord_P(dz/dt_P). Summing over all P: Σ_P ord_P(W) = 0 (W ∈ F^*),
Σ_P ord_P(dz) = 2γ − 2, Σ_P A_P = 0 (product formula), −Σ_P μ_P = H(v). Hence
0 ≥ −n₀(2γ − 2) + H(v) − n₀|𝒮|. ∎

### 3.2 Monomials on the curve of Section 1

From now on X = E, γ = 1, F = k̄(E), and the characteristic of K is 0 or p > 3.

**Lemma 3.2 (divisors).** div(ξ₁) = (P₀) + (−P₀) − 2(O) with x(±P₀) = 0, and
div(ξ₂) = (R − P₀) + (R + P₀) − 2(R). Let 𝒮 = {O, R, ±P₀, R ± P₀} (at most six points). Every
monomial ξ₁^a ξ₂^b is an 𝒮-unit. If **t ≠ 0** then ord_O(ξ₁^a ξ₂^b) = −2a and ord_R(ξ₁^a ξ₂^b) = −2b.

*Proof.* x has a double pole at O and its zeros are the points with x = 0. ξ₂ is x composed with
P ↦ R − P. At P = O, ξ₂ = x(R) = t, and at P = R, ξ₁ = x(R) = t; if t ≠ 0 (and t ≠ ∞, as R ≠ O) both
are units there. ∎

Write ‖n‖_p for the distance from the integer n to pZ if p > 0, and ‖n‖₀ = |n|.

**Theorem 3.3 (circuits are local).** Assume t ≠ 0. Let J be a finite set of exponent pairs
(a, b) ∈ Z²_{≥0}, m = |J|, such that the monomials ξ^e (e ∈ J) are linearly dependent over F₀ and
every proper subset is linearly independent over F₀ (a **circuit**). Then for all (a, b), (a′, b′) ∈ J

    ‖2(a − a′)‖_p + ‖2(b − b′)‖_p ≤ 6(m − 1)(m − 2).

*Proof.* There are φ_e ∈ F₀^* with Σ_{e∈J} φ_e ξ^e = 0 (all coefficients are nonzero by
minimality). Put v_e = φ_e ξ^e. Hypotheses 1–3 of Theorem 3.1 hold with 𝒮 of Lemma 3.2:
independence over F₀ of any m − 1 of the v_e is that of the ξ^e; outside 𝒮 the order of v_e is
that of φ_e, which is 0 for p = 0 (φ_e a nonzero constant) and divisible by p for p > 0 (φ_e a
p-th power). Hence H(v) ≤ ((m−1)(m−2)/2)·(0 + 6) = 3(m−1)(m−2). By (3.1),
2H(v) ≥ |ord_O(v_e/v_{e′})| + |ord_R(v_e/v_{e′})|, and by Lemma 3.2 these orders are congruent to
−2(a − a′) and −2(b − b′) modulo p (equal, for p = 0). ∎

For m = 2 the bound is 0; so two different monomials with exponents below p/2 are never
proportional over F₀.

**Lemma 3.4 (splitting along circuits).** Let M be a finite set of exponent pairs,
M = M₁ ⊔ M₂, such that no circuit J ⊂ M meets both M₁ and M₂. If f = Σ_{e∈M} c_e x^e ∈ 𝔭 with
c_e ∈ K, then Σ_{e∈M₁} c_e x^e ∈ 𝔭.

*Proof.* K ⊂ k̄ ⊂ F₀. In F, y := Σ_{M₁} c_e ξ^e = −Σ_{M₂} c_e ξ^e lies in
span_{F₀}(ξ^e : e ∈ M₁) ∩ span_{F₀}(ξ^e : e ∈ M₂). If y ≠ 0 choose B_i ⊂ M_i with {ξ^e : e ∈ B_i}
independent over F₀ and y in its span; then B₁ ∪ B₂ is dependent and contains a circuit, which
lies in neither B₁ nor B₂, so it meets both M₁ and M₂: contradiction. Hence y = 0. ∎

The reason for working over F₀ instead of the constants: in characteristic p the statement
"a minimal vanishing sum over the constants is local" is **false** (S_3(x₁^p, x₂^p, t) = S_3^p ∈ 𝔭
for a, b, t ∈ F_p is a nine-term vanishing sum with exponent differences p). Over F₀ = F^p these
nine monomials are pairwise proportional, and Theorem 3.3 says only that their exponents agree
modulo p (twice the difference is 0 modulo p).

---

## 4. Counting monomials of bounded weight in a window

For integers n ≥ 0, w ≥ 0 let C(n, ≤w) = Σ_{i≤w} binom(n, i).

**Definition.** Λ(w) is the largest integer ℓ ≥ 1 with 2^(ℓ−1) < 7776·C(2ℓ, ≤w)².

**Lemma 4.1.** Λ(w) is finite, Λ(w) ≥ 4w, the inequality 2^(ℓ−1) < 7776·C(2ℓ, ≤w)² fails for
every ℓ > Λ(w), and Λ(w) ≤ 13w + 23 for every w ≥ 1.

| D | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Λ(D−1) | 25 | 36 | 48 | 59 | 71 | 83 | 95 | 106 | 118 | 142 | 189 | 237 |
| k_D = D + Λ(D−1) | 27 | 39 | 52 | 64 | 77 | 90 | 103 | 115 | 128 | 154 | 205 | 257 |

*Proof.* binom(n+1, i)/binom(n, i) = (n+1)/(n+1−i) ≤ (n+1)/(n+1−w) for i ≤ w, hence
C(n+1, ≤w)/C(n, ≤w) ≤ (n+1)/(n+1−w). Applying this with n = 2ℓ and n = 2ℓ + 1,

    C(2ℓ+2, ≤w)²/C(2ℓ, ≤w)² ≤ [ (2ℓ+1)/(2ℓ+1−w) · (2ℓ+2)/(2ℓ+2−w) ]².

For ℓ ≥ 4w each of the two factors is < 8/7 (7(2ℓ+1) < 8(2ℓ+1−w) ⇔ 8w < 2ℓ+1, and the second factor
is smaller), so the right side is < (8/7)⁴ < 2. So ρ(ℓ) := 2^(ℓ−1)/(7776·C(2ℓ, ≤w)²) is strictly
increasing on ℓ ≥ 4w, and
ρ(ℓ+1)/ρ(ℓ) ≥ 2·(7/8)⁴ > 1.17 there, so ρ → ∞: the set of ℓ with ρ(ℓ) < 1 is finite, and its
part in [4w, ∞) is an initial segment. ρ(4w) < 1 (as C(8w, ≤w) ≥ binom(8w, w) ≥ 8^w, the right side
is ≥ 7776·2^(6w) > 2^(4w−1)), so Λ(w) ≥ 4w and no ℓ > Λ(w) satisfies the inequality. Finally
ρ(13w + 24) ≥ 1: with C(n, ≤w) ≤ (en/w)^w (1 ≤ w ≤ n) it suffices that
13w + 23 ≥ log₂ 7776 + 2w·log₂(2e(13 + 24/w)), i.e. 13w + 10.07 ≥ 2w·log₂(5.437·(13 + 24/w));
this holds for 1 ≤ w ≤ 7 by evaluation (the differences are 7.7, 7.7, 8.0, 8.5, 9.1, 9.7, 10.3)
and for w ≥ 8 because then 2·log₂(5.437·16) < 12.89 < 13. The table and the bound for w ≤ 120 are
recomputed by `ic_bounds_growth_lemmas.py lambda` (**exact** integer arithmetic). ∎

**Lemma 4.2 (window).** Let N = 2^k, p = 0 or p ≥ 2N + 1, σ ≥ 0 an integer, ℓ the least integer
with 2^ℓ ≥ 2σ + 1, and a₀ ∈ [0, N). The set { a ∈ [0, N) : ‖2(a − a₀)‖_p ≤ 2σ } is contained in the
union of at most 6 aligned blocks [c·2^ℓ, (c+1)·2^ℓ).

*Proof.* |2(a − a₀)| < 2N < p, so ‖2(a − a₀)‖_p ≤ 2σ means |a − a₀| ≤ σ or |a − a₀| ≥ (p − 2σ)/2
(only the first for p = 0). The first set is an interval of 2σ + 1 integers. The second is
contained in [0, σ) ∪ (N − σ, N): if a − a₀ ≥ (p − 2σ)/2 ≥ N + ½ − σ then a > N − σ, and
symmetrically. Each of the three intervals has at most 2σ + 1 ≤ 2^ℓ elements, so it meets at most
two aligned blocks. ∎

**Lemma 4.3 (a circuit of low weight is small).** Assume t ≠ 0 and p = 0 or p ≥ 2N + 1. Let J be
a circuit (Theorem 3.3) of exponent pairs in [0, N)² of weight ≤ w, w ≥ 1. Then
σ := 3(|J| − 1)(|J| − 2) satisfies **σ < 2^(Λ(w) − 1)**, and for (a, b), (a′, b′) ∈ J:
|a − a′| ≤ σ or |a − a′| ≥ (p − 2σ)/2 (the second only if p > 0); the same for b.

*Proof.* Let m = |J| and ℓ be the least integer with 2^ℓ ≥ 2σ + 1. Fix (a₀, b₀) ∈ J. By Theorem 3.3
every (a, b) ∈ J has ‖2(a − a₀)‖_p ≤ 2σ and ‖2(b − b₀)‖_p ≤ 2σ, so by Lemma 4.2 a and b lie in at
most 6 aligned blocks each. Inside one pair of blocks a pair (a, b) is determined by the low ℓ
bits of a and of b, and these 2ℓ bits contain at most w ones. Hence m ≤ 36·C(2ℓ, ≤w). If ℓ ≥ 1
then 2^(ℓ−1) < 2σ + 1 = 6(m−1)(m−2) + 1 ≤ 6m² ≤ 7776·C(2ℓ, ≤w)², so ℓ ≤ Λ(w) and
2σ + 1 ≤ 2^ℓ ≤ 2^Λ(w). If ℓ = 0 then σ = 0. The dichotomy for |a − a′| is Theorem 3.3 for the pair
and |2(a − a′)| < p. ∎

**Lemma 4.4 (no straddling).** Let N = 2^k, 1 ≤ λ ≤ k, σ < 2^(λ−1), λ − 1 ≤ j < k, p = 0 or
p ≥ 2N + 1, and a, a′ ∈ [0, N) with a ≥ N − 2^j > a′.

* If a − a′ ≤ σ then wt(a′) ≥ k − λ.
* If a − a′ ≥ (p − 2σ)/2 (p > 0) then wt(a) ≥ k − λ + 1.

*Proof.* First case: a′ = N − 2^j − e with 1 ≤ e ≤ σ < 2^(λ−1) ≤ 2^j. Then
a′ = (N − 2^(j+1)) + (2^j − e), the first summand has the k − j − 1 bits j+1, …, k−1, and
2^j − e ∈ (2^j − 2^(λ−1), 2^j) has all bits λ−1, …, j−1 equal to one: wt(a′) ≥ k − λ. Second case:
a ≥ a′ + N + ½ − σ > N − 2^(λ−1), so the bits λ−1, …, k−1 of a are one. ∎

Lemmas 4.2 and 4.4 are checked by brute force (`ic_bounds_growth_lemmas.py window`, `straddle`).

---

## 5. The growth theorem

**Theorem 5.1.** Let K have characteristic 0 or p > 2^(k+1), E nonsingular, t ≠ 0, w ≥ 1 and
k ≥ w + Λ(w) + 1. Then FEI_k(w) holds. Consequently, for D ≥ 2 and k ≥ k_D = D + Λ(D − 1):

    V_D(k) ⊂ S_D(k)    and    1 ∉ V_D(k).

*Proof.* Let f̃ ∈ 𝔭 have exponents in [0, N)² and weight ≤ w, and fix i (say i = 1; the argument
does not use the symmetry of S_3) and j < k. Let M = supp f̃, M₁ = {(a, b) ∈ M : a ≥ N − 2^j},
M₂ = M ∖ M₁. If M₁ = ∅ there is nothing to prove. Otherwise some a ≥ N − 2^j has weight ≤ w; its
bits j, …, k−1 are one, so k − j ≤ w and j ≥ k − w ≥ Λ(w) + 1.

Suppose a circuit J ⊂ M meets M₁ in (a, b) and M₂ in (a′, b′); so a ≥ N − 2^j > a′. By Lemma 4.3,
with λ = Λ(w): σ < 2^(λ−1), and a − a′ ≤ σ or a − a′ ≥ (p − 2σ)/2. By Lemma 4.4 (λ ≤ k, λ − 1 ≤ j)
wt(a′) ≥ k − λ or wt(a) ≥ k − λ + 1; but both weights are at most w < k − λ. So no circuit meets
both parts, and Lemma 3.4 gives f̃_hi = Σ_{M₁} c_e x^e ∈ 𝔭. This is FEI_k(w). The consequence is
Proposition 2.2 with w = D − 1 ≥ 1. ∎

**Corollary 5.2 (growth of the last fall degree).** Let k ≥ 2, K a field of characteristic 0 or
p > 2^(k+1), E: y² = x³ + ax + b with 4a³ + 27b² ≠ 0 over K, t ∈ K, t ≠ 0, and let F be the normalized
two-summand chain presentation of depth k (`PrimeM2General.inputs k depth a b t`). No condition
2^k ∣ char − 1 and no root of unity in K is assumed (in the repository 2^k ∣ p − 1 is the natural
setting, where the chain variables take values in F_p; it is not used here).

(a) For every D ≥ 2 with k_D ≤ k: 1 is not in the ordinary degree-D closure of F.

(b) If moreover F generates the unit ideal, then no e ≤ D is a last fall degree bound
(`¬ HasLFDBound F e`) for every D ≥ 2 with k_D ≤ k. With LFD(k) the least e such that
`HasLFDBound F e`:

    1 + max{D ≥ 2 : D + Λ(D − 1) ≤ k}  ≤  LFD(k)  ≤  k + 2            (k ≥ 27),
    ⌊(k − 10)/14⌋ + 1                   ≤  LFD(k)  ≤  k + 2            (k ≥ 38).

(c) *Existence.* For every k ≥ 2 there are instances satisfying all hypotheses of (b) at depth k:
let p be a prime factor of 2^(2^k) + 1, X ∈ F_p a non-square, n ∈ F_p a non-square, a = n − X²,
b = 0, t = X. Consequently the last fall degree of the two-summand prime chain systems is
unbounded: sup_k LFD(k) = ∞ over these instances, with LFD(k) ≥ ⌊(k − 10)/14⌋ + 1.

*Proof.* (a) Theorem 5.1 gives 1 ∉ V_D(k) and a subspace S_D(k) ⊂ M_≤D containing g, closed under
the degree-D rule and not containing 1. Under the identification of reduced polynomials with A_k
(U1 §1; `ChainTranslation.normalForm_X_mul_monomial`) its image W satisfies the hypotheses `gen`,
`products`, `one` of `ChainClosureBridge.prime_m2_not_one`, whose only field hypothesis is
`(2 : K) ≠ 0` (true: the characteristic is 0 or > 2^(k+1) ≥ 8). The translation of S_D(k) into such
a W is formalized: `GrowthReduction.not_one_of_fei` applies the bridge with W = `idealSpace` (= S_D(k)), and
`closure_normalForm_mem_idealSpace` gives V_D ⊂ S_D.

(b) `ChainClosureBridge.prime_m2_lfd_gt` (same hypotheses plus the unit ideal) gives
`∀ e ≤ D, ¬ HasLFDBound F e`. The upper bound is the kernel theorem
`SelectorExtension.lfd_le_depth_add_two_of_two_ne` (2 ≠ 0, k ≥ 2, 4a³ + 27b² ≠ 0, any t), so the
least bound exists and is ≤ k + 2. k_2 = 27 makes the first line non-empty from k = 27. For
k ≥ 38 put D₀ = ⌊(k − 10)/14⌋ ≥ 2; then k_{D₀} ≤ 14D₀ + 10 ≤ k (Lemma 4.1), so LFD(k) ≥ D₀ + 1.

(c) 2^(2^k) ≡ −1 mod p, so 2 has order exactly 2^(k+1) modulo p and 2^(k+1) ∣ p − 1 (Dirichlet's
theorem is not needed). Hence p ≥ 2^(k+1) + 1 > 2^(k+1): the characteristic hypothesis of Theorem 5.1
holds automatically. p is odd, so non-squares exist. X ≠ 0 (0 is a square), so t ≠ 0. a ≠ 0, because
a = 0 would make n = X² a square; so 4a³ + 27b² = 4a³ ≠ 0. X³ + aX = X(X² + a) = X·n is a product of
two non-squares of F_p^*, a nonzero square: there is Y ≠ 0 with Y² = t³ + at, and t is not a square.
These are the hypotheses of the kernel lemma `ICFallDegrees.UnitInstance.unit_zmod`
(p prime, 2 ≤ k, 2^(k+1) ∣ p − 1, Y² = t³ + at, Y ≠ 0, ¬IsSquare t), which gives
1 ∈ ⟨inputs k depth a 0 t⟩. There are ((p − 1)/2)² such pairs (X, n). ∎

Numerical check of (c) (`ic_bounds_growth_unit.py`, `unit-existence.json`; **exact**): for
p = 17, 97, 193, 257, 641, every non-square X, a sample of non-squares n and every k ≥ 2 with
2^(k+1) ∣ p − 1 (37,248 instance–depth pairs), the hypotheses hold and S_3(x₁, x₂, X) has no zero
with x₁, x₂ in the group of 2^k-th roots of unity.

Remarks. (i) The earlier version excluded p = 2^k + 1 under the convention 2^k ∣ p − 1; with the
hypothesis stated as p > 2^(k+1) that case needs no mention (and for 27 ≤ k < 2^33 no such prime
exists: k would be a power of two ≥ 32 and the Fermat numbers F_5, …, F_32 are composite). Whether
T1 holds for 2^k < p < 2^(k+1) is open: Lemma 4.4 fails there. (ii) Other sources of unit instances:
the generic system over K(a, b, t), and large primes with 2^k ∣ p − 1 by counting (review, F4).
(iii) Numerically the lower bound is far out: LFD ≥ 8 from k = 90 on. It exceeds the measured
values (5 from k = 6, 6 from 10, 7 at 14) nowhere in the computed range. Its content is the
unboundedness and the linear rate. The slope of the lower bound tends to 1/(1 + c) ≈ 1/13 with
c ≈ 11.99 the root of 4·H₂(1/(2c)) = 1 (k_D/D = 12.9 near D = 60).

**Corollary 5.3 (the ideal spaces are free extensions).** Let k ≥ D + Λ(D), characteristic 0 or
p > 2^(k+2), t ≠ 0. Then

    S_D(k+1) = ⊕_{q ∈ {0,1}²} x^(2^k q)·( S_D(k) ∩ M_≤D−wt q )

(lifts identified), hence dim(S_D(k+1) ∩ M_≤j) = s_j(k) + 2s_{j−1}(k) + s_{j−2}(k) with
s_j(k) = dim(S_D(k) ∩ M_≤j): the profile recursion measured for V_D holds for S_D, **proved**,
from the depth D + Λ(D) on.

*Proof.* ⊃: multiplication by a monomial preserves 𝔭 and adds wt q to the weight. ⊂: let
f̃ ∈ 𝔭 have exponents in [0, 2N)², weight ≤ D, and split its support by the bit k of the first
exponent. A circuit meeting both parts contains a ≥ N > a′ with (Lemma 4.3 at depth k + 1, weight
D, λ = Λ(D)) a − a′ ≤ σ < 2^(λ−1) or a − a′ ≥ (p − 2σ)/2 ≥ 2N + ½ − σ. In the first case
a′ ≥ N − σ > N − 2^(λ−1) has the bits λ−1, …, k−1 equal to one, weight ≥ k − λ + 1 > D. In the second
a > 2N − 2^(λ−1) has weight ≥ k − λ + 2 > D. So the two parts are in 𝔭 (Lemma 3.4); repeat with the
second exponent and divide by x^(Nq) (𝔭 is prime, x_i ∉ 𝔭). ∎

**Remarks.**

1. *Unequal depths.* With depths k₁ ≠ k₂ the reduction needs FEI in both coordinates, hence
   min(k₁, k₂) ≥ k_D, and the characteristic condition becomes characteristic 0 or p > 2^(max(k₁,k₂)+1). This is consistent with LFD ≤ min(k₁, k₂) + 2 (U1, Theorem 1): no growth
   comes from one block alone.
2. *Coset bases* x_i = s_i·u_{i,0}: replace ξ_i by ξ_i/s_i; nothing changes.
3. *t = 0.* Then ξ₂ vanishes at O and ξ₁ at R, ord_O(ξ₁^a ξ₂^b) = −2a + b (or −2a + 2b), and
   exponents far apart in the plane can have equal orders at O and at R. The proof does not
   apply; the orders at the remaining points of 𝒮 could replace them but this was not worked out.
4. *Why the constant is poor.* Theorem 3.1 is used with the bound 3(m−1)(m−2) for a circuit of m
   terms, and m is bounded by the number of all monomials of weight ≤ w in a window. On the data
   the circuits are much smaller and the splitting property holds much earlier (Section 6).
5. *More summands.* For m summands the curve is replaced by E^(m−1) and Theorem 3.1 by a
   higher-dimensional abc theorem (Hsia–Wang); not attempted.

---

## 6. Exact tests

All computations are over F_p with dense exact linear algebra (`scripts/ic_bounds_growth_lib.py`, own
closure engine, checked against the certified dimensions and profiles of the repository:
dim V_4(k) = 163, 386 (= |M_≤4|, 1 ∈ V_4) for k = 4, 5 and 605, 829, 1061, 1301 with low profile
(0, 0, 2, 4k+87) for k = 6, …, 9, and
dim V_3(k) = 2k + 1; every closure is re-checked by a full sweep). Instance unless said otherwise:
p = 12289, y² = x³ + 4021x + 9917, t = 4408 (record `D4-p12289-s12` of `ic-bounds-freeext`).
S_D(k) is the kernel of the evaluation matrix of the monomials on more than 4(2^k − 1) affine
F_p-points of C: a polynomial of bidegree < 2^k vanishing there lies in 𝔭 (Bézout on P¹ × P¹
against the (2,2)-curve; the kernel is the same for every such point set, which the scripts
check). All runs: `ulimit -v 1500000`, one thread, except the two marked "heavy lock".

### 6.1 The theorem cannot be tested in its range; what was tested instead

Theorem 5.1 starts at k = 27 (D = 2) and k = 39 (D = 3); no closure or kernel can be computed
there. The tests below are of (i) every lemma that is a finite statement, (ii) the reduction of
Section 2 on small cases, (iii) the analytic input (Theorem 3.3) in the forms that small cases
allow.

| test | script, output | result |
|---|---|---|
| Semaev identity S_3(x(P), x(R−P), x(R)) = 0 for the repository's formula | `ic_bounds_growth_lemmas.py semaev`, `semaev.json` | 2000 random pairs on 4 curves, all zero |
| Lemma 4.1: table of Λ, Λ(w) ≤ 13w + 23, initial-segment property, w ≤ 120 | `… lemmas.py lambda`, `lemmas.json` | ok (exact integers) |
| Lemma 4.2 (window): all a₀, k ≤ 8, 12 primes p ≥ 2N+1 per k, σ ≤ 60 | `… lemmas.py window` | 60,920 cases, never more than 2 blocks (bound 6) |
| Lemma 4.4 (no straddling): all pairs a ≥ N − 2^j > a′, k ≤ 8, all λ, j | `… lemmas.py straddle` | 74,856 pairs in the two cases, no failure |
| Proposition 2.1: "S_D(k) closed" ⇔ FEI_k(D−1), both computed independently (shift sweep against the echelon basis; top-run parts evaluated on the points) | `ic_bounds_growth_ideal.py`, `ideal-D3-*.json`, `ideal-D4-p12289-s12.json` (main instance, k = 4..9; regenerated by the reviewer, review item 5), `ideal-D4-p12289-unit-s10.json` (unit instance) | equal verdicts in all 33 computed cases (D = 3: k = 2, …, 9 and D = 4: k = 4, …, 9 on the main instance; D = 3: k ≤ 8 and D = 4: k ≤ 8 on the unit instance p = 12289, y² = x³ + x, t = 9362; D = 3: k ≤ 8 at t = 0) |
| Proposition 2.2: whenever S_D(k) is closed, V_D(k) ⊂ S_D(k) and 1 ∉ V_D(k) | same | holds in all cases |
| Theorem 3.3, **over the constants only**, on computed circuits: every reduced-echelon kernel vector of the evaluation matrix is a minimal vanishing sum over the constants; spread ≤ 3(m−1)(m−2)? (four column orders) | `ic_bounds_growth_circuits.py`, `circuits-w2-*.json`, `circuits-w3-*.json` | w = 2, k ≤ 9 and w = 3, k ≤ 8: 1,992 + 56 circuits, no violation; in fact spread ≤ m/2 in every circuit found (the proved bound is quadratic in m) |
| Theorem 3.3, contrapositive, **over the constants only**: m ≤ 7 monomials, exponents ≤ 400, pairwise spread > 3(m−1)(m−2) ⇒ independent on C (random, lattice, powers of two) | `… lemmas.py sparse` | 1,102 sets on 3 curves, none dependent |
| Theorem 3.3 in characteristic p, **dependence over the constants only**, exponents up to 4p (p = 101, 103, 127; points over F_{p²}): m ≤ 4 monomials pairwise far modulo p ⇒ independent over the constants; control: the nine monomials of S_3(x₁^p, x₂^p) are dependent over the constants | `ic_bounds_growth_charp.py`, `charp.json` | 803 sets, none dependent; control detected (rank 8 of 9) |
| Theorem 3.3 **over F^p itself** (added after review, F7): F^p-dependence decided by the Hasse–Wronskian criterion of Theorem 3.1 (a) at points of C(F_{p²}) through the local expansion x₁ = α + s, x₂ = y(s). Controls (p = 131): the 9 terms of S_3 dependent and every 8 independent; the **twisted** relation with exponents e + p·u_e dependent over F^p while its rank over the constants is 9 of 9; a pair e, e + p·u dependent. Tests: all three-term F^p-circuits {0, f₁, f₂}, exhaustive over all residues modulo p (p = 61, 131, and 103 with b = 0); random sets B of m − 1 residues pairwise farther than 6(m−1)(m−2), every residue f scanned as an extension (m = 5, p = 131, 40 sets; m = 6, p = 311, 6 sets) | `ic_bounds_growth_fpdep.py`, `fpdep.json` (18 s) | controls as stated (4 of 4 twisted relations: dependent over F^p, rank 9 over the constants); no three-term F^p-circuit at all; no dependent extension among 686,316 + 580,302 scanned residues |
| Corollary 5.2 (c): the construction of unit instances | `ic_bounds_growth_unit.py`, `unit-existence.json` | 37,248 instance–depth pairs on 5 primes: hypotheses hold, no zero on the grid |

### 6.2 The spaces V_D(k) and S_D(k) at small depth (**exact**)

`ic_bounds_growth_ideal.py 12289 4021 9917 4408 D kmin kmax` and
`ic_bounds_growth_lowpart.py` (same arguments). "low part" = V_D(k) ∩ M_≤D−1. S^♯ is defined in
Proposition 7.1.

**D = 3** (|M_≤3| = 15, 42, 93, 176, 299, 470, 697, 988):

| k | dim V_3 | 1 ∈ V_3 | dim S_3 | V_3 ⊂ S_3 | S_3 closed (FEI_k(2)) | dim low part of V_3 / dim S_2 | Crit_3(k) |
|---|---|---|---|---|---|---|---|
| 2 | 5 | no | 3 | no | no | 1 / 1 (equal) | yes |
| 3 | 7 | no | 16 | yes | no | 1 / 2 | yes |
| 4 | 9 | no | 41 | yes | **yes** | 1 / 2 | yes |
| 5 | 11 | no | 74 | yes | yes | 1 / 2 | yes |
| 6 | 13 | no | 111 | yes | yes | 1 / 2 | yes |
| 7 | 15 | no | 126 | yes | yes | 1 / 2 | yes |
| 8 | 17 | no | 130 | yes | yes | | |
| 9 | 19 | no | 134 | yes | yes | | |

**D = 4** (|M_≤4| = 163, 386, 794, 1471, 2517, 4048):

| k | dim V_4 | 1 ∈ V_4 | dim S_4 | V_4 ⊂ S_4 | S_4 closed (FEI_k(3)) | dim low part of V_4 / dim S_3 | dim S^♯ | Crit_4(k) |
|---|---|---|---|---|---|---|---|---|
| 4 | 163 | yes | 105 | no | no | 93 / 41 | 133 | no |
| 5 | 386 | yes | 270 | no | no | 176 / 74 | 305 | no |
| 6 | 605 | no | 564 | no | no | 111 / 111 (equal spaces) | 605 | **yes** |
| 7 | 829 | no | 1025 | **yes** | no (top-run map of rank 11 under each u_{i,6}) | 115 / 126 | 994 | yes |
| 8 | 1061 | no | 1657 | yes | **yes** | 119 / 130 | 1248 | yes |
| 9 | 1301 | no | 2403 | yes | yes | 123 / 134 (profiles) | | |

How far FEI_k(D−1) fails, independently of the basis (`ic_bounds_growth_obstruction.py`,
`obstruction-D*-p12289-s12.json`): the rank of the top-run map S_{D−1}(k) → F, f ↦ f̃_hi restricted to
the curve, for each (i, j). D = 3: rank 1 at (i, k−1) for k = 2, 3, zero from k = 4. D = 4:
k = 4: ranks 2 and 14 at levels 2 and 3; k = 5: 18 at level 4; k = 6: 21 at level 5; k = 7: 11 at
level 6; k = 8: zero (the same in both coordinates). The counts of failing basis rows printed by
`ic_bounds_growth_ideal.py` and `ic_bounds_growth_fei.py` (for instance 12 at D = 4, k = 7, or 165 at
weight 4, k = 11) depend on the basis and are not invariants; only the verdicts are.

The same numbers are obtained on the unit instance p = 12289, y² = x³ + x, t = 9362 (k ≤ 8;
`ideal-D*-p12289-unit-s10.json`, `lowpart-D*-p12289-unit-s10.json`).

**D = 5: the splitting property FEI_k(4) alone** (`ic_bounds_growth_fei.py 40961 21223 9887 4748 4 k`,
`fei-w4-k*-p40961.json`; k = 11, 12 under the heavy lock, 244 s and 665 s):

| k | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|
| monomials of weight ≤ 4 | 2517 | 4048 | 6196 | 9109 | 12951 |
| dim S_4(k) | 1657 | 2403 | 3163 | 3569 | 3855 |
| FEI_k(4) | no | no | no | no (165 rows fail, top level only) | **yes** |

So S_5(12) is closed under the degree-5 rule: by Proposition 2.2 this is a certificate of
1 ∉ V_5(12) (known from the closures) obtained from a kernel of 12,951 monomials on the curve,
with no closure in M_≤5 (55,455 monomials).

Two further observations on these dimensions (**exact** dimensions; the identifications are
**empirical**, and the numbers compared come from **different instances**: S_3 from p = 12289,
S_4 from p = 40961, the V_5 profiles from p = 12289, 40961, 65537, V_6(14) from p = 65537 and
163841):

* dim S_3(k) = 4k + 98 for k = 7, 8, 9 and dim S_4(12) − dim S_4(11) = 286 = 8·11 + 198 =
  2·dim S_3(11) + dim S_2(11) (taking dim S_3(11) = 142 from the formula): the free-extension
  recursion of Corollary 5.3 is consistent with being in force from k = 7 for weight 3 and from
  k = 11 for weight 4 (one step observed for weight 4), i.e. from 4(w+1) − 9, the base depth of V_{w+1}.
* Comparison with certified profiles. dim S_3(k) = 4k + 98 is the entry p_3 of every certified
  V_5(k) (138, 142, …, 162 for k = 10, …, 16). dim S_4(10) = 3163 is the entry p_4 of V_5(10) **only
  at k = 10**: for k = 11, 12 the profiles of V_5 have p_4 = 3441, 3727, less than
  dim S_4 = 3569, 3855. Continuing the recursion for S_4 from k = 12 gives
  3855 + 294 + 302 = **4451** at k = 14 (an extrapolation, not a computation), which is the entry
  p_4 of the profile (0, 0, 2, 154, 4451, 73587, 448922) of V_6(14)
  (`experiments/ic-bounds-k14/RESULTS.md`); and 4·14 + 98 = 154 is its entry p_3. So in these
  closures the part of weight ≤ D − 2 has the dimension of the no-wrap ideal space S_{D−2}(k),
  and at the first deficient depth also the part of weight D − 1.

Read-off (D = 3 and D = 4, one instance; the D = 5 column is the splitting property alone, on the
instance p = 40961, y² = x³ + 21223x + 9887, t = 4748):

| property | D = 3 | D = 4 | D = 5 | pattern |
|---|---|---|---|---|
| 1 ∉ V_D(k) (computed range: k ≤ 9) | k ≥ 2 | k ≥ 6 | k ≥ 10 (certified in the repository) | from 4D − 10 |
| low part of V_D(k) inside the ideal space S_{D−1}(k) (computed: D = 3, k ≤ 7; D = 4, k ≤ 8) | k ≥ 2 | k ≥ 6 | | from 4D − 10 |
| low part of V_D(k) **equal** to S_{D−1}(k) | k = 2 only | k = 6 only | dimensions agree at k = 10 (3163) | exactly at 4D − 10 |
| criterion of Proposition 7.1 (S^♯ has no new low elements; same range) | k ≥ 2 | k ≥ 6 (fails at 4, 5) | | from 4D − 10 |
| V_D(k) ⊂ S_D(k) (all of V_D lifts into the ideal; k ≤ 9) | k ≥ 3 | k ≥ 7 | | from 4D − 9 |
| S_D(k) closed, i.e. FEI_k(D − 1) | 4 ≤ k ≤ 9 (not at 2, 3) | k = 8, 9 (not at 4, …, 7) | k = 12 (not at 8, …, 11) | from 4D − 8 |

The last row is what Theorem 5.1 proves from k = 39, 52, 64 on.

### 6.3 Attack on the proof: the three weakest steps

1. **Theorem 3.1 in characteristic p (the step (3.3) and the use of F₀ = F^p).** This is the
   step that was written out here without the published statement at hand (Section 3.1). Checks:
   (a) the proof of (3.3) is a one-line computation (δ_t(t^p) = 0), re-derived by the reviewer;
   (b) the statement over the constants is false (Frobenius, control of `charp.py`), so the
   F₀-formulation is necessary; (c) **what the first version did not test** (review, F7): every
   original script decides dependence over the constants only, and the Frobenius control is itself
   a constant-coefficient relation, so dependence with non-constant coefficients in F^p was never
   exercised. `ic_bounds_growth_fpdep.py` now tests F^p-dependence itself (Section 6.1): its
   twisted control is dependent over F^p and has full rank 9 over the constants, i.e. it is
   invisible to the older tests. Result: no three-term F^p-circuit exists at all for the three
   curves (exhaustive over residues, rigorous: every pair is eliminated by a point with nonzero
   Wronskian), and no pairwise-far set has a dependent extension (m = 5, 6; random sets, a search,
   not a proof). Not tested here: circuits with 4 ≤ m exhaustively; m ≥ 7; curve coefficients
   outside F_p (the reviewer's own runs cover m = 7, 8 and coefficients in F_{p²}, with no
   violation). The "three" scan has no positive control of its own (no three-term circuit is
   known); its determinant code is the one validated by the controls.
   What survives if this step were wrong: Theorem 5.1 in characteristic 0 (Brownawell–Masser),
   hence for the generic system over Q(a, b, t), and, for fixed (D, k), for the **generic** system
   over F_p(a, b, t) for all but finitely many p (the closure is a finite sequence of rank
   decisions over Z[a, b, t]; not written out). This fallback says nothing about a specialized
   instance (a, b, t) ∈ F_p³: specialization can change rank decisions on a proper closed set of
   parameters (review, F3).
2. **Lemma 3.2 / Theorem 3.3: the orders at O and R, hypothesis t ≠ 0.** Attack: the same
   computations at t = 0 (`ideal-D3-p12289-t0.json`). The ideal is different there (three
   elements of weight 2 from k = 4 instead of two; S_3(k) closed from k = 5 instead of 4), so the
   hypothesis is not decorative; the conclusion (1 ∉ V_3(k)) still holds on that instance. No
   counterexample to the theorem, and no proof for t = 0.
3. **Lemma 4.3 (counting) and Lemma 4.4 (digits).** Elementary; brute-force checked as listed.
   The constant 36 = 6² in Lemma 4.3 is wasteful (the test never needs more than 2 blocks per
   coordinate) but nothing depends on it except the size of Λ.

A fourth point, not a step of the proof but of its use: Corollary 5.2 passes from the model to
the ordinary closure through the kernel bridge, which needs only 2 ≠ 0 in K (not 2^k ∣ p − 1, as
the first version said). The translation of S_D(k) into a submodule W of reduced polynomials
satisfying the bridge's hypotheses is formalized (`GrowthReduction.not_one_of_fei`, W = `idealSpace`).

---

## 7. What is not proved, refuted candidates, reductions

### 7.1 A weaker sufficient condition that is sharp on the data

**Proposition 7.1 (low-part criterion; proof).** Let D ≥ 3, S′ = S_{D−1}(k) and
S^♯ = S′ + Σ_{i,j} u_{i,j}·S′ ⊂ M_≤D(A_k) (products in A_k, with wrap-around). If

    S^♯ ∩ M_≤D−1 = S′                                                              (Crit_D(k))

then S^♯ is closed under the degree-D rule, V_D(k) ⊂ S^♯, V_D(k) ∩ M_≤D−1 ⊂ S_{D−1}(k), and
1 ∉ V_D(k). FEI_k(D − 1) implies Crit_D(k).

*Proof.* g ∈ S′. If f ∈ S^♯ has weight ≤ D − 1 then f ∈ S′ and u f ∈ S^♯ by definition; so S^♯ is
closed and contains V_D(k). 1 has weight 0 and is not in S′. If FEI_k(D−1) holds then
u S′ ⊂ S_D(k) (Proposition 2.1), so S^♯ ⊂ S_D(k) and S^♯ ∩ M_≤D−1 ⊂ S_D(k) ∩ M_≤D−1 = S′. ∎

In terms of the curve: for f_{i,j} ∈ S′ the element Σ u_{i,j} f_{i,j} has, by (1.1), the lift
−Σ (x_i^N − 1)·x_i^(2^j − N)·(f̃_{i,j})_hi modulo 𝔭. Crit_D(k) says: **if the weight-D monomials of
Σ u_{i,j} f_{i,j} cancel, then Σ_{i,j} (ξ_i^N − 1)·ξ_i^(2^j − N)·(f̃_{i,j})_hi = 0 on the curve.** It
involves only polynomials of weight ≤ D − 1 vanishing on C; no closure.

**Exact (D = 3 for k ≤ 7, D = 4 for k ≤ 8; Section 6.2): in this range Crit_D(k) holds if and only if
k ≥ 4D − 10**, i.e. exactly when 1 ∉ V_D(k), and at k = 4D − 10 the low part of V_D(k) is all of S_{D−1}(k) (so S^♯ = V_D(k) there).
As a general statement beyond this range the equivalence is **empirical**.
For D = 5 the dimensions agree: dim S_4(10) = 3163 = dim(V_5(10) ∩ M_≤4) (the certified profile).
If the pattern persists, dim S_5(14) = 73,587 (the profile of V_6(14)); not computed.

**Reduction A (the sharp lower half).** L_D: 1 ∉ V_D(4D − 10) follows from Crit_D(4D − 10), a
statement about the weight-(D−1) part of the ideal of one curve at one depth. **Conjecture:**
Crit_D(k) holds for every k ≥ 4D − 10 and FEI_k(w) for every k ≥ 4w − 4 (generic instance).
FEI_k(w) for k ≥ 4w − 4 alone would give 1 ∉ V_D(k) for k ≥ 4D − 8, i.e. LFD(k) ≥ ⌊(k + 12)/4⌋ for
unit ideals: the conjectured value ⌊(k + 14)/4⌋ when k ≡ 0, 1 mod 4 and one less otherwise.

What the present proof would need for that: Lemma 4.3 bounds a circuit by all monomials of
weight ≤ w in its window. The sharp statement needs the true size of circuits that contain a
monomial with a top run. On the data (Section 6.1) the spread of a circuit is at most half its
number of terms; a proved bound "spread ≤ c·m" in place of 3(m−1)(m−2) would give
k_D ≈ (1 + c′)·D with c′ the root of c′ = log₂(2e·c′) (about 4.7), still not 4. So the
constant 4 is not reachable by counting monomials in a window; it needs the structure of the
ideal elements of weight ≤ w themselves (their number at depth 4w − 4 against the 4·2^k
conditions imposed by the curve is the natural place to look: dim S_3(k) = 16, 41, 74, 111, 126,
130 for k = 3, …, 8 and 4k + 98 afterwards; dim S_4(k) = 105, 270, 564, 1025, 1657, 2403, 3163,
3569 for k = 4, …, 11).

### 7.2 The wrap-free closure

Let T_D ⊂ K[x₁, x₂] be the least subspace of polynomials of weight ≤ D (integer exponents, digit
sums) containing S_3 and closed under f ↦ x_i^(2^j) f, j ≥ 0, for f of weight < D. Let ρ_k be the
reduction of exponents modulo 2^k, and for a cut N = 2^k write f = Σ_h x^(N h) f_h with f_h of
exponents < N.

**Proposition 7.2 (proof).**
(a) T_D ⊂ (S_3) ⊂ 𝔭, so 1 ∉ T_D. (b) ρ_k(T_D) ⊂ V_D(k).
(c) Say T_D **splits at k** if f ∈ T_D implies f_h ∈ T_D for all h. If T_D splits at k then it
splits at k + 1, the space {f ∈ T_D : exponents < 2^k} is V_D(k) (lifts), and for every k′ ≥ k

    V_D(k′) = ⊕_h x^(2^k h)·( V_D(k) ∩ M_≤D−wt h )      (h over exponent pairs < 2^(k′−k)),

which is the free extension with its profile recursion; in particular 1 ∉ V_D(k′).

*Proof.* (a) the generators and the rule stay in the ideal. (b) ρ_k is a ring homomorphism to A_k
that does not raise weights and sends x_i^(2^j) to u_{i,j} or to 1. (c) *k → k+1:* the slice of f at
the cut 2N with index h is Σ_{q∈{0,1}²} x^(Nq) f_{2h+q}; wt(f_{2h+q}) ≤ D − wt(h) − wt(q), so the at
most two multiplications by x_i^N are allowed by the rule. *Low part:* let f ∈ T_D have exponents
< N and weight < D, and v = x_i^(2^j), j < k. Then v f ∈ T_D and v f = (vf)_0 + x_i^N (vf)_{e_i}; both
slices are in T_D, and the product in A_k is their sum. So the low part of T_D is closed in A_k
and contains g; by (b) it is contained in V_D(k). The displayed formula is the splitting at k for
elements of exponents < 2^(k′) together with the same argument at k′. ∎

**Exact (D = 4, k = 7).** In A_L a shift at level j ≤ L − D never wraps on an element of weight
< D, so the closure T^−(L) of g under those shifts is part of T_D. For L = 10, 11, 12 its part on
exponents < 2^7 has dimension 829 = dim V_4(7) and lies in V_4(7): **every element of V_4(7) is
derived without wrap-around** (`ic_bounds_growth_nowrap.py 12289 4021 9917 4408 4 7 12`,
`nowrap-D4-k7.json`; for L = 7, 8, 9 the dimension is 269, 551, 657). For D = 3, k = 3 the same
holds from L = 5.

**Reduction B.** Proposition 7.2 (c) shows: "T_D splits at k" implies the free-extension law from
k on. The converse needs one more hypothesis (review, F11): **the free-extension law from k
together with V_D(k) ⊂ T_D** (lifts) implies that T_D splits at k (the law then gives
V_D(k′) ⊂ T_D for all k′ ≥ k, every element of T_D lies in some V_D(k′), and its slices are in
V_D(k) ⊂ T_D). The inclusion V_D(4D − 9) ⊂ T_D is known only for D = 4 (E3) and D = 3. With it,
V_D(k) is, above the base, the part below 2^k of one fixed space T_D, which explains the
recursion and the independence of the instance. Splitting is a statement in the polynomial ring
with no wrap-around: T_D is closed under the coefficient maps ∂/∂u_{i,j} for j ≥ k (in the
multilinear coordinates u_{i,j} = x_i^(2^j)). The base depth itself (why 1 enters below it) is not
addressed by Reduction B; it is addressed by Crit_D.

### 7.3 Refuted candidates

| candidate | counterexample (exact) |
|---|---|
| V_D(k) = S_D(k) above the base depth ("the closure is the ideal part of the no-wrap model") | dim V_3(4) = 9 < 41 = dim S_3(4); dim V_4(8) = 1061 < 1657 = dim S_4(8) |
| S_D(k) is closed from the base depth 4D − 9 | D = 4, k = 7: the top-run map of S_3(7) has rank 11 modulo 𝔭 under each top shift u_{i,6} (`obstruction-D4-p12289-s12.json`) |
| V_D(k) ⊂ S_D(k) whenever 1 ∉ V_D(k) | D = 3, k = 2 and D = 4, k = 6: V_D uses the wrap-around (not inside S_D) and still omits 1; only its low part is inside the ideal |
| a minimal vanishing sum over the constants is local, in characteristic p | S_3(x₁^p, x₂^p) (control of `ic_bounds_growth_charp.py`) |
| a two-part support with a gap in one exponent splits an element of 𝔭 (no hypothesis on the other exponent) | the norm G(x₁^(2^j), x₂^(2^j)) = Π_{ζ,ζ′} S_3(ζx₁, ζ′x₂) ∈ 𝔭 (ζ, ζ′ over the 2^j-th roots of unity): its part of highest x₁-exponent is c(x₂)·x₁^e with c ≠ 0, which is not in 𝔭, although every other x₁-exponent is at least 2^j lower. The true statement is Theorem 3.3: far apart forces many terms (**proof**, not a computation) |

### 7.4 Open

* G2/G3, the sharp statements: 1 ∉ V_D(4D − 10) (reduced to Crit_D(4D − 10)) and 1 ∈ V_D(4D − 11).
  Nothing here bears on upper bounds.
* t = 0; more than two summands; proper ideals (Theorem 5.1 gives 1 ∉ V_D(k), which is no
  information when the ideal is proper; a lower bound there needs V_D ≠ I ∩ M_≤D).
* A linear bound for the spread of a circuit in its number of terms (true on all computed
  circuits).

---

## 8. Lean formalization: what is ready

| statement | needs | size |
|---|---|---|
| Propositions 2.1, 2.2, 7.1 (reduction to the ideal; the criterion) for an arbitrary prime ideal 𝔭 of K[x₁, x₂] containing the lift of g and not x_i, x_i^N − 1 | the exponent model of `DepthReduction`/`ChainModelMap`, the lift identity (1.1) | small (a few hundred lines), no computation |
| Lemmas 4.2, 4.4 | digit-sum lemmas on ℕ (`Nat.bitIndices` or a popcount of the repository) | small |
| Lemma 4.1 for one fixed w | `decide`/`norm_num` on integers of about 100 bits | small; the general bound Λ(w) ≤ 13w + 23 needs a real-number estimate |
| Lemma 3.4 | linear algebra over a subfield | small |
| Theorem 5.1 **conditional on Theorem 3.3** (stated as a hypothesis on 𝔭: exponent spread of circuits over a subfield F₀ ⊃ K) | the four items above | medium; this is the realistic kernel target |
| Theorem 3.1/3.3 | function fields of curves, places, the degree of the canonical divisor, Wronskians. Not in the pinned Mathlib | large; not recommended now |

The conditional development exists: `formal/ICFallDegrees/Growth*.lean` (six modules, imported, 35 audit lines; see
`docs/proofs/N8-lean-plan.md`).

## 9. Files

* `docs/proofs/N8-growth.md` (this file); review: `docs/reviews-N8-growth-redteam-20261112.md`.
* `scripts/ic_bounds_growth_lib.py` (engine), `ic_bounds_growth_ideal.py` (S_D(k), closedness, FEI, V ⊂ S),
  `ic_bounds_growth_lowpart.py` (low part, Crit_D), `ic_bounds_growth_fei.py` (FEI alone, blocked elimination),
  `ic_bounds_growth_nowrap.py` (wrap-free closure), `ic_bounds_growth_circuits.py`, `ic_bounds_growth_lemmas.py`,
  `ic_bounds_growth_charp.py` (constants only), and, added after the review, `ic_bounds_growth_fpdep.py`
  (dependence over F^p), `ic_bounds_growth_obstruction.py` (rank of the top-run map),
  `ic_bounds_growth_unit.py` (unit instances).
* `experiments/ic-bounds-growth/*.json` (outputs named above; new: `fpdep.json`, `unit-existence.json`,
  `obstruction-D3-p12289-s12.json`, `obstruction-D4-p12289-s12.json`).

## 10. Changes after review (2026-11-12)

Each claim of the review was checked before the change; none was found wrong.

| finding | what was checked | change |
|---|---|---|
| F1 (Lemma 4.1, ratio bound) | the old inequality was not implied by the quoted binomial ratio; the repair [(2ℓ+1)(2ℓ+2)/((2ℓ+1−w)(2ℓ+2−w))]² < (8/7)⁴ for ℓ ≥ 4w was verified by hand | proof step replaced; `lemmas.py lambda` rerun: Λ, the k_D table and k_D ≤ 14D + 10 unchanged |
| F2 (attribution) | derivation of the pairwise bound from the quoted form of Brownawell–Masser Theorem B redone (m − 1 units, g = 1, \|S\| ≤ 6 gives 3(m−1)(m−2)) | Section 3.1 and row T2: exact quoted statement, its hypothesis and constant; Voloch (1985) for characteristic p; no novelty claimed; stated that Voloch's statement was not retrieved and that Brownawell–Masser is quoted second-hand |
| F3 (fallback for all but finitely many p) | agreed: true for the generic system only | reworded in Section 6.3 |
| F4 (hypotheses of C1) | read `prime_m2_lfd_gt`, `prime_m2_not_one`, `lfd_le_depth_add_two_of_two_ne`, `unit_zmod`; existence construction checked on 5 primes | C1 and Corollary 5.2 restated with the hypotheses used; 2^k ∣ p − 1 removed; remark on p = 2^k + 1 corrected; existence (c) and the two-sided statement added |
| F5 (dim V_4 at k = 4, 5) | `ideal-D4-p12289-s12.json`: profiles (1, 9, 37, 93, 163), (1, 11, 56, 176, 386) | claim restricted to k ≥ 6 |
| F6 (E4 wording) | recorded V_5 profiles: p_4 = 3163, 3441, 3727 at k = 10, 11, 12 | E4 and the bullets of Section 6.2 reworded; instance of each number named |
| F7 (tests were over the constants) | agreed | new script `ic_bounds_growth_fpdep.py` with the twisted control; rows of Section 6.1 relabelled "over the constants only"; Section 6.3 item 1 rewritten |
| F8 (failing-row counts depend on the basis) | agreed (counts are of basis rows) | new script `ic_bounds_growth_obstruction.py`; ranks reported (11 at D = 4, k = 7; 21 at k = 6) |
| F9 (record t = 0, b = 0) | S_3(x₁, x₂, 0) = (x₁x₂ − 1)² for a = 1, b = 0 (`s3_coeffs`) | Section 0 paragraph rewritten; the t = 0, b ≠ 0 run cited as the evidence |
| F10 (slope 1/12.9) | root of 4·H₂(1/(2c)) = 1 is c = 11.99, slope 1/12.99 | "→ 1/13 asymptotically" |
| F11 (Reduction B) | agreed; converse argument written | hypothesis V_D(4D − 9) ⊂ T_D added |
| F12 (T3 label) | Proposition 2.1 is a proof for every instance | row T3 relabelled **proof** |
| F13 (docstring numbering) | grep of the scripts | `lemmas.py`, `charp.py`, `lowpart.py` renumbered |
| F14 (no frozen commit) | — | header states that the review applies to blob `6c98619e…` and that this revision would need a re-review |

## 11. Changes after the second review (2026-11-12)

Findings of `docs/reviews-N8-growth-redteam2-20261112.md` applied here and in the paper: stale statements that nothing
is kernel-checked and that the bridge with W = S_D is not formalized (R2-1) corrected; labels per claim as recommended
(header). In the paper: the abstract and introduction no longer call the hypothesis a published bound (R2-4), state
the characteristic hypothesis (R2-5), give the right role of t ≠ 0 (each coordinate is a unit at the other's pole,
R2-6), state the bound for twice the exponent differences (R2-7), and say where "nonsingular" enters the Lean
statements (R2-8). Also applied in this file (2026-11-12): the scope of E2 is restricted to the computed range D = 3, k ≤ 7 and
D = 4, k ≤ 8, where it is **exact**, and the general equivalence is labelled **empirical** (row E2 and Section 7.1;
R2-9); T1 now says that the lifts lie in the prime ideal 𝔭, the form the Lean statement uses, instead of "vanish on
the curve" (R2-10); Remark 1 on unequal depths states the condition p > 2^(max(k₁,k₂)+1) (R2-11).
