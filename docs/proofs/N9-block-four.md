# N9: block-of-four lemmas toward LFD ≤ ⌊(k+14)/4⌋ + O(1) — conjectural

**Status (2026-11-12): conjectural work. No bound with slope below 1 is proved here.** The measured
slope 1/4 (exact LFD 4, 5, 6, 7 at k ≤ 5, 6–9, 10–13, 14) is **not** established by anything in this file.
Evidence labels as in the repository: **kernel** (Lean-checked), **proof** (written out in full here),
**exact** (certified computation, command and output named), **empirical** (pattern in data),
**conjecture**. Nothing in this file is kernel-checked; no Lean file was written.

Scripts: `scripts/ic_bounds_block4_lib.py`, `scripts/ic_bounds_block4_mixed.py`,
`scripts/ic_bounds_block4_maps.py`. Outputs: `experiments/ic-bounds-block4/` (README there lists the
commands). Every closure below was followed by a deterministic full sweep (every shift of every final
basis row of weight < D reduces to zero): field `sweep_failures`, 0 in every record.

## 0. Results

| # | statement | label |
|---|---|---|
| R1 | Reduction: the block lemma BF(D): 1 ∈ V_D(4D−11) ⇒ 1 ∈ V_{D+1}(4D−7), for all D ≥ 4, together with U_4, gives LFD(k) ≤ ⌊(k+14)/4⌋ for every unit-ideal instance; the variant with an additive constant gives ⌊(k+14)/4⌋ + c. BF(3), BF(4) hold; BF(5) holds on three instances (V_6(13) full). BF(6) needs 1 ∈ V_7(17): not computed. | **proof** (Section 1) / **exact** (existing records) |
| M1 | **Min law.** For unequal depths (x₁ ∈ μ_{2^k₁}, x₂ ∈ μ_{2^k₂}), 1 ∈ V_D(k₁, k₂) iff min(k₁, k₂) ≤ 4D − 11, on all computed pairs: D = 3 (k₁ ≤ 8, k₂ ≤ 14), D = 4 (k₁ ≤ 8, k₂ ≤ 15), D = 5 (5 ≤ k₁ ≤ 9 up to (9, 11); (10, 10), (10, 11) fall), two unit instances (the second for D = 3, 4). The sum law k₁ + k₂ ≤ 8D − 22 is **refuted** (1 ∈ V_4(5, 15), sum 20). | **exact** per pair; the law **empirical** |
| M2 | **Coordinatewise free extension.** For min(k₁, k₂) ≥ 4D − 10 the dimension of V_D(k₁, k₂) depends only on s = k₁ + k₂: s + 1 (D = 3, 70 pairs), s² + 86s − 571 (D = 4, all pairs with min ≥ 6); one more bit in either coordinate multiplies the profile by (1 + z). | **exact** per pair; the law **empirical** |
| M3 | Consequence for the shape of a block lemma: under M1 the lemma can be stated in **one** coordinate, with the other one arbitrarily deep: BF1(D): 1 ∈ V_D(k, K) ⇒ 1 ∈ V_{D+1}(k+4, K) (k = 4D−11, K ≥ k+4). Instances D = 3 (K ≤ 15) and D = 4 (K = 9, 10, 11) hold. | **proof** of the equivalence under M1; instances **exact** |
| A1 | Candidate (a), doubling: Ψ_{Q,r}(w) = Q·w(x₁^{2^r}, x₂^{2^r}) maps V_D(k) into V_{D+s}(k+r) iff Q·g^{(r)} ∈ V_{D+s}(k+r) (wt Q ≤ s). | **proof** (Section 3.1) |
| A2 | On the data, g^{(r)} ∉ V_D(k) for every r ≥ 1 and k ≥ 4D − 10 (D = 3, 4, 5, 6): the pure doubling lemma is equivalent to the conclusion it is meant to prove. | **proof** from exact facts (Section 3.2) |
| A3 | **Refutation of doubling as an all-depth lemma.** No weight-≤1 multiplier works beyond the full depth: L_r(D,k) = {Q : Q g^{(r)} ∈ V_D(k)} is 0 for D = 3 (k = 2..5) and D = 4 (k = 6..9), r ≤ 5. In general, for k + r ≥ k_{D+s} (N8 range) every Ψ_{Q,r} fails unless Q·S_3(x^{2^r}, y^{2^r}, t) vanishes on the curve. | **exact** (Test A) / **proof** conditional on N8 T1 (Section 3.3) |
| A4 | Curve-compatible variant (norm-Frobenius): with h_r = Π_{ζ ∈ μ_{2^r}} g(ζx₁, x₂) (a polynomial in x₁^{2^r}), LFD_g(k₁+r, k₂) ≤ max(e_r, LFD_{h_r}(k₁, k₂)), e_r = entry degree of the norm in V[g]. On the data wt(h_r) = r+2, e₁ = e₂ = 4, and LFD_{h_r}(k₁, 6) = max(r+3, LFD_g(k₁+r, 6)) (r ≤ 3, k₁ ≤ 6): the route is lossless but trades one bit for one degree of generator weight. | **proof** (reduction) / **exact** (values) / **empirical** (the max-formula) |
| C1 | Candidate (c), top-block selectors in one coordinate: LFD_g(k₁+r, k₂) ≤ r + max_s LFD_{g_s}(k₁, k₂) over the 2^r coset instances g_s(y, x₂) = g(sy, x₂). The coset instances obey the same thresholds as g (Test C), so the route costs r degrees per r bits; the selector of a coset of index 2^r has weight exactly r, so no idempotent localization can do better. | **proof** (Section 4) / **exact** (Test C) |
| B1 | Candidate (b), lift charging one degree per carry block: stated precisely as BF1 composed with free deepening of the other coordinate (M1); no proof; the known obstruction (U1 T2′: low-degree parts do not lift) is restated. | **conjecture** (Section 5) |

**The exact gap.** LFD ≤ ⌊(k+14)/4⌋ + O(1) for unit ideals is equivalent to U_D (1 ∈ V_D(4D−11)) for
all D, up to a constant shift; proved is U_D for D ≤ 6 (exact, D = 6 on three instances) and the linear
bound LFD ≤ k+2 (kernel). Nothing here proves a single new U_D. For proper ideals the bound needs in
addition V_D = I_D (U1, Proposition 4), which no statement here addresses.

---

## 1. Setting and the target

A_{k₁,k₂} = F_p[Z/2^{k₁} × Z/2^{k₂}], monomial m_a, a = (a₁, a₂), weight wt(a) = wt(a₁) + wt(a₂) (binary
digit sums); u_{i,j} acts as a ↦ a + 2^j e_i (mod 2^{k_i}). V_D(k₁, k₂) is the least subspace of M_≤D
containing g = NF(S_3(x₁, x₂, t)) and closed under the u_{i,j} on elements of weight < D;
V_D(k) := V_D(k, k). By the kernel bridge (`ChainClosureBridge.chain_mem_closure_iff`) and Lemma 1.1 of
U1, for a unit ideal LFD(k) ≤ D ⟺ 1 ∈ V_D(k) (D ≥ 2). For unequal depths the bridge is written (U1,
Section 1) but not formalized; Section 2 uses unequal depths only in the model.

**Proposition 1.1 (proof).** Let D₀ ≥ 2 and c ≥ 0. Suppose 1 ∈ V_{D₀}(k₀) and, for every D ≥ D₀ and
every k, 1 ∈ V_D(k) ⇒ 1 ∈ V_{D+1}(k+4). Then 1 ∈ V_D(k) whenever k ≤ k₀ + 4(D − D₀), so for unit
ideals LFD(k) ≤ max(D₀, D₀ + ⌈(k − k₀)/4⌉). With k₀ = 4D₀ − 11 this is ⌊(k+14)/4⌋.

*Proof.* Induction on D gives 1 ∈ V_D(k₀ + 4(D − D₀)); depth restriction (U1 P3, **kernel**
`DepthRestrictionM2.not_one_of_le` in contrapositive form) gives every smaller depth. ∎

The hypothesis need only be checked at k = 4D − 11 (by depth restriction, 1 ∈ V_D(k) for k ≤ 4D−11
follows from the top depth, and the conclusion at k + 4 then follows from the one at 4D−7). Call that
instance **BF(D)**. It is (Step_D) of U1 Section 5.3. Status: BF(3) (1 ∈ V_3(1), 1 ∈ V_4(5)), BF(4)
(1 ∈ V_5(9)) and BF(5) (1 ∈ V_6(13), unit instances over 65537 and 163841) are **exact** from existing records
(`experiments/ic-bounds-k14/RESULTS.md`). BF(6) is the statement 1 ∈ V_7(17) (|M_≤7| at k = 17 is about
5.4·10⁶ columns): not computed.

## 2. Test M: unequal depths (the block structure is per coordinate)

`scripts/ic_bounds_block4_mixed.py OUT D k1min k1max k2max [maxcols] [unit|unit2]`. Instances: the unit
instance of `ic-bounds-k14` (y² = x³ + x over F_65537, t = 4271) and the second unit instance
(y² = x³ + 3x over F_163841, t = 73730). The solution count on μ_{2^k₁} × μ_{2^k₂} is recomputed for every
pair (field `solutions`, 0 everywhere: unit ideals).

**Engine checks.** The fast engine (block RREF of `csd/ic_bounds_d5_verify.Part`, unequal-depth layout)
reproduces the slow engine of `scripts/ic_bounds_growth_lib.py` (layout generalized) on all 84 D = 3 pairs
(`--slow-check`, `engines_agree: true`) and on the D = 4 pairs (4, 9), (5, 5), (3, 6); at equal depths it
reproduces dim V_4(6) = 605, V_4(7) = 829, V_4(8) = 1061 and the full spaces V_5(k) = M_≤5, 5 ≤ k ≤ 9.

**Results (exact per pair).**

| D | 1 ∈ V_D(k₁, k₂) for | 1 ∉ V_D(k₁, k₂) for | instances |
|---|---|---|---|
| 3 | k₁ = 1, k₂ = 1..14 | 2 ≤ k₁ ≤ 8, k₁ ≤ k₂ ≤ 14 | 65537 (84 pairs); 163841 (k₁ ≤ 4, k₂ ≤ 12) |
| 4 | k₁ ≤ 5, k₁ ≤ k₂ ≤ 15 | (6, 6..15), (7, 7..14), (8, 8..13) | 65537 (89 pairs); 163841 (4 ≤ k₁ ≤ 7, k₂ ≤ 12) |
| 5 | 5 ≤ k₁ ≤ 9 (pairs with ≤ 9000 columns, and (9, 9), (9, 10), (9, 11)) | (10, 10), (10, 11) | 65537 |

D = 5 rows computed: all pairs 5 ≤ k₁ ≤ 8, k₁ ≤ k₂ with |M_≤5| ≤ 9000 (full, 1 ∈ V_5), and the pairs
listed in Section 2.1 under the heavy lock.

**M1 (empirical).** 1 ∈ V_D(k₁, k₂) ⟺ min(k₁, k₂) ≤ 4D − 11. The deeper coordinate is free: 1 ∈ V_4(5, 15)
and 1 ∈ V_4(1, 15), while 1 ∉ V_4(6, 6). The **sum law** "1 ∈ V_D iff k₁ + k₂ ≤ 8D − 22" (which would
have made "four bits" mean two bits per coordinate) is **refuted** by (5, 15) for D = 4 and by (1, 14) for
D = 3 (exact). The law is consistent with U1 Theorem 1, LFD ≤ min(k₁, k₂) + 2, which is sharp at
min = 1, 2.

**M2 (empirical).** For min(k₁, k₂) ≥ 4D − 10, dim V_D(k₁, k₂) depends only on s = k₁ + k₂: s + 1 for D = 3
(70 pairs, 0 mismatches), s² + 86s − 571 for D = 4 (all pairs with min ≥ 6, 0 mismatches; at s = 2k this is
the known 4k² + 172k − 571). Profiles: (6, 6) → (6, 7) is (0,0,2,111,605) → (0,0,2,113,716), i.e.
p_j ↦ p_j + p_{j−1}: one bit in one coordinate multiplies the profile polynomial by (1 + z), so the known
recursion p_j(k+1) = p_j + 2p_{j−1} + p_{j−2} is two such steps. Both instances give the same dimensions and
profiles.

**Proposition 2.1 (proof).** Assume M1 for the instance. Then BF(D) holds iff 1 ∈ V_{D+1}(4D − 7, K) for
one (equivalently every) K ≥ 4D − 7, and BF(D) is implied by the one-coordinate statement

    BF1(D):  1 ∈ V_D(k, K)  ⇒  1 ∈ V_{D+1}(k + 4, K)      (k = 4D − 11, all K ≥ k + 4).

*Proof.* Under M1 membership of 1 depends only on min(k₁, k₂). ∎ (Without M1 the implication
BF1 ⇒ BF still holds when K = k + 4, which is BF itself.)

So the "block of four" is four bits of **one** coordinate per degree, with the other coordinate as deep as
one likes. Any proof has to explain why the second coordinate costs nothing. The only mechanism in the
repository with that property is the univariate fibre lemma (U1 Lemma 2.3: a quadratic fibre in x₂ reaches
1 with budget 2 at every depth of x₂); it is used in U1 together with a full selector in x₁ (weight k₁).
BF1 asks for a replacement of that selector costing one degree per four bits of x₁.

### 2.1 D = 5 near the threshold (heavy lock)

Unit instance over 65537, `ic_bounds_block4_mixed.py OUT 5 9 9 11 30000` and `... 5 10 10 11 30000` under
`flock /tmp/lfd-locks/heavy.lock` (`ulimit -v 7000000`, two threads). Every row closed (full sweep clean),
0 solutions.

| (k₁, k₂) | columns | dim V_5 | 1 ∈ V_5 | seconds |
|---|---|---|---|---|
| (9, 9) | 12,616 | 12,616 | yes | 219 |
| (9, 10) | 16,664 | 16,664 | yes | 437 |
| (9, 11) | 21,700 | 21,700 | yes | 856 |
| (10, 10) | 21,700 | 18,569 | no | 923 |

(10, 10) reproduces the certified dim 18,569 and profile (0, 0, 2, 138, 3163, 18569) of V_5(10) (third
engine). Prediction frozen here before the run finished, from M2 (one bit in one coordinate multiplies the
profile by 1 + z): V_5(10, 11) has profile (0, 0, 2, 140, 3301, 21732) and does not contain 1. Result:
| (10, 11) | 27,896 | 21,732 | no | 1,759 |

profile (0, 0, 2, 140, 3301, 21732), equal to the prediction (**exact**; one confirmed prediction of M2 at D = 5).

So 1 ∈ V_5(9, 11) although 9 + 11 = 20 > 18 = 2·9: the min law holds at the D = 5 threshold on the
computed pairs, and BF1(4) holds for K = 9, 10, 11 (1 ∈ V_4(5, K) from Test M and 1 ∈ V_5(9, K) here).

## 3. Candidate (a): doubling / Frobenius maps

### 3.1 The map lemma (proof)

For r ≥ 1 let φ_r : A_{k} → A_{k+r}, m_a ↦ m_{2^r a}. It is injective, preserves weight, and intertwines
u_{i,j} with u_{i,j+r} (2^r(a + 2^j) = 2^r a + 2^{j+r}, and reduction mod 2^k corresponds to reduction
mod 2^{k+r}). For Q ∈ A_{k+r} of weight ≤ s put Ψ_{Q,r}(w) = Q·φ_r(w); by (1.1) of U1 it raises weight by
at most s and it intertwines u_{i,j} with u_{i,j+r} (multiplication commutes). Write g^{(r)} = φ_r(g),
the generator with exponents multiplied by 2^r (the lift of S_3(x₁^{2^r}, x₂^{2^r}, t)).

**Lemma 3.1 (proof).** For D ≥ 2: Ψ_{Q,r}(V_D(k)) ⊆ V_{D+s}(k+r) ⟺ Q·g^{(r)} ∈ V_{D+s}(k+r). In that case
1 ∈ V_D(k) ⇒ Q ∈ V_{D+s}(k+r). If moreover s = 0, Q = 1, or s = 1 and Q = x₁ − β with β^{2^{k+r}} ≠ 1,
then 1 ∈ V_D(k) ⇒ 1 ∈ V_{D+s}(k+r).

*Proof.* (⇒) g ∈ V_D(k). (⇐) {w : Ψ(w) ∈ V_{D+s}(k+r)} is a subspace containing g; if w is in it and
wt(w) < D then wt Ψ(w) < D + s and Ψ(u_{i,j} w) = u_{i,j+r} Ψ(w) ∈ V_{D+s}(k+r); so it contains V_D(k).
For Q = x₁ − β: U1 Lemma 2.3 (degree-1 case, budget 2 ≤ D + 1) applied to x₁ − β gives 1 − β^{2^{k+r}}
in the closure. ∎

So BF(D) would follow from **A(D): g^{(4)} ∈ V_{D+1}(4D − 7)**, or from Q·g^{(4)} ∈ V_{D+1}(4D−7) for a
weight-1 Q as above.

### 3.2 Forcing: the doubling lemma at fixed degree is the conclusion itself (proof from exact facts)

**Proposition 3.2.** Let D be such that 1 ∈ V_D(k) for k ≤ 4D − 11 and 1 ∉ V_D(4D − 10) (true for
D = 3, 4, 5, 6 on the unit instance over 65537: **exact**, existing records and Test M). Then
g^{(r)} ∉ V_D(k) for every r ≥ 1 and every k ≥ 4D − 10.

*Proof.* At k = 4D − 10 and 1 ≤ r ≤ k − 1: if g^{(r)} ∈ V_D(k), Lemma 3.1 (s = 0) and 1 ∈ V_D(k − r) give
1 ∈ V_D(k), contradiction. For r ≥ k, g^{(r)} = g(1, 1)·1 with g(1, 1) ≠ 0 (unit ideal), not in V_D(k).
For larger k use the depth restriction π (U1 P3): π(V_D(k+1)) ⊆ V_D(k) and π(g^{(r)}) = g^{(r)}. ∎

Hence, on the data, A(D) holds exactly when 1 ∈ V_{D+1}(4D−7), and is no easier. The same argument with
weight-1 multipliers needs 1 ∉ V_{D+1}(k+r) and 1 ∈ V_D(k), which happens at k + r = 4D − 6, so
Q·g^{(4)} ∈ V_{D+1}(4D − 6) is not forced false by logic alone; Test A refutes it on the data where
computable (D + 1 = 4: Q·g^{(4)} ∉ V_4(6..9) for every weight-1 Q ≠ 0).

### 3.3 Obstruction for all-depth map lemmas (proof, conditional on N8 T1)

**Theorem 3.3.** Assume the hypotheses of N8 Theorem T1 (char K = 0 or p > 2^{k+r+1}; E nonsingular;
t ≠ 0; C the curve S_3(x₁, x₂, t) = 0), D + s ≥ 2 and k + r ≥ k_{D+s} (the N8 range, k_D ≤ 14D + 10).
Let Q have weight ≤ s and exponents < 2^{k+r} − 2^{r+1}. If Q·S_3(x₁^{2^r}, x₂^{2^r}, t) does not vanish
identically on C, then Ψ_{Q,r}(V_D(k)) ⊄ V_{D+s}(k+r).

*Proof.* Otherwise Ψ(g) = Q·g^{(r)} ∈ V_{D+s}(k+r) by Lemma 3.1. Its exponents are < 2^{k+r}, so it is read
as the polynomial Q·S_3(x₁^{2^r}, x₂^{2^r}, t), and N8 T1 says every element of V_{D+s}(k+r) so read
vanishes on C. ∎

So a depth-raising **linear map** lemma that is to hold at all depths must send g into the curve ideal
(vanish on C). Doubling does not: on the unit instance S_3(x₁^{2^r}, x₂^{2^r}, t) is nonzero at all 200
sampled points of C(F_p) for r = 1..6, and no nonzero weight-≤1 Q with exponents < 2^15 vanishes at all
of them (evaluation rank 31 of 31) (`curve.jsonl`, **exact** for the sampled points; with exponents < 2^16
the rank is 32 of 33, the one relation being x₁^{2^15} + x₂^{2^15} = 0 on C(F_p), a Legendre-symbol relation
on F_p-points that encodes the unit property of this instance, not a relation on C).
Since BF's hypothesis 1 ∈ V_D(k) holds only for k ≤ 4D − 11, a proof of BF by a map has to use the depth;
uniform maps of the doubling type are excluded.

### 3.4 Test A (exact)

`ic_bounds_block4_maps.py doubling OUT D k rmax` (records `doubling.jsonl`): membership of g^{(r)} and
dim L_r with L_r = {Q ∈ span(1, u_{i,j}) : Q g^{(r)} ∈ V_D(k)}, and membership of (x₁ − 3)g^{(r)}
(3 has order 2^16).

| D | k | dim V_D(k) | 1 ∈ V | g^{(r)} ∈ V (r = 1..5) | dim L_r |
|---|---|---|---|---|---|
| 3 | 2, 3, 4, 5 | 5, 7, 9, 11 | no | no | 0 |
| 4 | 4, 5 | 163, 386 | yes | yes | 2k+1 (all) |
| 4 | 6, 7, 8, 9 | 605, 829, 1061, 1301 | no | no | 0 |

(For 2^r ≥ 2^k the element g^{(r)} wraps to lower weight, field `g_r_weight`; the conclusion is the same.)
**Counterexample to the all-depth doubling lemma** (D = 3 → 4, r = 4): Q·g^{(4)} ∉ V_4(k+4) for every
weight-≤1 Q ≠ 0 at k + 4 = 6, 7, 8, 9, although Lemma 3.1 would need it there.

### 3.5 The curve-compatible variant: norm-Frobenius (proof + exact)

Let h_r(X, x₂) = Π_{ζ ∈ μ_{2^r}} g(ζx₁, x₂). It is invariant under x₁ ↦ ζx₁, hence a polynomial in
X = x₁^{2^r}, of X-degree ≤ 2 and x₂-degree ≤ 2^{r+1}; its X-exponents have weight ≤ 1 and its x₂-exponents
weight ≤ r + 1, so wt(h_r) ≤ r + 2. Let φ : A_{k₁,k₂} → A_{k₁+r,k₂}, X^a x₂^b ↦ x₁^{2^r a} x₂^b (weight
preserving, intertwines u_{1,j} ↦ u_{1,j+r}, u_{2,j} ↦ u_{2,j}); φ(h_r) = N_r := Π_ζ g(ζx₁, x₂) ∈ gA.

**Proposition 3.4 (proof).** Let e_r = min{E : N_r ∈ V_E[g](k₁+r, k₂)} (≤ 2^{r+1} by budgeted products).
For D ≥ e_r, φ(V_D[h_r](k₁, k₂)) ⊆ V_D[g](k₁+r, k₂). h_r is a unit at (k₁, k₂) iff g is a unit at
(k₁+r, k₂) (h_r(α^{2^r}, β) = Π_ζ g(ζα, β)). Hence, for unit ideals,
LFD_g(k₁+r, k₂) ≤ max(e_r, LFD_{h_r}(k₁, k₂)).

*Proof.* As Lemma 3.1, with g replaced by h_r on the left and the generator g on the right; the input
condition is φ(h_r) = N_r ∈ V_D[g]. ∎

Test N (`norm.jsonl`, `normentry.jsonl`, **exact**): wt(h_r) = 3, 4, 5, 6 for r = 1, 2, 3, 4; e₁ = 4 and
e₂ = 4 at (6, 6), (7, 7), (6, 8) (e₁ ≥ 4: not in V_3 at (4, 4), (6, 6)); LFD_{h_r}(k₁, 6) for k₁ = 1..6:

| r | wt h_r | 1 ∈ V_D[h_r](k₁, 6) | equals max(r+3, LFD_g(k₁+r, 6))? |
|---|---|---|---|
| 1 | 3 | D = 4: k₁ ≤ 4; D = 5: all k₁ ≤ 6 | yes (LFD_g(k₁+1, 6) = 4 for k₁ ≤ 4, 5 for k₁ = 5, 6) |
| 2 | 4 | D = 4: none; D = 5: all k₁ ≤ 6 | yes |
| 3 | 5 | D = 5: none; D = 6: all k₁ ≤ 6 | yes |
| 4 | 6 | D = 6: none (D = 7 not run) | consistent |

**Conjecture N (empirical on these cells).** LFD_{h_r}(k₁, k₂) = max(wt(h_r) + 1, LFD_g(k₁ + r, k₂)).
The norm route is then lossless and circular: it moves r bits of x₁ into the generator (weight r + 2, x₂-degree
2^{r+1}), one degree of generator weight per bit. For r = 4 it gives LFD_g(k₁+4, k₂) ≤ max(e₄, LFD_{h₄}(k₁, k₂))
with LFD_{h₄} ≥ 7 (a closure at D = wt(h₄) is span(h₄)), no gain for D ≤ 6 and no mechanism for one degree
per four bits.

## 4. Candidate (c): localization on the top block of x₁

For γ ∈ μ_{16} the selector of the coset {x₁^{2^{k₁}} = γ} at x₁-depth k₁ + r (r = 4 for the block) is

    e_γ = 2^{−r} Σ_{c < 2^r} γ^{−c} x₁^{c·2^{k₁}} = Π_{j=k₁}^{k₁+r−1} (1 + γ_j^{−1} u_{1,j}) / 2.

**Lemma 4.1 (proof).** e_γ is the unique idempotent supported on the coset; its weight is exactly r (the
monomial x₁^{(2^r − 1)2^{k₁}} has a nonzero coefficient). Hence any localization of x₁ to a coset of index
2^r by an idempotent costs r degrees.

**Proposition 4.2 (proof).** For s with s^{2^{k₁}} = γ put g_s(y, x₂) = g(sy, x₂) (a coset instance at
depths (k₁, k₂)). For D′ ≥ 2: if 1 ∈ V_{D′}[g_s](k₁, k₂) for all 2^r cosets, then 1 ∈ V_{D′+r}[g](k₁+r, k₂).
So LFD_g(k₁+r, k₂) ≤ r + max_s LFD_{g_s}(k₁, k₂) (unit ideals).

*Proof.* ψ(y^a x₂^b) = s^{−a} x₁^a x₂^b e_γ (a < 2^{k₁}) intertwines u_{1,j} (j < k₁) with
s^{2^j}·u_{1,j}, including the wrap (x₁^{2^{k₁}} e_γ = γ e_γ = s^{2^{k₁}} e_γ), and u_{2,j} with u_{2,j};
it raises weight by at most r, and ψ(g_s) = g·e_γ ∈ V_{2+r}[g] (budgeted product). As in Lemma 3.1,
ψ(V_{D′}[g_s]) ⊆ V_{D′+r}[g], so e_γ ∈ V_{D′+r}[g] for every γ, and Σ_γ e_γ = 1. ∎ (This is U1 Lemma 2.4
with a partial selector; U1 Section 5.3 item 1 records the two-coordinate version, +2 per level.)

Test C (`coset.jsonl`, **exact**, all 2^r cosets, every coset instance a unit ideal):

| D | (k₁, k₂) | r | 1 ∈ V_D[g_s] for all s | reading |
|---|---|---|---|---|
| 3 | (1, 6), (1, 8) | 4 | yes | bound LFD_g(5, 6) ≤ 7, LFD_g(5, 8) ≤ 7; true value 4 |
| 3 | (2, 6) | 4 | no | |
| 4 | (2, 6) | 4 | yes | bound LFD_g(6, 6) ≤ 8; true value 5 |
| 4 | (5, 9) | 4 | yes | bound LFD_g(9, 9) ≤ 8; true value 5 |
| 4 | (5, 7) | 1 | yes | |
| 4 | (6, 7), (6, 8) | 1, 4 | no | coset instances fall at x₁-depth 6, like g |

The coset instances satisfy the same min law as g (dimensions equal to those of g at the same depths:
794, 716, 1471, 829). So the route gives r + f(k₁) where the truth is f(k₁ + r) ≈ f(k₁) + r/4: it costs r
per r bits whatever r is, i.e. slope 1. **Candidate (c) with weight-4 selectors is refuted as a route to
"+1 per four bits"**: by Lemma 4.1 the selector alone costs 4, and the fibre problems are not easier.

## 5. Candidate (b): lifting a derivation with one degree for the carry block

Precise form. Write the depth-(k+4) exponent of x₁ as a₁ = l + 2^k h, h ∈ [0, 16) (carry block). A
derivation of 1 at depth k with budget D, replayed at depth k + 4 without the wrap at bit k, produces an
element whose components at h ≠ 0 record the carries. The candidate is that one extra degree suffices to
cancel them, uniformly in the other coordinate:

    (B)  1 ∈ V_D(k, K)  ⇒  1 ∈ V_{D+1}(k + 4, K)      (k = 4D − 11, every K ≥ k + 4),

which is BF1(D). Status:

* **exact instances**: D = 3 (1 ∈ V_3(1, K) and 1 ∈ V_4(5, K), K ≤ 15), D = 4 (1 ∈ V_4(5, K) and
  1 ∈ V_5(9, K) for K = 9, 10, 11; Section 2.1);
* **sharpness, exact**: 1 ∉ V_4(6, K) (K ≤ 15), 1 ∉ V_3(2, K): the step cannot be five bits;
* **no proof**. The known obstruction stands (U1 T2′): at the transition, elements of V_D(k) of low degree
  lift to V_D(k+1) only in full degree (the least degree of a lift of 1 from depth 5 to depth 6 is 4 = D),
  so the extra degree cannot be charged to a low-degree part of the derivation; it has to be charged to the
  whole space, as the coordinatewise free extension (M2) does above the base.
* The free deepening of the second coordinate (M1) is part of (B) and is itself unproved: the inclusion
  1 ∈ V_D(k, k) ⇒ 1 ∈ V_D(k, k + 1) has no proof here (depth restriction goes the other way; a top-bit
  selector in x₂ costs one degree; the norm route in x₂ changes the generator).

## 6. Summary: proved, tested, refuted, open

* **Proved (written).** Proposition 1.1 (reduction of the slope to BF), Proposition 2.1 (one-coordinate
  form under M1), Lemma 3.1 and Proposition 3.2 (doubling is equivalent to the conclusion on the data),
  Theorem 3.3 (all-depth map lemmas must preserve the curve ideal; conditional on N8 T1), Proposition 3.4
  (norm-Frobenius reduction), Lemma 4.1 and Proposition 4.2 (selector route, cost r per r bits).
* **Tested (exact).** Test M (min law and coordinatewise free extension, two unit instances, D = 3, 4, and
  D = 5 near k₁ = 9, 10), Test A (doubling: refuted beyond the full depth), Test C (coset instances obey the
  same thresholds), Test N (norm route lossless on the computed cells).
* **Refuted.** The all-depth doubling lemma (Section 3.4, counterexample D = 3 → 4 at depths 6..9, every
  weight-≤1 multiplier); the selector route as a +1-per-four-bits mechanism (Section 4); the sum law for
  unequal depths (Section 2).
* **Open.** BF(D) for D ≥ 6 (first case 1 ∈ V_7(17)); a proof of the min law M1 (even of
  1 ∈ V_D(k, k) ⇒ 1 ∈ V_D(k, k+1)); a mechanism charging one degree per four bits of one coordinate.
  Concrete next step: a selector-free proof of 1 ∈ V_4(5, K) for all K (the smallest nontrivial instance of
  BF1 with the second coordinate free; the selector proof gives only 1 ∈ V_7(5, K)), then the same for
  1 ∈ V_5(9, K).
