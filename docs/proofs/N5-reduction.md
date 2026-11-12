# N5: the depth reduction for the degree-4 closure, stated in the ordinary ring

Status (2026-11-08). Written proof, self-contained, of the reduction theorem N5
(`docs/IC_BOUNDS_PROGRAM.md`, Section 4) for the closure `InClosure` of
`formal/ICFallDegrees/Closure.lean`, that is, the ordinary iterated closure in the
polynomial ring, not only the chain-quotient group algebra. The parts that are kernel-checked
are named in Section 9 with their Lean theorems; the parts that remain written-only are named
there too. No `sorry`, `axiom` or `native_decide` is used; every theorem listed depends on
`propext`, `Classical.choice` and `Quot.sound` only.

The one-paragraph summary:

1. **Bridge (proved here, formalized).** For `F = C ∪ R` with `C` the chain equations and
   `R` any further inputs of degree `≤ D`, `D ≥ 2`, ordinary closure membership is
   `InClosure F D p ⟺ deg p ≤ D ∧ ρ(p) ∈ V_D^q`, where `ρ` is the chain normal form and
   `V_D^q` is the closure computed in the reduced (squarefree) space with the rule
   `w ↦ ρ(u_{i,j} w)` on elements of degree `< D`. In particular `1 ∈ V_D(F) ⟺ 1 ∈ V_D^q`,
   and a reduced containing space closed under `ρ(u_{i,j}·)` and excluding `1` refutes
   `1 ∈ V_D(F)` in the ordinary ring. This replaces the previously unproved premise
   "closure modulo the chains equals the ordinary closure for `D ≥ 2`".
2. **Model identification (proved here, formalized).** On squarefree monomials
   `ρ(u_{i,j} · m_a) = m_{a + 2^j e_i mod 2^k}` (binary addition with carries, wrap at the top
   level), so `V_D^q` is the group-algebra closure that `csd/ic_bounds_fall.py` computes
   (`ChainTranslation.normalForm_X_mul_monomial`).
3. **Reduction theorem (proved here, formalized in a two-level exponent model).** For any
   subspace `S` of the depth-`k0` algebra inside `M_{≤D}`, the free high extension
   `W_k(S) = ⊕_h τ_h (S ∩ M_{≤ D − wt h})` satisfies `1 ∈ W_k ⟺ 1 ∈ S`, and if `W_{k0+1}` is
   closed under all admissible variable products at depth `k0+1` then `W_k` is closed at every
   depth `k ≥ k0+1`. Nothing about `S` beyond `S ⊆ M_{≤D}` is used. With `k0 = 7`, `D = 4`
   this is N5 (iii) and (ii); N5 (i) (generator rows) is also proved but is not needed for
   the lower bound.
4. **Consequence (formalized, `ChainModelMap.n5_lfd_ge_five`).** On an instance whose base
   `S = V_4(7)` excludes `1`, contains the generator, and whose `W_8` passes the closedness
   certificate, `1 ∉ V_4(k)` in the ordinary ring for every `k ≥ 8`; a unit-ideal instance
   therefore has `LFD ≥ 5` at every depth `k ≥ 8`. The three base facts are the only
   hypotheses besides `2 ≠ 0` and `1 ∈ ⟨F_k⟩`. The dimension statement
   `dim V_4^q(k) ≤ dim W_k = 4k^2 + 172k − 571` holds for the quotient closure and transfers
   to the ordinary closure as `dim V_4(k) = dim(C ∩ M_{≤4}) + dim V_4^q(k)` (written only).

The written argument of `experiments/ic-bounds-fall/README.md` Section 4.2 is correct as a
statement about the group algebra; the gaps relative to a statement in the ordinary ring, and
two places where the earlier text is instance-specific or imprecise, are listed in Section 8.

---

## 1. Setting: the ordinary ring and the closure

Fix a field `K` with `2 ≠ 0` in `K` (in the application `K = F_p`, `p ≥ 5`), a number of
summands `m ≥ 1` (the application has `m = 2`) and a depth `k ≥ 2`. The ordinary ring is

    R_k = K[u_{i,j} : i ∈ {0,…,m−1}, j ∈ {0,…,k−1}],

`m k` variables; in Lean `MvPolynomial (Fin (m*k)) K`, with `u_{i,j} = X (var m k i j)`,
`var m k i j = j + k·i` (`ChainQuotient.var`, `row`, `depth`). The **chain equations** are

    c_{i,j} = u_{i,j}^2 − u_{i,j+1}  (j < k−1),      c_{i,k−1} = u_{i,k−1}^2 − 1,

the set `C = C_k = {c_{i,j}}` (`ChainQuotient.chainInputs K m k`); `⟨C⟩` denotes the ideal they
generate. For `m = 2` the **normalized generator** is

    g = NF(S_3(u_{0,0}, u_{1,0}, t)) = u_{0,1} u_{1,1} − 2(t^2+a) u_{0,0} u_{1,0}
        − 2t (u_{0,0} u_{1,1} + u_{0,1} u_{1,0}) + t^2 (u_{0,1} + u_{1,1})
        − (2at + 4b)(u_{0,0} + u_{1,0}) + (a^2 − 4bt)

(`PrimeM2General.normalized k depth a b t`; the raw `S_3` differs from `g` by an explicit
multiple of two chain equations, `raw_normalization_identity`). The input list is
`F_k = C_k ∪ {g}` (`Set.range (PrimeM2General.inputs k depth a b t)`; the lemma
`ChainClosureBridge.inputs_castSucc_mem_chainInputs` identifies the first `2k` entries with
the chain equations, and `ChainQuotient.chainInputs_two_subset` the converse).

**The closure.** For `F ⊆ R_k` and `D ∈ ℕ`, `V_D(F) = InClosure F D` is the least
`K`-subspace of `R_k` that contains every `f ∈ F` with `deg f ≤ D` and is closed under
`p ↦ a·p` whenever `deg a + deg p ≤ D` (`Closure.lean`; the multiplier rule applies to derived
rows). Its least-space property is `closure_minimal`; by `VariableClosure.closure_minimal_by_
variables` it suffices to check products by single variables `X v` on elements of degree
`< D`. `HasLFDBound F d` means `∀ f ∈ ⟨F⟩, f ∈ V_{max(d, deg f)}(F)`; for a unit ideal,
`HasLFDBound F e → 1 ∈ V_e(F)` (`one_of_hasLFDBound`), so `1 ∉ V_D(F)` gives
`∀ e ≤ D, ¬HasLFDBound F e`, that is `LFD(F) ≥ D+1`.

Degrees: `deg` is total degree; every element of `V_D(F)` has degree `≤ D`
(`InClosure.degree_le`), and `V_D(F) ⊆ ⟨F⟩` (`InClosure.sound`).

## 2. The chain normal form and the quotient model

### 2.1 Reduced space and normal form

A monomial `u^e` (`e : Fin (mk) →₀ ℕ`) is **standard** for `C` if no leading monomial
`u_{i,j}^2` of a chain equation divides it, that is, iff `e ≤ 1` coordinatewise: the chain
equations are squarefree-leading (`ChainQuotient.leading_chainPolynomial`:
`degLex.degree c_v = single v 2`; `standard_le_one`). The **reduced space** `Red` is the span
of the squarefree monomials (`FiniteQuotient.reducedSpace (chainInputs K m k)`). Two facts
from `FiniteQuotient.lean` / `ChainQuotient.lean` are used:

* (Division with budget, `exists_reduced_remainder`.) For every `p` there is `z ∈ Red` with
  `p − z ∈ InClosure C (deg p)`; in particular `p − z ∈ ⟨C⟩` and `deg z ≤ deg p`
  (`remainder_degree_le`).
* (Reduced-ideal triviality, `reduced_chain_ideal_trivial`, needs `2 ≠ 0`.) If `p ∈ ⟨C⟩`
  and `p ∈ Red` then `p = 0`. Hence the `z` above is unique.

Define `ρ(p)` to be this unique `z` (`ChainClosureBridge.normalForm`). Then
(`normalForm_*` in `ChainClosureBridge.lean`):

    ρ is K-linear;  ρ(z) = z for z ∈ Red;  ρ(c) = 0 for c ∈ ⟨C⟩;  deg ρ(p) ≤ deg p;
    p − ρ(p) ∈ InClosure C (deg p);   ρ(a·p) = ρ(a·ρ(p));   ρ(1) = 1.

The last holds because `1 = u^0` is squarefree (`chain_one_reduced`). The identity
`ρ(a p) = ρ(a ρ(p))` is uniqueness applied to `a p − ρ(a ρ(p)) = a(p − ρ(p)) + (a ρ(p) − ρ(a ρ(p)))
∈ ⟨C⟩`.

### 2.2 Squarefree monomials as exponent vectors; the translation identity

For a squarefree exponent `e` put `bitExp(e)_i = Σ_j 2^j e(var i j) ∈ [0, 2^k)`
(`ChainQuotient.bitExp`); this is a bijection between squarefree exponents and exponent
vectors `a ∈ [0,2^k)^m` (`bitExp_injective`, `bitExp_lt`), and `deg u^e = Σ_i wt(a_i)` where
`wt` is the binary digit sum (popcount). Write `m_a` for the squarefree monomial with
`bitExp = a`, and `M_{≤d}` for the span of the `m_a` with `wt(a) := Σ_i wt(a_i) ≤ d`; so
`Red ∩ {deg ≤ d} = M_{≤d}`.

**Lemma 2.2 (translation identity).** For every `a ∈ [0,2^k)^m`, `i < m`, `j < k`,

    ρ(u_{i,j} · m_a) = m_{a'},   a' = a + 2^j e_i  (mod 2^k, coordinatewise).

*Proof.* By induction on `k − j ≥ 1`. Let `e` be the exponent of `m_a` and `v = var i j`.

*Bit `j` of `a_i` is `0`.* Then `e(v) = 0`, so `u_{i,j} m_a = u^{e + δ_v}` is squarefree, hence
reduced, with `bitExp = a + 2^j e_i`; and `a_i + 2^j < 2^k` because the bits of `a_i + 2^j`
are those of `a_i` plus bit `j`, all below `k`. So `ρ(u_{i,j} m_a) = u^{e+δ_v} = m_{a'}`.

*Bit `j` of `a_i` is `1`.* Then `m_a = u_{i,j} · m_b` with `b = a − 2^j e_i` (exponent
`e − δ_v`, squarefree), and `u_{i,j} m_a = u_{i,j}^2 m_b ≡ tail(u_{i,j}) · m_b (mod ⟨C⟩)`.
If `j + 1 < k`, `tail = u_{i,j+1}` and by the induction hypothesis (at `j+1`, with `b`)
`ρ(u_{i,j+1} m_b) = m_{b''}`, `b'' = b + 2^{j+1} e_i mod 2^k = a + 2^j e_i mod 2^k` (since
`−2^j + 2^{j+1} = 2^j`); congruent polynomials have equal normal forms, so
`ρ(u_{i,j} m_a) = m_{a'}`. If `j = k−1`, `tail = 1` and `ρ(u_{i,j} m_a) = ρ(m_b) = m_b`, and
`b_i = a_i − 2^{k−1} = (a_i + 2^{k−1}) mod 2^k` because bit `k−1` of `a_i` is set and
`a_i < 2^k`. ∎

Thus multiplication by `u_{i,j}` followed by reduction is the translation
`τ_{i,j} : m_a ↦ m_{a + 2^j e_i mod 2^k}` of the group algebra `A_k = K[(Z/2^k)^m]`, with the
identification `m_a ↔ δ_a`, and `deg m_a = wt(a)`. This is exactly the model of
`csd/ic_bounds_fall.py` (`make_grid`, `perms`) and of `README.md` Section 1. In Lean,
Lemma 2.2 is `ChainTranslation.translate_monomial_mem_span` (the congruence, by induction on
`k − 1 − j`, with the bit lemmas `digits_update_of_bit_zero` / `digits_of_bit_one` for the
two cases) and `ChainTranslation.normalForm_X_mul_monomial` (the normal-form statement).

### 2.3 The quotient closure

For `F ⊆ R_k` and `D`, the **quotient closure** `V_D^q(F) = InQuotClosure C F D` is the least
subspace of `Red` containing `ρ(f)` for every `f ∈ F` with `deg f ≤ D`, and closed under
`w ↦ ρ(X_v · w)` for its elements `w` of degree `< D` and all variables `v`
(`ChainClosureBridge.InQuotClosure`). For `F = C ∪ {g}` the inputs reduce to `ρ(c) = 0` and
`ρ(g) = g` (`g` is squarefree-supported), and by Lemma 2.2 the product rule is the
translation rule; so `V_D^q(F_k)` is the group-algebra closure `V_D` of `README.md` Section 1
("the least subspace of `M_{≤D}` containing `g` and closed under `τ_{i,j}` on elements of
weight `< D`"; the generator rows `τ_μ g`, `wt μ ≤ D−2`, are then automatic, see Section 6).

Elementary properties (`InQuotClosure.reduced`, `.degree_le`): `V_D^q(F) ⊆ Red ∩ M_{≤D}`.

## 3. The bridge: ordinary closure = chain part ⊕ quotient closure

**Theorem 3.1** (`ChainClosureBridge.mem_closure_iff`, `chain_mem_closure_iff`). Let `G` be
reduced-ideal trivial with `0 ∉ G`, `G ⊆ F`, and `deg g ≤ D` for all `g ∈ G`. Then for every
`p ∈ R`:

    InClosure F D p  ⟺  deg p ≤ D  ∧  InQuotClosure G F D (ρ p).

For `G = C_k` the degree hypothesis is `D ≥ 2`.

*Proof.* (⇒) `deg p ≤ D` is `InClosure.degree_le`. For `ρ(p) ∈ V_D^q` we prove more
generally (`closure_normalForm_mem`, which needs neither `G ⊆ F` nor the degree bound):
`InClosure F D p → InQuotClosure G F D (ρ p)`, by induction on the derivation of `p`.
Inputs: `ρ(f)` is a quotient input. Zero, sums, scalars: `ρ` is linear and `V_D^q` is a
subspace. Products `a·p` with `deg a + deg p ≤ D`: `ρ(a p) = ρ(a ρ(p))`, `ρ(p) ∈ V_D^q` by
induction and `deg a + deg ρ(p) ≤ D`; it remains to show

  (3.2) `w ∈ V_D^q`, `deg a + deg w ≤ D` ⟹ `ρ(a w) ∈ V_D^q`   (`quot_mul_mem`).

Expand `a = Σ_e coeff_e(a) u^e` with `deg e ≤ deg a`; by linearity of `ρ` it suffices to treat
monomials `u^e` with `deg e + deg w ≤ D`, by induction on `e` one variable at a time: for
`e = δ_v + e'`, `ρ(u^e w) = ρ(X_v · u^{e'} w) = ρ(X_v · ρ(u^{e'} w))`, where `ρ(u^{e'} w) ∈ V_D^q`
by induction and has degree `≤ deg e' + deg w < D`, so the quotient product rule applies.
(This is the same induction as `VariableClosure.variable_closed_mul_mem`, with `ρ` inserted.)

(⇐) (`InQuotClosure.lift` and the transport of the chain part.) First, every element of
`V_D^q` is in `V_D(F)`: inputs `ρ(f) = f − (f − ρ f)`, with `f ∈ V_D(F)` an input and
`f − ρ(f) ∈ InClosure G (deg f) ⊆ InClosure F D` (`InClosure.transport`: the derivation from
`G` within budget `deg f ≤ D` is a derivation from `F ⊇ G` within budget `D`, because the
elements of `G` are inputs of degree `≤ D`); products `ρ(X_v w) = X_v w − (X_v w − ρ(X_v w))`
with `X_v w ∈ V_D(F)` (degree `1 + deg w ≤ D`) and the correction transported as before.
Then for `p` with `deg p ≤ D` and `ρ(p) ∈ V_D^q`: `p = (p − ρ p) + ρ(p)`, the first summand in
`InClosure G (deg p)` transported into `InClosure F D`, the second in `V_D(F)` by the lift. ∎

**Corollary 3.3.** With the hypotheses of 3.1 and `1 ∈ Red` (true for `G = C`):

* `1 ∈ V_D(F) ⟺ 1 ∈ V_D^q(F)` (`one_mem_iff`);
* for reduced `p`: `p ∈ V_D(F) ⟺ p ∈ V_D^q(F)` (`mem_closure_iff_of_reduced`);
* `V_D(F) = (⟨G⟩ ∩ {deg ≤ D}) ⊕ V_D^q(F)` as `K`-spaces: `p ↦ (p − ρ p, ρ p)`; the two
  summands meet in `0` because a reduced element of `⟨G⟩` is `0`. Hence
  `dim V_D(F) = dim(⟨C⟩ ∩ M_{≤D}^{ord}) + dim V_D^q(F)` where `M^{ord}_{≤D}` is the full
  degree-`≤D` part of `R`. (Stated here; the Lean file proves the iff, from which this is
  immediate linear algebra.)

**Theorem 3.4 (lower-bound witness)** (`not_one_of_witness`, `not_one_of_chain_witness`,
`prime_m2_not_one`). Let `G` be reduced-ideal trivial with `0 ∉ G` and `1 ∈ Red`. Let
`W ⊆ R` be a `K`-subspace such that `ρ(f) ∈ W` for every `f ∈ F` of degree `≤ D`,
`ρ(X_v w) ∈ W` for every `w ∈ W` of degree `< D` and every variable `v`, and `1 ∉ W`. Then
`1 ∉ V_D(F)`. If moreover `1 ∈ ⟨F⟩` then `¬HasLFDBound F e` for every `e ≤ D`
(`lfd_gt_of_witness`, `prime_m2_lfd_gt`). *Proof.* By `closure_normalForm_mem`,
`1 ∈ V_D(F)` gives `ρ(1) = 1 ∈ V_D^q(F)`, and `V_D^q(F) ⊆ W` by the least-space property of
the quotient closure (`quot_closure_minimal`). ∎

This direction needs neither `G ⊆ F` nor a degree bound on `G`; for `F_k = C_k ∪ {g}` the
input condition reduces to `ρ(g) ∈ W` because `ρ(c) = 0` for chain equations
(`normalForm_chain_eq_zero`), which is how `prime_m2_not_one` is stated.

**Remark (why `D ≥ 2` in the README's premise).** For `D ≤ 1` no input has degree `≤ D`, so
`V_D(F_k) = 0` while the formula of 3.3 would need the chain equations as inputs. The
decomposition of 3.3 is the exact content of "closure modulo the chains equals the ordinary
closure for `D ≥ 2`": the quotient closure is the ordinary closure with its chain part
`⟨C⟩ ∩ {deg ≤ D}` removed, and the two agree on reduced elements.

## 4. The two-level exponent model and the carry lemmas

Fix a base depth `k0 ≥ 1` and consider depths `k ≥ k0`. Every exponent `a_i ∈ [0, 2^k)`
is uniquely `a_i = l_i + 2^{k0} q_i` with the **low part** `l_i ∈ [0,2^{k0})` and the **high
pattern** `q_i ∈ [0, 2^{k−k0})`. The weight splits:

**Lemma 4.1** (`BinaryWeight.wt_add_mul_two_pow`). For `l < 2^n`: `wt(l + 2^n q) = wt l + wt q`.

Lean works with pairs `(l, q) ∈ (Fin m → ℕ) × (Fin m → ℕ)` (`DepthReduction.Ex m`), the
depth-`k` **box** `InBox k0 k (l,q) ⟺ ∀ i, l_i < 2^{k0} ∧ q_i < 2^{k−k0}`, the weight
`pwt(l,q) = Σ_i wt(l_i) + Σ_i wt(q_i)` (`lowWt`, `pwt`), and the algebra `Al m K = Ex m →₀ K`
(the depth-`k` algebra is the subspace of elements supported in the box, `Boxed k0 k`). The
translation by `u_{i,j}` at depth `k` is the linear map `translate k0 k i j = lmapDomain
(step k0 k i j)` where, writing `N = 2^{k−k0}`,

    step (l,q) = (l + 2^j e_i, q)                                    if j < k0 and l_i + 2^j < 2^{k0}
               = (l + 2^j e_i − 2^{k0} e_i, q + e_i mod N in block i)   if j < k0 and l_i + 2^j ≥ 2^{k0}
               = (l, q + 2^{j−k0} e_i mod N in block i)                 if j ≥ k0.

**Lemma 4.2 (model soundness).** Under `a_i = l_i + 2^{k0} q_i`, `step` is addition of
`2^j e_i` modulo `2^k`: `encode(step(l,q))_i = (a_i + 2^j) mod 2^k`, other coordinates
unchanged. *Proof.* Case `j < k0`, no carry: `a_i + 2^j = (l_i + 2^j) + 2^{k0} q_i ≤
(2^{k0} − 1) + 2^{k0}(2^{k−k0} − 1) = 2^k − 1`, no reduction needed. Case `j < k0`, carry:
`a_i + 2^j = (l_i + 2^j − 2^{k0}) + 2^{k0}(q_i + 1)`; if `q_i + 1 < N` this is `< 2^k` and is
the stated pair; if `q_i + 1 = N` it equals `(l_i + 2^j − 2^{k0}) + 2^k`, whose residue is the
low part with pattern `0 = (q_i + 1) mod N`. Case `j ≥ k0`: `a_i + 2^j = l_i + 2^{k0}(q_i +
2^{j−k0})`, and `2^{k0}(x mod N) ≡ 2^{k0} x (mod 2^k)` with `l_i + 2^{k0}(x mod N) < 2^k`. ∎
(Formalized as `ChainModelMap.encode_step`, with `encode k0 (l,q) i = l_i + 2^{k0} q_i`.)

On the box, `step` is a bijection with inverse `unstep` (`step_unstep`, `unstep_step`,
`step_inBox`, `unstep_inBox`), so `(translate w)(b) = w(unstep b)` for boxed `w` and `b` in
the box, and `= 0` for `b` outside it (`translate_apply`, `translate_apply_of_not_inBox`).

**Carry lemmas** (`BinaryWeight.lean`, all kernel-checked):

* (C1) `wt(a + 2^j) ≤ wt(a) + 1` (`wt_add_two_pow_le`): one added bit, any carries.
* (C2) If bit `j` of `a` is set, `wt(a + 2^j) ≤ wt(a)` (`wt_add_two_pow_le_of_odd`).
* (C3) For `r < 2^n`: `wt(2^n + r) = wt(r) + 1` (`wt_two_pow_add`).
* (C4) Carry out of an `n`-bit window: if `a < 2^n`, `j < n`, `a + 2^j ≥ 2^n`, then bit `j` of
  `a` is set (indeed bits `j..n−1` are all set, `carry_odd`), and
  `wt(a + 2^j − 2^n) + 1 ≤ wt(a)` (`carry_wt`): the window remainder loses weight.
* (C5) Wrap-around in a cyclic `n`-bit block: for `q < 2^n`, `j < n`,
  `wt((q + 2^j) mod 2^n) ≤ wt(q) + 1` (`wt_mod_two_pow_add_le`): either no wrap and (C1),
  or `q + 2^j = 2^n + r` and `wt r = wt(q + 2^j) − 1 ≤ wt q` by (C3), (C1).

**Lemma 4.3** (`pwt_step_le`). For `(l,q)` in the depth-`k` box, `k > k0`, `j < k`:
`pwt(step(l,q)) ≤ pwt(l,q) + 1`. *Proof.* `lowWt(update l i y) + wt(l_i) = lowWt(l) + wt(y)`
(`lowWt_update`). No-carry low step: (C1). Carry low step: the low part loses at least one
by (C4) and the pattern gains at most one by (C5) with `j = 0`. High step: (C5). ∎

## 5. The free high extension `W_k(S)`

Let `S` be a `K`-subspace of the **low algebra** `Lo m K = (Fin m → ℕ) →₀ K` (elements of
depth `k0`). For `q ∈ Fin m → ℕ` the **`q`-component** of `w ∈ Al m K` is the low element
`piece q w : l ↦ w(l, q)` (`piece`, linear), and `embed q v : (l,q) ↦ v(l)` places a low
element at pattern `q` (`embed`, linear), with `piece q (embed q v) = v` and
`piece q' (embed q v) = 0` for `q' ≠ q`. In the notation of the README, `embed q = τ_h` with
`h = 2^{k0} q`, and `piece q` is the projection onto the `h`-component followed by `τ_{−h}`.

**Definition 5.1** (`DepthReduction.W k0 D S k`).

    W_k(S) = { w ∈ Al m K : Boxed k0 k w  ∧  ∀ q, piece q w ∈ S  ∧
                 ∀ l ∈ supp(piece q w), wt(l) + wt(q) ≤ D }.

This is a subspace (supports of sums and scalar multiples shrink; `piece` is linear). It is
the component-wise description of `⊕_{q} embed q (S ∩ M_{≤ D − wt q})`: an element `w` of
the depth-`k` box is the sum of its components `embed q (piece q w)` over the finitely many
patterns `q`, the components have disjoint supports, and the weight of a monomial `(l,q)` is
`wt l + wt q` (Lemma 4.1), so "`w ∈ W_k`" says exactly that the `q`-component is `τ_q v_q`
with `v_q ∈ S` and `wt v_q ≤ D − wt q`. Patterns with `wt q > D` have zero component. For
`k = k0` only `q = 0` occurs and `W_{k0}(S) = S ∩ M_{≤D}`.

Two immediate facts:

* **(ii)** `δ_0 = 1 ∈ W_k(S) ⟺ δ_0 ∈ S` (`single_zero_mem_W_iff`): the `0`-component of
  `δ_{(0,0)}` is `δ_0`, every other component is `0`, and `wt 0 + wt 0 = 0 ≤ D`.
* `embed 0 v ∈ W_k(S)` for `v ∈ S` supported in the low box with `wt ≤ D`
  (`embed_zero_mem_W`).

**Definition 5.2.** `W` is **closed at depth `k`** (`Closed k0 k D W`) if for every `w ∈ W`
all of whose monomials have weight `< D` (`WtLT D w`), and every `i < m`, `j < k`,
`translate k0 k i j w ∈ W`. By Lemma 2.2 and Lemma 4.2 this is, for the image in `R_k`,
exactly the hypothesis "`ρ(X_v w) ∈ W` for `w ∈ W` of degree `< D`" of Theorem 3.4, and it is
what `check_containing_space` verifies (all `2k` shifts of every row of pivot weight `< D`;
the RREF with weight-sorted columns makes this equivalent to the condition on all elements of
weight `< D`, as in `FilteredSpan.rowSpace_variable_stable`).

## 6. The reduction theorem

**Theorem 6.1 (N5, general form).** Let `k0 ≥ 1`, `D ≥ 0`, `m ≥ 1`, and let `S ⊆ Lo m K`
be any subspace. Then:

* **(ii)** for every `k ≥ k0`: `1 ∈ W_k(S) ⟺ 1 ∈ S` (`single_zero_mem_W_iff`);
* **(iii)** if `W_{k0+1}(S)` is closed at depth `k0+1`, then `W_k(S)` is closed at depth `k`
  for every `k ≥ k0+1` (`W_closed`);
* **(i)** if `g ∈ S` for a low element `g` of weight `≤ D`, and `W_{k0+1}(S)` is closed, then
  for every `k ≥ k0+1` and every exponent vector `μ ∈ [0,2^k)^m` with `wt μ ≤ D − wt g`,
  the generator row `τ_μ g ∈ W_k(S)`.

*Proof of (iii).* Let `w ∈ W_k(S)` with all monomials of weight `< D`, and fix a variable
`(i, j)`, `j < k`. Write `w' = translate k0 k i j w`. We check the three clauses of 5.1.

*Box.* `w'` is supported on `step(supp w) ⊆ box` (`step_inBox`).

*Weight budget.* If `(l,q) ∈ supp w'` then `(l,q) = step(a)` for some `a ∈ supp w`
(support of a `mapDomain` lies in the image), and `wt l + wt q = pwt(step a) ≤ pwt(a) + 1 ≤ D`
by Lemma 4.3. This clause does not look at the components at all.

*Components in `S`.* Fix a pattern `q`. If `q` is outside the high box, `piece q w' = 0 ∈ S`
because `w'` is boxed. Otherwise:

*High variable, `j ≥ k0`* (`piece_translate_high`). `step` fixes the low part and permutes
the patterns by `σ : q ↦ q + 2^{j−k0} e_i mod N`, a bijection of the high box with inverse
`σ^{-1} : q ↦ q + N − 2^{j−k0} mod N` in block `i`. For low `l`:
`(piece q w')(l) = w'(l,q) = w(unstep(l,q)) = w(l, σ^{-1} q) = (piece (σ^{-1} q) w)(l)`;
for `l` outside the low box both sides vanish. So `piece q w' = piece (σ^{-1} q) w ∈ S`.
(The weight bound `wt l + wt(σ q) ≤ wt l + wt q + 1 ≤ D` for `wt l + wt q < D` is (C5); in
the formal proof it is subsumed by the budget clause above.)

*Low variable, `j < k0`* (`piece_translate_low`). Let `q⁻ = q − e_i mod N` in block `i`
(`qPrev`). Claim:

    piece q w' = piece 0 (T (embed 0 (piece q w))) + piece e_i (T (embed 0 (piece q⁻ w))),

where `T = translate k0 (k0+1) i j` is the same variable acting at depth `k0+1`, whose high
box is `{0,1}^m`. Evaluate both sides at a low `l`. The left side is `w(unstep_k(l,q))`.
If `2^j ≤ l_i` (no carry on the way back): `unstep_k(l,q) = (l − 2^j e_i, q)`, the first term
is `(embed 0 (piece q w))(unstep_{k0+1}(l,0)) = (embed 0 (piece q w))(l − 2^j e_i, 0) =
w(l − 2^j e_i, q)`, and the second term is `(embed 0 …)(l − 2^j e_i, e_i) = 0` because
`e_i ≠ 0`. If `l_i < 2^j` (a carry occurred): `unstep_k(l,q) = (l + (2^{k0} − 2^j) e_i, q⁻)`;
the first term is `(embed 0 …)(l + (2^{k0} − 2^j) e_i, (0 + 2 − 1) mod 2 · e_i) =
(embed 0 …)(…, e_i) = 0`, and the second is `(embed 0 (piece q⁻ w))(l + (2^{k0}−2^j) e_i,
(1 + 2 − 1) mod 2 · e_i) = (embed 0 (piece q⁻ w))(…, 0) = w(l + (2^{k0}−2^j) e_i, q⁻)`.
For `l` outside the low box all three terms vanish (boxedness at depths `k` and `k0+1`).
This proves the claim.

Now `v := piece q w ∈ S` with support in the low box (as `w` is boxed) and every monomial
`l` of `v` satisfies `wt l ≤ wt l + wt q < D`; so `embed 0 v ∈ W_{k0+1}(S)`
(`embed_piece_mem_W`) and all its monomials have weight `< D` (`embed_piece_wtLT`). By the
closedness hypothesis at depth `k0+1`, `T(embed 0 v) ∈ W_{k0+1}(S)`, whose `0`- and
`e_i`-components lie in `S`. The same holds for `q⁻`. Hence `piece q w'` is a sum of two
elements of `S`. ∎

*Proof of (i).* `τ_μ g` is obtained from `g` by `wt μ` single-variable translations (one per
set bit of `μ`, in any order); after `s` of them the element has weight `≤ wt g + s`, which
is `< D` as long as `s < wt μ ≤ D − wt g`, so each step stays in `W_k(S)` by (iii), starting
from `embed 0 g ∈ W_k(S)` (`embed_zero_mem_W`). ∎ (This is the README's statement (i);
note that it is not needed below. The README proves (i) without the closedness hypothesis
by the no-overflow bound `μ_low + s ≤ 96 + 2 < 128`; that argument is correct for `k0 = 7`,
`D = 4` and `g` supported on `{0,1,2}^2`, and the general form needs exactly the condition
that `supp g + {μ_low : wt μ_low ≤ D − wt g}` stays inside `[0,2^{k0})^m`.)

**Where the hypotheses enter.** (iii) uses only: `S` is a subspace (components are added),
and the closedness of `W_{k0+1}(S)`. It does not use that `S` is a closure, that `S` is
closed under anything, its dimension, or `S ⊆ M_{≤D}` (the budget is part of the definition of
`W`). The depth `k0+1` is forced by the low-variable case: a low variable can carry out of the
low window into exactly one new high bit, so the depth-`(k0+1)` algebra (patterns in `{0,1}^m`)
sees every effect a low variable can have on the components; a high variable only permutes
components and needs no certificate at all. The base depth `k0` itself is a choice; on the
primary instance `k0 = 6` fails because `W_7(V_4(6))` is not closed (README 4.3), and
`k0 = 7` passes.

## 7. Consequences in the ordinary ring

Fix `m = 2`, `k0 = 7`, `D = 4`, and an instance `(K, a, b, t)`, `2 ≠ 0` in `K`. Let
`S = V_4^q(F_7) ⊂ A_7` (the Python `V_4(7)`), read as a subspace of the low algebra via
`m_a ↔ δ_a`, and `W_k = W_k(S)` as in 5.1. Assume the two **base facts**:

* (B1) `1 ∉ S` (the RREF of `V_4(7)` has no pivot at the constant column; `829 × 3473`);
* (B2) `W_8` is closed at depth `8` (`check_containing_space(W_8)`: all `16` translates of
  every weight-`≤3` row of `W_8` reduce to `0` against `W_8`; `1061 × 2517`).

The third fact used, `g ∈ S`, is automatic (`g` is a generator row of `V_4^q(F_7)`;
`check_containing_space` re-checks that `W_8` contains all generator rows).

**Theorem 7.1.** For every `k ≥ 8`: `1 ∉ V_4(F_k)` in the ordinary ring `R_k`. Hence, if
`1 ∈ ⟨F_k⟩`, then `¬HasLFDBound F_k e` for `e = 0,…,4`, i.e. `LFD(F_k) ≥ 5`.

*Proof.* By Theorem 6.1 (iii) with (B2), `W_k` is closed at depth `k`; by (ii) and (B1),
`1 ∉ W_k`. Let `Φ_k : Al 2 K → R_k` be the linear map `δ_{(l,q)} ↦ m_{l + 2^7 q}` (the
squarefree monomial of the exponent vector), injective on the depth-`k` box. By Lemmas 2.2
and 4.2, `ρ(u_{i,j} Φ_k w) = Φ_k(translate 7 k i j w)` for boxed `w`, and
`deg Φ_k(w) = max pwt(supp w)` (Lemma 4.1). So `W := Φ_k(W_k)` is a subspace of `Red` with:
`ρ(g) = g = Φ_k(embed 0 g) ∈ W`; for `w' = Φ_k(w) ∈ W` of degree `< 4`, every monomial of
`w` has weight `< 4` and `ρ(X_v w') = Φ_k(translate w) ∈ W`; and `1 = Φ_k(δ_{(0,0)}) ∉ W`
since `Φ_k` is injective on the box and `δ_{(0,0)} ∉ W_k`. Theorem 3.4 (`prime_m2_not_one`)
gives `1 ∉ V_4(F_k)`, and `prime_m2_lfd_gt` the LFD statement. ∎

In Lean: `DepthReductionBridge.not_one_of_base` proves this for any `Φ_k` with the three
properties (`ModelMap`: reflects `1`, reflects degree `< D`, intertwines `translate` with
`ρ(u_{i,j}·)`); `ChainModelMap.Φ` is the concrete map `δ_{(l,q)} ↦ m_{l + 2^{k0} q}`
(`lmapDomain` of `digits k ∘ encode k0`), `ChainModelMap.modelMap` proves the three
properties (`Φ_one`, `Φ_degree`, `Φ_step`), and `ChainModelMap.n5_lfd_ge_five` /
`prime_field_n5_lfd_ge_five` are Theorem 7.1 with only the base facts as hypotheses; the
generator-in-base hypothesis is stated with `g` itself because `g` is reduced
(`normalized_reduced`, `normalForm_normalized`). See Section 9.

**Theorem 7.2 (dimension).** For every `k ≥ 8`, `V_4^q(F_k) ⊆ W_k`, hence

    dim V_4^q(F_k) ≤ dim W_k = Σ_q dim(S ∩ M_{≤ 4 − wt q}) = s_4 + 2(k−7) s_3 + C(2k−14, 2) s_2
                                   + C(2k−14, 3) s_1 + C(2k−14, 4) s_0,

`s_d = dim(S ∩ M_{≤d})`; on the primary instance `(s_0,…,s_4) = (0,0,2,115,829)`, giving
`829 + 230(k−7) + 2·C(2k−14,2) = 4k^2 + 172k − 571`, and `V_4^q(F_k) = W_k` is verified for
`k ≤ 12`. In the ordinary ring, `dim V_4(F_k) = dim(⟨C_k⟩ ∩ M^{ord}_{≤4}) + dim V_4^q(F_k)`
(Corollary 3.3), so the quotient-side number is the one that carries the instance
information; the chain part is instance-independent. *Proof of the inclusion.* `W_k` is
closed and contains `ρ(g) = g`, so it contains `V_4^q(F_k)` by `quot_closure_minimal`. The sum
is direct because the components have disjoint supports; the count is the number of patterns
`q` of weight `w` in the `2(k−7)` high positions, `C(2k−14, w)`. ∎

**Corollary 7.3 (proper ideals).** Let `Z ⊂ (μ_{2^k})^2` be the solution set of the instance
(finite, `|Z| ≤ 2^{k+1}` by N4a). The ideal `⟨F_k⟩` contains every reduced polynomial
vanishing on `Z` (the chain algebra is `K^{Z}`-valued: `R_k/⟨C_k⟩ ≅ K[(Z/2^k)^2]` is reduced
with all `4^k` characters rational when `2^k | |K^×|`), so `dim(⟨F_k⟩ ∩ M_{≤4}) ≥ |M_{≤4}| − |Z|`.
If `|M_{≤4}| − |Z| > dim W_k` there is a reduced `f ∈ ⟨F_k⟩` of degree `≤ 4` with `f ∉ W_k ⊇
V_4^q(F_k)`, hence `f ∉ V_4(F_k)` (Corollary 3.3), so `HasLFDBound F_k 4` fails and
`LFD(F_k) ≥ 5`. (The evaluation fact in this corollary is not formalized; the rest is.)

## 8. Review of the earlier written argument

The README Section 4.2 argument, read as a statement in the group algebra, is correct; the
second-agent check (`VERIFY.md`, claim 3) stands. Relative to a statement about the ordinary
closure, and in generality, the following points needed settling:

1. **The ordinary/quotient premise.** "Working modulo the chains equals the ordinary closure
   for `D ≥ 2`" was a numerically cross-checked premise (`IC_BOUNDS_PROGRAM.md` Section 5,
   last bullet). Section 3 proves the precise form: `V_D(F) = (⟨C⟩ ∩ {deg ≤ D}) ⊕ V_D^q(F)`
   for `D ≥ 2`, and `1 ∈ V_D(F) ⟺ 1 ∈ V_D^q(F)`. The proof uses only that graded-lex
   division by the chains stays within the degree budget and that the chain ideal has no
   nonzero reduced element; no Gröbner-basis statement about top forms is needed (the
   "top forms of chain-ideal elements of degree `≤ D` have no squarefree monomial" route is
   avoided by pulling back along `ρ` instead of adding `⟨C⟩ ∩ M_{≤D}`).
2. **The translation model.** That `ρ(u_{i,j} m_a) = m_{a + 2^j e_i mod 2^k}` (Lemma 2.2) was
   implicit in the README's Section 1 and in the Python grid. Its proof is the carry induction
   written out in 2.2, now kernel-checked (`ChainTranslation`), together with the two-level
   bookkeeping of Lemma 4.2 (`ChainModelMap.encode_step`).
3. **"`wt h ≤ 2`" versus "`wt h ≤ 4`".** The task statement and `IC_BOUNDS_PROGRAM.md` N5
   restrict `h` to weight `≤ 2`; the README sums over all `h` of weight `≤ 4`. They agree on
   the primary instance because `s_0 = s_1 = 0`, but the general theorem must take all `h`
   (equivalently, all patterns `q`, with patterns of weight `> D` contributing `0`): if
   `S ∩ M_{≤1} ≠ 0` for some instance, the weight-`3` patterns carry nonzero components and
   closedness of the `wt h ≤ 2` truncation would not follow. Definition 5.1 is the general
   one.
4. **Weight bound in the carried part.** The README states `wt(v_1) ≤ wt(m) − (7−j) ≤
   wt(v) − 1`; this is right, and the reason is (C4): a carry out of the window requires bits
   `j..6` all set, so the window remainder loses `7 − j ≥ 1` bits while the pattern gains one.
   In the proof above this bound is not even needed for closedness: the budget clause follows
   from the uniform Lemma 4.3 (`+1` per step), because the budget of `W` is `wt l + wt q ≤ D`
   on monomials and the input had weight `< D`. The sharper bound matters only for the
   dimension bookkeeping, where it is still true.
5. **The top-level wrap.** For a high variable `j = k−1` acting on a pattern with bit `k−1−k0`
   set, `h + u` wraps to a smaller pattern; (C5) covers this and the weight drops. For a low
   variable whose carry reaches a pattern `q_i = N − 1` (all high bits set), the new pattern is
   `0`: Lemma 4.2 and the formula `qPrev` handle it; the component then lands at `q = 0`, that
   is in `S` itself, and the proof is unchanged. At depth `k0+1` the "wrap" of the pattern `1`
   to `0` is exactly the case `k = k0+1` of the general statement; the closedness hypothesis at
   depth `8` therefore already includes the top-level identity `u_{i,7}^2 = 1`.
6. **Generator rows are not needed.** The lower bound uses only `g ∈ S`: the generator rows
   `τ_μ g` with `wt μ ≤ 2` are in `V_4^q(F_k)` as products of `g` by variables, and
   `quot_closure_minimal` needs only the inputs. README (i) is a correct additional statement
   about `W_k`, proved here as 6.1 (i) through (iii), and unconditionally by the README's
   no-overflow argument for `k0 = 7`.
7. **Depth range.** (iii) is stated for `k ≥ 8 = k0+1`; `W_7 = S` and nothing is claimed about
   the closedness of `S` itself. The lower bound for `k = 7` is the base computation (`1 ∉
   V_4^q(F_7)`, then Corollary 3.3), and for `k ≤ 6` the separate per-depth computations.
8. **What is instance-specific.** The reduction (Theorem 6.1) is instance-free and
   `(k0, D, m)`-free. The base facts (B1), (B2) and the profile `(s_d)` are per instance
   `(p, a, b, t)`; "the same for all generic instances" is the open genericity statement N6,
   on which nothing here depends.

No counterexample to (iii) was found; the statement is confirmed as written, with the
clarification of item 3.

## 9. Formal status and what remains

Files (new, under `formal/ICFallDegrees/`; none of the existing files is modified). Every
theorem depends on `propext`, `Classical.choice`, `Quot.sound` only; no `sorry`, `axiom` or
`native_decide`.

| file | content |
|---|---|
| `ChainClosureBridge.lean` | `normalForm` and its calculus; `InQuotClosure`; `mem_closure_iff`, `mem_closure_iff_of_reduced`, `one_mem_iff`, `not_one_of_witness`, `lfd_gt_of_witness`; chain specialisations `chain_one_reduced`, `chain_mem_closure_iff`, `not_one_of_chain_witness`; two-summand specialisation `inputs_castSucc_mem_chainInputs`, `prime_m2_not_one`, `prime_m2_lfd_gt` (Section 3) |
| `BinaryWeight.lean` | `wt`, (C1)–(C5), Lemma 4.1 (`wt_add_mul_two_pow`), cyclic inverses `mod_add_sub_cancel`, `mod_sub_add_cancel` (Section 4) |
| `DepthReduction.lean` | the two-level model (`Ex`, `Al`, `InBox`, `Boxed`, `WtLT`, `step`/`unstep`, `translate`, `piece`, `embed`), Lemma 4.3 `pwt_step_le`, `W`, `single_zero_mem_W_iff` (N5 ii), `embed_zero_mem_W`, `piece_translate_high`, `piece_translate_low`, `W_closed` (N5 iii) (Sections 5, 6) |
| `DepthReductionBridge.lean` | `ModelMap` (abstract interface for `Φ_k`), `GeneratorInBase`, `not_one_of_base`, `lfd_gt_of_base`, `n5_lfd_ge_five`, `prime_field_n5_lfd_ge_five`, all relative to a `ModelMap` hypothesis (Section 7) |
| `ChainTranslation.lean` | bits (`bit`, `bit_add_two_pow_of_lt/_self/_of_gt`, `bit_sub_two_pow`, `add_two_pow_lt_of_bit_zero`, `add_top_bit_mod`), `digits`, `digits_reduced`, `digits_update_of_bit_zero`, `digits_of_bit_one`, the carry congruence `translate_monomial_mem_span`, Lemma 2.2 `normalForm_X_mul_monomial` (Section 2.2) |
| `ChainModelMap.lean` | `sum_bits`, `wt_eq_sum_bits`, `digits_degree`, `digits_injOn`, `encode`, `encode_lt`, `encode_injOn`, Lemma 4.2 `encode_step`, the concrete `Φ`, `Φ_coeff`, `Φ_one`, `Φ_degree`, `Φ_step`, `modelMap`; `normalized_reduced`, `normalForm_normalized`; the assembled `GeneratorInBase'`, `not_one_of_base`, `lfd_gt_of_base`, `n5_lfd_ge_five`, `prime_field_n5_lfd_ge_five` with no interface hypothesis (Section 7) |

To register the modules, add to `formal/ICFallDegrees.lean`:

    import ICFallDegrees.ChainClosureBridge
    import ICFallDegrees.BinaryWeight
    import ICFallDegrees.DepthReduction
    import ICFallDegrees.DepthReductionBridge
    import ICFallDegrees.ChainTranslation
    import ICFallDegrees.ChainModelMap

and to `Audit.lean`:

    #print axioms ICFallDegrees.ChainClosureBridge.mem_closure_iff
    #print axioms ICFallDegrees.ChainClosureBridge.one_mem_iff
    #print axioms ICFallDegrees.ChainClosureBridge.not_one_of_witness
    #print axioms ICFallDegrees.ChainClosureBridge.prime_m2_not_one
    #print axioms ICFallDegrees.ChainClosureBridge.prime_m2_lfd_gt
    #print axioms ICFallDegrees.BinaryWeight.carry_wt
    #print axioms ICFallDegrees.DepthReduction.pwt_step_le
    #print axioms ICFallDegrees.DepthReduction.single_zero_mem_W_iff
    #print axioms ICFallDegrees.DepthReduction.W_closed
    #print axioms ICFallDegrees.DepthReductionBridge.not_one_of_base
    #print axioms ICFallDegrees.ChainTranslation.normalForm_X_mul_monomial
    #print axioms ICFallDegrees.ChainModelMap.encode_step
    #print axioms ICFallDegrees.ChainModelMap.modelMap
    #print axioms ICFallDegrees.ChainModelMap.not_one_of_base
    #print axioms ICFallDegrees.ChainModelMap.n5_lfd_ge_five
    #print axioms ICFallDegrees.ChainModelMap.prime_field_n5_lfd_ge_five

**The assembled Lean theorem** (`ChainModelMap.n5_lfd_ge_five`): for every field `K` with
`2 ≠ 0`, every `S : Submodule K (Lo 2 K)`, every `k ≥ 8` and `(a, b, t)`,

    Closed 7 8 4 (W 7 4 S 8)  →  single 0 1 ∉ S  →  GeneratorInBase' 7 4 S k hk2 a b t  →
    1 ∈ ⟨F_k⟩  →  ∀ e < 5, ¬ HasLFDBound F_k e,

and `not_one_of_base` gives `¬ InClosure F_k 4 1` without the unit-ideal hypothesis. The
hypotheses, none of which is discharged in Lean:

* `base : Closed 7 8 4 (W 7 4 S 8)` — (B2), the depth-8 closedness certificate for the base
  `S` (the model-side image of `V_4(7)`: `1061 × 2517`);
* `one : single 0 1 ∉ S` — (B1): `1 ∉ V_4(7)` (`829 × 3473`);
* `gen : GeneratorInBase' 7 4 S k hk2 a b t` — `g ∈ S`: some `v ∈ S` supported in `[0,2^7)^2`
  with weight `≤ 4` has `Φ 7 k (embed 0 v) = normalized k hk2 a b t` (in words: the element of
  `S` whose monomials are those of `g`); this is immediate from how `S` is built (`g` is a
  generator row), but it is part of the finite data;
* `unit : 1 ∈ ⟨F_k⟩` (for the LFD form; `not_one_of_base` does without it).

(B1), (B2) and `gen` are per-instance finite facts, to be supplied as kernel-checked
certificates or as hypotheses; they are not kernel-replayable without `native_decide` at
present. Everything else — the bridge, the model identification, the depth reduction — is
proved for every field with `2 ≠ 0`, every base depth `k0`, degree `D`, and number of summands
`m` (the two-summand restriction enters only through `PrimeM2General.inputs`).

Not formalized: Theorem 7.2 (the dimension count; it needs the direct-sum decomposition of
`W_k` by components and `Corollary 3.3` as a dimension statement) and Corollary 7.3 (the
evaluation fact `dim(⟨F⟩ ∩ M_{≤4}) ≥ |M_{≤4}| − |Z|`); the README's unconditional form of
(i) (not needed).

**Next concrete step.** The base certificates. `S` is given by its RREF (829 rows of
3473 columns at depth 7); the three hypotheses are: `1 ∉ S` (no constant pivot), `g ∈ S`
(one reduction), and closedness of `W 7 4 S 8` (`16 × 119` reductions of weight-`≤3` rows
of `W_8` against its `1061 × 2517` RREF, by `FilteredSpan.rowSpace_variable_stable`-style
pivot bookkeeping). A kernel-checked version needs either a sparse certificate format for
these sizes in the spirit of `SparseLowerCertificate` (the depth-6 degree-4 certificate was
already beyond kernel replay without `native_decide`), or the genericity theorem N6 that
replaces the per-instance base by a proof over `Q(a,b,t)`. Independently, Theorem 7.2 can be
formalized from `W_closed`, `quot_closure_minimal` and a component decomposition of `W_k`.
