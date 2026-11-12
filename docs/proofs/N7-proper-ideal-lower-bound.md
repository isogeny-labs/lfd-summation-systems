# N7. The all-depth lower bound without the unit-ideal hypothesis

Status (2026-11-10). Labels: **kernel** (Lean-checked), **proof** (written out),
**exact** (certified computation, command and certificate named), **empirical**, **conjecture**.

## 1. What changed

N5 (`docs/proofs/N5-reduction.md`, `ChainModelMap.not_one_of_base`) proves `1 ∉ V_4(k)` for every
depth `k ≥ 8` of the two-summand presentation, from a depth-7 base. Its LFD consequence
`ChainModelMap.lfd_gt_of_base` / `N5CertInstance*.lfd_ge_five` assumed `1 ∈ ⟨F_k⟩`. For a
proper ideal (the primary instance `p = 40961`, `(a, b, t) = (21223, 9887, 4748)` is one) that
assumption is false and the statement was vacuous.

N7 removes the assumption. **kernel**, for both certified instances, with no hypothesis left:

```
N5CertWitnessInstanceP40961.lfd_ge_five : ∀ (k : ℕ),
  8 ≤ k → ∀ (hk2 : 2 ≤ k), ∀ e < 5, ¬HasLFDBound (Set.range (PrimeM2General.inputs k hk2 21223 9887 4748)) e

N5CertWitnessInstanceP11.lfd_ge_five : ∀ (k : ℕ),
  8 ≤ k → ∀ (hk2 : 2 ≤ k), ∀ e < 5, ¬HasLFDBound (Set.range (PrimeM2General.inputs k hk2 2 5 7)) e
```

(over `ZMod 40961` and `ZMod 11`; printed by `#check`; axioms `propext`, `Classical.choice`,
`Quot.sound`). `lfd_ge_five'` is the same with `2 ≤ k` derived from `8 ≤ k`.

These are lower bounds only: `LFD(k) ≥ 5` for `k ≥ 8`. Nothing here bounds the LFD from above.

## 2. The argument (**proof**; every step is also **kernel**)

Notation of N5: `A_k = K[(Z/2^k)²]`, exponent `(l, q)` with low part `l ∈ [0, 2^{k0})²` and high
pattern `q`, weight `wt l + wt q`, `Φ` the model map to squarefree chain monomials,
`W_k(S) = ⊕_q τ_q (S ∩ M_{≤ D − wt q})` the free extension of a base `S ⊆ K[[0,2^{k0})²]`.

1. `HasLFDBound F d` says: every `f ∈ ⟨F⟩` lies in `InClosure F (max d (deg f))`. So if some
   `h ∈ ⟨F⟩` has `deg h ≤ D` and `h ∉ V_D := InClosure F D`, then no `e ≤ D` is an LFD bound.
   The unit `1` is one possible `h`; any ideal element will do.
2. Containing space (`ShiftWitness.normalForm_mem_map_W`). Under the N5 base hypotheses
   (B2) `W_{k0+1}(S)` closed and (G) generator in `S`, the chain normal form of every element
   of `V_D(k)` lies in `Φ(W_k(S))`, for every `k ≥ k0+1`. This is the N5 proof with the last
   step (`1 ∉ W`) left out: `quot_closure_minimal` with the space `Φ(W_k(S))`, closed by
   `DepthReduction.W_closed` and `ChainModelMap.Φ_step`.
3. `Φ` is injective on boxed elements and its values are reduced (`Φ_injOn`, `Φ_reduced`).
   Hence for boxed `w ∉ W_k(S)`: `Φ w ∉ V_D(k)` (`not_inClosure_of_not_mem_W`).
4. Take `w = embed 0 v` with `v` in the low box, weights `≤ D`, `v ∉ S`. Its pattern-0 piece
   is `v ∉ S`, so `w ∉ W_k(S)`; and `deg Φ w ≤ D` (`Φ_totalDegree_le`). If also `Φ w ∈ ⟨F_k⟩`,
   step 1 applies (`lfd_gt_of_base_witness`). The hypothesis (B1) `1 ∉ S` of N5 is not used.
5. Ideal membership at every depth. Multiplication by `x_i` (the chain variable `u_{i,0}`) acts
   on `embed 0 v` as the exponent shift `l ↦ l + e_i` as long as no exponent leaves the low
   box (`translate_embed_zero`: no carry into the high pattern), and
   `Φ(translate w) = NF(u·Φ w) ≡ u·Φ w` modulo the chain equations. So the set `LowIdeal k0 g`
   of low elements generated from the generator vector `g` by linear combinations and
   in-box unit shifts maps into `⟨F_k⟩` for every `k > k0` (`lowIdeal_mem_ideal`). The
   description of `v` is the same at every depth.

Together (`ShiftWitness.lfd_gt_of_lowIdeal`): closed `W_{k0+1}(S)`, `g ∈ S`,
`v ∈ LowIdeal k0 g`, weights of `v` at most `D`, `v ∉ S` ⟹ `∀ k ≥ k0+1, ∀ e ≤ D, ¬HasLFDBound e`.

## 3. The witness for `k0 = 7`, `D = 4` (**exact**, then **kernel**)

Command: `python3 scripts/ic_bounds_cert_witness.py --tag p40961-s7 --suffix P40961`
(and `--tag p11-a2-b5-t7 --suffix P11`); record `experiments/ic-bounds-cert/<tag>/witness.json`.
`S = V_4(7)` is the certified echelon basis of the N5 certificate (829 rows, 642 free columns).
`U_E = span{x₀^{e₀} x₁^{e₁} g : e₀, e₁ < E}` (no wrap for `E ≤ 126`), `M_{≤4}` = weight `≤ 4`.

The recommended single shifted generator does **not** work, for either instance:

| candidates | count | outside `S` |
|---|---|---|
| single shifts `x₀^{e₀} x₁^{e₁} g`, `e ≤ 125`, all monomials of weight `≤ 4` (includes `e = (7,0)`) | 131 | 0 |
| `U_6 ∩ M_{≤4}` | dim 29 | none |
| `U_14 ∩ M_{≤4}` | dim 105 | none |
| `U_30 ∩ M_{≤4}` | dim 270 | none |
| `U_46 ∩ M_{≤4}` | dim 503 | none |
| `U_47 ∩ M_{≤4}` | dim 534 | 1 basis vector |
| `U_48`, `U_50`, `U_54`, `U_62` | dim 544, 554, 562, 564 | yes |

The numbers are identical for `p = 40961` and `p = 11` (the count of basis vectors outside `S`
at `E = 54` is 9 resp. 7; it depends on the basis, it is not a dimension). The least box with a
witness is `E = 47`: a product `u = h·g` with `deg_{x₀} h = deg_{x₁} h = 46`, leading monomial
`x₀^{48} x₁^{48}` (weight 4), whose echelon residual modulo `S` is nonzero in the free column
`(48, 48)`. So what was used is a **linear combination of shifted generators**, not a short
product of low variables:

| instance | cofactor terms | monomials of `u` | residual free columns | chunks |
|---|---|---|---|---|
| `p40961-s7` | 1291 | 161 | 289 | 14 |
| `p11-a2-b5-t7` | 1159 | 146 | 249 | 13 |

Kernel certificate (`formal/ICFallDegrees/N5CertWitness.lean`, generic in the data):

* `prodOK d F terms_i part_i` (one `decide +kernel` per chunk of 96 cofactor terms):
  every term keeps the generator in `[0,128)²`, and `part_i = (∑ terms_i)·gen mod p`, by sorted
  merges; `sumOK`: `u = ∑ part_i mod p`. Soundness `witness_lowIdeal`, `sum_lowIdeal`:
  `LowIdeal 7 (toLo gen) (toLo u)`.
* `outsideOK d u slot`: all entries of `u` are admissible columns of weight `≤ 4`; the echelon
  accumulator of `u` cannot overflow (`coefSum u · p < B`); free column `slot` is not repeated;
  slot `slot` of the accumulator is not divisible by `p`. Soundness `not_mem_of_outsideOK`:
  `toLo u ∉ S` (an element of `S` vanishing on all pivots is `0`, by `N5CertSound.repr`), and
  `support_of_outsideOK`: low box and weight `≤ 4`.

The 43 kernel chunks of the N5 base certificate are imported (`N5CertKernel<P>.valid`), not
re-run. The new checks take about 40 s and 2.2 GB peak per instance.

## 4. Growth schema (**kernel**, conditional)

```
GrowthSchema.growth_schema : ∀ {K : Type u_1} [inst : Field K], 2 ≠ 0 →
  ∀ (a b t : K) (k0 D : ℕ), 2 ≤ k0 →
    ∀ (S : Submodule K (DepthReduction.Lo 2 K)),
      DepthReduction.Closed k0 (k0 + 1) D (DepthReduction.W k0 D S (k0 + 1)) →
        GrowthSchema.genVec a b t ∈ S →
          ∀ (v : DepthReduction.Lo 2 K),
            ShiftWitness.LowIdeal k0 (GrowthSchema.genVec a b t) v →
              (∀ l ∈ v.support, DepthReduction.lowWt l ≤ D) →
                v ∉ S →
                  ∀ (k : ℕ) (hk2 : 2 ≤ k), k0 + 1 ≤ k →
                    ∀ e ≤ D, ¬HasLFDBound (Set.range (PrimeM2General.inputs k hk2 a b t)) e
```

`genVec a b t` is the nine-monomial exponent vector of the normalized generator
(`Φ_genVec`: its image is `normalized k hk2 a b t` at every depth), so all hypotheses are
finite and depth-independent. This is the formal statement behind "each degree `D` is
eventually exceeded": **for each `D` for which such a base `(k0, S, v)` exists**, `LFD > D` at
every depth `k ≥ k0 + 1`. The theorem constructs no base. Only `D = 4`, `k0 = 7` is certified
(Section 3, two instances). That bases exist for every `D`, and the conjectured law
`LFD(k) = ⌊(k+14)/4⌋`, remain **conjecture**.

## 5. Lean statements

`formal/ICFallDegrees/ShiftWitness.lean`:

```
ShiftWitness.not_inClosure_of_not_mem_W : … Closed k0 (k0+1) D (W k0 D S (k0+1)) → k0+1 ≤ k → 2 ≠ 0 →
  GeneratorInBase' k0 D S k hk2 a b t → Boxed k0 k w → w ∉ W k0 D S k →
  ¬InClosure (Set.range (inputs k hk2 a b t)) D (Φ k0 k w)

ShiftWitness.lfd_gt_of_witness : … Boxed k0 k w → w ∉ W k0 D S k →
  Φ k0 k w ∈ Ideal.span (Set.range (inputs k hk2 a b t)) → (Φ k0 k w).totalDegree ≤ D →
  ∀ e ≤ D, ¬HasLFDBound (Set.range (inputs k hk2 a b t)) e

ShiftWitness.lfd_gt_of_base_witness : … (∀ l ∈ v.support, LowBox k0 l) → (∀ l ∈ v.support, lowWt l ≤ D) →
  v ∉ S → Φ k0 k (embed 0 v) ∈ Ideal.span (Set.range (inputs k hk2 a b t)) →
  ∀ e ≤ D, ¬HasLFDBound (Set.range (inputs k hk2 a b t)) e

ShiftWitness.lowIdeal_mem_ideal : 2 ≠ 0 → k0 < k → ∀ (I : Ideal (P (m*k) K)), chainInputs K m k ⊆ ↑I →
  (∀ l ∈ g.support, LowBox k0 l) → Φ k0 k (embed 0 g) ∈ I → LowIdeal k0 g v → Φ k0 k (embed 0 v) ∈ I

ShiftWitness.lfd_gt_of_lowIdeal : ∀ (k0 D : ℕ) (S : Submodule K (Lo 2 K)),
  Closed k0 (k0+1) D (W k0 D S (k0+1)) → ∀ (k : ℕ), k0+1 ≤ k → ∀ (hk2 : 2 ≤ k), 2 ≠ 0 → ∀ (a b t : K),
  ∀ g ∈ S, (∀ l ∈ g.support, LowBox k0 l) → (∀ l ∈ g.support, lowWt l ≤ D) →
  Φ k0 k (embed 0 g) = normalized k hk2 a b t →
  ∀ (v : Lo 2 K), LowIdeal k0 g v → (∀ l ∈ v.support, lowWt l ≤ D) → v ∉ S →
  ∀ e ≤ D, ¬HasLFDBound (Set.range (inputs k hk2 a b t)) e
```

`LowIdeal k0 g` is the inductive predicate with constructors `gen`, `zero`, `add`, `smul` and
`shift i : LowIdeal k0 g v → (∀ l ∈ v.support, l i + 1 < 2^k0) → LowIdeal k0 g (lowShift (Pi.single i 1) v)`.

Modules: `ShiftWitness` (general bridge), `GrowthSchema` (generator vector, schema),
`N5CertWitness` (checker and soundness), `N5CertWitnessData{P40961,P11}` (generated data),
`N5CertWitnessInstance{P40961,P11}` (kernel checks, instance theorems).

## 6. Limits

* The two instance theorems hold whether or not the ideal is the unit ideal; that the
  `p = 40961` ideal is proper is not proved in Lean and is not needed.
* The search of Section 3 is exhaustive only inside the stated families (single shifts with
  `e ≤ 125`; boxes `E ≤ 62`). It does not describe all ideal elements of weight `≤ 4`.
* The witness is dense (about 1200 cofactor terms). No structural description of it was found:
  single shifts of the partial norms `N_{x₀}^s N_{x₁}^t g` (`s + t ≤ 5`) are all in `S`
  (exploratory, `p = 11`, not recorded).
