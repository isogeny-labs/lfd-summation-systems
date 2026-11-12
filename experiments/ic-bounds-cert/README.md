# Kernel-checked base certificates for N5 (`1 ∉ V_4(k)` for every depth `k ≥ 8`)

The conditional theorem `ICFallDegrees.ChainModelMap.n5_lfd_ge_five` /
`ChainModelMap.not_one_of_base` (`formal/ICFallDegrees/ChainModelMap.lean`,
`docs/proofs/N5-reduction.md`) has three per-instance finite hypotheses for a base
subspace `S ⊆ F_p[(Z/2^7)^2]`:

| | hypothesis | finite content |
|---|---|---|
| (B1) | `Finsupp.single 0 1 ∉ S` | the constant monomial is not a pivot of the echelon basis of `S` |
| (B2) | `Closed 7 8 4 (W 7 4 S 8)` | for every row `v` of `S` of weight `≤ 3`, every block `i ∈ {0,1}`, level `j < 7`: the no-carry part and the carried part of `u_{i,j} v` both lie in `S` (`2·7·2·115 = 3220` membership tests) |
| (G) | `GeneratorInBase' 7 4 S k hk2 a b t` | the normalized generator `g` lies in `S`, and its sparse list maps to `normalized k hk2 a b t` under `Φ` at every depth `k` |

This directory holds the data of kernel-checked certificates for these hypotheses with
`S = V_4(7)` (the depth-7 degree-4 closure, as a reduced row echelon basis), and the Lean
modules under `formal/ICFallDegrees/` discharge them **inside the Lean kernel** (`decide +kernel`;
no `native_decide`, no `implemented_by`, no `ofReduceBool`, no axioms beyond `propext`,
`Classical.choice`, `Quot.sound`).  The result is the unconditional instance theorem

    N5CertInstance<P>.not_one_in_closure :
      ∀ k, 8 ≤ k → ∀ hk2 : 2 ≤ k, ¬ InClosure (Set.range (inputs k hk2 a b t)) 4 1

over `K = ZMod p`, and `lfd_ge_five` (`∀ e < 5, ¬ HasLFDBound … e`) under the unit-ideal
hypothesis `1 ∈ ⟨F_k⟩`, which is *not* certified here.

## Instances

| directory | p | a | b | t | basis | rows × free cols | free entries | low rows | slot width |
|---|---|---|---|---|---|---|---|---|---|
| `p11-a2-b5-t7/` | 11 | 2 | 5 | 7 | recomputed (`csd.ic_bounds_fall.closure`) | 829 × 642 | 82272 | 115 | 16 |
| `p40961-s7/` | 40961 | 21223 | 9887 | 4748 | `experiments/ic-bounds-fall/p40961-s7/k7-D4.npz` (sha256 `1a18c003…`) | 829 × 642 | 92180 | 115 | 40 |

`p40961-s7` is the primary instance of `experiments/ic-bounds-fall/` (random curve, random
target, independently verified there).  `p11-a2-b5-t7` is a small-prime instance with the same
dimension profile `(0, 0, 2, 115, 829)`; it was chosen first because its packed numbers are
smaller (16-bit slots).  Both have `3t² + a ≠ 0`; the degenerate `p = 5`, `a = 2, b = 3, t = 4`
(`3t² + a = 0`) and `p = 3` instances do **not** satisfy (B2) (`W_8(V_4(7))` is not closed),
so they are not certifiable by this route.

## Files

* `scripts/ic_bounds_cert_gen.py` — builds the certificate from the basis, re-runs every
  check in Python (a reference implementation of the Lean checker, entry by entry), and emits
  the Lean modules.  `--params p,a,b,t` recomputes the basis; the default loads the exported
  `.npz` of `p40961-s7`.
* `<instance>/certificate.json` — the data (hex-encoded packed rows, pivot and free codes,
  the sparse low rows, the generator) with the sha256 of the basis file and of the emitted
  Lean files.
* `formal/ICFallDegrees/N5Cert.lean` — the checker (plain `Nat`/`List` code, structural
  recursion only, so the kernel can evaluate it).
* `formal/ICFallDegrees/N5CertSound.lean` — the meaning of the data and the soundness
  theorems: `toLo_mem_of_member` (echelon membership), `one_not_mem` (B1), `closed` (B2),
  `Φ_embed_eCode`/`digits_two_bits` (for G), all relative to `Valid d` (the conjunction of the
  Bool checks).
* `formal/ICFallDegrees/N5CertData<P>.lean` — the data as Lean literals (generated).
* `formal/ICFallDegrees/N5CertKernel<P>.lean` — the `decide +kernel` theorems, chunked so that
  each kernel evaluation stays within a few GB, and `valid : Valid data` (generated).
* `formal/ICFallDegrees/N5CertInstance<P>.lean` — the instance theorems (generated).

## Certificate design

A monomial of the depth-7 algebra is the code `m = l₀ + 128 l₁`.  The reduced echelon basis of
`V_4(7)` (columns sorted by weight, heaviest first) has 829 pivot columns and 642 free columns.

* Rows are stored **packed**: row `idx` is one natural number whose free column `f` occupies
  bits `[w f, w (f+1))`, `B = 2^w`; the 829 packed rows sit in a binary tree indexed by `idx`.
  `S` is *defined* in Lean as the span of `δ_{pivot idx} + (slots of the packed row idx)`, so
  no lemma about the tree or about how the data was produced is needed: whatever the data is,
  the theorems are about that `S`.
* Membership of a sparse vector `u` in `S` is the echelon test `u = ∑_{pivot entries (m,x)} x·row_m`.
  Its free-column part is checked by accumulating
  `R = ∑_{pivot entries} x · packedRow + ∑_{free entries} (p − x) · B^f`
  with GMP-accelerated `Nat` arithmetic and then testing every slot of `R` for divisibility by
  `p`.  The slot sums are `< |u| · p² ≤ 512 p² < B`, so no slot overflows (`N5CertSound`
  proves this from the checked bounds).  Pivot and free hints come from a 16-bit lookup table
  indexed by code; every hint is re-verified (`piv idx = m`, `freeCode f = m`).
* The 115 rows of pivot weight `≤ 3` are also given as sparse lists; the kernel checks that
  their packing equals the packed row, applies the 28 translations
  (`shiftEntry`: block `i`, level `j < 7`, no-carry / carry), and runs the membership test on
  each image.  `N5CertSound.T_toLo` identifies these shifts with the components
  `piece 0 ∘ translate 7 8 i j ∘ embed 0` and `piece (unit i) ∘ …` of `DepthReduction`, and
  `closed_of` turns the finite checks into `Closed 7 8 4 (W 7 4 S 8)` along the lines of
  `W_closed` (`piece_translate_low/high`).

## Kernel feasibility numbers (Lean 4.19.0, this machine, one core)

Measured with `decide +kernel` on `Nat` loops with structural recursion:

* throughput: about `3·10⁴` primitive kernel operations per second (`Nat.mul/add/mod/div/land/
  shiftRight` on literals, `match`/`if` steps), i.e. ~35 µs each; ~7.5·10³ loop iterations/s
  for a 5-op loop body; big-number GMP operations on 10–40 kbit literals cost no more than
  small ones;
* memory: the kernel's whnf cache grows by ~4 kB per primitive operation and is only freed at
  the end of a declaration, so a single `decide +kernel` must stay below ~1.5·10⁶ operations
  (~6 GB); the certificate is therefore split into 43 theorems (29 chunks of 4 low rows, 9
  chunks of 100 rows for the slot bounds, 5 small ones), checked sequentially
  (`set_option Elab.async false`; with the default asynchronous elaboration the kernel checks
  run in parallel and exhaust 7 GB);
* well-founded recursion (`for` loops, `Std.Range`, `Nat.digits`) does not reduce in the kernel;
  `Array` has no kernel support; only structural recursion and GMP `Nat` literals are usable;
* literals: hex literals are parsed in quadratic time, so they are kept below ~10 kchar each
  (one packed row = 642·w bits); a 5.3 MB data module elaborates in about a minute and a half;
* operation count of the certificate: 3220 membership tests × (~80 entries × ~50 ops + 642
  slots × 4 ops) ≈ 2·10⁷ ops for (B2), plus 829 × 642 × 4 ≈ 2·10⁶ for the slot bounds of the
  rows: about 15–20 minutes of kernel time per instance.

## Measured runs

| module | wall time | peak RSS | result |
|---|---|---|---|
| `N5CertDataP11` (2.3 MB of literals) | 33 s | — | built |
| `N5CertKernelP11` (43 kernel theorems, sequential) | 777 s | 4.65 GB | all checks pass; `valid` depends on `[propext, Quot.sound]` |
| `N5CertInstanceP11` | < 1 min | — | `not_one_in_closure`, `lfd_ge_five`: `[propext, Classical.choice, Quot.sound]` |
| `N5CertDataP40961` (5.5 MB of literals) | 83 s | — | built |
| `N5CertKernelP40961` (43 kernel theorems, sequential) | 830 s | 5.02 GB | all checks pass; `valid` depends on `[propext, Quot.sound]` |
| `N5CertInstanceP40961` | 7 s | — | `not_one_in_closure`, `lfd_ge_five`: `[propext, Classical.choice, Quot.sound]` |

## Reproduce

    cd formal
    lake build ICFallDegrees.N5CertSound
    lake build ICFallDegrees.N5CertKernelP11      # ~15 min, one kernel check at a time
    lake build ICFallDegrees.N5CertInstanceP11
    lake build ICFallDegrees.N5CertKernelP40961   # ~20 min
    lake build ICFallDegrees.N5CertInstanceP40961

and `#print axioms ICFallDegrees.N5CertInstanceP11.not_one_in_closure` (see `Audit.lean`
lines in the task report).  To regenerate the data:

    python3 scripts/ic_bounds_cert_gen.py --params 11,2,5,7 --tag p11-a2-b5-t7 --suffix P11 \
        --out-dir experiments/ic-bounds-cert/p11-a2-b5-t7
    python3 scripts/ic_bounds_cert_gen.py --tag p40961-s7 --suffix P40961 \
        --out-dir experiments/ic-bounds-cert/p40961-s7
