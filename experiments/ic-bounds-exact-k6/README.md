# A kernel-checked exact last fall degree above 4: `LFD = 5` at depth 6

**Result (kernel).** For the normalized two-summand prime chain presentation at depth `k = 6`
over `F_1153` with curve `y² = x³ + 139 x + 1094` and target abscissa `t = 185`,

    ICFallDegrees.ExactK6InstanceP1153.exact_lfd_five :
      IsLastFallDegree (Set.range (PrimeM2General.inputs 6 hk (139 : ZMod 1153) 1094 185)) 5

that is, `HasLFDBound F 5` and `∀ e < 5, ¬ HasLFDBound F e` for the 13 polynomials `F` (12 chain
equations and the normalized `S_3`) in 12 variables.  It is the conjunction of

    not_one_in_closure_four : ¬ InClosure F 4 1          -- lower half
    one_in_closure_five     :   InClosure F 5 1          -- upper half

through `Closure.exactLFD_of_first_one`.  Both halves are finite certificates evaluated by the
Lean kernel (`decide +kernel`; no `native_decide`, `implemented_by`, `ofReduceBool`, `sorry`
or new axiom); all theorems depend on `propext`, `Classical.choice`, `Quot.sound` only.
The unit-ideal property is not a hypothesis: it follows from the upper half.

Until now every kernel upper bound was the dimension bound `2k+1 = 13` or a small sparse
example, and the kernel lower bounds stopped at "no `e < 5` is a bound" under an uncertified
unit-ideal hypothesis.  Depth 6 is the first depth where the measured two-summand profile
(`LFD = 4` for `k ≤ 5`, `5` for `k = 6..9`) leaves 4; here both sides are in the kernel.

## The instance

| | |
|---|---|
| field | `F_p`, `p = 1153 = 18·64 + 1` (the 64th roots of unity are rational) |
| curve | `y² = x³ + 139x + 1094`, discriminant `≠ 0`, `j = 1148`, `#E(F_p) = 1138` |
| target | `t = 185`, the abscissa of a rational point; `3t² + a = 197 ≠ 0` |
| scalar chain solutions | `|Z| = 0` (exact count over the `64 × 64` grid): unit ideal |
| `dim V_4`, profile `dim(V_4 ∩ M_≤j)` | `605`, `(0, 0, 2, 111, 605)`: the generic depth-6 values of `experiments/ic-bounds-locus` |
| `dim V_5` | `1586 = |M_≤5|` (everything), profile `(1, 13, 79, 299, 794, 1586)` |

`p = 1153` was chosen because unit ideals are rare at the small primes of the sweeps (the
expected number of scalar solutions is `4096/p`: about 21 at `p = 193`; no target of the curve
`(46, 29)` over `F_193` gives a unit ideal) and typical once `p` is a few thousand; `(a, b, t)` is
the first unit instance of a seeded random search (`random.Random(1153)`, curve and rational
target drawn uniformly).  The primary instance `p40961-s7` of `experiments/ic-bounds-fall` has
`|Z| = 2` at depth 6, so it is a proper ideal there and not usable for the unit route.

Labels: the dimensions and `|Z|` above are **exact** (generator log, and independently
`crosscheck_engine.py` with the repository engine `csd.ic_bounds_fall.closure`, output in
`p1153-a139-b1094-t185/crosscheck-engine.json`); the theorem is **kernel**.  (The engine reports
`closed = false` at `D = 5`: it stops when `1` appears; its rank 1586 is all of `M_≤5`, and the
kernel upper half does not use any closedness of `V_5`.)

## Route chosen and why

The unit route: for a unit ideal, `exactLFD_of_first_one` reduces the exact value to
`InClosure F 5 1` and `¬ InClosure F 4 1`.  The alternative (`LocalReductionBound.
lfd_bound_of_local` for a proper ideal) needs a certified dimension of `V_5` against
`|M_≤5| − |Z|` and a kernel-certified solution set; the unit route needs neither.

Both halves live in the quotient by the chains, the group algebra `F_p[(Z/64)²]`: a squarefree
chain monomial is a pair of block exponents `(e₀, e₁)`, its degree is `wt(e₀) + wt(e₁)`, and
reduced multiplication by `u_{i,j}` is the translation `e_i ↦ e_i + 2^j mod 64`
(`ChainTranslation.normalForm_X_mul_monomial`); the bridge to the ordinary closure is
`ChainClosureBridge` (`prime_m2_not_one`, `one_mem_iff`, `quot_mul_mem`).

### Lower half: `1 ∉ V_4(6)` (containing space at depth 6 itself)

`S = V_4(6)` as a reduced echelon basis (605 rows, 189 free columns), in the data format and with
the membership test of the N5 certificate (`N5Cert.Data`, `N5Cert.member`: packed big-integer
rows, echelon accumulator, slotwise divisibility).  The depth-7 model of `N5CertSound` (base at
`k0 = 7`, free extension to `k0 + 1`) does not instantiate at `k0 = k = 6`, so the closure part
is a direct depth-6 soundness proof (`ExactK6LowerSound.not_one`): for each of the 111 rows of
pivot weight `≤ 3` the 12 exact translates (with the wrap `2^6 ≡ 0`, `ExactK6Cert.shiftCode`) are
tested for membership (1332 tests, plus the generator).  The echelon part of `N5CertSound` is
reused unchanged: a monomial is stored under the *complement code* `m = 4095 − (e₀ + 64 e₁)`, so
every column of weight `≤ 4` has popcount `≥ 8`, the depth-7 specific low-row list of the record
is empty, and `N5CertSound.Valid` holds with `low = []` (`toLo_mem_of_member`, `repr`, `S` apply
as they are).  Free column 0 is a dummy with code 0 (N5Cert's table cannot address free column 0
in an entry).

### Upper half: `1 ∈ V_5(6)` (derivation certificate)

A derivation in the quotient closure with 83 rows.  Row 0 is the generator; row `n ≥ 1` is

    row n  =  Σ c · τ_s(row i),      i < n,   wt(s) + deg(row i) ≤ 5,

and row 82 is `1`.  The rule is sound by `ChainClosureBridge.quot_mul_mem` (multiplying by a
monomial within the degree budget and reducing stays in the quotient closure) and the degree
bound of every row is part of the certificate.  The rows are the module generators only: the
1586-dimensional span of their translates is never stored.

* **Packing.**  A vector is one natural with `128 × 128` slots of width 32; slot `e₀ + 128 e₁`
  carries the monomial `(e₀ mod 64, e₁ mod 64)`.  Stored rows use `e₀, e₁ < 64`, so `τ_s` for
  `s = (s₀, s₁)`, `s_i < 64`, is the single shift `row <<< (32·(s₀ + 128 s₁))`: **three
  big-number operations per term, whatever the density** (39 126 terms: about `1.2·10⁵`
  big-number operations in all).  The wrap `x^64 = 1` is done once per row by the scan, which
  folds the four quadrants slot by slot, compares with the stored row modulo `p`, and checks the
  row's shape, coefficient range and degree bound (`ExactK6Cert.scanBlocks`, about `10⁵`
  small operations per row).
* **Soundness** (`ExactK6UpperSound`).  All arithmetic on packed numbers goes through one
  statement, uniqueness of digits in polynomial form: `slot_eval : (∀ i, f.coeff i < B) →
  slot B (f.eval B) c = f.coeff c` for `f : ℕ[X]`.  The accumulator is the value at `B` of an
  explicit polynomial (`accPoly`), whose coefficients are bounded by `(Σ c)·p < B` (checked),
  so its slots are the coefficients; `regroup` splits the window into quadrants and
  `nf_M_mul : normalForm (M s * M c) = M (c + s)` turns the shift into the reduced product.
  `derivation_sound` is the induction over the rows.
* **Search** (`derive_one` in the generator).  Rounds of exact linear algebra mod `p`: in the
  span of all admissible translates of the current generators, the elements of unexpectedly low
  degree that are not translates-combinations of lower budget become new generators; only the
  least new degree is kept per round (lower degree reaches further); the generators the final
  row does not depend on are dropped.  Six rounds, 102 generators found, 83 kept
  (degrees: 1 row of degree 0, 1 of degree 1, 2 of degree 2, 79 of degree 3).

## Files

| file | content |
|---|---|
| `scripts/ic_bounds_exact_k6_gen.py` | builds both certificates (self-contained exact linear algebra, no repository engine), re-runs every Lean check in Python entry by entry (`ref_check_lower`, `ref_check_upper`: the reference checker), emits the Lean data, kernel and instance modules |
| `p1153-a139-b1094-t185/certificate.json` | pivot/free codes, low rows, generator, degrees and term lists; sha256 of the packed rows |
| `p1153-a139-b1094-t185/summary.json`, `generator.log` | sizes, operation counts, sha256 of the emitted Lean files |
| `p1153-a139-b1094-t185/crosscheck-engine.json`, `crosscheck_engine.py` | `V_4`, `V_5` of the same instance with `csd.ic_bounds_fall.closure` |
| `formal/ICFallDegrees/ExactK6Cert.lean` | the checkers (plain `Nat`/`List`, structural recursion) |
| `formal/ICFallDegrees/ExactK6LowerSound.lean` | `not_one` (lower soundness) |
| `formal/ICFallDegrees/ExactK6UpperSound.lean` | `derivation_sound` (upper soundness) |
| `formal/ICFallDegrees/ExactK6DataP1153.lean` | data (generated, 4.5 MB) |
| `formal/ICFallDegrees/ExactK6KernelLowerP1153.lean`, `ExactK6KernelUpperP1153.lean` | the `decide +kernel` theorems (generated) |
| `formal/ICFallDegrees/ExactK6InstanceP1153.lean` | `not_one_in_closure_four`, `one_in_closure_five`, `exact_lfd_five` (generated) |

## Measured builds (Lean 4.19.0, this machine, one build at a time, `ulimit -v 7000000`)

| module | wall time | peak RSS | result |
|---|---|---|---|
| `ExactK6UpperSound` | 9 s | 2.1 GB | built |
| `ExactK6LowerSound` (with the kernel check `weight_fact`) | 18 s | 3.2 GB | built |
| `ExactK6DataP1153` (4.5 MB of literals) | 241 s (366 s in a second run under load) | 3.0 GB | built |
| `ExactK6KernelUpperP1153` (85 kernel theorems: 82 rows, generator, global, final) | 365 s (470 s under load) | 2.9 GB | all checks pass; `urows_all` depends on `[propext]` |
| `ExactK6KernelLowerP1153` (40 kernel theorems: 28 chunks of 4 low rows, 5 chunks of row bounds, 7 small) | 402 s under load (306 s with chunks of 8, at 4.3 GB) | 3.2 GB | all checks pass; `valid` depends on `[propext, Quot.sound]` |
| `ExactK6InstanceP1153` | 8 s | 2.0 GB | `exact_lfd_five`, `not_one_in_closure_four`, `one_in_closure_five`: `[propext, Classical.choice, Quot.sound]` |

Peak RSS is `ru_maxrss` of the `lake build` process tree.  The generator itself: 37 s, below 1 GB.

`#check` output:

    ExactK6InstanceP1153.exact_lfd_five : IsLastFallDegree
      (Set.range (PrimeM2General.inputs 6 ExactK6InstanceP1153.hk 139 1094 185)) 5
    ExactK6InstanceP1153.not_one_in_closure_four : ¬InClosure
        (Set.range (PrimeM2General.inputs 6 ExactK6InstanceP1153.hk 139 1094 185)) 4 1
    ExactK6InstanceP1153.one_in_closure_five : InClosure
      (Set.range (PrimeM2General.inputs 6 ExactK6InstanceP1153.hk 139 1094 185)) 5 1

(`hk : 2 ≤ 6`; the field is `K = ZMod 1153`.)

Kernel cost is linear in the certificate: the upper half is 82 row theorems (one
`decide +kernel` each, 2 to 9 s, the largest row has 1583 terms and peaks at 1.3 GB), the lower
half 28 chunks of 4 low rows (48 membership tests each).

## Reproduce

    # certificates + Python reference check + Lean modules (about 40 s, < 1 GB)
    python3 scripts/ic_bounds_exact_k6_gen.py --params 1153,139,1094,185 --suffix P1153 \
        --out-dir experiments/ic-bounds-exact-k6/p1153-a139-b1094-t185
    cd formal
    lake build ICFallDegrees.ExactK6LowerSound ICFallDegrees.ExactK6UpperSound
    lake build ICFallDegrees.ExactK6DataP1153            # ~4 min
    lake build ICFallDegrees.ExactK6KernelUpperP1153     # ~6 min
    lake build ICFallDegrees.ExactK6KernelLowerP1153     # ~5 min
    lake build ICFallDegrees.ExactK6InstanceP1153
    # independent engine cross-check of the dimensions
    python3 experiments/ic-bounds-exact-k6/crosscheck_engine.py

The generator is deterministic (the emitted files are byte-identical on re-runs; hashes in
`summary.json`).  Another unit instance: `--params p,a,b,t --suffix <S>` with `64 | p − 1`; the
script refuses proper ideals and instances where `1 ∉ V_5`.

## Scope

* One instance.  Nothing here is a statement about other curves, targets, primes or depths;
  the generic profile `LFD(6) = 5` remains **empirical** beyond this instance.
* The theorem concerns the specified chain encoding (`PrimeM2General.inputs`), not relation
  yield or any cost.
* Not attempted: depth `k ≥ 7` exact values (the lower half exists for `k ≥ 8` as
  `N5CertInstance*`; an upper certificate at depth 7 would need `256 × 256` slots, four times
  the row size, and a new derivation search), and proper ideals.
