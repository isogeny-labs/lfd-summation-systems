# Verification of ITEM 1, three summands (m = 3)

Verifier report, 2026-11-08. I checked the producer's claims in `README.md` and the result files
in this directory with code I wrote myself. I did not import `csd/ic_bounds_m3*.py` or any closure
engine in the repository. The one exception is the separate ordinary-ring spot check, which feeds
inputs I built myself to the allowed validator `experiments/ic-curve-validation/code/lfd_engine.py`.
I ran the producer's scripts only to see whether the README commands run and reproduce. Those runs
used a copy in the scratchpad, so no producer file was overwritten. The producer files still have
the hashes listed in `README.md`, and git state is unchanged.

**Bottom line.** No exact or kernel claim was refuted. The verifier recomputed 1092 exact LFD
values (direct 370, split 261, m2 461) with an independent engine and found 0 mismatches. They
match in LFD, unit or proper status, `dim R/I` and the quotient rank at every degree. The written
proofs (local-to-global criterion, Proposition 1, Corollary 2, Remark 3) check out step by step.
The symbolic identities were re-derived from scratch. Two statements in the claim summaries are
imprecise:

* claim 6 says "316 cases" where there are 318;
* claim 9 says the split minus m2 gap is +2 at k = 3, but it is +1 in 6 of 117 cases. The README
  table has the right breakdown.

## What I wrote (all new files)

| file | role | sha256 (16) |
|---|---|---|
| `scripts/ic_bounds_m3_verify_symbolic.py` | sympy: `S_4 = Res_X(S_3(x1,x2,X), S_3(x3,x4,X))`, identities, numeric decompositions with my own group law | b81bc22533c5dae4 |
| `scripts/ic_bounds_m3_verify_engine.py` | independent closure engine modulo the chains, plus a free `z`; exact float64 mod-p RREF; `R/I` by CRT; local-to-global certificate | e896c50473485a81 |
| `scripts/ic_bounds_m3_verify_run.py` | stratified sample (seed 99173) and resumable runner | 5870019164887cda |
| `scripts/ic_bounds_m3_verify_compare.py` | compare with `rows.jsonl` | 9a338ce1b5b912f0 |
| `scripts/ic_bounds_m3_verify_extra.py` | case data, planted points, Corollary 2 hypotheses per case | 155ccc6a0a526619 |
| `scripts/ic_bounds_m3_verify_decomp.py` | every rational decomposition lies in Z (all 461 cases) | 3eef594d49c3ab74 |
| `scripts/ic_bounds_m3_verify_newprimes.py` | refutation attempt off the grid (seeds 424242, 424243) | 361d57e9326f311e |
| `scripts/ic_bounds_m3_verify_ordinary.py` | ordinary-ring spot check (validator engine, my inputs; seed 31337) | 38d0ccd9a6517a9a |
| `verify/s4_coeffs.json` | my S_4 over Z[a,b,x1..x4], 540 terms (sha of the term list a0e248de11eb31da) | 5ac0a6a994a38bc0 |
| `verify/symbolic_verify.json`, `verify/extra.json`, `verify/decomp.json` | check outputs | 08cd48b9 / 53b13287 / 84d651a3 |
| `verify/sample.json`, `verify/recheck.jsonl`, `verify/compare.json` | sample, 1092 recomputed profiles, comparison | 86c7443c / 34072baf / c5ef71d9 |
| `verify/newprimes.jsonl`, `verify/ordinary.jsonl` | off-grid instances (27), ordinary-ring checks (38) | e9a75b6f / d976a8f6 |

### How the engine works, and why it is sound

* **Space.** The monomials are chain masks `x_i^{e_i}` with `e_i mod 2^k`, of degree
  `sum popcount(e_i)`, times `z^e` for the split. Columns are sorted by degree, highest first. So
  in the reduced row echelon form (pivot = leftmost nonzero), `V cap M_<=e` is exactly the span of
  the rows whose pivot has degree `<= e`.
* **Closure.** A row with degree `<= D-1` is multiplied by every variable; for `u_{i,j}` that adds
  `2^j` to `e_i` mod `2^k`, and for `z` it adds 1 to the z exponent. This is the definition of
  `V_D` modulo the chains, which equals the ordinary closure for `D >= 2`. Every inserted row with
  low degree has its products queued, which closes the span by linearity. Exactness: all entries
  are below p <= 7681, so a product sum is at most n·p² < 2^53.
* **R/I by CRT.** The chain ideal is radical with zero set `H^m`. Hence `R/I = F_p^Z` for direct
  and m2, and `(+)_x F_p[z]/gcd(g1(x,z), g2(x,z))` for the split. For `I_field`, each factor is
  replaced by its gcd with `z^p - z`. The code stops if any fibre is infinite; none was.
* **Early stop.** The saturation stops early only when its rank reaches
  `|M_<=D| - rank(ev(M_<=D) -> R/I)` = `dim(I cap M_<=D)`. That number is an upper bound for
  `dim V_D`, and the code asserts it is never exceeded. Otherwise the closure runs to a fixed point.
  At k <= 2, and on 14 regression cases, it also re-checks closedness: every shift of every low row
  reduces to 0.
* **Certificate.** A case is certified at the first `D` where `M_<D` maps onto `R/I` and
  `V_D = I cap M_<=D`. The LFD is then the last strict fall. For a unit ideal it is the least `D`
  with 1 in `V_D`.

## Runs (all with `ulimit -v 3000000`, single-threaded BLAS)

```
python3 scripts/ic_bounds_m3_verify_symbolic.py                     # 3 s
python3 scripts/ic_bounds_m3_verify_run.py select                    # sample.json (seed 99173; 4 k=4 split + 1 k=4 direct added by hand, recorded in sample.json)
python3 scripts/ic_bounds_m3_verify_run.py run --k 1 2               # 16 s
python3 scripts/ic_bounds_m3_verify_run.py run --k 3                 # ~12 min
python3 scripts/ic_bounds_m3_verify_run.py run --k 4 --budget-min 23 # m2 all + direct; resumable, rerun until "done"
python3 scripts/ic_bounds_m3_verify_run.py run --k 4 --enc split --budget-min 23
python3 scripts/ic_bounds_m3_verify_run.py run --k 4 --enc direct --primes 257 641 7681 --budget-min 23
python3 scripts/ic_bounds_m3_verify_run.py run --k 4 --enc direct --budget-min 22
python3 scripts/ic_bounds_m3_verify_compare.py                       # compare.json
python3 scripts/ic_bounds_m3_verify_extra.py; python3 scripts/ic_bounds_m3_verify_decomp.py
python3 scripts/ic_bounds_m3_verify_newprimes.py; python3 scripts/ic_bounds_m3_verify_newprimes.py --k4
python3 scripts/ic_bounds_m3_verify_ordinary.py
```

Times on the shared, loaded host: k = 4 direct took 50 to 172 s per case and k = 4 split 166 to
282 s. Peak memory stayed well below the 3 GB cap. The `verify/` directory is 872 KB.

## Recomputation coverage and agreement (`verify/compare.json`)

| k | direct | split (I_free; I_field dim and status) | m2 |
|---|---|---|---|
| 1 | 111/111 agree (all) | 111/111 (all) | 111/111 (all) |
| 2 | 116/116 (all) | 116/116 (all) | 116/116 (all) |
| 3 | 117/117 (all) | 29/117: one random case per (prime, curve class), plus the 6 LFD-5 cases | 117/117 (all) |
| 4 | 26/117: at least 4 per prime, all 11 deg-5 cases | 5/117: one per prime (p = 17, 97, 257, 641, 7681) | 117/117 (all) |

For each recomputed case, the LFD, the status, `dim R/I` (`|Z|` or `sum deg gcd`) and the
quotient rank `q_D` at every computed `D` are identical to `rows.jsonl`. For split cases the
`I_field` dimension and the `>= p` status are identical too. There are 0 mismatches. This sample
includes all 6 LFD-0 and all 5 LFD-3 cases at k = 1, all 11, 10 and 11 deg-5 cases at k = 2, 3, 4,
and the README's representative case `4898841d8890a6d4` (direct `|Z| = 24`, LFD 7; split
`dim A = 36`, LFD 6).

## Further independent checks

* **S_4** (`verify/symbolic_verify.json`). I computed the resultant myself. It is symmetric in all
  four arguments, has degree 4 in each and total degree 12, and has 439 monomials in x1..x4, or 125
  after setting x4 = t. Both identities hold over Z[a,b,t]:
  - coefficient of x1³x2³x3³ = `-8(5t³+3at+6b)`;
  - coefficient of x3⁴ = `S_3(x1,x2,t)²`.

  I tested 100 random decompositions over p = 10007, 65537, 999983 and 1000003 with my own group
  law; S_4 vanishes on all 100. As a negative control, all 100 perturbed targets give nonzero
  values.
* **Case data** (`verify/extra.json`). For all 461 cases, p is prime, `2^k | p-1`, the curve is
  nonsingular and the target lies on it. All 113 planted triples have x_i in H, lie in Z, and sum
  to R under my group law.
* **Rational decompositions** (`verify/decomp.json`). I enumerated all of them by brute force with
  P_i in E(F_p) and x_i in H: 8331 x-triples over 461 cases. All 8331 lie in Z, and the per-case
  counts equal the producer's `rational_ordered` in 461/461 cases.
* **Corollary 2 hypotheses** (evaluated by direct computation, not taken from the producer). For
  every k >= 2 case, deg NF(g) = 6 exactly when c333 != 0 mod p. Also, max deg NF(x1·g) = 6, so
  (a) holds. For deg-6 cases at k >= 3, g is nonzero on at least 8 distinct slices {x1 = c}, so (b)
  holds. x1·g is never proportional to g on H³. The hypotheses hold in 105 / 107 / 106 cases at
  k = 2 / 3 / 4, which equals the producer's `carry_lower_bound_7_applies`.
* **Refutation attempt off the grid** (`verify/newprimes.jsonl`). I used new primes, random
  curves and random rational targets with c333 != 0:
  - k = 3: p = 41, 73, 89, 137, 233 and 409, four instances each;
  - k = 4: p = 113, 193 and 241.

  All 27 have exact LFD = 7. None is below the proven bound of 7, and the value 7 also persists
  off the grid.
* **Ordinary ring** (`verify/ordinary.jsonl`). I built the chains and NF inputs myself as explicit
  polynomials in the full monomial space and ran them through `lfd_engine.exact_lfd`. It agrees in
  38/38 cases (k = 1, 2), including all four LFD-0 cases it sampled at k = 1, LFD 3 and 6 cases,
  and 10 split cases.
* **Producer commands, run on a copy in the scratchpad.**
  - `ic_bounds_m3_symbolic.py` and `ic_bounds_m3_summarize.py` reproduce `symbolic.json`,
    `summary.json`, `rows.jsonl` and `table.md` byte for byte (hashes e7450717, be165905,
    2fed8fa1, 6f00ca92).
  - `ic_bounds_m3_validate.py` (734 s) reproduces `validation.json` byte for byte (bcfd3c79):
    140/140 direct and 122/122 split.
  - `ic_bounds_m3_selfcheck.py` (147 s) reproduces `selfcheck.json` byte for byte (2a9b78ab).
  - `ic_bounds_m3_run.py generate` reproduces `cases.jsonl` (342cc287) and `skips.json`
    (bee99170).
  - `ic_bounds_m3_run.py run --k 1` from empty results gives 333 lines. 322 are identical to the
    stored `results.jsonl` apart from timing. The 11 k = 1 split cases with `I_field != I_free`
    differ: the stored lines came from an older code path, with the field ideal "capped" after
    closures up to D = 9, while the current code records the Remark 3 `>= p` certificate. The
    free results are identical, and `summarize` turns both into the same `rows.jsonl`. This is a
    provenance note only: no number changes.

## Review of the written proofs

* **Local-to-global criterion** (README): correct as written. It needs only that `A = R/I` is
  finite dimensional, that `M_<D` spans A, and that `dim V_D = |M_<=D| - dim A`. The induction
  `V_{D'} = I cap M_<=D'` for `D' >= D` is valid, and "LFD = last strict fall" follows. The
  validator engine's docstring uses the same argument.
* **Proposition 1** (carry lemma): correct.
  - The chains are a Groebner basis for a degree order, so every `c` in C has a standard
    representation. Its top form therefore lies in the monomial ideal `(u_{i,j}^2)`.
  - `g`'s top form is squarefree, so `W = (C cap M_<=d) + F_p g` satisfies `W cap M_<d = C cap M_<d`.
    So W is closed under the rule, contains the inputs of degree at most d, and contains `V_d`.
  - `h = NF(u g)` lies in W only if `u g ≡ λ g` mod C. C is radical, so this means g is supported
    in one level set (λ = 0 gives g = 0 on H^m).
  - `V_e ⊆ V_{e+1}` then gives `LFD >= d+1`.
* **Corollary 2**: correct.
  - For k >= 2 the only exponent in 0..4 with popcount 2 is 3, so (3,3,3) is the unique possible
    degree-6 mask.
  - Multiplying by u_{1,0} maps the block-1 exponents 0,1,2,3,4 to 1,2,3,4,5, so the degree rises
    by at most 1, and (3,3,3) drops to weight 5. This is (a).
  - Each slice x1 ≠ t carries at most 2·2^k + 4·2^k = 6·2^k < 4^k zeros when k >= 3. There are
    at least 2^k - 1 >= 7 such slices. This is (b).
  - The hypothesis `p > 5` is automatic when 2^k | p-1 and k >= 3.

  The proof is labelled "proof", not kernel, which is accurate.
* **Remark 3**: correct. `z^p - z` has degree p, so `V_D(F_field) = V_D(F_free)` for D < p. If
  `dim A_free > dim A_field`, then `z^p - z ∉ I_free`. Its remainder `rho` modulo `I_free` has
  degree < D0 < p, lies in `I_field \ I_free`, and so lies in no `V_e` with e < p. If the
  dimensions are equal, the ideals and the LFDs coincide, because the free certificate degree is
  below p.

## Claim-by-claim notes

1. **Direct histogram (exact): confirmed.** I reproduce k = 1, 2, 3 in full and k = 4 on 26
   cases, all 7. The closure_check provenance is present on every profile row of the result files
   (`check_closed_exclusion_certificate` or unit for direct k <= 3; `all-low-shifts-reduce`,
   `rank-equals-dim` or unit for k = 4).
2. **k = 2 attains mk+1 = 7 (exact): confirmed.** All 105 deg-6 cases at k = 2 recomputed give 7.
   Corollary 2's k = 2 instance check (no proportionality) holds in all of them.
3. **Proposition 1 (proof): confirmed** (proof checked above).
4. **Corollary 2 (proof): confirmed.** The proof is valid and the identities were re-derived. The
   hypotheses hold in exactly 107 and 106 cases, no exact value contradicts the bound, and 27
   off-grid instances all give 7.
5. **Exceptional locus (exact): confirmed.** All deg-5 cases were recomputed (11, 10, 11). Each
   satisfies c333 ≡ 0, and their LFDs are 6, 7, 7.
6. **Empirical LFD = deg NF + 1: confirmed, with a miscount.** The deg-6 cases at k = 2, 3, 4
   number 105 + 107 + 106 = **318**, not 316. The pattern holds in all of them: all were
   recomputed at k = 2, 3, and 15 deg-6 cases at k = 4, plus the rows. It also holds in the 100
   deg-3 cases at k = 1. The sentence that the degree-7 certificate "must stop working at some
   depth" is a heuristic, not a result. It assumes |Z| grows like 8^k/p, which fails for unit
   ideals (27/117 at k = 4 have |Z| = 0) and for large p. The claim correctly says this is not a
   bound beyond 7.
7. **Split I_free histogram (exact): confirmed** on the sample: k = 1, 2 in full; k = 3 on 29
   cases, including all 6 LFD-5 cases (all y_R = 0); k = 4 on 5 cases, all 6. Coverage of the 117
   k = 4 split cases across the three result files is complete, with no duplicates.
8. **Field-equation split (proof): confirmed.** Remark 3 is valid, and the dimensions and `>= p`
   statuses match on all 261 recomputed split cases. The counts 100/87/80/69 and 11/29/37/48
   recount correctly from `rows.jsonl`.
9. **m2 reference and paired gaps: overstated, but only in wording.**
   - m2 is 3 at k = 1 and 4 at k = 2, 3, 4 in all 461 cases (all recomputed).
   - Direct minus m2 is +3 at k = 2 (105 cases), k = 3 (117) and k = 4 (117), as stated.
   - The summary "split minus m2 is +2 at k = 3, 4" is false for 6 of 117 cases at k = 3, which
     are +1. The README table states this correctly.
10. **Validator agreement (exact): confirmed.** Rerunning the producer's script reproduces it byte
    for byte. My own 38 ordinary-ring checks, with inputs I built myself, also agree.
11. **S_4 correctness (exact): confirmed** by my own resultant, the identities and the numeric
    decompositions. In addition, all 8331 sweep decompositions lie in Z. One clarification: the
    README's "125 terms generically" counts S_4(x1,x2,x3,t); the full S_4 in four variables has
    439 x-monomials.

## Not checked or limits

* At k = 4, only a sample was recomputed: 26/117 direct and 5/117 split. Each split case cost
  4 to 5 minutes on this loaded host, so the rest of k = 4 rests on the producer's runs, the rows
  and the structural checks.
* k = 3 split: 29/117 recomputed.
* Nothing here is Lean-checked; the task did not call for it.
* The heavy k = 5 runs were not attempted. The monomial counts they quote are correct:
  `|M_<=7| = 16384` for r = 15, and 17548 mixed monomials of degree at most 6.
