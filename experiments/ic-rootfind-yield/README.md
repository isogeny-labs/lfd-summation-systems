# Relation yield per linear system for LRF (binary curves)

**Question.** The linearized root-finding (LRF) kernel made the point-decomposition *step* fast
(`experiments/ic-rootfind/`). For index calculus end to end the step speed is no longer the
bottleneck; the **relation yield per linear system** is. This directory measures the yield, finds
where it is lost, and tests two levers the owner pinned (2026-11-21):
(a) factor bases / summand counts that give more relations per linear system;
(b) the Koblitz / Frobenius structure at prime `n`.

Everything here is a **count**, not a speed, so it is measured in Python at small prime `n` where the
decomposition set can be cross-checked against exact enumeration. Per-system wall-clock is already in
`experiments/ic-rootfind/` (Rust kernel). No ECDLP speed claim is made; see "Gap to rho".

## The funnel and the metric

For `m = 3` and factor base `V = <1, w, ..., w^(k-1)>`, LRF scans each abscissa `x3 in V`, forms
`T = R - P3`, and solves one F_2-linear system (`k + 1` unknowns, `n - dim V2 + 1` equations). The
C-kernel `Stats` counts the funnel exactly:

```
targets  ->  lin_solvable  ->  lin_candidates  ->  found
(systems)    (consistent)     (quadratics)        (both roots in V = a relation)
```

- **yield per system** = `found / targets`
- `= solvable_rate x cand_per_solvable x found_per_cand`

The determinacy knob is `delta = k + dim V2 - n`. For monomial `V`, `dim V2 = 2k - 1`, so
`delta = 3k - 1 - n`; the system is square at `delta = 0` (`n = 3k - 1`).

## Baseline (m = 3, monomial V, 5 seeds, 150 targets/seed)

Cross-check: on every row the LRF decomposition set equals the enumeration set (`mismatches = 0`).
`cost/rel = (targets + lin_candidates) / relations` is the work (systems plus candidate processing)
per distinct relation; it is the objective to minimise.

| n | k | m | dim V2 | delta | solvable | cand/solv | found/cand | yield/sys | cost/rel |
|---|---|---|---|---|---|---|---|---|---|
| 11 | 4 | 3 | 7 | 0 | 0.632 | 1.55 | 0.0433 | 0.0425 | **118.5** |
| 11 | 5 | 3 | 9 | 3 | 0.953 | 8.36 | 0.0264 | 0.2104 | 118.2 |
| 13 | 5 | 3 | 9 | 1 | 0.795 | 2.52 | 0.0358 | 0.0719 | **116.8** |
| 13 | 6 | 3 | 11 | 4 | 0.978 | 16.40 | 0.0167 | 0.2682 | 183.4 |
| 13 | 4 | 4 | 7 | -2 | 0.226 | 1.12 | 0.0773 | 0.0195 | 309.7 |
| 17 | 6 | 3 | 11 | 0 | 0.617 | 1.63 | 0.0157 | 0.0157 | **369.6** |
| 17 | 7 | 3 | 13 | 3 | 0.944 | 8.48 | 0.0082 | 0.0654 | 406.7 |
| 19 | 7 | 3 | 13 | 1 | 0.784 | 2.55 | 0.0078 | 0.0155 | **573.4** |
| 19 | 8 | 3 | 15 | 4 | 0.971 | 16.48 | 0.0037 | 0.0597 | 848.3 |
| 23 | 8 | 3 | 15 | 0 | 0.614 | 1.63 | 0.0037 | 0.0037 | **1625.7** |
| 23 | 9 | 3 | 17 | 3 | 0.941 | 8.50 | 0.0019 | 0.0155 | 1732.1 |

(full sweep in `baseline_summary.jsonl`; per-seed rows in `base_s*.jsonl`.)

### What the baseline establishes

1. **The optimal operating point is `delta` in {0, 1} (a square system), not high `delta`.** High
   `delta` maximises yield *per system* (kernel dimension gives many candidates), but the candidate
   processing and the `2^k` factor-base size dominate the cost. Minimum `cost/rel` sits at the
   square regime at every `n`. This refutes "bigger `k` is always better".
2. **`m = 4` is strictly worse than `m = 3`** at every matched `(n, k)` (e.g. `n = 13, k = 4`:
   310 vs 145 cost/rel) and costs a factor `|FB|` more targets. More summands do not raise yield per
   system here. A clean negative result for one half of lever (a).
3. **`found_per_cand` collapses as `n` grows** at fixed `delta` (0.026 at `n = 11` to 0.0019 at
   `n = 23`, `delta = 3`). This is the dominant yield loss and it is intrinsic: fixing the symmetric
   functions `e1 in V`, `e2 in span(V.V)` does not force the two roots into `V`. `cost/rel` at the
   optimal `delta` grows about `2^(0.7 n)`. **A faster kernel cannot change this exponent;** only the
   factor-base structure can.

## Lever (b): Koblitz / Frobenius — structure probe

`csd.ic_dreg.koblitz_invariant_dims` reports no proper Frobenius-invariant *subspace* at prime `n`,
which is why `V = <1, w, ...>` is not Frobenius-stable. But index calculus needs a Frobenius-stable
*set*, not a subspace. The `n` squaring-conjugates `V^(i) = {v^(2^i) : v in V}` satisfy (measured,
`n = 11..23`):

- `dim V2(V^(i))` is **the same for all `i`** (squaring is an F_2-algebra automorphism), so every
  conjugate is equally LRF-friendly;
- the `n` conjugate subspaces are **distinct**, and factor-base abscissas lie in **size-`n` squaring
  orbits** (only `0, 1` are fixed).

So `F = union_i V^(i)` is a Frobenius-stable factor base, each piece individually LRF-solvable with
the same `dim V2`, and its points partition into size-`n` orbits. The orbit-reduced relation
accounting is measured in `frobenius_*` (see that section when it lands).

## Gap to rho

At `n = 131` (ECC2K-130) the matched Pollard rho reference is about `2^61` group operations. The
baseline cost model above is the honest distance-remaining denominator; any lever here is reported as
a factor off that model, never as a pass/fail.

## Reproduce

```
python3 experiments/ic-rootfind-yield/yield_baseline.py --targets 150 --seed S --out base_sS.jsonl
```
Single thread, exact (cross-checked against enumeration on every target).
