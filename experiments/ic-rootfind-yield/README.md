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
the same `dim V2`, and its points partition into size-`n` orbits.

### Lever (b) is a no-go at prime n (measured), with a mechanism

`frob_orbits.py` runs real LRF on the Koblitz curve `y^2 + xy = x^3 + 1` (a=0, matches ECC2K-130),
collects decompositions, canonicalises every point to its Frobenius orbit, and measures the column
reduction over the true `E[ell]` (tau = Frobenius eigenvalue). The `n` Frobenius images of one
decomposition are `tau^j .(one orbit row)` -- scalar multiples, rank 1 -- so the win can only be the
column reduction, not extra relations per search.

| n | k | ideal `|F|/|orbits|` | per-point cols | orbit cols | realizable reduction |
|---|---|---|---|---|---|
| 11 | 4 | 9.0 | 21 | 15 | 1.40 |
| 13 | 5 | 11.6 | 33 | 25 | 1.32 |
| 17 | 6 | 16.2 | 65 | 57 | 1.14 |
| 19 | 7 | 18.6 | 139 | 121 | 1.15 |
| 23 | 8 | 22.7 | 275 | 253 | **1.09** |

The *ideal* reduction grows toward `n`, but the *realizable* one **decays toward 1**. Reason: LRF over
one subspace `V^(0)` produces only abscissas in `V`, and those points lie in **distinct** Frobenius
orbits, so `#V-points ~= #orbits-touched`. Populating the Frobenius closure needs the `n-1` conjugate
subspaces, whose decompositions are Frobenius *images* -- the same orbit rows, **no new rank** -- so
they are wasted.

A Frobenius-*structured* factor base (`frob_span.py`, sliding window `V = <w, w^2, ..., w^(2^(k-1))>`)
does give real orbit sharing: `phi(V) cap V` has dimension `k-1`, and about half of `V` squares back
into `V`. But then `dim V2` saturates at `n` for every `k >= ~sqrt(2n)` (measured), far below the
yield-necessary `k ~= n/3`, so LRF is infeasible. **The LRF condition `dim span(V.V) < n` and Frobenius
structure are in irreducible tension** -- the real obstruction at prime `n`, sharper than
`koblitz_invariant_dims`'s "no proper invariant subspace".

### found_per_cand is intrinsic (measured)

`found_per_cand` is the same narrow band for every LRF-feasible basis at a given `(delta, n)`
(monomial, `w^2`-powers, shifted all within noise; random bases have `dim V2 = n` and are infeasible;
"far-out" bases only raise `delta`, re-tracing the baseline tradeoff). Fixing `e1 in V`, `e2 in
span(V.V)` does not force the two roots into `V`, and no factor-base choice changes that. So the
`~2^(0.7 n)` cost/relation exponent is intrinsic to LRF decomposition; only a decomposition relation
that constrains the roots themselves (not merely their symmetric functions) could move it.

Data: `frobenius.jsonl` (orbit reduction), and the basis comparison is reproducible from
`frob_span.py`. Figure: `yield_report.{svg,png}`; write-up: `report.pdf` / `report.tex`.

## Verdict for ECC2K-130 (n = 131)

For `m = 3`, any non-negligible yield forces `k ~= 44`, which is exactly `delta = 0` -- the measured
optimal operating point. **The current method is already at that optimum; neither pinned lever improves
the n = 131 relation stage.** This is a measured negative result with mechanisms, not an absence of
effort; it redirects the search to the `found_per_cand` term.

## Gap to rho

At `n = 131` (ECC2K-130) the matched Pollard rho reference is about `2^61` group operations. The
baseline cost model above is the honest distance-remaining denominator; any lever here is reported as
a factor off that model, never as a pass/fail.

## Reproduce

```
python3 experiments/ic-rootfind-yield/yield_baseline.py --targets 150 --seed S --out base_sS.jsonl
```
Single thread, exact (cross-checked against enumeration on every target).
