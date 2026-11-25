# Relation yield per linear system for LRF (binary curves)

**Interpretation correction (2026-11-25):** high ANF degree of the projected
candidate-success indicator does not prove that valid candidates require
exhaustive kernel search. The exact quadratic root-coordinate lift and its
cost gate are in [the new screen](../quadratic-root-lift-screen-20261125/SCREEN.md).
The measured counts below are unchanged; statements about optimality and
prime-degree scaling are limited to the tested configurations.

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

1. **The lowest measured cost per relation occurs at `delta` in {0, 1} (a square system).** High
   `delta` maximises yield *per system* (kernel dimension gives many candidates), but the candidate
   processing and the `2^k` factor-base size dominate the cost. Minimum `cost/rel` sits at the
   square regime at every `n`. This refutes "bigger `k` is always better".
2. **`m = 4` is worse than `m = 3` in the matched cells tested** (e.g. `n = 13, k = 4`:
   310 vs 145 cost/rel) and costs a factor `|FB|` more targets. More summands do not raise yield per
   system here. A clean negative result for one half of lever (a).
3. **`found_per_cand` falls across the tested `n` values** at fixed `delta` (0.026 at `n = 11` to
   0.0019 at `n = 23`, `delta = 3`). This is the dominant measured yield loss: fixing the symmetric
   functions `e1 in V`, `e2 in span(V.V)` does not force the two roots into `V`. `cost/rel` at the
   tested `delta` fits about `2^(n/3)` (fit `2^(0.333 n)`, consistent with the heuristic
   systems/relation = `2^(n-2k)` at `delta=0`). A faster residual solver alone does not change
   the relation count, but its full cost still needs measurement.

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

### The tested Frobenius variant did not improve the measured prime-degree cells

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

The *ideal* reduction grows with `n`, but the *measured* reduction approaches
1 across these cells. LRF over one subspace `V^(0)` produces only abscissas
in `V`; on these fixtures they occupy almost distinct Frobenius orbits.
The Frobenius images of any one decomposition are scalar multiples of its
orbit row, so repeating those images alone adds no rank.

A Frobenius-*structured* factor base (`frob_span.py`, sliding window `V = <w, w^2, ..., w^(2^(k-1))>`)
does give real orbit sharing: `phi(V) cap V` has dimension `k-1`, and about half of `V` squares back
into `V`. In the measured sliding-window cases, `dim V2` saturates at
`n` near `k ~ sqrt(2n)`, making those bases infeasible for this LRF
implementation at the tested yield range. This is a concrete tradeoff for
the tested construction, not a theorem excluding every Frobenius-stable
factor-base set.

### Root-membership loss in the tested bases

`found_per_cand` is in the same narrow band for the LRF-feasible bases tested at a given `(delta, n)`
(monomial, `w^2`-powers, shifted all within noise; random bases have `dim V2 = n` and are infeasible;
"far-out" bases only raise `delta`, re-tracing the baseline tradeoff). Fixing `e1 in V`, `e2 in
span(V.V)` does not force the two roots into `V`. These measurements do not establish
a lower bound for every basis or solver. Constraining roots directly is a
separate formulation whose full cost needs testing.

Data: `frobenius.jsonl` (orbit reduction), and the basis comparison is reproducible from
`frob_span.py`. Figure: `yield_report.{svg,png}`; write-up: `report.pdf` / `report.tex`.

### Projected success-indicator structure

At high `delta` the 2-summand system has a `delta`-dim affine solution space (`2^delta` candidates);
only the ones whose quadratic has both roots in `V` are relations. If that "good" set were an
algebraically structured subset of the kernel (low `F_2` ANF-degree / affine), it could be enumerated
output-sensitively and we could run at large `k` (relation density `~1` per solvable system) without the
`2^delta` candidate blow-up. `goodset_structure.py` builds the Boolean good-indicator over the `delta`
kernel coordinates and measures its ANF degree and density.

| n | k | delta | ANF degree (mean/max) | density | good set |
|---|---|---|---|---|---|
| 23 | 10 | 6 | 5.98 / 6 | 0.0160 ~ 2^-6 | singleton |
| 23 | 11 | 9 | 8.84 / 9 | 0.0023 ~ 2^-9 | singleton |
| 29 | 12 | 6 | 6.00 / 6 | 0.0156 = 2^-6 | singleton |
| 29 | 13 | 9 | 8.96 / 9 | 0.0020 ~ 2^-9 | singleton |
| 29 | 14 | 12 | 11.88 / 12 | 0.00028 ~ 2^-12 | singleton |

The indicator has near-maximal ANF degree and density about `2^(-delta)`
on these cells: there is roughly one good candidate per solvable system.
This describes the projected Boolean function, not the complexity of finding
its support. Even a unique solution of a linear system has a degree-`delta`
indicator. The exact quadratic root-coordinate lift in the linked screen
shows why these data do not rule out a solver that avoids enumerating all
`2^delta` candidates. Whether such a solver beats the native scan, including
its extra variables and setup, remains unmeasured.

## Verdict for ECC2K-130 (n = 131)

For `m = 3`, the simple yield model places the useful range near `k ~= 44`,
which is near `delta = 0` for a monomial base at `n = 131`. The measurements
at smaller fields motivate that operating point; they do not verify an
`n = 131` optimum or rule out a different root-membership solver.

## Gap to rho

At `n = 131` (ECC2K-130) the matched Pollard rho reference is about `2^61` group operations. The
baseline cost model above is the honest distance-remaining denominator; any lever here is reported as
a factor off that model, never as a pass/fail.

## Reproduce

```
python3 experiments/ic-rootfind-yield/yield_baseline.py --targets 150 --seed S --out base_sS.jsonl
```
Single thread, exact (cross-checked against enumeration on every target).
