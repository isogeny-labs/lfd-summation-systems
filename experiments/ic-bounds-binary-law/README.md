# What governs the LFD of two-summand binary/Koblitz systems

Status 2026-11-10. Systems: Weil descent of `S_3(x1, x2, x_R)` with `x1, x2` in an F_2-linear coordinate
space `V` of dimension `k` in `F_(2^n)`, Boolean ring (field equations), exact ordinary LFD. The data
are the existing certified values (`experiments/ic-curve-sweep`, k = 2..7; `experiments/ic-bounds-binary`,
k = 8..10), deduplicated to **2277 distinct proper-subspace systems** (`k < n`), plus the full-field inverse
family (`experiments/ic-binary-inverse`, k = n = 1..12), plus **80 new preregistered cases (60 distinct
systems, all certified)** and 4 subfield checks. No cost or speed claim is made.

Labels: **proof** (written out here), **exact** (certified computation named), **empirical** (a pattern
in the data, with counts), **conjecture**.

## 1. Structure of every system (proof; exact on all 2277)

`S_3 = (x1 x2)^2 + x_R^2 (x1 + x2)^2 + x_R x1 x2 + b` in characteristic 2.

* **S1 (proof).** Squaring is F_2-linear, and `x1 x2 = sum a_i b_j v_i v_j` is bilinear in the
  coordinates `a` of `x1` and `b` of `x2`. So every Boolean equation is *affine bilinear*: its only
  quadratic monomials are `a_i b_j`, and the system is symmetric under `a <-> b`. The quadratic parts lie
  in the space of symmetric bilinear forms, so their rank `m_eff <= min(n, k(k+1)/2)`. `a` (the curve
  coefficient) does not occur at all; the system depends on `(n, V, b, x_R)` only. (exact: `bilinear` is
  true for all 2277 systems in `structure.jsonl`.)
* **S2 (proof): the trace equation.** For `x_R != 0` the quadratic part is `phi(x1 x2)`,
  `phi(u) = u^2 + x_R u`, whose image is `x_R^2 {w^2 + w} = x_R^2 ker Tr`. Applying
  `lambda(y) = Tr(y / x_R^2)` to the field equation kills `phi(x1 x2)` and leaves the affine-linear equation
  `Tr(x1) + Tr(x2) = Tr(b / x_R^2)` in the span of the Boolean equations. Hence: if `V` is not inside
  `ker Tr`, the span contains one linear equation (`r_eff = 2k - 1`, `m_eff <= n - 1`); if `V` is inside
  `ker Tr` (e.g. the trace-zero hyperplane `{z^2 + z}`), the equation is the constant `Tr(b/x_R^2)`:
  either `1` is in the input span (unit, LFD <= 2, `Z` empty) or there is no linear constraint and `|Z|`
  roughly doubles. This is the mechanism behind "trace-zero hyperplane doubles or empties `|Z|`".
  (exact consistency: among 1439 affine non-const systems `r_eff = 2k - 1` in 1047, `2k` in 155, of which
  154 on Frobenius-invariant bases; the remaining 237 have `m_eff < n - 1` at small `k`.)
* **S3 (proof).** For the order-two target `x_R = 0` the system is `x1 x2 = sqrt(b)`: homogeneous bilinear
  plus a constant, no linear part. If `V` is a scaled subfield `F_(2^k)` it is the inverse-family ideal of
  degree `k`.

Classes used below: `const` (1 in the input span; 502 systems), `order-two` (336), `affine` (`x_R != 0`,
1439). `e = r_eff - m_eff` is the log2 of the expected `|Z|` (`r_eff` = variables left after the linear
equations, `m_eff` = rank of the quadratic parts; `structure.py`).

## 2. Candidate laws on the existing data (empirical; counts of systems)

Full tables: `laws.md`, `laws.json` (`analysis.py`). A "lookup" law fits the modal LFD in each cell of
its variables; singleton cells fit trivially and are listed.

| law (affine class, 1439 systems) | cells (singleton) | satisfy | violate |
|---|---|---|---|
| LFD = f(\|Z\|) | 162 (74) | 1061 | 378 |
| LFD = f(k, \|Z\|) | 339 (165) | 1268 | 171 |
| LFD = f(n − k) | 10 (0) | 800 | 639 |
| LFD = f(2k − n) | 18 (0) | 1214 | 225 |
| LFD = f(n, k) | 54 (0) | 1307 | 132 |
| LFD = f(e) | 15 (0) | 1208 | 231 |
| LFD = f(k + e) | 21 (0) | 1315 | 124 |
| LFD = f(k, e) | 53 (1) | 1335 | 104 |
| LFD = f(n, k, base kind, curve set) (most specific) | 177 (7) | 1355 | 84 |

The 84 violations of the most specific cell law are target-to-target variation inside one
`(n, k, base, curve)` cell: no tested variable removes them.

Closed forms (all 1775 non-const systems unless stated):

| closed-form law | satisfy | violate |
|---|---|---|
| A (frozen before the check): affine, LFD = 2 / 3 / 4 / 5 for e ≤ −1 / ≤ 3 / ≤ 8 / ≥ 9 | 1208 | 231 |
| B (frozen before the check): affine, LFD = 2 + [k+e ≥ 3] + [k+e ≥ 9] + [k+e ≥ 19] | 1315 | 124 (93 below, 31 above) |
| B on order-two systems | 200 / 336 | 136 |
| semi-regular crossover: least D with Σ_{i≤D} c_i(r_eff, m_eff) ≤ \|Z\| (Boolean series (1+t)^r/(1+t²)^m) | 35 | 1740 |
| LFD = Boolean semi-regular d_reg(r_eff, m_eff) | 598 | 1177 (1173 below) |
| k = 10 grid rule "LFD 5 iff \|Z\| > 470" (proper, k ≥ 8) | 489 / 489 | 0 |

Observations (empirical):

* **|Z| alone does not determine the LFD.** Unit ideals (`|Z| = 0`, not const) have LFD 4 at
  `(n, k) = (12, 7)` (3 systems), while proper systems with `|Z| = 4..14` at the same `k` have LFD 3; at
  `e = 2..3` the LFD is 3 for `k <= 5` and mostly 4 for `k >= 7` (table by `e` in `laws.json`).
* **The closure does not follow the semi-regular model.** The crossover law fits 35 of 1775; the
  semi-regular `d_reg(r_eff, m_eff)` is an upper bound in 1771 of 1775 (4 violations, all `k <= 3`) but
  exceeds the LFD in 1173. A random symmetric affine-bilinear model with the same shape
  (`bilinear_model.py`, 16-run pilot `bilinear-pilot.jsonl`) gives LFD 4 at `(k, m) = (6, 8)` and `(6, 9)`
  where the curve systems at `(n, k) = (9, 6)`, `(10, 6)` have 3 (18 of 24) and 3 (24 of 24): the field
  structure lowers the LFD below a generic bilinear system (pilot only, not a controlled comparison).
* **Pre-certificate closure dimensions depend on (n, k, class, r_eff, m_eff) and not on curve or base.**
  Grouping by `(n, k, r_eff, m_eff, class)`, the whole profile `dim V_D` for `D < D*` equals the group's
  modal profile in 1678 of 1775 systems (97 differ: 60 affine/random, 32 affine/invariant,
  5 order-two/invariant; 40 of the 187 groups are singletons). Example: every affine random-base system
  at `(12, 10)` has codims `(1, 21, 180, 845)` before the certificate (the `180` is exactly
  `|M_<=2| − (n − 1) − r`: `n − 1` quadratics plus `x_i * (trace equation)`). So the binary LFD *is* "the
  first degree where a closure whose dimension depends only on (n, k, class) catches the ideal" — but
  that catch-up degree is not predicted by any |Z| crossover of the recorded dimensions: at `(12, 7)` the
  generic codim is 81 at both `D = 2, 3` and the catch-up happens at 3 for some `|Z| = 4` systems and at 4
  for others with `|Z| = 0..10`.
* `D* − LFD` is 0 in 1476 and 1 in 8 of 1486 proper systems (2 have 3).

## 3. Preregistered check (exact)

`predictions.json` was written by `prereg.py` (2026-11-10T02:29:54Z) **before any closure was run**:
structural class, `e` and `|Z|` (field enumeration only), and predictions of the two top closed laws
A (e only) and B (k + e), plus the k = 10 grid rule A_Z. Seed 20261110, runner
`scripts/ic_bounds_binary_run.py` (commands in `predictions.json`), results `results-prereg.jsonl`,
scored by `score.py` → `score.json`. All 80 cases certified (59 affine distinct systems scored, 1 const,
order-two recorded unscored); every row's `system_hash` and `|Z|` equal the preregistered ones.

| cell | systems | e | k+e | A predicts | B predicts | measured LFD | A | B | A_Z |
|---|---|---|---|---|---|---|---|---|---|
| P1 n=10, k=9, trace-zero hyperplane | 11 affine (\|Z\| 482–499) | 9 | 18 | 5 | 4 | **5** (11/11) | 11 hit | **11 miss** | 11 hit |
| P2 n=19, k=10, random base | 12 (\|Z\| 0–6) | 1 | 11 | 3 | 4 | **4** (12/12, incl. a unit) | 12 miss | 12 hit | 12 hit |
| P3 n=16, k=9, random base | 12 (\|Z\| 0–8) | 2 | 11 | 3 | 4 | **4** (12/12, incl. 2 units) | 12 miss | 12 hit | 12 hit |
| P4 n=19, k=8, random base | 12 (\|Z\| 0–2) | −3 | 5 | 2 | 3 | **3** (12/12) | 12 miss | 12 hit | 12 miss |
| P5 control n=11, k=9 | 12 (\|Z\| 110–152) | 7 | 16 | 4 | 4 | 4 (12/12) | hit | hit | hit |

Totals on the discriminating cells P1–P4 (47 systems): **A 11 hit / 36 miss, B 36 hit / 11 miss.**
Const case (P1, `Tr(b/x_R^2) = 1`): LFD 2 as predicted by S2. Order-two (unscored): P1 4 (3 systems,
`|Z|` 241–262, certificate at 5), P2–P4 3, P5 4.

**Outcome: both preregistered single-variable laws are refuted.** B fails exactly where `|Z|` is large
(P1: 4 predicted, 5 measured at `k = 9`, so the LFD-5 regime is not tied to `k = 10` or to `k + e >= 19`);
A fails wherever `|Z|` is small and `k` is large (P2–P4: unit and near-unit systems at `k = 8..10` have
LFD 3–4, not 2–3).

**Post-hoc (conjecture, not preregistered):** the data are consistent with two separate mechanisms,
`LFD = max(2 + [k+e >= 3] + [k+e >= 9],  A(e))`: a `k`-driven floor at small `|Z|` and an `e`-driven
(`|Z|`-driven) rise at large `|Z|`. It matches 1317 of 1439 existing affine systems (103 below, 19 above
the prediction) and all 59 new ones — the new ones were used to choose it, so it needs its own
preregistered check (next step: `k = 11, n = 14` where it predicts 4 and B predicts 5; `k = 8, n = 9`
trace-zero, `e = 8`, predicts 4; a `k = 9` trace-zero system at `n = 10` with `|Z|` forced below 470
cannot be made, `|Z|` concentrates near `2^9`).

## 4. Bounds and what a proof would need

The local criterion in Lean (`LocalReductionBound.lfd_bound_of_local`): LFD ≤ D if (i)
`dim V_D = |M_<=D| − |Z|` and (ii) `M_<D` evaluates onto `F_2^Z`.

* **Proof (condition (ii) from |Z| alone).** For any `Z ⊂ F_2^r` with `|Z| = N >= 1`, `M_<=floor(log2 N)`
  evaluates onto `F_2^Z`. Take a degree-compatible order; the leading monomials of `I(Z)` in the Boolean
  ring are closed under multiplication by variables (for `f ∈ I(Z)`, `x_i ∤ LM(f)`, `LM(x_i f) = x_i LM(f)`
  because terms containing `x_i` stay put and the others move up to `x_i t < x_i LM(f)`), so the `N`
  standard monomials form a down-closed family of subsets of `[r]`; a down-closed family containing a set
  of size `s` contains its `2^s` subsets, so `s <= log2 N`. The standard monomials span `R/I(Z) = F_2^Z`.
  So (ii) holds at every `D >= floor(log2 |Z|) + 1`, and the whole difficulty is (i).
* **No bound in |Z| alone exists (exact counterexamples).** LFD 4 with `|Z| = 0` (non-const units at
  `(12, 7)`; P2, P3); LFD 3 with `|Z| = 1` (28 proper systems violate `LFD <= max(2, floor(log2|Z|) + 2)`);
  LFD 4 with `|Z| = 2` (`k = 7`). Any valid bound must grow with `k` (or with the inputs) at fixed `|Z|`.
  So there is no provable bound in terms of the solution count that sits far below `min(k, l) + 2` for
  small `|Z|` without a `k`-dependent term, and condition (i) is a degree-of-regularity statement about the
  specific inputs, not about the point set.
* **"LFD ≤ 4 for |Z| ≤ 470" is false on proper subspaces (exact).** `subfield_check.py`: `V = F_(2^8)`
  inside `F_(2^16)` (proper, `k = 8 < n = 16`), `b = 1`, target `(0, 1)`: `|Z| = 255`, certified LFD 5
  (falls 3, 5), for `a = 0` and `a = 1`; `V = F_(2^6) ⊂ F_(2^12)`: `|Z| = 63`, LFD 4. Both equal the written
  inverse-family value `k/2 + 1` (S3), which grows without bound in `k` with `|Z| = 2^k − 1`. The empirical
  statements "LFD ≤ 4 whenever `|Z| <= 470`" (1458 of 1458 proper grid systems) and "LFD ≤ 4 whenever
  `k + e <= 18`" (1747 of 1747 non-const) hold only on the grid, which contains no scaled subfield of
  degree ≥ 8; the full-field inverse family violates both (n = 8: `|Z| = 255`, LFD 5).
* **What a proof of `LFD <= f(k, e)` for the affine class would need:** (ii) is free above
  `floor(log2|Z|) + 1`; (i) needs an explicit spanning set of `I(Z) ∩ M_<=D` inside the closure. The
  measured fact that the pre-certificate profile is the same for every curve and base with the same
  `(n, k, r_eff, m_eff)` suggests these spanning elements come from field identities (S1–S2 structure:
  `phi`, the trace equation, products `x1^(2^i) x2^(2^j)` that are all bilinear), as in the inverse-family
  proof, not from genericity. No such argument is written here. The proved bounds remain the kernel
  `k + l + 1` and the written `min(k, l) + 2`; empirically every system has LFD ≤ 5 at `k ≤ 10`.

## 5. Unexplained / undecidable at these sizes

* The within-cell variation (84 systems with identical `(n, k, base, curve)` but different LFD; e.g.
  `(12, 7)`, `|Z| = 4`: 7 systems LFD 3, 5 systems LFD 4) is not explained by any recorded variable.
* The order-two class has its own behaviour (LFD 3 up to `|Z|` ≈ 45, mostly 4 above; subfield cases follow
  `k/2 + 1`); no closed law is tested for it.
* Whether the LFD-5 threshold is in `e` (|Z|) alone at every `k`, and where LFD 6 first occurs for proper
  non-subfield spaces, needs `k >= 11` at `D = 5..6` (35,443+ columns, ≥ 3 GB per system): not run.
* No law explains why the curve systems sit below a random bilinear model; the pilot is 16 runs.

## Files and reproduction

`build_table.py` → `table.jsonl` (2481 case rows / 2289 systems incl. inverse) · `structure.py` →
`structure.jsonl` (equations regenerated for k = 8..10 with stored hashes checked) · `analysis.py` →
`laws.json`, `laws.md` · `prereg.py` → `predictions.json` (write-once) · runner → `results-prereg.jsonl` ·
`score.py` → `score.json` · `subfield_check.py` → `subfield-check.jsonl` · `bilinear_model.py` →
`bilinear-pilot.jsonl` · `plot.py` → `lfd_vs_kplus_e.{svg,png}`.
All runs with `ulimit -v 1500000`, one thread; peak 446 MB (P1, `D* = 5`, 66–95 s per system); no lock used.
