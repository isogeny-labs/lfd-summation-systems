# Preregistration: growth of the last fall degree with the depth for m >= 3 summands

Written 2026-11-20, before any closure of this study was run. Every number that this
study will produce is `[pcc:unverified]` (no pcc claim, no S3 record, one host).

## Question

For the Petit-Kosters-Messeng two-power subgroup factor base `H = mu_{2^k}` of `F_p^*`
and `m >= 3` summands, how does the last fall degree (LFD; on a unit ideal equal to the
Caminata-Gorla solving degree for every degree-compatible order, paper Section "The
solving degree") grow with the depth `k`? The paper proves only `LFD >= 7` for `m = 3`
at every `k >= 3` when `5t^3 + 3at + 6b != 0` (Theorem `thm:carry`) and `LFD <= mk + 1`
(Corollary `cor:allm`). Its data: exact `7` at `k = 2..5` on twelve instances (proper
and unit) and on the 234 seeded cases at `k = 3, 4`; at `k = 6` one proper instance
(`p7681-r`, `|Z| = 48`) has `V_7` closed and deficient, so `LFD >= 8` there; nothing
beyond (`experiments/ic-bounds-m3-depth`, paper Section `sec:msummands`). For `m = 4`:
exact `13 = deg g + 1` at `k = 4` on six instances; nothing at `k >= 5`.

## Presentation, and why

The **direct normalized chain presentation**
`F = {NF(S_{m+1}(x_1, ..., x_m, t))} ∪ G_{m,k}`, with `G_{m,k}` the `mk` chain equations
`u_{i,j}^2 - u_{i,j+1}`, `u_{i,k-1}^2 - 1`, `x_i = u_{i,0}`, and `NF` the chain normal
form (squarefree in the `mk` chain variables). This is the system to which the paper's
`m >= 3` statements apply: the carry lower bound (`thm:carry`), the dimension bound
`mk + 1` (`cor:allm`), Table `tab:msum` and the solving-degree statement
(`thm:sd` item 3). The split presentation (free auxiliary abscissae, `S_3` chain) has
different values (5 or 6 for `m = 3` at `k <= 4`, `tab:splitm`) and other bounds, and
PKM's System (4) is the raw split tower with auxiliaries; neither carries the paper's
`m >= 3` bounds, so neither is used here. Closures are computed modulo the chains; by
the chain bridge (`thm:bridge`) this is the ordinary closure for `D >= 2`.

`S_4 = Res_X(S_3(x_1, x_2, X), S_3(x_3, t, X))` and
`S_5 = Res_X(S_3(x_1, x_2, X), S_4(x_3, x_4, t, X))`. The generator `g = NF(S_{m+1})` is
computed as in the paper (Section `sec:msummands`, "We computed g exactly"): the values
of `S_{m+1}` on `H^m` from the resultant recursion with numbers (formal-degree Sylvester
determinants, so the values are those of the polynomial), then the inverse discrete
Fourier transform over `(Z/2^k)^m`; the coefficient of `x^e` (`e_i < 2^k`) is the
coefficient of the chain monomial whose block-`i` bits are the binary digits of `e_i`.

## Engine

The extracted second engine of the paper, `second-engine/engine`
(`hky_fp_tower_fast2`, "hky3", files unchanged), driven by a new Rust example
(`experiments/msum-growth/rust`, this study) that builds the ring of `mk` Kummer chain
variables (`y_{i,j}^2 = y_{i,j+1}`, top `y_{i,k-1}^2 = 1`, no free variable), the
generator above, and calls `closure_fast2` with the swap of the first two blocks as the
candidate symmetry (the engine verifies it on the generator before use; `S_{m+1}` is
symmetric in its summands) and `verify_closed = on`.

Checks before any new cell is read (reproductions of recorded values, not new data):
1. Generator: the values vanish at `x(P_1) + ... ` decompositions (planted points),
   `deg g = 6` for `m = 3` and the `(3,3,3)` coefficient is proportional to
   `5t^3 + 3at + 6b` (ratio to the `-8(5t^3+3at+6b)` of `thm:carry` constant across
   targets), `deg g = 12` for `m = 4` at `k >= 4`.
2. Engine on recorded cells: the profiles of `experiments/ic-bounds-m3-depth/closures.jsonl`
   for `p7681-r` at `k = 5` (`D = 7`) and `k = 6` (`D = 6, 7`; recorded
   `(0,0,0,0,0,3,2375,19958)`, closed, no `1`) must be reproduced exactly, and one
   `m = 4, k = 4` instance of `experiments/ic-bounds-msummands` at `D = 13` if it fits.
   A mismatch stops the study until explained.

## Instances

Primes with `2^k | p - 1` and `p > 2^{k+1}`:
* unit-ideal candidates (`p` large against `N^m`): `2013265921 = 15·2^27+1`,
  `998244353 = 119·2^23+1`, `469762049 = 7·2^26+1`;
* proper instances: the recorded `p7681-r` (`p = 7681`, `a = 1876`, `b = 2678`,
  `t = 6223`, `k <= 9`) and one seeded instance over `65537`.

Per prime one curve and target from a seeded generator (seed `msum-<p>-<i>`, `i = 0, 1, ...`):
random `a, b`, nonsingular; target a random rational point `(t, y)`; rejected if
`5t^3 + 3at + 6b = 0` (`m = 3`) so that `deg g = 6`. The same `(p, a, b, t)` is used at
every depth. `|Z|` (zeros of `g` on `H^m`, ordered tuples) is counted exactly at every
`(instance, k)`; the instance is a **unit** instance at that `k` iff `|Z| = 0`, else
**proper** and reported separately. At least two unit instances per `k` (different
primes), as far as resources allow.

## Quantity measured

* Unit ideals: `LFD = least D with 1 ∈ V_D`. For every smaller `D` run, the closure
  ended without `1` and passed the engine's closedness check (`verify_closed`: every
  generator and every variable times every final basis row of degree `< D` reduces to
  zero), so `1 ∉ V_D` is certified for that `D` in the sense of the paper's lower-bound
  certificate. `V_D = 0` for `D < deg g` and `V_{deg g} = span(g)` need no run.
* The profile `p_j = dim(V_D ∩ M_{<=j})`, the column count `|M_{<=D}|` (squarefree
  monomials of degree `<= D` in `mk` variables), wall time, peak RSS (`VmHWM`), the
  engine counters (`muladds`, rows reduced), and the SHA-256 of the binary.
* Proper ideals: the profile, `|Z|`, `rank ev(M_{<=j} -> F_p^Z)` for every `j`, so that
  `dim(I ∩ M_{<=j}) = |M_{<=j}| - rank ev(M_{<=j})`; the value is exact when the
  local-to-global criterion holds at `D` (`rank ev(M_{<D}) = |Z|` and
  `dim V_D = |M_{<=D}| - |Z|`), and then is the last strict fall; a closed `V_D` with
  `dim V_D < dim(I ∩ M_{<=D})` gives `LFD > D`.
* Caps: a single closure stops at about 2 h wall or about 12 GB RSS; a stopped closure
  is reported as censored and bounds nothing.

## Predictions (written now)

* **P1** (`m = 3`, `k = 2..5`, unit): `LFD = 7` on every instance. Falsified by any other
  value.
* **P2** (`m = 3`, `k = 6`, unit): `1 ∉ V_7`, i.e. `LFD >= 8`, as on the proper
  `p7681-r`. Falsified if `1 ∈ V_7` on any unit instance.
* **P3** (`m = 3`, `k = 6`, unit): `LFD = 8` exactly (`1 ∈ V_8`). Falsified if a closed
  `V_8` misses `1`.
* **P4** (growth law, the main hypothesis): the three-summand value is the two-summand
  value plus three, `LFD(3, k) = floor((k + 14)/4) + 3`: `7` at `k <= 5`, `8` at
  `k = 6..9`, `9` at `k = 10..13`; i.e. it grows with `k`, by about one per four
  depths, with the same step depths as `m = 2`. Consequence testable here: `LFD(3, 7) = 8`
  (`1 ∉ V_7`, `1 ∈ V_8` at `k = 7`). Falsified by `LFD(3, 7) >= 9` (faster growth, e.g.
  the linear alternative `LFD = k + 2`, which predicts `8` at `k = 6` and `9` at `k = 7`)
  or by `LFD(3, k) = 7` at some `k >= 6`. `k = 6` alone cannot separate P4 from the
  linear alternative; `k = 7` can, and `k = 10` would test the second step.
* **P5** (proper instances): the proper instances have the same value as the unit ones
  at the same `k` (in particular `p7681-r` has `LFD = 8` at `k = 6`). Falsified by a
  different exact value.
* **P6** (`m = 4`): `LFD = 13 = deg g + 1` at `k = 4` (reproduction) and at `k = 5`
  (no growth yet; by analogy with `m = 3`, where the first step comes at `k = 2m = 6`,
  the first step for `m = 4` is expected only at `k >= 8`). Falsified by `LFD(4, 5) >= 14`.
  This cell may not fit the host (`|M_{<=13}|` about `9.9·10^5` columns at 20
  variables); if it does not run, P6 at `k = 5` stays untested.

## Reading rules

* An exact value needs a closure with `1` (unit) or the local-to-global criterion
  (proper) at `D`, and certified closures without `1` at every `D' < D` down to
  `deg g + 1` (unit).
* A censored closure gives no bound. "`>= D`" is reported only from a certified closed
  closure without `1` at `D - 1` on a unit ideal, or a closed deficient closure on a
  proper ideal.
* Costs reported are those of this engine on this host (4 vCPU, 15 GB); the
  extrapolation to a 32 vCPU / 256 GB host is an estimate, labelled as such.
