# `hky.basis/v1`: export of a closed degree-D closure V_D

Status: specification, frozen 2026-11-18 (`wt/paper-engine`). Writer: `pkm_tower_pilot --instance … --export-basis <path>`
(engine `hky3`, also `hky2`). Readers: the round-trip checker `hky_basis_check` in this directory's harness, and the
degree-lift worker. Any change to the layout below is a new version (`hky.basis/v2`), never an edit of v1.

## What is exported

The F_p-subspace V_D ⊆ M_≤D of the two-summand prime chain system of the paper
(`curve-lfd-bounds.tex`, normalized two-summand chain presentation), written in **original monomial coordinates**:
no σ-orbit (summand-swap) coordinates appear in the file; the engine's split parts are lifted back losslessly
(see "Lifting from the σ-split", below). A file is written only for a **closed** closure (the engine's queue ran
empty, or 1 became a pivot); a censored closure (budget, memory, max-entries) writes no basis.

If 1 ∈ V_D then V_D = M_≤D (1 has degree 0 < D; multiplying by distinct chain variables reaches every square-free
monomial of degree ≤ D), and the file holds the identity basis (every row is a pivot with no other entry).

## Monomials

Chain variables u_{i,j} = X_i^{2^j}, i ∈ {1, …, m} (summand), j ∈ {0, …, k−1}; ring
F_p[X_1, …, X_m]/(X_i^{2^k} − 1). Modulo the chain every element is a combination of **square-free** monomials in the
m·k chain variables. A monomial is a `u64` mask: **bit (i−1)·k + j is u_{i,j}**. For m = 2 this is the paper's mask
(block 1 in bits 0..k−1, block 2 in bits k..2k−1; `csd/ic_bounds_d5.py`, `ic_dreg.pkm_generator`): the mask of
X_1^{e_1} X_2^{e_2} is `e_1 | e_2 << k`. Degree = popcount. m·k ≤ 56 (engine limit); free (non-chain) variables do not
occur for m = 2 and are not supported by v1 (`n_free` must be 0).

## Column orders

The file carries the full column list, so a reader never has to reconstruct an order; `column_order` names the rule:

- `0` — **degree descending, then mask descending** (as unsigned integers). This is the engine's unsplit order and
  equals the paper's `monomials(k, D)` order (`csd/ic_bounds_d5.py`: weight descending, then mask descending).
- `1` — **σ-orbit grouped**. σ swaps the two summands: σ(mask) = (mask >> k) | ((mask & (2^k − 1)) << k) for m = 2.
  Orbits {s, σs} are listed by their representative r = max(s, σs), in degree-descending-then-mask-descending order of
  r; each orbit contributes r, then σr if σr ≠ r. Degrees are non-increasing along this order too (σ preserves degree).
  This is the order written when the engine ran σ-split (the default for m = 2).

In both orders column 0 is the highest-degree monomial with the largest mask and the last column is the monomial 1
(mask 0).

## Rows: a pivoted row echelon form

`rank` rows, **sorted by pivot column ascending**. Row i has pivot column π_i with coefficient exactly 1 (implicit,
not stored) and further entries only at columns **> π_i**; pivots are distinct (π strictly increasing). The form is
**not** required to be fully reduced: a row may have nonzero entries at other rows' pivot columns. These properties
make the rows linearly independent and let a vector be reduced against the basis by one left-to-right pass:
for c in 0..n_cols, if the vector has x ≠ 0 at c and c is some π_i, subtract x · row_i.

Since degrees are non-increasing along the column order, a row's degree is the degree of its pivot column, and the
filtration profile is recomputed from pivots alone:
profile[d] = dim(V_D ∩ M_≤d) = #{i : deg(column π_i) ≤ d}, d = 0..D (the paper's `profile`; profile[D] = rank).
contains_one ⇔ the last column (monomial 1) is a pivot.

## Binary layout (little-endian throughout, no padding beyond what is listed)

| offset | type | field |
|:--|:--|:--|
| 0 | `[u8; 8]` | magic `HKYBAS01` (ASCII) |
| 8 | `u32` | format version = 1 |
| 12 | `u32` | m (summands; 2) |
| 16 | `u64` | p (prime, < 2^31) |
| 24 | `u32` | k (chain depth) |
| 28 | `u32` | D (closure degree) |
| 32 | `u32` | n_vars = m·k (chain variables) |
| 36 | `u32` | mask_encoding = 1 (u64 square-free mask, bit (i−1)·k + j = u_{i,j}) |
| 40 | `u32` | column_order (0 or 1, above) |
| 44 | `u32` | flags: bit 0 = contains_one; other bits 0 |
| 48 | `u64` | n_cols = Σ_{d ≤ D} C(m·k, d) |
| 56 | `u64` | rank |
| 64 | `u64 × (D+1)` | profile[0..=D] |
| … | `u64 × n_cols` | column monomials (masks), in column order |
| … | rows × rank | each row: `u32` pivot column, `u32` t (number of further entries), then t × (`u32` column, `u32` coefficient), columns strictly ascending and > pivot, coefficients in [1, p) |
| end − 32 | `[u8; 32]` | SHA-256 of every preceding byte of the file |

A reader must check: magic, version, mask_encoding = 1, n_cols matches the binomial sum and the column list (distinct
masks, popcount ≤ D, order rule), pivots strictly increasing, entries strictly ascending and > pivot, coefficients in
[1, p), the profile against the pivots, and the trailing SHA-256.

## Lifting from the σ-split (how the writer converts, losslessly)

With σ-split the engine stores V = V⁺ ⊕ V⁻ on orbit coordinates: a `+` row with coefficient x on representative r
stands for x·(e_r + e_{σr}) (x·e_r when σr = r); a `−` row stands for x·(e_r − e_{σr}) (fixed orbits do not occur in
`−`). Both parts are row echelon in representative order with monic leading coefficients. In column order 1 the lifted
rows are already echelon except where a `+` row h⁺ and a `−` row h⁻ lead at the same orbit r; the writer replaces that
pair by (h⁺ + h⁻)/2 (pivot r, coefficient 1, zero at σr) and (h⁺ − h⁻)/2 (zero at r, pivot σr, coefficient 1). This is
an invertible change of basis of the same space, so the span, the rank and the profile are unchanged. Without the
split the engine's rows are written as they are, in column order 0.

## Sidecar JSON (`<basis>.json`)

```json
{
  "schema": "hky.basis.sidecar/v1",
  "basis_file": "<file name>",
  "file_sha256": "<sha256 of the whole .bin file>",
  "content_sha256": "<the trailing SHA-256 stored in the file>",
  "bytes": 0,
  "p": 65537, "m": 2, "k": 14, "D": 6, "n_vars": 28,
  "n_cols": 499178, "rank": 448922, "contains_one": false,
  "profile": [0, 0, 2, 154, 4451, 73587, 448922],
  "column_order": 1, "mask_encoding": 1,
  "instance": {"path": "…", "sha256": "…", "prime": 65537, "a": 1, "b": 0, "target": [4271, 11720], "k": 14},
  "engine": {"name": "hky_fp_tower_fast2", "commit": "…", "binary_sha256": "…", "options": {}},
  "closure": {"muladds": 0, "rows_reduced": 0, "seconds": 0.0, "resumed": false},
  "created_utc": "…"
}
```

`file_sha256` is the checksum to compare across hosts; `content_sha256` is self-verifying inside the file.
