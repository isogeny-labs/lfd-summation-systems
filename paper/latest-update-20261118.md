# Paper update, 18 November 2026

Source edits to `paper/curve-lfd-bounds.tex` only. The PDF, the figures, the
ePrint abstract and the digests in `eprint-curve-lfd-bounds/` were not
rebuilt or changed (no Tectonic in the editing environment); the owner
rebuilds. Until then `scripts/check_publication_artifacts.py` reports the
manuscript digest as stale. The abstract itself is unchanged.

## Results folded in

All from the repository `isogeny-labs/crypto` at commit `89a2a70`
(branch `claude/lfd-ffd-paper-eprint-status-qq8fit`), labelled
`[pcc:unverified]` there.

| Item | Record | Where in the paper |
|---|---|---|
| Second engine (Rust `hky2`, `src/cryptanalysis/hky_fp_tower_fast.rs`, generator `examples/pkm_tower_pilot.rs`), k = 14, random curve and coset base over p = 2013265921: degree-6 closure without 1, dim 448,922 of 499,178, profile (0,0,2,154,4451,73587,448922), part ranks 224,688 + 224,234, 28,747 s on 32 vCPU | `research/hky_pod_runs_20261117/results/t14_hky2.jsonl` (sha256 `093245bb…`), `t14_progress.jsonl`, `t14_rss.log` | `sec:depth14`, new paragraph "A second engine on a fourth prime"; `tab:profile` rows 13 and 14; Limitations |
| Same engine, k = 13, two targets: V_5 dim 39,657 without 1, V_6 = M_{≤6} (313,912) contains 1 | `results/t13_hky2.jsonl` (sha256 `aeba21e3…`) | same |
| Split presentation, later single-engine readings: (6,2) 7, (7,2) 8, (8,2) 8, (9,2) no 1 through 8, (6,3), (7,3), (8,3) no 1 through 7; (5,3) 8 and the calibration cells (4,2), (4,3), (5,2) reported in the summary only | `results/prog_m6-T2.jsonl`, `b_m7-T2.jsonl`, `b_m8-T2.jsonl`, `m8pod_m8-T2*.jsonl`, `prog_m9-T2-cols12M.jsonl`, `prog_m6-T3-d3-long.jsonl`; `RESULTS.md` | `sec:msummands`, new paragraph "Later readings of one engine (not certified)"; kept out of `tab:splitm` |
| Depth reach: size of M_{≤6} at k = 15, 16, 17, 18, 20 and of M_{≤7}, M_{≤8} at k = 18 | computed as Σ_{w≤c} C(2k, w) | Limitations, item "Depth" |

## Checks made for this update

- The full degree profile of the k = 14 record equals the paper's: c = 3: 29;
  c = 4: (0,0,2,143,2621); c = 5: (0,0,2,154,4323,47849); c = 6:
  (0,0,2,154,4451,73587,448922), with split columns 249,824 / 249,354 and part
  ranks 224,688 / 224,234. At k = 13 the c = 5 profile (0,0,2,150,4021,39657),
  part ranks 19,864 / 19,793 and V_6 = M_{≤6} equal the paper's.
- The instance is the paper's coset presentation: in
  `examples/pkm_tower_pilot.rs` the Kummer tower has steps y_{j+1} = y_j^2,
  top constant 1 and abscissa x = g·y_0 (also so at the engine README's base
  commit `94479f7e3e`). The weighted homogeneity
  S_3(cx_1, cx_2, t; a, b) = c^4 S_3(x_1, x_2, t/c; a/c^2, b/c^3) was checked
  numerically mod 2013265921.
- The engine has no closedness sweep; the paper says so where the numbers
  appear and does not use them for any bound.
- No per-closure record for m = 4, 5 exists in `results/`; the (5,3) value is
  cited as reported in `RESULTS.md` only.
- A draft-mode pdflatex compile of a scratch copy (figures replaced by
  boxes) gave no errors, no undefined references and no overfull boxes.
  This is not the Tectonic build.

## Pending

Every place where the runs `paper-runs-20261118` (second engine on the unit
instances over p = 65537 and 163841 at k = 14; depth 15 on the unit
instance over p = 65537; depth 16 on a new unit instance over a prime with
2^17 | p − 1) will go carries a LaTeX comment `% PENDING(paper-runs-20261118)`.
No pending result is in the visible text.
