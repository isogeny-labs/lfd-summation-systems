# Paper update, 18 November 2026

The first pass (commit `431bd79`) edited `paper/curve-lfd-bounds.tex` only.
The second pass (below, "Second pass") folds in the second engine's run on
the unit instance over p = 65537, shortens the abstract, adds a significance
paragraph to the introduction, and rebuilds the PDF, the figures, the ePrint
abstract and the digests with Tectonic.

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

## Second pass, 18 November 2026

### Result folded in (run J1 of `paper-runs-20261118`)

Records in `isogeny-labs/crypto`, branch `wt/paper-engine`, commit `98b7d87056`
(local in the session that produced them; it must be pushed before the
citation in the data appendix resolves), directory
`research/hky_paper_runs_20261118/`: `README.md` (system correspondence,
gate), `pod/LAUNCH.md`, `results/J1/record.json`, `job.json`,
`basis.bin.json`; engine notes `research/hky_engine_opt_20261118/README.md`.
Labelled `[pcc:unverified]` there.

| Item | Value |
|---|---|
| Engine | `src/cryptanalysis/hky_fp_tower_fast2.rs` (hky3), driver `examples/pkm_tower_pilot.rs` (`--instance` mode); engine commit `450d5b295b`; binary sha256 `39acdeac81689a49fa98a9b71e9924bdc29418d80247bece72e372aa41cac74f` |
| Instance | `unit-p65537-a1-s10-k14` (y² = x³ + x over p = 65537, R = (4271, 11720)), file sha256 `fa06148e0276dd459a45ad45e08b751a699ceb2a9a440a8ee07770a9ebeaa1ab`; generator checked against the paper's formula and the file (`generator_matches_pkm_generator`, `generator_matches_file`) |
| Gate before the run | same binary reproduced 11 recorded closures: D = 5 at k = 10..15 (both unit instances at k = 13, 14), D = 6 at k = 10, 11, 12 |
| Closure, D = 6, k = 14 | 499,178 columns; profile (0,0,2,154,4451,73587,448922); 1 ∉ V_6; rank 448,922; parts 249,824 / 249,354 columns, ranks 224,688 + 224,234 |
| Closedness check (engine-internal, `--hky3-verify-closed on`) | closed; 0 non-zero over 4,120,873 reductions; 18,495 s |
| Resources | wall 26,723 s on 32 shared vCPU (EPYC 7713); peak RSS 73,882,436 KiB (`peak_rss_kb_this_attempt`, about 70 GiB; with the check) |
| Basis | `hky.basis/v1`, file sha256 `06fdd4e4d00095ad02d9a9cea735ca59eed62843552ee14b43726f76d1ad9d1c` |

Where in the paper: `sec:depth14`, new paragraph "The second engine on the
unit instance over p = 65537"; `tab:profile` row 14 and caption; Limitations
items "The value 7 at k=14", "Second engines", "Depth", "The conjecture";
Future directions, "Larger computations"; data appendix, new row in the
second data table.

Label: a corroboration, not a certificate. The check has the form of the
paper's closedness sweep but runs inside the engine; the paper's checker has
not read the exported basis; one run. Not done: the canonical-hash
comparison of this basis with the paper's stored W, and the second lift on
this basis (running). So W = V_6(14) and the lift are not settled. The paper
states one consequence that needs no hash: since the checker gives
W ⊇ V_6(14) with dim W = 448,922, W ⊆ V_6(14) also follows if this engine adds
only admissible products (neither engine is certified to).

### Abstract and significance

- Abstract shortened from 432 to 180 words (1291 characters). It keeps the
  kernel-checked upper bounds, the at-every-depth lower bounds, the certified
  values 4, 5, 6 through k = 13, the depth-14 split (lower bound certified,
  upper bound by the lift), the conjecture, the written growth and
  solving-degree results, and the no-cost statement; it drops the sweep size,
  the m-summand and binary details and the second-lift sentence (all still in
  the body).
- Introduction: new paragraph "Why this matters" (relation search versus
  Pollard rho; what the proved bounds cap and what they rule out; no cost
  claim).
- Date of the manuscript set to 18 November 2026.

### Build and layout

- Tectonic 0.17.0 (conda-forge), default bundle v33 (`tlextras-2022.0r0`,
  fetched from data1.fullyjustified.net), through
  `python3 scripts/build_ic_curve_paper.py`: no warnings, no overfull boxes,
  no unresolved references. The regenerated figures changed bytes
  (matplotlib 3.11.2 here).
- The data table overflowed its page ("Float too large" by 73 pt) once the
  new row was added; it is now split into `tab:data` (formal development,
  proofs, reviews) and `tab:data2` (computations).
- 52 pages (the source of `431bd79` alone builds to 51 pages with this
  toolchain; the committed PDF before this pass had 49 and predated
  `431bd79`).

### Still pending

The `% PENDING(paper-runs-20261118)` markers remain (13, rewritten where J1
filled part of them): J2 (p = 163841), the canonical-hash comparison and the
second lift on the second engine's bases, k = 15 on unit-p65537 and k = 16 on
the new unit instance over p = 786433.
