# Paper update, 19 November 2026

Folds the remaining finished runs of `paper-runs-20261118` into
`paper/curve-lfd-bounds.tex` on top of the manuscript finalized in PR #32
(`4f1084b`). The abstract is unchanged (180 words, 1291 characters); the
date of the manuscript is 19 November 2026.

## Records

`isogeny-labs/crypto`, branch `claude/lfd-ffd-paper-eprint-status-qq8fit`
(PR #4735), commit `12308b4877` (merge of `wt/paper-engine` `b4aba21adf`),
directory `research/hky_paper_runs_20261118/`: `results/J1/lift/`,
`results/J2/`, `results/J3/`, `results/J4s/`, `pod/LAUNCH.md`,
`PREDICTIONS.md`, `README.md` (section 6: the depth-16 unit instance);
skip criterion `research/hky_engine_skip_20261118/` (`PROOF.md`,
`README.md`); lift tool `tools/hky_lift/`. Labelled `[pcc:unverified]` there.
The J2 lift landed later, at commit `6fc6ccf847` (merge of `wt/paper-engine`
`752424ef08`), `results/J2/lift/` (`lift_ref.rust.json`,
`lift_ref.verify.json`, `lift_ref.stdout`, `lift_ref.inputs.sha256`).

| Run | Result | Label in the paper |
|---|---|---|
| J1 lift (unit p = 65537, k = 14) | Rust `hky_lift` (binary sha256 `6824a8bd…`) on the second engine's exported basis (`06fdd4e4…`): 1 ∈ V_7, dim V_7 = 1,683,218, dim F = 50,256, dim(V_7 ∩ M_≤6) = 499,178, profile equal field by field to `runs/lift-unit-p65537-k14-D6to7.json` (`lift_ref.verify.json`, `all_equal: true`; the regular-top count is informational, a different notion); input checks: 4 generators in V_6, 84 shifted random combinations stay in V_6; 4 rounds (551, 9166, 35,167, 5372 new falls); 5.60e14 multiply-adds, 13,374 s (about 3.7 h) on 32 vCPU, `maxrss_mb` 234,285 | corroboration; with the J1 closure, a second implementation of both halves of LFD(14) = 7 on the unit instance; not a certificate |
| J2 (unit p = 163841, k = 14) | 499,178 columns, profile (0,0,2,154,4451,73587,448922), 1 ∉ V_6, part ranks 224,688 + 224,234, closedness check 0 non-zero over 4,120,873 products, wall 49,797 s, peak 73,998,080 KiB, basis sha256 `3d1f2996…` | corroboration, not a certificate |
| J2 lift (unit p = 163841, k = 14) | same Rust `hky_lift` (binary sha256 `6824a8bd…`) on the J2 basis (`3d1f2996…`): 1 ∈ V_7, dim V_6 = 448,922, dim V_7 = 1,683,218, dim F = 50,256, dim(V_7 ∩ M_≤6) = 499,178, fixed point, profile equal field by field to `runs/lift-unit-p163841-k14-D6to7.json` (`lift_ref.verify.json`, `all_equal: true`; regular-top count 1,166,824 against 1,173,002 informational, a different notion); input checks: 4 generators in V_6, 84 shifted random combinations stay in V_6; 4 rounds (551, 9166, 35,167, 5372 new falls); 5.60e14 multiply-adds, 21,199 s (about 5.9 h) on 32 threads, `maxrss_mb` 234,299; records at `isogeny-labs/crypto` commit `6fc6ccf847` (PR #4735), `results/J2/lift/` | corroboration; the second implementation of both halves of LFD(14) = 7 now holds on both unit instances; not a certificate; hash comparison with W not run |
| J3 (unit p = 65537, k = 15) | 768,212 columns, profile (0,0,2,158,4761,82643,600547), 1 ∉ V_6, rank 600,547 (parts 300,535 + 300,012), closedness check 0 non-zero over 4,958,581 products, peak 79,614,980 KiB; resumed twice from checkpoints after host restarts | single run, not certified; adds no lower bound (Theorem `thm:restrict`) |
| J4s (unit p = 786433, k = 16) | instance `unit-p786433-a1-s10-k16` (y² = x³ + x, seed 10, R = (485051, 497175)); 1,149,017 columns, profile (0,0,2,162,5079,92323,770594), 1 ∉ V_6, rank 770,594 (parts 385,594 + 385,000); lead-cover skip criterion on, no closedness check; wall 10,238 s, peak 51,529,076 KiB | single run, weaker (skip, no check), not certified |

## Where in the paper

- `sec:depth14` (now titled "Depths eleven to sixteen"): intro sentence; the
  J1 paragraph ends with the run's qualifications; new paragraphs "A third
  lift, on the second engine's space" (J1 lift, both-halves statement, and
  the conditional input containment: if the second engine adds only
  admissible products then S ⊆ V_6 ⊆ W with dim S = dim W, so S = V_6 = W;
  if the third engine does, W = V_6; neither certified, hashes not
  compared), "The second engine on the unit instance over p = 163841" (J2),
  followed by the J2 lift (the same lift program on the J2 basis, equal
  field by field; both halves now hold on both unit instances), and
  "Degree six at depths fifteen and sixteen (one engine, not
  certified)" (J3, J4s, the depth-16 instance, the recursion check).
- `sec:recursion`: one sentence on the D = 6 recursion at k = 15, 16.
- `tab:profile`: row 14 evidence (both unit instances; third lift on both unit
  instances, one run each), row 15
  (dim V_6 = 600,547; "≤ 7 open"), new row 16; caption. Table set in
  `\footnotesize` with a wider evidence column so that it fits its page.
- `sec:conjecture`: prediction scoring (below).
- Introduction, Contributions item 6: scoring.
- Limitations: "The value 7 at k = 14" (third implementation on both unit
  instances, two chains of separately written programs on both), "Second
  engines", "Depth" (memory of the k = 15/16 runs, lift estimates), "The
  conjecture".
- Future directions, "Larger computations".
- Data appendix: `tab:data2` new rows (commit `12308b4`; J2 lift at commit
  `6fc6ccf`); prediction ids
  listed in the free-extension row.

## Prediction scoring

The test is of prediction `d6-base-15` in
`experiments/ic-bounds-freeext/predictions-2026-11-10.json` (frozen
2026-11-10): "W_15(V_6(14)) is NOT closed; W_16(V_6(15)) is closed and
equals V_6(16); the recursion holds for V_6 from k = 14", and of the
conjecture's "recursion" field for D = 6 (k ≥ 4D − 10 = 14). Only the
recursion clause is tested. The new text in `sec:conjecture`:

> The sixth (`d6-base-15`) states that the free extension of V_6(14) to
> depth 15 is not closed, that the free extension of V_6(15) to depth 16 is
> closed and equals V_6(16), and that the recursion holds for V_6 from
> k = 14. Its recursion clause held on the instances computed: the profiles
> of V_6(15) on the unit instance over p = 65537 and of V_6(16) on the unit
> instance over p = 786433 equal, entry by entry, the recursion applied to
> the certified profile of V_6(14) (Section 5.5). These are single runs of
> one engine, the second without a closedness check, and not certified. Its
> two closedness clauses are untested.

Contributions item 6 now reads: "three held in full, two agree with every
value computed, and of the sixth the recursion of the degree-6 profile held
at k = 15 and 16 in single runs of one engine, not certified, while its
closedness clauses are untested". The `lfd-table` values 7 at k = 15, 16
remain untested on the upper side (only the lower bound 7 is known).

## Checks made for this update

- Recursion arithmetic: (0,0,2,154,4451,73587,448922) → (0,0,2,158,4761,
  82643,600547) → (0,0,2,162,5079,92323,770594), equal to the J3 and J4s
  record profiles entry by entry. |M_≤6| = 768,212 (k = 15), 1,149,017
  (k = 16); |M_≤7| = 2,804,012, 4,514,873.
- Record fields read directly: `record.json` (`columns`, `profile`, `rank`,
  `contains_one`, `engine.part_rank`, `engine.verify_products`,
  `engine.verify_nonzero`, `engine.options.skip`, `engine.options.verify_closed`,
  `io.export.file_sha256`, `system.generator_matches_file`,
  `generator_matches_pkm_generator`), `job.json` (`wall_s_this_attempt`,
  `peak_rss_kb_this_attempt`, `binary_sha256`), `lift_ref.rust.json`,
  `lift_ref.verify.json` (J1 and J2; for J2 also `lift_ref.inputs.sha256`:
  binary `6824a8bd…`, basis `3d1f2996…`, expected record `4b5bb99e…`).
- Theorem `thm:restrict` carries the certified exclusion 1 ∉ V_6(14) to
  k = 15 on the same coefficients over p = 65537, so J3 adds no lower bound.
  It does **not** carry it to the k = 16 instance, which is over a different
  prime: there 1 ∉ V_6(16), and with the unit ideal LFD ≥ 7, rest on the
  single J4s run (the paper says so).
- J4s ran in 10,238 s against a prediction of 230,000 s (band 100,000 –
  950,000 s) written before the run; it is below the band and the records do
  not explain it. Noted in a LaTeX comment, not in the visible text.
- Not done: the canonical-hash comparison of either second-engine basis
  with W; any degree-7 step at k = 15, 16 (the
  lift is estimated at about 540 GB and 1.2 TB of memory).

## Build

Tectonic 0.17.0 (conda-forge), bundle `tlextras-2022.0r0` from
data1.fullyjustified.net, through `python3 scripts/build_ic_curve_paper.py`:
no warnings, no overfull boxes, no unresolved references; 59 pages (57
before); `curve-sweep.png` regenerated with different bytes (matplotlib
3.11.2). ePrint sheet: page count, digests and a revision paragraph
updated; the abstract file is unchanged.
`python3 scripts/check_publication_artifacts.py`: passed.

Rebuilt the same way after the J2 lift was folded in: no warnings, no
overfull boxes, no unresolved references; 59 pages; figures unchanged; PDF
SHA-256 `d0a446b7…`; ePrint sheet digests and revision paragraph updated;
`check_publication_artifacts.py`: passed.

## Still pending

No `% PENDING(paper-runs-20261118)` marker remains: the last one, for the
Rust lift on the J2 basis (`3d1f2996…`), is filled with the finished lift.
Not run: the canonical-hash comparison of either second-engine basis with W.

## Consolidation into this repository (later on 19 November)

Owner decision: the paper material of `isogeny-labs/crypto` PR #4735 moves here.
`second-engine/` (see its `README.md` and `PROVENANCE.md`) now holds the records
cited above (`records/hky_paper_runs_20261118/`, `records/hky_engine_opt_20261118/`,
`records/hky_engine_skip_20261118/`, the cited files of
`records/hky_pod_runs_20261117/`, `records/hky_engine_speed_20261117/README.md`),
the Rust lift (`lift/`, unchanged), and a standalone Rust extract of the second
engine (`engine/`: the engine files unchanged from crypto commit `7014d5dd2a`, the
paper-instance loader and the items of the pilot it uses). The runs were made with
binaries built from crypto commits `450d5b295b` (J1–J3, sha256 `39acdeac…`) and
`43d5c40d1f` (J4s, `36a047ff…`); the lift binary built here has the recorded
hash `6824a8bd…`.

Local gate of the extract (`second-engine/gate-20261119/`): the eleven recorded
closures that fit in 16 GB (degree 5 at k = 10..15 and on both unit instances at
k = 13, 14; degree 6 at k = 10, 11, 12) give equal columns, profile and
membership of 1, closed, with and without the skip criterion, and the same
multiply-add and row counts as the recorded gates of the two run binaries; the
engine's closedness re-check reduced 1,637,168 products to zero; the independent
basis checker accepts the exported bases at k = 10, D = 5 and 6 (closed, generator
in span); a SIGKILL and a SIGTERM resume reproduce the reference closure and its
basis file byte for byte, and that basis file is the one the run binary exported
in its own test; the three negative controls at k = 10, D = 6 (planted target
X = 11932, perturbed X, wrong a) equal the record; the engine's 20 unit tests
pass; the lift's validation suite gives 41 of 41 equal.

Paper: the second-engine paragraphs of `sec:depth14` and the `tab:data2` rows cite
`second-engine/…` instead of crypto commits; the data appendix gains one
provenance sentence, and the run-J1 paragraph one sentence on the extract's gate.
No value changes.
