# verify/: independent check of item 4 (see ../VERIFY.md for the verdicts)

Every run used `ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
The host was shared (load about 5 on 4 cores). Peak RSS was 2.0 GB, in the k = 8, D = 6 run.

## Commands

```sh
python3 scripts/ic_bounds_exact10_verify.py selftest --dmax 5 --per-case 2          # selftest.json, 12/12 agree, ~10 s
python3 scripts/ic_bounds_exact10_verify.py grid --ks 5 6 --dmax 6 --tag=-k56        # grid-k56.jsonl, 16 pairs, ~85 s
python3 scripts/ic_bounds_exact10_verify.py grid --ks 7 --dmax 6 --instances j1728-p7681-x936 rand-p257-x58 --tag=-k7
python3 scripts/ic_bounds_exact10_verify.py grid --ks 7 --dmax 6 --instances rand-p7681-x2244 rand-p12289-x11294 \
    j0-p7681-seed41 coset5-p12289-seed42 proper-p12289-seed43 proper-coset7-p7681-seed44 --tag=-k7   # ~135 s per instance
python3 scripts/ic_bounds_exact10_verify_compare.py                                    # compare.json: all producer rows vs mine
python3 scripts/ic_bounds_exact10_verify_block.py                                      # block-rerun.json (existing block engine, 5 rows)
python3 scripts/ic_bounds_exact10_verify.py midscale --system 0 --dmax 6 --mix 320 --mix-from 5 --tag=-k8s0   # 1.2 ks
python3 scripts/ic_bounds_exact10_verify.py midscale --system 14 --dmax 5 --mix 320 --mix-from 5 --tag=-k8
python3 scripts/ic_bounds_exact10_verify.py midscale --system 6 --dmax 5 --mix 320 --mix-from 5 --tag=-k8
python3 scripts/ic_bounds_exact10_verify.py midscale --depth-dir experiments/ic-chain-depth9 --system 0 --dmax 5 --mix 360 --mix-from 5 --tag=-k9   # 760 s
python3 scripts/ic_bounds_exact10_verify.py dryrun --cols                              # dryrun-k10.json
```

`k10-onto.json` came from an inline call to `eval_ranks` on `Space(10, 5)`, for the two proper
depth-10 instances (p = 13313, x = 7412; p = 15361, a = 4021, b = 9917, x = 7776).

Seeds:
- The mixing seeds are `numpy.default_rng(1000 + D)` (grid) and `2000 + D` (midscale).
- The selftest targets come from `random.Random(31 p + k)`.
- The grid targets are the producer's recorded `target_x`.

Mixing affects only speed. Every reported dimension is either `closed` (deterministic worklist
exhausted) or `reached` (rank = |M_<=D| - rank ev(M_<=D)).

## Results

| check | result |
|---|---|
| grid k = 5, 6, 7, D = 2..6, 8 instances | all 72 producer rows (and 48 v2 rows) match: dims, low dims, 1 in V_D, targets, statuses |
| profiles | 24/24: LFD 4 (k = 5), LFD 5 (k = 6, 7), falls [4] / [4,5] |
| k = 8 systems 0 / 14 / 6 | V_4 = 1061, closed. V_5 = 6633 / 6885 (unit) / 6803, reached |
| k = 8 system 0, D = 6 | V_6 = 14641 = 14893 - 252, reached |
| k = 9 system 0 | V_4 = 1301, closed. V_5 = 12592, reached (low 4024) |
| k = 10 instance construction | x = 9362 / 7412 / 7776 with 0 / 98 / 60 points |

## Hashes (sha256)

Scripts:
```
a6f29e0c5ce8695fe20f90fb33327fbaea9c8c702e521d1cf12becb69fc22dd9  scripts/ic_bounds_exact10_verify.py
b788abedd9c02ef51374fd2db2e8d7b702d50516d138f0827dd8491d38adab04  scripts/ic_bounds_exact10_verify_block.py
b31367c72291181430d180ee8eb8e05b8cb46e0a4b640ac000741c543aa9113b  scripts/ic_bounds_exact10_verify_compare.py
```

Outputs:
```
97081670137e0c1b26770a452c398a6c021e285add3eaa5ef64955e5eb171360  grid-k56.jsonl
e9089de6aac2b9f48aa53ae75e27aac47397526c6179dac45059eca4e252da5d  grid-k7.jsonl
45231d2ac2ff3872eba449a3681c31c5b71d2dafe8d2a1f27536362f8ce488d0  midscale-k8s0.jsonl
e52faf1501d99f67446fcb4fc0f56e7a8c2905bb320ccf1e1653cd3ac05ed25b  midscale-k8.jsonl
2e2232a3c7ac685853741f3da4970890898add0a7c63f2a65a954cfa0d71f830  midscale-k9.jsonl
bd47449dc17917032458257be96f862707b0655ae975cf14bccad71287ca0a60  compare.json
66dd9efc4ac9af23c53159bc3e691f2dcbef0f635d95a8d29a1b49e2deb78578  block-rerun.json
7afa4a2f604e083d87dcba1dbb6cdff8794ea1c8e5d0041e6022682fd1d248c2  selftest.json
28d15b8c359b5ca4a261681ca6a191b36aef99f74ed97d45c379a1607777bff2  dryrun-k10.json
c1cfb44f2fe8b4efda8ced0d181550efb1502d365bc03a248215cd32bbf1ad62  k10-onto.json
```
