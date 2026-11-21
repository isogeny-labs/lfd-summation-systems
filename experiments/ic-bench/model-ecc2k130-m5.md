n=131 m=5  t_lin=1849.72 ns  t_cand=260 ns  t_add=150.05 ns  t_la=10 ns (assumed)  rho: 2^60.81 iterations

| k | log2 mu | log2 relations (s) | log2 algebra (s) | log2 total (s) | total core-years | algebra memory | gap to rho (log2) |
|---|---|---|---|---|---|---|---|
| 42 | 72.1 | 75.9 | 58.7 | 75.9 | 2.29e+15 | 172 TiB | 37.8 |
| 43 | 77.1 | 74.9 | 60.7 | 74.9 | 1.14e+15 | 344 TiB | 36.8 |
| 44 | 82.1 | 73.9 | 62.7 | 73.9 | 5.72e+14 | 688 TiB | 35.8 |
| 45 **←** | 87.1 | 73.9 | 64.7 | 73.9 | 5.68e+14 | 1.38e+03 TiB | 35.8 |
| 46 | 92.1 | 75.2 | 66.7 | 75.2 | 1.41e+15 | 2.75e+03 TiB | 37.1 |
| 47 | 97.1 | 77.1 | 68.7 | 77.1 | 5.22e+15 | 5.5e+03 TiB | 39.0 |
| 48 | 102.1 | 79.1 | 70.7 | 79.1 | 2.07e+16 | 1.1e+04 TiB | 41.0 |

rho: 3.03e+11 core-seconds = 9.6e+03 core-years
floor with a FREE relation finder (t_lin = 0): 2^37.5 core-seconds at k = 30 (-0.7 vs rho in log2; algebra memory 2^-4.8 TiB)
required speed-up of the relation finder to reach rho (k re-optimized): 2^49.5 = 7.87e+14x, i.e. t_lin <= 2.35e-12 ns, break-even at k = 31
per-target budget for ANY relation finder (all decompositions of one R, one core), best k = 31: 2.33e+07 s  (mu = 2^17.1, 2^30.1 relations, 167 s per relation)
Assumptions: |FB| = 2^k, mu = |FB|^m/(m! 2^n), t_la is not measured, kernel limits (e.g. the Rust kernel needs k >= 35 at n = 131) are not modelled.
