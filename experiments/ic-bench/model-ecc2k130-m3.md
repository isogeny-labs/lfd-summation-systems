n=131 m=3  t_lin=1849.72 ns  t_cand=260 ns  t_add=150.05 ns  t_la=10 ns (assumed)  rho: 2^60.81 iterations

| k | log2 mu | log2 relations (s) | log2 algebra (s) | log2 total (s) | total core-years | algebra memory | gap to rho (log2) |
|---|---|---|---|---|---|---|---|
| 41 | -10.6 | 72.6 | 56.0 | 72.6 | 2.29e+14 | 72 TiB | 34.5 |
| 42 | -7.6 | 71.6 | 58.0 | 71.6 | 1.14e+14 | 144 TiB | 33.5 |
| 43 | -4.6 | 70.6 | 60.0 | 70.6 | 5.72e+13 | 288 TiB | 32.5 |
| 44 **←** | -1.6 | 69.6 | 62.0 | 69.6 | 2.87e+13 | 576 TiB | 31.5 |
| 45 | 1.4 | 69.6 | 64.0 | 69.6 | 2.89e+13 | 1.15e+03 TiB | 31.5 |
| 46 | 4.4 | 70.9 | 66.0 | 71.0 | 7.27e+13 | 2.3e+03 TiB | 32.8 |
| 47 | 7.4 | 72.8 | 68.0 | 72.8 | 2.7e+14 | 4.61e+03 TiB | 34.7 |

rho: 3.03e+11 core-seconds = 9.6e+03 core-years
floor with a FREE relation finder (t_lin = 0): 2^43.0 core-seconds at k = 34 (+4.9 vs rho in log2; algebra memory 2^-0.9 TiB)
required speed-up of the relation finder to reach rho: none suffices at this m (the floor is above rho)
per-target budget for ANY relation finder (all decompositions of one R, one core), best k = 31: 0.126 ns  (mu = 2^-40.6, 2^30.1 relations, 207 s per relation)  -- shorter than forming the target (150 ns): infeasible
Assumptions: |FB| = 2^k, mu = |FB|^m/(m! 2^n), t_la is not measured, kernel limits (e.g. the Rust kernel needs k >= 35 at n = 131) are not modelled.
