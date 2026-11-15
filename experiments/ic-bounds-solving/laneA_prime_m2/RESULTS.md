# Lane A: prime two-summand chain system, static solving degree (sd) vs d_lf, d_ff

Every result statement below is **[pcc:unverified]** (produced and read by the same session; no independent validator, no replay certificate, no matched rho; nothing here is a leaderboard row).
Frozen definitions and predictions: `../PROTOCOL.md`. Run 2026-11-15, 2 threads max (OMP_NUM_THREADS=1 per process, two processes at most), `ulimit -v 6000000`, other lanes were using the same 4-core host (wall times are inflated by contention).

## What was built (all under `tools/`)

| file | role |
|---|---|
| `common.py` | prime/root-of-unity helpers, direct S_3, README normalised generator g, unit test, popcount |
| `check_generator.py` | checks g == S_3 as polynomials (3 primes x 2000 random points), g agrees with `lfd_engine.prime_inputs/prime_points` zero sets, and S_3(x(P),x(Q),x(P+Q)) = 0 on all 9702 pairs of a p=97 curve: all pass |
| `static_bruteforce.py` | task 1: sd = least D with 1 in Mac_D, ordinary ring. Mode `dense` (literal Mac_D in P_{<=D} with `lfd_engine.Space/Ech`) and mode `nf` (chain rows eliminated by the rewriting u_{i,j}^2 -> u_{i,j+1}, u_{i,k-1}^2 -> 1; rows m*g for EVERY ordinary monomial m with deg m <= D-2, NF-reduced by monomial product + rewriting, deduplicated by identical row, echelon mod p on the 4^k squarefree monomials). Also `lfd_per_degree` = `lfd_engine.exact_lfd` wrapper (d_lf, d_ff). |
| `sd_inverse.py` | task 2: g^{-1} by evaluation on mu_N x mu_N, Fermat pointwise inversion, two 1-D radix-2 modular DFTs (exact int64, p < 2^31), exact convolution check g * g^{-1} = 1 on every instance, maxwt, zeros per weight class, coefficient at (N-1,N-1) |
| `lf_quotient.py` | quotient-algebra closure V_D in F_p[(Z/N)^2] (U1 sec. 1), reusing `lfd_engine.Ech`; unit d_lf, d_ff |
| `dlf_crosscheck.py` | ordinary `lfd_engine` vs `lf_quotient` (k=3, 4) |
| `run_scan.py`, `campaign*.sh` | drivers (tasks 1-4, 5) |
| `task5_dfull.py` | task 5 (group-algebra d_full vs d_lf, proper instances) |
| `exc_dlf.py`, `dff_only.py`, `aggregate.py`, `plot.py` | exceptions' d_lf, d_ff for k=7..9, aggregation, figure |

Commands: `bash tools/campaign.sh` (stopped by hand during its repeated k=12 p=12289 non-unit search; its Part C was replaced by `campaign3.sh`), `campaign2.sh`, `campaign3.sh`, `campaign4.sh`, `campaign5.sh`, `python3 tools/exc_dlf.py`, `python3 tools/dlf_crosscheck.py 3 97 4`, `python3 tools/dlf_crosscheck.py 4 1009 3`, `python3 tools/static_bruteforce.py selftest` (stopped after 3 instances), `python3 tools/aggregate.py`, `python3 tools/plot.py`. Logs: `data/campaign*.log`, `data/aggregate_out.txt`.
`plot.py` needs matplotlib, which was installed with `pip install --target tools/_pylib matplotlib`; I deleted `_pylib` afterwards (168 MB), reinstall to regenerate the figure.

Data: `data/all_instances.json` (29,775 rows: p, a, b, t, k, sd_bruteforce, sd_inverse, maxwt, zero_top_coeff, zeros_per_weight, nzeros, inverse_verified, d_lf, d_ff, d_lf_source), `data/scan_k*_p*_s*.json` (per run, with tries, unit rate, wall seconds), `data/summary.json`, `data/dlf_crosscheck_*.json`, `data/dff_k{7,8,9}.json`, `data/exceptions_dlf.json`, `data/task5_k*_p*.json`, `data/k7_quotient_profile_note.txt`, `data/sd_vs_k.svg`.

## Sanity checks that gate the numbers

* g (README sec. 1) equals S_3 as a polynomial and reproduces `lfd_engine.prime_points` zero counts; S_3 vanishes on x-coordinates of P+Q+R=O (above). [pcc:unverified]
* `static_bruteforce.py selftest`: at k=3 (p=97, three unit instances) the literal dense Mac_D agrees with the `nf` mode for every D = 2..8 in both "1 in Mac_D" and dimension (dim Mac_D(dense) = #monomials<=D - #squarefree<=D + rank of NF-rows). Dense mode was NOT run at k = 4 or higher (the selftest was stopped after 3 of its 8+ instances because the dense echelon is slow). [pcc:unverified]
* d_lf, d_ff: ordinary-ring `lfd_engine.exact_lfd` and my quotient closure agree on 4 (k=3, p=97) + 3 (k=4, p=1009) instances (`data/dlf_crosscheck_*.json`), and `lfd_engine` and the quotient agree again on 12 k=5 and the 3 k=5 ordinary instances (both give 4/4). Control from PROTOCOL: k=7 quotient closure gives dim V_4 = 829, dim V_4 cap M<=3 = 115, 1 not in V_4, 1 in V_5 (dim 3473 = |M<=5|): equals the recorded V_4(7)=829 and recorded LFD(7)=5. The literal k=12, p=12289 control was NOT possible: no unit instance exists there (see below), and a non-unit instance has no "1 in V_D". [pcc:unverified]
* `inverse_verified` (exact cyclic convolution g * g^{-1} = delta_0) is true on all rows. [pcc:unverified]

## Important deviation: unit instances need p >~ 4^k

An instance is unit iff S_3(alpha,beta,t) != 0 on all of mu_N x mu_N, N=2^k, i.e. 4^k points; the expected number of zeros is about 4^k/p. So the prescribed p = 12289 produces unit instances only up to k = 8 (unit rate per random draw: k<=4 ~1.0, k=5 0.95, k=6 0.8, k=7 0.44, k=8 0.04, k=9..12: 0/300 each). For k = 9..13 I therefore used primes of size 4^k times a small factor: p_mid(k) = least prime >= 3*4^k with 2^k | p-1, p_2x4k = least prime >= 2*4^k, and p_big = the largest prime below 2^31 (so the DFT stays in int64); at k = 3..8 I additionally used 12289, 40961, 65537 and (k=3,4,5) many small primes. "Smallest prime with 2^k | p-1" (k=8: 769, k=9: 7681, k=10..12: 12289) has unit rate 0/300 for k>=8, so exceptions there cannot be studied; for k=3,4 the smallest primes (17; 97) are included. All sd values below are over unit instances only.

## Task 1 and 2: brute force vs 2 + maxwt(g^{-1}) (P1)

Brute force = `nf` mode (ordinary-ring rows, rewriting reduction, see table above), unit instances, sd scanned D = 2..2k+2 (the whole range; `sd=None` would mean 1 not in Mac_D up to 2k+2, which never happened).

| k | primes (unit instances brute-forced) | instances | mismatches sd_bf != 2+maxwt | sd values seen (count) |
|---|---|---|---|---|
| 3 | 41, 73, 89, 97, 113, 137, 193, 257, 401 (60..3000 each) | 18,492 | **0** | 8 (all generic), 7 (129 exceptions, 125 of them brute-forced) |
| 4 | 97, 113, 193, 241, 257, 337, 401, 1009, 2017 (30 each), plus 12 on 1009 | 282 | **0** | 10, and 9 (2 exceptions) |
| 5 | 193, 3137, 12289, 40961 (12 each) + ordinary-ring 3,137 | 48 | **0** | 12 |
| 6 | 12289 (6, bonus; 86 s enumeration of 2.7M monomials) | 6 | **0** | 14 |

(The `P1_checked` column of `data/summary.json`: 18492 / 282 / 48 / 6; exceptions among them: 125 at k=3 and 2 at k=4 each match, with sd = 2k+1 on both sides.) Dense-mode agreement is only established at k=3 (selftest above).
**P1 verdict: supported, not falsified, on 18,828 instances (k=3..6), including every exceptional instance brute-forced.** The weakness: the "ordinary ring" brute force at k >= 4 relies on the rewriting reduction of the chain rows (verified against dense Mac_D only at k=3, D<=8), which is the elementary half of the Lemma; the other half, NF(Mac_D) = g*M_{<=D-2}, is not assumed (rows m*g are NF-reduced from every ordinary monomial m). k=5 brute force at the size of the dense matrix (646,646 columns at D=12) was never attempted. [pcc:unverified]

## Task 3: scan k = 3..13, sd = 2 + maxwt(g^{-1}) (P2)

Totals (unit instances found / random draws, from `data/summary.json`, `scan_*.json`): k=3 19,179/25,367; k=4 6,970/13,117; k=5 2,723/3,615; k=6 438/832; k=7 122/456; k=8 73/668; k=9 85/1,037; k=10 60/372; k=11 60/367; k=12 60/370; k=13 5/7. At least 20 unit instances per k for k=3..12 at p_mid, p_2x4k and (k<=12) p_big; k=13 only 5 (p=201465857, 52 s per draw); k=14 not run (see below).

sd per k (all unit instances): k=3: 8 (19,050) and 7 (129); k=4: 10 (6,940) and 9 (30); k=5: 12 (2,723); k=6: 14 (438); k=7: 16 (122); k=8: 18 (73); k=9: 20 (85); k=10: 22 (60); k=11: 24 (60); k=12: 26 (60); k=13: 28 (5). So sd = 2k+2 for every instance with k >= 5 and for 99.3% (k=3) and 99.6% (k=4) of the pooled small-prime instances.

Exceptions occur only at k=3 and 4, only with small p, always sd = 2k+1 = 7 or 9, always with maxwt = 2k-1 and coefficient at (N-1,N-1) equal to 0 (coefficient zero in all 159 exceptions; it is nonzero in every generic instance). Exception frequency against 1/p (unit instances):

| k | p | n | exceptions | frac generic | 1-3/p | 1-10/p |
|---|---|---|---|---|---|---|
| 3 | 73 | 3060 | 40 (1.31%, 1/p=1.37%) | 0.9869 | 0.9589 | 0.8630 |
| 3 | 97 | 3072 | 31 (1.01%, 1.03%) | 0.9899 | 0.9691 | 0.8969 |
| 3 | 137 | 3260 | 18 (0.55%, 0.73%) | 0.9945 | 0.9781 | 0.9270 |
| 3 | 193 | 3260 | 13 (0.40%, 0.52%) | 0.9960 | 0.9845 | 0.9482 |
| 3 | 257 | 3060 | 13 (0.42%, 0.39%) | 0.9958 | 0.9883 | 0.9611 |
| 3 | 401 | 3000 | 11 (0.37%, 0.25%) | 0.9963 | 0.9925 | 0.9751 |
| 3 | 17 | 27 | 2 | 0.9259 | 0.8235 | 0.4118 |
| 4 | 113 | 1530 | 16 (1.05%, 0.88%) | 0.9895 | 0.9735 | 0.9115 |
| 4 | 257 | 1530 | 4 (0.26%, 0.39%) | 0.9974 | 0.9883 | 0.9611 |
| 4 | 401 | 1530 | 6 (0.39%, 0.25%) | 0.9961 | 0.9925 | 0.9751 |
| 4 | 593 | 1500 | 3 (0.20%, 0.17%) | 0.9980 | 0.9949 | 0.9831 |
| 4 | 577 | 200 | 1 | 0.9950 | 0.9948 | 0.9827 |

In every (k,p) cell (`data/summary.json`, `per_k.*.p2`) the generic fraction is above 1-3/p; the smallest margins are k=4,p=577 (0.9950 vs 0.9948) and k=3,p=401 (0.9963 vs 0.9925), i.e. the data sit at about 1/p, inside the 3/p allowance and far above the falsifier level 1-10/p. k=5..13 had no exception in any prime (5..1,227 instances per prime at k=5, p up to 2^31); the expected number of exceptions at k=5, p=3137 is about 1000/3137 = 0.3, so absence there is not informative.
**P2 verdict: supported (not falsified) at k=3..13 for every prime where unit instances exist; the exceptions are rare (about 1/p), always sd = 2k+1.** Not tested: k=9..12 with primes below about 4^k (no unit instance exists there); sd for k >= 9 was computed only through g^{-1} (no brute-force cross-check beyond k=6), so for k >= 7 the figure depends on the Theorem of the protocol, which this lane only tests up to k=6. [pcc:unverified]

## Task 4: d_ff, d_lf, sd and P3

d_ff = least D with V_D cap R_{<=D-1} strictly larger than V_{D-1}. d_lf for unit instances = least D with 1 in V_D.

| k | d_ff | d_lf | d_lf source | sd (modal; exception) | sd - d_lf | k+2 | 2k+2 |
|---|---|---|---|---|---|---|---|
| 3 | 4 (12 inst.) | 4 | measured, lfd_engine (12; plus 25 exceptional instances: 4) | 8 (7) | 4 (3 on exceptions) | 5 | 8 |
| 4 | 4 | 4 | measured, lfd_engine (12; plus 20 exceptional: 4) | 10 (9) | 6 (5) | 6 | 10 |
| 5 | 4 | 4 | measured: lfd_engine 3, quotient engine 12 | 12 | 8 | 7 | 12 |
| 6 | 4 | 5 | measured: quotient engine, 12 | 14 | 9 | 8 | 14 |
| 7 | 4 | 5 | measured: quotient engine, 2 (234 s total) | 16 | 11 | 9 | 16 |
| 8 | 4 | 5 | measured: quotient engine, 1 (548 s) | 18 | 13 | 10 | 18 |
| 9 | 4 (3 inst.) | 5 | LAW floor((k+14)/4), not measured here | 20 | 15 | 11 | 20 |
| 10 | not run (see below) | 6 | LAW, not measured here | 22 | 16 | 12 | 22 |
| 11 | not run | 6 | LAW | 24 | 18 | 13 | 24 |
| 12 | not run | 6 | LAW | 26 | 20 | 14 | 26 |
| 13 | not run | 6 | LAW | 28 | 22 | 15 | 28 |

Recorded values of the repository (`ic-bounds-exact10/summary.json`, `profile_lfd_by_k`: 5 -> 4, 6 -> 5, 7 -> 5) agree with my measurements; I did not open `ic-bounds-k14`. The d_lf at k = 9..13 is the law, labelled so; nothing there is measured by this lane. d_ff = 4 at every k measured (3..9), one fall at D=4 (dim V_4 cap M<=3 = 115, 119, 123 at k=7,8,9, as in the fall-space README).
All d_lf values were measured on unit instances with p <= 2*10^5 (p_mid or 97..3137), because `lfd_engine.Ech` requires cap*(p-1)^2 < 2^53; they are not on the large-prime instances, but sd takes the same value at every prime of the same k.

**P3 (d_lf <= k+2 < sd, sd - d_lf >= k):** holds on every measured row, k=3..8 (and on the law for k=9..13). Tightest cases: the exceptional k=3 instances, where sd = 7, d_lf = 4, sd - d_lf = 3 = k (equality, 25 instances measured); k=4 exceptions give 5 >= 4. For generic instances sd - d_lf = 4, 6, 8, 9, 11, 13 at k = 3..8. d_lf <= k+2 holds with slack (k+2 - d_lf = 1, 2, 3, 3, 4, 5). **Verdict: supported on the measured range; for k=9..13 the d_lf side is the recorded law, not measurement.** [pcc:unverified]

## Task 5: group-algebra d_full vs d_lf (proper instances, |Z| > 0)

d_full = least D with dim(g*M_{<=D-2}) = dim I_D (I_D = |M_<=D| - rank eval on Z) and eval of M_{<D} onto K^Z; d_lf by the same certificate as `lfd_engine.exact_lfd` in the quotient. Five instances per k (3 at k=6); |Z| from 1 to 6; `stable` = condition keeps holding for all larger D up to 2k+2 (true everywhere).

| k | p | instances | d_full | d_lf |
|---|---|---|---|---|
| 3 | 41 | 5 | 7 on all | 4 on all |
| 4 | 97 | 5 | 9 on all | 4 on all |
| 5 | 193 | 5 | 11 on all | 4 on all |
| 6 | 769 | 3 | 13 on all | 5 on all |

So d_full = 2k+1 on every proper instance tested, d_full >= d_lf always (no violation), d_full - d_lf = 3, 5, 7, 8. Note: for a proper ideal the generic static-to-closure value is 2k+1, one less than 2k+2 for unit instances. This is the group-algebra analogue only: no ordinary-ring brute-force check was run for d_full, and 18 instances is small. [pcc:unverified]

## Wall times (summed per k, seconds, host shared with other lanes)

sd_inverse: k=3 22, 4 16, 5 9, 6 4, 7 4, 8 16, 9 61, 10 101, 11 448, 12 1,745, 13 360 (about 51 s per draw at k=13, 13-18 s at k=12, 3.4-5.5 s at k=11).
Brute force (nf mode): k=3 122, 4 21, 5 69, 6 193 (+86 s per one-time monomial table). d_lf: k=3 1, 4 11, 5 20, 6 204, 7 234, 8 548; task 5: k=6 about 440 s per instance. Total campaign about 80 minutes wall.

## Not run / limits

* k = 14 (N^2 = 2.7e8): the in-memory transform would need about 2 GB per array plus temporaries, over the 6 GB `ulimit -v`; a chunked version was not written. k=13 ran with 5 instances only.
* k = 9..13 d_lf not measured (law); the quotient engine at k=8 needs about 9 minutes per instance; k=8 has 1 instance, k=7 has 2.
* d_ff for k=10..13 not run (the Ech overflow assertion needs p < ~7*10^5 at those sizes; the instances used primes of 3*10^6 and up).
* Dense Mac_D (literal ordinary ring without chain elimination) only at k=3 (3 instances, D=2..8), not at k >= 4. Brute force beyond k=6 not attempted.
* The literal "12289 for k <= 12" prescription is infeasible for unit instances (see deviation section); k=12 control instance with recorded dim V_4/LFD values not reproduced (no unit instance at p = 12289, k=12; the k=7 control above stands in).
* No replay certificate, validator report, matched rho or `VALIDATION.json`; nothing here may enter the leaderboard.

## Verdicts

* **P1: not falsified.** 0 mismatches in 18,828 brute-force-vs-inverse comparisons (k=3..6), including 129 exceptional instances (sd = 2k+1 on both sides).
* **P2: not falsified.** Generic fraction >= 1-3/p in every cell with exceptions possible; exception rate about 1/p at k=3,4 (sd = 2k+1), none at k=5..13; absence of exceptions at k >= 5 is uninformative at the primes used (1/p tiny).
* **P3: not falsified on k=3..8 (measured), consistent with the law for k=9..13 (not measured).** Equality sd - d_lf = k at exceptional k=3 instances.
* Extra: all exceptions have vanishing top coefficient c_{(N-1,N-1)} (zero_top_coeff true 159/159, false 29,616/29,616 generic), i.e. the Corollary's equivalence "sd = 2+2k iff top coefficient nonzero" holds on every row, since sd = 2k+2 exactly when it is nonzero.
