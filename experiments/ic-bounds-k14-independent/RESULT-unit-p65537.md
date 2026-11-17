# Result: independent lift at k = 14, unit instance over p = 65537

Reading by the rule of `PLAN.md` section 2: **CONFIRM**. Records: `pod/`. Run of 2026-11-16 on
one H200 (launch record `LAUNCH-20261116.md`).

## Preconditions

| check | record | result |
|---|---|---|
| arithmetic on the card | `pod/indep-arith-gpu.json` | passed |
| k = 10, D = 5 -> 6, unit and proper, three settings each | `pod/indep-lift-*-k10-D5to6-*.json` | six runs equal to the predicted values |
| closure V_6(14) | `pod/unit-p65537-k14-D6.json` | closed, rank 448,922, 1 not in; both basis hashes equal to the recorded ones (`eeac52ef…`, `a58c6c7d…`); 6,707 s |
| independent checker on that stored basis | `pod/certcheck-unit-p65537-k14-D6.json` | `all_true` |

## The independent lift (`pod/indep-lift-unit-p65537-k14-D6to7.json`)

`status: done`, `fixed_point: true`, `sigma_stable: true`, `saturated: true`,
`closed_check.nonzero: 0`; 699.6 s; host peak 46.7 GB. The driver verified every file of the
stored basis and recomputed its canonical hashes before lifting.

| field | predicted (`PLAN.md`, before the run) | observed |
|---|---|---|
| `contains_one` | true | true |
| `dim_VD` | 448,922 | 448,922 |
| `dim_low` = dim(V_7 n M_<=6) | 499,178 | 499,178 |
| `dim_low_parts` | [249,824, 249,354] | [249,824, 249,354] |
| `dim_VD1` = dim V_7 | 1,683,218 | 1,683,218 |
| `n_special_tops` | 11,038 | 11,038 |
| `n_regular_tops` | 1,173,002 | 1,173,002 |
| `n_N_weight_D` | 1,406 | 1,406 |
| `n_T - n_single` | 61,294 | 61,356 - 62 = 61,294 |
| `rank_K - n_single` | 61,294 | 61,356 - 62 = 61,294 |
| `falls_new` | 50,256 | 50,256 |
| `n_single` | not predicted | 62 |

The space saturated (K = T) after the first chunk of regular tops and three rounds of falls
(`counts`: 1 chunk, 3 fall rounds).

## The first lift on the same files (`pod/lift-unit-p65537-k14-D6to7-again.json`)

`contains_one: true`, `dim_VD1` 1,683,218, `dim_T` = `dim_K` = [30,655, 30,639], falls
[25,136, 25,120], low part ranks [249,824, 249,354]; 673.4 s. Equal to its record of
2026-11-10/11 and to the independent lift.

## What this shows and what it does not

Two implementations that share no code, in different coordinates, give the same verdict and the
same dimensions for the step from the stored degree-6 space to degree 7: 1 is in the lifted
space and the lifted space is all of M_<=7.

It does not show that the stored degree-6 space is contained in V_6 (`PLAN.md` section 8, item 1;
the manuscript's limitations state this). That inclusion still rests on the degree-6 engine, here
reproduced a third time with the same hashes.

## Cost

Pod `yu1co05fj8x398`, US$5.29 per hour, created 19:33 UTC, deleted 22:27:36 UTC (HTTP 204, a
following listing no longer shows it): about 2.9 hours, about US$15.4. The chain itself ended at
22:14 UTC; the pod idled about 13 minutes until the next scheduled check. The stored basis
(about 45 GB) was not archived and went with the pod.
