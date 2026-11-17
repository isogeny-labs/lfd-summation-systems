# Result: independent lift at k = 14, unit instance over p = 163841

Reading by the rule of `PLAN.md` section 2: **CONFIRM**. Records: `pod-p163841/`. Run of
2026-11-17 on one H200 (launch record and the interruption: `LAUNCH-20261116.md`).

## Preconditions

| check | record | result |
|---|---|---|
| arithmetic on the card, and k = 10, D = 5 -> 6, six runs | first attempt on this pod, 22:29 to 22:35 UTC (`chain-independent-status.json` of that attempt was overwritten by the restart; the six runs are those of the same code and card type recorded under `pod/`) | passed on this pod before the interruption; not repeated in the run that counts |
| closure V_6(14) | `pod-p163841/unit-p163841-k14-D6.json` | closed, rank 448,922, 1 not in; both basis hashes equal to the recorded ones (`e0c4afd7…`, `44a89923…`); 5,802 s |
| independent checker on that stored basis | `pod-p163841/certcheck-unit-p163841-k14-D6.json` | `all_true` |

The run that counts is job `chain3` (`--stage k14`), started after a provider restart of the
container killed the first attempt during the closure. It skipped the small checks.

## The independent lift (`pod-p163841/indep-lift-unit-p163841-k14-D6to7.json`)

`status: done`, `fixed_point: true`, `sigma_stable: true`, `saturated: true`,
`closed_check.nonzero: 0`; 698.1 s; host peak 46.8 GB; basis files verified and canonical
hashes equal to the closure record.

| field | predicted (`PLAN.md`, before any run) | observed |
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

## The first lift on the same files (`pod-p163841/lift-unit-p163841-k14-D6to7-again.json`)

`contains_one: true`, `dim_VD1` 1,683,218, `dim_T` = `dim_K` = [30,655, 30,639], falls
[25,136, 25,120], low part ranks [249,824, 249,354]; 671.2 s.

## Both unit instances

The two instances (p = 65537 and p = 163841) give the same verdict and the same dimensions in
both implementations. The observed fields are identical between the instances, including the
unpredicted `n_single` = 62.

The limit of `RESULT-unit-p65537.md` applies unchanged: both lifts start from a degree-6 space
produced by one engine, and neither run shows that this space is contained in V_6.

## Cost

Pod `hilmllrqee4wp2`, US$5.29 per hour, created 22:28 UTC on 2026-11-16, deleted 03:14:43 UTC
on 2026-11-17 (HTTP 204; a following listing shows none of this experiment's pods): about 4.8
hours, about US$25.3. Of that, about 1.6 hours went to the interrupted first attempt and about
1 hour to the pod sitting idle after the restart until the next check. The stored basis was not
archived.

Both instances together: about US$40.7, against the planned ceiling of US$106.
