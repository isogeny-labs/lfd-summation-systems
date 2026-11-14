# ic-bounds-block4: tests of block-of-four lemmas (conjectural)

Experiments for `docs/proofs/N9-block-four.md` (2026-11-12). Exact arithmetic modulo p in the chain
quotient A_{k1,k2} = F_p[Z/2^k1 x Z/2^k2] (weight = binary digit sum). Every closure is followed by a
deterministic full sweep (`sweep_failures` = number of (low row, shift) products that do not reduce to zero;
0 in every record of this directory). Instances: unit instance of `experiments/ic-bounds-k14`
(y^2 = x^3 + x over F_65537, t = 4271) and, for `*-unit2`, y^2 = x^3 + 3x over F_163841, t = 73730. The
solution count on mu_{2^k1} x mu_{2^k2} is recomputed per record (`solutions`, 0 = unit ideal).

Engines: `scripts/ic_bounds_block4_lib.py` (unequal-depth layout; fast closure on the block-RREF class `Part`
of `csd/ic_bounds_d5_verify.py`, imported unchanged; slow closure from `scripts/ic_bounds_growth_lib.py`,
imported unchanged, used as a cross-check).

Run from `scripts/` with `ulimit -v 1500000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` (the D = 5 rows with
k1 >= 9 under `flock /tmp/lfd-locks/heavy.lock`, `ulimit -v 7000000`, two threads).

| file | command | content | time |
|---|---|---|---|
| `mixed-D3-unit.jsonl` | `ic_bounds_block4_mixed.py OUT 3 1 8 14 2000 --slow-check` | Test M, D = 3, 84 pairs, both engines (`engines_agree`) | 2 min |
| `mixed-D4-unit.jsonl` | `ic_bounds_block4_mixed.py OUT 4 1 8 15 9000` | Test M, D = 4, 89 pairs | 9 min |
| `mixed-D5-unit.jsonl` | `ic_bounds_block4_mixed.py OUT 5 5 8 12 9000`, then `... 5 9 9 11 30000` and `... 5 10 10 11 30000` (heavy lock) | Test M, D = 5 | see records |
| `mixed-D3-unit2.jsonl`, `mixed-D4-unit2.jsonl` | `... 3 1 4 12 3000 unit2`, `... 4 4 7 12 5000 unit2` | Test M on the second instance | 3 min |
| `doubling.jsonl` | `ic_bounds_block4_maps.py doubling OUT D k rmax` for D = 3, k = 2..5, rmax = 4 and D = 4, k = 4..9, rmax = 5 | Test A: g^(r) in V_D(k)? dim of weight-<=1 multipliers Q with Q g^(r) in V_D(k) | 1 min |
| `curve.jsonl` | `ic_bounds_block4_maps.py curve OUT 6 200` | Test A': S_3(x^(2^r), y^(2^r), t) on 200 points of C(F_p); evaluation rank of weight-<=1 functions | seconds |
| `coset.jsonl` | `ic_bounds_block4_maps.py coset OUT D k1 k2 r` for (3,2,6,4) (4,2,6,4) (3,1,6,4) (4,5,7,1) (4,6,7,1) (4,5,9,4) (4,6,8,4) (3,1,8,4) | Test C: coset instances g(s y, x2), all 2^r cosets | 1 min |
| `norm.jsonl` | `ic_bounds_block4_maps.py norm OUT D k1 6 r`, r = 1..4, D = 3..6, k1 = 1..6 (cells with \|M_<=D\| <= 5000) | Test N: closure of the norm h_r = prod_zeta g(zeta x1, x2) | 2 min |
| `normentry.jsonl` | `ic_bounds_block4_maps.py normentry OUT Dmax k1 k2 r` for (4,6,6,1) (4,6,6,2) (4,7,7,1) (4,7,7,2) (4,6,8,2) (3,4,4,1) (3,6,6,1) | Test N: least E with the lifted norm in V_E[g] | seconds |

Labels: every number is **exact** for the instance and depths of its record; the laws read off them
(min law, coordinatewise free extension, Conjecture N) are **empirical**.
