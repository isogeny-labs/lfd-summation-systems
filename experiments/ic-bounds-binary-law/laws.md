| law | systems | cells (singleton) | satisfy | violate |
|---|---|---|---|---|
| `|Z|` lookup, all_nonconst | 1775 | 183 (76) | 1335 | 440 |
| `|Z|` lookup, affine | 1439 | 162 (74) | 1061 | 378 |
| `|Z|` lookup, order_two | 336 | 91 (50) | 334 | 2 |
| `(k,|Z|)` lookup, all_nonconst | 1775 | 384 (173) | 1495 | 280 |
| `(k,|Z|)` lookup, affine | 1439 | 339 (165) | 1268 | 171 |
| `(k,|Z|)` lookup, order_two | 336 | 165 (103) | 335 | 1 |
| `n-k` lookup, all_nonconst | 1775 | 10 (0) | 984 | 791 |
| `n-k` lookup, affine | 1439 | 10 (0) | 800 | 639 |
| `n-k` lookup, order_two | 336 | 9 (0) | 280 | 56 |
| `2k-n` lookup, all_nonconst | 1775 | 18 (0) | 1429 | 346 |
| `2k-n` lookup, affine | 1439 | 18 (0) | 1214 | 225 |
| `2k-n` lookup, order_two | 336 | 16 (1) | 325 | 11 |
| `(n,k)` lookup, all_nonconst | 1775 | 54 (0) | 1509 | 266 |
| `(n,k)` lookup, affine | 1439 | 54 (0) | 1307 | 132 |
| `(n,k)` lookup, order_two | 336 | 47 (2) | 325 | 11 |
| `(n,k,class)` lookup, all_nonconst | 1775 | 101 (2) | 1632 | 143 |
| `(n,k,class)` lookup, affine | 1439 | 54 (0) | 1307 | 132 |
| `(n,k,class)` lookup, order_two | 336 | 47 (2) | 325 | 11 |
| `(n,k,class,base_kind,curve)` lookup, all_nonconst | 1775 | 314 (35) | 1686 | 89 |
| `(n,k,class,base_kind,curve)` lookup, affine | 1439 | 177 (7) | 1355 | 84 |
| `(n,k,class,base_kind,curve)` lookup, order_two | 336 | 137 (28) | 331 | 5 |
| `e,class` lookup, all_nonconst | 1775 | 27 (0) | 1534 | 241 |
| `e,class` lookup, affine | 1439 | 15 (0) | 1208 | 231 |
| `e,class` lookup, order_two | 336 | 12 (0) | 326 | 10 |
| `k+e,class` lookup, all_nonconst | 1775 | 39 (0) | 1633 | 142 |
| `k+e,class` lookup, affine | 1439 | 21 (0) | 1315 | 124 |
| `k+e,class` lookup, order_two | 336 | 18 (0) | 318 | 18 |
| `(k,e,class)` lookup, all_nonconst | 1775 | 98 (3) | 1661 | 114 |
| `(k,e,class)` lookup, affine | 1439 | 53 (1) | 1335 | 104 |
| `(k,e,class)` lookup, order_two | 336 | 45 (2) | 326 | 10 |

| closed-form law | systems | satisfy | violate | LFD below / above prediction |
|---|---|---|---|---|
| B: affine, LFD = 2+[k+e>=3]+[k+e>=9]+[k+e>=19] | 1439 | 1315 | 124 | 93 / 31 |
| A: affine, LFD = 2/3/4/5 for e<=-1 / <=3 / <=8 / >=9 | 1439 | 1208 | 231 | 29 / 202 |
| B on order-two (same thresholds) | 336 | 200 | 136 | 131 / 5 |
| semi-regular crossover min{D: sum_{i<=D} c_i(r_eff,m_eff) <= |Z|} | 1775 | 35 | 1740 | 256 / 1484 |
| semi-regular d_reg(r_eff, m_eff) | 1775 | 598 | 1177 | 1173 / 4 |
| k=10 grid rule: LFD 5 iff |Z|>470 else <=4 (proper, k>=8; equality to 5 / <=4 checked) | 489 | 489 | 0 | - / - |

| bound | systems | satisfy | violate |
|---|---|---|---|
| LFD <= min(k,l)+2 = k+2 | 1775 | 1775 | 0 |
| LFD <= semi-regular d_reg(r_eff, m_eff) | 1775 | 1771 | 4 |
| LFD <= max(2, floor(log2|Z|)+1) (proper) | 1486 | 1197 | 289 |
| LFD <= max(2, floor(log2|Z|)+2) (proper) | 1486 | 1458 | 28 |
| LFD <= 4 when |Z| <= 470 (proper) | 1458 | 1458 | 0 |
| LFD <= 4 when k+e <= 18 (nonconst) | 1747 | 1747 | 0 |
| LFD <= 3 when k+e <= 7 (affine) | 690 | 690 | 0 |
