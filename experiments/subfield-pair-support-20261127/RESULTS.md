# Exact subfield pair-support screen

## Frozen question

Does a subfield-coordinate factor base improve two-summand relation support
and complete single-target cost, compared with an equal-size coordinate
subspace on the **same** non-descending binary curve and target subgroup?
This screen uses exact point addition to verify the relations behind the
binary third summation polynomial. It does not use a polynomial-system
solver, so its pair-table timings are a favorable oracle for relation
finding, not a measured S3/LRF solver time.

## Fixture and validity

`replay.py` fixes `K=F_(2^16)` with modulus
`z^16+z^5+z^3+z^2+1` and the curve
`E: y^2+xy=x^3+z^11*x^2+z`. Its group order is `65986=2*32993`,
with prime `32993`. Since `j(E)=1/z` is outside `F_(2^8)`, this curve
cannot descend to that subfield even after a `K`-isomorphism.

The candidate x-space is the embedded `F_(2^8)`. The control x-space is
`<1,z,...,z^7>`. Both lie in the absolute trace-zero hyperplane and
their nonzero-x points occupy the same non-prime 2-coset. The candidate
has 128 nonzero x-values with two points over each; the control has 129,
and the last sorted x-class is dropped to match **128 unknown factor-base
logs and 256 signed points**. The product-space dimensions are 8 and 15,
respectively. The order-two point at x=0 is excluded from both bases.

Each table entry is an actual equality `P_i+P_j=Q` on `E(K)`.
For every nonzero `Q`, the replay checks the corresponding binary S3
identity

    (x_i*x_j+x_i*x_Q+x_j*x_Q)^2+x_i*x_j*x_Q+z = 0.

It projects the equality into the prime-order subgroup by doubling and
solves the resulting linear system for all base logs and the target log.
The recovered scalar is checked independently by scalar multiplication.

## Results

Exact enumeration over all 32,993 subgroup elements gives:

| Equal-size base | Distinct two-summand sums, including zero | Nonzero query coverage | Maximum multiplicity away from zero |
|---|---:|---:|---:|
| Subfield | 24,605 | 24,604 / 32,992 = 74.58% | 2 |
| Monomial control | 20,959 | 20,958 / 32,992 = 63.52% | 6 |

Both bases have exactly 32,896 unordered signed pairs. The subfield
base therefore has **17.4% greater distinct-query coverage** at fixed
base size and pair count. This is a measured component benefit, not a
lower total solve cost.

The following are cold single-target seconds from the frozen
SageMath 10.10.rc0 replay on this host. Each relation run charges base
construction, all 32,896 pair additions, queries through first full
rank, and linear elimination. The scalar is planted for verification;
the algorithms receive only its public target point. Query counts are
the first full-rank prefix of a deterministic 1,024-query stream.

| Target scalar | Subfield queries / seconds | Control queries / seconds | Cold BSGS seconds | Cold rho seconds |
|---:|---:|---:|---:|---:|
| 3425 | 321 / 1.561 | 298 / 1.163 | 0.00213 | 0.00105 |
| 17123 | 203 / 1.171 | 353 / 1.100 | 0.00185 | 0.01573 |
| 29113 | 257 / 1.170 | 626 / 1.112 | 0.00219 | 0.00165 |

The subfield base raises the chance that one query yields a relation,
but this does not translate to a consistent first-full-rank advantage.
The pair table alone costs roughly one second here, while cold BSGS
uses 182 baby steps and finishes in roughly two milliseconds. Even this
favorable exact-group relation oracle is far behind the matched generic
methods on these three targets. The result is a bounded negative cost
screen for this fixture; it is not a universal complexity bound.

`cold_total_seconds` is optimistic: it uses the known first-full-rank
prefix and excludes the cost of detecting rank during live collection.
Rho timing varies with the deterministic walk seed. All three targets,
query-hit counts, row counts and rank gates remain in `results.json`;
`replay.py` regenerates the individual queries and relation rows.

## Replay

```bash
sage -python experiments/subfield-pair-support-20261127/replay.py \
  > /tmp/subfield-pair-support-replay.json
python3 experiments/subfield-pair-support-20261127/verify.py \
  /tmp/subfield-pair-support-replay.json
```

The verifier compares every discrete field, support, query, rank, and
recovered-scalar record with `results.json`; wall-time fields are omitted
from equality because they vary by host and load. The committed result
SHA-256 is `db52aec9b1236292aff9559c59063568715fb74af5d8604a2a6d9e251a359998`.
