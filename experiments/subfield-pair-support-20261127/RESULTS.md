# Exact subfield pair-support screen

## Frozen question

Does a subfield-coordinate factor base improve two-summand relation support
and complete single-target cost, compared with an equal-size coordinate
subspace on the **same** non-descending binary curve and target subgroup?
This screen uses exact point addition to verify the relations behind the
binary third summation polynomial. It benchmarks both an exhaustive pair
table and a direct subfield-quadratic S3 oracle. The latter is a specialized
two-summand construction; it does not measure a general S3/LRF solver.

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

## Frobenius-eliminated S3 oracle

Here `q=256` and `K/F_q` is quadratic. For a nonzero query point `Q`, put
`t=x(Q)`, and for candidate abscissas `x1,x2 in F_q`, put
`s=x1+x2` and `p=x1*x2`. The binary S3 equation becomes

    p^2 + t^2*s^2 + t*p + b6 = 0.

The exact base enumeration proves every point with `x in F_q` lies outside
the odd-order subgroup, apart from the excluded order-two point. Hence
every nonzero subgroup query has `t not in F_q`. Write
`A=t+t^q`, `B=b6+b6^q`, `N=t*t^q`, and `C=B/A`. Taking the
`q`-Frobenius difference of S3 yields `p=A*s^2+C`. Substituting
`u=s^2` back into S3 gives the **quadratic over F_q**

    A^2*u^2 + N*u + D = 0,    D=C^2+t*C+b6 in F_q.

For each of its at most two roots, `s` is uniquely determined because
squaring is bijective in `F_q`; `p=A*u+C`. The two abscissas are the roots
of the second quadratic `X^2+s*X+p`. Each unordered abscissa pair gives
at most one signed representation of a fixed nonzero `Q`, because the
base excludes order-two summands. Thus this oracle returns at most two
signed relations per query. It precomputes the Artin-Schreier map
`v -> v^2+v` on `F_q` to solve the two quadratics, and verifies every
returned point relation. A separate exhaustive pair-table audit confirms
exact set equality on all 1,024 queries for each target.

This lowers the **univariate degree of this two-summand relation query**
to two. It is not a bound on the solving degree of other arities or a
generic binary S3 presentation.
Subfield-coordinate factor bases are established index-calculus objects
([Diem, Sections 1–2](https://msp.org/ant/2013/7-6/ant-v7-n6-p01-s.pdf));
this record makes no novelty claim for the construction or elimination.

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
SageMath 10.10.rc0 replay on this host. Pair-table runs charge base
construction, all 32,896 pair additions, queries through first full
rank, and linear elimination. Direct-oracle runs charge base
construction, its Artin-Schreier lookup, queries, and elimination; the
pair table used only to audit completeness is excluded. The scalar is
planted for verification;
the algorithms receive only its public target point. Query counts are
the first full-rank prefix of a deterministic 1,024-query stream.

| Target scalar | Direct subfield queries / seconds | Subfield pair-table seconds | Control pair-table queries / seconds | Cold BSGS seconds | Cold rho seconds |
|---:|---:|---:|---:|---:|---:|
| 3425 | 321 / 1.041 | 1.570 | 298 / 1.270 | 0.00139 | 0.00098 |
| 17123 | 203 / 0.628 | 1.529 | 353 / 1.138 | 0.00184 | 0.02325 |
| 29113 | 257 / 0.718 | 1.270 | 626 / 1.155 | 0.00219 | 0.00168 |

The subfield base raises the chance that one query yields a relation,
but this does not translate to a consistent first-full-rank advantage.
The direct quadratic oracle saves the pair-table work, but its cold
method still takes 0.63–1.04 seconds. Cold BSGS uses 182 baby steps and
finishes in roughly two milliseconds. Even this specialized oracle is far
behind the matched generic methods on these three targets. The result is
a bounded negative cost screen for this fixture, not a universal
complexity bound. The specialization requires a proper degree-eight
subfield of the degree-16 field, so it does not transfer directly to
prime extension degree 131.

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
SHA-256 is `84ec881a8d4480c4514092b11f1d30b17f1a36854b4de1bafbf688cdf330e6a4`.
