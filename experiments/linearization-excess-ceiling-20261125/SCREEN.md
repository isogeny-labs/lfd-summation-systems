# Linearization-excess ceiling on the frozen LRF bases

## Question

An existing candidate lead proposes choosing a factor-base subspace with
larger **linearization excess** `e = n - k - dim(V²)`, where
`V² = span_F2{uv : u,v in V}`. A fitted two-summand PDP rule associates
larger `e` with smaller solving degree. Could that basis-selection lever
improve the current prime-degree binary LRF campaign at fixed `k`?

## Ceiling and proof

Let `K=F_(2^n)` with `n` prime and let `V` be any nonzero `k`-dimensional
`F_2`-subspace. The finite-field linear Kneser inequality gives

    dim(V²) >= min(n, 2k-1).

One direct source is [Hou, Leung and Xiang, Corollary 2.5](https://bpb-us-w2.wpmucdn.com/sites.udel.edu/dist/a/6387/files/2017/08/Kneser-167tfeu.pdf);
the prime-extension specialization is stated as
[Corollary 3.3 here](https://eccc.weizmann.ac.il/report/2009/063/revision/1/download/).
For clarity, the stabilizer `H={a in K : a V² subset V²}` in the linear
Kneser inequality is a subfield. Prime `n` leaves only `F_2` and `K`.
If `H=K`, then `V²=K` and its dimension is `n`. Otherwise `dim H=1`
and `dim V² >= 2k-1`. Thus

    e <= n-k-min(n,2k-1).

When `2k-1 <= n`, the monomial space
`V=<1,t,...,t^(k-1)>` attains equality: its product space is
`<1,t,...,t^(2k-2)>`, and these powers are independent because `t`
has minimal polynomial of degree `n`. Hence `e_max=n-3k+1`.

## Frozen campaign result

`verify.py` recomputes all pair products with the exact field modulus in
`../subspace-remainder-cost-20261125/manifest.json`, checks an independent
binary rank against `csd.ic_rootfind.product_space`, and writes the frozen
`verification.json`. The manifest SHA-256 is pinned in the verifier.

| (n,k) | `dim(V²)` | `e` | maximal possible `e` at this `k` |
|---|---:|---:|---:|
| (23,8) | 15 | 0 | 0 |
| (23,10) | 19 | -6 | -6 |
| (23,11) | 21 | -9 | -9 |
| (29,12) | 23 | -6 | -6 |
| (29,13) | 25 | -9 | -9 |
| (29,14) | 27 | -12 | -12 |

**Decision.** At each frozen `(n,k)` the existing monomial base already
maximizes `e`. Selecting another `k`-dimensional subspace cannot increase
this predictor or lower the LRF kernel dimension through a smaller
`dim(V²)`. This closes the *maximize-linearization-excess* lever for these
cells. It does **not** prove that every other subspace has the same
three- or four-summand solving degree or total cost: the fitted rule was
measured for a two-summand PDP presentation, and other structural features
could matter. A new basis proposal needs a measured advantage beyond `e`,
with matched target and base controls and the full collection cost charged.

## Composite-degree scope boundary and cost gate

The prime-degree hypothesis matters. If `k` divides a composite `n`, then
`V=F_(2^k)` is a subfield of `K=F_(2^n)` and `V²=V`, so
`dim(V²)=k` and `e=n-2k`. At `n=2k`, its nominal LRF kernel excess is
`k+dim(V²)-n=0`, whereas a monomial `k`-space has product-space
dimension `2k-1` and nominal excess `k-1`. Actual affine-system rank
and kernel dimension still require measurement.
This is a genuine **component-level** opening beyond the prime-degree
ceiling; subfield-coordinate factor bases themselves are established
index-calculus constructions (see [Diem, Sections 1 and 2](https://msp.org/ant/2013/7-6/ant-v7-n6-p01-s.pdf)).

There are two separate gates before claiming a relation advantage.
First, if the curve descends to `F_(2^k)`, the `x`-subfield base has a
strong Frobenius structure. For the `2^k`-Frobenius `pi`, a point has
`x(P) in F_(2^k)` exactly when `pi(P)=P` or `pi(P)=-P`: the two points
above one abscissa are negatives. On the odd-order component, the
projectors `(1+pi)/2` and `(1-pi)/2` split any target into these two
parts. The group order also factors as
`#E(F_((2^k)^2))=(2^k+1-t)(2^k+1+t)`, where
`t=2^k+1-#E(F_(2^k))`. Thus the same structure that helps decomposition
gives a much stronger generic reference on the factors, whose orders
are on the order of `2^k` rather than `2^(2k)`. The factor base is **not**
confined to `E(F_(2^k))`; it also contains the Frobenius `-1` eigenset.
For this case, the matched rho reference must use the target's actual
prime-order factor after projection, not the full curve order.
For a curve with no subfield descent, including through a `K`-isomorphic
model, this projector shortcut is unavailable; target support must still
be measured rather than inferred from a coefficient choice.

The following exact small-field replay checks the minus eigenset and the
order factorization. It finds an `F_16`-abscissa point on `E/F_256`
which is not `F_16`-rational:

```bash
sage -python - <<'PY'
from sage.all import GF, EllipticCurve
K = GF(2**8, 'a')
F = K.subfield(4, 'b')
E0 = EllipticCurve(F, [1, 0, 0, 0, 1])
E = E0.change_ring(K)
t = 17 - E0.cardinality()
assert E0.cardinality() == 16
assert E.cardinality() == 288 == (17-t)*(17+t)
P = next(P for P in E.points() if not P.is_zero()
         and P[0]**16 == P[0] and P[1]**16 != P[1])
assert E([P[0]**16, P[1]**16, 1]) == -P
print('subfield-x minus eigenset verified')
PY
```

Second, at `n=2k` a base with `B` signed points has `B` on the order of
`sqrt(N)` for `N=#E(K)`. On a non-descending curve, if the uniform-sum
heuristic passes an exact target-count check, a two-summand
query has `Theta(B²/N)=Theta(1)` relations. A method that determines
`Theta(B)` factor-base logarithms therefore needs `Omega(B)` such queries
even with a constant-cost relation oracle, and constructing the base
already touches `Omega(B)` points. This matches the generic square-root
*exponent* when the attacked subgroup order is comparable to `N`, before
sparse linear algebra, failed queries, and target
resolution are charged; it does not rule out a constant crossover.
At arity three, one query can return `Theta(B)` relations, but must
output and verify them, and their independence is not guaranteed.

The next useful experiment, if this route is pursued, must freeze a
non-subfield-curve target panel and the exact subfield factor base, compare
the complete signed relations and ranks with a size-matched ordinary base,
and charge base construction, failed queries, relation collection, sparse
elimination, target resolution, and a cold same-target rho reference. A
kernel-dimension or single-query win alone does not pass that gate.
