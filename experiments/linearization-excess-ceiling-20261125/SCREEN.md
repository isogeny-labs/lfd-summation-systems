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
