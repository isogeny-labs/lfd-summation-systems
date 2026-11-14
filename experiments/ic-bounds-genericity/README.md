# Genericity of the depth-7 base (N6): cheap experiments

Companion to `docs/proofs/N6-genericity.md`. Everything here is **exact** (closures verified
closed by re-reducing all shifts of the low rows, or dimensions/pivots compared between closed
bases) unless stated otherwise. Nothing is kernel-checked. Outputs total under 0.3 MB.

Model: the group algebra `A_k^{(c)} = F_p[x1, x2]/(x1^(2^k) - c, x2^(2^k) - c)`, monomials =
exponent pairs, degree = Hamming weight, chain variable `u_{i,j}` = translation by `2^j` in block
`i`, with the factor `c` on wrap-around. `c = 1` is the chain algebra; `c = mu^(-2^k)` is the
coset `mu^(-1) H` as factor base, equivalently the parameters `(mu^2 a, mu^3 b, mu t)`
(N6 document, Proposition 3.1). `c = 0` is the no-wrap truncation.

Engine: `genericity_engine.py` (grid and `fp_echelon` from `csd/`; the twisted shift is new).
It reproduces the item-5 rank histories `[79, 123, 212, 272, 363, 422, 491, 557, 598, 605]`
(k = 6) and `[106, 160, 282, 360, 483, 566, 661, 759, 818, 829]` (k = 7) of
`experiments/ic-bounds-fall/p40961-s7/` entry by entry.

## Commands (repository root; every run under 20 min, peak RSS < 0.5 GB)

    export PYTHONPATH=$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ulimit -v 3000000
    python3 experiments/ic-bounds-genericity/twist_test.py            # twist_test.json      (~4 min)
    python3 experiments/ic-bounds-genericity/coset_test.py            # coset_test.json      (~9 min)
    python3 experiments/ic-bounds-genericity/extension_twist_test.py  # extension_twist_test.json (~2 min)
    python3 experiments/ic-bounds-genericity/line_degree.py --k 6 --N 1200          # line_degree_k6.json (12 min + 9 s)
    python3 experiments/ic-bounds-genericity/line_degree_chunks.py sample --k 7 --s-start 1 --s-end 700      # ~17 min, scratch npy
    python3 experiments/ic-bounds-genericity/line_degree_chunks.py sample --k 7 --s-start 701 --s-end 1400   # ~17 min
    python3 experiments/ic-bounds-genericity/line_degree_chunks.py reconstruct --k 7                         # line_degree_k7.json

The line scripts keep their sampled rows (int32, about 250 MB per 700 points at k = 7) in the
session scratch directory, not in the repository.

## Results

### twist_test.json: c-dependence of V_4 in A_k^{(c)}

Instances: `p40961-s7` (a=21223, b=9887, t=4748), `p12289-r3` (3899, 9710, 9895),
`p12289-j1728` (1, 0, 9362), and at p = 193: generic (22, 23, 1), order two (22, 170, 1),
tangent (190, 23, 1). For each, `c ∈ {1, 0, 2, 3, p-1, four random}`.

* **k = 7:** for every instance and every c the closure is the **same subspace** as for c = 1
  (equal row spaces), dimension 829 (835 order two, 828 tangent), same profile, closed, 1 not in.
* **k = 6:** for every instance and every c the dimension (605 / 607 / 603), profile and pivot set
  equal those of c = 1, but for c != 1 the subspace differs (neither inclusion holds).

### coset_test.json: all cosets of H at p = 12289

Instance `p12289-r3`: `(mu^2 a, mu^3 b, mu t)` for one mu in each coset of H.
k = 6 (H of order 64, 192 cosets): all 192 give dim 605, profile (0,0,2,111,605), the pivot set
of mu = 1, closed. k = 7 (96 cosets): all 96 give 829, (0,0,2,115,829), same pivots.
(Item 2 could test only the 3 cosets that exist at p = 193.)

### extension_twist_test.json

* Free extension W_8 of V_4(7) on `p40961-s7` (1061 rows): closed under the twisted shifts for
  c in {1, 0, 2, 777}, 1 not in W_8, and equal to the twisted closure V_4^{(c)}(8) for each c.
* Pivot sets of V_4(6) and V_4(7) at random points `(A, B, 1)` over p = 193, 12289, 40961,
  33554393 (two points each): identical to the reference (symmetric difference 0), dims 605 / 829.
* Density of the RREF basis at a random point: k = 6: 68,803 nonzero of 605 × 794 (14.3 %),
  8,841 in the 111 fall rows; k = 7: 93,011 of 829 × 1471 (7.6 %), 8,893 in the 115 fall rows.
  The counts agree at p = 40961 and p = 33554393; a third random point (p = 12289) gives 68,801
  at k = 6 and 93,011 at k = 7, and the nonzero patterns are not identical: a few numerators
  vanish at individual points, so the support is a property of the generic basis, not of a point.

### line_degree_k6.json: degree of the pivot denominator along a line (k = 6)

Line `(A, B) = (A0 + s A1, B0 + s B1)`, t = 1, p = 33554393, s = 1..1200: the pivot set is the
generic one at every s. 60 rows stored (2 of weight 2, 28 of weight 3, 30 of weight 4), 6
nonzero entries sampled per row, each reconstructed as a reduced rational function in s
(Newton interpolation, extended Euclid, degree budget 1200 - 12).

* 360 / 360 reconstructed; max numerator degree 538; max denominator degree 480; **lcm of
  denominators: degree 480**.
* By pivot weight: weight 2: denominators of degree <= 1 (`A + 3`); weight 3: quartiles
  7 / 18 / 79, max 115; weight 4: quartiles 30 / 82 / 361, max 480 (21 of 180 entries >= 400).
* The reconstructed lcm vanishes at the s where the line meets `1 + A + B = 0` and `A = -3`,
  as the specialisation theorem predicts (both are exceptional: larger closure, resp. the
  denominator of the weight-2 fall element w).

Label: exact lower bound `deg D >= 480` for the restriction to this line (an unsampled entry
could have a larger denominator); for a generic line this is the total degree of D(A, B).

### line_degree_k7.json (k = 7)

Same line (seed 1, so the same A0, B0, A1, B1), p = 33554393, s = 1..1400 in two 700-point
chunks (916 s and 926 s), generic pivot set at every s, 60 rows stored, 360 entries reconstructed
(degree budget 1400 - 12).

* 360 / 360 reconstructed; max numerator degree 533; max denominator degree 480; **lcm of
  denominators: degree 480**, the same as at k = 6.
* By pivot weight: weight 2: degree <= 1; weight 3: quartiles 4 / 18 / 78, max 111; weight 4:
  quartiles 15 / 73 / 99, max 480.
* The lcm vanishes at the order-two and tangent intersection points of the line.

Label: exact lower bound `deg D >= 480` on this line at k = 7. Whether the k = 6 and k = 7 lcm
polynomials coincide was not checked (only degrees were recorded).

### Checks of the two counterexamples quoted in the document (one-off, Section 3.2)

    python3 - <<'EOF'
    import sys; sys.path.insert(0,'experiments/ic-bounds-genericity'); import genericity_engine as G
    p,k=193,6; tg=G.TwistGrid(k,4)
    for a in (13, 52):   # 13 in H (64th powers), 52 = mu^2 * 13 with mu^2 not in H
        r=G.closure(G.generator(a,0,0,p,k),tg,p,1); print(a, r['rank'], r['contains_one'], r['profile'])
    EOF

gives `13: 730, False, [0,1,35,235,730]` and `52: 794, True, [1,13,79,299,794]`. The C2 example
`(A,B) = (20, 0)`: `plane/p193-k6-B0.jsonl` (t = 1: |Z| = 10, V4 = 784) against
`data/p193-k6-a5-b0.jsonl`, t = 96 (|Z| = 32, V4 = 762), both saturated.
