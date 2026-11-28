# Unramified affine degree-two quotient: excluded

## Claim boundary

Under the [frozen hypotheses](PROTOCOL.md), there is no genus-131 curve
`C/F_2` with a degree-131 map to
`E: y²+xy=x³+1` whose geometric Galois closure `D/E` has group
`C_131 semidirect C_2` and whose normal `C_131`-cover `D/E'` is
geometrically connected, unramified, and constant over `F_2`.
This excludes one subcase of the non-Galois cover residue in the
[cover-of-E review](../ecc2k130-weil-factor/round2/cover-of-E-review/REPORT.md).
It makes no assertion about a ramified normal cover, a larger complement,
another monodromy group, or a Jacobian without such a map.

## Proof

Let `G=C_131 semidirect C_2` act on its 131 affine points. A reflection
fixes one point and swaps the other 130 in 65 pairs. Since `D/E'` is
unramified, all inertia in `D/E` lies, up to conjugacy, in the order-two
complement. The permutation representation for `C/E`, restricted to this
inertia, is `66*trivial + 65*sign`. The permutation representation for
`E'/E` is `trivial + sign`. Artin conductors therefore give, at every
place, including wild places and with residue degrees counted,

    deg Diff(C/E) = 65 deg Diff(E'/E).

Both `C` and `E'` map to the genus-one curve `E`. Riemann–Hurwitz and
`g(C)=131` give `deg Diff(C/E)=2*131-2=260`. Thus
`deg Diff(E'/E)=260/65=4`, and Riemann–Hurwitz for the degree-two
map gives `2g(E')-2=4`, hence `g(E')=3`.

The four `F_2`-points of `E` are the point at infinity,
`(0,1)`, `(1,0)`, and `(1,1)`. The map `E' -> E` induces a nonzero
norm homomorphism `J(E') -> E`. Equivalently, the Frobenius polynomial of
`E` is a factor of that of `J(E')`, so an `F_2`-isogeny decomposition
`J(E') ~ E x B` exists with `dim B=2`. In particular,

    4 = #E(F_2) divides #J(E')(F_2).

A geometrically connected unramified *constant* `C_131`-cover of `E'`
requires a quotient of `Pic^0(E')(F_2)=J(E')(F_2)` of order 131.
Therefore `131` also divides `#J(E')(F_2)`, forcing
`#J(E')(F_2) >= 4*131 = 524`.

For any genus-three curve over `F_2`, each of its six Frobenius
eigenvalues has complex absolute value `sqrt(2)`. Therefore

    #J(E')(F_2) = product_i (1-alpha_i)
                  <= (1+sqrt(2))^6
                   = 99+70 sqrt(2) < 198.

The last inequality is exact: `70²*2=9800 < 99²=9801`.
The required lower and upper bounds contradict each other.

## Checks and interpretation

`sage -python experiments/ecc2k130-agl2-gate-20261127/control.py` replays
the finite controls. Direct enumeration in SageMath 10.10.rc0 returned precisely the four
listed points on `E(F_2)`. A separate integer check counted one fixed
sheet and 65 transpositions for `x -> -x` on `F_131`, and verified the
strict radical inequality. These checks confirm the small arithmetic
inputs; the exclusion rests on the written conductor, class-field, and
Weil-bound arguments.

```bash
sage -python - <<'PY'
from sage.all import EllipticCurve, GF
E = EllipticCurve(GF(2), [1, 0, 0, 0, 1])
assert sorted(str(P) for P in E.points()) == [
    '(0 : 1 : 0)', '(0 : 1 : 1)', '(1 : 0 : 1)', '(1 : 1 : 1)']
assert E.cardinality() == 4
p = [(-x) % 131 for x in range(131)]
assert sum(p[x] == x for x in range(131)) == 1
assert sum(x < p[x] for x in range(131)) == 65
assert 70**2 * 2 < 99**2
print('four base points; affine involution and strict bound checked')
PY
```

As a positive control for the genus and divisibility calculations, take
the Artin–Schreier double cover `E': w²+w=y` over `E`. The function `y`
has a single pole of order three at infinity, so the cover has different
degree four and genus three. Directly counting the two possible `w` values
over each finite base point gives
`#E'(F_2)=3`, `#E'(F_4)=11`, and `#E'(F_8)=3`. Newton identities then give
the zeta numerator

    1 + 3T² - 2T³ + 6T⁴ + 8T⁶,

and `#J(E')(F_2)=16`, divisible by four as required. This control shows
that the genus-three quotient class is nonempty; the proof excludes the
additional unramified `C_131` layer. The three point counts are
computed by enumerating `E(F_(2^m))` for `m=1,2,3`, adding one point above
infinity and, at each finite point `P`, counting solutions of
`w²+w=y(P)` in `F_(2^m)`.

The initially planned enumeration of genus-three Artin–Schreier models is
unnecessary for this subcase: every possible `E'` already violates the
order bound. The next mathematical gate for the affine route is the
**ramified** `C_131`-cover case and complements `C_d` with `d>2`.
The existing review's degree-ten, genus-eleven example lies outside this
screen. A complete transfer would still require an explicit curve,
a map carrying the actual target subgroup, and a charged comparison with
same-target generic solving.
