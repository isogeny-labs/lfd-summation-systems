# ECC2K-130 fixed-base two-torsion recipe screen — 2026-11-25

## Fixed sources and scope

Read-only inputs from the `isogeny-labs/crypto` checkout:

- `docs/ic/params/ecc2k130-fixed.json`, SHA-256
  `f41972681d260ef4eeb2fea7441f31859f4a121fdd5f2ecc4ff3dee54082aeae`:
  `K_0/F_(2^131)`, cofactor 4, and a Hamming-weight-two factor-base recipe.
- `ecc2k130/codegen/indexcalc_fixed.py` at the inspected crypto commit:
  `subgroupBase` tests `toCoords(x)` in the type-II normal basis, even though
  public points in the JSON are stated in polynomial coordinates.
- `ecc2k130/codegen/field.py` at that commit: ONB basis elements are
  `gamma_i=zeta^i+zeta^-i`, with `zeta` of order 263, represented as symmetric
  vectors modulo the all-ones vector.

The output is a **recipe-level necessary compatibility test**, not a run of
the exact subgroup point-selection algorithm. Enumerate every nonzero ONB
coordinate of Hamming weight at most two: 131 singleton coordinates and
`C(131,2)=8515` pairs. For each, compute `x -> 1/x`, the abscissa map under
translation by the rational order-two point `(0,1)`, in the exact type-II ONB.
Count how many reciprocals still have weight at most two, by input weight
and overall. Verify multiplication, identity and sampled reciprocal products;
batch inversion is allowed if the result is checked.

This test can show that the weight-two **recipe** is not torsion-stable. It
cannot by itself show whether the *selected subgroup factor-base point set*
is closed, because some recipe candidates may fail to lift or may be
excluded by the subgroup condition. It does not determine solver degree,
useful relation yield or a complete cost ratio. Keep the native
`koblitz-symmetrised` pipeline and its archived measurements distinct from
this fixed ECC2K-130 base.

Also verify `ord_131(2)=130`. The additive frame `u=1/(x+1)` needs a proper
Frobenius-invariant `F_2`-subspace containing 1; primitive Frobenius action
on `F_(2^131)/F_2` makes the only such proper subspace `F_2` itself. This is
a structural size obstruction for directly transplanting the native
subspace construction to ECC2K-130, separate from the weight-two recipe
screen.
