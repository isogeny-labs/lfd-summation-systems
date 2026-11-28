# Geometric dihedral degree-131 cover: rational-fiber exclusion

## Claim boundary

Under the [frozen hypotheses](PROTOCOL.md), a genus-131 curve `C/F_2`
with a degree-131 map to `E: y²+xy=x³+1` and geometric monodromy
`C_131 semidirect C_2` cannot have `J(C)` in the isogeny class
`Res_(F_(2^131)/F_2)(E) ~ E x A`. This closes the geometric **dihedral
degree-131** subcase of the cover route to that precise Weil factor.
It says nothing about larger geometric affine complements, alternating
or symmetric monodromy, other map degrees, or other Jacobian classes.

## Required point count

Direct enumeration gives `E(F_2)={O,(0,1),(1,0),(1,1)}`, of size four.
The Weil restriction `W=Res_(F_(2^131)/F_2)(E)` has Frobenius acting
cyclically on 131 geometric copies of `E`. Its degree-one Frobenius
trace is zero. If `J(C)` is isogenous over `F_2` to `W`, then

    #C(F_2) = 2+1-Tr(Frob | H¹(C)) = 3.

Equivalently, the existing [cover review](../ecc2k130-weil-factor/round2/cover-of-E-review/REPORT.md)
computes `Tr(E)=-1`, `Tr(A)=1`, hence `Tr(E x A)=0`.

## Arithmetic group and local fibers

The normalizer in `S_131` of the geometric dihedral group is
`AGL(1,131)`. Indeed, its translation subgroup `C_131` is the unique
Sylow-131 subgroup, so every normalizer preserves it; the normalizer of
that regular cyclic subgroup is `AGL(1,131)`. Conversely, every affine
map conjugates translations to translations and reflections to
reflections. Write the multiplier of arithmetic Frobenius as `u`.
The geometric multiplier subgroup is `{+1,-1}`.

Suppose `u` lies outside `{+1,-1}`. At any rational place `v` of `E`,
the multiplier of a local Frobenius lift is `u` or `-u`, so it also lies
outside `{+1,-1}`. Inertia lies in the geometric dihedral group because
constant-field extensions are unramified. Its only possibilities are
`1`, a reflection `C_2`, or the translation group `C_131`. Full dihedral
inertia is impossible: wild inertia would be a normal 2-subgroup, which
the dihedral group lacks, while a wholly tame inertia group is cyclic.

Rational points of `C` above `v` correspond to inertia orbits on the
131 sheets that are fixed by local Frobenius:

- If inertia is trivial, the affine Frobenius `x -> ax+b` has exactly one
  fixed sheet because `a != 1`.
- If inertia is `C_131`, its regular action has one orbit, hence one
  rational point above `v`.
- If inertia is a reflection `x -> 2c-x`, its orbits are `{c}` and the
  65 pairs `{c+a,c-a}`. A Frobenius lift normalizes the reflection and
  therefore fixes its center `c`; it acts as `x -> c+a(x-c)` with
  multiplier `a` outside `{+1,-1}`. It fixes `{c}` and no pair, since
  `a z = +z` or `-z` for nonzero `z` would force `a=+1` or `-1`.

Each of the four rational base points therefore has exactly one
rational preimage. Thus `#C(F_2)=4`, contradicting the required count
`3`.

If instead `u` lies in `{+1,-1}`, arithmetic Frobenius adds no multiplier
beyond the geometric group. The full arithmetic and geometric groups
coincide, and the Galois closure has constant field `F_2`. The
[ramification gate](../ecc2k130-affine-ramification-20261127/SCREEN.md)
then makes its normal `C_131` layer unramified, and the
[genus-three quotient gate](../ecc2k130-agl2-gate-20261127/SCREEN.md)
excludes that case. Both gates' geometric connectedness and constant
deck-group conditions follow here from the full constant dihedral
Galois closure. The two arithmetic cases are exhaustive.

## Reproduction and limits

`python3 experiments/ecc2k130-dihedral-fiber-count-20261127/control.py`
enumerates the four rational points of `E` and checks every affine
multiplier outside `{+1,-1}` and every possible reflection center on
131 sheets. It independently checks that the 131-cycle has trace zero
as a permutation matrix and that the constant dihedral arithmetic
alternative is the only omitted multiplier case. These finite checks
verify the permutation arithmetic; the local-fiber correspondence and
isogeny-class trace are the mathematical proof.

This is an exclusion of a proposed transfer construction. It does not
measure a summation-polynomial solver, relation yield, or total ECDLP
cost. The larger affine-complement residue, including the degree-ten
conditional configuration in the review, remains open.
