# Frozen rational-fiber count gate for geometric dihedral monodromy

## Exact question

Let `E/F_2` be `y²+xy=x³+1`. Let `C/F_2` be a smooth,
geometrically connected genus-131 curve with a separable degree-131 map
`C -> E`. Suppose the geometric monodromy of its Galois closure is the
nontrivial `C_131 semidirect C_2` in the natural 131-point action.
The arithmetic Galois group may be larger than the geometric group.
Can `J(C)` have the Weil-restriction isogeny class
`Res_(F_(2^131)/F_2)(E) ~ E x A`?

## Predeclared gates

1. Derive the required `#C(F_2)` from the Frobenius trace of this exact
   Weil-restriction class. Check `#E(F_2)` directly.
2. Show the arithmetic group normalizes the geometric dihedral group,
   hence embeds in `AGL(1,131)`; write its Frobenius multiplier `u`.
3. For `u` outside `{+1,-1}`, classify possible inertia groups at every
   rational base place. Count Frobenius-fixed inertia orbits in the
   131-point action, including wild order-two inertia and tame
   order-131 inertia. Compare the total `#C(F_2)` with gate 1.
4. If `u` lies in `{+1,-1}`, check whether arithmetic and geometric
   monodromy coincide; then apply the independently proved constant
   dihedral [ramification](../ecc2k130-affine-ramification-20261127/SCREEN.md)
   and [unramified](../ecc2k130-agl2-gate-20261127/SCREEN.md) gates,
   preserving all their hypotheses.

If all implications hold, record an exclusion only for this geometric
dihedral, degree-131, genus-131 Weil-restriction target class. Otherwise
retain the failing case. Do not extend to complements of order greater
than two, other monodromy, or maps of other degree.

## Verification

Use a written permutation and local-fiber proof. An exhaustive finite
check over affine multipliers and reflection centers may verify the
permutation counts, and direct enumeration checks the four base points.
Finite checks do not replace the monodromy and Frobenius argument.
