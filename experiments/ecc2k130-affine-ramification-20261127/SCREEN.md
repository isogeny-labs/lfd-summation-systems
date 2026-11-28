# Constant dihedral degree-131 cover: ramification gate

## Exact claim

Under the [frozen hypotheses](PROTOCOL.md), the normal `C_131` layer
`D -> E'=D/C_131` is unramified. Combined with the separately proved
[unramified quotient gate](../ecc2k130-agl2-gate-20261127/SCREEN.md), this
excludes a genus-131 intermediate curve `C=D/C_2` in this **constant
arithmetic-and-geometric dihedral** case. The argument does not exclude
larger affine arithmetic monodromy, a complement of order greater than
two, non-dihedral geometric monodromy, or a Jacobian transfer without a
map `C -> E`.

## Ramification proof

Suppose that the normal `C_131` layer ramifies over a place `v` of `E`.
In the dihedral group `G=C_131 semidirect C_2`, any inertia subgroup
containing a nonidentity element of `C_131` is either `C_131` or `G`.
Wild inertia is a normal 2-subgroup of inertia. The full dihedral group
has no normal order-two subgroup, so `G` cannot be inertia. Thus inertia
is exactly `C_131`, and the ramification is tame.

Let `d=deg(v)`. A lift of arithmetic Frobenius at `v` conjugates a
generator `tau` of tame inertia by `tau -> tau^(2^d)` (up to the harmless
choice of Frobenius convention). Every element of the full arithmetic
group `G` acts on `C_131` by either `+1` or `-1`. Therefore

    2^d = +1 or -1 (mod 131).

The order of `2` modulo `131` is `130`, and `2^65=-1 (mod 131)`. Hence
`65 | d`. In the 131-sheet permutation action for `C/E`, the subgroup
`C_131` is regular. Its fixed subspace has dimension one, so the tame
Artin conductor at `v` is `131-1=130`. By the conductor-discriminant
identity, `v` contributes at least

    deg(v) * 130 >= 65 * 130 = 8450

to `deg Diff(C/E)`. But Riemann-Hurwitz for `g(C)=131`, `g(E)=1` gives

    deg Diff(C/E) = 2g(C)-2 = 260.

All local different contributions are nonnegative, so this is impossible.
The normal `C_131` layer is therefore unramified. Because the assumed
full Galois group acts over `F_2` and `D` is geometrically connected,
`D -> E'` is a geometrically connected unramified constant `C_131`
cover. The prior quotient gate applies and excludes this case.

## Arithmetic control

`python3 experiments/ecc2k130-affine-ramification-20261127/control.py`
checks the order of `2` modulo `131`, its half-order value, the
normal-subgroup fact, the regular-action conductor, and the strict
budget contradiction. This finite computation checks the numerical
inputs, while the written local-field and conductor implications carry
the proof.

## Boundary and next test

If the geometric monodromy is dihedral but the arithmetic group is
larger, arithmetic Frobenius can act on `C_131` by a multiplier other
than `+1` or `-1`; this screen does not apply. The next gate is to
classify those arithmetic normalizers and their smallest possible
branch degrees, then compare their conductor costs with `260`.

The tame Frobenius relation and conductor-discriminant identity used
above are stated in [Local fields, tame ramification](https://openmaths.org/courses/NT-LOC/NT-LOC-08.html)
and [Class field theory, conductors and discriminants](https://openmaths.org/courses/NT-CFT/artin-l-functions-conductors-and-discriminants.html).
