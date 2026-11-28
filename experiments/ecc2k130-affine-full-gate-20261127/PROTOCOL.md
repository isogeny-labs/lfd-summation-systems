# Frozen gate: all degree-131 affine complements

## Exact scope

Let `E/F_2` be `y²+xy=x³+1`. Let `C/F_2` be a smooth,
geometrically connected curve of genus 131 with a separable degree-131
map `C -> E`. Assume `J(C)` is `F_2`-isogenous to
`W=Res_(F_(2^131)/F_2)(E)`, and that the geometric Galois closure of
`C -> E` has affine group `C_131 semidirect C_d` in its natural action,
where `d` is a nontrivial divisor of 130. The arithmetic group may be
larger. Test whether any such cover exists.

This is a cover-transfer gate, not a statement about all summation
polynomials or all Jacobian factors.

## Predeclared gates

1. Check `#E(F_2)=4` and the required `#C(F_2)=3` from the precise Weil
   class. Generalize the [dihedral rational-fiber proof in PR #83](https://github.com/isogeny-labs/lfd-summation-systems/pull/83)
   to every `d`: if arithmetic Frobenius has multiplier outside the
   geometric `C_d`, count fixed inertia orbits at each rational base
   place, including all possible inertia groups. Decide whether this
   forces arithmetic and geometric monodromy to coincide.
2. With constant affine group, classify whether inertia containing
   `C_131` can also contain nontrivial complement elements in
   characteristic two. For inertia whose complement has order divisible
   by `ell=5` or `13`, use its affine fixed center, local tame Frobenius,
   and the 131-sheet permutation conductor to bound the degree-weighted
   different. Compare with `2g(C)-2=260`.
3. If no such odd-complement inertia is possible, test whether the
   corresponding `C_ell` quotient of the Galois closure is a geometrically
   connected, constant, unramified cover of `E`. Check its compatibility
   with `#E(F_2)=4`, using either function-field class field theory or
   the elliptic-isogeny argument.
4. Apply the previously proved constant-dihedral gates only if the
   surviving complement is `C_2`. Re-examine the review's conditional
   degree-ten ramification example against the local Frobenius constraint.

If all gates hold, record an exclusion of this exact affine monodromy
route. Otherwise retain the surviving `d` and the failed implication.
Leave non-affine monodromy, map degrees other than 131, and Jacobians not
in the stated Weil isogeny class open.

## Verification

An independent finite script should enumerate the divisors of 130,
the affine inertia orbit counts, `ord_5(2)`, `ord_13(2)`, and minimum
conductor contributions. The local-field and unramified-cover arguments
must be written; small arithmetic does not substitute for them.
