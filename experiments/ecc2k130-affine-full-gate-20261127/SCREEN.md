# All affine degree-131 complements: excluded for the Weil class

## Exact claim and boundary

Under the [frozen hypotheses](PROTOCOL.md), no genus-131 degree-131
cover `C -> E` with geometric affine monodromy
`C_131 semidirect C_d`, `d>1`, has `J(C)` in the isogeny class
`W=Res_(F_(2^131)/F_2)(E)`. This closes the **affine degree-131**
part of the cover route to that precise Weil factor. It does not treat
`A_131` or `S_131` monodromy, maps of other degree, or Jacobians outside
this isogeny class. It does not produce a relation-finding speedup.

## 1. Arithmetic monodromy must equal geometric monodromy

The geometric group `G=C_131 semidirect C_d` lies in `AGL(1,131)`;
`d` divides `130`. Its regular translation subgroup is characteristic
in `G`, so the arithmetic normalizer also lies in `AGL(1,131)`.
Write `H=C_d` for the geometric multiplier group. The arithmetic
Frobenius multiplier determines a coset `uH`.

Direct enumeration gives `#E(F_2)=4`. Frobenius on the 131 geometric
copies making up `W` is a 131-cycle, so its degree-one trace is zero.
Thus `J(C) ~ W` requires `#C(F_2)=2+1-0=3`.

Suppose `u` is outside `H`. At any rational base place, local
Frobenius has multiplier `a` in `uH`, hence `a` is outside `H`.
Inertia in the geometric affine group is either the regular
`C_131` or a subgroup of a conjugate of `H`:

- A subgroup containing `C_131` and a nontrivial complement is
  noncyclic. If its order is odd, tame inertia cannot have that form.
  If its order is even, it has no nontrivial normal 2-subgroup, so it
  cannot have wild inertia either.
- A nontrivial subgroup of order prime to 131 fixes one affine center `c` and is
  conjugate into `H`.

For trivial inertia, `x -> ax+b` has one fixed sheet because `a != 1`.
For `C_131` inertia, its 131 sheets form one orbit. For inertia
`I<=H` fixing `c`, the inertia orbits are `{c}` and cosets of `I`
on the nonzero offsets. The local Frobenius normalizes `I`, hence fixes
its unique center and acts on offsets by multiplication by `a`.
A nonzero orbit could be fixed only if `a in I`, impossible because
`a` is outside `H`. Rational points above a rational base place are
the Frobenius-fixed inertia orbits, so every one of the four rational
base points has exactly one rational preimage. This gives
`#C(F_2)=4`, a contradiction. Hence `u in H` and the full arithmetic
group equals `G`; its Galois closure has constant field `F_2`.

This extends the [degree-two fiber count in PR #83](https://github.com/isogeny-labs/lfd-summation-systems/pull/83)
to every affine complement. The proof here is self-contained so that
the conclusion does not depend on the PR's merge order.

## 2. Odd complement inertia exceeds the different budget

Riemann-Hurwitz for `g(C)=131`, `g(E)=1` gives
`deg Diff(C/E)=260`. Let `ell` be `5` or `13`, and suppose some
inertia has complement order divisible by `ell`. By the classification
above, it is a subgroup of an affine complement, with a unique fixed
center. Every element of its decomposition-group normalizer fixes
that center and is a scalar map about it, so it **centralizes** the
order-`ell` tame inertia. The local tame Frobenius relation therefore
forces `2^deg(v)=1 (mod ell)`.

| `ell` | `ord_ell(2)` | Minimum local conductor on 131 sheets | Minimum degree-weighted contribution |
|---:|---:|---:|---:|
| 5 | 4 | `130(1-1/5)=104` | `4*104=416` |
| 13 | 12 | `130(1-1/13)=120` | `12*120=1440` |

The conductor lower bounds are codimensions of the inertia-fixed
subspace of the permutation representation. Wild terms, if present,
are nonnegative. Either minimum already exceeds the *entire* different
budget of 260. Thus no inertia can have a 5- or 13-component.

## 3. Unramified odd quotients are impossible

If `ell=5` or `13` divides `d`, project the constant affine group to
its `C_ell` multiplier quotient. Every inertia group maps trivially to
this quotient: translation inertia maps trivially, and complement
inertia has no `ell`-component by section 2. The quotient is therefore
a geometrically connected, constant, unramified cyclic `ell`-cover
`E_ell -> E` over `F_2`.

Riemann-Hurwitz makes `E_ell` a genus-one curve. The Hasse lower bound
`3-2 sqrt(2)>0` supplies an `F_2`-point. After translating origins,
the cover is a separable isogeny to `E`; hence
`#E_ell(F_2)=#E(F_2)=4`. Its constant deck group acts by rational
translations in an order-`ell` kernel, forcing `ell` to divide 4.
Neither 5 nor 13 does. This contradiction excludes every `d`
divisible by 5 or 13.

The nontrivial divisors of 130 are `2,5,10,13,26,65,130`, so only
`d=2` remains. That constant dihedral case is excluded by the
[ramification](../ecc2k130-affine-ramification-20261127/SCREEN.md)
and [unramified quotient](../ecc2k130-agl2-gate-20261127/SCREEN.md)
gates. These hypotheses follow from the constant geometric Galois
closure established in section 1.

## Correction to the earlier degree-ten example

The [cover review](../ecc2k130-weil-factor/round2/cover-of-E-review/REPORT.md)
called a degree-ten configuration locally consistent, with rational
order-five and order-ten inertia. In a constant affine closure, a
decomposition-group normalizer centralizes the order-five inertia,
whereas tame Frobenius at a rational place must square it. Since
`2 != 1 (mod 5)`, those rational branch data cannot occur. Moving
such inertia to a degree-four place already costs at least 416,
exceeding 260. This corrects that specific configuration; the proof
above excludes all affine `d`, not merely the example.

## Reproduction

`python3 experiments/ecc2k130-affine-full-gate-20261127/control.py`
checks all divisors of 130, all outside-complement multipliers and
their inertia-orbit counts, the orders of 2 modulo 5 and 13, the
conductor bounds, and the four rational base points. These checks
verify finite arithmetic. The local Frobenius relation, ramification
classification, and elliptic quotient argument are written above.

The local tame Frobenius relation is also stated in
[Local fields: tame ramification](https://openmaths.org/courses/NT-LOC/NT-LOC-08.html).
