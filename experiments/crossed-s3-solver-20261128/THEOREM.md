# Projection of O-aware binary S3 trees

Let `K` be a finite field of characteristic two and
`E/K: y^2+xy=x^3+a*x^2+b`, with `b != 0`. Let `X` be any set of
abscissas of `K`-rational points. An auxiliary label is either an element of
`K` or the distinguished symbol `O`; `O` is **not** a field element.
For finite `x,y,z`, set

`S3(x,y,z)=(xy+xz+yz)^2+xyz+b`.

Define `Pair(x,y,u)` to mean `S3(x,y,u)=0` for finite `u`, or `x=y`
when `u=O`. Define `Root(u,v,t)` for finite target `t` as follows: use
`S3(u,v,t)=0` when both children are finite; require the finite child to
equal `t` when exactly one is `O`; reject when both are `O`.

**Pair lemma.** For rationally liftable `x,y` and any label `u`,
`Pair(x,y,u)` holds exactly when there exist `P,Q in E(K)` with those
abscissas and label `u` for `P+Q`. To see this, choose one lift of each
abscissa. If `x != y`, all sums `P+Q` and `P-Q` are finite. Substitution in
the group law shows their abscissas are the two roots, counted with
multiplicity, of `S3(x,y,z)`; its leading coefficient in `z` is
`(x+y)^2 != 0`. Negating either input exchanges these two orientations.
If `x=y != 0`, opposite lifts sum to `O`, while equal lifts double to the
finite abscissa `x^2+b/x^2`, the unique root of
`S3(x,x,z)=x^4+x^2*z+b`. If `x=y=0`, the unique point above zero has order
two, so only `O` is attained; `S3(0,0,z)=b != 0`. This covers all branches.

**Tree theorem.** Fix any four labels `x1,...,x4 in X` and finite target
abscissa `t`. The existential system

`Pair(x1,x2,u) and Pair(x3,x4,v) and Root(u,v,t)`

has a `K`-rational witness, including possible `O` children, exactly when
there exist rational point lifts `P1,...,P4` whose sum has abscissa `t`.
The pair lemma supplies point sums for each child. Each finite child can be
replaced by its negative while retaining the same pair abscissas: negate
both leaves of that pair. Thus the orientations required by `Root` can be
chosen consistently. The converse follows by taking the two actual partial
sums. The same proof holds for the crossed pairing `(1,3)|(2,4)`.

Consequently the **intersection** of the two complete O-aware tree systems
has the same projection onto `(x1,...,x4)` as either tree alone. It adds
constraints and auxiliary variables, but cannot add relations. The theorem
does not assert lower solving degree or cost. It also does not justify an
affine-only encoding that drops the `O` branches, or a polynomial solve that
allows nonrational leaf abscissas without membership constraints.

The finite-field generator and the independent Sage replay in this directory
check the theorem's implementation on the frozen small-field panel; those
checks are separate from the general argument above.
