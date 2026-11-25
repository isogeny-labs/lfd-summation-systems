# Complete PDP correspondence after membership and an infinity branch

**Status:** written proof, not kernel-checked. The exhaustive audit is a
separate finite check; it is not the proof of the general statement below.

Let k be a finite field of characteristic two, E/k the nonsingular curve
y²+xy=x³+ax²+b, R an affine k-point and V any subset of k. Let X be the set
of x-coordinates of affine E(k)-points. Write r=x(R). All tuple sets below
are ordered x-coordinate triples; they do not count signed point triples.

Define

- D = {(x1,x2,x3) in (V intersect X)^3 : S4(x1,x2,x3,r)=0};
- F = {(x1,x2,x3) in (V intersect X)^3 : there is u in k with
  S3(x1,x2,u)=S3(u,x3,r)=0};
- H = {(x,x,r) : x in V intersect X}, when r is in V, and the empty set otherwise;
- P = {(x(P1),x(P2),x(P3)) : Pi affine in E(k), x(Pi) in V,
  P1+P2+P3=R}.

**Proposition.** D = F union H = P. No degree or gcd restriction is needed
for this PDP statement. Membership in X and the infinity branch are essential
parts of the statement; it asserts no equality of the raw polynomial zero sets.

**Proof.** Use the defining property of Semaev's summation polynomials over
an algebraic closure: an S_j zero corresponds to affine points with those
x-coordinates whose sum is O. Since every xi is in X and r=x(R), all points
with these coordinates are k-rational: each corresponding quadratic in y
already splits over k. The S4 property therefore gives D=P, allowing a
simultaneous sign change if the target point is -R.

For a decomposition P1+P2+P3=R, put Q=P1+P2. If Q is not O, its coordinate
u is in k and the two S3 equations hold, so the tuple is in F. If Q=O,
then x1=x2 and P3=R, so it is in H. Thus P is contained in F union H.

Conversely, an S3 zero for x1,x2,u supplies points A1,A2,U with sum O.
A1 and A2 are rational by membership, so U is rational too. The second
S3 equation supplies points U',A3,T with sum O, where x(U')=x(U),
x(A3)=x3 and T is R or -R. U' is U or -U. Simultaneously change the
signs of A1,A2,U if needed to make U'=-U. Adding the two identities
then gives A1+A2+A3+T=O, and a simultaneous sign change if necessary
gives a decomposition of R. Every tuple in F is consequently in P.
Every tuple (x,x,r) in H is realized by A,-A,R for any point A over x,
including the point of order two. Hence H is contained in P as well. ∎

The standard summation-polynomial property is stated in
[Kosters–Yeo, Proposition 2.1](https://arxiv.org/pdf/1503.08001).
The split construction is described in
[Karabina, Section 3.1](https://arxiv.org/pdf/1504.02347).

## Small explicit counterexample to the uncorrected affine split

Over GF(4)=GF(2)[w]/(w²+w+1), take E:y²+xy=x³+1 and R=(0,1).
The triple (0,0,0) is a true decomposition: R has order two and R+R+R=R.
S4(0,0,0,0)=0, but each S3(0,0,u)=1, so there is no finite auxiliary u.
The intermediate sum R+R is O. This witness is checked against both the
Sylvester resultant and the enumerated curve group in the saved audit/tests.

A different phenomenon appears for b=w on y²+xy=x³+w: target (w,1),
triple (1,w,1). Both intermediate quadratics are w*u²+w*u+1,
with no root in GF(4), but their resultant is zero. The triple is not a
rational-point decomposition. Membership constraints exclude it.
