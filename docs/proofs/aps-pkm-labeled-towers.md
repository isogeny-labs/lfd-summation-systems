# Labeled factor bases and rational-map towers

This note gives a probability bound for a random APS-style factor base, an
explicit composite-order construction, and its cost and prime-order limit.
It does not claim a faster prime-order ECDLP algorithm.

## A random APS base and a short tower

Let $G=\langle P\rangle\subseteq E(\mathbb F_p)$ have odd prime order $n$,
let $Q=wP$, and set

$$
X_G=\{x(T):T\in G\setminus\{O\}\},\qquad |X_G|=(n-1)/2.
$$

Independently sample $u_j,v_j\in\mathbb Z/n\mathbb Z$ uniformly and put
$R_j=u_jP+v_jQ$. For each fixed $w$, the scalar $u_j+v_jw$ is uniform because
$u_j$ is uniform. Conditioned on all $R_j\ne O$ and all their $x$-coordinates
being distinct, the unordered coordinate set $S$ is uniform among the
$\binom{(n-1)/2}{s}$ subsets of $X_G$ of size $s$.

**Proposition.** Fix $D\ge2$ and $k\ge1$. Let $\mathcal T_{D,k}$ be the
rational functions over $\mathbb F_p$ expressible as a composition of $k$
rational functions of degree at most $D$. Then

$$
\Pr\left[\exists L\in\mathcal T_{D,k}:
 S=\{x\in X_G:L(x)=0\}\right]
\le
\min\left(1,\frac{p^{2k(D+1)}}{\binom{(n-1)/2}{s}}\right).
$$

*Proof.* A rational function of degree at most $D$ has a presentation $A/B$
with both polynomials of degree at most $D$. Counting all coefficient pairs,
including invalid and duplicate presentations, gives at most $p^{2(D+1)}$
choices. There are at most $p^{2k(D+1)}$ tuples of $k$ stages. Each tuple
determines at most one set $\{x\in X_G:L(x)=0\}$. The random set $S$ is
uniform among $\binom{(n-1)/2}{s}$ possibilities. Apply the union bound. □

If $n\asymp p$, $s\asymp p^{1/m}$ with fixed $m\ge2$, $D$ is fixed, and
$k=O(\log s)$, then

$$
\log_p\binom{(n-1)/2}{s}=(1-1/m+o(1))s.
$$

The numerator of the probability bound has exponent $O(\log s)$, so the
probability is at most $p^{-(1-1/m+o(1))s}$. This bound uses the specified
random-sampling model and stages defined over $\mathbb F_p$. It does not
exclude a deliberately constructed base.

## A labeled tower from a smooth kernel

Suppose $E(\mathbb F_p)=\langle P\rangle$ has order $N$. Choose a divisor
$d\mid N$ with $\gcd(d,p)=1$, put $M=N/d$ and $K=MP$, and let
$H=\langle K\rangle$, of order $d$. Choose $t$ points
$R_i=u_iP+v_iQ$ whose images under the quotient isogeny
$\phi:E\to E/H$ are distinct up to sign and are not 2-torsion. Define

$$
F=\{R_i+jK:1\le i\le t,\quad 0\le j<d\}.
$$

Every point has an affine expression in the unknown log:

$$
R_i+jK=(u_i+jM)P+v_iQ.
$$

Write $x(\phi(T))=A(x(T))/B(x(T))$ and $c_i=x(\phi(R_i))$. Define

$$
L(X)=\prod_{i=1}^{t}\bigl(A(X)-c_iB(X)\bigr).
$$

For a normalized separable isogeny of degree $d$, $\deg A=d$ and
$\deg B=d-1$, so $\deg L=td$. Its $\mathbb F_p$-roots are exactly the
distinct $x$-coordinates of $F$:

1. If $T\in E(\mathbb F_p)$, then $L(x(T))=0$ exactly when
   $\phi(T)=\pm\phi(R_i)$ for some $i$, equivalently
   $T\in\pm R_i+H$. Negation does not change $x$.
2. A root $x\in\mathbb F_p$ that does not lift to $E(\mathbb F_p)$ lifts
   to a point $T\in E(\mathbb F_{p^2})$ with Frobenius $T\mapsto-T$.
   Since $\phi$ is defined over $\mathbb F_p$, its image would also
   satisfy Frobenius $\phi(T)\mapsto-\phi(T)$. But
   $x(\phi(T))=c_i$ makes $\phi(T)=\pm\phi(R_i)$ rational. This would force
   $\phi(R_i)$ to be 2-torsion, excluded by assumption.

If the prime factors of $d$ are small, a subgroup chain factors $\phi$
into small-degree isogenies. The outer polynomial has degree $t$, so the
description is a bounded-degree tower when $t$ is also bounded.

For a relation $\sum_a\epsilon_a(R_{i_a}+j_aK)=O$, put

$$
A_0=\sum_a\epsilon_a(u_{i_a}+j_aM),\qquad
B_0=\sum_a\epsilon_a v_{i_a}.
$$

Then $A_0P+B_0Q=O$. If $\gcd(B_0,N)=1$, one relation gives
$w=-A_0B_0^{-1}\bmod N$.

### Checked instance

On $E:y^2=x^3+x+2$ over $\mathbb F_{101}$, the group has order $100$.
Take $P=(6,27)$, $Q=(9,29)=17P$ for verification, $d=4$, and
$K=25P=(1,99)$. Choose:

| $i$ | $(u_i,v_i)$ | $R_i$ | $x(\phi(R_i))$ |
|---:|---:|---:|---:|
| 1 | $(30,38)$ | $(18,14)$ | 100 |
| 2 | $(13,92)$ | $(41,53)$ | 73 |
| 3 | $(50,61)$ | $(85,43)$ | 0 |
| 4 | $(19,11)$ | $(20,12)$ | 7 |

The two degree-two $x$-maps in the isogeny chain are

$$
h_1(X)=\frac{X^2+X+4}{X+1},\qquad
h_2(Z)=\frac{Z^2-3Z+8}{Z-3}.
$$

Their composition is $A/B$, where

$$
A=X^4-X^3+11X^2+9X+12,\qquad B=X^3-X^2-X+1.
$$

The polynomial $L=(A-100B)(A-73B)A(A-7B)$ has degree $16$ and exactly
16 roots in $\mathbb F_{101}$, matching the factor-base coordinates.
The factor-base relation

$$
(20,12)+(26,38)=(85,43)
$$

is $R_4+(R_4+3K)=R_3$. Its labels give $63P+61Q=O\pmod{100}$, hence
$w=-63\cdot61^{-1}=17\pmod{100}$. The [verification script](../../scripts/verify_aps_pkm_labeled_towers.py)
enumerates the curve and checks both rational maps, all roots, the relation
and the recovered log. Run it with
`python3 scripts/verify_aps_pkm_labeled_towers.py`.

### Orbit recurrence is a weaker property

There is an easy way to generate labeled points in a prime-order group:
take $\{[2^j]P,[2^j]Q:0\le j<4\}$. Their $x$-coordinates follow the
degree-four doubling map, but the polynomial having those coordinates as
roots need not decompose.

For $E:y^2=x^3+x+32$ over $\mathbb F_{101}$, the group has prime order
101. With $P=(4,10)$ and $Q=3P$, the eight coordinates are
$4,23,49,27,43,34,24,28$. Their root polynomial is

$$
f=x^8+71x^7+66x^6+38x^5+25x^4+97x^3+29x^2+34x+3
\quad\text{in }\mathbb F_{101}[x].
$$

It has no nontrivial polynomial decomposition. A degree-two right
component can be normalized to $h=x^2+cx$; the $x^7$ coefficient
forces $c=71/4=43$. Then $f(-43-x)=f(x)$ is necessary, but
$f(-43)=93\ne3=f(0)$. A degree-four right component can be normalized
to $h=x^4+86x^3+72x^2+89x$. The top coefficients force this choice.
If $f=g(h)$ with $\deg g=2$, the residual $f-h^2$ must equal
$ch+c_0$. Its $x^4$ coefficient is 36, forcing $c=36$; its $x^3$
coefficient would then be $36\cdot86=66$, but is 7. Since the only
proper degree factorizations of 8 are $2\cdot4$ and $4\cdot2$, $f$
is indecomposable. This is a counterexample to the inference from a
short point recurrence to a compositional root polynomial, not a
general impossibility theorem.

## Cost and prime-order limit

The quotient $E(\mathbb F_p)/H$ has order $M=N/d$. For independent random
representative images, a nondegenerate $m$-tuple has probability about $1/M$
of summing to zero. Thus $t^m$ candidate tuples give about $t^m/M$
relations in expectation. A hash search over $m-1$ choices costs about
$t^{m-1}$ quotient-group operations. At $t\asymp M^{1/m}$, this is
$M^{1-1/m}\ge\sqrt M$ for $m\ge2$, before setup, factor-base enumeration
and the smooth-kernel discrete logarithm. This is a heuristic cost estimate
for this search, not a bound on every algebraic solver.

If the target subgroup $G=\langle P\rangle$ has prime order $n$ and
$n\nmid |H|$, then $H\cap G=\{O\}$, so the isogeny is injective on $G$.
Every fiber contains at most one labeled point of $G$. If $n\mid|H|$,
the kernel degree is at least $n$, and its subgroup chain must include a
stage of degree at least $n$. This construction does not enlarge a factor
base inside a large prime-order subgroup.

## Sources and scope

- Petit, Kosters and Messeng,
  [Algebraic Approaches for the Elliptic Curve Discrete Logarithm Problem
  over Prime Fields](https://people.maths.ox.ac.uk/petit/files/16PKC_primeECDLP.pdf):
  rational-map composition for a factor base.
- Amadori, Pintore and Sala,
  [On the Discrete Logarithm Problem for Prime-Field Elliptic Curves](https://ora.ox.ac.uk/objects/uuid%3Ac99fdb39-f5b1-428d-999c-770b63363e25):
  randomly generated points with affine expressions in the unknown log and
  one-relation extraction.
- Vélu isogeny formulas:
  [Robert, Isogenies and endomorphism rings of elliptic curves](https://www.math.u-bordeaux.fr/~damienrobert/pro/teaching/slides/2011-09-Nancy-ECCSummerSchool.pdf).

The probability bound is an elementary parameter count, with no novelty
claim. The worked construction is a composite-order example, with no
prime-order ECDLP breakthrough claim.
