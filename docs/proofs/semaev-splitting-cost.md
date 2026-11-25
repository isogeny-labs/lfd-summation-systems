# Empirical cost reduction from splitting a Semaev system

**Status (2026-11-24):** literature benchmark and written algebraic lemma;
not kernel-checked. The claim below concerns reported average Gröbner-basis
costs for a specified benchmark. It does not assert a per-instance or
asymptotic complexity separation.

## Formulations

Let $k=\mathbb F_{2^n}$ and let
$E/k: y^2+xy=x^3+ax^2+b$ be nonsingular. Fix an
$\ell$-dimensional $\mathbb F_2$-subspace $V\subseteq k$ and a target
coordinate $r=x(R)$. Put

$$
S_3(x,y,z)=(xy+xz+yz)^2+xyz+b.
$$

The original formulation is

$$
S_4(x_1,x_2,x_3,r)=0,\qquad (x_1,x_2,x_3)\in V^3.
$$

The split formulation introduces $u\in k$ and uses

$$
S_3(x_1,x_2,u)=0,\qquad S_3(u,x_3,r)=0.
$$

Here $\ell$ is the source's $n^{\prime}$. Before preprocessing, coordinate
expansion over $\mathbb F_2$ gives respectively
$3\ell$ and $3\ell+n$ Boolean variables, and $n$ and $2n$ coordinate
equations, in addition to the Boolean field equations. In the split formulation,
the first block has degree at most three and the second at most two: squaring
in $k$ is $\mathbb F_2$-linear, and $r$ is fixed.
See [Karabina, Section 3.1](https://arxiv.org/pdf/1504.02347).

## Empirical proposition

Let $\widehat T_A(n,\ell,+)$ and $\widehat M_A(n,\ell,+)$ denote the
reported average Gröbner-basis time and memory for formulation $A$ in the
benchmark category with solutions, using the implementation and platform of
Huang's experiment. Table 5.1, printed page 54, reports the following for
$(n,m,\ell)=(23,3,6)$:

| Formulation | Reported degree of regularity | GB time (seconds) | Memory (MB) |
| --- | ---: | ---: | ---: |
| Original FPPR | 7 | 3888.70 | 6656.13 |
| Split | 4 | 168.12 | 794.30 |

Hence

$$
\widehat T_{\mathrm{split}}<\widehat T_{\mathrm{orig}},\qquad
\widehat M_{\mathrm{split}}<\widehat M_{\mathrm{orig}},\qquad
N_{\mathrm{split}}=41>18=N_{\mathrm{orig}}.
$$

The reported speedup is approximately $23.13$, and the ratio of original to
split memory is approximately $8.38$. These are published aggregate
measurements, not independently reproduced results or guarantees for every
instance. The benchmark uses Magma V2.18-9's `GroebnerBasis` function.
Source: [Huang, Chapter 5 and Table 5.1](https://catalog.lib.kyushu-u.ac.jp/opac_download_md/1500510/math0178.pdf#page=67).

## Algebraic correspondence and its scope

For a fixed tuple $(x_1,x_2,x_3)$, write
$A(t)=S_3(x_1,x_2,t)$ and $B(t)=S_3(t,x_3,r)$. By definition,

$$
S_4(x_1,x_2,x_3,r)=\operatorname{Res}_t(A,B).
$$

**Lemma.** Suppose $A$ and $B$ both have degree two and
$\deg\gcd(A,B)\leq 1$. Then the original equation holds if and only if
there is a unique $u\in k$ satisfying the split equations.

**Proof.** The resultant vanishes if and only if the gcd is nonconstant.
Under the hypotheses, this is equivalent to the gcd being linear over $k$,
which is equivalent to the existence of a unique common root in $k$.
Conversely, any common root makes the resultant zero. Coordinate expansion
in fixed bases is bijective, so this correspondence survives Weil descent
with Boolean field equations. $\square$

The degree-two hypotheses require $x_1\ne x_2$ and $x_3\ne r$,
since the leading coefficients are $(x_1+x_2)^2$ and $(x_3+r)^2$.
The lemma is a pointwise statement on this restricted locus, not an audit
of the benchmark's treatment of exceptional cases.

This lemma does not establish equivalence on the exceptional locus. Degree
drops and quadratic gcds require separate analysis; a common root over an
algebraic closure need not lie in $k$. Relating affine coordinates to the full
point-decomposition problem also requires handling intermediate points at
infinity and verifying recovered points and their sum.

## Supported conclusion

For the reported binary-field point-decomposition benchmark at
$(n,m,\ell)=(23,3,6)$, splitting the Semaev equation before Weil descent
strictly reduces the measured average Gröbner-basis time and memory,
despite increasing the number of Boolean variables from 18 to 41.

The numerical inequalities establish this empirical conclusion. Smaller
input degrees or first-fall-degree bounds alone do not establish it, and do
not prove an unconditional asymptotic separation. For the distinction between
first-fall degree and Gröbner-basis computation degree, see
[Kosters and Yeo, Notes on summation polynomials](https://arxiv.org/abs/1503.08001).

## Relation to the LFD paper

This comparison concerns splitting $S_4$ into two $S_3$ equations over a
binary extension field, followed by coordinate Weil descent. It supplies no
new bound on the last fall degree or solving degree of the recursive
presentations studied in the LFD paper. The reported degree of regularity
is retained as a source label, not identified with the paper's precisely
defined solving-degree invariant.
