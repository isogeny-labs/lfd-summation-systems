# Quadratic lift of the binary root-membership condition

## Exact relation condition

Work over K=F_(2^n) on E: y²+xy=x³+a x²+b with b≠0. Let V be a k-dimensional F₂-subspace and T=(u,v) a point with u≠0. The existing LRF derivation in `csd/ic_rootfind.py` gives an affine solution space for the symmetric coefficients of a two-point decomposition. Write its free kernel coordinates as g∈F₂^δ. On every LRF solution, e₁(g), e₂(g), and the line slope λ(g) are affine F₂-maps into K; e₁(g) lies in V. Set

    x(z) = Σ_i z_i v_i,  z_i∈F₂,  and
    F(g,z) = x(z)² + e₁(g)x(z) + e₂(g).

The n coordinate equations F(g,z)=0 have Boolean degree at most **two** in (g,z): squaring is F₂-linear, the middle term is bilinear, and e₂ is affine. If they hold, x(z) and x(z)+e₁(g) are both in V and are exactly the two roots of Z²+e₁(g)Z+e₂(g). This includes e₁=0 and the repeated-root branch. With μ=λ(g)u+v+u, the points

    P₁=(x(z), λ(g)x(z)+μ),
    P₂=(x(z)+e₁(g), λ(g)(x(z)+e₁(g))+μ)

are the two curve intersections on the line through −T, with multiplicity. The cubic factorization in the LRF derivation therefore gives P₁+P₂=T, and each abscissa is in V. Conversely, every two-point relation summing to a nondegenerate T lies on such a line, yields an LRF solution g, and satisfies F(g,z)=0 for either root. Exact signed point recovery and a group-law check remain mandatory. Targets T=O or u=0 use the existing degenerate-case path; this screen makes no generic-lift claim for them.

This is a **quadratic presentation of the same relation family**, closely related to the existing LRF root recovery. It is neither a new summation law nor evidence that the *solving degree* is two. The extra k root-coordinate variables can make elimination expensive even though every input equation is quadratic.

## What the previous ANF screen proves

`experiments/ic-rootfind-yield/goodset_structure.py` measured the algebraic normal form (ANF) of the **projected success indicator** on g. Near-maximal ANF degree and an approximately singleton success set were observed at δ=6,9,12. Those measurements describe the projected Boolean function. They do not bound the degree or cost of an equation system with auxiliary root variables. A direct counterexample is the indicator of the unique solution g=0:

    1_{g=0} = ∏_(i=1)^δ (1+g_i).

Its ANF has degree δ, yet δ linear equations g_i=0 recover the point. Thus the earlier inference that a maximal-degree indicator forces kernel scanning was too strong. This correction does **not** imply that the quadratic lift is fast.

The normal-basis Hamming-ideal experiment in `research/hamming_ideal_pdp_20260930/RESULT.md` of the crypto repository is a relevant negative control: bounded-degree input presentations were sound, but the measured Gröbner oracle was far slower than direct enumeration on its frozen cases. Input degree, projected-indicator degree, observed solving degree, and full relation cost are distinct measurements.

## Cost gate before a solver prototype

The current m=3 kernel still scans one factor-base point for each residual two-point problem. The lifted system could at most improve its residual kernel search; it does not itself raise relation yield or remove that outer scan. On the frozen matched cells, a candidate solver must:

1. produce the same canonical signed point relations as native LRF, the image filter, and direct enumeration, including empty and diagonal cases;
2. record actual maximal solver degree and all setup, failed systems, root/sign recovery, verification, wall, CPU, and memory cost;
3. beat the fastest complete native relation-query arm on identical targets and factor bases, then include relation collection, linear algebra, and target resolution in a cold same-target comparison against the measured generic reference.

At (23,11), (29,13), and (29,14), the first-query BSGS medians are 10.2, 57.5, and 56.9 ms respectively; the existing image first-query medians are 22.6, 85.0, and 1,024.1 ms. A reduced input degree alone cannot pass this gate. The cheapest next falsifier is to build the lifted equations for the **already frozen** high-δ two-root systems, verify their solution projection against LRF, and measure one solver's complete query cost. Do not infer a large-field improvement from that screen.
