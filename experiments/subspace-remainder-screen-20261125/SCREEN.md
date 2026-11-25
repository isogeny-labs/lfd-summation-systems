# Next summation-polynomial formulation screen

## Scope and prior-art gate

This is a mathematical screen, not a relation-cost measurement or a speedup
claim. The existing split S4 presentation already has a corrected
point-decomposition correspondence and a matched 33-case small-field solver
comparison in `experiments/semaev-matched-20261124/`. The one-T quotient has
no three-repeat total-cost crossover on its frozen prime-field range in
`experiments/torsion-crossover-20261125/`. The full rational E[2] invariant
on E: y²=x³+x is also not a new relation map: the duplication formula gives

    J(x) = (x²-1)²/[x(x²+1)] = 4 x([2]P).

Its poles are nonidentity two-torsion points. A doubled-point relation only
determines an original relation modulo E[2], so exact point recovery and group
verification remain necessary. A direct group baseline can use the same
doubling quotient. On the previous μ_N factor base, the extra two-torsion
translations generally leave the base, so quotient size alone is not a
same-base yield gain. The degree-two map classification and torsion-equivariant
map degree constraint are in [Faugère et al., Eurocrypt 2014, Props. 3 and 5](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf).

## Binary subspace-polynomial remainder

Let V be a k-dimensional F₂-subspace of F_{2^n}, and let
L_V(Z)=∏_{v∈V}(Z-v)=Σ_{j=0}^k c_j Z^{2^j}. It is linearized and squarefree.
For Q(Z)=Z²+e₁Z+e₂, reduce L_V modulo Q:

    L_V(Z) ≡ A(e₁,e₂) Z + B(e₁,e₂) (mod Q).

Starting with A₀=1, B₀=0, the residues of Z^{2^j} obey

    A_{j+1}=e₁ A_j²,       B_{j+1}=e₂ A_j²+B_j²,
    A=Σ c_j A_j,           B=Σ c_j B_j.

For distinct roots x,y of Q, Q divides L_V iff both x,y lie in V. If
e₁=0 then Q=(Z+x)² and cannot divide squarefree L_V, even when x∈V. The
following symmetric equations give the exact condition in both cases:

    e₁ A=0,
    e₂ A²+e₁ A B+B²=0.

Indeed L_V(x)=Ax+B and L_V(y)=Ay+B; the displayed equations are their sum
and product. Over a field, both vanish exactly when both images vanish. For
e₁≠0, they reduce to A=B=0. LRF already enforces e₁∈V, which implies
e₁A=L_V(e₁)=0, so the remaining test is the second equation. For e₁=0,
the second equation retains the diagonal x∈V. Signed points must still be
recovered and group-checked; the equations alone concern abscissae.

`verify.py` independently constructs L_V by adjoining basis vectors and
exhaustively checks all 1,024 monic quadratics over F₃₂ (k=2) and all
16,384 over F₁₂₈ (k=3), comparing the remainder equations with direct root
enumeration. `verification.json` records the passed counts.

## Cost gate and next experiment

This is the binary specialization of a known symmetric factor-base membership
idea already used in `experiments/torsion-symm-screen-20261125/`; it is not a
new relation family. The conditional Artin-Schreier image
W_{e₁}={x²+e₁x:x∈V} has dimension k−1 for nonzero e₁∈V, so one can
precompute a linear membership test for each e₁. That setup costs across
|V| values of e₁. Substituting LRF's affine e₁,e₂ kernel coordinates into
the symmetric equations can also produce a high-degree Boolean predicate.
The existing `experiments/ic-rootfind-yield/goodset_structure.py` reports
sparse, near-maximal-ANF-degree good-candidate sets at δ=6,9,12; those data
are a warning, not a lower bound against every solver.

Before a timed prototype, freeze a panel near both δ≈0 and a larger δ,
including random and planted targets on the identical binary curve and
factor-base points. Compare a remainder-circuit solver with native LRF's
2^δ candidate scan and native direct/MITM point methods. Charge subspace
polynomial construction, W_{e₁} preprocessing, every failed system, root and
sign recovery, and group verification; require equal canonical signed
relation sets and compare cost per independent rank. Reject this formulation
if solving the remainder equations costs at least as much as the scan on the
same systems. The measured LRF yield/cost ledger is in
`experiments/ic-rootfind-yield/README.md`; Gröbner-only or fall-degree
improvements are secondary diagnostics. [Kosters–Yeo](https://arxiv.org/abs/1503.08001)
give a relevant caution that low first-fall degree need not control the
observed regularity of a descended summation system.
