# N8: an explicit rational two-torsion *unit-ideal* family

Status: **written proof candidate, not kernel formalized**. This is a distinct route to
unbounded LFD that avoids constructing the genus-one function field for every
target. A pinned Lean build is still required even for the accompanying algebraic
lemmas. Do not promote the result to a kernel theorem on these notes alone.

## The new family

For each k >= 2 let N = 2^k. Choose a prime p dividing 2^(2^k)+1.
The multiplicative order of 2 modulo p is 2^(k+1), hence
2N divides p-1, so p > 2N and N divides p-1.
Let B = F_p and K = F_(p^2). Fix a in K \\ B and put

    t = 1,  b = -1-a,  E_a: y^2 = x^3 + ax - 1-a.

Then (1,0) is a K-rational 2-torsion point, and

    4a^3 + 27b^2 = (a+3)^2 (4a+3) != 0.

Both exceptional values a=-3 and a=-3/4 belong to B. (p is odd.)

For every x,z in a field of characteristic p,

    S3(x,z,1; a,-1-a) = ((x-1)z -(x+2+a))^2.            (A)

The associated reduced curve is the rational graph

    z = (x+2+a)/(x-1).                                 (B)

This is **not** the genus-one curve E: it is the reduced S3 curve for
the target (1,0). It is sufficient for the ideal-kernel / no-wrap
method since S3 is killed by this evaluation.

## Unit ideal, without a GPU closure

The N-th roots of unity all lie in B because N | p-1. Every solution
of the chain equations over an algebraic closure has x,z in mu_N(B).
If S3=0, (A) would force

    a = (x-1)z -x-2  in B,

contradicting a not in B. Therefore the normalized two-summand
chain presentation has **no geometric common zero**. By the weak
Nullstellensatz its generated ideal is the unit ideal.

This argument works for EVERY a in K \\ B. It does not assume a
degree-D closure certificate or any LFD values.

The algebraic step is isolated in
`GrowthRationalUnit.s3_ne_zero_of_parameter_not_in_subfield`.
An explicit quotient-ring certificate is given below; the chain quotient bridge is not yet in Lean.

## Constructive Bézout certificate (stronger than Nullstellensatz)

The preceding unit-ideal argument can be made entirely algebraic,
without any algebraic closure or geometric zero criterion. In the
quotient by `X^N-1` and `Z^N-1`, write

    q = (X-1)Z-X-2,        u = q-a.

Because N divides p-1, the quotient has X^p=X and Z^p=Z, and hence
q^p=q by the characteristic-p Frobenius homomorphism. Since a is
outside F_p in F_(p²), c=a-a^p is nonzero. Thus

    u^p-u = (q^p-a^p)-(q-a) = c.

Consequently the explicit polynomial

    v = (u^(p-1)-1) / c

satisfies u*v=1 in this quotient. As S3=u², it follows that

    S3 * v² = 1  modulo (X^N-1, Z^N-1).

This is a **constructive unit-ideal certificate** with an explicit
high-degree inverse, not a closure degree bound. It is valid for any
field of characteristic p containing a with a^p != a, provided
N | p-1. The finite degree-two extension is simply one convenient
way of choosing such an a.

`GrowthUnitCertificate.power_difference_unit_certificate`,
`square_power_difference_unit_certificate` and
`square_unit_of_frobenius_mismatch` are Lean proof candidates for
the generic algebraic bridge. The remaining formal input is
`q^p=q` in the specific chain quotient (and Freshman's dream);
they must not be assumed in the final theorem. The numerical
regression also evaluates the exact inverse v at every tested pair.

---

## Concrete function field and two pole orders

Put F=K(X) and

    xi_1 = X,   xi_2 = (X+2+a)/(X-1) in F.

The first function has divisor (0)-(infinity). The second has
divisor (-(2+a))-(1). As a not in B, 0, infinity, 1 and -(2+a)
are four distinct K-rational places. In particular:

                  at infinity       at X=1
    ord(xi_1)        -1                0
    ord(xi_2)         0               -1.

This supplies the concrete `on_curve` and `ne_zero` data:
(A) shows S3(xi_1,xi_2,1)=0 and each xi_i is nonzero in F.
It also supplies precisely the two distinguished valuations.

## Characteristic-p Wronskian spread

Since K is perfect, F^p=K(X^p), and [F:F^p]=p. Let J be a
minimal dependent set over F^p of distinct monomials
xi_1^u xi_2^v, with m=|J|. Any m-1 of them are F^p-independent,
so m-1 <= p. This explicitly includes m=p+1; derivatives of
orders 0,...,p-1 suffice and **no p-th derivative is used**.

Choose a minimal relation sum_{e in J} phi_e xi^e=0 with
phi_e in (F^p)^*. Write v_e=phi_e xi^e. Their orders outside the
four marked places are multiples of p, and the ordinary derivative
d/dX annihilates F^p. The genus-zero Wronskian proof with r=m-1
columns gives, for the projective height,

    H(v_e : e in J) <= binom(m-1,2) * (-2 + 4)
                     = (m-1)(m-2).                    (C)

(The canonical divisor of dX is -2[infinity].) The elementary
projective-height inequality at the two distinguished places yields

    |ord_infinity(v_e/v_e')|+|ord_1(v_e/v_e')|
         <= 2H <= 2(m-1)(m-2).

The orders are congruent to -(u-u') and -(v-v') modulo p because
phi_e/phi_e' is a p-th power. Doubling and choosing the same
multiples of p gives

    ||2(u-u')||_p + ||2(v-v')||_p
         <= 4(m-1)(m-2)
         <= 6(m-1)(m-2),                              (D)

exactly `GrowthSpread.SpreadBound` (with a stronger constant).

The critical remaining Lean task is to construct F, F^p, their scalar
tower and the genus-zero valuation/Wronskian machinery, and prove (C)
and (D) without any analytic axiom. The currently committed
`GrowthRationalPoleArithmetic.spread_from_simple_poles` supplies
only the final integer step and **assumes** the height estimate.

## LFD consequence (written, not kernel)

Once (D) and the concrete evaluation field are instantiated in
`GrowthMain.CurveModel`, the existing Lean FEI and closure bridge
give 1 not in V_D(k) whenever k >= D+Lambda(D-1), and the unit-ideal
argument makes LFD(k)>D. Together with the existing universal k+2
upper bound, for k>=38,

    floor((k-10)/14)+1 <= LFD(k) <= k+2.

Thus this **constructed family** exhibits LFD(k)=Theta(k).
This does not prove the sharp measured slope 1/4, or the same lower
bound for every curve and target, nor say anything about ECDLP cost.

## Independent finite-field test

`python scripts/check_n8_rational_unit.py` uses a quadratic extension
K=B[alpha] with alpha^2 a nonsquare in B, sets a=alpha, and exactly
enumerates every (x,z) in mu_N(B)^2 for k=2,...,8. It checks 87,376
pairs, the curve discriminant and nonvanishing of S3. These tests
confirm finite instances but are not a replacement for (C) or Lean.

## Suggested formalization order

1. Build the pinned algebraic Lean candidates and repair real errors.
2. Formalize the mu_N-in-B root statement, no-common-zero -> unit ideal.
3. Define K(X), xi_1, xi_2 and prove all divisor/valuation statements.
4. Prove the characteristic-p P1 Wronskian estimate (C), retaining
   nonconstant F^p coefficients and the r=p boundary.
5. Instantiate SpreadBound, CurveModel, growth theorem, and LFD
   bounds. Verify `#print axioms` and the project audit on the exact
   commit before promoting the theorem to kernel.

Independent mathematical re-review should focus on the Wronskian
height convention, the r=p case, the rational-place valuations and
the exact normalization of the chain ideal.
