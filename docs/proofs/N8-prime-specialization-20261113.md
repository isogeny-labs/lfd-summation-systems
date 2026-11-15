# N8: prime-field specialization of the unconditional complex E6 theorem

**Status (2026-11-13):** written proof target plus exact k=2 regression. This is **not a Lean kernel theorem**. It does not claim that all primes work or that characteristic-p rational Wronskian bounds are formalized. It builds on the verified fixed-complex E6 final theorem of #226.

## A concrete obstruction to uniform reduction

Fix `E:y²=x³+6x−7`, `t=1`. The two-torsion Semaev relation is the square of

```
h(x,z)=(x−1)z−x−8.
```

At depth `k=2`, `N=2^k=4`, and `p=17`, take `x=z=4`. Then `4^4=1 (mod 17)` and `h(4,4)=0` in `F_17`; in fact the displayed integer value vanishes exactly. The curve remains nonsingular because its discriminant factor `4*6^3+27*(-7)^2=2187=3^7` is nonzero modulo 17; moreover `2^(k+1)=8 | (17−1)`. Hence the root grid does contain a zero of the Semaev relation, so the corresponding chain system **cannot be a unit ideal** for that specific `(p,k)`. The grid point extends to a zero of the full chain system by taking successive squares for each chain variable. Consequently **there is no all-primes reduction** of the fixed E6 complex unit-ideal statement.

In the quotient `A=ℤ[x,z]/(x^4−1,z^4−1)`, multiplication by `h` has integral determinant

```
279791370354825 = 3^14 * 5^2 * 7^2 * 17 * 53^2.
```

Over `ℂ` this multiplication is invertible. Over fields with characteristic outside the displayed finite set it is invertible as well. Since `S3=h²`, invertibility of `h` and `S3` are equivalent. The exact independent regression is `scripts/check_n8_prime_specialization.py`.

## A different prime-field theorem from rational specialization

For every **fixed depth k≥24**, the complex theorem #226 proves that the chain ideal is the unit ideal and that the degree-D closure does **not** contain 1 for `D=floor((k−10)/14)`. The polynomial system is defined over `ℤ`. Both properties can be specialized to all but finitely many rational primes, *without* proving any positive-characteristic Wronskian theorem.

### Step A. Unit ideal after reduction

Because the chain inputs have integer coefficients, a complex polynomial Bezout representation `1=Σ_i H_i F_i` implies a representation with `H_i∈ℚ[x_1,…,x_(2k)]`. Reason: take a finite set of monomial coefficients supporting the complex witnesses and solve the resulting finite rational linear system over ℚ; consistency after tensoring to ℂ implies consistency over ℚ by rank invariance. Clear all denominators to obtain an integer `M_k≠0` and integer-coefficient polynomials `G_i` such that

```
M_k = Σ_i G_i F_i.
```

For every prime `p∤M_k`, reduction gives `1∈(F_i mod p)` in `𝔽_p[x_1,…,x_(2k)]`.

### Step B. Nonmembership in bounded closure after reduction

Let `D=floor((k−10)/14)` and `W_D` the finite-dimensional monomial span of total degree ≤D. The degree-D closure `V_D` is the least input-containing vector subspace of `W_D` stable under allowed multiplications. It may be computed by iterating finite-dimensional linear-span extension using variable multiplication, restricted to inputs whose degree is at most D. The iteration stabilizes in at most `dim W_D` strict enlargements. **Important technical requirement:** formally prove this iteration agrees with the repository's *inductive* `InClosure` definition, including degree cancellation; use existing `VariableClosure` and `closure_minimal_by_variables` as the starting point.

At fixed k,D all input coefficient matrices and multiplication maps are integral. Gaussian elimination over ℚ uses only finitely many nonzero rational minors. Exclude primes dividing the finitely many chosen pivot/minor numerators or denominators at every iteration. The ranks and membership test of the vector `1` then agree over ℚ and `𝔽_p`. The complex nonmembership of 1 implies rational nonmembership by faithful scalar extension, so for all primes outside another finite set `S_{k,D}`:

```
1 ∉ V_D(𝔽_p).
```

One cannot merely assert that failure to solve a linear system over ℚ automatically means failure modulo *every* prime; the excluded-minor condition is essential.

### Step C. Prime-field LFD bounds

For any `k≥24` and all primes outside the finite exceptional set from A+B, also exclude `p=2,3` to ensure characteristic ≠2 and nonsingularity. The existing selector upper theorem works in `𝔽_p`; the specialized unit ideal and closure nonmembership then give for the actual last fall degree `d_{p,k}`:

```
floor((k−10)/14) < d_{p,k} ≤ k+2.
```

Therefore **at each fixed depth k there are all but finitely many good primes**, including infinitely many admissible primes `p≡1 mod 2^(k+1)` by Dirichlet's theorem. Choosing one such good prime for each k gives a growing-prime family with `d_{p_k,k}=Θ(k)`. This is **not** a statement for one fixed prime with unbounded k (indeed `2^k|p−1` eventually fails), nor a uniform bound valid for *every* admissible prime, and it makes no ECDLP runtime claim.

## Required Lean formalization and tests

1. Formalize `InClosure` as a stabilizing finite-dimensional `ℚ`-linear iteration, with a characteristic-free base-change result at fixed k,D.
2. Prove a finite nonzero integer bad-prime certificate for each fixed k,D, obtained from the rational row space and rational Bezout certificate. Its statement must not assume the desired finite-field unit ideal or closure nonmembership.
3. Apply the existing #226 complex theorem over coefficients defined in ℤ, and transfer it to rational inputs by faithful scalar extension.
4. Construct a finite exceptional set `S_k` and prove the specialized prime-field LFD window for all primes outside it. Add an explicit existential prime corollary (Dirichlet/infinitude of primes in the progression if available in pinned Mathlib).
5. Run full pinned Lean 4.19.0/Mathlib v4.19.0 CI and `#print axioms` on every final theorem. Never label this document kernel verified.

**Separate stronger research goal:** a uniform positive-characteristic `GrowthSpread.SpreadBound` over `F_p(X^p)`, or an all-eligible-primes theorem with a different unit family, still requires characteristic-p geometry/linear algebra and distinct unit-ideal certificates. NIST P-192/P-224 applications are not consequences of the specialization proof.
