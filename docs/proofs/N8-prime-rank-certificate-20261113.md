# N8 prime-field LFD: exact rank-specialization certificate

**Evidence label:** mathematical linear-algebra lemma + exact executable regression. This file is **not** a Lean formalization of the N8 prime-field theorem.

## Basic exact specialization lemma

Let `M∈ℤ^(n×m)`, `b∈ℤ^n`, and `r=rank_Q(M)`.

**Outside case.** If `b∉colspan_Q M`, then `rank_Q[M|b]=r+1`. Choose one nonzero `r×r` minor `Δ` of `M` (take `Δ=1` if `r=0`) and one nonzero `(r+1)×(r+1)` minor `Γ` of the augmented matrix. Then for any prime `p∤ΔΓ`,

```
rank_Fp(M mod p) = r,
rank_Fp([M|b] mod p) = r+1,
b mod p ∉ colspan_Fp(M mod p).
```

The upper rank of `M mod p` is at most `r` because every `(r+1)×(r+1)` minor vanishes **as an integer**, and thus after every reduction. The nonzero `Δ mod p` proves rank at least `r`. The nonzero `Γ mod p` proves augmented rank at least `r+1`. A single chosen pair of integer minors therefore yields an explicit finite bad-prime set.

**Inside case.** If `b∈colspan_Q M`, rational Gaussian elimination produces `x∈ℚ^m` with `Mx=b`. Clearing denominators gives `z∈ℤ^m` and `d∈ℤ_{>0}` with `Mz=db`. For every prime `p∤d`, the residue of `d` is invertible and `b mod p=d^{-1}M z mod p`, so membership persists.

Together these two cases guarantee that the membership decision of **any fixed integer matrix** and right-hand side is constant under reduction outside a finite set, with completely explicit certificates.

## Application to the N8 closure (not yet formalized)

For a fixed depth k and degree D, choose a monomial basis of `W_D`. Each stage of the linear closure algorithm expands a bounded-degree space by multiplying its suitable lower-degree part by variables and adjoining the input generators. All linear operations are defined by integer matrices. The crucial next theorem must establish that this finite algorithm exactly computes the repository's inductive `InClosure`, **including cancellation where the support degree of a linear combination falls**.

A finite sequence of stages and their rank certificates determines a finite exceptional set. Outside it, all pivot and target-membership decisions match their rational counterparts. Because the complex fixed-E6 theorem provides bounded-degree nonmembership and unit ideal, a rational-base-change argument can then transfer those two properties and obtain an unconditional prime-field LFD window for each fixed k≥24, all but finitely many primes.

**Do not bypass these missing interfaces:** (1) full equivalence of finite closure iteration with `InClosure`, (2) faithful extension from ℚ to ℂ for unit/closure statements, (3) assembly of a common finite set of primes preserving all iteration stages, and (4) kernel checking of the resulting prime-field LFD theorem. The determinant certificate alone does not imply any N8 growth theorem.

## Regression

`python scripts/check_n8_prime_rank_certificates.py` exercises exact `Inside` and `Outside` certificates, checks exact determinants and integer identities, and verifies the classification against finite-field matrix ranks at multiple primes outside each computed divisor. The script is a diagnostic/algorithmic helper, not a substitute for Lean.
