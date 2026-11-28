# Lead 3: an explicit A' with zeta not in End(A'), as a Deligne module

The lead-3 agent did this analysis; the harness blocked it from saving the report file. The parent session
re-checked the numbers below in gp: lambda^2 = 2 mod p and ord(lambda) for all four (ell, p), and the lattice
counts. It also checked the result against lead 1 (see the last section).

Reproduce, from this folder:
- `python3 deligne_chars.py`
- `python3 group_structure.py`
- `gp -q stable_lattices.gp` (about 2 min)
- `gp -q group_end_gor.gp` (about 7 min)

Two independent methods were used:
- **Method A (PARI):** lattice linear algebra, Hermite and Smith normal forms.
- **Method B (pure Python):** characters over F_p and F_{p^2}, and local p-adic group structure.

## 1. Minimal degree (proved, and attained)
**Lower bound.** An F_2-isogeny A -> A' with zeta not in End(A') has degree divisible by a prime dividing U_ell.
Away from U_ell, R = Z[pi, pi-bar] is locally O_K, so a lattice whose index avoids U_ell is O_K-stable.

**Minimum attained.** The minimum is 23 at ell = 11, 457 at ell = 19 and 263 at ell = 131. Among isogenies of
p2-power degree at ell = 131, the minimum is p2 itself. At ell = 5 and 13 no such A' exists.

**Count of R-stable lattices of index p:**

| ell, p | index-p R-stable lattices | O_K-stable | not O_K-stable |
|---|---|---|---|
| 11, 23 | 218 | 20 | 198 |
| 19, 457 | 7788 | 36 | 7752 |
| 131, 263 | 34058 | 260 | 33798 = 129 * 262 |
| 131, p2 | p2 + 1 | 0 | p2 + 1 |

**Controls.** Primes not dividing U_ell give 0 non-O_K-stable lattices: ell = 5 at 11 and 31, ell = 13 at 53
and 79, ell = 11 at 29 and 67, ell = 19 at 191, ell = 131 at 787.

## 2. The constructed A'
Use the conjugation-stable merged character pair (chi2 = chi1 o conj).

- **End(A') = O, an order of index p in O_K.** For split p (263), O = {x : chi1(x) = chi2(x) mod p}.
  For p2, O = {x : chi1(x) mod P in F_p}.
- **Aut(A') = {+-1}.** There is no automorphism of order ell.
- **O is Gorenstein**, by the conductor test and, at ell = 11 and 19, by invertibility of the trace dual.
- **Frobenius on the kernel.** The kernel is cyclic of order p, etale, and Frobenius acts on it by lambda,
  with lambda^2 = 2.

| ell, p | lambda | kernel points defined over F_{2^m}, m = ord(lambda) |
|---|---|---|
| 11, 23 | 18 | 11 |
| 19, 457 | 377 | 152 |
| 131, 263 | 110 | 262 |
| 131, p2 | 14975904886999252 | (p2 - 1)/12 = 12208813656794060 |

## 3. Group structures
- **Over F_2.** A'(F_2) = A(F_2) = Z/N, so the ECDLP is unchanged.
- **Over extensions.** The structures differ from A's only at p-power torsion over particular extensions,
  e.g. the 263-part at k = 0 mod 262. The two methods agree on every case computed.

## 4. Isogeny invariance
**Invariant.** These are shared by every member of the isogeny class, so every test already run on A
holds for A':
- the Weil polynomial P_A;
- virtual point counts and place counts;
- simplicity, ordinarity, and End^0 = K.

**Not invariant, and testable on the Deligne module:**
- End and Aut.
- **Rosati stability.** A principal polarization forces End to be conjugation-stable. So only lattices from
  the conjugation-stable pair can be principally polarizable. At ell = 131, p = 263 that leaves 262 of 33798;
  the other 33536 are not principally polarizable, hence not Jacobians.

## Combined with lead 1
Lead 1's parity rule: eps = (-7/131) = -1, and n(T) counts the self-glued inert places where End is
non-maximal.
- **Index-263 models.** 263 is of split type, so n = 0. The rule then forbids a principal polarization, so
  even the conjugation-stable 263-model is not principally polarizable.
- **Every ppav without zeta** is non-maximal at v0 above p2.
- **Smallest isogeny to a ppav model.** So the smallest isogeny from A to a possible Jacobian model has
  p2-power degree. Its kernel points live over F_{2^m} with m = (p2-1)/12, about 1.2e16.

**Practical consequence (an observation, not a theorem about every route).** If such a model were a
Jacobian, carrying a DLP instance across the minimal-degree isogeny would mean evaluating an isogeny whose
kernel generators lie in a field of degree about 1.2e16 over F_2. That rules out kernel-point (Velu-type)
methods. Other ways to evaluate an isogeny exist in principle, but none is known in dimension 130.

## Scope
- A' is constructed as an explicit lattice. No equations or curve were produced.
- Lattices were counted, not isomorphism classes: Pic(O) was not computed.
- Primality of p2 rests on PARI.
