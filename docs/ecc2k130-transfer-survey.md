# Changing the group for ECC2K-130: what the program has already settled

Date: 2026-11-23. Context: PR #44 showed that nothing inside linearized root finding (factor base,
summand count, Frobenius orbits, candidate structure) moves the relation cost below about
`2^(n/3)` per relation, so any improvement must transfer the discrete logarithm to another group.
This note collects what `isogeny-labs/crypto` already proved or measured about such transfers for
`E: y^2 + xy = x^3 + 1` over `F_(2^131)`. Paths are relative to the crypto repository. The main
source is `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` (HYP). Spot-checked in this
session: HYP lines 81, 87, 200, 302 and 316 say what is quoted below, and the order of 2 mod 131 was
recomputed (it is 130). The rest is as reported by a read-only survey and has not been re-read line
by line. Results in the n = 7 study directories are tagged `[pcc:unverified]` there.

## Closed for n = 131

- **GHS/Hess descent to F_2.** 2 is a primitive root mod 131, so the magic number is 1 or 131 and
  the genus is 1, `2^130 - 1` or `2^130` (HYP §3). For b = 1 the m = 1 transfer is exactly the zero
  map (HYP §3.2). Isogeny walks do not change this (HYP §3).
- **Any correspondence to a curve over F_2** lands on genus at least 130 (HYP line 302): the
  restriction of scalars splits as `E x A` with `A` simple of dimension 130.
- **Covers over F_(2^131) itself** cost at least `2^131`, above rho.
- **Base change, Joux-Vitse, MOV**: closed (`RESEARCH_ECC2K130_EXTENSION.md`,
  `RESEARCH_ECC2K130_ROUTES.md`).
- **Cyclic covers**: genus at least 1300 over a base of genus at most 1, and no carrier of `A` below
  genus 394 has an automorphism of order 131 (`RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md`,
  `research/cover_decomposition_gap_search_20261018/`).

## Open

- **The window question** (HYP §6): is `A` (simple, dimension 130 over F_2) isogenous to the
  Jacobian of a curve of genus 130 to about 300, by an isogeny that can be evaluated? HYP models
  index calculus at genus 130 at `2^37.17` against rho at `2^60.81`. Every known construction lands
  outside the window; HYP says this is "not 'impossible'". Theory work, not a measurement.
- Another commutative algebraic group over F_2 with a point of order 131 (`LEADS.md`, L14):
  dimension at most 4 excluded, beyond that open.
- Whether the class group of `Z[tau]` or the Koblitz CM structure has been used for index calculus
  (`RESEARCH_ECC2K130_ROUTES.md`, route 5): no evidence found, which is not evidence of absence.

## Bearing on PR #44

The `2^(n/3)` relation-cost wall measured in PR #44 is a property of point decomposition on the
elliptic curve itself. It does not bound index calculus on a genus-130 carrier, which uses a
different relation (smooth divisors). So the two results are consistent: the only remaining route to
a sub-rho cost is the window question, and it is a construction problem, not a solver problem.
