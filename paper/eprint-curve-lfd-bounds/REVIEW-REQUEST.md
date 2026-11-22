# Draft: request for an outside review

A ready-to-send draft for the author. The recipient is not chosen here: fill
in the bracketed fields, attach the PDF, and send it from the author's own
address. Nothing in this file has been sent. The status statements below
follow `../proof-status.json` and the manuscript as of 22 November 2026;
check them again if the paper changes before you send.

---

**To:** [Recipient name] <[recipient email]>
**Subject:** Review request: last fall degrees of elliptic curve summation systems (preprint, about 60 pages)

Dear [Recipient name],

I am writing to ask whether you would be willing to read and comment on a
preprint before I post it on the IACR ePrint Archive. I am asking you because
of [your work on ...].

The paper is "Last fall degrees of elliptic curve summation systems:
machine-checked bounds and exact values" (PDF attached, 63 pages; the
sources, Lean development and records are public at
https://github.com/isogeny-labs/lfd-summation-systems). It asks whether the
degree a Gröbner-basis or linear-algebra solver needs on Semaev
summation-polynomial systems stays bounded as the factor base grows, and
answers it for two-summand prime-field systems with two-power subgroup
factor bases: the last fall degree (LFD) grows with the depth k, though
slowly.

So that you know what rests on what, here is the status of each claim:

- **Checked by the Lean 4 kernel** (Mathlib, only `propext`,
  `Classical.choice` and `Quot.sound`; no `sorry` or project axioms):
  LFD ≤ k+2 for the two-summand chain presentation on every nonsingular
  cubic at every depth k ≥ 2 (Theorem 12); LFD ≥ 4 for unit-ideal
  two-summand systems; LFD ≥ 7 for three summands under a stated
  nondegeneracy condition; the first fall degree is 3 or 4 at every depth;
  depth restriction of excluded constants; the core of the m-summand carry
  bound for an arbitrary generator; and growth for a target of order two
  (Theorem 23).
- **Exact computation with certificates** (not Lean): LFD = 4 for k ≤ 5,
  5 for 6 ≤ k ≤ 9, 6 for 10 ≤ k ≤ 13 on a unit instance over F_65537, and
  LFD ≥ 7 at k = 14.
- **Written proofs only, not formalized:** the growth theorem for unit-ideal
  instances with t ≠ 0 and p > 2^(k+1) (Theorem 22, full proof in
  Appendix C; its Lean reduction keeps the curve model as a hypothesis); the
  degree lift behind LFD ≤ 7 at k = 14 (no final closure certificate); the
  steps of the m-summand bound outside its kernel core (Theorem 21); the
  equality of solving degree and LFD for unit ideals (Theorem 30); and the
  relation-search cost consequence in Section 8.
- **Empirical:** the law LFD(k) = ⌊(k+14)/4⌋, fitted on k ≤ 13.

No independent human review of the paper exists yet. The internal proof
reviews reported in Appendix B.4 were made by AI sessions, and the paper
says so; AI assistants were also used to prepare the Lean development, the
programs and the text, under my direction (see "Use of AI tools" at the end
of Section 1).

The parts where an expert reading would help most:

1. The written proof of the growth theorem (Appendix C), in particular its
   use of the vanishing-sum bound for function fields in characteristic 0
   and p > 2^(k+1).
2. The degree lift at depth 14 (Section 5.5): whether its written proof
   supports LFD ≤ 7 as stated.
3. The steps of Theorem 21 outside the kernel core: the grid step of part
   (1), the degree bound of S_(m+1) for m ≥ 4, and the genericity statement.
4. The Frobenius argument of Section 8 that an attempt on a target of order
   greater than two finds no relation exactly when its ideal is the unit
   ideal, and the cost consequence drawn from it.
5. Whether the Lean statements say what the paper says they do (Appendix A,
   Table 12, "Statements and Lean declarations").

Any part of this list, or a short opinion on whether the main claims hold,
would be valuable. I would of course acknowledge your help in the paper,
or not, as you prefer. If you cannot take this on, I would be grateful for
the name of someone who might.

[Optional: the date by which comments would be most useful.]

With thanks and best regards,

Adam Buran
adam@isogenylabs.io
