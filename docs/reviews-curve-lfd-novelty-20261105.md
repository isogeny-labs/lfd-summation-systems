# Prior-art and novelty check: `paper/curve-lfd-bounds.tex`

Date: 2026-11-05. Scope: candidate contributions (a), (b), (c). Read first: abstract, Sections 1-4, 7 (related work) of the manuscript and `docs/IC_NOVELTY_AUDIT.md`. Items already settled in that audit (HKYY, Caminata-Gorla, Kousidis-Wiemers, Guo et al. Lean GB, CompPoly, Huang semilocal, Preparata) are not re-argued here except where this check adds something.

## Headline findings

1. **Contribution (b) has a direct empirical precedent that the manuscript does not cite as such.** Petit-Kosters-Messeng (PKM, PKC 2016), Table 4 ("Prime case, p-1 subgroups", n1 = 1..11), reports the Magma corrected degree of regularity D_av^corr of the two-summand subgroup-chain systems as **4 for n1 <= 5, 5 for n1 = 6..9, 6 for n1 = 10, 11**. The step positions (5->6 and 9->10) are exactly the manuscript's LFD breakpoints (4 for k<=5, 5 for k=6..9, >=6 at k=10). PKM report it as an observation with no explanation. The manuscript's new content is (i) exact certified LFD rather than a solver step degree, and (ii) the dimension-count account. The "who measured the steps first" credit belongs to PKM and must be stated.
2. **Contribution (a) is, as a mathematical statement, the Groebner-basis property of a system with pairwise coprime pure-power leading terms** (Buchberger's first criterion), applied to the PKM chain construction. The manuscript's own Remark concedes this. The grid-of-roots-of-unity proof is a different route to a textbook fact; the bound `max(mk+1, delta)` is the generic finite-quotient/Macaulay-type bound. Novelty is limited to the Lean formalization and the explicit statement. I found no prior statement of the exact `LFD <= max(mk+1, delta)` form, but a failed search is not proof.
3. **Contribution (c) is partly anticipated qualitatively** (Kosters-Yeo Section 5: a structured vector space changes observed dreg; GGMP SAC 2020: Frobenius-invariant bases reduce the number of systems and linear algebra, not claimed to change degree). No source I could access states the specific "no change in exact LFD beyond solution count" finding. It is a small-grid empirical statement, not a theorem.
4. **Two citation/characterization errors in the manuscript's related-work section** (Section "Related work"): (i) Kosters-Yeo is described as giving "LFD bounds for Weil descent of nonprime fields". It is not: it is about the *first* fall degree (= 2 for ordinary binary S3, via the trace morphism), NP-completeness of summation-polynomial evaluation, and growth of dreg with n. The LFD-for-Weil-descent results are HKYY, Huang-Kosters-Yeo CRYPTO 2015, Huang-Kosters, Gorla-Mueller-Petit, Huang 2103.07282. (ii) PKM is said to "bound its first fall degree heuristically". In the PKM text I extracted, the first fall degree assumption is mentioned only for the binary case (Section 2.4), and the prime-case evidence is Magma degree-of-regularity experiments (Section 4) compared with Bardet's semi-generic value (3+n)/2. I found no FFD bound for the prime systems in PKM. Re-check against the published version before keeping the sentence.

## Method and access limits

Searched (web search tool, standard and extended): PKM prime-field subgroup factor base + FFD; HKY CRYPTO 2015; LFD bounds with field equations; Groebner basis formalization Coq/Isabelle/Lean; Amadori-Pintore-Sala; GGMP subfield curves; "Kosters Pagano" FFD Weil descent; Semaev/Karabina binary symmetric subspaces; "formalization last fall degree Lean/Coq/Isabelle"; follow-ups on PKM degree growth.

Full text read (downloaded, converted with pdftotext): PKM PKC 2016 author PDF; Kosters-Yeo arXiv:1503.08001 (Sections 4-5); Caminata-Gorla arXiv:2112.05579 (Thm 2.8, 2.9, 5.1, 5.2, Examples 4.2-4.4); HKYY arXiv:1505.02532 (intro, Thm 1.1, Section 6); Gorla-Mueller-Petit arXiv:2011.01050 (Thm 11, Cor 1); Kousidis-Wiemers arXiv:1906.05594 (Remark 3.3, Thm 3.2; grep only); Steiner arXiv:2310.03637 (grep of structure and the coprime-leading-monomial argument; not fully read).

Abstract/metadata only: GGMP SAC 2020 (ePrint 2020/1315 returned HTTP 403; Birmingham record read, no Groebner content visible); Amadori-Pintore-Sala; Semaev ePrint 2015/310; Karabina ePrint 2015/319; Huang arXiv:2201.00653 (product decomposition: not LFD-relevant).

Could not access: all `eprint.iacr.org` PDFs (403), Springer/Elsevier full texts, Google Scholar, DBLP query results, zbMATH/MathSciNet, Mathlib/AFP/Coq-contrib source trees, and the Huang-Kosters (n-independent bound) paper by title; I could not identify any "Kosters-Pagano" paper on FFD bounds (searched; nothing found; the name may be a misremembering). Galbraith-Gebregiyorgis INDOCRYPT 2014 was not opened; the manuscript's description of it ("measured Groebner behaviour for subfield and random coordinate spaces") is unverified by me. Rocq/Isabelle/Lean claims below rest on search-result descriptions, not on reading the libraries.

## (a) Chain equations as a reduced-ideal-trivial subsystem; LFD <= max(mk+1, delta); Lean

**Construction is PKM's.** PKM Section 3.3: for V the subgroup of order prod p_i in F_p^*, they set L_j(x) = x^{p_j} for j < n0 and L_{n0}(x) = 1 - x^{p_{n0}}, and L = L_{n0} o ... o L_1; the polynomial system is the "SRA-like" decomposition of L(x_i) = 0 into chain variables (their experiments use p_j = 2). That is exactly the squaring chain u_{j}^2 = u_{j+1}, u_{k-1}^2 = 1. Coset variant mentioned there as well ("we could also choose V as a coset aG"). The manuscript's "prime presentation due to PKM" attribution is right; cite PKM Section 3.3 and Section 4.1 (set-up) and Table 4.

**The mathematical content is the Groebner-basis property.**
- Leading terms u_v^2 are pairwise coprime for any degree-compatible order, so G_{m,k} is a Groebner basis (Buchberger's first criterion; Cox-Little-O'Shea, *Ideals, Varieties, and Algorithms*, Ch. 2 Section 9, Prop. 4). The manuscript's Remark states this and says it holds in every characteristic.
- "Every polynomial in <G> all of whose monomials are standard is zero" is the standard characterization of a Groebner basis (CLO Ch. 2 Section 6, Thm 6). The manuscript's Theorem 2 therefore restates this fact, proved by a different argument (grid vanishing over the 2^k-th roots of unity, char != 2 only).
- Same coprime-leading-monomial structure is used for iterated systems by Steiner (ToSC 2024(1)), who cites CLO15 Ch. 2 Section 9 for chains of the form x_{i+1} - f(x_i) and field equations y^q - y and builds LFD/solving-degree bounds on it (Section 3, Theorem 3.2/3.7; Section 8 constructs systems with given degree falls). This is the closest published method to (a).
- Caminata-Gorla Thm 2.9 (with Remark 2.10): for a Groebner basis {g_i} w.r.t. a degree-compatible order, d_F = max{d_{g_i} : deg g_i < d_{g_i}}, and d_F = 0 when all d_{g_i} = deg g_i. Their Example 4.3 (F = {x1^2 - 1, x2^2 - 1}, which is the m = 2, k = 1 chain) has LFD 0. So "the chain subsystem alone has LFD 0, and any enlargement's LFD is controlled by the new Groebner elements" is available from the general theory, and the paper's Theorem 3 is an instance of "LFD <= (max degree bound on GB elements of the enlarged ideal) in a finite quotient with dim bound".
- Bounds of the type "LFD / solving degree is bounded when the system contains field equations or other pure powers": CG Thm 5.1 ([Dubois-Gama; Lazard-type] sd(F) <= reg(F^h)); CG Thm 5.2 (Semaev-Tenti) sd <= 2 d_reg - 2; HKYY Thm 1.1 / Thm 4.4 (n-independent LFD bound of Weil descent, depending on number of solutions s); Gorla-Mueller-Petit Thm 11 / Cor 1 (LFD <= (q-1)ceil(log_q(d)+1)+1 for univariate systems with field equation). None states `max(s+1, delta)` for a general GB-with-standard-monomials-of-degree-s. I did not find the exact statement, but the finite-quotient argument is routine and I would expect it to be folklore. Do not describe it as a new bound; the manuscript already says "not new as a bound".

**Meaningfulness caveat.** For the normalized two-summand system the proved bound 2k+1 is weak against exact values 4/5 (Table 2 of the manuscript). Reviewers will see (a) as a correct but loose dimension bound. HKYY Prop 2.8 says that solving a zero-dimensional system needs V_{max(d_F, e)} with e the number of solutions; LFD alone does not bound solving cost. The abstract/introduction sentence "LFD ... governs the linear-algebra solving" should carry this qualification.

**Lean / formalization status (what a reviewer can cite against you).**
- Coq: Persson (c. 1999) and Thery, "A machine-checked implementation of Buchberger's algorithm", J. Automated Reasoning 26 (2001) 107-137 (Coq; extraction to OCaml). Isabelle/HOL: Immler-Maletzky AFP entry "Groebner_Bases" (reduction, S-polynomials, Buchberger criteria, reduced GBs), plus Maletzky-Immler, "Groebner bases of modules and Faugere's F4 algorithm in Isabelle/HOL", CICM 2018, arXiv:1805.00304. Mizar: Schwarzweller (2006). Lean 4: Guo-Shen-Liu-Zhi, arXiv:2602.12772 (2026); Lean 3 earlier work by Dermitzakis (as described in that paper); Shen et al. arXiv:2604.13514 (tactics, certificates). Mathlib `MonomialOrder.div` and the combinatorial Nullstellensatz file are what your proof uses.
- No formalization of LFD, FFD or solving degree in any proof assistant turned up in my searches (WebSearch, standard and extended). That is a negative keyword result only; CompPoly, WuProver GroebnerTactic, and Lean repositories were not searched at source level by me.
- Therefore: "first Lean formalization of Groebner bases" must not be implied; "first (that we could find) formalization of the last-fall-degree definition and of an LFD upper bound for an ECDLP-derived family" is the defensible form, with the search described.

**What to cite / soften for (a).**
- Cite PKM Section 3.3 for the chain system; cite CLO Ch. 2 Section 9 and CG Thm 2.9/Example 4.3 for the Groebner-basis and LFD-0 facts; cite Steiner 2024 for the coprime-leading-monomial technique in LFD bounds for iterated systems; cite GMP Thm 11 and HKYY Thm 1.1 for n-independent LFD bounds; cite the Coq/Isabelle/Lean GB work for formalization context.
- Rewrite "The new ingredient is a proof that ... form a reduced-ideal-trivial subsystem" as: "This is the Groebner-basis property of the PKM chain system (coprime leading terms); we give a short root-of-unity proof, kernel-checked, and use it with the finite-quotient bound." Drop "new ingredient".
- State that `max(mk+1, delta)` is a generic Macaulay-type dimension bound and is far from the exact values; say the characteristic-two case is covered classically (the Remark), not only by the formalized route.

## (b) Exact LFD of two-summand prime-field subgroup chains: 4 (k<=5), 5 (k=6..9), >=6 (k=10); dim V_4 crossover

**Directly overlapping prior art.** PKM Table 4 (n1 = 1..11; p the smallest prime with 2^{n1} | p-1 and p >= 2^{2 n1}; random reduced Weierstrass curve and random point, 100 runs per row; Magma Groebner routine; average number of solutions about 1):

| n1 | 1-5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|
| D_av^corr | 4.00 | 5.00 | 5.00 | 5.00 | 5.00 | 6.00 | 6.00 |
| D_av (uncorrected) | 4.00-4.48 | 5.00 | 5.36 | 5.36 | 5.78 | 6.36 | 6.29 |

PKM Table 5 (isogeny-kernel factor base) shows 5 for 2^2..2^9-order kernels and 6 for 2^10, 2^11. PKM observe "a much smaller dependency in n of the degree of regularity" than the semi-generic (3+n)/2 (Bardet) and say the systems are "easier to solve than semi-generic systems" but give no explanation of the steps.

Differences the manuscript can legitimately claim: (i) PKM's number is the Magma-reported maximal step degree (a solver-dependent dreg proxy, with a correction that drops steps with no reduced pair), averaged, with a free choice of p and curve; the manuscript's is the exact iterated-closure LFD certified by rank, over fixed primes/curves/targets, including cosets; (ii) the mechanism (dim V_4 = 4k^2 + 172k - 571 versus binom(2k, <=4)); (iii) the all-depth Lean value 4 for one target family. The agreement of the breakpoints strongly suggests that LFD, not just a solver artefact, is behind PKM's table, and the manuscript should say so explicitly as corroboration.

Arithmetic check of the manuscript's dimension claim (done here): 4k^2+172k-571 gives 389 (k=5), 605 (k=6), 829 (k=7), 1061, 1301, 1549 (k=10); binom(2k,<=4) = 386 (k=5), 794 (k=6): consistent with the text. Note the formula is fitted at k = 6..10 (five points for a three-parameter quadratic) and is called empirical in the manuscript. The k=5 value 389 does not follow the formula's regime (V_4 is the full degree-4 ideal part there), so "the formula crosses the binomial between k=5 and k=6" is a statement about the extrapolated formula, not a derived threshold.

**Other prior art relevant to (b).**
- Kosters-Yeo Section 4-5 and Remark 4.12: in the binary setting, the (Magma) degree of regularity of S3 Weil-descent systems grows with n (n <= 40) and the FFD is 2, so FFD is a poor predictor. Same phenomenon (observed dreg rises while a bound is flat) in a different family. Also HKY CRYPTO 2015 provide evidence against the FFD assumption in char 2.
- Caminata-Gorla Examples 4.1-4.4: separation of FFD, LFD, solving degree and dreg for small systems; Example 4.3 is the k=1 chain.
- Yokoyama-Yasuda-Takahashi-Kogure, J. Math. Cryptol. 14 (2020), complexity bounds on Semaev's naive index calculus for prime fields, discuss degrees of fall in that setting (per the earlier audit).
- Amadori-Pintore-Sala, FFA 51 (2018): single-relation variant of the prime-field Semaev approach; compares with PKM; I did not read their experiments for degree data.
- Huang-Kosters-Yeo / HKYY Thm 1.1 and Section 6: LFD bounds depend on number of solutions s and are not directly available for ECDLP systems without field equations. Relevant to the "solution-count effects" language.
- I found no paper computing exact LFD of PKM systems as k grows, and no published proof or explanation of the 4 -> 5 -> 6 pattern. Absence in search is not proof.

**What to cite / soften for (b).**
- Cite PKM Table 4 and Table 5 as the original observation of the 4/5/6 pattern in the corrected degree of regularity, and say the manuscript's data reproduce its breakpoints at the level of exact LFD.
- Replace "explains" with "gives an empirical account"; the abstract already says "predicted"; keep the claim that the identity for dim V_4 is empirical and that no universal lower bound is proved.
- The k=10 result rests on one unit-ideal instance (the abstract says so). Say "LFD >= 6 for one instance"; PKM's n1 = 10 row (6.00 averaged over 100 runs) is wider evidence for solver behaviour but not for LFD. Do not claim the pattern "grows by one each time" beyond the two observed steps (the manuscript is careful here; keep it).
- Note unit-ideal/solution-count confounds: HKYY Prop 2.8 / Thm 1.1 make the solution count e part of the cost; PKM's runs average about one solution.

## (c) Frobenius-invariant (symmetric) factor bases do not change exact LFD beyond solution-count effects

- GGMP (Galbraith-Granger-Merz-Petit, SAC 2020, LNCS 12804, pp. 115-138): Frobenius-invariant factor bases cut the number of polynomial systems (best case by 1/n), speed up the linear algebra by n^2 for Koblitz, and improve symmetry breaking. Their claims concern relation counts and cost; nothing I could see claims a change in the degree of the individual system. Full text not read (ePrint 403).
- Faugere-Gaudry-Huot-Renault, J. Cryptology 27 (2014): symmetries in the index calculus (automorphism/torsion based). Faugere-Huot-Joux-Renault-Vitse (EUROCRYPT 2014) symmetrized summation polynomials: these change the polynomial (degree and variable count), a different mechanism from Frobenius-stable base.
- Kosters-Yeo Section 5 (reply to Semaev 2015/310 and Karabina 2015/319): with n=45, m=2, t=2 a *specific* (non-random) sub vector space gave dreg 5 where Semaev claimed constant 4, and "the choice of the specific vector space is a good one" at n=40 (dreg 4). So structured versus random spaces were already observed to change solver degree in this family. This is the closest qualitative precedent for the manuscript's "random versus symmetric" comparison and conflicts in spirit with a blanket "no change" reading; the manuscript's finding is restricted to LFD, to its seeded grid, and to strata of equal k and solution count.
- Galbraith-Gebregiyorgis (INDOCRYPT 2014) compare subfield and other factor bases experimentally (the manuscript says so; I did not verify). Semaev 2015 and Karabina 2015 use subspaces of F_2^n with specific structure.
- HKYY Thm 1.1 and Prop 2.8: LFD upper bounds and solving cost depend on the number of solutions; since Frobenius-stable solution sets have more solutions, "effects through solution count" is in line with theory; this is a plausible mechanism, not a new observation.

What to cite / soften: cite GGMP and KY Section 5 for the structured-versus-random context; state (c) as an empirical finding on the sampled grid (the manuscript's statistical language "no residual effect detectable" is appropriate); do not generalize to Koblitz curves beyond the tested k <= 7; make explicit that GGMP's benefit (fewer systems, faster linear algebra) is untouched by this null result on degree.

## Specific edits recommended in the manuscript

1. Related work: correct the Kosters-Yeo description (first fall degree = 2 for ordinary binary S3; dreg growth; NP-completeness) and move "LFD for Weil descent" to HKYY, HKY, Huang-Kosters, Huang 2103.07282, GMP. Correct the PKM sentence about FFD (see headline 4).
2. Add PKM Tables 4-5 as precursor of the exact data (new subsection in Results or related work).
3. Change "The new ingredient" (abstract) and item 1 of the introduction to describe the Groebner-basis property and say the contribution is the formalization plus the explicit LFD statement.
4. Add Steiner (ToSC 2024), GMP, Kousidis-Wiemers, and Caminata-Gorla Thm 2.9/Example 4.3 where the LFD of Groebner-basis subsystems is discussed.
5. Add a sentence on formalization precedents (Coq/Isabelle/Lean GB libraries; Guo et al. 2026) and state that no LFD formalization was found by search.
6. Qualify "LFD governs the linear-algebra solving" with the number-of-solutions requirement (HKYY Prop 2.8).
7. Describe the dim V_4 identity as an empirical fit at k = 6..10 (five data points), and avoid "explains".

## Unverified / needs author check

- Whether Galbraith-Gebregiyorgis compare subfield versus random bases the way the manuscript says.
- Whether the published PKM (Springer) text differs from the author PDF with respect to first fall degree for the prime case.
- Exact bibliographic data flagged "verify" below were taken from search-result metadata rather than the publisher page.
- Nothing here establishes that a statement identical to (a) or (b) is absent from the literature.

## BibTeX-ready additions

```bibtex
@inproceedings{hky15,
  author    = {Ming-Deh A. Huang and Michiel Kosters and Sze Ling Yeo},
  title     = {Last Fall Degree, {HFE}, and {Weil} Descent Attacks on {ECDLP}},
  booktitle = {Advances in Cryptology -- CRYPTO 2015},
  series    = {LNCS},
  volume    = {9215},
  pages     = {581--600},
  year      = {2015},
  note      = {IACR ePrint 2015/573, https://eprint.iacr.org/2015/573; volume/pages/DOI to verify}
}
@article{kw19,
  author  = {Timo Kousidis and Andreas Wiemers},
  title   = {On the first fall degree of summation polynomials},
  journal = {arXiv preprint},
  year    = {2019},
  note    = {arXiv:1906.05594; Thm 3.2, Remark 3.3; check journal version}
}
@article{gmp21,
  author  = {Elisa Gorla and Daniela Mueller and Christophe Petit},
  title   = {Stronger bounds on the cost of computing {Gr\"obner} bases for {HFE} systems},
  journal = {Journal of Symbolic Computation},
  year    = {2021},
  doi     = {10.1016/j.jsc.2020.07.011},
  note    = {arXiv:2011.01050; Theorem 11, Corollary 1; volume/pages to verify}
}
@article{steiner24,
  author  = {Matthias Johann Steiner},
  title   = {Solving Degree Bounds For Iterated Polynomial Systems},
  journal = {IACR Transactions on Symmetric Cryptology},
  volume  = {2024},
  number  = {1},
  pages   = {357--411},
  year    = {2024},
  note    = {arXiv:2310.03637}
}
@inproceedings{ggmp21,
  author    = {Steven D. Galbraith and Robert Granger and Simon-Philipp Merz and Christophe Petit},
  title     = {On Index Calculus Algorithms for Subfield Curves},
  booktitle = {Selected Areas in Cryptography -- SAC 2020},
  series    = {LNCS},
  volume    = {12804},
  pages     = {115--138},
  year      = {2021},
  doi       = {10.1007/978-3-030-81652-0_5},
  note      = {ePrint 2020/1315}
}
@article{aps18,
  author  = {Alessandro Amadori and Federico Pintore and Massimiliano Sala},
  title   = {On the discrete logarithm problem for prime-field elliptic curves},
  journal = {Finite Fields and Their Applications},
  volume  = {51},
  pages   = {168--182},
  year    = {2018},
  doi     = {10.1016/j.ffa.2018.01.009},
  note    = {ePrint 2017/609}
}
@misc{semaev15,
  author = {Igor A. Semaev},
  title  = {New algorithm for the discrete logarithm problem on elliptic curves},
  year   = {2015},
  note   = {IACR ePrint 2015/310; arXiv:1504.01175}
}
@inproceedings{karabina15,
  author    = {Koray Karabina},
  title     = {Point Decomposition Problem in Binary Elliptic Curves},
  booktitle = {ICISC 2015},
  series    = {LNCS},
  volume    = {9558},
  year      = {2016},
  doi       = {10.1007/978-3-319-30840-1_10},
  note      = {ePrint 2015/319}
}
@article{fghr14,
  author  = {Jean-Charles Faug{\`e}re and Pierrick Gaudry and Louise Huot and Guénaël Renault},
  title   = {Using Symmetries in the Index Calculus for Elliptic Curves Discrete Logarithm},
  journal = {Journal of Cryptology},
  volume  = {27},
  pages   = {595--635},
  year    = {2014},
  doi     = {10.1007/s00145-013-9158-5}
}
@article{hps14,
  author  = {Timothy J. Hodges and Christophe Petit and Jacob Schlather},
  title   = {First fall degree and {Weil} descent},
  journal = {Finite Fields and Their Applications},
  year    = {2014},
  note    = {ScienceDirect S1071579714000872; volume/pages to verify}
}
@misc{huang21,
  author = {Ming-Deh A. Huang},
  title  = {On the last fall degree of Weil descent polynomial systems},
  year   = {2021},
  note   = {arXiv:2103.07282}
}
@article{huang-semilocal,
  author  = {Ming-Deh A. Huang},
  title   = {Last fall degree of semi-local polynomial systems},
  year    = {2026},
  doi     = {10.1137/21M1466293},
  note    = {arXiv:2311.02804; journal name/volume to verify}
}
@article{yytk20,
  author  = {Kazuhiro Yokoyama and Masaya Yasuda and Yasushi Takahashi and Jun Kogure},
  title   = {Complexity bounds on {Semaev}'s naive index calculus method for {ECDLP}},
  journal = {Journal of Mathematical Cryptology},
  year    = {2020},
  doi     = {10.1515/jmc-2019-0029}
}
@article{semaev-tenti21,
  author  = {Igor Semaev and Andrea Tenti},
  title   = {Probabilistic analysis on {Macaulay} matrices over finite fields and complexity of constructing {Gr\"obner} bases},
  journal = {Journal of Algebra},
  volume  = {565},
  pages   = {651--674},
  year    = {2021},
  note    = {cited as [20, Thm 2.1] in Caminata-Gorla; verify}
}
@inproceedings{pq12,
  author    = {Christophe Petit and Jean-Jacques Quisquater},
  title     = {On Polynomial Systems Arising from a {Weil} Descent},
  booktitle = {ASIACRYPT 2012},
  series    = {LNCS},
  volume    = {7658},
  year      = {2012}
}
@book{clo15,
  author    = {David A. Cox and John Little and Donal O'Shea},
  title     = {Ideals, Varieties, and Algorithms},
  edition   = {4},
  publisher = {Springer},
  year      = {2015},
  note      = {Ch. 2 Sec. 6 (Groebner basis characterization), Ch. 2 Sec. 9 (coprime leading terms)}
}
@article{thery01,
  author  = {Laurent Th{\'e}ry},
  title   = {A Machine-Checked Implementation of {Buchberger}'s Algorithm},
  journal = {Journal of Automated Reasoning},
  volume  = {26},
  number  = {2},
  pages   = {107--137},
  year    = {2001},
  note    = {Coq; volume/pages to verify}
}
@inproceedings{maletzky-immler18,
  author    = {Alexander Maletzky and Fabian Immler},
  title     = {Gr{\"o}bner Bases of Modules and {Faug{\`e}re}'s {F4} Algorithm in {Isabelle/HOL}},
  booktitle = {Intelligent Computer Mathematics (CICM 2018)},
  year      = {2018},
  note      = {arXiv:1805.00304; AFP entries Groebner\_Bases and Signature\_Groebner}
}
@misc{guo26,
  author = {Guo and Shen and Liu and Zhi},
  title  = {Formalizing {Gr\"obner} Basis Theory in {Lean}},
  year   = {2026},
  note   = {arXiv:2602.12772; first-name initials to fill from the arXiv record}
}
@misc{shen26,
  author = {Shen and Guo and Liu and Zhi},
  title  = {Automated Tactics for Polynomial Reasoning in {Lean} 4},
  year   = {2026},
  note   = {arXiv:2604.13514; first-name initials to fill}
}
```

Already in the manuscript bibliography (no change needed except the descriptions above): hkyy, cg, pkm, fppr, semaev (2004), fhjrv, gg, ky.
