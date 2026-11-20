# Style audit of `curve-lfd-bounds.tex` (20 November 2026)

Scope: the manuscript `paper/curve-lfd-bounds.tex` as on `origin/main` at
`1ff74a7` (4283 lines, 35,161 words of source, 62 pages). Line numbers below
refer to that version. Counts exclude LaTeX comment lines. The abstract
(lines 41-60) was approved separately and is out of scope.

The mathematics, the numbers and the evidence distinctions are not in
question here. The problem is the prose: the paper reads as a log of runs
and review rounds layered on top of a mathematics paper. The patterns below
are listed in order of how much they cost the reader.

## 1. Run logistics in the body

The body carries the operational detail of individual runs, which a reader
of the mathematics does not need and which belongs in the data appendix.

| Pattern | Count | Examples |
|---|---|---|
| Wall times (`~s`, `~h`, CPU-hours) | 12 on 10 lines | 1697 "took $28{,}747$~s on $32$ shared virtual CPUs"; 1746 "The run took $26{,}723$~s ... $18{,}495$~s of it in the check"; 1779 "about $3.7$~h ($13{,}374$~s)"; 1831; 1870; 1434 "about $52$ CPU-hours" |
| Memory (KiB, GiB, MB, GB, TB, "resident set") | 27 on 17 lines | 1698 "peak resident set of $23{,}129{,}176$~KiB (about $22$~GiB)"; 1747; 1780 "$234{,}285$~MB"; 1811; 1833; 1856; 1871; 1883 "about $540$~GB ... $1.2$~TB" |
| Hosts, CPUs, GPUs, restarts | 15 lines | 1701 "also on a second host"; 1811 "on $32$ shared virtual CPUs of a loaded host"; 1857 "after two host restarts the run was resumed from checkpoints"; 3674 "a provider restart interrupted the first attempt" |
| Hashes and SHA-256 prefixes | 33 on 32 lines | 1748 "(SHA-256 \code{06fdd4e4\ldots})"; 1765 "binary SHA-256 \code{6824a8bd\ldots}"; 1770; 1812 "\code{3d1f2996\ldots}"; 3523 "basis hash \code{d4fe1334\ldots}" |
| File and directory paths in running text (`\code{...}`) | 308 on 102 lines (most in the appendix tables, where they belong) | body: 1668, 1721-1724, 1735-1737, 1766, 1833, 1874, 2854, 2859, 2972, 2988 |
| Engine names, binary names, job names | — | 1721 "\code{hky\_fp\_tower\_fast2.rs}, with the instance loader of its generator \code{pkm\_tower\_pilot.rs}"; job names J1-J4s appear in comments 1708-1848 and in the data table |
| Machine-checked counts of internal operations | — | 1745 "all $4{,}120{,}873$ reductions gave zero"; 1778 "$5.6\cdot10^{14}$ multiply-adds" |

## 2. Process narration

The text tells the history of how results were obtained rather than stating
them.

| Pattern | Count | Examples |
|---|---|---|
| "later", "first version", "had previously", "since kernel-checked", "after the snapshot" | 19 lines | 1037 "an independent check of $292$ random generators made later"; 1106 "which had previously been checked only numerically"; 1720 "A later version of the second engine"; 2031 "Later readings of one engine"; 2369 "where the upper bound $7$ was later obtained by the lift"; 3630 "(since kernel-checked, Theorem~\ref{thm:bridge})"; 3717 "Its first version failed on 29 of the 40" |
| Prediction/preregistration bookkeeping in the body | 13 lines | 1890, 1912, 1986, 2189-2191 (a disclosed exception about smoke tests), 2223, 2236, 2380 "That the predictions preceded the runs rests on the dates recorded in the files" |
| Process sequences ("Before the run, the same binary reproduced ...", "It first checked that ...", "The same lift program was then run") | — | 1731-1737, 1770-1772, 1824-1826 |
| Review history inside results | 6 lines | 1030-1040 "A review of the written proof by a separate AI session ... made before the formalization"; 2475 |

## 3. The same evidence hedge restated sentence by sentence

The paper defines evidence labels in the introduction (lines 82-93), then
restates the status of the same result in almost every sentence that touches
it.

| Hedge | Count |
|---|---|
| "not certified", "not a certificate", "corroborat-" | 29 on 23 lines |
| "single run", "one run", "run once", "one engine", "single-engine" | 37 on 32 lines |
| "unverified" | 3 |
| "kernel-checked" | 77 lines |
| "written proof" | 31 lines |

The $k=14$ upper bound alone ("lifted", "no independent certificate",
"rests on the lift") is qualified in the abstract, the evidence labels, the
contributions list, Figures `fig:fall` and `fig:pipeline`, Section `sec:computations` (six times), the caption of Table `tab:profile`,
the conjecture paragraph, two Limitations bullets (lines 3067-3104 and
3105-3122, together 56 lines), Future directions and the conclusion. The
word "lift" occurs 45 times in Section `sec:computations` and 17 times in Limitations;
"second engine" 19 times in Section `sec:computations` and 14 times in Limitations.

## 4. Redundant restatement across sections

Each of these results is stated in full in at least four places:

- $\mathrm{LFD}\le k+2$: introduction (contributions and "Why this
  matters"), Section `sec:bounds` roadmap, Theorem `thm:selector` and its afterword,
  Table `tab:twosided`, Section `sec:solving`, Conclusion.
- The profile 4, 5, 6, (7): contributions, the box of Figure `fig:degrees`,
  the last paragraph of `sec:lfd`, the comparison in `sec:related`,
  `sec:depth10`, the summary paragraph of `sec:depth14`, Table `tab:profile`,
  the conjecture paragraph, `sec:solving`, Conclusion.
- The growth slope $1/14$: contributions, "Why this matters", Theorem `thm:growth`
  afterword, "Scope" paragraph, the end of `app:growthproof`, Conclusion.
- Second-engine runs at $k=14$-$16$: `sec:depth14` (five paragraphs), the $k=14$ row of `tab:profile`,
  Limitations (three bullets), Future directions.

## 5. Sentence and paragraph shape

- Semicolons: 537 in the body. Typical chain (line 1040): "The bound is
  attained at $k=2$, where ... ; for unit ideals ... is kernel-checked. For
  $k\ge3$ no computed system attains it (...), so it is not claimed to be
  sharp. For singular ..."
- Stacked parentheticals: 1479 opening parentheses in the body. Example
  (line 1115): "(one per block) on the levels ... (in the notation of ...)".
- 143 sentences longer than 45 words (math counted as one word).
- Paragraphs above 250 words: 27. The longest prose blocks are the
  Limitations list (1762 words as one list), Future directions (935), the
  "second engine on a fourth prime" paragraph (706), the "Observations"
  paragraph of the sweep (374) and the depth-step paragraph (350).
- Bullet-like prose: the "Status" paragraph of `sec:growth` and the Limitations
  bullets are lists of clauses joined by semicolons rather than arguments.

## 6. Meta-commentary and boilerplate

The usual filler words ("we note that", "importantly", "crucially",
"notably", "it is worth") do not occur. The equivalent here is defensive
meta-text:

- "This paper neither claims nor verifies these statements, and we do not
  list them as bounds" (line 983).
- "This search does not establish that none exists" (line 588).
- "These data refute the tested solution-count law; they do not rule out
  other bounds in the solution count" (line 2229).
- "These observations distinguish coordinate-space classes on the computed
  instances; they do not exclude other universal upper bounds" (caption of Figure `fig:binarydeep`).
- "Within these limits, ..." closing paragraph of Limitations, which
  repeats the evidence legend.

## 7. Terminology drift

- The same programs are called "the reference engine", "the third engine",
  "a closure engine of another repository", "the second engine", "a later
  version of the second engine", "the faster closure engine of the same
  repository", "the same binary".
- "lift", "degree lift", "first lift", "second implementation of the
  lift", "third lift", "the same lift program".
- "unit instance", "the unit instance over $p=65537$", "both unit
  instances", "the paper's own unit instance".
- "certified", "exact", "certified exact", "exact computation", "exact
  (one instance, one engine)".

## 8. Narrative

The introduction lists results but does not state the question they answer.
The result sections often open with the machinery ("The dense chain closure
does not reach $k=8$; a shift engine ... does") and give the finding at the
end. The conclusion summarises rather than saying what has changed and what
the sharpest open question is.

## Plan

1. Introduction: open with the question (does the solver degree stay
   bounded as the factor base grows?), what was believed, what we find,
   and what it does and does not mean. Contributions list kept, shortened.
2. One evidence legend, defined once; each status stated once where the
   result is introduced; the $k=14$ qualifications collected in one place
   (`sec:depth14`) and referenced elsewhere.
3. Run logistics (times, memory, hosts, hashes, binaries, file paths, job
   names) moved to a new "Computational record" subsection of the data
   appendix; the body keeps instance, computation, result and status.
4. Process narration removed from the body.
5. Sections lead with the finding; bridging sentences at section starts.
6. Limitations and Future directions rewritten as short paragraphs without
   repeating `sec:computations`.
7. No change to any theorem-environment body, table, figure, label,
   citation, number or evidence status. Abstract untouched.

## Result of the rewrite

Before: 62 pages, 34,311 words of non-comment source, 41,960 words of
extracted PDF text. After: 60 pages, 33,082 words of source, 40,521 words of
PDF text. The page target (15-25% fewer) was not reached: about 45% of the
pages are theorem and proof environments, tables, figures, the generated
sweep tables, the Lean declaration table, the growth proof and the
bibliography, none of which may change, and the hard rule that every number
of the old body survives leaves little prose that can be deleted rather than
rewritten.

Pattern counts in the body (before the appendices, comments excluded):

| Pattern (body, before the appendices) | Before | After |
|---|---|---|
| wall times/memory/hosts in body | 46 | 0 |
| SHA-256 / hash prefixes in body | 8 | 0 |
| file paths (\code) in body | 22 | 3 |
| "not certified"/"not a certificate"/"corroborat-" | 29 | 19 |
| "single run"/"one run"/"run once"/"one engine" | 37 | 13 |
| "later"/"before the run"/"first version" | 11 | 1 |
| semicolons | 451 | 407 |

The three remaining `\code{}` uses in the body are inside the proof of
Theorem `thm:sd` (unchanged by rule), the prediction name `d6-base-15`, and
the caption of Table `tab:fermatlfd` (unchanged by rule). The remaining
"one run"/"not certified" occurrences are mostly inside Table `tab:profile`
(unchanged by rule) and in the single place where each status is stated.

New section outline: Introduction (question, prior belief, findings,
evidence labels, contributions, organization); Background; Systems and
notation; Proved bounds (the ceiling $k+2$ and the floors); Exact
computations (sweep; depths 8-9; the step at depth 6; depth 10; depths
11-16 and the frontier at $k=14$, with the degree lift and the
corroboration collected in one place; the recursion; the profile and the
conjecture, moved here from the end of the section; three to eight summands;
binary and Koblitz curves); Growth in the depth; The solving degree; What
the results mean, and what they do not (new); Limitations; Future
directions; Conclusion; appendices, with a new subsection "Computational
record" (`app:record`) holding the run logistics moved out of the body.

## Invariants

Mechanical comparison of `origin/main` (`1ff74a7`) against this branch
(scratch script, not committed). Theorem-type environments, proofs, tables,
figures and equations are compared as byte-identical multisets of
`\begin{...}...\end{...}` blocks; numbers are all digit tokens of the old
non-comment source (thousands separators normalised, bibliography
excluded), searched in the whole new file.

```
labels old 111 new 113 removed [] added ['app:record', 'sec:meaning']
cite keys old 28 new 28 removed []
bibitems identical True
refs to undefined labels in new []
env theorem: old 16 new 16 identical-as-multiset True
env lemma: old 11 new 11 identical-as-multiset True
env corollary: old 3 new 3 identical-as-multiset True
env proposition: old 2 new 2 identical-as-multiset True
env conjecture: old 1 new 1 identical-as-multiset True
env definition: old 3 new 3 identical-as-multiset True
env example: old 1 new 1 identical-as-multiset True
env remark: old 3 new 3 identical-as-multiset True
env proof: old 28 new 28 identical-as-multiset True
env table: old 10 new 10 identical-as-multiset True
env figure: old 10 new 10 identical-as-multiset True
env equation: old 4 new 4 identical-as-multiset True
env thebibliography: old 1 new 1 identical-as-multiset True
abstract identical True
AI-use paragraph identical True
distinct numeric tokens in old body 518 missing from new file 0 []
of which missing from new body (non-comment) []
```

Exceptions: none. No label or citation was removed; two labels were added
(`sec:meaning`, `app:record`). Every theorem, lemma, corollary,
proposition, conjecture, definition, example and remark body, every proof,
every table and every figure is byte-identical (the subsection
`sec:conjecture`, with Table `tab:profile` and Figure `fig:lfddepth`, moved
as a block). The abstract and the "Use of AI tools" paragraph are
byte-identical. All 518 distinct numeric tokens of the old body occur in
the new body (outside comments).

Build: Tectonic 0.17, producer xdvipdfmx, 0 warnings, 0 overfull boxes,
no "??" in the text. `scripts/check_publication_artifacts.py`: passed,
439 checks.

## Second pass (20 November 2026)

A sharper pass on top of the first rewrite (branch head `945ea83`), with the
same hard rules: no theorem/lemma/definition/proposition/conjecture body,
proof, table value, figure, `\label` or `\cite` changed; abstract and the
"Use of AI tools" paragraph byte-identical; every numeric value of the
pre-pass body still present somewhere. The first pass had removed the worst
slop and added the story arc and `sec:meaning` but only reached 60 pages and
left some mechanical reporting prose.

### What was condensed or cut (prose only; all numbers preserved)

- **Larger split-cell readings m=6–9** (`sec:msummands`, "Larger cells
  (corroboration)"): multi-clause narration of the not-certified readings
  replaced by one compact table (`tab:splitlarge`, 7 cells in a 3-wide
  layout) plus two sentences; all column counts and readings moved into the
  table/caption. ~20 prose lines → table + 2 sentences.
- **Depth-11-to-16 corroboration** (`sec:depth14`): the intro, "A fourth
  prime", "The two unit instances", "A third lift", "Degree six at depths
  fifteen and sixteen", "The input of the lift", "A random curve",
  "Two more primes" paragraphs tightened sentence by sentence (dropped the
  engine-hardware aside; merged restated "same as the unit instance"
  clauses); every dimension, profile, part rank and the evidence-status
  reasoning kept. ~16 lines.
- **Decomposition without a Gröbner basis** (`sec:future`): derivation prose
  tightened to its essential claim; all measurements (293×, 94 s, 1.85 µs,
  one core-year, 110→6) kept. ~5 lines.
- **Exploratory items**: "Towards genericity" (`sec:depthsix`), "Generic
  closure dimensions" and "What remains" (`sec:blockfour`,
  `sec:growthnontorsion`), "Exact closures with four and five summands"
  (`sec:msummands`) — connective/analysis prose trimmed, numbers kept. ~8
  lines.
- **Cross-section de-duplication**: the depth-ten profile restatement
  "4 for k=2..5, 5 for k=6..9, 6 at k=10" replaced by a reference to
  `tab:profile`; the conclusion and the depth-14 engine-roadmap sentence
  tightened (each result now stated fully once in the body and referenced).
  ~12 lines.
- **Appendix padding** (`app:engines`, `app:checks`, `app:data`,
  `app:record`): the third-engine and reference-engine validation prose, the
  proof-review and second-lift notes, the data-intro, and the per-job record
  bullets (factored out the repeated "32 shared virtual CPUs"/"peak resident
  set" framing) tightened; all times, memory figures and hashes kept. ~12
  lines.

### Before / after

| | original `1ff74a7` | first pass `945ea83` | second pass |
|---|---|---|---|
| pages | 62 | 60 | **59** |
| PDF-text words | 41,960 | 40,521 | **39,835** |
| source words (non-comment) | 34,311 | 33,082 | **32,440** |
| non-comment semicolons | — | 503 | 497 |

The page count floor is 59: ~45% of the pages are theorem/proof
environments, the 10 figures, the 11 tables (incl. the generated sweep and
Lean-declaration tables), the growth proof and the bibliography, none of
which may change, and the rule that every number of the pre-pass body
survives caps prose-only cuts. A further page reduction needs a decision to
remove a *result* (candidates below), which was not taken.

### Result-removal candidates NOT acted on (owner decision)

These would each save roughly a third to a half page but remove reported
content, so they were left in place:

1. The not-certified larger split-cell readings m=6–9 (`tab:splitlarge`) —
   external corroboration, not a result of this paper; could be dropped
   rather than tabulated.
2. The "decomposition without a Gröbner basis" item (`sec:future`) — an
   exploratory benchmark unrelated to the LFD results.
3. The "Blocks of four: conjectural remarks" subsection (`sec:blockfour`,
   "nothing in it is proved") — could move to a companion note.
4. The depth-six "Towards genericity" written-analysis enumerate and the
   binary/Koblitz statistical "Observations" (associations on small strata,
   not causal) — could be cut to their headline claims.
5. One of the four depth-14 corroborating instances (the certified
   `LFD ≥ 7` is established on four instances over three primes; reporting
   fewer would thin the evidence and was not done).

### Invariants (scratch script, not committed; `945ea83` vs this pass)

```
env theorem:   old 16 new 16 identical True
env lemma:     old 11 new 11 identical True
env corollary: old 3  new 3  identical True
env proposition: old 2 new 2 identical True
env conjecture: old 1 new 1  identical True
env definition: old 3 new 3  identical True
env example:   old 1  new 1  identical True
env remark:    old 3  new 3  identical True
env proof:     old 28 new 28 identical True
env table:     old 10 new 11 (added tab:splitlarge; none removed or altered)
env figure:    old 10 new 10 identical True
env equation:  old 4  new 4  identical True
env thebibliography: identical True
labels old 113 new 114  removed []  added ['tab:splitlarge']
cite keys old 26 new 26  removed []  added []
refs to undefined labels in new: []
abstract identical: True
AI-use paragraph identical: True
distinct numeric token-values in old body: 543; missing from new file: 0
```

Exceptions: one table added (`tab:splitlarge`), holding the relocated m=6–9
readings; no existing table, figure, theorem, proof or equation changed. No
label or citation removed. Abstract and AI-use paragraph byte-identical.
All 543 distinct numeric tokens of the pre-pass body occur in the new file.

Build: Tectonic 0.17.0, producer xdvipdfmx, 0 warnings, 0 overfull boxes,
no "??" in the text; 59 pages A4. `scripts/check_publication_artifacts.py`:
passed, 439 checks. `EPRINT.md` page count and PDF/tex digests updated.
