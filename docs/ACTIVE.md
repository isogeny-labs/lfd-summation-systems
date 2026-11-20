# Active goal and items

Owner-pinned work queue. `tools/goal_brief.py` prints this file at the start
of every Claude and Codex session (SessionStart hooks) and for `/goal`
(AGENTS.md, "Goal: what to work on next"). When this file has no unchecked
item, sessions fall back to the newest "Open threads" in `docs/session-logs/`.

Keep it short: one goal, a few ordered items. Sessions add their task here as
the first step of the task (AGENTS.md, "First step of every task"), keep the
line current, and tick it when done; history belongs in the session log or the
PR, not in this file.

## Current goal

_(none pinned — set one with `/goal <text>`)_

## Active items (top = next)
- [ ] Paper: style pass on curve-lfd-bounds.tex (remove agent-output patterns, move run logistics to the appendix, cut repetition, narrative arc, "What the results mean" section; no mathematical, numeric or evidence-status change), with a style audit and invariant check (Claude, started 2026-11-20, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #40; draft, awaiting owner review)
- [x] Paper: replace the abstract with the owner-approved text, resync the ePrint abstract and sheet, list the kernel-checked bounds in proof-status.json, rebuild with Tectonic (Claude, started 2026-11-20, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #39)
- [x] Improve the LFD manuscript for clarity and readability: expanded background on the FFD/LFD/solving-degree problem, three new conceptual figures (index calculus at a glance; the three degrees at a glance; the degree-fall timeline), an intuitive §2.1 opening and a decompressed introduction; extended 2026-11-20 at the owner's request: abstract rewritten for readability (LFD gloss and the first-fall-degree contrast, both already proved in the paper) with the submission sheet and abstract file resynced, roadmaps for Sections 4-6 and the depth-11-to-16 evidence, and sentence-level splits in the limitations, future directions and conclusion; no bound, label or citation change; build gate and 439 publication checks pass (Claude, started 2026-11-19, branch claude/lfd-paper-clarity-visuals, PR #38)
- [x] Consolidate the second-engine paper material from crypto PR #4735 here: run records, lift crate, standalone engine extract gated locally, paper cites lfd paths (Claude, started 2026-11-19, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #37)
- [x] Merge the open PRs (#33, #35), fold the E6 proper-ideal theorem into the manuscript, rebuild and reverify the submission artifacts (opencode, started 2026-11-19, branch paper/finalize-submission-20261119, PR #36)
- [x] Paper: second engine at depth 14 on both unit instances (closures + lifts), V_6 at depth 15 and 16 against the frozen recursion (Claude, started 2026-11-19, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #34)

- [ ] Complete depth-14 input certificate, split degree-8 computation and depth-18 breakpoint, local-only (Codex, started 2026-11-18, branch research/lfd-completion-20261118, PR #35; merged as a record 2026-11-19, the three obligations remain open, see Blocked/questions)
- [x] Prove E6 proper-ideal growth via polynomial first-derivative invariant (Codex, started 2026-11-18, branch formal/e6-proper-ideal-growth-20261118, PR #33)
- [x] Reconcile and merge the open submission PRs (#27, #28, #30, #31); rebuild and verify the manuscript (opencode, started 2026-11-18, branch paper/finalize-submission-20261118, PR #32)
- [x] Paper update: second-engine depth-14 closures, depth-15/16 values, k>=17 out-of-reach statement, cross-engine evidence (claude, started 2026-11-18, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #30)
- [x] GPU path for split-presentation closures with m >= 4 summands: dense device engine built and gated on CPU; verdict: it does not beat the sparse CPU engine; no GPU run made (claude, started 2026-11-17, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #30)
- [x] Execute five LFD follow-ups: derivation certificates, lift review, split closure, exceptional primes, depth-18 test (Codex, started 2026-11-17, branch research/lfd-followups-20261117, PR #31; local-only results merged, remaining obligations stay in the PR)

- [x] Finalize and locally validate both depth-14 replications and submission provenance; final CI and merge tracked in PR #24 (Codex, started 2026-11-16, branch paper/finalize-pr24-20261116, PR #24)
- [x] Build and benchmark a root-finding solver for index-calculus summation systems to replace F4/F5 (resultant/GCD elimination vs Groebner baseline, measured) (Claude, started 2026-11-17, branch claude/index-calculus-root-finding-e95kd6, PR #28)
- [x] Multi-word root-finding kernel in Rust (no Python in the execution path) and a measurement at n = 131 on the ECC2K-130 curve (Claude, started 2026-11-17, branch claude/index-calculus-root-finding-e95kd6, PR #28)
- [x] Update the LFD paper: independent review and fixes of the 16 Nov sections (rational prime-field growth, solving degree, function field and pole orders), a root-finding remark, and refreshed status/provenance (Claude, started 2026-11-18, branch claude/index-calculus-root-finding-e95kd6, PR #28)
- [x] Standardize the repository CI check contract (Codex, started 2026-11-17, branch codex/standard-ci-check, PR #29)

- [x] Manuscript cleanup before ePrint: full growth proof in an appendix, shorter abstract, AI-assistance disclosure, no revision-history wording (Claude, started 2026-11-17, branch claude/ecstatic-dirac-7f5qyt, PR #27)
- [ ] Verify the theorem, definition, table and section numbers cited in the LFD paper against the primary sources, and fix any that are wrong (Claude, started 2026-11-20, branch claude/index-calculus-root-finding-e95kd6, PR #42)

<!-- One line per task, added as the task's first step with its draft PR:
- [ ] <task> (<runtime>, started YYYY-MM-DD, branch <branch>, PR #<n>) -->

## Blocked / questions for the owner


- PR #27: the author must confirm the wording of the "Use of AI tools" paragraph (it names Claude and Codex) before submission.
- PR #31 / PR #35: full depth-14 certificate, split D8, and depth-18 breakpoint remain open; local probes finished with explicit caps (RAM/scratch limits). PR #35 was merged on 2026-11-19 as a record of the capped searches; the merge does not claim any of the three obligations complete.
