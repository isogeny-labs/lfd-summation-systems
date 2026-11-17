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

- [x] Finalize and locally validate both depth-14 replications and submission provenance; final CI and merge tracked in PR #24 (Codex, started 2026-11-16, branch paper/finalize-pr24-20261116, PR #24)

- [x] Manuscript cleanup before ePrint: full growth proof in an appendix, shorter abstract, AI-assistance disclosure, no revision-history wording (Claude, started 2026-11-17, branch claude/ecstatic-dirac-7f5qyt, PR #27)

<!-- One line per task, added as the task's first step with its draft PR:
- [ ] <task> (<runtime>, started YYYY-MM-DD, branch <branch>, PR #<n>) -->

## Blocked / questions for the owner

- PR #27: the author must confirm the wording of the "Use of AI tools" paragraph (it names Claude and Codex) before submission.
