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

- [ ] Paper update: second-engine depth-14 closures, depth-15/16 values, k>=17 out-of-reach statement, cross-engine evidence (claude, started 2026-11-18, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #30)
- [x] GPU path for split-presentation closures with m >= 4 summands: dense device engine built and gated on CPU; verdict: it does not beat the sparse CPU engine; no GPU run made (claude, started 2026-11-17, branch claude/lfd-ffd-paper-eprint-status-qq8fit, PR #30)
- [x] Finalize and locally validate both depth-14 replications and submission provenance; final CI and merge tracked in PR #24 (Codex, started 2026-11-16, branch paper/finalize-pr24-20261116, PR #24)

<!-- One line per task, added as the task's first step with its draft PR:
- [ ] <task> (<runtime>, started YYYY-MM-DD, branch <branch>, PR #<n>) -->

## Blocked / questions for the owner
