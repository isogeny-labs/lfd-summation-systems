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

- [ ] Execute five LFD follow-ups: derivation certificates, lift review, split closure, exceptional primes, depth-18 test (Codex, started 2026-11-17, branch research/lfd-followups-20261117, PR #31; in progress 2026-11-18: polynomial/FEI first-derivative bridge)

- [x] Finalize and locally validate both depth-14 replications and submission provenance; final CI and merge tracked in PR #24 (Codex, started 2026-11-16, branch paper/finalize-pr24-20261116, PR #24)

<!-- One line per task, added as the task's first step with its draft PR:
- [ ] <task> (<runtime>, started YYYY-MM-DD, branch <branch>, PR #<n>) -->

## Blocked / questions for the owner

- Local-only constraint: full depth-14 certificate and direct depth-18 D7/D8 runs exceed current RAM/scratch capacity; local probes finished with explicit caps; full depth-14 certificate, split D8, and depth-18 breakpoint remain open (PR #31).
