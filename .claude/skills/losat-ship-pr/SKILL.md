---
name: losat-ship-pr
description: Open, watch, and finish a LOSAT pull request - the base and head branches of each track, the body, the staged-size check, which CI runs, one watcher per session, reading failing logs, and when a session may merge. Use when a branch is ready for a pull request, while its CI runs, and when it merges.
---

# Ship a LOSAT pull request

This skill covers the procedure; the session instruction decides when a PR is due and what it
contains.

## Branches

- LOSAT Web campaign: the engine branch `feature/losat-web-gui` goes to `main` when an
  instruction asks for "`main` への PR"; the app branch `feature/losat-web-gui-app` is merged
  into `feature/losat-web-gui` at the end of each app session (README rule 1), not into `main`.
- Other work: a branch from `origin/main` (`fix/`, `feature/`, `docs/`, `audit/`, `ci/`) with a PR
  to `main`.

## Before opening

1. Run the full tier of skill `losat-gates` for what the branch touches, and record it in the gate
   record.
2. Check staged and committed sizes: no file over about 1 MB without the Owner's consent (skill
   `losat-worktree`); `git diff --stat origin/main...HEAD | tail -1` for the size of the PR.
3. Write the body to a file in `$TASK_DIR`: what changed and why, the stages and DW decisions it
   covers, the gate record paths and results, the evidence hashes, the residuals, and
   Owner-delegated choices. Then `gh pr create --base main --head <branch> --body-file <file>`.

## Watch CI

- `ci.yml` (`rust`, `fast output regressions`) runs on every PR; `web.yml` runs when `web/**`,
  `LOSAT/src/**`, `LOSAT/Cargo.*`, or the WASI test files changed. Expect 5-20 minutes.
- Run one watcher per session, not one per agent:
  `python3 .claude/skills/losat-ship-pr/scripts/watch_prs.py --prs <n>` under Monitor or in the
  background, and act on its event lines (merged, closed, conflict, failed, cancelled, green).
  For long runs (nightly, dispatch workflows) check back at most every 30 minutes.
- Before acting on an event (merge, update, report done), re-read the PR:
  `gh pr view <n> --json state,mergeStateStatus,statusCheckRollup`. The Owner or another session
  may have merged or changed it.
- Read only failing logs, into a file: `gh run view <id> --log-failed > "$TASK_DIR/logs/<id>.log"`,
  and quote the failing lines. A run whose jobs passed but that ended `cancelled` is recovered
  with `gh run rerun <id> --failed`.

## Merge

- Merge only when the session instruction says the session merges, with CI green:
  `gh pr merge <n> --merge --match-head-commit <full sha>`.
- A PR that changes agent instructions, skills, settings, product decisions (`docs/product_decisions/`),
  or `AGENTS.md` waits for the Owner's merge.

## After the merge

Record the merge commit in the gate record and `STATE.md` (skill `losat-campaign`); clean up
task worktrees (skill `losat-worktree`).
