---
name: losat-ship-pr
description: Open, watch, and finish a LOSAT pull request - the base and head branches of each track, the body, the staged-size check, which CI runs, one watcher per session, reading failing logs, and when a session may merge. Use when a branch is ready for a pull request, while its CI runs, and when it merges.
---

# Ship a LOSAT pull request

This skill covers the procedure; the session instruction decides when a PR is due and what it
contains.

## Branches

`dev` is the integration branch (Owner, 2026-10-17, as in gbdraw). Work pull requests target
`dev`, run the fast checks and merge there; `dev` goes to `main` in one pull request.

- Work: a branch from `origin/dev` (`fix/`, `feature/`, `perf/`, `docs/`, `audit/`, `ci/`) with a
  PR to `dev`.
- LOSAT Web campaign: the engine branch `feature/losat-web-gui` goes to `dev` when an
  instruction asks for "`main` への PR"; the app branch `feature/losat-web-gui-app` is merged
  into `feature/losat-web-gui` at the end of each app session (README rule 1), not into `dev`.
- Promotion: when the Owner asks, a branch `promote/<date>` at the `dev` commit to promote and
  one PR from it to `main`. It contains no change of its own (not `dev` itself as the head: each
  merge into `dev` would restart its Platform checks); they run on it once for the batch, and
  auto-merge merges it when they pass.

## Before opening

1. Run the full tier of skill `losat-gates` for what the branch touches, and record it in the gate
   record.
2. Check staged and committed sizes: no file over about 1 MB without the Owner's consent (skill
   `losat-worktree`); `git diff --stat origin/dev...HEAD | tail -1` for the size of the PR.
3. Write the body to a file in `$TASK_DIR`: what changed and why, the stages and DW decisions it
   covers, the gate record paths and results, the evidence hashes, the residuals, and
   Owner-delegated choices. Then `gh pr create --base dev --head <branch> --body-file <file>`.
4. Open the PR as a draft (`--draft`) while its content may still change: marking it ready turns
   on auto-merge into `dev`.

## Watch CI

- `ci.yml` (`rust`, `fast output regressions`) runs on every PR; `web.yml` runs when `web/**`,
  `LOSAT/src/**`, `LOSAT/Cargo.*`, or the WASI test files changed. Expect 5-20 minutes. A PR into
  `dev` needs these (`rust`, `fast output regressions`, `web result`) and merges by auto-merge.
- `platform-checks.yml` (`platform checks result`) runs its platform jobs only on a PR into
  `main` (the promotion PR), when `LOSAT/src/**`, `LOSAT/Cargo.*`, `LOSAT/build.rs` or its own
  files changed and the PR is not a draft: outputs on Linux x64, Linux arm64, Windows x64 and
  macOS arm64 and `cargo test` on Linux arm64 and Windows. Expect 1-2 hours, most of it the one
  macOS job waiting for a runner (GitHub runs five macOS jobs at a time for the whole account; the
  full five-platform set runs on `main` and weekly). On a PR into `dev` it passes without jobs.
- A PR into `dev` that another merge made `DIRTY` merges `origin/dev`, runs its gates again and
  pushes (no force); its checks take minutes, so several PRs of one program may be ready at once.
  Do not merge `origin/dev` into a PR that has no conflict.
- The promotion PR: do not push to it unless it is `DIRTY` or a check failed. A push cancels its
  Platform checks run and queues it again behind every other PR (2026-10-16: eight cancelled runs
  cost 41 hours of queue). Fix a failure through a work PR into `dev`, then merge `origin/dev`
  into the promotion branch once.
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

- A PR into `dev` merges by auto-merge when its checks pass (`auto-merge.yml` turns it on when the
  PR is opened as ready or marked ready). A session merges by hand only when its instruction says
  so, with CI green: `gh pr merge <n> --merge --match-head-commit <full sha>`.
- The promotion PR (`promote/<date>` to `main`) is opened only on the Owner's request; it merges
  by auto-merge unless the Owner turns that off.
- A PR that changes agent instructions, skills, settings, product decisions (`docs/product_decisions/`),
  or `AGENTS.md` waits for the Owner's merge.

## After the merge

Record the merge commit in the gate record and `STATE.md` (skill `losat-campaign`); clean up
task worktrees (skill `losat-worktree`).
