---
name: losat-worktree
description: Set up, use, and clean up LOSAT git worktrees on the Linux clone - which branch each track uses, where builds and temporary files go, the Rust toolchain and target directories, committing and pushing each verified step, the staged-size check, and cleanup after a merge. Use when starting or resuming work that edits the repository, and when finishing it.
---

# LOSAT worktree lifecycle

This skill covers mechanics only; the session instruction decides what the task includes.
The machine paths come from `CLAUDE.local.md` at the clone root: `$WORK_ROOT` (the clone),
`$TASK_DIR` (the task folder outside the repository), `$BUILD_ROOT` (build and run output).

## Where to work

- Never work in a checkout on `/mnt/c` (the Windows drive). Its 9p mount failed eight times
  between 2026-09-28 and 2026-10-06 and stranded uncommitted work. The Owner's
  `/mnt/c/Users/genom/GitHub/LOSAT` checkout is read only.
- The clone root stays a detached `origin/main`. Do not edit or commit there.
- Long-lived worktrees of the LOSAT Web campaign (`docs/losat_web_gui_sessions/README.md`,
  rule 1):

| Track | Worktree | Branch |
| --- | --- | --- |
| Engine (`LOSAT/`, `web/adapter/`) | `$WORK_ROOT/.worktrees/web-gui` | `feature/losat-web-gui` |
| App (`web/app/` only) | `$WORK_ROOT/.worktrees/web-gui-app` | `feature/losat-web-gui-app` |

- Any other task: one worktree per task.

```bash
git -C "$WORK_ROOT" fetch origin
git -C "$WORK_ROOT" worktree add --detach ".worktrees/<topic>" origin/main
cd "$WORK_ROOT/.worktrees/<topic>"
git switch --no-track -c <prefix>/<topic> origin/main   # fix, feature, docs, audit, ci
```

- At the start of a session: `git status` must be clean, then `git pull --ff-only`. Merge
  `origin/main` only when the task needs a change from it, and resolve conflicts in their own
  commit (README rule 2).

## Builds and temporary files

- Rust: `RUSTUP_TOOLCHAIN=1.92.0` (or `cargo +1.92.0`) for every cargo or rustc call, including
  scripts that call cargo themselves (for example `LOSAT/tests/build_wasi_artifacts.py`).
- Build output goes outside the worktree and is reused: `--target-dir "$BUILD_ROOT/<purpose>"`
  with the names `native`, `test`, `wasi`, `adapter`, and `gate-<stage>` for a gate's frozen
  build. Do not create a new target directory per attempt; each one starts cold and costs
  gigabytes (168 such directories filled 197 GiB by 2026-10-06).
- `cargo test --all-features` needs `LOSAT_BLASTX_WORKER_LOG` set to a writable file, as in CI.
- Logs, probes, sweeps, and captures go to `$TASK_DIR` or `$BUILD_ROOT/<task>/`, never into the
  worktree. Set `TMPDIR` to the task folder's `tmp/`; `/tmp` is wiped by a WSL restart.
- Do not run `git stash`. Every worktree of the clone shares one stash list.
- Track background processes by PID (`$!`, a PID file); never `pkill -f` or `pgrep -f` (a hook
  denies them).

## Commit and push

- Commit each coherent, verified change and push it at once: `git push -u origin <branch>`.
  A push to a work branch is the off-machine backup. Keep application commits apart from engine
  commits (`web/AGENTS.md` rule 9). For app work the Owner asked for frequent commits and pushes
  ("随時どんどんコミットプッシュ", 2026-10-06).
- Stage named paths only; never `git add -A` on an evidence tree.
- Before every commit, list the staged sizes and stop and ask the Owner if any file is over
  about 1 MB or the commit is unexpectedly large ("くれぐれもこういうクソでかいファイルはコミットとか
  プッシュしないでね", 2026-10-04):

```bash
git diff --cached --name-only -z | xargs -0 -I{} sh -c 'printf "%s\t%s\n" "$(git cat-file -s ":{}" 2>/dev/null || echo deleted)" "{}"' | sort -rn | head
```

- Bulky run output stays under `$BUILD_ROOT`; `docs/evidence/` gets summaries, logs, scripts,
  and hashes only.

## Finish

- After the PR merges (or the branch is otherwise done), run `git status` in the worktree, then
  `git worktree remove <path>` and `git branch -D <branch>`. Keep the two campaign worktrees.
  Remove only worktrees you created; ask before touching another session's.
- Delete the task's `$BUILD_ROOT/<task>/` scratch when the evidence it backs is committed.
  Keep target directories named in an open gate record until that stage is merged.
