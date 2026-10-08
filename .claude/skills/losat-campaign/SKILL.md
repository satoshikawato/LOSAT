---
name: losat-campaign
description: Keep the state of a multi-session LOSAT stage in files - which tracked documents are the authority, the volatile STATE.md in the task folder that is reloaded after a compaction, when to end a session, how to hand off to the next session's instruction file, and when to split work into a parallel session. Use when starting or resuming a session from docs/losat_web_gui_sessions/ (or any work that spans sessions), before compacting or ending such a session, and when asked how far the work has progressed.
---

# LOSAT campaign state

## Where the state lives

| What | Where | Rule |
| --- | --- | --- |
| Plan, completion conditions, maintainer decisions (DW-n) | `docs/losat_web_gui_plan.md` | Authority; change only as the instruction says |
| Session order and status | `docs/losat_web_gui_sessions/README.md` table | Updated by the engine track at the end of a session |
| This session's task | its instruction file in `docs/losat_web_gui_sessions/` | Pasted by the Owner; do not loosen its conditions |
| Results, commits, residuals | gate record `docs/evidence/losat_web_<stage>/README.md` | Written at the end, with run directories that are never rewritten |
| Volatile state (worktree, head SHA, pushed or not, running jobs and PIDs, logs, next step, open questions) | `$TASK_DIR/STATE.md` | Overwrite; at most 150 lines; read first after a resume or compaction |
| Owner-delegated choices taken with the recommended option | `$TASK_DIR/DECISIONS.md`, then the gate record | Append; cite them in the final report |

`$TASK_DIR` is `/home/kawato/losat-baselines/<stage>-<yyyymmdd>/` (see `CLAUDE.local.md`). Do not
keep resume notes inside build directories.

## Start or resume

1. Create `$TASK_DIR` (with `tmp/` and `logs/`) from `templates/STATE.md` if it is new.
2. Bind the session so that the SessionStart hook prints `STATE.md` after a compaction or resume:
   `python3 .claude/skills/losat-campaign/scripts/bind_session.py "$TASK_DIR"`.
3. Read `STATE.md`, the instruction file, and the previous gate record. Check README rule 2
   (branch, clean status, `git pull --ff-only`).

## While working

- Update `STATE.md` in the same turn as a commit, a push, a started or finished background job,
  an abandoned approach, or a decision.
- Answer "how far along is it?" from the instruction's numbered steps and `STATE.md`, not from
  memory.

## End the session

End it, and continue in a new session, when any of these holds:

- a PR merged and the next work touches other files;
- the context is past about 300k tokens, or a second compaction would be needed;
- the instruction's completion conditions are met;
- the session sat idle for more than an hour (the prompt cache is gone).

Before ending (README rule 8): commit and push; write the gate record; update the next session's
instruction file from the measurements and results, and the README table; commit and push those.
The final answer names the commit SHAs, the push result, the verification results, the
residuals, the next instruction's path, and its full text. If a completion condition remains,
do not call the stage complete; make it the first step of the next instruction.

## Parallel sessions

At most one engine-track and one app-track session at a time (plan DW-7). When a session from the
table can run alongside the current one (its own worktree, entry condition met or waived), do not
do it in the current session: update its instruction file, commit and push, and give the Owner its
full prompt (Owner, 2026-10-03). Heavy oracle runs are serialised across both sessions by the lock
in skill `losat-oracle-runs`.
