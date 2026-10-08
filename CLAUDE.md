# CLAUDE.md

@AGENTS.md

`AGENTS.md`, imported above, is the authoritative guidance for LOSAT, including current parity
status, accepted exceptions, build commands, and verification rules. If this file conflicts with
`AGENTS.md`, follow `AGENTS.md`. Code under `web/` also follows `web/AGENTS.md`.

Use repository plans, release documents, session instructions, and gate records for
task-specific state. Do not add session-resume notes, old investigations, or recently fixed lists
to this file.

## Claude Code tooling

Machine paths (`$WORK_ROOT`, `$TASK_DIR`, `$BUILD_ROOT`, `$NCBI_SRC`, `$NCBI_BIN`,
`$ORACLE_JOBS`) come from the untracked `CLAUDE.local.md` at the clone root.

- Skills (`.claude/skills/`): `losat-worktree` (worktrees, branches, builds, commits),
  `losat-gates` (verification tiers quick, standard, full), `losat-oracle-runs` (NCBI BLAST+ as
  an oracle: concurrency cap, lock, work dirs), `losat-campaign` (state files and session splits
  for multi-session stages), `losat-ship-pr` (pull requests and CI), and
  `verify-ncbi-parity-and-speed` (parity and performance method).
- Agents (`.claude/agents/`): `losat-implementer` (one step of a port or fix),
  `losat-mechanic` (merges, records, CI fixes with a clear log), `losat-test-runner` (runs named
  checks, reports only failures), `losat-inventory` (one range of an NCBI call-path inventory),
  and `losat-reviewer` (plan review, code review, independent parity audit, screen review;
  read-only).
- The main session orchestrates. It hands reading and porting to `losat-implementer` and long
  runs to `losat-test-runner`, and keeps logs in files rather than in the conversation.

# Compact instructions

When compacting, keep: the stage and its instruction file path; the worktree, branch, head SHA,
and whether it is pushed; running background jobs with their PIDs and output paths; the task
folder and its `STATE.md`; open Owner decisions and DW numbers; the next step.
