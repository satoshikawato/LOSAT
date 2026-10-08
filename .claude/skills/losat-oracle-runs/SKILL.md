---
name: losat-oracle-runs
description: Run NCBI BLAST+ and LOSAT searches for comparison (oracle runs, sweeps, captures, gate scripts, audits that execute searches) without overloading the machine - the concurrency cap, the machine-wide lock, work directories, inputs that reach the network, Gate A's /tmp fixture, process handling, and how to wait for long runs. Use before starting any command or agent that runs many searches.
---

# Running the NCBI oracle and LOSAT searches

NCBI BLAST+ is a comparison oracle only (`AGENTS.md` rule 2). Paths come from
`CLAUDE.local.md`: `$NCBI_BIN` (BLAST+ 2.17.0, the build the gate records cite), `$NCBI_SRC`
(source at the pinned commit), `$BUILD_ROOT`, `$ORACLE_JOBS`.

## Why there is a cap

The machine has 32 logical CPUs and 31 GB of memory.

- 2026-10-05: the TBLASTX option sweep, V-ABI full, and Gate A overlapped; one
  `-threshold +inf` case peaks at about 8.6 GB. WSL restarted.
- 2026-10-06: five inventory agents ran oracles at once. The load average reached about 21, the
  `/mnt/c` mount failed, and two sessions lost their uncommitted work.

The Owner had allowed many parallel Sonnet agents (2026-10-02). That still holds for agents that
only read; agents that run searches follow the cap below (Owner approval, 2026-10-08).

## The cap

- At most `$ORACLE_JOBS` (8) NCBI or LOSAT search processes at once, machine-wide, and a load
  average under 12. Check `uptime` before starting a pool.
- Always pass `--jobs` explicitly. Several scripts default to `os.cpu_count()` (32):
  `ci_fast_regressions.py`, `range_sweep.py`, `range_regression_fixtures.py`,
  `tblastx_regression_fixtures.py`, `ctoolkit_compare.py`; `run_v_abi_parallel.py` defaults to 12
  and `option_sweep.py` to 8. Use `--jobs 4` for V-ABI full and capture, `--jobs 3` for sweeps,
  and `--jobs 2` for TBLASTN/TBLASTX option sweeps.
- Run memory-heavy sweeps one at a time. Run the stages of a full gate in sequence (quick checks,
  then capture, then V-ABI full alone), not as parallel background streams.
- At most two agents that run searches at the same time (audit angles, inventory ranges). Other
  angles wait or work from source only.
- Take the machine-wide lock for any heavy run (capture, sweeps, V-ABI full, Gate A, audits that
  run searches), so the engine and app sessions do not overlap:

```bash
flock "$BUILD_ROOT/oracle.lock" python3 docs/evidence/losat_web_e1a/capture_outputs.py run ... --jobs 4
```

- V-PERF needs a quiet machine (load under 1) and the existing lock script
  (`~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`; the app track waits while
  `~/.cache/losat-web-gui-target/vperf.lock` exists, README rule 5).

## Work directories and inputs

- Each run gets its own directory under `$BUILD_ROOT/<task>/` (or `$TASK_DIR`), and each agent
  gets its own output path and target directory. Append results as they arrive, so that a stop
  or a usage limit loses nothing.
- Read inputs from the worktree on the Linux clone, never from `/mnt/c`.
- Gate A reads its lexical fixtures from `/tmp/losat-pr5-runtime-cert-*`, created by
  `stage_lexical_fixtures()` in `LOSAT/tests/ci_fast_regressions.py`. Never delete it during a
  gate, and recreate it after a WSL restart. When cleaning `/tmp`, first check that no process
  uses an entry (`/proc/*/cwd`, `/proc/*/fd`) and keep system and tool entries.
- `~/.ncbirc` must not exist (parity runs assume it is absent).
- FASTA lines that NCBI reads as a Seq-id make NCBI contact the network. Keep such cases few.

## Processes and waiting

- Start long runs in the background with output to a log, keep the PID
  (`cmd > "$LOG" 2>&1 & echo $! > "$TASK_DIR/<name>.pid"`), and stop them with
  `kill "$(cat <pid file>)"` or `pkill -P <pid>`. `pkill -f` and `pgrep -f` are denied by a hook:
  they matched their own shell in four sessions.
- Wait with one watcher: a background `until ! kill -0 "$PID"; do sleep 60; done` or Monitor. For
  runs longer than half an hour, check back at most every 30 minutes (the Owner's prompt-cache
  concern, 2026-10-02). Do not busy-poll.
- If `/mnt/c` returns `Input/output error`, it does not recover without a WSL restart: stop
  agents and polls, save anything held only in `/tmp` or in context to the task folder, write a
  resume note, and tell the Owner it is safe to restart WSL.
