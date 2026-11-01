# Threaded Wasm in browsers: the shared-memory guard at build time (2026-11-01)

Release v0.2.0 GO item (`release-v020` DECISIONS 52). gbdraw's browser runs of LOSAT's threaded
command trapped with `RuntimeError: memory access out of bounds` in a `memory.copy` or
`memory.fill` of one thread after another thread had grown the shared memory: V8 checks bulk
memory instructions against bounds that can be stale in that isolate. The Node runner has
guarded the module at load time since S09 (`LOSAT/tests/wasi_shared_memory.js`), and so has
LOSAT Web's Vite build (`web/app/build/reactors.ts`); the threaded artifacts that CI and the
release publish (`losat-threaded-command.wasm`, `losat-threaded-reactor.wasm`, the
wasm32-wasip1-threads tarball) were not guarded, and a browser host that loads them as they
are has no guard.

This branch guards those artifacts at build time with the same function, marks a guarded
module with a custom section `losat.shared_memory_guard` (`{version, fill, copy}`), and lets
every host run a marked module as it is.

## Step 1: traps with and without the guard, in one browser

### gbdraw-type harness, threaded command (`scripts/serve.py`, `scripts/run.py`, `scripts/rounds.sh`)

gbdraw 900992caf's browser host (`web/js`, `wasm/losat`, `vendor/browser_wasi_shim`, exported
read-only with `git archive`) and its `harness.html`, served on 127.0.0.1 with COOP/COEP and
driven by Python Playwright 1.64.1 with Chromium 156.0.8078.4 (headless shell 1248). One
"pass" runs every job of a job list, 4 at a time, each with 4 threads, as gbdraw does (gbdraw
OV-2369). The arms of a round run in one browser process, in alternating order.

| Arm | Module | Passes | Jobs | Jobs that trapped |
| --- | --- | ---: | ---: | ---: |
| `devraw` | `dev` b94e9844f + this branch, threaded command as cargo writes it (sha256 `cc6ac44e5109fc90...`) | 13 | 164 | 5 |
| `devguarded` | the same, guarded by `build_wasi_artifacts.py` (`5771656804c45c2b...`, fill 1312, copy 6734) | 13 | 164 | 0 |
| `gshipped` | the module gbdraw ships at 900992caf (LOSAT build of 2026-05-25, no guard) | 7 | 56 | 6 |

- `devraw` against `devguarded`: one-sided Fisher exact p = 0.030; at the `devraw` rate
  (5/164), 0 traps in 164 jobs has probability 0.006.
- Job lists: `jobs8` (8 self-comparisons, rounds a and b) and `jobs14` (gbdraw's 14 jobs,
  round c, 10 passes), both from gbdraw's harness data.
- Every trap is a bulk memory instruction (`results/command_harness_traps.tsv`, offsets in the
  module V8 compiled): `devraw` 4 `memory.fill` (in `calloc`, `BlastIntervalTree::add_hsp`,
  `BlastpSubjectScratch::new`, `x_hits_flat::x_lists_by_query`) and 1 `memory.copy`
  (`encode_query_frame`), 66-359 ms after the WASI thread started; `gshipped` all in the
  `memory.copy` of `blastp::blast_engine::fasta_id` (the parallel FASTA id parse).
- Outputs: one SHA-256 per job across all passes and both `dev` arms
  (`results/command_harness_outputs.tsv`). gbdraw's host runs a trapped job again (up to 2 retries),
  so the output of a trapped job comes from its retry.
- The `dev` command declares 43 initial pages (gbdraw's host passes 21) and takes NCBI-style
  options (`-query`, `-num_threads`; gbdraw's host writes `--query`, `--num-threads`): the
  harness takes the initial size from the module and rewrites the options for the `dev` arms
  only.
- System load 16-96 when a pass started (other sessions); per-pass rows, with that load, in
  `results/command_harness_passes.tsv`.

### LOSAT Web, reactor ABI v2 (`web/app/tests/e2e/guard-measure.spec.ts`)

The opt-in measurement spec (`LOSAT_WEB_MEASURE=guard`, wrapper
`scripts/losat_web_guard_measure.sh`) builds the contract harness twice in one browser: once
with the threaded module as the Vite build serves it (guarded) and once with the reactor
artifact as it is. Each round runs the same BLASTP searches in both arms (ABBA), each search in
a new runtime so that the shared memory grows from its initial size. Reactors: `dev` 5ce40a9e0
(2026-10-30 build).

| Run | Chromium | Searches at a time x threads | Rounds | Guarded stopped | Unguarded stopped |
| --- | --- | --- | ---: | ---: | ---: |
| r1-n8 | 149.0.7827.55 (Playwright 1.61.1) | 1 x 8 (AP027078+AP027131 vs AP027132+AP027133) | 10 | 0/10 | 0/10 |
| r2-c4n4 | 149.0.7827.55 | 4 x 4 (gbdraw jobs 1, 3, 5, 7) | 4 | 0/16 | 0/16 |
| r3-c4n4-c156 | 156.0.8078.4 | 4 x 4 (gbdraw jobs 8-11, 1.7-2.0 MB pairs) | 10 | 0/40 | 0/40 |

The reactor path never trapped without the guard (0/66). That does not show it cannot: at a
3 % rate, 0 in 66 has probability 0.13. The shipped LOSAT Web serves the guarded module, and
`web/app/tests/e2e/engine-guard.spec.ts` now checks that in CI. Concurrent searches first
failed with `TypeError: Failed to fetch` (Chromium `net::ERR_CACHE_WRITE_FAILURE` when two
workers write the same module into the HTTP cache); the harness server sends
`Cache-Control: no-store` when several searches run at a time.

## Time cost of the guard (Node, in-process)

One Node 26.8.2 process runs both arms of a case in turn (`scripts/wasm_guard_ab.js`,
`scripts/run_ab.sh`), through the host's `createThreadHost`: the 9 daily cases at 8 threads, 5
rounds after one discarded search per arm, the order alternating by round. Arms: the guarded
command against the same module as cargo writes it. The null arm is the guarded command against
itself with one more custom section (other bytes, so V8 shares no compiled code; same code).
No `LOSAT_*` variables except the guard cache switched off. Local machine (WSL2, 32 logical
CPUs), load 2.9-11.1.

| Case | Median guarded / raw, s | Guarded / raw, median of rounds | Rounds | Null, median of rounds | Rounds |
| --- | --- | ---: | ---: | ---: | ---: |
| blastp | 7.84 / 7.74 | x1.014 | x1.009-1.038 | x0.999 | x0.944-1.033 |
| blastp-fast | 10.28 / 10.66 | x1.019 | x0.949-1.044 | x0.967 | x0.922-1.040 |
| megablast | 7.87 / 7.64 | x1.019 | x0.947-1.051 | x1.008 | x0.975-1.129 |
| blastn | 9.83 / 10.15 | x0.985 | x0.927-1.067 | x0.966 | x0.878-0.993 |
| dc-megablast | 2.94 / 2.80 | x1.089 | x0.932-1.108 | x1.088 | x0.950-1.131 |
| blastn-short | 3.15 / 3.30 | x0.960 | x0.928-1.024 | x0.973 | x0.791-1.049 |
| tblastn -db | 2.69 / 2.55 | x1.055 | x0.965-1.137 | x1.019 | x0.967-1.098 |
| blastx -db | 1.51 / 1.47 | x1.041 | x0.962-1.114 | x1.026 | x0.989-1.140 |
| tblastx | 9.58 / 9.50 | x0.993 | x0.823-1.009 | x1.017 | x0.914-1.121 |
| geometric mean of the medians | | x1.019 | | x1.006 | |

The cost is about 1 % over the 9 cases and this run does not separate it from noise: every
case's median lies inside its null arm's range, and dc-megablast's x1.089 has a null of x1.088.
blastp is the one case with all 5 rounds above 1 (x1.009-1.038), so a small cost there is
possible. Each case has one output hash across both runs and all arms (`results/node_ab_summary.tsv`).

This is the cost for a host that had no guard, such as a browser host. The Node runner has run
guarded code since S09: its module is the same as before plus the 64-byte section, and it now
skips the guard and its cache at load.

### One trap of the raw arm in Node

The first attempt of this A/B ran at load 8.6-69.7 and is not used for timing
(`results/node_ab_first_run_trap.*`). In its third round the raw arm of blastn trapped with
`memory access out of bounds` at 0x57b240, a `memory.fill` (`fc 0b 00`) in
`blastn::alignment::gapped::x_prepare_nt_tables`, in a rayon worker: the same fault as in
browsers, and the reason the Node runner guards. Counting every search of all the runs, with the
discarded ones: unguarded 1 trap in 82, guarded 0 in 190 (null arms included). The rerun, at low
load, had none in either arm. The process then spun for hours, because this script caught
SIGTERM, which is how a worker ends the process on a trap (`failWorker`); the main thread stays
blocked in Wasm. The shipped runner does not catch it. The script now lets the process end,
records the search as failed and goes on with the next round in a new process.


## What the guard does not cover

- It wraps `memory.fill` and `memory.copy` only. Each threaded module keeps 4 `memory.init`
  instructions (`__wasm_init_memory` 3, `__wasm_init_tls` 1), and none of the 12 traps above was
  one. LOSAT's hosts refresh the bounds before a thread starts (`wasi_thread_worker.js`,
  `web/app/src/infra/reactor/instance.ts`); gbdraw's host at 900992caf does not.
- Atomic accesses are not guarded. Line WO (`release-v020`, 2026-11-01) reports a trap of the same
  kind at Liftoff's atomics in V8 and plans to widen the guard on top of this change.

## Not checked

- Firefox and WebKit (the traps are a V8 behaviour; the guard runs there, but no trap counts).
- Node 24.21.0 (the release toolchain): local runs used Node 26.8.2; the release assembly was
  run locally with the toolchain capture mocked and `--allow-unfilled-lineage` (guarded: `WASM_ARTIFACT_READY`, smoke 5/5;
  unguarded: refused with "no verified shared-memory guard"; `verify_release_archive_v020.py`
  passes the guarded archive).
- `release-readiness-v020.yml` on GitHub. `wasm-threading.yml` passed on the branch at aaf15116a
  (run 41303426863); both threaded artifacts of that run verify (`losat-threaded-command.wasm`
  `eb75215caf76f336...`, `losat-threaded-reactor.wasm` `5d113fa43ffd2d5f...`).
- LOSAT Web reactors built from this branch's head (the step-1 runs used `dev` 5ce40a9e0;
  the guard is applied by the Vite build in both).
- gbdraw's own port (`gbdraw/web/js/services/losat-shared-memory-guard.js`) does not know the
  marker and wraps the helpers of a guarded module once more until gbdraw drops it.

## Files

- `results/command_harness_passes.tsv`, `results/command_harness_traps.tsv`: the gbdraw-type
  harness, one row per pass and per trap.
- `results/command_harness_outputs.tsv`: the same harness, the output hash of each job.
- `results/losat_web_rounds.tsv`: LOSAT Web, one row per arm and round.
- `results/node_ab_summary.tsv`, `results/node_ab_guarded_vs_raw.jsonl`, `results/node_ab_null.jsonl`,
  `results/node_ab_tblastn_guarded_vs_raw.jsonl`, `results/node_ab_tblastn_null.jsonl`: the timing
  (tblastn ran separately, its inputs copied from the workstation).
- `results/node_ab_first_run_trap.jsonl`, `results/node_ab_first_run_trap.log`: the first attempt,
  with the trap.
- `scripts/`: the harnesses as run (`serve.py`, `run.py`, `rounds.sh`, `harness.html` from
  gbdraw's evidence, `losat_web_guard_measure.sh`, `command_harness_tables.py`) and the timing (`wasm_guard_ab.js`,
  `run_ab.sh`, `summarize_ab.py`); paths point at the session's scratch directories.
