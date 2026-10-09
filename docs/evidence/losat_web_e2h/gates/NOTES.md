# SFb (E2h) gate scripts: notes

**Status (2026-10-09):** drafts written in SFb and checked with `bash -n` / `py_compile` only. SFc changed `sf_gates.sh` before the first run: the Stage G namespace mounts a tmpfs on `/mnt` and binds the worktree at the path inside it (`STAGE_G_BIND` defaults to 1; `guard` checks the path only inside the namespace, so nothing outside reads `/mnt/c`), and the wasm32 test filter is `web_api::` (S10 moved the `bio` tests to `web_api::v1_bio`). SFc tries `STAGES=fast-all` first.

Drafts written 2026-10-08, not run (checked with `bash -n` and `py_compile` only). They go to
`$WT/docs/evidence/losat_web_e2h/gates/` (all files in this directory; `verify_added.py` is the one with E2h changes).
Paths are relative to the worktree root as in S11. Env defaults inline (`WT`, `BUILD_ROOT`, `NCBI_BIN`, `NCBI_SRC`).

## Run

```bash
source /home/kawato/losat-baselines/sfb-e2h-20261008/env.sh; cd $WT
setsid nohup docs/evidence/losat_web_e2h/gates/sf_gates.sh > $TASK_DIR/logs/gate.log 2>&1 &   # keep the PID; end the turn, wait once
# after a stop:   .../sf_gates_resume.sh      (skips stages with a marker, checks HEAD and artifacts.sha256)
# then alone:     .../sf_gate_a.sh            (3-5.6 h)      then:   .../sf_perf.sh   (quiet machine)
# then:           FORCE=1 STAGES=collect .../sf_gates_resume.sh      (the run record picks up Gate A and V-PERF rows)
```

State: `$BUILD_ROOT/sfb-e2h/gate-sf.run` (run id), run output `$BUILD_ROOT/sfb-e2h/gate-<UTC>/` (`status.tsv`, `stage-times.tsv`,
`.done/<stage>`), tracked record `docs/evidence/losat_web_e2h/run-<UTC>/` (written by stage `collect`: `status.tsv`, `status-failed.txt`,
small logs and tsv, `failed/` full logs of failures, hashes, capture `hashes.tsv`, reactors json; files over 300 KB get a sha256 and a tail or gzip).
`GATE=<name>` changes the build prefix (`$BUILD_ROOT/<name>/...`, default `gate-sf`). `ALLOW_DIRTY=1` builds a dirty tree (the record names it).
Every stage restages `/tmp/losat-pr5-runtime-cert-*` and checks HEAD against `head.txt`. Gates that fail are recorded, not fatal; a
build failure, a dirty tree or a missing frozen input stops the run. Heavy stages run as `flock "$BUILD_ROOT/oracle.lock" sf_gates.sh --stage X`
(all but `prep`, `lint`, `collect`), one at a time, `--jobs` explicit.

## Stages in order (expected time; source: skill `losat-gates` unless marked)

| # | Stage | What | Compared with | Time |
| --- | --- | --- | --- | --- |
| 1 | prep | head, tree clean check | | <1 min |
| 2 | lint | verify_refs + `verify_added.py` (base f3048ffde), protein tables `--check`, fmt (LOSAT, adapter), pure-Rust boundary + 2 unittests | | ~2 min |
| 3 | clippy | 4 LOSAT configs + 3 adapter configs `-D warnings` | | ~5 min warm; cold unmeasured (10-20) |
| 4 | tests | `cargo test --all-features` (includes the reader unit tests), reader filter `blastinput::fasta_reader` alone, adapter `cargo test` + `--test scan_properties` + `--test scan_ncbi_properties` (99 if the file is missing), wasm32 `web_api::tests` | | 5-10 min |
| 5 | build | native, native-serial, wasi (+serial), reactors, build identity, `artifacts.sha256` (one frozen build) | | cold unmeasured (10-20) |
| 6 | quick-fixtures | `check_losat.py` threads 1/2/4, `run_oracle.py` | frozen fixtures | ~2.5 min |
| 7 | regression-fixtures | range 35, TBLASTX 84, BLASTN 183 (`--jobs 4`) | frozen NCBI hashes | ~3 min |
| 8 | sf-fixtures | `LOSAT/tests/fasta_input_fixtures.py check --jobs 3` (all programs; expected `$BUILD_ROOT/sf-e2h/fixtures/ncbi`) | NCBI 2.17.0 frozen (3932 rows; before: same 1131, differs 13, rejects 2764, pending 24; after: all `same`/listed rejections) | unmeasured (est. 5-15 min) |
| 9 | sf-sweeps | `fasta_sweep.py check` (2280) and `check_inputs.py check` (1032), `--jobs 3`, run from `$BUILD_ROOT/sfb-e2h/sweeps`; generate/freeze only if missing; `summary` | frozen NCBI in `sweeps/ncbi`; `unlisted`, `differs`, `timeout` must be 0 | check seconds; freeze 100 s + 43 s |
| 10 | input-sweeps | `ctoolkit_compare`, `punct_defline`, `title_sweep`, `protein_title_sweep`, `range_sweep` (`--jobs 3`) | NCBI live | ~5 min |
| 11 | fast-all | `ci_fast_regressions.py --all-cases --jobs 4` (S02 baseline, Gate A + Stage G frozen hashes, fixtures of every program) | frozen hashes (1 known mismatch `Sakai.MG1655.megablast`) | 30-40 min (est.; same 236 cases as capture) |
| 12 | v1-wasi | `check_wasm_threading.py` + `v1_requests.js`, `cmp` with S11's `v1-requests-after.jsonl` | S11 run (ABI v1 frozen) | 2-3 min |
| 13 | vabi-quick | `v_abi_cases.py --suite quick` + `v_abi.js` | | 4-10 min |
| 14 | option-blastp | `option_sweep.py` (`--jobs 3`, timeout 1200) + class comparison | E2e run `sweeps/after-blastp.tsv` (1196 rows); changed rows listed in `sweeps/compare-blastp.txt` | unmeasured (est. 10-20 min) |
| 15 | option-tblastn | same, `--jobs 2` | E2e `after-tblastn.tsv` (1427) | 30-34 min |
| 16 | option-tblastx | same, `--jobs 2`, + `gencode_api_check.py` (S11 skipped this sweep) | E2e `after-tblastx.tsv` (641); up to 8.6 GB per case, so alone | 30-34 min |
| 17 | capture | `capture_outputs.py run --jobs 4` (236 cases), `compare` | S02 `docs/evidence/losat_web_e1a/baseline/hashes.tsv` and `~/.cache/losat-web-gui-target/sf/capture-before/hashes.tsv` (S11 build) | median 37 min |
| 18-21 | vabi-blastn, -blastp, -tblastn, -tblastx | `run_v_abi_parallel.py --jobs 4` on the full suite split by program (each alone) | native vs WASI serial/threaded inside the tool; split count must add up (`v-abi-full-cases-split`) | median 167 min together |
| 22 | collect | tracked run record | | <1 min |
| | `sf_gate_a.sh` | Gate A, 20 pairs, alone, under flock; restages the `/tmp` fixture | v0.1.0 outfmt 6 parity (NCBI oracle) | 3-5.6 h |
| | `sf_perf.sh` | V-PERF, before = `sf/bin/` (S11), after = gate build, alternating; 11 standard + 4 read-heavy cases; over +5% re-measured `--repeat 5` | before binaries | 2-3 min standard; read-heavy unmeasured |

Whole `sf_gates.sh` about 6-7 h (V-ABI 2.8 h, option sweeps 1.2 h, capture + fast-all 1.2 h, builds/lint/tests 0.5-1 h), plus Gate A and V-PERF: 10-13 h wall.

## What changed from S11 (`s11_gates*.sh`, `s11_gate_a.sh`, `s11_perf.sh`)

- Added: SF fixtures (stage 8), the two SF sweeps (9), reader unit tests and the two adapter property tests named in the log (4), `GATE`/run-id pointer,
  `status.tsv`, `collect`, markers per stage, per-program V-ABI stages, `guard` (restage fixture, HEAD check), `ALLOW_DIRTY`, load wait (`need_load`).
- Run now (S11 skipped): TBLASTX option sweep, capture 236 with both comparisons, Gate A 20 pairs, V-ABI full for all 4 programs.
- Dropped or replaced: V-ABI full as a background job beside the fixtures (now alone, sequential); `--jobs 6` everywhere (now 4 / 3 / 2); the E2g
  `check_inputs.py` run (E2h's `check_inputs.py` contains E2g's 300 BLASTN rows, frozen NCBI, plus 244 each for TBLASTX, TBLASTN, BLASTP);
  S11's resume list (TBLASTX `e2d.` filter) is replaced by marker-based resume; E2d `rundir` file; `/mnt/c` worktree paths -> `$WT`; `$S/ncbi/c++` -> `$NCBI_SRC/c++`
  (same pinned commit, same line numbers); `/tmp/claude-1000` Gate A output -> `$BUILD_ROOT/sfb-e2h/gate-<UTC>/`; builds `s11-gate-*` -> `$BUILD_ROOT/gate-sf/*`;
  the capture-vs-S11-gate comparison now uses `sf/capture-before` (the S11 build captured at SF's start). Sweep results of S11 for BLASTP/TBLASTN are not used; E2e's classes are the reference for all three.
- Kept from `~/.cache/losat-web-gui-target/`: `s07p-resume/vperf_lock.sh` (its `APP_WORKTREE` still names the old `/mnt/c` path, so `sf_perf.sh` resets it to
  `$WORK_ROOT/.worktrees/web-gui-app` after sourcing), `s08p/verify_refs.py` (as S11; the `s07p-resume/` copy hard-codes `/mnt/c/.../ncbi-blast`), `s08p/api-oracle/` (TBLASTN gencode oracle), `sf/bin/`, `sf/capture-before/`.
- `verify_added.py`: base is an argument (default `f3048ffde`, was `4fab73fdb`).
- V-PERF: new `perf_cases.py` (wraps E2b's) with `blastn-q100k`, `blastn-genome-1line`, `blastn-genome-80col`, `blastp-many`; inputs from `gen_perf_inputs.py`
  (seeded, `$BUILD_ROOT/sfb-e2h/perf-inputs/`, not committed); standard list adds `blastn`, `blastn-large-fmt0`, `blastn-many`; automatic `--repeat 5` for cases over +5%.

## Open decisions (recommendation first)

1. TBLASTN Stage G cases name `/mnt/c/Users/genom/GitHub/LOSAT/...` (capture, fast-all, V-ABI; `guard` checks the directory exists). Reading the Owner's Windows
   checkout over 9p breaks the "never read /mnt/c" rule and risks EIO. SFc (2026-10-09): `STAGE_G_BIND=1` is the default: those stages run under `unshare -rm` with a tmpfs on `/mnt`,
   the path created in it and `$WT` bound there (private namespace, no sudo, nothing read from /mnt/c; probe 2026-10-09 listed and wrote through the bind).
2. Capture (37 min) repeats the 236 cases that fast-all already hashes. Recommend keeping both (instruction names capture; cheap against the rest).
3. `v_abi_cases.py --suite full` has no SF fasta-input cases, so V-ABI full does not exercise the new reader's rows through the ABI beyond existing cases. Recommend the implementer add a `fasta_input` group before the gate (or the verification cell cites the scan property tests only).
4. Read-heavy perf sizes (100k x 300 nt, 5 Mb, 20k proteins) are not timed; the 100k query in WASI modes may run minutes per sample. Recommend one native timing first, then `gen_perf_inputs.py --queries/--proteins` smaller or `PERF_READ_MODES=native,serial-wasi`.
5. The tracked run record keeps small files only; sf-fixtures/sweep TSVs over 300 KB are gzipped (or sha + tail). Recommend this.
