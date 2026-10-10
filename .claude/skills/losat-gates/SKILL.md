---
name: losat-gates
description: Choose and run LOSAT verification by what a change touches - the quick, standard, and full tiers with their commands, run times, and when each is required for engine, adapter, app, and docs changes; the conditional full-gate items (Gate A, V-ABI full, option sweeps, capture, V-PERF, independent audit); and GitHub CI. Use before committing, before pushing, before a pull request to main, and when writing a session instruction's gate section.
---

# LOSAT verification tiers

Pick the tier from the paths the diff touches, not from the kind of session. A higher tier
includes the lower ones. A session instruction (`docs/losat_web_gui_sessions/*.md`) that
explicitly requires a gate wins over this table; when you write an instruction, choose its gates
from this table (Owner approval, 2026-10-08). Run times were measured on this machine in
2026-09/10 (sources: gate logs, file times, `gh`, transcripts). Every search-running step follows
the cap and lock in skill `losat-oracle-runs`. Hand runs longer than a few minutes to the
`losat-test-runner` agent, or run them through `.claude/scripts/run_quiet.py`, so the full output
stays in a log file.

`AGENTS.md` rule 10 still applies to parity work: do not run "maybe helpful" tests while known
NCBI divergences remain open. App commits follow `web/AGENTS.md`.

Notation: `$N` = the LOSAT release binary under test; all cargo calls use
`RUSTUP_TOOLCHAIN=1.92.0` and `--target-dir "$BUILD_ROOT/<purpose>"`; run from the worktree root.

## Which tier

| Diff touches | Per commit | Before push | Before a PR to `main` or a parity/performance claim |
| --- | --- | --- | --- |
| Docs, evidence records, instructions only | none | none | none (CI `rust` runs; fast regressions select no program) |
| `web/app` only | quick | standard (app) | `run_gate.sh` once at the end of the stage |
| `web/adapter` (ABI, reactors) | quick | standard + V-ABI quick + v1 WASI matrix | + V-ABI full for the touched programs |
| One program's directory under `LOSAT/src/algorithm/` | quick | standard, restricted to that program | full, restricted to that program |
| Shared engine code (`cli`, `report`, `blastinput`, `core`, `stats`, `utils`, readers, option parsing) | quick | standard, all programs | full, all programs |

## Quick (each commit during a step; under 5 minutes, no NCBI runs)

Engine:

```bash
cd LOSAT && cargo fmt --check                                                   # ~12 s
cargo build --release --locked --target-dir "$BUILD_ROOT/native"                # ~1 min incremental
export LOSAT_BLASTX_WORKER_LOG="$BUILD_ROOT/blastx-worker.log"
cargo test --locked --all-features --lib --target-dir "$BUILD_ROOT/test" <module path>   # inline tests
cargo test --locked --all-features --test <binary> --target-dir "$BUILD_ROOT/test" <filter>   # integration
cd .. && python3 docs/evidence/losat_web_e2a/check_losat.py --losat "$N" --threads 1 --programs <touched>   # ~26 s
```

plus the one focused fixture from the stage's `AUTHORITY.md`. Adapter: `cargo fmt --check` and
`cargo test --locked` in `web/adapter` (1-2 min). App: `cd web/app && npm run check` (median
77 s); `npm ci` only when `package-lock.json` changed.

- A bare `cargo test <filter>` builds and links all ten test binaries (the lib tests and the nine under `LOSAT/tests/`) and the bin before it filters.
- `LOSAT/Cargo.toml` sets `[profile.test] opt-level = 1`, as CI does, so `CARGO_PROFILE_TEST_OPT_LEVEL=1` is no longer needed.

## Standard (end of a step and before every push; 10-25 minutes)

Engine, in addition to quick:

- clippy, the four LOSAT configurations of the gate script (`--all-targets --all-features -D warnings`,
  `--all-targets --no-default-features`, `--lib --target wasm32-wasip1 --no-default-features`,
  `--lib --target wasm32-wasip1-threads --features wasm-threads`): about 5 min warm;
- the whole `cargo test --locked --all-features` (920 tests, 1-2 min warm);
- `python3 LOSAT/tests/ci_fast_regressions.py --losat "$N" --out "$BUILD_ROOT/<task>/fast" --changed-from <base> --jobs 4`
  (S02 baseline hashes, Gate A and TLOSAN Stage G frozen hashes, outfmt 0 fixtures, and the
  BLASTN/TBLASTX frozen regressions of the selected programs; 1-2 min);
- `check_losat.py --threads 1`, `2`, and `4` (1.5 min); `python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir "$NCBI_BIN" --out <dir>` when fixtures changed (49 s);
- `python3 LOSAT/tests/{blastn,tblastx,range}_regression_fixtures.py check --losat "$N" --jobs 4 --out <tsv>` for the touched programs (under 1 min).

When `web/adapter`, `web_api`, WASI, or threading changed: build the reactors
(`web/adapter/tools/build_reactors.py`), `check_build_identity.py`, V-ABI quick
(`v_abi_cases.py --suite quick` + `node web/adapter/tests/v_abi.js ...`, 4-10 min, run alone), and
the v1 WASI matrix (`LOSAT/tests/check_wasm_threading.py`, 2-3 min).

App: `cd web/app && npm ci && npm run check && npm run e2e` (FakeEngine build, about 3.5 min).
Add the engine-build E2E (`LOSAT_WEB_REACTORS=... LOSAT_WEB_NATIVE=... npm run e2e`, 2.5-4.4 min)
only when engine-worker, data-worker, or adapter-facing code changed.

## Full (before a PR to `main` or a release-facing claim; one gate at a time)

The stage's gate script (copy of `docs/evidence/losat_web_e2d/gates/s11_gates.sh` with a new
`GATE=` prefix), restricted to the program groups the diff touches, run stage after stage under
`flock "$BUILD_ROOT/oracle.lock"`. Today's whole gate takes 112-273 min (median 193). Add each item
below only when its condition holds:

| Item | Time | Required when |
| --- | --- | --- |
| Capture, all 236 cases (`capture_outputs.py run ... --jobs 4`, then `compare` with the S02 baseline and the previous gate) | median 37 min | TBLASTX or shared code changed; otherwise the standard tier's fast regressions cover the same hashes |
| Gate A, 20 pairs (`LOSAT/tests/audit_tblastx_v010.py`) | 3-5.6 h when overlapped; run alone | TBLASTX or shared scoring, statistics, or reader code changed, unless the capture of the same commit matches the S02 baseline and the Gate A frozen hashes for all 20 TBLASTX pairs (a change proven byte-identical; the NCBI side is the fixed oracle; Owner 2026-10-10) |
| V-ABI full (`run_v_abi_parallel.py ... --jobs 4`) | median 167 min | ABI, adapter, `run_local`, or report layer changed; only the touched programs' cases |
| Option sweeps (`docs/evidence/losat_web_e2e/option_sweep.py`, TBLASTN/TBLASTX `--jobs 2`) | 30-34 min, up to 8.6 GB per case | option parsing, validation, or a program's defaults changed |
| Range and title sweeps (`docs/evidence/losat_web_e2d/range_sweep.py`, the latest stage's `title_sweep.py` and `check_inputs.py` under `docs/evidence/losat_web_e2*/`) | under 1 min each | range handling, FASTA reading, or defline reporting changed |
| V-PERF (`perf_cases.py run ... --repeat 3`; settle cases over +5% with `--repeat 5`, not 10, Owner 2026-10-01) | 2-3 min, quiet machine | hot loops, threading, or output writers changed (plan §6.2 limits the +5% rule to the S02-S04 and SX refactors) |
| Independent audit (`losat-reviewer`, one angle per agent, at most two running searches) | first round about 73 min | new NCBI-ported behavior or a release-facing claim; not for a change proven byte-identical by capture |
| App gate `docs/evidence/losat_web_w3/run_gate.sh` | 19.5 min; `LOSAT_WEB_GATE_STEPS=after-review` 11.8 min after a review fix | end of an app stage |
| Nightly (`gh workflow run nightly.yml`; do not rerun locally) | about 2 h on GitHub | release candidate, or `LOSAT/src` changed since the last green nightly |

## GitHub CI

A PR runs `ci.yml` (`rust` 2.9 min; `fast output regressions` 4.5 min with all programs) and,
when `web/**` or `LOSAT/src/**` changed, `web.yml` (about 20 min; the `adapter` job reruns V-ABI
quick). Nightly runs `--all-cases` and the WASI matrix (median 120 min). Do not rerun locally
what CI runs for the same commit.
