# TBLASTN Stage F native and command-WASI parallel gate — 2026-09-25

Status: **PASS — fixture-scoped TBLASTN Stage F native and command-WASI parity**. The independent [read-only audit](stage_f_independent_audit_20260925.md) supports this scope. Stage F begins at Stage E pass commit `911d91f5cfc31eae3be48cf5ca9efac9970474bd` on `feature/tlosan-tblastn-v0.2.0`. The only behavior authority is fixed NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. The local comparison executable is BLAST+ 2.17.0+ `/home/kawato/micromamba/bin/tblastn`, SHA-256 `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`. The selected-code C++ API oracle has SHA-256 `d8e0f100143031db5121e13907ec34c31a84e258adbf3f92241976916f1536bb` and is comparison-only.

## Program, source, and scheduling boundary

Program/task: TBLASTN local file `-query`/`-subject`, outfmt 0/6/7. The pinned `c++/src/algo/blast/api/prelim_stage.cpp:145-188` starts and joins preliminary search threads. `c++/src/algo/blast/core/blast_engine.c:1410-1476` obtains a subject OID and calls `s_BlastSearchEngineCore`; its frame/chunk search, HSP append, link, E-value calculation and preliminary reap are at lines 804-905. `blast_hspstream.c:289-319` consumes retained query lists in source order. `blast_kappa.c:3383-3427,3525-3736`, `compo_heap.c:252-275,330-391,439-466`, and `blast_hits.c:3243-3297,3383-3437,3420-3437` define redo, containment, heap and final hitlist operations.

LOSAT owner: [stage_d_pipeline.rs](../../../LOSAT/src/algorithm/tblastn/stage_d_pipeline.rs) and [args.rs](../../../LOSAT/src/algorithm/tblastn/args.rs). A search-scoped pool runs the existing six-frame preliminary search, link/direct E-value calculation and preliminary reap within each independent subject job. At most one pool-sized batch of subject results is retained. Rayon's indexed collection restores input OID order before the serial preliminary collector and the unchanged Kappa, containment, heap, hitlist and report paths. One-subject work remains serial. Diagnostics record worker slots and the maximum number of simultaneously active subject calculations. No NCBI binary, library, FFI or subprocess is in LOSAT's runtime, build or fallback path.

First observed Stage F divergence: **none**. The initial 112-subject outfmt-6 check matched NCBI for native threads 1, 2, 4 and 8. No accepted exception was used for code 1. For non-default subject codes, only `PD-TLOSAN-LOCAL-GENCODE-32` applies: the selected code drives local subject translation/search/reporting. All 27 IDs, including 32, use the pinned `FindGeneticCode(code)` C++ API search/formatting oracle; each code-1 API control is calibrated against the local CLI on the same input. The API's optional pairwise query-cover column is removed only under the Stage E calibrated procedure. Code 32 is never passed to the rejecting NCBI CLI.

## Exact output and deterministic repetition

| Scope | Exact result | Evidence |
| --- | ---: | --- |
| 73 physical Stage E commands × 0/6/7; NCBI CLI freshly rerun and checked against its Stage E registered SHA | **1,794/1,794** native output comparisons across 1/2/4/8 and repeated order-sensitive cases | [rows](native_physical/comparison.jsonl), [summary](native_physical/summary.json), [environment](native_physical/environment.json) |
| All 27 Stage E code fixtures × 0/6/7; same-input CLI/API code-1 calibration and selected-code API oracle freshly rerun | **396/396** native comparisons; 81/81 code-1 calibrations | [rows](native_codes/comparison.jsonl), [summary](native_codes/summary.json) |
| All 27 codes, each source subject duplicated into eight independently searched, tied subjects with unique IDs; CLI/API calibration rerun | **810/810** native comparisons; 81/81 code-1 calibrations; all 27 codes exercise the parallel path | [rows](native_codes_multi/comparison.jsonl), [summary](native_codes_multi/summary.json), [generated inputs](native_codes_multi/generated_inputs) |
| Stage E four long generated chunk families and three-query long subject | **27/27** complete serial report bytes against fresh NCBI | [rows](stage_e_generated/comparison.jsonl), [log](stage_e_generated.log) |
| `-out` on 112 subjects × 0/6/7 × 1/2/4/8 | **12/12** file byte matches and empty stdout | [rows](output_file/comparison.jsonl), [log](output_file.log) |
| Stage E unsupported options with thread support updated; zero threads replaces the former two-thread rejection | **14/14** fail explicitly, no stdout or output file | [rows](negative/comparison.jsonl), [log](negative.log) |

The native physical, single-code and eight-subject matrices contain **3,000/3,000** complete output matches. The physical 219 case/format conditions and selected-code 81 conditions retain the Stage E **300/300** share of its 327/327 pass; the fresh generated long reports give the other **27/27**. The long-input runner regenerates its subjects deterministically and verifies their SHA-256 against Stage C manifests. The code-32 eight-subject fixture additionally tests selected translation and tied subject ordering under real parallel scheduling. Every positive row records the exact command, query/subject SHA-256, complete stdout SHA-256, executable SHA-256, thread count, repeat index, exit status and equality. First mismatches would save both raw streams and stop the runner. The 112-subject case's four initial outfmt-6 hashes were `003259739787e61a4dca66ee5179c69ad59e062187f5fe958e3442bc6848e054`.

Order-sensitive physical fixtures (multi-query, multi-HSP, SEG, heap/containment/result order, 112 subjects) repeat native 2/4/8 runs three times. All 27 eight-subject genetic-code fixtures repeat 2/4/8 runs three times. The diagnostic logs record actual maximum concurrent preliminary work of **2, 4 and 8** respectively, not merely accepted thread counts. The source-ordered collector and Kappa stages remain serial; intermediate HSP retention is bounded by the requested worker count and the capped collector.

## Command-WASI target boundary

| Target and feature | Thread counts | Result | Evidence |
| --- | --- | ---: | --- |
| `wasm32-wasip1`, `--no-default-features` command | 1 | **102/102** bytes equal native and its NCBI contract | [rows](wasm/comparison.jsonl), [environment](wasm/environment.json) |
| `wasm32-wasip1-threads`, `wasm-threads` command | 1/2/4/8 | **462/462** bytes equal native and its NCBI contract, including repeated heap/112-subject/code-32 cases | [rows](wasm/comparison.jsonl), [summary](wasm/summary.json) |
| Plain `wasm32-wasip1` command | 2 | Explicit nonzero rejection; no stdout | [rejection](wasm/plain_wasi_threads_rejected.json) |

The **564/564** WASI comparisons cover seven physical cases (multi-query, multi-HSP, SEG, alternate matrix, containment, heap replacement, 112 subjects) and all 27 eight-subject genetic-code cases in outfmt 0/6/7. The real threaded-WASI diagnostics show maximum simultaneous subject work **2, 4 and 8**. The compiled release SHA-256 values are [recorded](build_artifacts.log): native `87fc615c1c1e78afc52f081a4744f1b01ceafd76270c9e1b37063f510c237339`, plain command-WASI `b68278ee51b6d4e5765d357b6ce07e02a2aae3f40a90ffb1aec0c765b19a57d7`, threaded command-WASI `ec76f820345d6a55310e973d271c985bc621782431aad8e594b02cae24e3fe8f`. The runner uses Node 26.8.2 and repository WASI hosts without altering BLAST output.

No command target was unavailable: both required WASI targets were installed and compiled. `wasm32-unknown-unknown` and reactor/browser APIs are outside this command-WASI Stage F declaration; no certification or unavailability claim is made for them. Stage G owns broader platform and performance certification. NCBI local `-subject -num_threads N` is an output oracle only and is never described as NCBI N-thread performance.

## Reproduction and regression gates

Run from the repository root. Every output directory below must be new. The Stage E comparison-only API oracle may be rebuilt with `bash docs/evidence/tlosan_stage_e/run_stage_e_codes.sh /tmp/tlosan-stagee-api-rebuild`; that helper verifies the fixed NCBI source and copies the compiled oracle to its output directory. The matrix runners require the recorded oracle SHA-256 and verify source and executable hashes before search.

```bash
cargo build --manifest-path LOSAT/Cargo.toml --release
cargo build --manifest-path LOSAT/Cargo.toml --release --bin LOSAT --target wasm32-wasip1 --no-default-features --target-dir LOSAT/target/serial-command
cargo build --manifest-path LOSAT/Cargo.toml --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads --target-dir LOSAT/target/threaded-command
python3 docs/evidence/tlosan_stage_f/run_stage_f_matrix.py physical /tmp/tlosan-f-physical-new
python3 docs/evidence/tlosan_stage_f/run_stage_f_matrix.py codes /tmp/tlosan-f-codes-new --api-oracle /tmp/tlosan-stagee-api-rebuild/tblastn_stage_e_local_oracle
python3 docs/evidence/tlosan_stage_f/run_stage_f_matrix.py codes_multi /tmp/tlosan-f-codes-multi-new --api-oracle /tmp/tlosan-stagee-api-rebuild/tblastn_stage_e_local_oracle
python3 docs/evidence/tlosan_stage_f/run_stage_f_wasm.py /tmp/tlosan-f-physical-new /tmp/tlosan-f-codes-multi-new /tmp/tlosan-f-wasm-new
python3 docs/evidence/tlosan_stage_e/run_stage_e_generated.py /tmp/tlosan-f-generated-new
python3 docs/evidence/tlosan_stage_f/run_stage_f_output_file.py /tmp/tlosan-f-output-file-new
python3 docs/evidence/tlosan_stage_f/run_stage_f_negative.py /tmp/tlosan-f-negative-new
python3 /mnt/c/users/genom/github/losat/docs/evidence/tlosan_stage_d/replay_stage_d_new_oracles.py
cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn -- --test-threads=1
cargo test --manifest-path LOSAT/Cargo.toml --lib utils::threading -- --test-threads=1
cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings
cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check
git diff --check
```

The [Stage D oracle replay](stage_d_replay.log) reproduced **130/130** files; the [14 checksum manifests](stage_d_checksum.log) validated **615/615** input/output files. The [focused suite](focused_test.log) passed **106/106** tests, and the [search-pool lifecycle tests](threading_test.log) passed **3/3**. [Clippy](clippy.log), [fmt/diff/script syntax](static_checks.log) and the supported `-out` and rejection matrices passed. The [Stage F checksum manifest](evidence.sha256) covers every retained evidence file except itself and can be checked with `cd docs/evidence/tlosan_stage_f && sha256sum -c evidence.sha256`. The original Stage E 327/327 byte gate and Stage D numerical/call-order gates are retained, not redefined by a native fingerprint. Baseline/candidate performance samples: **not applicable in Stage F**, which claims observed concurrent work and byte parity, not a speedup. Formal warmup and three timed repetitions, NCBI prebuilt-`-db` speed protocol, and cross-program certification belong to Stage G.

Remaining scope for Stage G: existing BLASTN/BLASTP/TBLASTX regressions, full certification matrix, formal native/WASI performance and resource measurement, any additional product-supported platform, and the final independent release-facing audit. Arbitrary untested option combinations and unsupported composition modes 1/3 remain outside this fixture-scoped gate.
