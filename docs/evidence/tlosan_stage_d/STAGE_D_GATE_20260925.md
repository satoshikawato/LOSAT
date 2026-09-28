# TBLASTN Stage D local-subject gate — 2026-09-25

Status: **PASS for the scoped internal, native, single-thread local-subject Stage D fixture gate**. This record covers the internal Rust local `-subject` Stage C→D path on `feature/tlosan-tblastn-v0.2.0`, based on `267efede89bbdccad90a93a29b3eac0491d28868`. The sole behavior authority is NCBI source commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; the comparison executable is NCBI BLAST+ 2.17.0+ `/home/kawato/micromamba/bin/tblastn`, SHA-256 `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

The Rust path calculates the initial local parameters, runs the natural Stage C search, links and reaps preliminary HSPs, consumes them in the NCBI stream order, performs composition redo or ordinary traceback, and retains score, exact `double` bit score/E-value, coordinates, edit script, raw matrix-adjustment rule, internal and report identities/positives, alignment length, mismatch and gap counts through containment, link/reap, composition heap, hitlist and result sorting. The public TBLASTN CLI still explicitly returns an unimplemented error. Stage E outfmt 0/6/7 byte comparison is outside this gate; tabular output below is a numeric comparison oracle.

## Scope and comparison result

| Case family | Saved comparison-only authority | Rust assertion |
| --- | --- | --- |
| Natural 112-subject search, 13-subject chain, multiple queries, hard SEG | [Prior call-state and owned payload checkpoints](CONTINUATION_20260924.md), [12 redo events and retained results](kappa_heap_rejection_20260925/result_order_20260925), [natural sort, containment and heap fixtures](natural_positive_20260925/run_20260925) | Initial C→D and Kappa incoming HSP order; all 12 redo edit scripts transferred by event to retained HSPs in the broad run; accepted result scores, exact bit/E values, coordinates and report counts; natural E-value sort, postredo containment, and heap replacement. |
| 27 `gc.prt` codes (IDs 1–6, 9–16, 21–29, 30–33), each selected code and code-1 control | [54 calibrated API D/Kappa traces and reports](all_codes_20260925/run_20260925), [27 per-input code-1 CLI controls](all_codes_20260925/cli_controls_20260925), [custom report fields and raw matrix rules](all_codes_20260925/fields_rules_20260925) | Exact C→D input/order, effective length, Kappa incoming order, retained raw score/full-precision E/bit/coordinates/link count, each owned script and rule, all report numeric fields. Code 32 uses `FindGeneticCode(32)` through comparison-only API search/formatting; its code-1 API path is calibrated against the CLI. |
| Six SEG/lowercase hard/soft mask inputs, both composition controls | [12 numeric reports and D/Kappa traces](masking_options_20260925/run_20260925) | C→D input, Kappa input where applicable, exact retained numerics, edit-script ownership and report counts. |
| Five remaining physical Stage C local FASTA families, both composition controls | [10 numeric reports and D/Kappa traces](remaining_local_20260925/run_20260925) | C→D input/order, Kappa input, exact retained numerics, coordinates and report order. Includes ambiguity, lowercase subject, multi-HSP, no-fence x-drop and six-frame merged inputs. |
| Four generated Stage C translated chunk boundaries, both composition controls | [8 numeric reports and D/Kappa traces](extended_chunks_20260925/run_20260925) | Input SHA-256 matches the Stage C manifest, C→D and Kappa input/order, exact retained result values for masked boundary, long chunk and two no-range middle cases. |
| Long 15,000,962-nt subject and three queries including no-hit context | [NCBI ordinary call trace/report](long_subject_20260925/run_20260925) | Query-context effective search space, ordinary batch order, full-translation retry statistical length, exact retained double E/bit, coordinates and report order. |
| Independent `-comp_based_stats`/`-sum_stats` controls, plus alternate protein scoring | [Mode 0/sum true and mode 2/sum false](option_cross_20260925/run_20260925); [BLOSUM45, gap 14/2, word size 2](alternate_matrix_20260925/run_20260925) | Selected branch, C→D input/order, link or direct E-value length, Kappa input where applicable, exact numeric/report results. Default mode 2/sum true and control mode 0/sum false are also exercised in the families above. |

This is a finite fixture gate for the internal single-thread local-subject path. Composition modes 1 and 3 and arbitrary untested option combinations are outside this Stage D fixture gate. ADAPTIVE_CBS stream close returns an explicit internal error; the public TBLASTN CLI remains explicitly unimplemented. No NCBI executable, library, source, FFI or subprocess is in the LOSAT runtime, build or fallback path. NCBI `-db` statistics and headers are not used as local `-subject` expected results. The approved `PD-TLOSAN-LOCAL-GENCODE-32` exception applies only to selected subject translation for non-default codes.

## First differences found and corrected

1. Internal HSP identity on masked query bytes differed from NCBI tabular report counts computed from original aligned sequence strings. Separate internal and report counts now travel with each HSP (`tabular.cpp:971-1021`; `blast_kappa.c:515-526`).
2. Ordinary posttraceback E-value input required the final translated window length, including full-subject retry state (`blast_traceback.c:294,425-433,717-719`; `blast_hits.c:1231-1235`).
3. Initial effective length used actual local subject count, and the four composition/sum-statistics branches use NCBI's independent switch points (`blast_setup.c:964-985`; `blast_traceback.c:1481-1501`; `blast_kappa.c:411-427`).
4. BLOSUM45 gap-14/2 requires its pinned NCBI Gumbel coefficients. Stage D query-context validity also had to use the selected scoring matrix and word length (`blast_stat.c:2778-2803`; `lookup_wrap.c:91-100`).
5. Ordinary multi-query traceback consumed the same-OID HSP lists in reverse query order from NCBI's OID-sorted batch stream (`blast_hspstream.c:91-96,158-204,569-610`). Kappa uses its distinct score-sorted stream.
6. NCBI calls preliminary link/direct E-value on an allocated empty HSP list for each subject before reaping it. The Rust path now enters the same no-op call for all 112 natural subjects, including the 82 without retained preliminary HSPs (`blast_engine.c:870-905`; `link_hsps.c:1774-1777`).
7. Initial direct E-value for translated subjects uses the last nonempty translated frame length. NCBI uses 119 aa, rather than 120 aa, for a 360-nt subject before ordinary traceback; the Rust path now follows that exact initial input and tests short-frame boundaries (`blast_engine.c:728-729,804-813,882-887`; `blast_util.c:508-527,1070-1101`).

## Reproduction and checksums

Run from the repository root. Every runner checks the pinned NCBI source/executable and writes exact commands, input SHA-256 values, diagnostic traces, reports and an output checksum manifest. The code-32 runner builds a temporary **comparison-only** C++ API oracle from the pinned source; no generated binary is committed or shipped.

```bash
python3 docs/evidence/tlosan_stage_d/replay_stage_d_new_oracles.py
bash docs/evidence/tlosan_stage_d/all_codes_20260925/run_ncbi_all_codes.sh /tmp/tlosan-d-all-codes-replay
bash docs/evidence/tlosan_stage_d/all_codes_20260925/run_fields_and_rules.sh /tmp/tlosan-d-fields-rules-replay
python3 docs/evidence/tlosan_stage_d/extended_chunks_20260925/run_ncbi_extended_chunks.py /tmp/tlosan-d-extended-replay
gcc -std=c11 -O0 -shared -fPIC docs/evidence/tlosan_stage_d/ncbi_kappa_traceback_trace.c -ldl -o /tmp/tlosan-d-kappa-probe.so
gcc -std=c11 -O0 -shared -fPIC docs/evidence/tlosan_stage_d/ncbi_d_call_trace.c -ldl -o /tmp/tlosan-d-call-probe.so
python3 docs/evidence/tlosan_stage_d/all_codes_20260925/verify_cli_controls.py /tmp/tlosan-d-cli-controls-replay /home/kawato/micromamba/bin/tblastn /tmp/tlosan-d-kappa-probe.so /tmp/tlosan-d-call-probe.so
find docs/evidence/tlosan_stage_d -name outputs.sha256 -path '*20260925*' -print0 | while IFS= read -r -d '' manifest; do (cd "$(dirname "$manifest")" && sha256sum -c outputs.sha256); done
```

The `all_codes_20260925/fixtures.sha256` file pins every code input; `cli_controls_20260925/outputs.sha256` pins all 27 CLI code-1 output/D/Kappa controls; `fields_rules_20260925/fields_outputs.sha256` pins the exact raw-rule and custom positive/gap reports. The derived custom-format oracle differs from the calibrated API source only in the `CBlastFormat` field-list argument. Its 27 code-1 custom outputs and D call traces were compared byte for byte to the local NCBI CLI for each input. The four generated chunk inputs are reconstructed by [the runner](extended_chunks_20260925/run_ncbi_extended_chunks.py), then checked against their Stage C query and subject SHA-256 values before any oracle run.

## Final local verification

- [Focused TBLASTN test log](stage_d_tblastn_gate_test_20260925.log): `cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn -- --test-threads=1` — **105 passed, 0 failed**. This includes the complete 27-code selected/control run, 112-subject redo ownership, all retained Stage C fixture families, natural sort/containment/heap positives, 15,000,962-nt subject, SEG, scoring and option controls, and short-frame boundaries.
- [Fresh oracle replay](stage_d_oracle_replay_20260925.log): six NCBI runner families, respectively 44, 9, 5, 5, 37 and 30 files, **130 files byte-identical** to the saved comparisons.
- [Checksum gate](stage_d_checksum_gate_20260925.log): **14 manifests and 615 input/output files OK**, including 54 code FASTA inputs, 27 CLI calibration controls and all custom report/matrix-rule traces.
- [Clippy log](stage_d_clippy_20260925.log): `cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings` completed with zero warnings. `cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check` and the tracked `git diff --check` both exited zero. The [staged source/document diff check](stage_d_scoped_diff_check_20260925.log) covers every staged `.rs`, `.py`, `.sh`, `.c`, `.cpp` and `.md` file. The full staged check reports only trailing whitespace inside raw NCBI stderr/build logs and the Cargo test log; the NCBI output bytes are retained to preserve checksum and replay validity, and the Cargo log is retained as emitted by the test command.

## Independent read-only audit and decision

The independent `ncbi_parity_auditor` reviewed the pinned source, staged Rust/source evidence, boundary inputs and call order, owned edit scripts and matrix rules, oracle replay, checksums and validation logs. Its verdict is **SUPPORTED** for the internal native single-thread local `-subject` Stage D fixture gate above; it found no remaining concrete source, input, ordering or payload-ownership discrepancy. It independently checked the final translated-frame-length correction and accepted the documented raw-oracle whitespace exception. The positive comparison is exact for every fixture and supported option combination listed in the scope table; no residual difference is waived.

The decision does not certify the public TBLASTN CLI, Stage E outfmt 0/6/7 bytes, composition modes 1/3, arbitrary untested option cross-products, Wasm, threading or performance. The public CLI continues to return its explicit unimplemented error. The final response records the commit SHA and remote tip.
