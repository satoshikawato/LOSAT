# TBLASTN Stage E local-subject gate — 2026-09-25

Status: **PASS for the scoped Stage E native local-subject gate**. Work resumes from Stage D gate commit `d05503f6b415af213fc66ec0d095afe7cbad5bfd` on `feature/tlosan-tblastn-v0.2.0`. Behavior authority is NCBI C/C++ source commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison executable `/home/kawato/micromamba/bin/tblastn` has SHA-256 `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`. The [Stage C gate](../tlosan_stage_c/STAGE_C_GATE_20260924.md) and [Stage D gate](../tlosan_stage_d/STAGE_D_GATE_20260925.md) supply the preceding search, candidate, statistics and result-order boundaries.

The public Rust `tblastn -query ... -subject ...` path now runs the Stage D search and writes complete outfmt 0, 6 or 7 bytes. Pairwise output uses the retained raw score, full-precision E/bit values, edit script, identity/positive/gap counts, matrix-adjustment method and link count from the same HSP that Stage D retained. The formatter follows the pinned NCBI report path for query and subject identifiers, translation and frame coordinates, subject and HSP order, no-hit query statistics, SEG display case, report headers, spacing and final newline. It buffers the complete result before writing the requested output.

## Exact byte comparisons

| Scope | Input and output record | Result |
| --- | --- | --- |
| 73 physical Stage D local commands, including 112 subjects, multi-query, six-frame, SEG/lowercase hard and soft masks, BLOSUM45, independent composition/sum controls, natural sort, containment and heap replacement | [CLI matrix](cli_matrix_20260925/comparison.jsonl), [environment](cli_matrix_20260925/environment.json) | **219/219** outfmt 0/6/7 stdout pairs byte-identical |
| Four reconstructed 15/30-MB chunk families in both composition modes and the 15,000,962-nt three-query subject | [Generated matrix](generated_20260925/comparison.jsonl), [environment](generated_20260925/environment.json) | **27/27** stdout pairs byte-identical; every generated query/subject SHA-256 checked against Stage C/D records before search |
| All 27 `gc.prt` subject codes, each selected code and same-input code-1 control | [API comparison](all_codes_20260925/comparison.jsonl), [environment](all_codes_20260925/environment.json), [source verification](all_codes_20260925/source_verified.log) | **81/81** code-1 CLI/API calibrations and selected-code API/LOSAT stdout pairs byte-identical, including ID 32 |
| Supported `-out` file path | [File comparison](output_file_20260925/comparison.jsonl) | **3/3** outfmt 0/6/7 files byte-identical, zero stdout bytes |
| Unsupported public option paths | [Negative matrix](negative_20260925/comparison.jsonl) | **14/14** nonzero exit, zero stdout bytes, no output file |

The three search matrices total **327/327** byte-identical complete reports. The physical matrix includes both `uneven_gap_20260924` per-fixture manifests and six earlier `run_20260924` manifests, in addition to the consolidated Stage D manifests.

Each positive comparison row records exact commands, input and output SHA-256 values, exit status and equality; the negative matrix records the failed command, exit status, stdout checksum and absence of an output file. The [evidence checksum manifest](evidence.sha256) covers the 37 other files in this Stage E directory and verifies with `cd docs/evidence/tlosan_stage_e && sha256sum -c evidence.sha256`. The generated-input runner reconstructs subjects from the Stage D recipe and checks the Stage C manifest hashes before comparison. Outfmt 0/6/7 matching includes values, coordinates, translated amino-acid strings, gaps, IDs, statistics, headers, ordering, empty lines and trailing newline because the comparison is over the complete byte stream.

For non-default subject codes, the approved `PD-TLOSAN-LOCAL-GENCODE-32` difference is restricted to applying the selected code during local subject translation/search/reporting. NCBI `-db` statistics and headers are never used as the local `-subject` oracle. The [comparison-only C++ API oracle](tblastn_stage_e_local_oracle.cpp) uses `FindGeneticCode(code)` and the pinned search and `CBlastFormat` paths through all three output formats. Its code-1 output is checked against the CLI for every fixture and format before selected-code comparison. The separately distributed API formatter emits an optional query-coverage column in pairwise summaries despite the requested default sort flag; [the calibration helper](run_stage_e_codes.py) removes only that column, retains raw and calibrated SHA-256 values, and requires complete same-input code-1 CLI equality. Code 32 is never sent to the rejecting NCBI CLI with `-db_gencode 32`.

## First differences corrected

1. NCBI `blast_format.cpp:1541-1557` renders SEG and lowercase query mask locations in lower case in pairwise alignments. The Rust report now carries the search mask positions to display.
2. NCBI `align_format_util.cpp:986-993` uses `%4.1lf` for tabular bit score, retaining the leading space for scores below ten.
3. NCBI `blast_kappa.c:331-342` and `showalign.cpp:3599-3604` choose the displayed composition method per HSP matrix-adjustment rule.
4. NCBI `blast_results.cpp:82-115` returns before initializing Karlin blocks for an invalid query context; `blast_format.cpp:445-477` then writes only its zero search space and blank lines.
5. NCBI `showalign.cpp:3595-3598` writes `Expect(n)` for a linked HSP set; Stage D's retained link count now reaches pairwise output.
6. NCBI `align_format_util.cpp:3152-3161` computes pairwise identity/positive/gap percentages with `int(0.5 + 100*n/d)`, caps non-perfect values at 99, and renders a perfect match as 100.

## Supported public boundary and explicit rejections

Supported: native, single-thread, local file `-query` and `-subject`; outfmt 0/6/7; BLOSUM62 with word size 3, threshold 13, window 40, gap 11/1 and composition 0 or 2; BLOSUM45 with word size 2, threshold 16, window 60, gap 14/2 and composition 0. The finite gate exercises sum-statistics true/false, default/custom SEG, lowercase masks, hard/soft masking, E-value and max-target controls, both supported x-drop parameters, multiple queries and all 27 subject genetic codes. This is a fixture-scoped certification, not an assertion about arbitrary untested cross-products.

The CLI explicitly rejects unported database search, remote search, PSI checkpoints, `tblastn-fast`, subject ranges, ungapped search, multiple threads, nonzero intron length, composition 1/3, unsupported scoring/lookup combinations and matrices, other output formats, and custom tabular fields. Stdin query input, Wasm, multithreaded execution and performance certification are outside this gate. No NCBI binary, library, source, FFI or subprocess enters LOSAT runtime, build or fallback code.

## Reproduction

Run from the repository root. Each output directory must be new.

```bash
cargo build --manifest-path LOSAT/Cargo.toml --release
python3 docs/evidence/tlosan_stage_e/run_stage_e_cli_matrix.py /tmp/tlosan-e-cli-replay
python3 docs/evidence/tlosan_stage_e/run_stage_e_generated.py /tmp/tlosan-e-generated-replay
bash docs/evidence/tlosan_stage_e/run_stage_e_codes.sh /tmp/tlosan-e-codes-replay
python3 docs/evidence/tlosan_stage_e/run_stage_e_negative.py /tmp/tlosan-e-negative-replay
python3 docs/evidence/tlosan_stage_e/run_stage_e_output_file.py /tmp/tlosan-e-output-file-replay
cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn -- --test-threads=1
python3 /mnt/c/users/genom/github/losat/docs/evidence/tlosan_stage_d/replay_stage_d_new_oracles.py
cargo clippy --manifest-path LOSAT/Cargo.toml --lib --tests -- -D warnings
cargo fmt --manifest-path LOSAT/Cargo.toml --all -- --check
git diff --check
```

The Stage D [checksum gate](../tlosan_stage_d/stage_d_checksum_gate_20260925.log) records the exact 14 retained manifests and 615 files. The Stage E repeat [passed all 615](stage_e_stage_d_checksum_gate_20260925.log). The [focused Rust suite](stage_e_tblastn_focused_test_20260925.log) passed **106/106** tests: the Stage D 105-test baseline plus one complete Stage E outfmt 0/6/7 oracle-byte test. The [fresh Stage D replay](stage_e_stage_d_oracle_replay_20260925.log) reproduced **130/130 files byte-identically** across six NCBI runner families. [Clippy](stage_e_clippy_20260925.log), [rustfmt](stage_e_fmt_20260925.log) and [tracked source diff whitespace](stage_e_diff_check_20260925.log) exited zero. The four trailing spaces in the two saved NCBI outfmt 0 oracle files are original output bytes and are intentionally retained; the staged source, scripts, documentation and logs pass `git diff --cached --check` when those two raw oracle files are excluded.

The initial replay call through the differently capitalized `/mnt/c/Users/genom/GitHub/LOSAT` path failed on `manifest.txt` and its checksum file only. Its report and all D/Kappa traces for the masking family were byte-identical; [the captured failure](stage_e_stage_d_oracle_replay_path_case_20260925.log) records that diagnostic. Re-invoking the unchanged runner with the saved `/mnt/c/users/genom/github/losat` path spelling reproduced all 130 files, including manifests and checksum files. This spelling is observable metadata in the saved commands on this case-insensitive mount.

## Final build and independent audit

The final native release executable `LOSAT/target/release/losat` has SHA-256 `da85b0a87a0bc92fcfa1bcc92fb8aa76716424286468a60a92994d2ddfca5e87`. All 327 positive rows in the retained comparison manifests record this same binary SHA-256. The [release build log](stage_e_release_build_20260925.log) and the validation logs above refer to this candidate. The comparison-only NCBI C++ oracle is compiled and run solely by the evidence runner, never by LOSAT's runtime or build.

The independent, read-only `ncbi_parity_auditor` returned **SUPPORTED** for the stated native Linux x86_64, single-thread, local `-subject` fixture and option boundary. The auditor checked the fixed NCBI source and the Rust search/report timing, values, ordering, masks and frames; independently reran an uneven-gap outfmt 0 case, a natural heap-order outfmt 7 case and a genetic-code-32 outfmt 0 case; and checked the Stage D regression and explicit-rejection records. Its three findings were resolved before the final verdict: every matrix row now identifies the final release-binary checksum, eight per-fixture manifests including uneven-gap are included, and the `Expect(n)` Rust branch has the pinned NCBI `blast_seqalign.cpp:1180-1182` source snippet immediately above it. The auditor found no remaining concrete discrepancy. See the [audit record](stage_e_independent_audit_20260925.md).

**Gate decision: PASS**, bounded by the supported options and finite fixtures stated above. Wasm, multithreading, performance and arbitrary untested option combinations have no Stage E certification claim.
