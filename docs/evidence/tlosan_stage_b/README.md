# TLOSAN v0.2.0 Stage B — CLI and genetic code

Status: Stage B input and translation-table boundary only. Every TBLASTN search
still exits with `TBLASTN local search is unimplemented (Stages C-E)` after
argument validation. This record makes no TBLASTN search/output parity or
release claim.

## Authority and Stage A gate

- LOSAT branch: `feature/tlosan-tblastn-v0.2.0`; starting tree was clean.
- NCBI authority: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`.
  The checked-out NCBI worktree is modified. A fresh `git archive` of that
  commit passed every checksum in
  [Stage A ncbi_source.sha256](../tlosan_stage_a/ncbi_source.sha256).
- All retained Stage A CLI and C++ API output checksums passed; logs are
  [source](stage_a_source_checksum.log), [CLI](stage_a_cli_checksum.log),
  and [API](stage_a_api_checksum.log). The starting LOSAT parent, platform,
  toolchain, and tested release-binary hash are in
  [build_provenance.txt](build_provenance.txt). Stage A's
  code-32 `FindGeneticCode(32)` API oracle produced 100.000% identity,
  251 bits and `3.39e-93` for the TAG-rich subject; its code-1 control
  produced 95.833%, 231 bits and `1.75e-85`. The same pinned API outputs
  are asserted by the Rust unit test and the comparison-only checker here.
- The approved local-subject code-32 exception remains limited to TBLASTN;
  TBLASTX's CLI still accepts only NCBI's 26 CLI IDs.

## Reproduce

From the repository root, use a fresh output directory:

```bash
cd LOSAT
cargo fmt --all -- --check
cargo test --lib algorithm::tblastn::args::tests
cargo test --lib utils::genetic_code::tests -- --nocapture
cargo test --test cli_v2
cargo test --test unit_tests unit::tblastx::translation
cargo test --test unit_tests unit::tblastx::args
cargo build --release
cargo clippy --all-targets -- -D warnings
cd ..
python3 docs/evidence/tlosan_stage_b/verify_stage_b.py \
  --output /tmp/tlosan-stage-b-fresh
python3 LOSAT/tests/audit_tblastx_v010.py \
  --manifest docs/evidence/tlosan_stage_b/tblastx_genetic_manifest.tsv \
  --output-dir /tmp/tlosan-stage-b-tblastx-fresh --repeatability-runs 3
```

The checker reads pinned `gc.prt` and `blast_stat.c` blobs with `git show`,
compares them with [gc_prt_27.tsv](gc_prt_27.tsv) and the Rust matrix gap
tables, asserts the Stage A C++ API code-32 control, and runs the built LOSAT
binary through the CLI matrix. NCBI binaries and
libraries are never called from LOSAT's build/runtime paths.

## Genetic-code result

[genetic_code_comparison.tsv](run_20260923/genetic_code_comparison.tsv)
records all 27 IDs and the number of stop codons in each 64-codon NCBI table.
The Rust unit test checks each of the 1,728 codons in TCAG order against that
pinned fixture; [rust_gc_test.log](rust_gc_test.log) records one passing line
per ID. Representative ambiguous codons and every explicit stop codon are
checked for each ID, along with code 32 TAG/TAR and every invalid `u8` ID. The previous duplicate tables in `utils/genetic_code.rs` and
`core/gencode_singleton.rs` disagreed with `gc.prt` at ATA for IDs 2, 3, 5,
13, and 21, and at TGA for ID 31. ID 32 was absent. The core module now
re-exports the single table owner; invalid IDs return an error at the fallible
boundary, and internal validated callers fail rather than silently use code 1.

## CLI acceptance and rejection

The complete 27-ID and invalid-ID run is in
[cli_accept_reject.tsv](run_20260923/cli_accept_reject.tsv).
[matrix_gap_comparison.tsv](run_20260923/matrix_gap_comparison.tsv) checks
all nine NCBI protein matrix tables, their allowed gap pairs, and the
preferred defaults used when a cost is omitted. NCBI also derives the word
threshold and window size from the selected matrix. The Rust tests check
BLOSUM62 `13/40`, PAM30 `18/15`, BLOSUM45 `16/60`, IDENTITY `29/40`,
and word-size-5/6/7 threshold overrides `19.3/21.0/20.25`.

| Input | Rust CLI result | NCBI source / product boundary |
| --- | --- | --- |
| `-task tblastn`, the 26 NCBI CLI IDs from `gc.prt` | Validated; search explicitly unimplemented | `tblastn_args.cpp:55-62`, `blast_args.cpp:997-1056` |
| `-db_gencode 32` | Validated; search explicitly unimplemented | `gc.prt:340-347`, PD-TLOSAN-LOCAL-GENCODE-32; pinned NCBI CLI rejects 32 |
| Invalid code 0, 7, 8, 17, 20, 34, 255, 256, −1, text | Parse error naming genetic code; no code-1 substitution | `blast_args.cpp:997-1013`, `blast_aux.cpp:588-600` |
| `tblastn-fast`, PSI checkpoint, DB search, remote | Explicit unsupported error | `tblastn_args.cpp:55-62,125-152`, `blast_args.cpp:2364-2385` |
| `-db` with `-subject`, `-in_pssm` with `-query`, remote with subject range | CLI conflict error | `blast_args.cpp:2364-2385`, `tblastn_args.cpp:125-152` |
| `-ungapped` with composition adjustment | NCBI-sourced incompatibility error | `blast_args.cpp:872-876` |
| Word size 8, threshold 0, E-value 0 | NCBI-sourced value error | `blast_options.c:1303-1348,1518-1523` |
| Invalid matrix, unsupported gap pair, IDENTITY with word size 6 | NCBI-sourced matrix/score error | `blast_stat.c:183-428,577-586,2952-2994`; `blast_options.c:913-943,1778-1796` |
| PAM30, IDENTITY, and valid explicit gap pair | Validated; search explicitly unimplemented | `blast_args.cpp:258-273`; `blast_stat.c:3374-3399` |
| Matrix or word-size default adjustments | Effective threshold/window and gap costs follow NCBI | `blast_args.cpp:289-304,486-497,599-618`; `blast_options.c:1174-1236` |
| Other registered NCBI options without Rust behavior | Explicit unsupported option error | `tblastn_args.cpp:64-129`, pinned NCBI `-help` |
| `-outfmt 0/6/7` | Accepted as future output profile; search always fails | `blast_args.cpp:2657-2660` |

Pinned NCBI CLI additionally accepted `-comp_based_stats 2u`, `2xyz`, and
`bad`: its source switch examines the first byte and has no default error.
Rust accepts these tokens and applies the same first-byte rule for the
`-ungapped` incompatibility. NCBI accepted the tested Boolean forms
`true/false`, `T/F`, `yes/no`, `y/n`, and `1/0`; Rust uses the same grammar.

## Rust ↔ NCBI mapping

| Rust owner | Pinned NCBI owner |
| --- | --- |
| `src/algorithm/tblastn/args.rs`: task, options, defaults, rejects | `tblastn_args.cpp:45-152`; `tblastn_options.cpp:53-85`; `blast_args.cpp:258-273,997-1056,2364-2385,2657-2660,3425-3427`; `blast_options.c:913-943,1174-1236,1303-1348,1518-1523,1778-1796` |
| `src/algorithm/tblastn/scoring.rs`: allowed gap pairs, preferred costs, suggested threshold/window | `blast_stat.c:183-428,577-586,2952-2994,3374-3399,3578-3638`; `blast_options.c:1174-1236` |
| `src/cli.rs`, `src/main.rs`, `src/algorithm/mod.rs`: dispatch and unported option error | `tblastn_args.cpp:45-129`; `tblastn_app.cpp:288-301` |
| `src/utils/genetic_code.rs`: 27 tables and codon/ambiguity translation | `gc.prt:105-357`; `blast_aux.cpp:588-613`; `blast_util.c:369-424`; `blast_encoding.c:80-103` |
| `src/core/gencode_singleton.rs`: shared table adapter | `gencode_singleton.c:65-69` |

## Existing TBLASTX regression and remaining work

The focused existing TBLASTX audit uses the committed v0.1.0 manifest cases
`p12`, `p14`, `d04`, and `d06`, with a fresh temporary output directory and
three repeats where the runner requires them. Its classifications and output
hashes are retained in [tblastx_classifications.tsv](tblastx_classifications.tsv).
The default and query-code-4 cases are `EXACT_TEXT` (4,319/4,319 and
14,871/14,871 rows). Both local subject-code-4 cases pass the approved
TBLASTX subject genetic-code exception; their LOSAT hashes equal the frozen
canonical hashes. These are focused genetic-code regressions, not full-suite
TBLASTX certification. The unit test logs for TBLASTX translation and
argument parsing are also retained here: [translation](rust_tblastx_translation_test.log)
and [arguments](rust_tblastx_args_test.log).

Stage C must implement the NCBI TBLASTN local protein query and six-frame
subject search, exact frame/context ordering, candidate/HSP construction,
re-evaluation, and first-stage trace diagnostics. It must carry the selected
subject code through that entire search path, including 32. Stages D and E
remain responsible for linking/statistics/composition adjustment and
outfmt 0/6/7 reporting. Database, PSI, fast, remote, stdin, and unported
option paths remain explicit errors until separately implemented and verified.
