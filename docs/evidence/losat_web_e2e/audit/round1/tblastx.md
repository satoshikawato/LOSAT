# TBLASTX application flow and option checks (auditor: TBLASTX)

Setup: `D=/home/kawato/.cache/losat-web-gui-target/s08p/audit`, `W=$D/work/tblastx`; inputs and the harness are in `$W`
(`q1.fna q2.fna qm.fna qbig.fna` queries, `s1.fna s3.fna sm.fna sbig.fna` subjects, `cmp.sh` runs NCBI and LOSAT and compares
stdout, stderr and exit status; run results are kept in `$W/r/`). NCBI = `/home/kawato/micromamba/bin/tblastx`, LOSAT = `$D/LOSAT tblastx`.
About 1,000 option/value probes were run (every comparison below is stdout+stderr+exit unless said otherwise).

Counts: high 2, medium 3, low 5 (10 findings).

## Findings

### TX-1 (high, CONFIRMED) -window_size >= 2^31 - (translated query length) gives hits where NCBI gives "No hits found"
- LOSAT: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:2021` (`while diag_array_size < (query_length + window)`, i32) and the
  word finder that follows (`diag_offset: window`, run_impl.rs:3980-4090).
- NCBI: `src/algo/blast/core/blast_extend.c:52-61` (`while (diag_array_length < (qlen+window_size))`, `diag_table->offset = window_size`);
  `Int4 qlen+window_size` overflows, the loop is skipped, diag_array_length = 1, mask 0, and the later `offset + subject_offset` arithmetic wraps.
  NCBI's result is deterministic on the shipped binary (3 identical runs): no hits.
- `-window_size` is in AUTHORITY.md §K as "supported (>= 1)"; the value parses (`<= 2147483647`) and is not rejected.
- Repro (in `$W`):
  ```
  /home/kawato/micromamba/bin/tblastx -query q1.fna -subject s1.fna -window_size 2147483640 -outfmt 6   # empty stdout, exit 0
  $D/LOSAT tblastx          -query q1.fna -subject s1.fna -window_size 2147483640 -outfmt 6   # q1 s1 40.000 15 9 0 290 246 20 64 0.46 18.9 , exit 0
  /home/kawato/micromamba/bin/tblastx -query q2.fna -subject s3.fna -window_size 2147483647 | tail   # ***** No hits found *****  (outfmt 0, 919 bytes)
  $D/LOSAT tblastx          -query q2.fna -subject s3.fna -window_size 2147483647                     # 1955 bytes: s2 335 bits 5e-96 alignment
  ```
  Same for `0x7fffffff`, 2147483600, 2147483620, 2147483630, 2147483640, 2147483646, 2147483000 (q2/s3: the boundary is INT_MAX - qlen, about 2147481400
  for q2; 2147480000 hangs, see TX-2). Outputs are in `$W/r/` (ids 8599f497 960a7f28 138741e6 ce409afd 52bfd4f5 bd7c62ae 695341bf 7183abe2).
- Under the NCBI defect policy (deterministic NCBI behaviour is reproduced; otherwise reject) LOSAT must either reproduce "no hits" or reject
  `qlen + window > INT_MAX`.

### TX-2 (medium, CONFIRMED, known "D8" not implemented) accepted -window_size in (2^30 - qlen, 2^31 - qlen) loops forever, no explicit rejection
- LOSAT: run_impl.rs:2020-2022 (`diag_array_size <<= 1` on i32 never reaches `query_length + window`). NCBI: blast_extend.c:54-57 (same infinite shift loop; Int4 wraps 2^31 -> INT_MIN -> 0).
- AUTHORITY.md §M D8 says these values are to be rejected ("未実装"), but the item is not in §N and no code rejects them.
- Repro: `timeout 60 /home/kawato/micromamba/bin/tblastx -query q1.fna -subject s1.fna -window_size 1500000000` and
  `timeout 60 $D/LOSAT tblastx -query q1.fna -subject s1.fna -window_size 1500000000` both run at 100% CPU (RSS 6 MB) until killed (exit 124).
  Also 1073741823 (q2), 1073741824, 1073741825, 2000000000, 2147480000. LOSAT should fail with "... not supported by LOSAT's TBLASTX".
- Side note: `-window_size 536870912..1000000000` makes `vec![DiagStruct::default(); 2^30]` per worker (8-16 GB); it matched NCBI here but an allocation failure aborts (approved exception).

### TX-3 (high, CONFIRMED) `-query - -subject -` (and `-subject -` with the default query `-`) from a seekable stdin: LOSAT prints nothing, NCBI prints the report
- LOSAT: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:907` `let seekable = std::io::Seek::stream_position(&mut query_file).is_ok();`
  BLASTN (`blastn/blast_engine/run.rs:5343`), TBLASTN (`tblastn/args.rs:699`) and BLASTP (`blastp/blast_engine.rs:4454`) all carry
  `!(args.query == "-" && args.subject == "-") && ...`; TBLASTX is the only one without it.
- NCBI: `src/app/blast/blast_app_util.cpp:856-860` (`orig_p = in.tellg(); if(orig_p < 0) return false;`): after the subject has been read from `cin`,
  `cin` is at EOF, `tellg()` fails, the query is "not empty", and NCBI prints the prolog/epilog of a run with no query (`tblastx_app.cpp:132-137`).
- Repro (in `$W`, stdin from a regular file):
  ```
  /home/kawato/micromamba/bin/tblastx -query - -subject - < sm.fna   # exit 0, stdout 617 bytes ("TBLASTX 2.17.0+ ... Database: User specified sequence set (Input: -).  63 sequences; 55,000 total letters ...")
  $D/LOSAT tblastx                    -query - -subject - < sm.fna   # exit 0, stdout 0 bytes, stderr "Warning: [tblastx] Query is Empty!"
  /home/kawato/micromamba/bin/tblastx -subject - < qm.fna            # exit 0, 613 bytes (5 sequences; 3,100 total letters); LOSAT: 0 bytes + "Query is Empty!"
  cat sm.fna | $D/LOSAT tblastx -subject -                           # exit 1 explicit rejection ("an empty query from a stream without a position ... not supported by LOSAT's TBLASTX"): acceptable
  ```
  Silent: exit 0 with a different stdout.

### TX-4 (medium, CONFIRMED) -seg locut/hicut with a malformed number: LOSAT says "not supported by LOSAT's TBLASTX" where NCBI rejects with "Invalid input for filtering parameters"
- LOSAT: `LOSAT/src/blastinput/app.rs:417-427` (`Err(NcbiDoubleError::Unsupported)` arm bails with the LOSAT text) and `value_parsers.rs:260-285`
  (`ncbi_string_to_double` returns `Unsupported` for every token that starts with a digit/sign/point but is not a Rust `f64`).
- NCBI: `blast_args.cpp:423-428` (CStringException eConvert -> `Invalid input for filtering parameters`, exit 1): `NStr::StringToDouble` fails when `strtod` does not consume the whole token.
- Same exit code (1), different text; also the NCBI FASTA-Reader warnings are dropped. Tokens: `2.5x`, `2,2`, `+`, `-`, `.`, `.e1`, `1e`, `1e+`, `1e-`, `0x`, `1_0`, `2.5,`, `2.2x`.
- Repro: `for s in '12 2.2 2.5x' '12 2,2 2.5' '12 1e 2.5' '12 0x 2.5'; do $W/cmp.sh -q q2.fna -s s3.fna -- -seg "$s"; done`
  NCBI: `BLAST query/options error: Invalid input for filtering parameters` / `Please refer to the BLAST+ user manual.` (exit 1);
  LOSAT: `Error: the SEG locut or hicut "2.5x" (not a finite decimal number) is not supported by LOSAT's TBLASTX` (exit 1).
- Correct classification: only tokens that `strtod` reads fully but LOSAT cannot (hex `0x10`, `0X1p3`, `+inf`, `-infinity`, `-nan`, `1e400`, `1e99999`) are "not supported"
  (these were verified: NCBI runs them, LOSAT rejects explicitly, which is allowed).

### TX-5 (medium-low, CONFIRMED) LOSAT's explicit rejections come before NCBI's own error (text or exit code differs, both fail)
- LOSAT: `run_impl.rs:696-700` (`read_nucleotide_subjects` hard errors before `check_ncbi_options`), `run_impl.rs:631-650` (`report_format` rejection before the option checks).
- Examples, all `cmp.sh` runs in `$W` (subject/query files `xs.fna` = `>a\nACGTNNXX`, `nohdr2.fna` = `ACGT`, `empty.fna`):
  - `-s xs.fna -- -threshold 0`: NCBI stderr `FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 7-8` + `BLAST query/options error: Non-zero threshold required`, exit 1;
    LOSAT `Error: subject record 1 (a) has 'X' at residue 7, which is not an IUPAC nucleotide letter; ... not supported by LOSAT's TBLASTX`, exit 1. Same with `-seg abc`, `-evalue 0`, `-word_size 5`.
  - `-s nohdr2.fna -- -threshold 0` (also -seg abc / -evalue 0 / -word_size 5): NCBI option error (exit 1) vs LOSAT `failed to read subject FASTA` rejection (exit 1).
  - `-outfmt 5 -threshold 0` and `-outfmt 5 -seg abc`: NCBI `Non-zero threshold required` / `Invalid number of arguments to filtering option`; LOSAT `output format 5 is not supported by LOSAT's TBLASTX`.
  - `-q empty.fna -outfmt 5`: NCBI exit 0 (`Query is Empty!`), LOSAT exit 1 (fine, value is a listed rejection). `-s empty.fna -outfmt 5`: NCBI exit 3 (`Empty CBlastQueryVector`), LOSAT exit 1 (exit code differs).
- Order of the NCBI-reproduced checks (`-seg`, `-outfmt`, hit-list warning, `-threshold`, `-word_size`, `-evalue`, file errors, `Query is Empty!`) was checked with 50+ multi-error combinations and matched.

### TX-6 (low, CONFIRMED) -num_threads rejections do not use the required phrase
- LOSAT: `LOSAT/src/utils/threading.rs:67-82` (`validate_threads`), `with_search_pool`.
- NCBI accepts any int >= 1 (warnings only; ignored with -subject). LOSAT:
  - `-num_threads 2147483647` -> `Error: requested 2147483647 threads exceeds Rayon maximum 65535, which is not supported by LOSAT` (exit 1);
  - `-num_threads 1000` under `ulimit -v 12000000` -> `failed to build tblastx pool with 1000 threads; a thread count that the system cannot start is not supported by LOSAT` (exit 1);
  - builds without parallel/wasm threads: `unsupported num_threads=N: this build does not support parallel search` (no "not supported by LOSAT" at all) - SUSPECTED, not built.
  The contract text is "not supported by LOSAT's TBLASTX". (`-num_threads 2..200` matched NCBI's stdout; the NCBI thread warnings are an approved difference.)

### TX-7 (low, CONFIRMED) wrong NCBI line range in a comment
- `LOSAT/src/stats/protein_options.rs:492` cites `blast_options.c:1515-1520` for the "expect value or cutoff score must be greater than zero" check; the check is at `blast_options.c:1518-1523`
  (1515 is `return BLASTERR_OPTION_VALUE_INVALID;` of the previous check, the message string is at 1521). (Shared by BLASTP/TBLASTN; the TBLASTX copy at `run_impl.rs:~795` cites 1518-1523 correctly.)

### TX-8 (low, SUSPECTED; ABI v1 is frozen, plan TD-1) web ABI v1 `-evalue` for TBLASTX uses Rust `f64::parse`
- `LOSAT/src/web_api.rs:403-411` (`parse_tblastx_args`): `"inf"`, `"nan"`, `"infinity"` (no sign) are accepted and run; NCBI rejects them at argument parsing (first character must be a digit, point or sign, `ncbistr.cpp:1313-1318`).
  Also `-threshold`/`-window_size`/`-seg` are not accepted by v1 at all ("unsupported tblastx argument for web API"). Not built/run here.

### TX-9 (low, CONFIRMED by code reading + CLI runs) adapter `validate` is stricter than the CLI for an empty query
- `web/adapter/src/run.rs:88-90` calls `check_options` = `check_ncbi_options` + `check_losat_limits`; the CLI only applies `check_losat_limits` in `search()` after the `Query is Empty!` return (`run_impl.rs:907-925`, `search` 1021+).
  `-q empty.fna -s s3.fna -- -window_size 0` / `-word_size 2` / `-word_size 4 -threshold 15`: NCBI and LOSAT CLI exit 0 with `Warning: [tblastx] Query is Empty!` (SAME); `validate` rejects them.
  `validate` also omits `validate_threads`, `check_unsupported_environment`, `query_batch_size` (environment/threads; `num_threads` is adapter-owned, `describe.rs:14`).
  For every other option value validate and the CLI share `try_parse_from` + `check_ncbi_options` + `check_losat_limits`, so texts and order are identical (static comparison; the adapter has no Cargo.toml in the snapshot and was not built).

### TX-10 (low, CONFIRMED) AUTHORITY.md §K / M inconsistency
- §K lists "`-window_size` of 1 or more" as supported without the upper limit, while §M D8 says the huge values are (to be) rejected and §N does not list the open item (see TX-1, TX-2).

## Checked OK (identical stdout, stderr and exit status to NCBI, or an allowed rejection)
- `-threshold`: 0, -0, +0, -0.0, 0e5, 00.0, 0x0 (rejected as not supported), 0.0001, 0.5, 0.9, 11, 12, 12.5, 13, 13.9999, 14, 100, 1000, 1e9, 1e10, 2147483647, 2147483648 (INT_MIN path), 4294967296, 1e19, 1e400, +inf/+INF/+Infinity, -1/-inf/+nan/nan (parse errors, exit 2 vs 1 approved), 123456.7, 1e15, 1e-300, 5e-324, 1e308 ... 35 more spellings with the outfmt 0 epilog `Neighboring words threshold` (`%g`) identical.
- `-word_size` 1/0/-1 (parse), 2 and 4 (LOSAT rejection), 3, 5, 6 (NCBI `Word-size must be less than 6 for protein comparison` text and exit 1 identical).
- `-evalue`: 0, -1, 1e-400 (NCBI error text identical), 1e-20, 1e-300, 5e-324, 1e300, 1e400, 1e5, +inf, +nan, -nan, +nan(abc), 0x10 (LOSAT explicit rejection), inf/nan (parse error), with outfmt 0/6/7.
- `-seg`: no, yes, default, `12 2.2`, `12 2.2 x`, abc, `yes extra`, 4 tokens, empty/space values, `YES/NO/true/T/L/n/y/1/0` (all NCBI error), leading/trailing/double spaces and tabs, window `+12 -0 012 1000 2147483647 -2147483648`, 2147483648 and 4294967296 (error), `0 0 0`, `-5 2 2`, `12 -1 -1`, `12 2.5 2.2`, `12 1e5 2.5`, `12 5e-324 2.5`, `12 1e-400 2.5`, `12 nan 2.5`, `12 inf 2.5`, `12 infinity 2.5` (NCBI error identical), `-seg=...` form.
- `-window_size` 1, 2, 3, 10, 39, 40, 41, 100, 1000, 100000 ... 1000000000 (and 536870911/2, 715827882/3), 0/-0/+0/00/0x0 (LOSAT rejection), hex/octal/`+5`, -1, 2147483648, `0x`, `1e2`, ` 5` (parse errors).
- `-culling_limit` 0, 1, 2, 3, 4, 5, 6, 7, 10, 25, 100, 100000, 2147483643..2147483647 (the +3 Int4 wrap), hex, `+2`; -1, 2147483648 (parse). 125+ combinations with `-max_target_seqs`, `-evalue`, `-window_size`, `-threshold`, `-seg`, outfmt 0/6/7 on `qm/sm` and on the repeat-rich `qbig/sbig` all identical.
- `-max_target_seqs` 1..10, 20, 1000, 100000, 1000000000, 1073741799..1073741825, 1431655765/6, 2000000000, 2147483597..2147483647 (>10 subjects per query, so the prelim list size wrap would show), hex, `+3`, `03`; 0, -1, 2147483648, `0x`, `-0` (parse). The `Examining 5 or more matches is recommended` warning and its order before the threshold error match.
- `-outfmt`: 0, 6, 7, ` 6`, `+6`, `06`, `006`, `+0`, `00`, `6\t`, `\t6`, `6 ` , `0 delim=,`, `0 std`, `0 qseqid`, `6 delim=`, `7 delim=`, `6 `; 19, 21 (exit 1), 22, 99, -1 (exit 255) identical text; abc, `0x6`, `6x`, `6,`, `1e1`, `6.0`, empty, space (exit 1 identical); 13 and 14 (identical `Please provide a file name for outfmt 13.`); 17 (identical); 1-5, 8-12, 15, 16, 18, 20, `6 std`, `6 delim=,`, custom fields: explicit LOSAT rejection (NCBI succeeds).
- `-num_threads` 1, 2, 4, 8, 64, 100, 200 (stdout identical; NCBI's two warnings are an approved difference), 0/-1 (parse), hex/`+2`.
- Unported options: every `-name` that `tblastx -help` lists (57) plus all 120 names of `cmdline_flags.cpp` and `-help-full -xmlhelp -logfile -conffile -version-full*`: options NCBI accepts (`-dbsize -export_search_strategy -import_search_strategy -line_length -matrix -max_hsps -max_intron_length -qcov_hsp_perc -query_loc -searchsp -soft_masking -subject_loc -sum_stats -xdrop_ungap -lcase_masking -strand -h -version -db ...`) give "not supported by LOSAT's TBLASTX"; options NCBI rejects as unknown give a parse error (exit 2 vs 1, approved).
- Check order with 2+ wrong options (50+ combinations of `-seg -threshold -word_size -evalue -outfmt`, missing/unreadable subject/query/-out, empty subject/query, `-outfmt 13/14/99`): first error text and exit code identical (exceptions in TX-5).
- Syntax forms: `-evalue=1e-5`, `-outfmt=6`, `-seg=12 2.2 2.5`, duplicates, missing values, extra positional (parse errors exit 2 vs 1, approved), `--evalue` (rejected by both).
- Inputs with default and non-default options: lowercase, mixed case, IUPAC, N-only, CRLF, `U`, 1/2/3/5 nt queries, multi-record, long deflines (`Title ends with at least 50 valid amino acid characters` warning identical), query `-`/stdin file (identical), `-out file`, `-out -`, `/dev/null`, `/dev/full` (exit 6 identical), 100- and 250-character file names.
- Random option combinations: 80 on `qm/sm` and 45 culling-centred on `qbig/sbig` (outfmt 0/6/7): stdout identical in all; only the approved `num_threads` warning differs.
- Adapter: `validate` for TBLASTX = `parse` + `check_options`; the CLI parser/`check_ncbi_options`/`check_losat_limits` and their order are shared (TX-9 for the differences).
- NCBI line references: 100+ NCBI file:line + snippet comments in `tblastx/*`, `blastinput/app.rs`, `value_parsers.rs`, `cli.rs`, `blastn/input.rs`, `stats/protein_options.rs`, `report/*` were machine-checked (every snippet line occurs inside the cited range; script `$W/refcheck.py`, results `$W/refs1.txt refs2.txt`); the following were also verified by hand against the NCBI source with line numbers:
  blast_hits.c:1996, blast_options.c:1303-1311 / 1358-1364 / 1518-1523, blast_aalookup.c:237 and 245, blast_args.cpp:2975-2978 / 3335-3340 / 3624-3627 / 2745-2748 / 3425-3427 / 578-583 / 396-406 / 375-384 / 423-428,
  tblastx_app.cpp:132-137, blast_app_util.cpp:856-860, aa_ungapped.c:214-224, blast_options_local_priv.hpp:1320-1341, setup_factory.cpp:330-341, blast_extend.c:52-61, blast_options.c:1749-1810 (order Extension, Scoring, Lookup, InitialWord, HitSaving). Only TX-7 was wrong.
