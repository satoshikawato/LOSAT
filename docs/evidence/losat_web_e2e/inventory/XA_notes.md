# Notes, range XA (TBLASTX application and argument layer)

Files: `XA.tsv` (110 rows), this file, scratch in `scratch_XA/`. Oracle/LOSAT runs: for a case `NAME`, `NAME.out/.err/.rc` are NCBI 2.17.0 `tblastx` and `NAME.lout/.lerr/.lrc` are `LOSAT-before tblastx` (made by `batch.sh`, `batch2.sh` (with a stdin file), `envcmp.sh`, `cmp.sh`; `show.sh NAME...` prints rc and the error text with the USAGE block collapsed; case lists are `b_*.txt`). Inputs: `q1.fa`/`s1.fa` (one 2400-nt query vs a mutated copy inside 600 nt of flanks), `q2.fa`/`s2.fa` (2 queries x 3 subjects), `*_lc.fa` (lower case), `*_pd*.fa` (ids). Sweeps: `sweeps.txt` (gencodes, matrix names), `sweeps2.txt` (seg, evalue), `epilog_runs.txt` (epilog lines), `ord_summary.txt` (error order).
Row counts: faithful 15, exception 9, rejected 4, unported 23, divergent 59 (n/a 0). Impact: no row is `high`; `medium` for stdin/missing-subject/real-valued -threshold/word sizes/most missing options.

## 1. Call path (NCBI)

`CTblastxApp::Init` (tblastx_app.cpp:79-89: `CTblastxAppArgs`, `HideStdArgs`, `SetupArgDescriptions`) -> NCBI argument parser (`CArgDescriptions::x_CreateArg` ncbiargs.cpp:2941, value classes `CArg_Int8/Integer/Double/Boolean` ncbiargs.cpp:351-500, constraints `CArgAllowValues*` blast_input_aux.hpp:100-140, dependencies `x_PostCheck` ncbiargs.cpp:3162) -> `CTblastxApp::Run` (tblastx_app.cpp:91): `SetOptions` (blast_args.cpp:3585-3641) -> `InitializeSubject` -> `IsIStreamEmpty` -> `CBlastFormat` -> `PrintProlog` -> batches -> `PrintEpilog`.

Arg groups of tblastx in registration order (tblastx_args.cpp:55-118; this is also the order of `ExtractAlgorithmOptions`): search-strategy, program description, `CBlastDatabaseArgs` (db, lists, `-subject`, `-subject_loc`, `-dbsize`), `CStdCmdLineArgs` (`-query`, `-out`), `CGenericSearchArgs(query_is_protein=true, is_tblastx=true)` (`-evalue`, `-word_size`, `-qcov_hsp_perc`, `-max_hsps`, `-xdrop_ungap`, `-searchsp`, `-sum_stats`; no gap options), `CLargestIntronSizeArgs`, `CFilteringArgs(protein)` (`-seg`, `-soft_masking`), `CMatrixNameArg`, `CWordThresholdArg`, `CHspFilteringArgs`, `CWindowSizeArg`, `CQueryOptionsArgs(nucleotide)` (`-lcase_masking`, `-query_loc`, `-strand`, `-parse_deflines`), `CGeneticCodeArgs` x2, `CFormattingArgs`, `CMTArgs`, `CRemoteArgs`, `CDebugArgs`.

### The order in which NCBI can fail (verified by `ord_*` runs)
1. Parser (CArgException, USAGE text, rc 1): value conversions and constraints, in command-line order; unknown/duplicate/dependent options. LOSAT: clap, rc 2 (approved, exception 1, as long as NCBI also fails at this stage).
2. `-outfmt` text parsed (`ArchiveFormatRequested`, blast_args.cpp:3624): `'X' is not a valid output format` (rc 1) or `Error: Formatting choice is out of range` (rc 255).
3. `ExtractAlgorithmOptions` loop: `-subject` missing (`Either a BLAST database or subject sequence(s) must be specified`), subject file not accessible / read / `Empty CBlastQueryVector` (rc 3); query file and `-out` opened (`Command line argument error: Argument "query"|"out". File is not accessible`); `-seg` text (`Invalid number of arguments to filtering option`, `Invalid input for filtering parameters`); `-query_loc`; the outfmt checks for SAM/AIRR/FASTA and `delim`.
4. `Validate` (`BLAST_ValidateOptions`): Extension, Scoring, Lookup (`Non-zero threshold required`, word size > 4 `Word-size must be less than 6 for protein comparison`), InitialWord (`x_dropoff must be greater than zero`), HitSaving (`expect value or cutoff score must be greater than zero`, `Uneven gap linking of HSPs is allowed for blastx, tblastn, and psitblastn only`). Printed as `BLAST query/options error: MSG` + `Please refer to the BLAST+ user manual.`, rc 1.
5. Run: `Query is Empty!` (seekable empty stream), search, formatting.

LOSAT today runs 1 (clap) for `-seg`, `-outfmt`, `-threshold`, `-word_size`, `-window_size`, `-evalue` and the subject/query requirement; 2-4 are interleaved in `run()` (run_impl.rs:621) in NCBI's order for the file steps only. Rows 13-19 describe the port.

## 2. What the port needs (data and code)

- Integer/real/boolean readers for every TBLASTX option: `ncbi_integer`, `ncbi_double`, `ncbi_constrained_integer`, `ncbi_string_to_bool` already exist (value_parsers.rs:190-290, ncbi_environment.rs). TBLASTX still uses `positive_usize`/`nonnegative_f64`/`positive_i32` (Rust `from_str`).
- C cast of a double to Int4 (`(Int4)threshold`, blast_aalookup.c:245): out of range -> INT_MIN on x86-64 (core/blast_util.rs:261 `ncbi_int4_from_double`).
- C++ `ostream << double` (6 significant digits, `%g` style) for the epilog threshold line.
- Matrices (rows 93, 110): NCBI `util/tables/raw_scoremat.c` (9 standard matrices BLOSUM45/50/62/80/90, PAM30/70/250, IDENTITY, selected by `NCBISM_GetStandardMatrix`), `blast_options.c:1174-1230` (suggested threshold 11 for BLOSUM62, 14 BLOSUM45, 12 BLOSUM80, 16 PAM30, 14 PAM70, 27 IDENTITY, +2 for a translated subject; window 40, BLOSUM45 60, BLOSUM80 25, PAM30 15, PAM70 20), Karlin blocks computed by `Blast_ScoreBlkKbpUngappedCalc` (blast_stat.c:2737). LOSAT has the matrix tables (utils/matrix.rs, stats/tables.rs) and the BLOSUM62 computation (stats/karlin_calc.rs).
- ASN.1 text/binary writers and parsers for strategies (rows 106-107) and outfmt 5-16/18/20 are not in LOSAT; recommended `reject`.
- `CSeq_id` parsing (`-parse_deflines`, row 98) is large.

## 3. Surprises

1. `-threshold` is a Real truncated to Int4 by the lookup table: 12.5 = 12, 0.5 = 0 (exact words only, no error), >= 2^31 -> INT_MIN (all neighbours), `-threshold 0`/`-0` pass the parser and fail in Validate; the epilog prints the double (`12.5`, `1e+10`, `inf`).
2. `+inf`/`+nan`/`1e400` are valid `-evalue` and `-threshold` values (strtod); `inf`, `nan`, ` 5` are not (first character rule). `-evalue 0x10` = 16 and `-evalue 1e` = 1 (NStr::StringToDoublePosix), but `-seg` reals use plain strtod (`1e` fails, `0x10`/`0x1p1` work).
3. The post-parse errors (-seg, -outfmt, -threshold 0, word size, evalue <= 0) are NOT parser errors, so exception 1 does not cover LOSAT's clap versions of them (rows 13-16, 33, 45-46, 62-69).
4. Word size: NCBI accepts 2, 3, 4 and fails for >= 5 with a misleading text ("less than 6"); `-word_size 0x3` works (constraint re-reads with strtod).
5. `-window_size 0` selects the one-hit finder and removes the `Window for multiple hits` line from the outfmt 0 epilog.
6. `-matrix` with an unknown name is not rejected by the parser: the search fails with rc 3 `BLAST engine error: Error: Unknown error code -1 Error: Unknown error code -1 ` (outfmt 0 has printed its prolog). The epilog prints the matrix name as typed.
7. `-outfmt` accepts `+6`, `06`, `6 std`, `7 std`, `0 std`; unknown field names in a custom list are silently dropped; formats 13/14 need `-out`; the out-of-range error has rc 255 and the `Error:` form, not `BLAST query/options error`.
8. An empty PIPE as `-query -` is not "empty" for NCBI (`tellg` < 0): it prints a report for no query; `-subject -` plus `-query -` makes the query an empty pipe.
9. `-num_threads` >= 65535 (or any count the OS cannot start) makes LOSAT fail although NCBI reduces the count and searches on one thread.
10. `-remote` with `-subject` really contacts NCBI ("RID: ..."). My oracle run `o_remote` sent one request; do not repeat it.
11. TBLASTX has no `is_unported_*` list: every missing option is a generic "unknown option or argument" (rc 2). For the options that stay rejected, add an explicit `is_unported_tblastx_arg` list like TBLASTN/BLASTN so that the message names LOSAT's TBLASTX (rows 3-5, 71-85, 105-107).
12. TBLASTX does not call `check_ncbi_application_settings` (main.rs:129; BLASTN does at main.rs:110): DIAG_*, NCBI_CONFIG*, `.ncbirc`/`tblastx.ini` and `BLAST_USAGE_REPORT=bogus` (NCBI SIGABRT, rc 134) are silently ignored.
13. The result_A note of a `-seg "3 -5 -5"` difference did not reproduce on the base binary (q1/s1 and q2/s2; `sweeps2.txt`).
14. The adapter (`web/adapter/src/describe.rs`) generates the parameter list from clap, so every option added to `TblastxArgs` appears in `describe`; options that the GUI must not offer (anything that reads files, `-import_search_strategy`, `-export_search_strategy`, `-remote`, `-html`) need to be filtered through `ADAPTER_OWNED` or kept out of clap. `validate` = clap + `check_options`, so every new post-parse check must live in `check_options` (not in `run`), or the web validate and the CLI diverge.

## 4. Decisions for the session

- Rows with a `reject:` proposal: -version (build date), hidden toolkit options (-logfile, -conffile, -version-full, -xmlhelp), strtod hexadecimal reals in `-seg` (BLASTN precedent E2c §G), custom outfmt field lists and outfmt 1-5, 8-12, 15, 16, 18, 20 (scope), `-db` family and `-entrez_query` (no database), `-remote` (network), `-query_loc`/`-subject_loc` (S11), search strategy import/export (ASN.1). The common rules prefer `port` where NCBI is deterministic; I chose `reject` only where reproduction needs an object model or a non-reproducible value. Please confirm `-import/-export_search_strategy`, `-html` (port:L) and `-parse_deflines` (port:L, BLASTN rejects it).
- strtod spelling `0x10` / `1e` for `-evalue`: BLASTN keeps an explicit rejection; TBLASTX could do the same or port `StringToDoublePosix`. Row 31 proposes `port:S` for the sign/inf/nan forms and leaves hex/exponent forms to the BLASTN precedent.
- Row 54: clamp the thread pool to the CPU count (exception 2 allows honouring the count but not failing).
- Engine-range work that the argument layer only passes through: word sizes 2 and 4 (row 38), one-hit finder (row 42), sum statistics off (row 91), soft masking (row 49), culling (XC), best hit (rows 94-95), `-strand` frames (row 97), lowercase masks (row 96), effective lengths for `-dbsize`/`-searchsp` (rows 86-87), `-qcov_hsp_perc`/`-max_hsps` (rows 88-89), `-xdrop_ungap` (row 90).
- Messages that NCBI prints only after the prolog (outfmt 0) or never (dryrun) need the order of rows 13-14 respected: the outfmt 0 prolog is written after `SetOptions`, so every Validate error comes before it; Empty CBlastQueryVector for the subject comes before it too (stdout empty), while a header-only/empty query error (rc 3) comes after the prolog (S08 row 69).
