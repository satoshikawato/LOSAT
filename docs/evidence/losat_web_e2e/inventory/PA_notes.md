# PA notes: BLASTP application and argument layer

Rows: `PA.tsv` (118). Evidence: `scratch_PA/` (`mx.py` runs NCBI blastp 2.17.0 and `LOSAT-before blastp` on the same
arguments and writes `r/<group>_<label>.{n,l}.{out,err,rc}`; groups are named in the `evidence` column; `envt*.py`
for environment variables; `findings.txt` is the running log; inputs `q1.faa` (one query, first 660 aa of
AvCLPV), `q3.faa` (3 queries), `s6.faa` (6 subjects), `s_full.faa` (120 subjects), `qbig_*.faa` (a 12000-residue first
query for batch tests)). Nothing in a repository was modified.

## 1. NCBI call path followed

`blastp_app.cpp` `NcbiSys_main -> CBlastpApp().AppMain` (`ncbiapp.cpp:831`):
1. `AppMain` pre-parses the command line for `-logfile`, `-conffile`, `-version`, `-version-full[-xml|-json]`,
   `-dryrun` (ncbiapp.cpp:845-1025), resets the environment (`DIAG_POST_LEVEL`, `ABORT_ON_THROW`, 1029-1045), loads
   the registry (`blastp.ini`, `.ncbirc`), then `Init` -> `SetupArgDescriptions(CBlastpAppArgs::SetCommandLine())`.
   Any `CArgException` prints USAGE plus `Error: ...` and exits 1 (ncbiapp.cpp:1087-1113; PD-1 exception).
2. `CBlastpAppArgs` (blastp_args.cpp:44-110) registers, in this order (it is also the order of
   `ExtractAlgorithmOptions`): search-strategy args, program description, task (blastp/blastp-fast/blastp-short),
   `CBlastDatabaseArgs` (db, filters, `-subject`, `-subject_loc`, `-dbsize`), `CStdCmdLineArgs` (`-query`, `-out`),
   `CGenericSearchArgs(protein, no tblastx, no sum stats)` (evalue, word_size, gapopen, gapextend, qcov_hsp_perc,
   max_hsps, xdrop_*, searchsp), `CFilteringArgs` (seg, soft_masking), `CMatrixNameArg`, `CWordThresholdArg`,
   `CHspFilteringArgs`, `CWindowSizeArg`, `CQueryOptionsArgs` (lcase_masking, query_loc, parse_deflines),
   `CFormattingArgs`, `CMTArgs`, `CGappedArgs`, `CRemoteArgs`, `CCompositionBasedStatsArgs`, `CDebugArgs` (empty in
   release builds). The 2.17.0 help text is in `scratch_PA/blastp_help.txt`.
3. `CBlastpApp::Run` (blastp_app.cpp:115): `SetDiagPostLevel(Warning)`, prefix `blastp` (messages read
   `Warning: [blastp] ...`), `RecoverSearchStrategy`, else `CBlastAppArgs::SetOptions` (blast_args.cpp:3585): first
   `m_FormattingArgs->ArchiveFormatRequested` = `ParseFormattingString` (so `-outfmt` errors precede everything),
   create the options handle from `-task`, call `ExtractAlgorithmOptions` of each arg class in the order above
   (subject read here; query and `-out` opened here; hit list and the `Examining 5 or more matches` warning;
   thread warnings; `-ungapped` + composition error), then `Validate()` (`BLAST_ValidateOptions`:
   evalue, threshold, word size, matrix and gap pair; PV range). `CException`s of `Validate` are re-thrown as
   `CInputException`.
4. `x_RunMTBySplitDB` (blastp_app.cpp:181; with `-subject` `CMTArgs` forces one thread, so
   `x_RunMTBySplitQuery` is unreachable): `InitializeSubject` (blast_app_util.cpp:163; `BL2SEQ_LEGACY`),
   `IsIStreamEmpty` (blast_app_util.cpp:845) -> `Query is Empty!` rc 0 before any prolog, `CBlastInput(batch =
   GetQueryBatchSize = 10000, env BATCH_SIZE)`, `CBlastFormat`, `PrintProlog`, per batch
   `GetNextSeqBatch -> CLocalBlast::Run -> PrintOneResultSet`, `PrintEpilog`; everything inside `CATCH_ALL`
   (blast_app_util.hpp:167-270: exception class -> message, exit 1/3/4/6/255).

LOSAT path: `cli.rs:63 try_parse_from` (clap; `is_unported_*_arg` exists for blastn/tblastn/blastx but **not for
blastp**, so every unimplemented blastp option is a generic `unknown option or argument`) -> `main.rs:117
blastp::run(args)?` -> `blast_engine.rs:4022 run`: `resolve()` (args.rs:517), `validate_requested_blastp_support`
(blast_engine.rs:1871), read query then subject with `bio::io::fasta`, build `ReportOutputs` (file created last),
`run_local`. Any error is an `anyhow` error: `Error: msg` plus `Caused by:` chain, exit 1 (never NCBI's classes 3,
6, 255). `NativeError`/`exit_on_native_error` (cli.rs:431-530) exist but blastp does not use them.

The web adapter (`web/adapter/src/run.rs`) parses the same `BlastpArgs` through `try_parse_from` and calls
`run_local_blastp`; `validate(words)` is a no-op for blastp (`_ => {}`, run.rs:98). Anything the session moves
from the clap layer to a post-parse check must also be reachable from `validate` (as `blastn::scoring::check_scoring`
and `tblastx::check_options` are), otherwise the web path loses the NCBI messages. The adapter refuses `-out`
and `-outfmt` itself, so rows about `-out`/`-outfmt` concern the CLI only; the `-outfmt` text parser should still
be a shared function.

## 2. Tables and data the port needs

No large NCBI data tables belong to this range. What the argument layer needs, all small:
- NCBI's integer grammar (`s_StringToInt8` + range, decimal/`0x`; `blastinput/value_parsers.rs:190 ncbi_integer`
  and `:280 ncbi_constrained_integer` already exist and are used by BLASTN) and real grammar (`strtod` /
  `StringToDouble(fDecimalPosixOrLocal)`, `:232 ncbi_double`). `ncbi_double` still rejects hexadecimal reals
  (`0x10`, `0x1p3`) and `1e`; BLASTN keeps these as explicit rejections. For BLASTP the oracle accepts them, see
  decision D2.
- Boolean grammar (`NStr::StringToBool`: `ncbi_environment.rs:109 ncbi_string_to_bool` exists).
- The NCBI message texts quoted in the rows (`BLAST query/options error: ...`, `Command line argument error:
  Argument "x". File is not accessible:  `f'`, `BLAST engine error: Empty CBlastQueryVector`, ...).
- Reusable LOSAT pieces: `blastn/input.rs` (`open_input`, `check_utf8_file_name`, `read_records`),
  `blastn/hsp.rs:78 parse_blastn_output_format`, `blastn/blast_engine/run.rs:5134-5360` (the order of a run, stdin,
  `-out`, `Query is Empty!`, write failure), `report/query_warnings.rs:122 few_matches_warning`,
  `blastinput/ncbi_environment.rs:461 check_ncbi_application_settings`, `blastinput/query_batch.rs`,
  `tblastx/blast_engine/run_impl.rs:1433-1530` (environment checks for `CHUNK_SIZE`, `BL2SEQ_LEGACY`,
  `PRE_FETCH_SEQS_LIMIT`, `BATCH_SIZE`).

## 3. Surprises (things the oracle showed that the sources alone do not)

1. **`-out -` writes a file named `-`** in LOSAT blastp (NCBI: standard output). Row 29.
2. `-evalue 0`, `-0`, `1e-400` are accepted by LOSAT (rc 0, empty report) where NCBI stops with `expect value or
   cutoff score must be greater than zero` (rc 1); negative values are rejected by the clap layer (rc 2). Row 35.
3. NCBI accepts hexadecimal everywhere (`-word_size 0x3`, `-gapopen 0x0B`, `-evalue 0x10`, `-threshold 0x10`,
   hexadecimal-float `-seg "12 0x1p1 2.5"`) and `+inf`, `1e999`, `+nan`, `-nan` for `-evalue`; LOSAT's blastp parsers
   (plain `usize`/`i32`/`f64` parse) accept none of these. Rows 18-20, 36, 48, 51.
4. `-comp_based_stats` has no error path in NCBI: an unknown first character (`4`, `x`, empty, a blank) means mode 0, extra
   characters after the first are ignored except a `u`/`U` second character (unified P). LOSAT rejects all of
   them. Row 75. With the `u` suffix the NCBI output is **not reproducible** in the q1/s6 test (4 runs of `1u` gave
   4 different reports; `2u` 3 different; q3/s_full gave identical output): recommend an explicit rejection of
   unified P (NCBI defect policy). Row 74.
5. `-seg`, `-threshold 0`, `-outfmt abc` and similar options are checked by LOSAT inside clap (exit 2) although NCBI raises
   them **after** parsing as `BLAST query/options error` (exit 1 or 255) and in a fixed order relative to
   the file errors. S08 AUTHORITY F.8 listed these as TD-13. Rows 46, 51, 57, 95.
6. Error order (row 95): NCBI reads and checks `-outfmt`, then the subject (missing -> `Either a BLAST database or
   subject...` rc 1; empty -> `Empty CBlastQueryVector` rc 3), then the query, then creates `-out`, then the
   options, finally `Validate`; warnings printed by an earlier handler stay on stderr when a later one fails. LOSAT:
   clap, `resolve`, engine support check, query, subject, `-out` last.
7. Environment variables that change NCBI blastp output and that LOSAT ignores: `BL2SEQ_LEGACY` (huge change),
   `OLD_FSC` (E-values; also set by `-searchsp`), `ADAPTIVE_CBS`; that stop NCBI: `BATCH_SIZE`,
   `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE`, `PRE_FETCH_SEQS_LIMIT` (non-integer -> rc 255), `BLAST_USAGE_REPORT=bogus`
   (SIGABRT), `ABORT_ON_THROW`; stderr-only: the `DIAG_*` family. `BATCH_SIZE=0` gives `Empty CBlastQueryVector` rc 3.
   Rows 8-13, 107-112. `main.rs` does not call `check_ncbi_application_settings("blastp")`.
8. An empty query file: NCBI prints only `Warning: [blastp] Query is Empty!` (rc 0, no prolog); LOSAT prints a
   complete outfmt 0 report. A query on a pipe is never "empty" for NCBI (`printf '' |` runs with zero queries).
9. Write failure: outfmt 0 `-out /dev/full` -> NCBI `BLAST failed to write output` rc 6, LOSAT rc 1; outfmt 6/7 -> LOSAT
   exits **0 silently** (PD-3 promises a reported error). Rows 102-103.
10. `-help`'s sibling `-h` fails in LOSAT blastp (`unknown option '-help'`), contrary to PD-1. `-version` is unknown.
11. `-num_threads 2147483647`: LOSAT rejects (`Rayon maximum 65535`); NCBI accepts and caps. Row 70.
12. `-remote` with `-subject` really contacts the NCBI servers (15 s in the oracle) and returns a different report.
13. `-searchsp N` (any N, even 0) sets the process environment `OLD_FSC=true` before the search (blast_args.cpp:303-307);
    `-searchsp 0` therefore changes the report (368 bytes against 548 in the test).
14. `-task blastp-fast|blastp-short` are valid in NCBI; LOSAT's parser accepts only `blastp` although `resolve()`
    already knows the three tasks (args.rs:299, 195).
15. NCBI ignores unknown tokens in a custom `-outfmt "6 ..."` list and understands `delim=`; LOSAT errors (out of the
    declared scope, rows 59-60).
16. Batching is invisible for BLASTP reports (BATCH_SIZE=1 and 1000000 equal the default) but a failure in a later batch leaves the
    prolog and the first batch on stdout with rc 3 in NCBI (rows 99, 100).

## 4. Decisions the session must take

- D1. Make `-query` and `-subject` optional in `BlastpArgs` (stdin default / NCBI missing-subject error). The web
  adapter shares `BlastpArgs`; its `parse(words)` and `describe` read the clap arguments, so a change of
  `required` shows in the `describe` JSON (check `web/adapter/src/describe.rs:36,85`).
- D2. Hexadecimal reals and `1e`: port (`ncbi_double` extended; the Rust std parser needs a hex-float reader) or keep
  BLASTN's explicit rejection ("not supported by LOSAT's BLASTP"). Integers in hexadecimal are cheap (`ncbi_integer`).
- D3. `-evalue +inf/+nan/1e999`: needs engine behaviour (infinite cutoff; NaN gives an empty report in the oracle).
- D4. Hidden and toolkit options (`-version-full`, `-logfile`, `-conffile`, `-xmlhelp`, `-dryrun`): proposal is to port
  `-version` and `-dryrun`, reject the others explicitly.
- D5. The `'num_threads and mt_mode' is currently ignored when 'subject' is specified.` warning (printed for
  `-num_threads N>1 -mt_mode 1`) is not literally covered by PD-2's wording; treat it like the `num_threads` warning.
- D6. Unified P (`-comp_based_stats Nu`): reject (non-reproducible in the oracle).
- D7. Search strategy import/export, `-remote`, `-html`: proposed explicit rejection.
- D8. The structure of `blastp::run`: restructure it like `blastn::run` (rows 95-97) so that the order of errors,
  file creation and `Query is Empty!` follow NCBI; this is the prerequisite of most `port:S` rows. Express the
  post-parse option checks in a function callable from the adapter's `validate`.
- D9. Add `is_unported_blastp_arg` (or implement the options) so that the remaining unimplemented options give an
  explicit "not supported by LOSAT's BLASTP" error instead of the generic text (status `divergent` in the table
  becomes `rejected`).

## 5. Statuses

See the reply of the agent for the counts; the `proposal` column uses `port:S|M|L` and `reject:<reason>`.
Rows tagged "PV"/"PS"/"IP" in their notes continue in the ranges PV (options layer and Validate), PS (matrices,
statistics, composition-based statistics) and IP (protein FASTA input); this range records only the argument and
application-layer part of them.
