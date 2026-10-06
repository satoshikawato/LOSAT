# PV notes: BLASTP options layer, task defaults, NCBI option checks

Table: `PV.tsv` (78 rows). Evidence: `scratch_PV/` (files `out_<tag>.{ncbi,losat}.{stdout,stderr,rc}` written by
`scratch_PV/cmp.sh <tag> <blastp args>`, which runs NCBI 2.17.0 and `LOSAT-before` on `q.faa`/`s.faa`
(4 queries, 20 subjects cut from AvCLPV/CoBV); `segrun.sh`/`segrun.out` (21 SEG parameter sets on `q2.faa`/`s2.faa`,
50 x 90 proteins); `cmp_tables.py` (NCBI vs LOSAT gap tables); `epilog.sh` + `epilog_{matrix,ws,task}.txt`
(NCBI epilog matrix/gaps/threshold/window for matrices, word sizes, tasks). Only BLASTP was run; the NCBI functions
named below are shared with TBLASTN (`BlastScoringOptionsValidate`, gap tables, `BLAST_GetSuggested*`), so the
matrix/gap rows apply there too, but I did not run TBLASTN. The machine was shared with other inventory agents.

## Call path followed (NCBI, in order)

1. `CBlastpApp::Run` (`src/app/blast/blastp_app.cpp:114-215`, `SetOptions` at 132): `CBlastAppArgs::SetOptions(args)` is the first
   thing after the diagnostics setup. The subject is opened (`InitializeSubject`) and `Query is Empty!` is decided
   only afterwards (`x_RunMTBySplitDB`:181-215), so every option error of this range precedes a missing/empty input.
2. `CBlastAppArgs::SetOptions` (`blast_args.cpp:3585-3640`): `x_CreateOptionsHandle` -> `CBlastOptionsFactory::CreateTask`
   (`blast_options_handle.cpp:381-402`; `CBlastAdvancedProteinOptionsHandle` over `CBlastProteinOptionsHandle`
   defaults), then each arg class's `ExtractAlgorithmOptions` in the order of `CBlastpAppArgs`'s constructor
   (`blastp_args.cpp:44-117`): program, task, db, std (query/out), **CGenericSearchArgs** (evalue; gap costs from
   `BLAST_GetProteinGapExistenceExtendParams` when `-matrix` is given and each explicit side overrides; word size >4 =>
   compressed lookup, threshold 19.3/21/20.25 for 5/6/7; xdrops; max_hsps), **CFilteringArgs** (soft masking flag,
   `-seg` tokens), **CMatrixNameArg**, **CWordThresholdArg** (explicit threshold, else suggested-by-matrix only when
   `(int)threshold == 11`), HSP filtering, **CWindowSizeArg** (explicit, else suggested-by-matrix), query options,
   formatting (`-outfmt`, "Examining 5 or more matches" warning), MT, **CGappedArgs** (`-ungapped`), remote,
   **CCompositionBasedStatsArgs** (`s_SetCompositionBasedStats`: first-character switch, ungapped error, `u` suffix,
   SW flag), debug. Then `m_OptsHandle->Validate()`; a `CBlastException` is rethrown as `CInputException`
   (`BLAST query/options error: <msg>` + `Please refer to the BLAST+ user manual.`, exit 1).
3. `BLAST_ValidateOptions` (`blast_options.c:1749-1812`) order: `BlastExtensionOptionsValidate`,
   `BlastScoringOptionsValidate` (matrix name via `Blast_KarlinBlkGappedLoadFromTables`, status 1 -> `BLAST_PrintMatrixMessage`,
   status 2 -> `BLAST_PrintAllowedValues`; skipped when `-ungapped`), `LookupTableOptionsValidate` (threshold > 0, then
   word size: <=0, >7), `BlastInitialWordOptionsValidate`, `BlastHitSavingOptionsValidate` (hitlist, evalue <= 0),
   `s_BlastExtensionScoringOptionsValidate`, then the IDENTITY word-size check.
4. Not option layer but option-dependent and reached before the search: `CSetupFactory::CreateScoreBlock`
   (`setup_factory.cpp:144-156`, IDENTITY resets composition statistics to 0 with a warning per query),
   `BlastSetup_ScoreBlkInit` -> `Blast_ScoreBlkMatrixInit` (matrix name uppercased in the score block only;
   unknown name with `-ungapped` gives `BLAST engine error: Error: Unknown error code -1 ` x queries, exit 3),
   `Blast_ScoreBlkKbpUngappedCalc` (fully masked query warning).

## What the port needs (data and code, with locations and sizes)

* Gap/Karlin tables: `src/algo/blast/core/blast_stat.c` lines 183-198 (BLOSUM45, 14 rows), 219-236 (BLOSUM50, 16),
  258-271 (BLOSUM62, 12), 290-301 (BLOSUM80, 10), 317-326 (BLOSUM90, 8), 340-357 (PAM250, 16), 379-391 (PAM30, 11),
  409-419 (PAM70, 9), 577-581 (IDENTITY, 2); 11 columns (open, extend, INT2_MAX, lambda, K, H, alpha, beta,
  C, a/Alpha, Sigma) plus the `*_prefs` arrays (BEST row). LOSAT (`stats/tables.rs`, duplicated in
  `core/blast_stat/lookup_tables.rs`) keeps 5 columns; PAM30 lacks rows 15/3, 14/2, 14/1, 13/3 and has alpha 0.34
  (NCBI 0.35) for 10/1; PAM70 lacks 11/2, 12/3 and has ungapped beta -0.9 (NCBI -0.7). The other six tables are equal.
  The sentinel row uses INT2_MAX=32767, LOSAT uses `i32::MAX`.
* Substitution matrices: `src/util/tables/sm_*.c`; LOSAT already has the 8 non-identity ones in `utils/matrix.rs`
  (compared number by number: equal). `sm_identity.c` (625 numbers, 9 diagonal, -5 elsewhere) is missing.
* Composition adjustment data for the other matrices: `composition_adjustment/matrix_frequency_data.c` (1418 lines,
  8 matrices) and `core/matrix_freq_ratios.c` (1753 lines); LOSAT has `blosum62_start_freq_ratios.inc.rs` only.
  `compo_mode_condition.c` (256), `unified_pvalues.c` (260), `smith_waterman.c` (615) for modes 3/1, `u` suffix and
  `-use_sw_tback` (details: range PS).
* Messages to reproduce verbatim: matrix list order BLOSUM80, BLOSUM62, BLOSUM50, BLOSUM45, PAM250, BLOSUM90,
  PAM30, PAM70, IDENTITY (each `NAME \n`, note the space), the "supported values are:" row lists, 1024-byte (matrix
  message) and 2048-byte (gap message) buffers with `snprintf` truncation, `Non-zero threshold required`,
  `Word-size must be less than 8 for a tblastn, blastp or blastx search`, `expect value or cutoff score must be
  greater than zero`, `Word size larger than 5 is not supported for the identity scoring matrix`,
  `Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a
  gapped search` (sic), `Invalid number of arguments to filtering option`, `Invalid input for filtering parameters`.
* Option errors must go out as `NativeError{exit:1, "BLAST query/options error: <msg>\nPlease refer to the BLAST+ user
  manual.\n"}` (model: `algorithm/blastn/scoring.rs:130-160`), engine-time matrix failures as
  `BLAST engine error: ` + `Error: Unknown error code -1 ` repeated per query (model: `karlin_error`).

## Surprises worth knowing

* LOSAT's `-task blastp-fast`/`blastp-short` exist in `args.rs` but are unreachable (`value_parser ["blastp"]`).
  `resolve()` is wrong for blastp-fast: NCBI's threshold is 20 (`BLAST_WORD_THRESHOLD_BLASTP_FAST`), also with
  `-word_size 3`; LOSAT gives 19.3 and 11. NCBI `-task blastp-fast` output equals `-word_size 5 -threshold 20`
  on the probe, and LOSAT already reproduces that run.
* Many LOSAT clap checks are stricter than NCBI's declared ranges and therefore are **not** covered by
  PD-LOSAT-CLI-NONSEARCH-DIFFERENCES item 1 (NCBI raises them after parsing as option errors, exit 1): `-evalue`
  sign/zero, `-threshold 0`, negative `-gapopen/-gapextend`, every `-seg` token error, `-comp_based_stats` string
  (NCBI has no constraint: unknown first character means mode 0, only `[1]` is read for `u`, case-insensitively,
  trailing characters ignored). Declared ranges (item 1 applies): `-word_size >= 2`, `-threshold >= 0`,
  `-window_size >= 0`, `-task` set, `-max_target_seqs >= 1`, `-culling_limit >= 0`.
* NCBI accepts `+inf`, `1e400`, `+nan` for `-evalue` and `+inf`/`1e999`/hex floats in `-seg`; bare `inf`/`nan` fail the
  first-character test. LOSAT (CLI v2 rule) rejects all non-finite values.
* `(Int4)threshold` of an out-of-range double: x86-64 gives INT_MIN (NCBI: every word is a neighbour, same output as
  `-threshold 1`); Rust `as i32` saturates to INT_MAX. Compressed lookup multiplies by 100 first and then **aborts**
  (heap corruption) for thresholds >= 21474836.48 with word size >= 5.
* `-matrix` keeps the typed spelling in NCBI (epilog `Matrix: pam30`, messages), comparisons are `strcasecmp`
  except `strcmp(matrix, "BLOSUM62")` for chaining (blastp-fast only).
* `-ungapped` skips the gap-pair check, keeps sum statistics off, and accepts unknown matrices until the engine
  fails (exit 3). `-ungapped` with any enabled composition mode is an option error, but `-comp_based_stats x`
  silently means mode 0 and therefore passes.
* `-use_sw_tback -ungapped -comp_based_stats 0` crashes NCBI (segmentation fault, deterministic here).
* IDENTITY is a valid matrix with special cases (threshold 27, gap 15/2, composition statistics reset with a
  warning per query, word size > 5 rejected, no compressed lookup).
* The extra HSPs LOSAT reports for huge `-evalue` (1e300, 1e308; the RP hand-off for 1000) are raw-score-0 HSPs
  (bit score 4.6): cause not isolated (zero-score drop in `blast_seqalign.cpp:672-674` vs cutoff clamp).
* A fully SEG-masked query makes NCBI print the Karlin-Altschul warning per query; LOSAT prints nothing.
* BLASTX in this tree (`algorithm/blastx/args.rs`) only accepts BLOSUM62 and composition 0/2 (also string-checked),
  so it is no model for the matrix/composition port; BLASTN's option chain (`blastn/scoring.rs`,
  `blastn/args.rs:309` for `-dust`) is the model for placement and message style.

## Decisions the session must take

1. Move the option checks to NCBI's order and message style (rows 2, 3), including `-seg`, `-comp_based_stats`,
   `-evalue`, `-threshold`, `-gapopen/-gapextend` out of clap (rows 31, 45, 52-54, 64, 65, 72-74). Needs a rule for
   which clap checks stay (declared NCBI ranges only).
2. Non-finite numbers: port NCBI's `strtod` behaviour (rows 48, 54, 74) or keep rejecting with an explicit
   "not supported by LOSAT" message.
3. NCBI crashes: row 49 (compressed threshold overflow, abort) and row 69 (`-use_sw_tback -ungapped`, SIGSEGV):
   approved exception with the intended result, deterministic reproduction, or explicit rejection.
4. Row 27 (custom matrix files through `-ungapped`) and the `.ncbirc`/`BLASTMAT` route: reject or ignore.
5. Order of work suggested by the table: option-layer ports (S rows) first; then mode 0 and `-ungapped` (rows 60, 67),
   matrices (row 11 pipeline, then 13-22), word sizes 2/4/6/7 (with range LK), unified P / modes 1, 3 / SW (range PS).
