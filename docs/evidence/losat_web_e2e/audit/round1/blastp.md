# BLASTP application flow and option checks - independent audit

Scope: LOSAT blastp (audited commit source in $D/src, binary $D/LOSAT) vs NCBI blastp 2.17.0 (/home/kawato/micromamba/bin/blastp).
Method: work/blastp/cmp.sh runs both programs with the same argv (-query/-subject from tests/fasta/outfmt0 unless noted), compares stdout, stderr and exit status; batch.sh runs lists (about 1,000 argv in total). All runs used `timeout 60`, at most 2 processes of mine at a time.
Parser-syntax errors (clap text, exit 2) are the approved exception PD-LOSAT-CLI-NONSEARCH-DIFFERENCES and are not findings. Severity: high = silent output difference or NCBI crash/behaviour not reproduced and not rejected; medium = wrong message/exit/check order; low = wording/reference/doc.

Totals: high 5 (BP-1..BP-5), medium 2 (BP-6, BP-7), low 6 (BP-8..BP-13). All CONFIRMED by running unless marked SUSPECTED.

## Findings

### BP-1 (high, CONFIRMED) All-empty subject set: outfmt 0 prints "Effective search space used: <n>"; NCBI prints 0
- LOSAT: report epilog / search-space statistics for a subject set whose records all have zero length (algorithm/blastp/blast_engine.rs report path); input not an option but part of the app flow (empty-subject handling, AUTHORITY §A).
- NCBI: database length 0 gives effective search space 0 for every query (blast_format.cpp epilog).
- Repro (F=$D/src/LOSAT/tests/fasta/outfmt0; printf '>x\n' > hdronly.faa):
  `blastp -query $F/e2e_protein_query.faa -subject hdronly.faa` and `$D/LOSAT blastp ...same...`: both exit 0, same stderr ("Warning: [blastp] Subject_1 x: Subject sequence contains no data"); stdout differs in 32 diff lines: NCBI "Effective search space used: 0", LOSAT "Effective search space used: 480" (386, 408, 436, 849, ... one per query). outfmt 6 and 7 identical (empty).

### BP-2 (high, CONFIRMED) -max_target_seqs in [2^30-25, 2^31-25]: NCBI SIGSEGV (exit 139), LOSAT runs and prints a full report (exit 0)
- LOSAT: algorithm/blastp/args.rs `max_target_seqs: blastn_count` (no upper check); algorithm/blastp/hsp.rs:300-326 `get_prelim_hitlist_size` (usize, saturating_mul/saturating_add, no Int4 wrap) used at blast_engine.rs:5059. BLASTN/TBLASTN use blastn/hsp.rs:508 (wrapping i32) and reject the crashing range (DW-15); BLASTP does neither.
- NCBI: blast_hits.c:43-70 `GetPrelimHitlistSize` in Int4: with cbs (default for blastp) and h > 500 the size is 2h+50 wrapped; from h = 2^30-25 it is negative (crash) until h = 2^31-25 (size 0, crash).
- Repro (Q=$F/e2e_protein_query.faa, S=$F/e2e_protein_subject.faa):
  `blastp -query $Q -subject $S -max_target_seqs 1073741799 -outfmt 6` -> "timeout: the monitored command dumped core", exit 139, no output. LOSAT: exit 0, report.
  Same for 1073741823, 1073741824, 2147483597, 2147483598, 2147483622, 2147483623 (outfmt 0 and 6). 1073741798: identical, exit 0.
- Policy: PD-LOSAT-NCBI-DEFECTS (AGENTS.md rule 3) needs an approved exception (with a checkable valid result) or an explicit rejection; neither exists for BLASTP.

### BP-3 (high, CONFIRMED) -max_target_seqs in [2^31-24, 2^31-1]: NCBI keeps only (2h+50) mod 2^32 = 2..48 subjects per query (deterministic), LOSAT keeps all
- LOSAT: same cause as BP-2 (hsp.rs:300 does not wrap). DW-15 says deterministic NCBI results must be reproduced.
- Repro (Q=$F/e2e_protein_query.faa, S=$F/e2e_many_subject.faa (300 subjects)), `-outfmt 6`; rows for query BDT62569.1 (NCBI vs LOSAT):
  -max_target_seqs 2147483624: 2 vs 306 ; 2147483630: 14 vs 306 ; 2147483640: 34 vs 306 ; 2147483646: 46 vs 306 ; 2147483647: 48 vs 306 (total rows 61 vs 319; outfmt 0 bytes 60666 vs 327860).
  With the small subject file (e2e_protein_subject.faa) 2147483647 is identical because fewer than 48 subjects hit (the difference only shows with enough hits).
  AUTHORITY §K lists "-max_target_seqs >= 1" as supported without limit.

### BP-4 (high, CONFIRMED) -evalue +inf / 1e999 (any infinite value): NCBI SIGSEGV when enough hits exist; LOSAT prints a full report
- LOSAT: accepted by blastp_real and validate_protein_options (only evalue <= 0 is checked); AUTHORITY §K: "-evalue ... including +inf and 1e999" supported.
- NCBI: crash cause not isolated (SUSPECTED: cutoff-score computation from an infinite e-value). Deterministic crash, no output.
- Repro: Q=$F/e2e_protein_query.faa S=$F/e2e_many_subject.faa: `blastp -query $Q -subject $S -evalue +inf -outfmt 6` -> exit 139, empty stdout; LOSAT exit 0, 1198734 bytes. -evalue 1e10, 1e100, 1e308 on the same input: identical.
  With the compressed lookup the small input already crashes NCBI: `-word_size 5 -evalue +inf -outfmt 6` (and 1e999, outfmt 0/7, `-task blastp-fast -evalue +inf`): NCBI exit 139 (outfmt 0: 17 bytes then crash); LOSAT exit 0, 2541 bytes (outfmt 6) / 30265 bytes (outfmt 0). With word size 3 and the small input both exit 0 identically.

### BP-5 (high, CONFIRMED) -outfmt "6 frames" / "6 sframe" (and 7): LOSAT prints 1 and 1/1, NCBI prints 0 and 0/0 when no other field makes NCBI compute the alignment
- LOSAT: blast_engine.rs:2295 `blastp_frame_value`, 2362-2363, 2406-2408, 2494-2496 (writer prints "1", "1", "1/1" unconditionally).
- NCBI: objtools/align_format/tabular.cpp:900-910 (the SetFields block runs only if one of qstart/qend/sstart/send/length/gaps/gapopen/qseq/sseq/nident(>0)/positive/mismatch/ppos/pident/qframe/btop/sstrand is requested), 1071-1094 (frames set there), 106 (m_QueryFrame = m_SubjectFrame = 0 otherwise).
- Repro: `blastp -query $Q -subject $S -outfmt "6 sframe"`: every row "0"; LOSAT "1". `"6 frames"`: "0/0" vs "1/1". Also differ: `"6 sseqid sframe"`, `"6 evalue frames"`, `"6 qseqid sseqid frames"`, `"6 stitle frames"`, `"7 frames"`, `"7 sseqid sframe"`, with qlen/slen/score/bitscore/qacc/stitle/evalue. Identical: `"6 qframe sframe"`, `"6 ppos frames"`, `"6 std frames"` (any list that contains a triggering field).

### BP-6 (medium, CONFIRMED) -seg "W LOCUT HICUT": malformed locut/hicut get LOSAT's "not supported" text, NCBI says "Invalid input for filtering parameters"
- LOSAT: blastinput/app.rs `parse_seg_option` (~l.407-428): any token that Rust `f64::from_str` rejects but whose first char is a digit/./+/- becomes NcbiDoubleError::Unsupported -> "the SEG locut or hicut "..." (not a finite decimal number) is not supported by LOSAT's BLASTP".
- NCBI: blast_args.cpp:396-428 (StringToDouble, CStringException -> "Invalid input for filtering parameters"), exit 1 (same exit code, different text).
- Repro: `-seg "12 1e 2.5"`, `"12 2.2 2.5e"`, `"12 2.2 2.5x"`, `"12 . 2.5"`, `"12 +. 2.5"`, `"12 -. 2.5"`, `"12 1e+ 2.5"`, `"12 1,5 2.5"`, `"12 2.5.5 2.5"`, `"12 --5 2.5"`, `"12 1_0 2.5"`, `"12 .e1 2.5"`: NCBI "BLAST query/options error: Invalid input for filtering parameters"; LOSAT the "not supported" text (both exit 1).
  The explicit rejection is only right for values NCBI accepts and LOSAT does not read: `+inf`, `-inf`, `1e400`, `0x10`, `0x1p3`, `-nan`, `+nan(5)` (NCBI exit 0 and runs; LOSAT exit 1 with the "not supported" text - fine).

### BP-7 (medium, CONFIRMED) LOSAT's own rejections that run before NCBI's later checks hide NCBI's error text and exit code
- LOSAT: blast_engine.rs run(): `report_format` (unsupported -outfmt 1-5, 8-12, 15, 16, 18, 20; 13/14 with -out), `validate_cli_outfmt` (custom fields NCBI writes and LOSAT does not, delim=), `validate_threads` all execute before the subject is read and before the option handlers.
  NCBI order: outfmt parse (blast_args.cpp:2801-2851), subject read (2553-2562), query/out open, handlers (3631-3639), Validate, "Query is Empty!".
- Repro (all exit 1 for both unless noted):
  `-outfmt 5 -matrix FOO`, `-outfmt 5 -threshold 0`, `-outfmt 5 -ungapped`, `-outfmt 5 -seg "1 2"`: NCBI prints its option error, LOSAT "output format 5 is not supported by LOSAT's BLASTP".
  `-subject empty.faa -outfmt 5` (and blank.faa): NCBI exit 3 "BLAST engine error: Empty CBlastQueryVector", LOSAT exit 1.
  `-query empty.faa -outfmt 5`: NCBI exit 0 ("Query is Empty!"), LOSAT exit 1.
  `-num_threads 100000 -matrix FOO`: NCBI "BLAST query/options error: FOO is not a supported matrix..." (after its thread-clamp warning), LOSAT "requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT" (both exit 1).
  Contract (3) holds (message has "not supported by LOSAT's BLASTP"); contract (2) is not met where NCBI itself fails first.
  Checked equal: the order among NCBI's own checks (see Checked OK).

### BP-8 (low, CONFIRMED) NCBI toolkit argument `-dryrun` (and toolkit args in option-value position)
- LOSAT: cli.rs `is_ncbi_toolkit_arg` (list lacks `dryrun`).
- NCBI: corelib/ncbiargs.cpp:88 `s_ArgDryRun`, ncbiapp.cpp:1000-1001 (stripped from argv anywhere, runs DryRun(): no output, no search).
- Repro: `blastp -query $Q -subject $S -dryrun` NCBI exit 0, empty stdout/stderr; LOSAT exit 2 "unknown option or argument '-dryrun'" (no "not supported by LOSAT's BLASTP").
  ncbiapp.cpp:955-1001 scans every argv word, including values: `-evalue -version` / `-task -version` NCBI prints the version, exit 0; LOSAT clap error exit 2. `-out -version` NCBI prints the version, LOSAT runs the search and creates a file named "-version" (silent). `-out -dryrun`: NCBI USAGE exit 1 (-out lost its value), LOSAT writes the report to a file "-dryrun", exit 0.
  `-outfmt 6 --`: NCBI exit 0 (empty trailing `--` accepted), LOSAT exit 2 "unknown option or argument '--'".

### BP-9 (low, CONFIRMED) Rejections whose text lacks "not supported by LOSAT's BLASTP"
- `-num_threads 100000`: "requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT" (utils/threading.rs:80; also "unsupported num_threads=N: this build does not support parallel search").
- Environment and registry rejections from blastinput/ncbi_environment.rs end with "is not supported by LOSAT" (e.g. `DIAG_POST_LEVEL=Info`, `NCBI_CONFIG__BLAST__X=1`, `.ncbirc` `[BLAST] BATCH_SIZE=1`); the .ncbirc message also reads "which LOSAT does not know to change no output of NCBI BLAST+" (garbled), and NCBI itself ignores `[BLAST] BATCH_SIZE` in a registry file (NCBI exit 0 with BATCH_SIZE=0 there, while the environment variable gives exit 3), so that rejection is conservative.

### BP-10 (low, CONFIRMED) `LOSAT_STARTUP_TRACE=1` adds "[startup] enter main" / "[startup] after clap parse" to stderr (main.rs) - a LOSAT-only stderr difference, and a feature NCBI does not have (AGENTS.md rule 5).

### BP-11 (low, CONFIRMED) `--help`: NCBI exit 1 (unknown argument), LOSAT prints help and exits 0 (cli.rs maps "--help" to clap's help). `-help` is the approved exception; `--help` is not NCBI syntax.

### BP-12 (low, CONFIRMED) NCBI file:line references that do not contain the cited snippet (checked with work/blastp/vref2.py: 69 snippet blocks in args.rs, app.rs, value_parsers.rs, protein_options.rs, cli.rs, main.rs, plus 395 in blast_engine.rs, none of the blast_engine.rs ones is session-added; spot-checked 30+ inline citations by hand)
| LOSAT | cited | actual NCBI | session-added |
|---|---|---|---|
| stats/protein_options.rs:390 | blast_options.c:913-943 | gapped matrix/gap check is at 910-940 (`if (options->gapped_calculation ...` line 910) | yes |
| stats/protein_options.rs:492 | blast_options.c:1515-1520 | `expect value or cutoff score must be greater than zero` is line 1521 (check at 1518-1522) | yes |
| blastinput/app.rs:638 | blast_input.hpp:258-259 | those lines are an enum-name switch; `TSeqPos m_BatchSize` is at 364, 445-451 | yes |
| algorithm/blastp/blast_engine.rs:7011 | blast_format.cpp:2266 | `options.GetMatrixName()` is line 2267 (2266 is `else {`) | yes |
| algorithm/blastp/args.rs:91 | blast_args.cpp:845-886 | `switch (comp_stat_string[0])` is line 838 | no |
| blastinput/value_parsers.rs:550 | cmdline_flags.cpp:46-50 | `kArgSubject` is line 51 | no |
| cli.rs:181 | cmdline_flags.cpp:46-94 | kArgWordSize is 107, kArgCompBasedStats 143 | no |
| cli.rs:209 | tblastn_args.cpp:64-129 | first `m_BlastDbArgs.Reset` is line 63 | no |
Also: protein_options.rs `validate_protein_options` `word_size <= 0` ("Word-size must be greater than zero") is an unreachable port (parser enforces >= 2) shown only as `...` in the snippet.

### BP-13 (low, SUSPECTED by code reading) Web adapter vs CLI
- web/adapter/src/run.rs:92 `validate` calls `blastp::blast_engine::check_options` = `args.resolve()` (the same `BlastpArgs::check_options` the CLI calls, with a sink for warnings) + `validate_requested_blastp_support`: the option verdicts (NCBI errors and LOSAT "not supported" rejections) are the same as the CLI's; it reuses the same clap parser (`LOSAT::cli::try_parse_from`), so parse-level rejections are identical. The adapter forbids -out/-outfmt itself.
- Not checked by `validate` (nor by the CLI's check_options): `validate_threads` (CLI `run` calls it first, blast_engine.rs:4312; the adapter run path reaches it in `with_search_pool`), so `-num_threads 100000` or `-num_threads 2` in a build without threads passes validate and fails at run.
- `run_local` (adapter/web path, blast_engine.rs:4537) does not call `check_protein_input_of`/`check_protein_residues_of` (CLI `search_cli`, blast_engine.rs:4344 and 4481 reject residues other than letters and `*`); records the host passes are not checked here. I could not execute the wasm adapter, so the behaviour for such records is unverified.

## Checked OK (identical stdout, stderr, exit status to NCBI unless noted)
Order of NCBI's own checks, several wrong options at once: `-matrix FOO` with `-outfmt 99/abc/17`, `-threshold 0`, `-word_size 8`, `-evalue 0`, `-seg "1 2"`, `-seg "x 2.2 2.5"`, `-comp_based_stats 2 -ungapped`; `-threshold 0` with `-word_size 8`/`-evalue 0`/`-seg`; `-seg bad` with `-outfmt 17/99`, `-ungapped -comp_based_stats 2`; `-outfmt 17/19/21/13/14` with `-ungapped`/`-matrix FOO`/`-threshold 0`; `-gapopen 10/5 -gapextend 5` with matrix/threshold/word_size/evalue; `-task blastp-short/fast` with the same; `-use_sw_tback` with `-matrix FOO`/`-threshold 0`; unreadable query, missing/empty subject, `-out` to a missing directory/directory/`/proc/x` with each of the above; Query is Empty with each of -outfmt 13/14/17/19/21/99.
Values (both run identically, or both reject with the same text and exit code): -matrix blosum62/Blosum62/BLOSUM62/-matrix=BLOSUM62; unknown matrix FOO, "", " BLOSUM62"; BLOSUM50 with bad gap pairs; IDENTITY with 5/2 and word size 6; negative/huge/2^31-1 gap costs (matrix error text and the allowed-values text); -threshold 0, -0, 0.0, 1e-400 (NCBI error), 0.1, 0.5, 0.999, 1, 5, 10.99, 11.5, 11.99, 12, 15, 24, 100, 1e-5, 1e1, 1e9, 1e10, 1e18, 1e300, 2147483647, 2147483648, 3e9, 4294967296, 1e400, +inf, +infinity, .5e1; -word_size 3, 5, 0x5, 8, 9, 100, 0x8, 2147483647, with -threshold 0.5/1/5/10/12/20/20.5/22/24/25/30/100/5000/2e7/2.1e7/2.14e7/2.147e7/2.1474e7; -window_size 0, 1, 2, 3, 15, 40, 1000, 2147483647, 0x10; -evalue -1, 0, -0, +0, 0.0, 1e-400, 1e-999, 2e-324 (error), 4.9e-324, 1e-320, 1e-300, 1e-100, 1e-30, 1e5, 1000, 1e10, 1e100, 1e308, 1e999, +inf, +infinity, +INF, -nan, +nan, +nan(1), +nan(abc), -nan(1), +nan(), +NaN, +NAN(x_1), .5, +5, 5., 5.e1, .5e1, 1E5, 1e+5 (BP-4 excepted: +inf crash cases); -inf/-infinity NCBI error text; -max_target_seqs 1-5, 10, 100, 1000000, 536870912, 600000000, 1000000000, 1073741798, 2147483647 (small input), 0x3, 0x7fffffff; -max_hsps 1, 2, 3, 100, 0x1, 2147483647; -seg no/yes/"12 2.2 2.5"/"10 1.8 2.1"/"+12 2.2 2.5"/"-1 ..."/"0 ..."/"2 ..."/"100000 ..."/"12 -1 2.5"/"12 .5 2.5"/"12 2. 2.5"/"12 1e1 2.5"/"12 1e-400 2.5"/"12 5 6"/"12 2.5 2.2", errors for ""/" "/"12"/"12 2.2"/"12  2.2 2.5"/" 12 2.2 2.5"/"12 2.2 2.5 "/"12 inf 2.5"/"x 2.2 2.5"/"0x10 2.2 2.5"/"2147483648 2.2 2.5"/YES/No/true/T; -comp_based_stats 2, 2x, 2xyz, 22, 2abc, D, d, T, t, Dx, Tx, true, TRUE, "2 " (all mode 2, identical); all other spellings (0, 1, 3, F, f, x, garbage, "", " 2", -1, 0x2, +2, 4, u, Fx and the u-suffixed 0u/1u/2u/2U/3u/Du/tu/2uu) are explicit rejections: "-comp_based_stats 0|1|3 is not supported by LOSAT's BLASTP" or "unified P-values ... not supported by LOSAT's BLASTP" (exit 1), and `-ungapped` with any cbs reproduces NCBI's "Composition-adjusted searched are not supported with an ungapped search..." for non-zero modes first.
-outfmt: 0, 6, 7, " 6", "6 ", "  6  ", +6, 06, "0 std", "0 qseqid", "0 foo", "0 delim=, qseqid", "6 std", "6 foo", "6 -qseqid", "6 delim= qseqid", all 50 NCBI field names singly (6 and 7; 28 identical, 20 rejected with "the output field X is not supported by LOSAT's BLASTP", except frames/sframe = BP-5), pair lists with sseqid/qlen/slen/evalue/stitle/score/bitscore/qacc x frames/sframe/qframe/nident/mismatch/pident/positive/ppos/gaps/gapopen/length/btop/qseq/sseq/qstart/qend/sstart/send, "6 qseqid sseqid pident length evalue bitscore score nident positive gaps ppos qframe sframe qlen slen", "7 ..." same; errors: "" / " " / abc / 0x6 / 6.0 / 1e1 / 6x / "6x std" / 2147483648 / 4294967302 / 99999999999 ("'X' is not a valid output format", exit 1), -1 / 22 / 23 / 99 ("Error: Formatting choice is out of range", exit 255), "6 delim" / "0 delim" (Delimiter format is invalid), 17 / 19 / 21 / "17 foo" ("... only applicable to blastn/igblastn/magicblast"), 13 and 14 to stdout ("Please provide a file name for outfmt 13/14."); LOSAT rejections (exit 1, text "output format N is not supported by LOSAT's BLASTP" for 1-5, 8-12, 15, 16, 18, 20, "a custom delimiter (delim=) ... not supported by LOSAT's BLASTP" for every delim= spelling).
Large-input runs (1 query x 300 subjects and 8 queries x 300 subjects, outfmt 0/6/7, 40+ combinations): defaults, -max_target_seqs 1/2/3/10, -max_hsps 1/2/3, -evalue 1e-300..1e308, -seg yes / three-value, -window_size 0/1/3/15/1000/2147483647, -threshold 0.001..+inf, -word_size 5 with thresholds, -task blastp-fast (+ seg, window_size 0, threshold, max_hsps), -num_threads 4 (only stderr differs: NCBI's two thread warnings, approved), custom field lists, -matrix BLOSUM62 -gapopen 11 -gapextend 1, -comp_based_stats 2/D/t: all identical.
Compressed-lookup overflow rejection: `-word_size 5 -threshold 21475000` / `2.14749e7` / `1e9` / `+inf` and blastp-fast: LOSAT "... not supported by LOSAT's BLASTP: its neighboring words fill NCBI BLAST+'s 1024 overflow banks ..." (NCBI hangs/times out >60 s); up to 2.1474e7 identical.
Unsupported-option rejections: all 43 names of is_unported_blastp_arg (each run as `-NAME x`) plus -h, -help-full, -xmlhelp, -version, -version-full(-xml/-json), -logfile, -conffile give "the NCBI BLAST+/C++ Toolkit option -X is not supported by LOSAT's BLASTP" (exit 2); blastn/tblastn/tblastx-only options and -bogus: "unknown option" (NCBI USAGE exit 1; approved parser exception). NCBI -help option list (63) is fully covered by LOSAT options + rejections (no NCBI blastp option missing).
App flow: -query/-subject from stdin (`-`), `-subject -` with empty stdin (exit 3 both), `-out -`, empty query via stdin, `-query -subject -` (LOSAT explicit rejection, D10); file name length 200-300, empty names; duplicate flags and missing values (approved parser exceptions); `-evalue=1`, `-outfmt=6`, `-matrix=BLOSUM62`.
Environment: BATCH_SIZE=0 (exit 3 both), 1, 100, -1; BATCH_SIZE=abc/""/0x10/4294967297 (explicit rejection), PRE_FETCH_SEQS_LIMIT 0/5, CHUNK_SIZE=3, OVERLAP_CHUNK_SIZE=5, ADAPTIVE_CBS, BLAST_USAGE_REPORT, BLASTDB, NCBI_DONT_USE_LOCAL_CONFIG identical; OLD_FSC, BL2SEQ_LEGACY, DIAG_*, NCBI_CONFIG__*, bad CHUNK_SIZE/PRE_FETCH values: explicit rejection with the "not supported by LOSAT" text (see BP-9).
NCBI file:line spot-check: about 60 snippet blocks verified programmatically in the six audited files, 8 wrong references listed in BP-12; manually confirmed correct: blast_args.cpp:317-318, 332-349, 475-488, 586-600, 879-883, 2657-2660, 2745-2748, 2910-2927, 2960-2977 (the 2976 warning), 3158-3163; blast_options_handle.cpp:395-399; blast_options.c:46-48, 62-64, 1301-1308, 1303-1395; blast_stat.c:3068-3070; ncbistr.cpp:1332-1340; blast_advprot_options.cpp:57; split_query_aux_priv.cpp:73; blast_app_util.hpp:252-255; blast_seg.c:2247-2250; tabular.cpp:1100-1108; blast_kappa.c:332-342; showalign.cpp:3595-3598; seqsrc_multiseq.cpp:226.
