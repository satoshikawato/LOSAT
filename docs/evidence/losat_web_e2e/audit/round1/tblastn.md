# TBLASTN application flow and option checks - audit notes (in progress)

Environment: D=/home/kawato/.cache/losat-web-gui-target/s08p/audit, NCBI tblastn 2.17.0+, LOSAT binary $D/LOSAT.
Harness scripts: $D/work/tblastn/{cmp2.sh,t4.sh,pairs.py,pairs2.py,sweep.py,verify_refs.py}. Generated inputs: $D/work/tblastn/in/gen_q.faa, gen_q_lc.faa, gen_s.fna (reverse-translated, mutated query fragments on both strands).

## Findings

### TN-1 (high, CONFIRMED) -xdrop_gap_final / -xdrop_gap values below ~16 bits with -comp_based_stats 0 change HSP e-values (silent difference)
- LOSAT: LOSAT/src/algorithm/tblastn/stage_d_pipeline.rs:1167-1181 (finish_local_mode0_no_sum_stats, `final_xdrop`), search_gapped.rs:1416-1454 (`stat_length` from last `translated_length`) - suspected area, not root-caused.
- NCBI: core/blast_parameters.c:455-463 (gap_x_dropoff / gap_x_dropoff_final), core/blast_traceback.c:425-433 (stat_length), :717-719 (Blast_HSPListGetEvalues).
- With -comp_based_stats 0 (mode 0 path) and an effective final X-drop (MAX(final, prelim) bits) of 15 bits or less, or an extreme one (>= 1e9 bits, +inf), LOSAT reports the same HSPs as NCBI but a few HSPs get a different E-value (LOSAT smaller by a constant factor ~0.37-0.40, as if the subject stat_length differs). Stdout differs, exit 0 both, no message.
- Default final X-drop 25 and any final >= 16 bits (and prelim <= final): identical. -comp_based_stats 2 (default): identical in all cases tried.
- Repro 1 (fixtures only, X-drop small):
  `cd $D/src/LOSAT/tests/fasta/outfmt0; Q=e2e_protein_query.faa S=e2e_many_subject.fna`
  `/home/kawato/micromamba/bin/tblastn -query $Q -subject $S -comp_based_stats 0 -xdrop_gap 5 -xdrop_gap_final 5 -outfmt 6 > n.txt; $D/LOSAT tblastn -query $Q -subject $S -comp_based_stats 0 -xdrop_gap 5 -xdrop_gap_final 5 -outfmt 6 > l.txt; diff n.txt l.txt`
  -> 4 rows differ in column 11 only, e.g. `BDT62562.1 n209 ... NCBI 3.0 vs LOSAT 1.2`, `BDT62620.1 n36 NCBI 0.91 vs LOSAT 0.36`, `n138 NCBI 1.8 vs LOSAT 0.65`, `n145 NCBI 0.16 vs LOSAT 0.059` (one row also moves position because of the sort by E-value).
- Repro 2 (fixtures only, huge final X-drop): same files, `-comp_based_stats 0 -xdrop_gap_final 1e10 -outfmt 6` (also 1e9, 4e9, 1e300, +inf): identical 4 rows differ. AUTHORITY.md section F claims `-xdrop_gap_final 1e9` was fixed; it holds for -comp_based_stats 2 only.
- Repro 3 (generated inputs): `-query $D/work/tblastn/in/gen_q.faa -subject $D/work/tblastn/in/gen_s.fna -comp_based_stats 0 -xdrop_gap_final 10 -outfmt 6` -> gq1/gs22 NCBI 2.5, LOSAT 0.93 (identical for -xdrop_gap_final 25).
- -sum_stats false does not remove the difference (same 4 rows), so it is not the linking code.
- Also found by random sweep ($D/work/tblastn/sweep.py): `-xdrop_gap_final -5 -comp_based_stats F` (effective final = prelim 15) differs on gen_q.faa.
- Likely area (not root-caused): LOSAT's traceback keeps one `TargetTranslation` cache across the two passes (`search_gapped.rs:1412` is created before `for pass in 0..2`) and takes `stat_length` from the last window translated; NCBI makes a new `SBlastTargetTranslation` per `Blast_TracebackFromHSPList` call (core/blast_traceback.c:337) and `stat_length = range[2*context+1]` (blast_hits.c:1147-1232, kMaxTranslation 99). The e-value change tracks the subject length used (constant scale per subject).

### TN-2 (high, CONFIRMED) with -comp_based_stats 0 and a hard-masked query: e-values differ when the preliminary hit list is truncated (-max_target_seqs <= 10) with a hard-masked query
- LOSAT: search_gapped.rs:1412-1470 / stage_d_pipeline.rs:1167-1220 (same stat_length area). NCBI: blast_hits.c:43-70 (prelim hit list size), blast_traceback.c:425-433.
- Repro (generated inputs; `D=/home/kawato/.cache/losat-web-gui-target/s08p/audit`):
  `Q=$D/work/tblastn/in/gen_q_lc.faa S=$D/work/tblastn/in/gen_s.fna`
  `/home/kawato/micromamba/bin/tblastn -query $Q -subject $S -comp_based_stats 0 -lcase_masking -seg no -max_target_seqs 1 -outfmt 6` vs `$D/LOSAT tblastn` with the same arguments.
  NCBI row 1: `gl0 gs20 100.000 250 0 0 33 282 753 4 9.33e-177 479`; LOSAT: `... 9.00e-177 479` (4 more rows differ the same way: 2.24e-04 vs 2.22e-04, ...). Same for -max_target_seqs 1,5,10 (grid.py: DIFF only for lcase+seg no+soft_masking false+cbs 0+mts<=10; -max_target_seqs 100/500, -comp_based_stats 2, soft_masking true: identical). A BLOSUM45 profile variant: `-matrix BLOSUM45 -word_size 2 -comp_based_stats 0 -lcase_masking -seg no -outfmt 6 -max_target_seqs 2` -> 1 row (4.31e-157 NCBI vs 4.21e-157 LOSAT).

### TN-3 (high, CONFIRMED) outfmt 0: a subject row that is all gaps prints coordinates; NCBI prints none
- LOSAT: LOSAT/src/report/pairwise.rs:3826-3831 (`write_tblastn_alignment`: `write!(writer, "Sbjct  {}", s_pos)` ... `"  {}\n", s_end` unconditionally; the BLASTX writer `write_blastx_alignment` (pairwise.rs:4374-4462, `if qn > 0` / `if sn > 0`) already omits the numbers when a row has no residues, TBLASTN's does not). NCBI: objtools/align_format/showalign.cpp (row coordinates are left out when the row has no residues).
- Repro (generated inputs): `Q=$D/work/tblastn/in/gen_q.faa S=$D/work/tblastn/in/gen_s.fna`; `tblastn -query $Q -subject $S -xdrop_gap_final 100 -seg no` (stdout, outfmt 0), diff NCBI vs LOSAT:
  `< Sbjct       ------------------------------------------------------------  `
  `> Sbjct  545  ------------------------------------------------------------  544` (two rows, lines 2201 and 2205; an HSP with 189/312 gap letters).
  Also `-xdrop_gap_final 100 -soft_masking true -sum_stats false`. Default options: identical (no all-gap rows occur).

### TN-4 (high, CONFIRMED) huge -xdrop_gap / -xdrop_gap_final: LOSAT runs for minutes / uses tens of GB / is killed; NCBI answers in under a second
- LOSAT: stage_d_pipeline.rs:563 (prelim `gap_xdrop`), :1167-1181 (final), gapped DP memory scales with the X-drop. NCBI: blast_parameters.c:455-463.
- Fixtures (e2e_protein_query.faa x e2e_tblastn_subject.fna), NCBI always 0.1-0.6 s, 70-220 MB:
  `-comp_based_stats 0 -xdrop_gap_final 1e7`: LOSAT 3.7 s 427 MB; `1e8`: 36 s 4.1 GB; `2e8`: >60 s timeout (rc 124); `5e8`: SIGKILL (rc 137, e2e_many_subject).
  `-xdrop_gap 1e7`: 1.2 s 414 MB; `1e8`: 7.1 s 4.1 GB; `5e8`: 30 s, 20.3 GB peak (default -comp_based_stats 2). Output identical when it finishes. `-xdrop_gap 1e9` and `-xdrop_gap_final 5e8` with cbs 2: fast (INT_MIN path).
- Commands: `/usr/bin/time -f '%es %MKB' $D/LOSAT tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_tblastn_subject.fna -comp_based_stats 0 -xdrop_gap_final 1e8 -outfmt 6` (F=$D/src/LOSAT/tests/fasta/outfmt0).
- AUTHORITY.md section F/K list "任意の実数" for -xdrop_gap/-xdrop_gap_final as supported; neither a rejection nor NCBI's speed is delivered above ~1e7 bits.

### TN-5 (high, CONFIRMED; extreme e-values) BLOSUM45 profile with large -evalue: equal-score HSPs placed on different frames/coordinates
- LOSAT: BLOSUM45 14/2 word 2 threshold 16 window 60 cbs 0 path (stage_d_pipeline.rs mode-0 traceback / report ordering). NCBI: blast_hits.c ordering/heap tie-breaks.
- Repro (fixtures only): `$D/LOSAT tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_tblastn_subject.fna -matrix BLOSUM45 -word_size 2 -comp_based_stats 0 -evalue 1e4 -outfmt 6` vs NCBI: 4 of 29,900 rows differ only in the subject coordinates, e.g. `BDT62620.1 LvMJNV_160001_165000 50.000 6 3 0 39 44 4712 4695 3879 10.3` (NCBI) vs `... 4710 4693 ...` (LOSAT); same score, e-value, length. -evalue 2000 and below: identical; 5000: 6 rows; 1e5: 12; 1e10/+inf also (many rows plus ordering). BLOSUM62 default profile is identical up to 1e10 in the same files.
- Practical impact is small (E-values of 1e4+) but the value is in the accepted set (AUTHORITY.md K: "-evalue").

### TN-6 (medium, CONFIRMED) NCBI crashes (SIGSEGV) on -evalue +inf / 1e999 for some inputs; LOSAT prints a result and exit 0, no recorded exception
- LOSAT: value_parsers.rs:676 (`tblastn_real`), args.rs:58, AGENTS.md PD-LOSAT-NCBI-DEFECTS lists no -evalue crash; docs/losat_web_gui_plan.md DW-15 says NCBI's +inf/1e999 are "ported" because NCBI is deterministic.
- Repro: `/home/kawato/micromamba/bin/tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_many_subject.fna -evalue +inf -outfmt 6` -> rc 139 (SIGSEGV, deterministic over 3 runs; also with -sum_stats false, -seg no, 1e999; -comp_based_stats 0 and -max_target_seqs 5 do not crash; e2e_many_query does not crash). LOSAT: rc 0, 7221 lines, equal to NCBI at -evalue 1e300/1.7e308 (which NCBI completes, same 7221 lines). Policy says an NCBI crash with a checkable intended result is an approved exception - it needs recording, or a rejection.

### TN-7 (medium, CONFIRMED) `-dryrun` (hidden toolkit flag of tblastn) is not handled
- LOSAT: LOSAT/src/cli.rs:362 (`is_ncbi_toolkit_arg` list lacks `dryrun`; falls to the generic "unknown option" path). NCBI: corelib/ncbiargs.cpp:88,2340-2347,4124 (`s_ArgDryRun`, hidden by fHideDryRun in the help only).
- Repro: `/home/kawato/micromamba/bin/tblastn -query $F/e2e_protein_query.faa -subject $F/e2e_tblastn_subject.fna -dryrun` -> rc 0, empty stdout/stderr. `$D/LOSAT tblastn ... -dryrun` -> rc 2, `error: unknown option or argument '-dryrun'; use -help for CLI v2 syntax` (no "not supported by LOSAT's TBLASTN"). With `-dryrun -matrix FOO` NCBI is also rc 0 (nothing checked). Same for the other three programs presumably.

### TN-8 (low, CONFIRMED) messages for rejected environment/threads lack the program name required by the contract
- LOSAT: LOSAT/src/utils/threading.rs:80 ("requested N threads exceeds Rayon maximum 65535, which is not supported by LOSAT"), :162, and the NCBI_CONFIG__* registry-override rejection in main.rs ("... is not supported by LOSAT"). Contract text: "not supported by LOSAT's TBLASTN".
- Repro: `$D/LOSAT tblastn -query $Q -subject $S -num_threads 100000` -> `Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT` rc 1 (NCBI: warning "Number of threads was reduced to 32...", rc 0). `NCBI_CONFIG__BLAST__BATCH_SIZE=3 $D/LOSAT tblastn ...` -> "...is not supported by LOSAT", rc 1 (NCBI honours the entry).

### TN-9 (low, CONFIRMED) a rejection of the input or of an unsupported -outfmt pre-empts NCBI's own option error
- LOSAT: args.rs:575 (`run` rejects outfmt 5/8/custom before `check_options`), blastn/input.rs:713-760 (`read_nucleotide_subjects` rejects non-IUPAC subject residues while reading, before `check_options`).
- Repro: `-evalue 0 -outfmt 5` (NCBI: `expect value or cutoff score must be greater than zero`, rc 1; LOSAT: `output format 5 is not supported by LOSAT's TBLASTN`, rc 1); subject file with `ACGT1ACGT...` and `-matrix FOO` or `-evalue 0` (NCBI: FASTA warning + the options error, rc 1; LOSAT: `subject record 1 (s) has '1' at residue 5 ... not supported by LOSAT's TBLASTN`, rc 1). Explicit rejection (contract 3) but not NCBI's text/order; cheap to reorder for the outfmt case (NCBI parses -outfmt first, then the options, then formats).

### TN-10 (low, CONFIRMED) NCBI line references off by a few lines (session-added comments)
- protein_options.rs:390 cites blast_options.c:913-943; the quoted `if (options->gapped_calculation && !Blast_ProgramIsRpsBlast(...` is at line 910 (block 910-940).
- protein_options.rs:492 cites blast_options.c:1515-1520; "expect value or cutoff score must be greater than zero" is at 1521 (`if` at 1519).
- stage_e_report.rs (matrix_name comment) cites blast_format.cpp:2266; `options.GetMatrixName()` is at 2267.
- args.rs:373 cites blast_args.cpp:2874-2886; the `hitlist_size < 5` snippet is in the second cited range (2960-2977, line 2975) - fine but the first range is unrelated to the snippet.

### TN-11 (low, CONFIRMED) comments that do not carry a snippet or that sit above unrelated code
- args.rs:440-444: `blast_args.cpp:838-866; blast_traceback.c:1481-1501` comment has no snippet and sits above the doc comment of `search_settings`, not above the `composition_mode2` match it describes (args.rs:469-479).
- args.rs:448-453: snippet `lookup->threshold = (Int4)(kMatrixScale * opt->threshold)` (blast_aalookup.c:1292) is placed above the "compressed lookup is not supported" rejection, which is a LOSAT limit, not that NCBI code.
- args.rs:467-471 (`blast_parameters.c:774-815`, longest_intron), args.rs:485-488 (`blast_args.cpp:3152-3187` num_threads) and args.rs:955-958 (`blast_format.cpp:1443-1451`): cite file:line with at most a paraphrase; ranges verified correct (longest_intron at 793-813, kArgNumThreads constraint at 3162).

## Checked OK (byte-identical stdout, stderr and exit status unless noted)
- Harness: ~1,150 two-option combinations ($D/work/tblastn/pairs.py: 26 NCBI-rejected options in all ordered pairs; pairs2.py: 18 LOSAT-limit options x 14 NCBI-rejected options, both orders). Result: every pair is byte-identical, or the NCBI USAGE parse error (approved exception 1), or an explicit LOSAT rejection; check order (outfmt parse, -seg, matrix/gap pair, threshold, word size, evalue) matches NCBI in every combination. NCBI crash cases (-max_target_seqs 1073741799 with cbs 2, etc.) are rejected by LOSAT as documented.
- Message text/exit code identical: -matrix FOO (rc 1, matrix list), -gapopen/-gapextend pairs (1/1, 5 alone, -1), -threshold 0, -word_size 8/9, -evalue 0/-1, -seg "12 2.2"/"a b c"/4 tokens, -outfmt 99 (rc 255) / x / -6 / 99999999999, -ungapped with cbs 2/1/3 (Composition-adjusted searched ...), IDENTITY with word size 6, -outfmt 13 to stdout, empty query/subject, missing files, -out to a directory, BATCH_SIZE 0 (rc 3).
- Parse errors where NCBI prints USAGE (LOSAT text, exit 2; approved exception): -threshold -1, -word_size 0/1/-3, hex/negative/overflow ints for typed options, -evalue inf/nan (no sign), -window_size -1, -max_intron_length -1, -soft_masking/-sum_stats garbage, -task Tblastn/blastp, -num_threads 0.
- Boolean spellings for -soft_masking and -sum_stats: true/false/t/f/T/F/yes/no/Yes/NO/y/n/on/off/ON/Off/1/0 identical outputs; maybe/garbage/""/" true"/2/-1/tru rejected by both.
- Doubles for -xdrop_gap (56 spellings: 5., .5, +5, 5e0, 5E+0, inf/+inf/-inf, +nan(1), 1e400, 4.9e-324 ...) identical; hex floats and `1e` forms are rejected by LOSAT with "not supported by LOSAT's TBLASTN" (NCBI accepts: contract 3, documented D6).
- Supported profiles: BLOSUM62 11/1 w3 t13 win40 with -comp_based_stats 0 and 2, BLOSUM45 14/2 w2 t16 win60 with cbs 0 (also with explicit -gapopen 14 -gapextend 2 -threshold 16 -window_size 60, matrix spellings blosum45/Blosum45/bLoSuM45), outfmt 0/6/7, on e2e_* fixtures and on generated sets (6 queries, 42 subjects, both strands, ambiguity codes, lowercase): identical with defaults and for: -evalue 1e-300..1e300 / +inf / +nan (BLOSUM62), -max_target_seqs 1..1000 (all 4 data sets), -seg no/yes/12 2.2 2.5/20 2.0 2.5/10 1.0 1.0/0 0 0/5 3 3, -sum_stats true/false, -soft_masking true/false, -lcase_masking (queries and subjects with lowercase) with/without -seg no, -xdrop_gap 0/1/1e10/1e300/+inf/-5/100/1000 and -xdrop_gap_final 0/1/1e10/+inf/-5/100/1000 with the default -comp_based_stats 2 (differences only with cbs 0, TN-1/TN-4), BATCH_SIZE 1..100000 on 6 queries, -num_threads 1..1000, -out FILE for outfmt 0/6/7, stdin query/subject (not both: both-from-pipe is an explicit rejection).
- Random sweeps ($D/work/tblastn/sweep*.py, 6 runs, 1 to 6 random options of 11 per run, 700+ runs): all identical except the TN-1..TN-5 classes.
- LOSAT rejections with the required phrase (checked): -word_size 2/4/5/6/7 (non-BLOSUM45), -task tblastn-fast (with word 3, 4, 5), -window_size 0/1/41, -threshold 11/13.5/13.9/14/1e999/+inf, -matrix PAM30/PAM70/PAM250/BLOSUM50/BLOSUM80/BLOSUM90/IDENTITY, BLOSUM62 with other gap pairs, BLOSUM45 with cbs 2 or other word size, -comp_based_stats 1/3, -ungapped, -max_intron_length != 0, -outfmt 1-5/8-12/15-16/18/20 and custom specs/delim, hex doubles, ADAPTIVE_CBS (cbs 2), OLD_FSC, BL2SEQ_LEGACY, BATCH_SIZE=x, PRE_FETCH_SEQS_LIMIT=x, header-only/non-IUPAC records. No panic or silence seen on any of them. Every option NCBI lists in `tblastn -help` (65) is either implemented or answered with "the NCBI BLAST+ option -X is not supported by LOSAT's TBLASTN"; toolkit options (-logfile, -conffile, -version-full*, -help-full, -xmlhelp) likewise, except -dryrun (TN-7).
- NCBI file:line references spot-checked: 66 fenced snippets in session-added comments of tblastn/args.rs, stage_d_pipeline.rs, stage_d_results.rs, stage_e_report.rs, blastinput/app.rs, cli.rs, value_parsers.rs, stats/protein_options.rs verified mechanically against the NCBI tree ($D/work/tblastn/verify_refs-style script): 63 exact; plus 22 unfenced references checked by hand (blast_args.cpp 3425-3427, 1029-1056, 997-1003, 220-229, 280-286, 258-273, 1938-1942, 2549-2557, 2745-2748, 2801-2851, 2894-2927; blast_parameters.c 455-463, 774-815; blast_traceback.c 286-300, 405-433; seqsrc_multiseq.cpp 226-240; blast_seqsrc.h 205; blast_setup.c 969-973; local_db_adapter.cpp 133-135; blast_seqalign.cpp 672-674, 1484-1490; prelim_stage.cpp 145-147): all correct. Exceptions are TN-10/TN-11.
- Adapter (web/adapter/src/run.rs:92-96): `validate` for Tblastn is `parse` + `tblastn::check_options` = `TblastnArgs::resolve` + `search_settings`, the same two functions the CLI runs (args.rs `search_cli` and `search`), so option-only parity holds by construction; it rejects -out/-outfmt before parsing (run writes all formats), drops the "Examining 5 or more matches" warning, and does not run the environment checks (CLI does). `run_local` order (subject records, title warnings, check_options, empty-query warning) equals the CLI's. The session added no TBLASTN case to the adapter tests (only blastp/tblastn dispatch lines).

## Counts
high 5 (TN-1, TN-2, TN-3, TN-4, TN-5), medium 2 (TN-6, TN-7), low 4 (TN-8..TN-11). All CONFIRMED by running both binaries; TN-1/TN-2 root cause is a suspected area only.
