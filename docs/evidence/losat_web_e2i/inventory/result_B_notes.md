# E2i result notes, range B (C++ options layer and option checks)

Final code `90c5f0181`; NCBI 2.17.0 as oracle; scratch `/home/kawato/.cache/losat-web-gui-target/sd-res-B` (script `cmp.sh` runs both programs and compares stdout, stderr, exit status; runs are named `cNNN`).

## Counts (rows 1-48 of B.tsv plus extras X1-X5)
ported 15, faithful 15 (+5 extras faithful), n/a 11, rejected 7, exception 0, GAP 0, UNSURE 0.

## GAP rows
None found.

## UNSURE rows
None.

## What was checked
- Every option value of the two tasks in the final `task_defaults` / `configure_task` against `blast_options_handle.cpp:344-380`, `blast_nucl_options.cpp` and `disc_nucl_options.cpp`:
  - dc-megablast: word 11, reward 2, penalty -3, gaps 5/2, e-value 10, DUST on, mask_at_hash on, DP extension and traceback, gap x-drop 30 / final 100, ungapped x-drop 20 bits (x_dropoff_init not zeroed), min diag separation 6 (megablast value, SetMBHitSavingOptionsDefaults is not overridden), window 40, scan_range 0, template coding/18, lookup table MB. All equal to NCBI.
  - blastn-short: word 7, reward 1, penalty -3, gaps 5/2, e-value 1000, no DUST by default (`-dust` can switch it on), DP, x-drop 30/100, min diag separation 50, window 0, no template, small/normal lookup table by BlastChooseNaLookupTable, chunk 1,000,000. All equal to NCBI.
  - Settings that are not LOSAT options (threshold, stride, hitlist 500, mask level 101, perc identity 0) are shared by all tasks.
- NCBI `-export_search_strategy` for the four tasks confirms MBTemplateType 0, MBTemplateLength 18, WordSize 11, WindowSize 40 for dc, and MatchReward 1, MismatchPenalty -3, EvalueThreshold 1000, WordSize 7, DustFiltering FALSE, MaskAtHash TRUE for blastn-short.
- Option checks: NCBI order is BlastExtensionOptionsValidate, BlastScoringOptionsValidate (penalty, gap extension 0), LookupTableOptionsValidate (word size limit 100, template word size 11/12, lookup table type), BlastInitialWordOptionsValidate (off_diagonal_range, unreachable: option rejected), BlastHitSavingOptionsValidate (e-value), s_BlastExtensionScoringOptionsValidate (zero gap costs). `check_scoring_options` has the same order and the same texts. Stdout, stderr and exit status equal to NCBI in 30 runs on both tasks (c8-c13, c16-c27, c601-c607), including combinations in which two checks fail at once and the `-dust` token errors (which NCBI raises during option extraction, before Validate) combined with check errors.
- Value matrices, all equal to NCBI (stdout, stderr, exit status): default outfmt 0/6/7 for both tasks; blastn-short -word_size 4-13, 16, 28; dc -word_size 11, 12; the 18 template combinations on dc; -evalue, -reward, -penalty, -gapopen, -gapextend singles and pairs on both tasks (including NCBI errors for unsupported gap costs, c501-c539); -dust yes/no/three numbers on both tasks with low-complexity inputs; -max_hsps, -perc_identity, -max_target_seqs, -subject_besthit, -lcase_masking; queries of 4 to 40 residues, IUPAC, lower case, reverse complement, all-N (c701-c712); 2.5 Mb and 22-record queries with the default and with CHUNK_SIZE=300000/120000 (c401, c402, d_*).
- Approved exceptions seen: clap text and exit 2 for `-task DC-megablast`, one-sided `-template_type`/`-template_length`, `-template_length 20`, positive `-penalty` (c14, c15, c28, c29, c119-c121); no `-num_threads` warning with `-subject` (t1, t2).
- Explicit rejections confirmed with both tasks (exit 2, "the NCBI BLAST+ option -X is not supported by LOSAT's BLASTN"): -window_size, -off_diagonal_range, -xdrop_ungap, -xdrop_gap, -xdrop_gap_final, -no_greedy, -ungapped, -soft_masking, -use_index, -index_name. `-task rmblastn`: "the task rmblastn is not supported by LOSAT's BLASTN ...".

## Remarks
- Row 26 (chunk size of dc, 5,000,000): ported (`task_uses_megablast_chunks`), but I found no input whose output depends on 1,000,000 vs 5,000,000: NCBI's own output for a 2.5 Mb dc query is identical with `CHUNK_SIZE=1000000` and the default, and LOSAT equals NCBI under `CHUNK_SIZE=300000/120000`. The row is checked by reading the code (`run.rs:5523`, `query_split.rs:23`) and the unit test only.
- Row 6: the final code keys the discontiguous lookup table on `config.mb_template_length > 0` (as NCBI) instead of `task == dc-megablast`; classified faithful because dc behaved correctly before; megablast with a template (c123, c124) and blastn / blastn-short with a template (c122, c23, c24) follow NCBI.
- Row 48: `zero_gap_extension_formula` is limited to megablast and blastn (run.rs:6106); for dc and blastn-short zero gap costs never reach the epilog because they stop in the greedy check, so no output differs.
- Not in my range, noted: NCBI blastn-short with `-reward 0 -penalty 0` gives output byte-equal to the default run (435 lines), dc with 0/0 gives 31 lines (default 22); LOSAT stops for 0/0 with "a reward of 0 ... is not supported by LOSAT's BLASTN" (`check_losat_limits`, explicit rejection from E2g).
- `min_ungapped_score` in `TaskConfig` (`MIN_UNGAPPED_SCORE_*`, marked DEPRECATED) is read only in debug `eprintln!` lines of `run.rs`; no effect on output for any task.

## Inventory corrections
None of the 48 rows was wrong. Notes that only need updating: row 6 (`discontig_template` is now keyed on the template), row 26 (ported, see remark).

## Extra rows (not in B.tsv)
X1 BLAST_ValidateOptions call order; X2 BlastScoringOptionsValidate gap-extension check; X3 BlastHitSavingOptionsValidate e-value check; X4 CBlastAppArgs::SetOptions / CFilteringArgs::ExtractAlgorithmOptions order (-dust errors before Validate; -dust enabling DUST for blastn-short after ClearFilterOptions); X5 LookupTableOptionsValidate word-size limit 100. All faithful for both tasks.
