# Range E result notes (gapped extension, traceback, hit saving; dc-megablast and blastn-short)

Final code 90c5f0181, binaries `LOSAT-90c5f0181` and NCBI 2.17.0. Scratch: `/home/kawato/.cache/losat-web-gui-target/sd-res-E` (inputs, `run.sh`, `sens.py`, outputs).

## GAP rows
None. No input was found on which LOSAT and NCBI differ for dc-megablast or blastn-short in this range.

## UNSURE rows
None.

## What was run (all outfmt 6 unless noted, stdout + stderr + exit status compared with `cmp`)
dc-megablast (DP, min_diag_separation 6):
- 28 self-comparisons: LC738874 14 kb windows `w00..w19` (step 14 kb) and synthetic tandem repeats `t0..t5` (unit 12..1000 bp, 5..20% substitutions and 1..8 bp indels), `tandem`, `tandem2`: all byte-identical.
- 20 kb windows `win0..win2`, reverse-complement subjects (`rc_t0`, `rc_t4`, `wrc`), nested window (`wsub` vs `win0`), overlapping windows (`ov1` vs `ov2`), 700x700 subjects (`msub` vs `msub`, 44376 lines), LC738874 vs itself (670 HSPs, `big1`) and vs LC738875 (`big2`): all identical.
- Options: -max_hsps/-evalue, -word_size 12 -template_length 21 -template_type optimal, -reward 1 -penalty -3 -gapopen 2 -gapextend 2: identical. outfmt 0 and 7 on `w15` identical.
- Sensitivity of the inputs to 6 versus 50: a script (`sens.py`) lists, in NCBI's own output, HSP pairs where the lower-scoring HSP is box-contained in the higher-scoring one and both end diagonals differ by >= 6 but one by < 50 (kept by c=6, deleted by c=50 at the final purge). Such pairs exist in `win0`, `w14`, `w15`, `win3`, `t0`, `t4` and `big1` (2 to 8 pairs each, e.g. win0: (2375-2524, 2262-2412, 151) holds (2414-2524, 2264-2375, 113) with diagonal differences 37 and 37); LOSAT reproduces them, so the final purge uses 6. The init-hit and the traceback pre-check cannot be isolated from the CLI; they read the same captured variable (see reader list) and the HSP sets of all those runs are equal.
blastn-short (DP, min_diag 50, e-value 1000):
- 60 queries of 15..80 bases from LC738874 vs 30 kb/40 kb/20 kb windows (`s1,s2,s6`), 30 queries 22..200 bases (`s3`, `s9` with -word_size 4: 42984 lines), 20 kb vs 30 kb and 20 kb vs 20 kb windows (`s4,s5`): identical.
- 5 queries x 700 subjects (`m1`, `m2` -max_target_seqs 20, `m3` -max_hsps 2 -num_threads 4, `m4` -word_size 5 -evalue 100000: 84510 lines; every query reaches the 500-subject hit list limit): identical (m3 stderr differs by NCBI's approved num_threads warning only).
- E-value: default (1000) gives 9716 lines on `s6`; -evalue 10 gives 3041 lines in both programs (`s7`); explicit -evalue 1000 equals the default (`ev3`); -evalue 20000, 1e-3, 1e-30, 0 (exit 1, same message) identical.
- Short subjects/queries: 61 queries of 4..20 bases at the end of a 5003-base subject with -word_size 4,5,7,8,10 and default (`tiny*`): identical (row 19 rollback).
- DUST off by default and on with `-dust yes` on a low-complexity pair (`d1,d2`); dc-megablast keeps DUST (`d3,d4`): identical. outfmt 0 and 7 on `w15` identical (83952 and 9035 lines).
- Rejections: -xdrop_gap, -xdrop_gap_final, -ungapped, -no_greedy with both tasks: exit 2 and "the NCBI BLAST+ option -X is not supported by LOSAT's BLASTN". 0/0 gap costs with both tasks: both programs exit 1 with the same message.

## Readers (grep over LOSAT/src at 90c5f0181)
- `min_diag_separation`: written only by `coordination.rs` `task_defaults`/`configure_task` (587, 618) into `TaskConfig`; read in `run.rs` at 7191 (one local) and passed to `containing_hsp` at 10328 (init hit vs gapped HSPs, BLAST_GetGappedScore), 11266 (speculative traceback prefilter), 11360 (traceback pre-check) and 12293 (final purge). `interval_tree.rs` takes it as a parameter. Not read: `BlastnArgs.min_diag_separation` (args.rs:155, `#[arg(skip)]`, set 0 in web_api.rs:256) and `blast_engine/mod.rs:257` (see below).
- `use_dp`: `coordination.rs:558` (`!defaults.greedy`), `run.rs` 6276 (speculative traceback only if DP), 7171, 10401, 10576, 11520, 11643, 11680, 11765, 11969; `gapped.rs` parameter; `blast_engine/mod.rs:67` `_use_dp` unused. All keyed on the flag, never on the task name.
- Gapped X-dropoffs: `configure_task` 563/569 -> `run.rs:7188` `gap_x_dropoffs` (min lambda scaling) -> 10425/10512/10576 (`x_drop_gapped`) and 11310/11499 (`x_drop_final`). Gap trigger 27 is task independent (`run.rs:7500`).
- E-value threshold: `determine_evalue` only (`run.rs:7194`, `scoring.rs:147`); consumers `run.rs` 7491 (cutoff_score_max), 10995 (preliminary reap), 12409 (final reap). No remaining reader of `args.evalue` with a fixed 10.

## Inventory remarks
- All LOSAT line numbers of the inventory are of a92fa902f; the result TSV has the final numbers.
- Rows 12 (`s_BlastSetUpAuxStructures`, status n/a "LOSAT has no discontiguous table") and 14 (CreateTask) are stale: SD ported the discontiguous word finder and the task. Row 15's "args.rs:49 evalue default 10.0" is stale (now `Option<f64>`).
- The E.tsv has 23 data rows (the file shows 24 lines with the header; there is no 24th row).
- Row 5 remark about `blast_engine/mod.rs:257` `filter_hsps` hard-coding 0: still true. `filter_hsps` is `pub` in the public module `blast_engine` but has no caller in the crate, tests or benches (grep), so it is not on the CLI, web ABI or `run` path. Not a gap; if the function is ever called it would ignore the task's separation.

## Extra rows
X1 `BlastHitSavingOptionsValidate` (faithful), X2 `BlastHspNumMax` / -max_hsps (faithful), X3 `Blast_HSPListReapByRawScore` branch for matrix scoring (rejected via reward 0 / rmblastn).
