# Result notes, range F (formatting)

Binaries: LOSAT `/home/kawato/.cache/losat-web-gui-target/sd/bin/LOSAT-90c5f0181` against NCBI `/home/kawato/micromamba/bin/blastn` 2.17.0. Harness and outputs: `/home/kawato/.cache/losat-web-gui-target/sd-res-F` (`cmp.sh`, `matrix*.log`, `out/`, inputs in `in/`).

## GAP rows
None.

## UNSURE rows
None.

## Comparisons run (stdout, stderr and exit status compared; all equal unless noted)
- matrix.log: 84 runs. Tasks dc-megablast and blastn-short x outfmt 0, 6, 7 x 14 inputs: multi-query (hit, all-N, no-hit, reverse strand, lowercase, mixed case), no hit, all-N, two all-N queries, queries of 30/5/12/20 bases, dc_query/dc_subject, short_query/short_subject, short_nohit_query, -max_target_seqs 3 and 1 (12 subjects, with NCBI's "Examining 5 or more matches" warning), -lcase_masking on a partly and an all-lowercase query and on the multi-query file, and the same without -lcase_masking.
- matrix2.log: 78 runs. dc-megablast with -template_type coding/optimal/coding_and_optimal x -template_length 16/18/21 x -word_size 11/12 x outfmt 0/6/7 on dc_query/dc_subject (54), plus templates 16 and 21 on the multi-query, lowercase, all-N and -max_target_seqs inputs (24).
- matrix3.log: 36 option runs for both tasks (9 option sets) x outfmt 0 and 7: -reward 1 -penalty -2, -reward 2 -penalty -3, -gapopen 2 -gapextend 2, -reward 1 -penalty -2 -gapopen 0 -gapextend 2, -gapopen 0 -gapextend 0 and -gapopen 5 -gapextend 0 (both fail identically in NCBI and LOSAT: same stderr, rc 1), -reward 5 -penalty -4 -gapopen 8 -gapextend 6, -evalue 1000 -word_size 7 -dust no (for dc-megablast also an identical template/word-size error, rc 1), -evalue 1e-50 -perc_identity 90 -max_hsps 1. All 36 equal. The remaining 12 runs in matrix3.log are -window_size 40, -window_size 0 and -use_index true with dc-megablast, blastn-short, megablast and blastn: all rejected by LOSAT (exit 2, "error: the NCBI BLAST+ option -window_size|-use_index is not supported by LOSAT's BLASTN"), NCBI runs them (rc 0). Known decision.
- matrix4.log: option runs outside the formatter path (custom -outfmt field lists, -num_alignments, -strand, -soft_masking) are rejected by LOSAT with explicit "not supported by LOSAT's BLASTN" messages for both tasks and are not task-specific; the runs with -lcase_masking -dust no, an IUPAC query, edge_short_invalid and edge_batch_allN equal NCBI for both tasks.

Observed footer/epilog facts (identical in NCBI and LOSAT): dc-megablast ends with `Matrix: blastn matrix 2 -3`, `Gap Penalties: Existence: 5, Extension: 2`, `Window for multiple hits: 40`; blastn-short ends with `Matrix: blastn matrix 1 -3`, `Gap Penalties: Existence: 5, Extension: 2` and no Window line, with `1.37 0.711 1.31` in both Karlin blocks; the prolog of both is the Altschul reference with `BLASTN 2.17.0+`.

## Inventory notes
- The inventory file has 20 rows (not 21; the line count includes the header and a trailing empty line).
- Rows 7, 10, 16 and 18 were divergent/unported before the port and are `ported` now; the other rows kept their status. Row 5 and 9 are `n/a` with their reasons verified by runs (the `-use_index` rejection; the identical NCBI and LOSAT option errors for a zero gap extension).
- Inventory LOSAT line numbers are of `a92fa902f`; result_F.tsv uses the final code (`90c5f0181`).

## Extra rows
- X1: `CBlastFormat::PrintEpilog` prints the Gap Penalties line only when `options.GetGappedMode()`; LOSAT rejects `-ungapped`, so gapped mode is always on (`rejected`).
- X2: `x_PrintOneQueryFooter` skips the Karlin-Altschul blocks of an invalid (all-N) query; LOSAT matches NCBI on all-N inputs (`faithful`).
