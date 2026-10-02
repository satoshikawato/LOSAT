# Range C result notes: the description table of `tblastx -outfmt 0` (final code `0d533ba76`)

Result TSV: `result_C.tsv` (41 rows, one per row of `C.tsv`). Scratch, inputs and every oracle output:
`rscratch_C/` (`cmp.sh NAME Q S [opts]` runs NCBI 2.17.0 `tblastx` and the S08 LOSAT binary with `-outfmt 0` on the same paths and
compares stdout byte for byte; files `runs/NAME.{ncbi,losat}.{out,err}`; sets in `rscratch_C/sets/`, the older sets of the inventory in `scratch_C/sets/`).

## 1. Result

No GAP and no UNSURE row for TBLASTX. The description table is written by the one shared writer
`write_blastn_description_table` (`report/pairwise.rs:1783`, new parameter `show_sum_n`), called with `show_sum_n = true`
from `write_tblastx_pairwise_report` (`pairwise.rs:2800-2810`). The N column, the first-HSP N, the best-bit HSP for Score/E,
the total-score width rule (last-row assignment), the `num_descriptions` slice, the title cleaning and the 68-byte cut all match NCBI
on every input I ran (about 100 full outfmt 0 comparisons, listed in section 3; all byte equal, stderr equal).

Row 38 (the "switch question") is `n/a` for TBLASTX but leaves a finding for other programs (section 4).
Incidental finding outside the table: `-culling_limit` (section 5).

## 2. How each group of rows was decided

- Rows 1, 3, 5, 6, 10, 28: the table gate, `eShowSumN`, the header and the blank lines. Header bytes, `cat -A` checked in `multi1` (N width 2, trailing space
  after one-digit N) and `link` (N width 3: `102`, `159`, `84 `), and in the score-width-9 fixture `tblastx.lc.0`.
- Rows 7-9 (`-max_target_seqs`): NCBI `blast_args.cpp:2924-2928` (outfmt <= 4 branch: both counts = N, hitlist = N) and the default 500/250; LOSAT
  `run_impl.rs:1154-1157` and `final_hit_order` (`report.rs:765`, shared `HitList` port) with the warning from `process_options` (`run_impl.rs:759`).
  Own oracle with 540 subjects (`many540b`): default 500 rows and 250 alignments; `-max_target_seqs` 1, 5, 510, 530, 2147483647 equal (stdout and stderr).
  Ties at the hit list border: 12 identical subjects (`ties_s.fa`, `-max_target_seqs` 1, 3, 5, 11, 12) and `dup_s.fa` (1 to 5) equal (higher oid first on ties).
- Rows 12-20 (`x_InitDeflineTable`, `GetSeqAlignSetCalcParams`): LOSAT builds, per described subject, (best-bit HSP with strict `>` in Seq-align order,
  total of all bit scores, N of the first HSP). I checked the NCBI source again (`showdefline.cpp:1051-1160`, `align_format_util.cpp:4214-4315`) for the last-row quirk:
  `if (m_MaxTotalScoreLen < total.size()) m_MaxScoreLen = total.size()` (score width assigned, total width not updated for the last row); it applies only when the subjects
  fit in the descriptions (last group non-empty). Oracle at the border: `big3_XY` with `-max_target_seqs 2` (subjects = descriptions: quirk applies) and
  `-max_target_seqs 1` (does not), `big3_XYY2` with 1, 2, 3: equal.
- Rows 24-26, 32-33: number formats. Oracle sets with E from `0.0` to `5939`, bits from ` 9.3` to `1.392e+05`: equal, also with `CTOOLKIT_COMPATIBLE=1` in the environment of that one command
  (header `(bits)`, same on both sides; the NCBI 2.17.0 binary does not apply the compile-time E-value shift).
- Rows 29-31: titles. Oracle `titles4` (16 MAG:/TPA_asm:/UNVERIFIED:/TLS: prefixes, a 118-byte title with space runs at the cut, `id7 x.,;~ .`, `id8;, ~` quirk, `multispecies:`),
  `titles5` (21: `lcl|`, `gnl|BL_ORD_ID|`, `ref|`, `gi|`, `Subject_5` prefixes, bracket modifiers `[organism=..]`, `[topology=..]`, `[gcode=4]`, `[lineage=x]`, `[`, `]`, `[]x`),
  `title_ascii` (66 to 120 bytes), `crlf_s/crlf_both` (CRLF files), `trail` (trailing blanks and a trailing tab): all equal. LOSAT prints no id column; NCBI hides every `lcl|Subject_N` id
  (`showdefline.cpp:893-909`), which is every `-subject` id, and the same text comes out.
- Row 36 (`s_EvalueCompareHSPLists`): subject order equal on `ties_s` (c11 first), `dup_s`, `big3_YY2`, the fixtures `tblastx.multi.0` and `tblastx.many.0`.
- Row 37 (non-ASCII titles): rejected by `check_report_titles` (`run_impl.rs:889-930`). Oracle: `scratch_C/sets/title2_s.fa` and `title_s.fa` (records with a 2-byte character): LOSAT exits 1 with
  "Error: subject record N has a defline that is empty, starts with white space or has a control character or a non-ASCII byte, ...". The check looks at every record of both files, hit or not.
- Rows 2, 11, 27, 39-41: no byte effect (inventory reasons kept).

## 3. Oracle runs (all `IDENTICAL`, stdout; stderr equal where it is not empty)

| name | input | options | what it covers |
|---|---|---|---|
| multi1, ev | multi1 (12 subjects, first HSP != best-bit HSP), ev (14 subjects) | -, `-evalue 100000` | N widths, best-bit row, E formats |
| big_A_only, big_A_small, big_A_B, big_B_only, small_B, big2_XY, big2_Y, big2_YY2, big2_XYY2, big3_XY, big3_YY2, big3_XYY2 | big_q vs big sets (up to 3539 HSPs, totals 128618 to 139175 bits) | - | last-row width rule, total > 99999, duplicate subjects |
| big3_XYY2_mts2/3/1, big3_XY_mts1/2 | same | `-max_target_seqs` | last-row rule at subjects == descriptions, hit list |
| sc1_* (24) | 24 window pairs of LC738874/LC738875 | - | first-HSP bits < best-bit (3 of 24), front HSP single while others linked |
| titles4, titles5, title_ascii, dup, crlf_s, crlf_both, trail | title/defline sets | - | title text, cut, ids, duplicate deflines |
| ties_mts1..12, dup_mts1..5 | 12 identical subjects; 5 subjects with duplicate deflines | `-max_target_seqs` | hit-list ties |
| many540b (+ mts1, 5, 510, 530, 2147483647) | q1500 vs 540 mutated copies | `-max_target_seqs` | more than 500 subjects: 500 rows, 250 alignments; stderr |
| link | colinear chunks with spacers (3 subjects) | - | N of 3 digits (102, 159, 84) |
| lowbits | 60 short subjects | `-evalue 1e9 -threshold 11` | bits < 10 (` 9.3`), large E |
| ctk_multi1, ctk_big2Y, ctk_lowbits | as above | env `CTOOLKIT_COMPATIBLE=1` | `(bits)` header |
| fr1..fr8, sw1..14 | see section 6 | | fresh genome pairs, option sweep |

## 4. Other programs (row 38)

The final code has three table writers: `write_blastn_description_table` (BLASTN, TBLASTN since S08 at `pairwise.rs:3378`, TBLASTX), and the older
`write_subject_summary_table_with_sum_n` (first HSP, plain score width; BLASTP `pairwise.rs:1576,2986`, BLASTX `pairwise.rs:3823`). NCBI's code has no program test, so
BLASTP and BLASTX can differ from NCBI for a subject whose first HSP has fewer bits than another HSP of the same subject, or whose last row total exceeds 99999 bits.
I ran 12 earlier BLASTX window sets (`rscratch_C/bx`, inputs from `scratch_C/scan_tn2`): all byte equal, so I have no concrete BLASTX case; the TBLASTN cases of the inventory (p417, p716) were the ones S08 fixed.
This is outside the TBLASTX scope and is not counted as a GAP.

## 5. Incidental findings outside the table

1. `-culling_limit N` (N >= 1) differs from NCBI for TBLASTX, also in outfmt 6, so it is an engine/option matter (E2e, `docs/losat_web_gui_sessions/session_s08p_e2e_protein_options.md` lists `-culling_limit` among the unverified TBLASTX options),
   not the table. Repro: `tblastx -query rscratch_C/sets/q1500.fa -subject rscratch_C/sets/many540b_s.fa -culling_limit 1 -outfmt 6`: NCBI 21 lines, LOSAT 0 lines (outfmt 0: LOSAT prints `***** No hits found *****`).
   With `-culling_limit 2`: NCBI 42 lines, LOSAT 11563; `-culling_limit 5`: NCBI 109, LOSAT 13183. On `multi1` with `-culling_limit 2` the outfmt 0 reports differ in the HSP set (LOSAT keeps/loses different HSPs). The option is accepted silently.
2. Inventory section 7 items now: `-max_target_seqs` applied, warning below 5 printed (rows 7-9); empty subject defline is rejected for outfmt 0/7 (row 37), while outfmt 6 still prints `unknown` where NCBI prints `Subject_N` (deferred: `Subject_N` ids).

## 6. Fresh runs and option sweep

Fresh genome pairs (full outfmt 0, byte equal): fr1 LC738868 vs LC738869 (1.1 MB), fr2 PemoMJNVA vs PemoMJNVB (23.8 MB), fr3 LC738884 vs LC741431, fr4 two queries vs three subject records,
fr5 `outfmt0/width_*`, fr6 `outfmt0/longdef_*` (long deflines), fr7 `outfmt0/many_*`, fr8 AP027131 vs AP027133 (4.3 MB).
Option sweep (`sweep2.sh`, two inputs each, all byte equal): `-seg no`, `-threshold 11`, `-threshold 15`, `-window_size 25`, `-evalue 10`, `-evalue 1e-3`, `-query_gencode 4`, `-num_threads 3`,
`-max_target_seqs 7`, `-seg yes`. `-word_size 2` and `-word_size 4` stop in LOSAT with "unsupported TBLASTX word_size: only 3 is implemented" (explicit rejection, exit 2; NCBI runs them): not a table matter.
`-seg no` on the 12-subject `multi1` pair was not finished (NCBI itself needs minutes because of the low-complexity runs); it was stopped, not a difference.

## 7. Inventory corrections

None of the 41 rows is wrong about NCBI behaviour; the notes of the inventory (best-bit HSP for Score/E, first HSP for N, last-row width quirk, higher oid first on ties, 68-byte cut)
were all reproduced. Two small points: (a) row 4 has `-sum_stats` as `rejected` in the inventory, but the S08 decisions list it as deferred; I used `deferred`;
(b) the inventory's statement "BLASTX is the only LOSAT caller that passes show_sum_n=true" is still true outside TBLASTX; TBLASTX now passes it as well (`pairwise.rs:2807`).

## 8. What the port does now, per NCBI reference (checked against the final code)

1. Table for every query with hits, between `Length=` and the alignments: header line 1, header line 2 with the `N`, empty line, rows, empty line, empty line (`blast_format.cpp:1539,606-613,1545`; `showdefline.cpp:757-849`) -> `pairwise.rs:1832-1856,1902,2809`.
2. Widths by `x_InitDeflineTable` (score 6, E 5, N 1; last-row total assigned to the score width) -> `pairwise.rs:1815-1831` (`showdefline.cpp:1055-1062,1100-1124,1133-1158`).
3. Per subject: best-bit HSP (strict `>`, Seq-align order) for Score and E, first HSP N (`num > 1 ? num : 1`), total of all bit scores -> `pairwise.rs:1791-1810` (`align_format_util.cpp:4214-4313`, `blast_seqalign.cpp:1180-1198`).
4. HSPs per subject in `Blast_HSPListSortByEvalue` order and subjects in `s_EvalueCompareHSPLists` order with the hit list of `-max_target_seqs` (default 500) -> `report.rs:765` (`blast_seqalign.cpp:1577`, `blast_hits.c:1396-1457,3078-3110,3243-3300`).
5. Number and title formats -> `outfmt6.rs:256,313`, `defline.rs:94` (`align_format_util.cpp:940-1004`, `showdefline.cpp:893-930`).
6. 500 descriptions and 250 alignments by default, N and N with `-max_target_seqs N`, warning below 5 -> `run_impl.rs:759-772,1154-1157` (`blast_args.cpp:2910-2977`).
7. Titles that NCBI reads differently (empty, leading white space, control or non-ASCII byte, HTML reference, punctuation past the end) are rejected -> `run_impl.rs:889-930`.
