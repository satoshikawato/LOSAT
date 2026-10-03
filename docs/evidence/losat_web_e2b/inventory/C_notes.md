# Range C notes: the description table of `tblastx -outfmt 0`

Source numbering is the CR-stripped numbering of `c++/` at commit 598d8ae6. Evidence files are under
`/home/kawato/.cache/losat-web-gui-target/s08/inventory/scratch_C/` (`sets/`, `scan1/`, `scan_tn2/`, helper scripts `lib.py`, `xl.py`, `mk_*.py`, `scan*.py`).
The TSV is `C.tsv` (41 rows). Oracle = `/home/kawato/micromamba/bin/tblastx` 2.17.0 (and `tblastn` for the switch question). I also ran the pre-built, unmodified LOSAT binary
`/home/kawato/.cache/losat-web-gui-target/s08-base/release/LOSAT` (read-only use, nothing was built) to compare TBLASTN/TBLASTX outfmt 6 and TBLASTN outfmt 0.

## 1. Call path and order of writes

1. `PrintOneResultSet` (`blast_format.cpp:1411`ff.): query preamble, then for a query with hits `aln_set = results.GetSeqAlign()`;
   `m_IsUngappedSearch` (tblastx: `tblastx_options.cpp:73` `SetGappedMode(false)`) -> `PrepareBlastUngappedSeqalign` (1519-1521) splits each subject's Seq-align
   (one Std-seg per HSP, `blast_seqalign.cpp:1397-1450`) into one Seq-align per HSP, in HSP order.
2. Table gate `(!m_IsBl2Seq || m_IsDbScan) && !(m_DisableKAStats || kIsGlobal)` (1539): `-subject` is bl2seq (121) but in DbScan mode (122; `blast_app_util.cpp:206-210`), so the table is printed.
3. `x_DisplayDeflines` (566-614): `CShowBlastDefline showdef(*aln_set, scope, 68, m_NumSummary + 0)`; `x_ConfigCShowBlastDefline` sets `eShowSumN` only (498-523);
   `showdef.DisplayBlastDefline(out)`; then `"\n"` (613). `PrintOneResultSet` writes another `"\n"` (1544) before the alignments.
4. `DisplayBlastDefline` (`showdefline.cpp:1004-1011`) = `x_InitDeflineTable()` then `x_DisplayDefline(out)`. **`x_InitDeflineTable` runs, never `x_InitDefline`** (that one is only reached through `Init()` with no templates, which blast_format never calls).
5. Nothing is written before `x_DisplayDefline`, so the header is written after all widths are known (all rows computed first). The table has no flushes of its own.

Stream order for one query: `...Length=N` line, header line 1, header line 2, empty line, rows, empty line (613), empty line (1544), `> subject` heading.

## 2. What is printed per subject (exact rule)

`x_InitDeflineTable` collects consecutive Seq-aligns with the same subject id into a set `hit` (loop 1076-1131; at most `m_NumToShow` subjects), then `x_GetScoreInfoForTable(hit)` (1518-1560)
-> `GetSeqAlignSetCalcParamsFromASN` finds no `seq_hspnum` (BLAST never writes it), so `GetSeqAlignSetCalcParams(hit, m_QueryLength, false)` (`align_format_util.cpp:4247-4313`):

| shown | value | source |
|---|---|---|
| Score | highest `bit_score` over ALL HSPs of the subject; strict `bits > highest_bits` from 0, scanning in Seq-align order (first of equal bit scores wins) | `afu:4285-4300` |
| E | the E-value of that same HSP (variable is named `lowest_evalue`; it is NOT the subject's lowest E-value) | `afu:4299-4302,4310` |
| N | `sum_n` of the FRONT HSP of the subject (`GetSeqAlignCalcParams(aln.Get().front())`, 4259, 4230); `-1` (no `sum_n` score) -> 1; never touched by the loop | `afu:4215-4245` |
| total (not shown) | sum of the bit scores of all HSPs; used only for the width rule | `afu:4285,4311` |

Which Seq-align scores are read (`GetAlnScores` 811-861 -> `s_GetBlastScore` 230-263): `score`, `bit_score`, `e_value` **or** `sum_e` (the same variable `evalue`; BLAST writes exactly one: `blast_seqalign.cpp:1186-1189`
`score_type = (hsp->num <= 1) ? "e_value" : "sum_e"`), `sum_n` (only written when `hsp->num > 1`, 1180), `num_ident`, `comp_adjustment_method`, `use_this_gi`. If the Seq-align has no `e_value/sum_e`, the Std-seg's scores are read (`seg.GetStd().front()`).
E is `0.0` when `hsp->evalue < 1e-180` (`SMALLEST_EVALUE`, `blast_seqalign.cpp:62,1186`); `GetScoreString` uses the same threshold.

Order of HSPs inside a subject (the Seq-align order, decides "first/front"): `Blast_HSPListSortByEvalue` (`blast_seqalign.cpp:1577`; `blast_hits.c:1396-1457`): E-value ascending (two E-values below 1e-180 count as equal),
then raw score descending, then subject offset ascending, subject end descending, query offset ascending, query end descending. The NCBI comment says why: "lower scoring HSPs might have lower e-values, if they are linked with sum statistics" (1396-1400).
Subject order: `s_EvalueCompareHSPLists` (`blast_hits.c:3078-3110`): best E-value (same epsilon), first HSP raw score descending, then **higher oid first** (identical subjects print in reverse input order: `sets/title3.out`, `sets/big3_YY2.out`).
Because each tblastx query context has its own ungapped Karlin block (`blast_stat.c:2764-2790`), equal raw scores can have different bit scores, so the highest bit score need not be the highest raw score; use `Hit.bit_score`.

## 3. Widths, header and row bytes

Start widths (`1055-1062`): score 6 (`(Bits)`), E 5 (`Value`), N 1, total 5 (`Total`). Per row (mid-loop 1100-1124): score = max(score, len(bit string)), total = max(total, len(total string)), E, N (digits of sum_n).
**Last row** (1133-1158, only when the last `hit` group is non-empty, i.e. subjects <= m_NumToShow, which is always true with in-scope options because the engine keeps at most `hitlist_size = max_target_seqs` subjects):
score = max(score, len(bit string)); then `if (m_MaxTotalScoreLen < total string length) m_MaxScoreLen = total string length` (assignment to the SCORE width, 1140-1142; the total width is not updated for the last row).
A total string is longer than 5 only above 99999 bits (`%5.3le`, 9 chars). So the score column is 9 wide iff the last row's total exceeds 99999 bits AND no earlier row's total did. Oracle:
`sets/big_A_only.out` and `scratch_C/full1.out` (one subject, total 139175) width 9; `sets/big2_XY.out`, `sets/big3_XY.out` (31596 then 110666) width 9; `sets/big3_XYY2.out`, `sets/big3_YY2.out` (earlier row also >99999) width 6;
`sets/big3_XYY2_mts2.out` (`-max_target_seqs 2`, second row is last) width 9; `sets/big_A_small.out` (first total 128618, last small) width 6.

Header (`x_DisplayDefline` 757-849), first-pass, `eShowSumN`:
```
line 1: AddSpace(68+2) "Score" AddSpace(Wscore-5) AddSpace(2) AddSpace(2) "E" "\n"                       (no trailing spaces)
line 2: "Sequences producing significant alignments:" AddSpace(68-43) AddSpace(1) "(Bits)" AddSpace(Wscore-6) AddSpace(2)
        "Value" AddSpace(Weval-5) AddSpace(2) "N" "\n"                                                 (no trailing spaces, no pad after N)
"\n"
```
Row (886-985): title (68 columns, see below), `"  "`, bit string padded to Wscore, `"  "`, E string padded to Weval, `"  "`, N padded to Wn, `"\n"`.
Trailing spaces exist only after N when its column is wider than that row's N.

Bytes (`cat -A`, title area shortened with ..):
```
full1.out  (LC738874 vs LC738875, one subject, Wscore 9):
                                                                      Score        E$
Sequences producing significant alignments:                          (Bits)     Value  N$
$
LC738875.1 MAG: Melicertus latisulcatus pemonivirus Okinawa2016 D...  776        0.0    4$
$
$
> LC738875.1 MAG: ...
```
Why `776        0.0    4`: Wscore = 9 because the subject's total (139175 bits, `1.392e+05`) is 9 chars; `776` + 6 pad + 2 margin = 8 spaces; E `0.0` padded to Weval 5 (+2) + 2 margin = 4 spaces; N width 1.
```
sets/multi1.out (12 subjects, 11 rows, Wscore 6, Weval 6, Wn 2):
                                                                      Score     E$
Sequences producing significant alignments:                          (Bits)  Value   N$
sub12 window ..    612     0.0     12$
sub6 window ..    148     8e-100  17$
sub11 window ..    138     3e-82   11$
sub7 window ..    116     9e-77   8 $          <- trailing space (Wn 2)
sub2 window ..    138     7e-65   4 $
sub10 window ..    145     1e-34   1 $
sub5 window ..    55.6    3e-12   3 $
sub3 window ..    38.6    0.017   3 $          <- first HSP 26.7 bits N=3; shown HSP is the 38.6-bit singleton
sub8 window ..    29.9    7.2     1 $
sets/ev.out (-evalue 100000): "0.021", "0.97", "1.3", "12", "157" as E strings; "27.2", "20.3", "16.6" bit strings
```
Title (`913-930`, `498`): `sdl->defline = CDeflineGenerator().GenerateDefline(handle, fLeavePrefixSuffix)` (cleaned title: trailing `.,;~ ` trimmed, space runs compressed, `MAG:`-style prefixes kept);
the id `lcl|Subject_N` is hidden (893-909), so `line_length` is 0 and the title is the whole defline. `size() > 68` bytes -> first 65 bytes + `...`; else padded to 68. Byte based (`sets/title.out`: 66/67/68 padded, 69/70/120 cut;
`sets/title2.out` M6 cut inside a 2-byte character). Title cleaning oracle `sets/title3.out`: `a  b   c`->`a b c`, `T. `->`T`, `x;,~`->`x`, `MAG: x` stays (the `> ` heading shows `x`), empty defline -> empty title and id `Subject_N` in outfmt 6.

Number formats (`GetScoreString`, `afu:940-1004`): E: `<1e-180` `0.0`; `<1e-99` `%2.0le`; `<0.0009` `%3.0le`; `<0.1` `%4.3lf`; `<1` `%3.2lf`; `<10` `%2.1lf`; else `%2.0lf`.
Bits: `>99999` `%5.3le`; `>99.9` `%3.0ld` of `(long)bits`; else `%4.1lf` (so `' 9.9'`). Total: same but `%2.1lf` at or below 99.9. LOSAT's `format_evalue_ncbi` / `format_bitscore_ncbi` implement these.

## 4. LOSAT today

- `write_subject_summary_table_with_sum_n` (`report/pairwise.rs:799`; BLASTP 1576, TBLASTN 2727, BLASTX 3162, 2347): header with N is faithful (incl. `max_evalue-5+2`), but it implements the dead `x_InitDefline` logic: first HSP of the subject (score, E, N) and a plain score width. It also builds the label as id + title without `ncbi_nucleotide_title` (TBLASTN: `a  b   c.` printed where NCBI prints `a b c`).
- `write_blastn_description_table` (1765): the full `x_InitDeflineTable` rule (highest bit score, that HSP's E, total of all HSPs, last-row width assignment with `last_row_counted`, title cleaning, 68-byte cut), no N column.
- TBLASTX has no outfmt 0/7 (`error: invalid value '0' for '-outfmt': unsupported TBLASTX outfmt: only 6`), `-sum_stats`, `-sorthits`, `-num_descriptions` are unknown options.
- TBLASTX's engine knows the linked count (`ExtendedHit.num`, `algorithm/tblastx/chaining.rs:131`) but the final `Hit` (`run_impl.rs:2825-2877`) and `PairwiseHit::from(Hit)` (`pairwise.rs:120`) drop it.
- TBLASTX outfmt 6 equals NCBI byte for byte on 3 of my 4 sets (ev, multi1, big); the 4th differs only in subject ids of empty-defline records (NCBI `Subject_N`, LOSAT `unknown`, see 7).

### What TBLASTX needs
One writer = `write_blastn_description_table` plus a `show_sum_n` parameter:
(a) header: append `AddSpace(Weval-5) + AddSpace(2) + "N"` to line 2 (as `write_subject_summary_table_with_sum_n` does);
(b) per row: N = `sum_n` of `hits.first()` (Seq-align order) as `num > 1 ? num : 1` (not Some(0)), printed as `"  " + N` padded to `Wn`, with `Wn = max(1, digits(N))` over the described rows including the last;
(c) keep everything else as is: score and E from the highest-bit HSP (strict >), the last-row total rule, title via `ncbi_nucleotide_title(.., true)`, 68-byte cut, two blank lines after;
(d) carry `num` from the TBLASTX engine into `PairwiseHit.sum_n` and keep each subject's HSPs in `Blast_HSPListSortByEvalue` order (`common.rs::ncbi_order_evalue_hsp_order`);
(e) apply `-max_target_seqs` to the number of subjects in the engine (default 500) before the writer; the writer's `described`/`last_row_counted` logic then holds.

## 5. Is the G.3 rule NCBI's rule for every program? Which BLASTP/TBLASTN/BLASTX outputs change?

Yes. `CShowBlastDefline::DisplayBlastDefline`, `x_InitDeflineTable`, `x_GetScoreInfoForTable`, `GetSeqAlignSetCalcParams` and `GetScoreString` contain no test on the program. The only program-dependent inputs are
`eShowSumN` (from `use_sum_statistics && m_IsUngappedSearch`, so only ungapped programs: tblastx, `blastx -ungapped`), the line length 68, the number of descriptions and the Seq-align data.
If BLASTP/TBLASTN/BLASTX moved to the G.3 table:
- **TBLASTN and BLASTX** (both set `SetSumStatisticsMode()` even when gapped: `tblastn_options.cpp:70`, `blastx_options.cpp:89`): a subject whose first HSP (lowest E) has fewer bits than another HSP of the same subject now shows that HSP's bit score and its (larger) E-value.
  Oracle, TBLASTN, 931 six-frame ORFs (>=60 aa) of LC738874 vs LC738875 (`scan_tn2/q_all.fa`, `scan_tn2/s_all.fa`, `scan_tn2/o_all_tn.txt`): 3 of 738 (query, subject) rows differ from the first-HSP rule:
  `p417` HSPs (18.1, Expect(2)=4.8), (17.7, Expect(2)=4.8), (20.0, Expect=7.3): NCBI row `20.0    7.3  `, LOSAT `18.1    4.8  `; `p716` (19.6 Expect(2)=0.76 ...; 20.4 Expect=1.8 ...): NCBI `20.4    1.8  `, LOSAT `19.6    0.76 `.
  Repro files: `sets/tn_p417_q.fa`, `sets/tn_p716_q.fa` against `scan_tn2/s_all.fa`; outputs `sets/tn_p417.out` (NCBI) and `sets/tn_p417.losat.out`; `diff` of the two reports from `Query=` on shows only the row. Smaller random runs (418 rows, 12 windows) showed none: rare, but real. BLASTX behaves the same by code (same linking), not run.
- **BLASTP** (no sum statistics; one Karlin block; E-value order = raw score order): the first HSP already has the highest bit score, so only the last-row width rule can change a report: a subject with more than 99999 total bits (many HSPs on one protein) widens the score column to 9 for the whole table.
- All three: the score column is also widened when the last row's total exceeds 99999 bits (TBLASTN genome subjects with many HSPs).
- Titles: NCBI cleans the title through `CDeflineGenerator` for protein subjects too; the shared writer does not (TBLASTN shows the raw title; for TBLASTN `ncbi_nucleotide_title(.., true)` is the right call, BLASTP/BLASTX need the protein title rule: not inventoried here).
- N column: no change for gapped runs (column off); `blastx -ungapped` would also get N from the first HSP with score/E from the best-bit HSP.
Recommended order: switch TBLASTN (nucleotide titles available) after TBLASTX proves the writer; BLASTP/BLASTX need the protein title cleaning first. Any switch must be checked against the frozen BLASTP/TBLASTN/BLASTX fixtures: only rows with first-HSP != best-bit HSP or totals >99999 can change.

## 6. Oracle input where the first HSP has fewer bits than the best one (tblastx), and the row NCBI prints

All use query `LC738874` windows and subject `LC738875` windows (`scan1/q_*.fa`, `scan1/s_*.fa`, outputs `scan1/o_*.txt`; HSP lists parsed from the alignment section):

| input | first HSP (Seq-align order) | best-bit HSP | row printed |
|---|---|---|---|
| `scan1/*_123796_183796_199027_219027` (query 123796-183796, subject 199027-219027, 22 HSPs) | 26.7 bits (52), `Expect(3)` = 3e-05 | 30.9 bits (61), `Expect` = 0.14 (no sum_n) | `s199027 ...  30.9    0.14   3` |
| `scan1/*_129875_189875_235662_255662` (470 HSPs) | 137 bits, Expect(19) = 0.0 | 140 bits, `Expect(7)` = 2e-117 | `  140     2e-117  19` |
| `scan1/*_182408_242408_233511_253511` (742 HSPs) | 138 bits, Expect(22) = 0.0 | 167 bits, Expect(3) = 8e-58 | `  167     8e-58  22` |
| `sets/multi1_q.fa` vs `sets/multi1_s.fa`, subject `sub3` | 26.7 bits, Expect(3) = 9e-04 | 38.6 bits, Expect = 0.017 | `sub3 ...  38.6    0.017   3 ` |
| `sets/big_q.fa` vs `sets/big2_Y_s.fa` (3539 HSPs) | 232 bits, Expect(80) = 0.0 | 274 bits, Expect(49) = 8e-133 | `  274        8e-133  80` |

Rows show: score and E from the best-bit HSP, N from the front HSP (even though the best-bit HSP has a different or no `sum_n`). When the front HSP has no `Expect(n)` (single) and other HSPs do, N prints 1
(`scan1/o_200199_260199_229579_249579.txt`: first HSP `Expect = 1e-04`, row `  40.9    1e-04  1`).

## 7. Incidental findings outside the table (for the owners of other ranges)

- LOSAT TBLASTX never applies `-max_target_seqs` (field `max_target_seqs` in `algorithm/tblastx/args.rs:57` is unused). `sets/big_q.fa` vs `sets/big3_XYY2_s.fa` with `-max_target_seqs 2`: NCBI outfmt 6 has 2 subjects (2769 lines), LOSAT 3 subjects. The default of 500 subjects is not applied either.
- NCBI prints `Warning: [tblastx] Examining 5 or more matches is recommended` for `-max_target_seqs` 1..4 (`blast_args.cpp:2975-2977`); LOSAT prints nothing.
- Empty defline (`>` or `> `): NCBI id `Subject_N` (outfmt 6 and heading), LOSAT TBLASTX outfmt 6 prints `unknown` (`sets/title3.ncbi6` vs `sets/title3.losat6`).
- A title that ends in a multibyte character loses its last byte in NCBI (`>E3 ü` -> `E3 \xC3`, `sets/title2.out`); LOSAT BLASTN rejects non-ASCII deflines (`algorithm/blastn/input.rs:105`), TBLASTX accepts them in outfmt 6; `ncbi_nucleotide_title` has `expect("ASCII deflines")`. A decision is needed for TBLASTX outfmt 0 (UNSURE row).
- `AddSpace(size_t)` wraps for negative widths; not reachable for `-subject` (ids hidden, widths never below a row's size).

## 8. What the port must do (with NCBI references)

1. Print the table for every query with hits, between the `Length=` line and the alignments: header line 1, header line 2 with the trailing `N`, empty line, one row per subject, empty line, empty line (`blast_format.cpp:1539,606-613,1544`; `showdefline.cpp:757-849,852-985`).
2. Use the x_InitDeflineTable rule for all widths, including the last-row assignment of the total-score width (`showdefline.cpp:1055-1062,1100-1124,1133-1158`).
3. For each subject choose the HSP with the highest `bit_score` (strict `>`, Seq-align order) for Score and E; take N from the first HSP in Seq-align order with `num > 1 ? num : 1`; total = sum of all HSP bit scores (`align_format_util.cpp:4247-4313`, `blast_seqalign.cpp:1180-1198`).
4. Keep HSPs per subject in `Blast_HSPListSortByEvalue` order and subjects in `s_EvalueCompareHSPLists` order (higher oid first on ties); this is already what the outfmt 6 path produces (`blast_seqalign.cpp:1577`, `blast_hits.c:1396-1457,3078-3110`).
5. Format E, bits and total with `GetScoreString` rules (`align_format_util.cpp:940-1004`), title with `ncbi_nucleotide_title(.., true)`, 68-byte cut with `...`, hidden `lcl|Subject_` ids (`showdefline.cpp:498,893-930`).
6. Carry `num` from the TBLASTX engine to the report and apply `-max_target_seqs` (default 500) to the subject count before writing, including the stderr warning below 5 (`blast_args.cpp:2910-2977`, `blast_hits.c:3243-3300`).
7. Do not copy `write_subject_summary_table_with_sum_n`'s first-HSP/plain-width rule for TBLASTX (it models the unused `x_InitDefline`).

## 9. Oracle commands run (all with outputs under `scratch_C/`)

- `tblastx -query LC738874.fasta -subject LC738875.fasta -outfmt 0 -out full1.out` (5.7 s): header/row bytes above.
- Window scans `scan1.py` (24 query x subject windows; 3 of 24 show first-bit < best-bit), `mk_multi1.py` (12 subjects), `mk_big*.py` (width rule, duplicates for tie order, `-max_target_seqs 2`), `mk_title*.py` (title cut, cleaning, multibyte), `ev` set with `-evalue 100000` (E formats).
- `tblastn` with 931 ORF queries (`scan_tn3.py`) and per-query repros; LOSAT `tblastn -outfmt 0` on the same repros; LOSAT vs NCBI `tblastx -outfmt 6` on four sets (`cmp`).
- Not run: sets above 500 subjects (range A), `-num_descriptions`, `BL2SEQ_LEGACY`, HTML, `-sorthits`.
