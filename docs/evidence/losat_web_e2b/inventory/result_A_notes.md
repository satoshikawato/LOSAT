# Result notes, range A (tblastx application layer and CBlastFormat orchestration)

Files: `result_A.tsv` (132 rows, one per row of `A.tsv`), this file, scratch in `rscratch_A/` (inputs `*.fa`; for every oracle comparison `NAME.out/.err/.rc` = NCBI 2.17.0 and `NAME.lout/.lerr/.lrc` = LOSAT, made by `rscratch_A/cmp.sh NAME args...`; the row appender is `add.py`, the row scripts `rows_*.py`).

## 0. What was compared with what

- NCBI: `/home/kawato/micromamba/bin/tblastx` 2.17.0 (no `.ncbirc`; BATCH_SIZE only where a case says so). Source: `ncbi-blast/c++` (598d8ae6).
- LOSAT "final code" = commit `0d533ba76`. All `losat_location` line numbers in `result_A.tsv` are those of that commit (I extracted it with `git archive` to `rscratch_A/committed/`).
- WARNING about the binary and the working tree: while I worked, the LOSAT working tree got uncommitted edits (`run_impl.rs`, `blastn/input.rs`, `tests/tblastx_regression_fixtures.py`) and `s08/native/release/LOSAT` was rebuilt from them at 18:53:31. The edits (a) reject `-culling_limit > 0`, (b) reject records without residues and residues that are not IUPAC nucleotide letters (`check_residues_of`, `check_records_have_residues_of` with program "TBLASTX"), (c) read `U` as `T` (`with_u_as_t`). All comparisons made BEFORE 18:53 (`cmp.sh`, outputs `NAME.lout/.lerr/.lrc`) used the binary of the commit. To be safe I then built the commit in scratch (`cargo +1.92.0 build --release --locked --offline`, target dir `rscratch_A/target`, source `rscratch_A/committed/LOSAT`; nothing in the repositories was touched) and re-ran every GAP, every observation of section 3 and a regression sweep of about 100 cases with `cmpc.sh` (outputs `NAME.clout/.clerr/.clrc`, prefixes `c_`, `k_`, `cm_`, `cf_`). All GAPs and observations reproduce on the pristine build; the sweep is byte-identical except `mix_empty.fa` (deferred row 69). The classification in `result_A.tsv` is about the commit; where the uncommitted edits already change a result I say so.

## 1. Counts (result_A.tsv)

ported 51, reused 25, faithful 4, rejected 3, exception 3, deferred 6, n/a 34, GAP 6 (rows 6, 26, 41, 52, 116, 126). UNSURE (kept as n/a, as in the inventory): rows 30 and 88.

## 2. GAP rows in detail

All commands are `tblastx -query Q -subject S ...` for NCBI and `LOSAT tblastx ...` for LOSAT, run in `rscratch_A/`.

### Row 6 (query input config): empty or blank query defline in outfmt 6 (outfmt 0/7 are rejected)
- Input: `emptytitle.fa` = a line `>` followed by the q1 sequence (also `>   `, and `> abc def` with a leading space), `s1.fa` subject, `-outfmt 6`.
- NCBI first column of every row: `Query_1` (for `> abc def`: `abc`), e.g. `Query_1\tswin0\t47.826\t23\t12\t0\t1682\t1750\t436\t504\t0.002\t25.4`.
- LOSAT: `unknown\tswin0\t47.826\t...` (same length of output, 3235 B each, bytes differ). Files `et_emptytitle_6.out/.lout`, `et_spacetitle_6.*`, `lt_lead_space_title.*`.
- Subject counterpart (the deferred `Subject_N` ids): `-subject` with `>` as defline: NCBI `Subject_1`, LOSAT `unknown` (`st_s_emptytitle.*`); `> abc def`: NCBI `abc`, LOSAT `unknown` (`st_s_spacetitle.*`).
- outfmt 0 and 7 reject these deflines with `query record 1 has a defline that is empty, starts with white space ...`.

### Row 26 (CObjReaderParseException; input outside the declared scope but neither reproduced nor rejected)
- `bin.fa`: `>a`, then the bytes `00 01 02` and `ACGT`. NCBI: stderr `BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` rc 1; outfmt 0 stdout is the 379-byte prolog. LOSAT: complete report (760 B), rc 0, stderr empty (`binq.*`, `bins.*`).
- Residues that CFastaReader ignores with a diagnostic (`digits.fa` = `ACGTACGT12ACGTTT-GACC...*`): NCBI stderr `CFastaReader: Hyphens are invalid and will be ignored around line 2` + `FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 9-10, 67`, and the sequence shrinks (outfmt 7 row `a swin0 46.154 13 7 0 60 22 ...`); LOSAT keeps the characters (`... 63 25 ...`), no stderr (`digits.*`, `in_badres.*`, `uq_x_q_*`, `uq_star_q_*`).
- Text before the first defline (leading blank line, `;` comment, UTF-8 BOM, no defline): NCBI reads the file; LOSAT stops with `Error: failed to read query FASTA ... Expected > at record start.` (rc 1) (`in_lead_blank.*`, `in_nodefline.*`, `in_semicolon.*`, `in_bom.*`). This is an explicit stop, not a gap, but the message is bio's.
- The uncommitted edit (b) above already rejects the non-IUPAC residues and records without residues.

### Row 41 (GetSubjectFile): `-subject` file name that is not valid UTF-8
- `cp s1.fa "$(printf 'caf\xe9.fa')"`; `-subject` that file, `-outfmt 7` (or 0): NCBI `# Database: User specified sequence set (Input: caf<E9>.fa)` (raw byte), LOSAT writes U+FFFD (EF BF BD): stdout 3421 vs 3423 B (`sp_nonutf.*`). UTF-8 names (relative, absolute, spaces, wrapped at 68 columns, `café.fa`) are identical.

### Row 52 (CStdCmdLineArgs::ExtractAlgorithmOptions): `-out` the same file as `-query`
- `cp q1.fa qcopy.fa; tblastx -query qcopy.fa -subject s1.fa -out qcopy.fa -outfmt 7`: NCBI stderr `Warning: [tblastx] Query is Empty!`, rc 0, the file is truncated to 0 bytes (it opens the query stream, truncates `-out`, reads later). LOSAT reads the whole query first and writes a 3419-byte report into the file, rc 0, no warning. With `-out` the same file as `-subject` both programs agree (the subject is read first).

### Row 116 (PrintEpilog threshold line): an integer `-threshold` of 1000000 or more
- NCBI prints the double with the C++ default precision (6 significant digits, `%g`): `-threshold 1000000` -> `Neighboring words threshold: 1e+06`; `-threshold 1234567` -> `1.23457e+06`. LOSAT (i32 printed with `{}`) -> `1000000`, `1234567` (stdout 2268 vs 2270 B, 2274 vs 2270 B; `thr_1000000.*`, `thr_1234567.*`). Up to 999999 identical. Non-integer thresholds are rejected by LOSAT's value parser (the deferred `ncbi_double` item); the epilog must print `%g` when that is ported.

### Row 126 (AcknowledgeBlastQuery / title text): first white space of a defline is not a space
- LOSAT rebuilds the title as `record.id() + " " + record.desc()` (`algorithm/tblastx/report.rs:485`). `bio` splits the header at the first white-space character (tab, VT, FF, NBSP and every Unicode white space included) and drops it, so `check_report_titles` (`run_impl.rs:889`, which rejects control characters and non-ASCII) checks the rebuilt string and never sees it.
- Query defline `abc<TAB>def`, `-outfmt 7`: NCBI `# Query: abc`, LOSAT `# Query: abc def`; `-outfmt 0`: `Query= abc` vs `Query= abc def` (stdout 13537 vs 13541 B; `tq_tab_0/7`, `tq_vt_*`, `tq_ff_*`, `tq_tabsp_*`, `tq_nbsp_*`). Subject defline with the same shape: description row and `> abc` heading are `abc` in NCBI, `abc def` in LOSAT (`ts_tab_0`, `ts_vt_0`, `ts_ff_0`, 13545 vs 13549 B).
- A tab AFTER the first white space (`abc def<TAB>ghi`) and a control character as first character (`abc\x01def`) are rejected as intended (`tq_soh_first_*`).

## 3. Differences on the in-scope path that belong to no row of range A (found while testing; for the other ranges or the parent)

- `-culling_limit N>0` (committed binary): outfmt 6 on `rq.fa` vs `rs.fa` (7 queries x 5 subjects): `-culling_limit 1` NCBI 5507 B vs LOSAT 0 B; `2`: 6523 vs 6082; `5`: 7403 vs 7210 (`r_cull*_6`, `r_culling_limit2_0/7`). The uncommitted edit (a) now rejects it.
- RNA letter `U` in a query or subject (committed binary): `u_q.fa` (q1 with every T as U): NCBI output equals the one with T; LOSAT finds other HSPs (`uq_u_q_7`: `# 58 hits found` vs `# 54`; minus-strand frames differ). The uncommitted edit (c) addresses it.
- `-seg "3 -5 -5"` (negative cut-offs), outfmt 7 on q1 vs s1: NCBI 150 B, LOSAT 25589 B (`seg_3_-5_-5`). Engine range; not analysed further.
- `-word_size 2` runs in NCBI (27965 B for q1 vs s1) and is an explicit rejection in LOSAT (`unsupported TBLASTX word_size: only 3 is implemented`); `-window_size 0` and BATCH_SIZE text are the listed rejections.
- Records without data (deferred item), committed build: a header-only query (`hdronly.fa`) is an invalid query of an unsearched batch in LOSAT (warning `Could not calculate ungapped Karlin-Altschul parameters`, report, rc 0) while NCBI exits 3 after the prolog with `BLAST engine error: Warning: Sequence contains no data `; a header-only subject: NCBI `Warning: [tblastx] Subject_1 hdr only: Subject sequence contains no data` then `BLAST engine error: The average subject length is too short` rc 3 (380-byte prolog for outfmt 0), LOSAT writes the same 380 bytes but `Error: invalid subject length` rc 1 and no warning; `mix_empty.fa` (valid, `>emptyseq`, all-N): NCBI warns `Query_2 emptyseq: Sequence contains no data `, LOSAT is silent (stdout identical); empty subject record inside a subject set (`s_withempty.fa`): NCBI warns `Subject_2 emptysubj: Subject sequence contains no data`, LOSAT silent (stdout identical).

## 4. Rows that are not plain `ported`/`reused`

- deferred: 24 (CInputException texts: -outfmt text, missing -subject, Validate), 25 (query/subject not accessible; `-out` is ported with the same text), 35 and 53 and 54 (-outfmt text and range: NCBI rc 1 / 255, LOSAT clap rc 2), 69 (records without data).
- rejected: 40 (`-query -`), 55 (custom field list), 56 (formats 17/19/21 and all others).
- exception: 2 (clap), 33 (outfmt 6/7 write failure), 129 (thread warnings). The `-db_gencode` exception shows in rows 97 and 108 (cmp `dbg2_*`, `dbg4_*` differ from NCBI by design).
- UNSURE kept as n/a: row 30 (CIOException other than eFlush leaves the status 0; no trigger found) and row 88 (`results.HasErrors()` error branch: the only per-query error-severity messages are the filtering failure `Failure at filtering` (blast_filter.c:1292) and internal failures; odd `-seg` windows (1 0 0, 2 1 1, 100000 2.2 2.5, 2147483647 ...) produce the same warnings in both programs, no error).

## 5. Oracle runs worth keeping

- Interleaving (`2>&1`) of stdout and stderr: identical for outfmt 0, 6, 7 on `b_seq2.fa`, `b_seq1.fa`, `two_allN.fa`, `mixed.fa -max_target_seqs 2`, `allN.fa -max_target_seqs 3`, BATCH_SIZE=1 and 300 (`mg_*`).
- Batches: query lengths 10001/10002/10003/5001+5001 plus a 3000-nt and a 700-nt query (`bd_*`), `b_seq1/2/3.fa`, BATCH_SIZE 0, 1, 500, 700, -1 and text (text: NCBI rc 255 `CStringException`, LOSAT rejects, rc 1).
- Hit lists: 300 identical subjects with default / 1 / 4 / 5 / 280 for outfmt 0, 6, 7 (`s300*`); q1 vs s2 with 1 and 2.
- 7 queries x 5 subjects (`rq.fa`/`rs.fa`) outfmt 0 and 7 with default, `-evalue 1e-5`, `-evalue 1000`, `-threshold 11`, `-threshold 15 -window_size 30`, `-seg no`, `-seg yes -window_size 60`, `-query_gencode 4`, `-max_target_seqs 3`, `-max_target_seqs 2 -evalue 100`, `-num_threads 4`: byte-identical (only the thread warning differs: approved exception). `-culling_limit 2` differs (section 3).
- Degenerate queries (`deg_*`: poly-A/T/C/G, `TAA`x100, `TAA`x100 + sequence, `TGA`x40 + sequence, `AC`x150, R x300, `ACGN`x80, stop-rich, `ATGNNN`x50): outfmt 0 and 7 identical, including warnings.
- Write failure: `-out /dev/full` outfmt 0 (q1, all-N, mixed): `BLAST failed to write output`, rc 6, stdout empty; empty query: `Query is Empty!` rc 0; outfmt 6/7: NCBI SIGABRT rc 134, LOSAT rc 1 (approved); stdout to `/dev/full` outfmt 0: rc 6 both.

## 6. What the port still has to do (open items for the parent)

1. Row 126: take the title from the raw defline (or reject a tab/VT/FF/non-ASCII white space as the first white space) in `report.rs:485` / `check_report_titles`.
2. Row 6: the generated `Query_<n>` / `Subject_<n>` ids for empty or space-leading deflines in outfmt 6 (or reject them in outfmt 6 as in 0/7).
3. Row 116: `%g` formatting of the threshold (needed for integer values >= 1e6 now, and for non-integers when `ncbi_double` is ported).
4. Row 52: open `-out` before reading the query file (NCBI order), or reject `-out` equal to `-query`.
5. Row 41: reject or reproduce non-UTF-8 `-subject` names.
6. Row 26: invalid residues / binary bytes (the uncommitted edit (b) rejects them).
