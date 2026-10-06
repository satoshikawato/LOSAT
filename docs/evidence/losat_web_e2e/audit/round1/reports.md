# Audit: BLASTP / TBLASTN reports (outfmt 0, 6, 7)

Harness: work/reports/cmp.sh (runs NCBI and LOSAT, compares stdout, stderr, combined 2>&1 and exit code byte for byte; outputs under work/reports/runs/<tag>/{n,l}.{out,err,both,rc}). Inputs under work/reports/in/.
About 1,000 NCBI/LOSAT command pairs were run for the items below.

## Findings

### RP-1 (medium, CONFIRMED) Order of the two query warnings is reversed when a query has O residues and is invalid

- LOSAT: LOSAT/src/algorithm/blastp/blast_engine.rs:6937-6948 and LOSAT/src/algorithm/tblastn/args.rs:961-973 (messages built as [O-replaced, then invalid-Karlin]); comment at tblastn/args.rs:959-960 states the O message comes before the Karlin one; LOSAT/src/report/query_warnings.rs:107-125 (`query_warning`) joins them in that order.
- NCBI: src/algo/blast/api/blast_setup_cxx.cpp:920-932 (O message is appended to `warnings`), src/algo/blast/format/blast_format.cpp:1450-1452 (prints `results.GetWarningStrings()`). The NCBI binary puts the Karlin-Altschul message FIRST, then the O message, on the same Warning line.
- Applies to every query that is invalid (all X / O / `*` after the O->X replacement) and contains at least one O: for blastp and tblastn, outfmt 0, 6, 7 (stderr and 2>&1).
- Repro:
  `printf '>oo\nOO\n' > q.faa; blastp -query q.faa -subject e2e_protein_subject.faa -outfmt 6` (NCBI vs `$D/LOSAT blastp ...`, compare stderr)
  (inputs: work/reports/in/o6.faa, o1.faa, o2.faa, o7.faa, o8.faa, o9.faa, mixE.faa; runs ov16, ov1, ov4, ov19, ov22, ov25, iv25)
- Diff (`diff runs/ov16/n.both runs/ov16/l.both`):
```
< Warning: [blastp] Query_1 oo: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options One or more O characters replaced by X for alignment score calculations at positions 0, 1 
---
> Warning: [blastp] Query_1 oo: One or more O characters replaced by X for alignment score calculations at positions 0, 1 Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options 
```
  Same for tblastn (iv26, ov1t...). Fix: put INVALID_QUERY_MESSAGE first (blastp/blast_engine.rs:6941-6947, tblastn/args.rs:966-972) and fix the comment at tblastn/args.rs:959-960. A valid query with O residues (e.g. 3 O among normal residues, or 25 O, mixF, o3, o4, o5) matches NCBI exactly.

### RP-5 (medium, CONFIRMED) Query-defline "Title ends with at least 50 valid amino acid characters" warning is written before the report prolog; NCBI writes it after the prolog (outfmt 0, 2>&1)

- LOSAT: LOSAT/src/algorithm/blastp/blast_engine.rs:4480 and LOSAT/src/algorithm/tblastn/args.rs:721 call `fasta_input::write_protein_title_warnings(&queries, outputs.diagnostics)` while the query file is read, before the prolog (`write_blastp_pairwise_prolog` / `write_tblastn_pairwise_prolog`) is written. Definition: LOSAT/src/algorithm/blastn/input.rs:223-236.
- NCBI: the warning is posted by the FASTA reader when the first query batch is read (src/objtools/readers/fasta.cpp:1650-1673 via blastp_app.cpp/tblastn_app.cpp reading the batch after CBlastFormat::PrintProlog, see the `BATCH_SIZE` note in AUTHORITY.md section H "prolog is written before the first query batch is read"). In the NCBI output the warning therefore appears after the "Database: ... N sequences; M total letters" block and before the first `Query=` line. (Subject-title warnings are read before the prolog and match; only query-title warnings are misplaced.)
- Only visible when stdout and stderr go to the same file/pipe; separate stdout/stderr and outfmt 6/7 are identical (outfmt 6 has no prolog, outfmt 7 prints its header after the read).
- Repro (qh11.faa: `>` + 80 `x` + sequence; qw1.faa: second record title `g1 AAAA...(60)`):
  `blastp -query in/qh11.faa -subject e2e_protein_subject.faa -outfmt 0 2>&1` and the same with `tblastn ... -subject e2e_tblastn_subject.fna`; runs qwqh110, qwtqh110, qwqw10, qwqw20, qwqh120.
- Diff (`diff runs/qwqw20/n.both runs/qwqw20/l.both`, two warning lines, query file in/qw2.faa):
```
0a1,2
> FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?
> FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?
23,24d24
< FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?
< FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?
```
  (NCBI line 23 is after the `Database:` block, i.e. between the database header and `Query=`.)

### RP-2 (low, CONFIRMED) Wrong NCBI line range in a session-added comment (tabular.cpp:179-188)

- LOSAT: LOSAT/src/algorithm/blastp/blast_engine.rs:2215 `NCBI reference: c++/src/objtools/align_format/tabular.cpp:179-188,1132-1149`; the comment says "eAccession uses GetLabel(...,0), eAccVersion uses fLabel_Version".
- NCBI: the GetLabel(..., 0) call is at tabular.cpp:178 and the fLabel_Version one at tabular.cpp:184 (case blocks are 175-186); 179-188 holds `break;`, `}`, `case eGi`, `id_str = ...FindGi`. Range should be 175-186 (or 177-185).
- Status CONFIRMED (comment only).

### RP-3 (low, CONFIRMED) Off-by-one start line in comment (blast_seqalign.cpp:1484-1490)

- LOSAT: LOSAT/src/algorithm/blastp/blast_engine.rs:6765-6773 quotes `seqalign = s_BlastHSP2SeqAlign(...); } if (seqalign.Empty()) continue;` under `blast_seqalign.cpp:1484-1490`.
- NCBI: `seqalign =` is at line 1485 (1484 is `} else {`), `if (seqalign.Empty()) continue;` at 1490. Should read 1485-1490. Comment only.

### RP-4 (low, SUSPECTED, outside report code) TBLASTN raw score/bit score differs for one long HSP under composition-based statistics

- Not a report-formatting defect; found while testing >100000-bit rows. Engine (compositional matrix adjustment) for a 4998-aa identical match.
- Repro: work/reports/in/bigq.faa (60000 aa random) vs work/reports/in/bigs.fna (5 back-translated records), `tblastn -query bigq.faa -subject bigs.fna -outfmt 6` -> record sbig4 first HSP: NCBI bitscore 10404 (raw 26998), LOSAT 10412 (raw 27020). With only sbig4 as the subject (in/bigs4.fna) NCBI gives 10383 and LOSAT 10404 (so NCBI itself depends on the other subjects), `-comp_based_stats 0` and `-seg no -evalue 1e-5` agree. BLASTP on the protein version of the same data agrees.
- Runs: bgt6, bgt4 (runs/bgt6/{n,l}.out). Listed for the engine auditors; the contract bucket is "silent output difference" but needs separate triage.

## Checked OK (byte-identical stdout, stderr, 2>&1, exit status, unless stated)

1. Description table / alignment headings for protein subjects (x_InitDeflineTable port, cleaned titles): hand-made stress file with 53 titles (multiple brackets, ". [" / ", [" , TPA:/TPA_exp:/TPA_inf:/TPA_asm:/TPA_reasm:/MAG/MULTISPECIES:/UNVERIFIED: prefixes in all positions, 70/71/200 char titles, word runs, trailing , ; . ~ spaces, HTML-like `<b>` `>` `"` `'`, pipes, backslashes, `%`, empty titles, titles that are only a prefix), with ids plain, lcl|, gb|..|, sp|..|..; plus 70 random-title files (60 records each, alphabet of punctuation/prefix tokens, 40 with random ids of 1-12 chars from punctuation) - outfmt 0 all identical; 35 titles of 55..89 chars (ellipsis boundary 68). Titles with HTML character references (&amp; ...), non-ASCII bytes, tab, CR, control characters, `>?` deflines, empty records and `-` residues are explicit LOSAT rejections.
2. Number of descriptions/alignments: e2e_many_subject (300 hits) and a 720-record mutated set (protein and nucleotide) for blastp and tblastn with default (500/250), -max_target_seqs 1,2,3,5,250,251,260,400,500,501,600,700,1000000,2147483647; -max_target_seqs 0/-1/abc/1.5 are rejected by both (syntax-error text differs: approved exception 1).
3. SEG masked query rows lowercase: blastp -seg default/yes/no/"12 2.2 2.5"/"10 1.5 2.0"/"20 3.0 3.3"/"8 1.5 2.0"/"15 2.5 3.0"/"10 1.0 1.2"/"20 3.5 3.7", tblastn (seg default on) same list, with 7 queries (mid, leading, trailing, whole-query, short, multi-region, both ends), subjects with full-length, truncated fragments, indels, both strands, lowercase/mixed-case query input, tblastn -soft_masking true/false/0, -lcase_masking combinations. -seg bad forms ("12 2.2", "x 1 2", double or leading/trailing spaces) fail in both with exit 1.
4. Epilog: -matrix blosum62 / BLOSUM62 / Blosum62 / bLoSuM62 (printed as typed) for blastp and tblastn; threshold 11.5, 12, 1e9, +inf, 1e-5, 100000, 1e6, 1e7, 123456789, 1e15, 1e20, 0.5, 1e300, 1e400, 2147483647, 2147483648, 4294967296 for blastp; tblastn threshold 13/13.0/0.1/1e-3/12345678 (others are explicit LOSAT rejections); window sizes 0,1,2,20,40,100,1000,2147483647 (blastp; tblastn non-40 rejected); gap penalties 11/1 and the explicit rejections for others. `-threshold inf`/`nan`/`-1` are rejected by both (text differs: approved exception 1). 0x10 is an explicit LOSAT rejection.
5. Zero-score HSPs: -evalue 1e300, +inf, 1e999, 1e308, 1e100, 1e50, 1e30, 1e10, 10000 x outfmt 0/6/7 x 6 query/subject pairs (blastp and tblastn, up to 43k lines of output) and extreme-evalue combos with -threshold 6/8/10, -window_size 0/20, -word_size 5, -task blastp-fast: identical.
6. Invalid queries and mixed batches: all-X, O-only (RP-1), `*`-only, tiny queries, lowercase x, U-only, B/Z/J-only, X-run inside a valid query, in all orders (valid/invalid/valid, invalid first, invalid last, all invalid, no-hit queries mixed), outfmt 0/6/7, with and without 2>&1 and with `2>&1 | cat`: footers, "Effective search space used: 0", -1 footers, outfmt 7 "# 0 hits found" presence, Karlin warning placement before each query's report are identical (only the RP-1 line differs). Empty-record and '-' queries are explicit rejections. Subject files with O residues, empty-record subjects and e2e_titles_subject combined with these queries: identical.
7. TBLASTN Sbjct re-translation with B/Z/J/X: e2e_amb_subject.fna and a generated set (36 records, ambiguity rates 3/15/50 %, both strands, ambiguity at codon positions 1-3) with default and `-seg no -comp_based_stats 0`; queries containing B/Z/J (pident in outfmt 6/7): identical (outfmt 0, 6, 7). -db_gencode non-1 differs from NCBI's -subject path (approved exception in AGENTS.md; 24 codes sweep shows the known score differences, not counted as findings).
8. outfmt 6/7 custom fields for BLASTP: 82+48 combinations (std expansion, `-field` deletion of present and absent fields, repeated fields, unknown names ignored, empty/blank spec, multiple spaces, leading/trailing space, case STD/QSEQID (ignored by both), `-std`, `qseqid -qseqid`, all 30 supported names, accession fields on gb|/sp|/lcl|/ref|/gi|/pdb|/emb| ids) against hit-bearing, no-hit and mixed query sets: identical. `delim=`, outfmt 10 and unported fields are explicit rejections.
9. Epilog database block with 1,278,723 letters / 3,200 sequences (comma grouping), long subject path wrapping with spaces: identical. -out file and stdin query: identical. Subjects with multiple HSPs (sum statistics for tblastn, 12 subjects with several domains in both strands/frames): identical for outfmt 0/6/7, -evalue 1000, -seg no, -comp_based_stats 0.
10. A 60000-aa query with total score above 100000 bits (the "Total"/max-score width quirk in the description table, last-row-counted case) for blastp outfmt 0/6/7: identical (tblastn: only the RP-4 score difference).
11. Additional checks, all identical: duplicate subject IDs (adjacent and interleaved, with -max_target_seqs 1..8) for blastp and tblastn outfmt 0/6/7; short queries/subjects (6-100 aa) with both "Method: Composition-based stats." and "Method: Compositional matrix adjust." rows (364 rows in one run), biased-composition subjects; -task blastp-fast, -word_size 5, -threshold 19/22/25.5, -window_size 0/30 epilogs; blastp -max_hsps 1/2/3/100 (tblastn -max_hsps is an explicit rejection); 12 random mixed query files (valid fragments, X runs, O, `*`, low-complexity runs, lowercase) x blastp/tblastn x outfmt 0/7 - the only differences are RP-1 and RP-5; query deflines with odd ids/titles (lcl|, gb||, 80-char id, 150-char title, trailing punctuation) and query-warning labels (`Query_N <id> ..` truncation) for O-containing queries; subject/query FASTA-reader warnings for subjects (protein 50-letter title, tblastn 20-nucleotide-letter title) match in placement.

## Spot-check of NCBI file:line references in session-added comments (item 9)

Verified OK (cited text found at the cited lines): blast_format.cpp:1547-1556, blast_aux.cpp:936-951, showalign.cpp:4305-4318 (range runs 2 lines past the quoted text, harmless), showalign.cpp:2518-2521, alnvec.cpp:116-121, alnvec.cpp:911-918, tabular.cpp:1277-1278, local_db_adapter.cpp:133-135, seqsrc_multiseq.cpp:226-240, blast_aalookup.c:245, blast_aalookup.c:1292, blast_seqalign.cpp:672-674, blast_format.cpp:1450-1452, showalign.cpp:3600-3603, blast_format.cpp:606-611 (606-608 for the quoted text), align_format_util.cpp:1014-1040 (PruneSeqalign), showalign.cpp:2273, local_blast.cpp:177-180,204-208, blast_format.cpp:445-478 and align_format_util.cpp:581-613 (footer/Karlin printing), blast_results.cpp:82-103, blast_setup_cxx.cpp:920-932, create_defline.cpp:3273, 4092-4095, 3431-3446, 3952-3960, 4050-4062, showdefline.cpp:498, format_flags.cpp:219,221, blast_args.cpp:2913-2927, format_flags.cpp:41-195, ncbistr.cpp:4223. Wrong: RP-2, RP-3.
