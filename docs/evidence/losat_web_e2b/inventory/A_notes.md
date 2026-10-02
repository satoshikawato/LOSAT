# Range A notes: tblastx application layer and CBlastFormat orchestration (outfmt 0 / 6 / 7)

Rows: `A.tsv` (one row per NCBI function/branch). Oracle inputs and outputs: `scratch_A/` (`NAME.out`, `NAME.err`, `NAME.rc` = NCBI 2.17.0+; `NAME.lout/.lerr/.lrc` = LOSAT-base). Line numbers of NCBI files are CR-stripped (`tr -d '\r' | cat -n`), pinned commit 598d8ae6.

## 0. Method and a warning about the LOSAT tree

- NCBI oracle: `/home/kawato/micromamba/bin/tblastx` 2.17.0+ (build Aug 11 2025). No `.ncbirc`; `BATCH_SIZE`, `CHUNK_SIZE`, `BL2SEQ_LEGACY`, `CTOOLKIT_COMPATIBLE`, `OLD_FSC` unset except in the single commands that say otherwise (BATCH_SIZE in four commands, to study batching).
- LOSAT comparison binary: `s08/bin/LOSAT-base`. LOSAT line numbers in the TSV are those of the committed tree `4fc67f9ab` (extracted read-only to `scratch_A/base/LOSAT/src/`; `scratch_A/checkloc.py` re-checks every cited location against it). While I worked, somebody else modified the LOSAT working tree (uncommitted changes to `report/pairwise.rs`, `common.rs`, `algorithm/tblastx/*`; line numbers there shifted by up to 57). Re-check the locations against the tree you edit.
- Inputs cut from LC738874/LC738875 (`mk.py`): `q1.fa` (3000 nt, id `qwin0`), `q2.fa`, `s1.fa` (3000 nt), `s2.fa` (`swin0`, `swin40k`), `allN.fa` (300 N), `short2.fa` (`AC`), `mixed.fa` (valid, all-N, valid), `two_allN.fa`, `mix_empty.fa`, `b_seq1.fa`, `b_seq2.fa`, `b_seq3.fa` (10003-nt queries to force several batches), `empty.fa`, `ws.fa`, `hdronly.fa`, `emptytitle.fa`, `big_tiny.fa`, ...
- Helper scripts: `run.sh NAME args` (NCBI), `runl.sh NAME args` (LOSAT), `flushmap.py args` (gdb `catch syscall write`: lists every write(2) of NCBI tblastx with its size and first/last bytes), `addrows.py` (row appender).

## 1. Order of everything NCBI does for `tblastx -query Q -subject S -outfmt F` (what the port must reproduce)

1. AppMain parses the command line (`tblastx_app.cpp:79-89`). Syntax errors: USAGE + `Error: ...` on stderr, exit 1 (LOSAT: approved exception, clap, exit 2).
2. `Run` sets the diagnostics: warnings are `Warning: [tblastx] <text>\n` on stderr (`:93-101`).
3. `SetOptions` (`:110`, `blast_args.cpp:3585-3641`) runs ExtractAlgorithmOptions of every argument group in registration order (`tblastx_args.cpp:44-119`): database/subject args (the `-subject` FASTA is read and its adapter built here), std args (`-query` and `-out` are opened here; `-out` is created/truncated), generic search args, ..., formatting args (outfmt parse errors; `Examining 5 or more matches is recommended`), MT args (`num_threads` warnings), then `Validate()` (`BLAST query/options error: ...`, exit 1). So stderr order is: formatting warning, then thread warnings(s), then `Query is Empty!`.
4. `InitializeSubject` (`:118`; for -subject the reading already happened) - subject problems (`Empty CBlastQueryVector`, missing file) beat query problems.
5. `IsIStreamEmpty(query)` (`:132-135`): empty or whitespace-only seekable file -> `Warning: [tblastx] Query is Empty!`, exit 0, nothing on stdout, no epilog, no `# BLAST processed` line. (A pipe is never "empty": the run goes on and ends with `# BLAST processed 0 queries` for outfmt 7.)
6. `CBlastFormat` constructed (`:145-163`): report stream gets `exceptions(badbit)` (`blast_format.cpp:119`); dbscan subject info (`:129-139`).
7. `PrintProlog()` (`:173`): outfmt 6/7 write nothing. outfmt 0 writes `TBLASTX 2.17.0+\n` + `endl` (flush, write #1 = 17 bytes `TBLASTX 2.17.0+\n\n`) + `endl` (flush, 1 byte), then the Gapped BLAST reference (4 wrapped lines + blank), `\n\n`, the `Database:` block - the last part stays in the stream buffer. This is before the first query is read.
8. Loop over query batches (`:176-207`): batch = queries added until 10002 nucleotides are reached (`blast_input_aux.cpp:130-134`, `blast_input.cpp:138-171`); FASTA-reader diagnostics for the batch appear while it is read; `CLocalBlast::Run` searches the batch (`local_blast.cpp:166-299`) or, when no query in the batch has any context with Karlin parameters (`BlastScoreBlkCheck`, `blast_stat.c:853-877`), returns the not-searched results (`:177-225`: NULL Seq-align, ancillary data -1/-1/-1, search space 0, warnings for every query of the batch). Then one `PrintOneResultSet` per query.
9. `PrintOneResultSet` (`blast_format.cpp:1410-1590`): `m_QueriesFormatted++` -> (errors) -> `ERR_POST(Warning << query warnings)` (stderr; cerr is tied to cout so pending stdout is flushed first) -> tabular branch (`x_PrintTabularReport`) or pairwise branch (preamble `\n\nQuery= ...\nLength=N\n`, no hits / description table + `\n` + alignments, footer `x_PrintOneQueryFooter`).
10. `PrintEpilog` (`:209`): outfmt 0 `\n\n` (endl endl) + Database block + `\n\nMatrix: ...` + threshold + window; outfmt 7 `# BLAST processed N queries`; outfmt 6 nothing.
11. `~CBlastFormat` flushes. Exit status = CATCH_ALL status (`blast_app_util.hpp:167-267`).

### Observed write(2) map (gdb, `flushmap.py`)

- outfmt 0, one query (q1 vs s1): `[17] "TBLASTX 2.17.0+\n\n"`, `[1] "\n"`, three `[4096]` buffer-full writes, `[920]` (up to and including the last alignment's blank lines), `[1] "\n"` (footer `NcbiEndl`), `[90]` (Lambda/K/H block + `Effective search space used: 948676\n` + first epilog `NcbiEndl`'s `\n`), `[1] "\n"`, `[240]` (database block + Matrix + threshold + window; written at exit). I could not find which statement flushes the 920 bytes before the footer (no flush in showalign.cpp/showdefline.cpp); the exact position only matters for the partial output left behind by a failing write.
- outfmt 7: one write per query block (flush in `~CBlastTabularInfo`, `tabular.cpp:160-163`) plus `[28]` for the `# BLAST processed` line: `mixed.fa` gives `[4096][1434] | [108] | [3025] | [28]`.
- outfmt 6: nothing is written until a query's tabular block ends (`~CBlastTabularInfo` flush) - one write per query.

### Interleaving seen with `2>&1` (saved: `merged_bseq2_0.txt`, `merged_bseq2_7.txt`)
`b_seq2.fa` = [all-N 10003 nt] [valid 3000 nt, all-N 10003 nt] [all-N 300 nt]: batches are {Q1} {Q2,Q3} {Q4}. Warnings only for Q1 and Q4 (their batches have no valid context); Q3 (all-N inside a searched batch) is silent. outfmt 0: the warning for Q1 comes after the first 13 lines of the prolog (the database block's trailing blank line is flushed by the cerr tie) and before `\n\nQuery= allN first 10003`; the warning for Q4 comes right after Q3's `Effective search space used: 0` line.

## 2. Oracle results (bytes and exit status)

| Case | stdout | stderr | rc |
|---|---|---|---|
| empty query (empty.fa, ws.fa), outfmt 0/6/7 | empty | `Warning: [tblastx] Query is Empty!\n` | 0 |
| `-max_target_seqs 1` + empty query | empty | `Warning: [tblastx] Examining 5 or more matches is recommended\nWarning: [tblastx] Query is Empty!\n` | 0 |
| query `>hdr only` | outfmt 0: whole prolog (379 B); 6/7: empty | `BLAST engine error: Warning: Sequence contains no data \n` | 3 |
| 2-nt query / 1-nt / 300 N alone, outfmt 0 | prolog, `Query= ...`, `Length=N`, `***** No hits found *****`, footer with `-1.00` blocks and `Gapped` heading, `Effective search space used: 0`, epilog | `Warning: [tblastx] Query_1 <title>: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n` | 0 |
| same, outfmt 7 | `# TBLASTX 2.17.0+`, `# Query: <title>`, `# Database: ...` and `# BLAST processed 1 queries` (no `# 0 hits found`) | same warning | 0 |
| same, outfmt 6 | empty | same warning | 0 |
| all-N / 2-nt query inside a batch with a valid query | outfmt 0: `No hits found`, 5 blank lines, `Effective search space used: 0` (no Karlin block); outfmt 7: header + `# 0 hits found` | none | 0 |
| 3-nt query (`ATG`) | valid search, `# 0 hits found`, Lambda 0.318 K 0.134 H 0.401, search space 1000 | none | 0 |
| `-max_target_seqs 3` / `1` | normal; subjects limited to N | `Warning: [tblastx] Examining 5 or more matches is recommended\n` (once, all formats) | 0 |
| `-evalue 1e-100` (no hits) | outfmt 0: `No hits found`, 3 blank lines, valid Lambda block, search space; outfmt 7: `# 0 hits found`, `# BLAST processed 1 queries`; outfmt 6: empty | none | 0 |
| `-out /dev/full`, outfmt 0 | empty | `BLAST failed to write output\n` (nothing else, even with an all-N query) | 6 |
| `-out /dev/full`, outfmt 6 / 7 (also stdout to /dev/full) | empty | `terminate called after throwing an instance of 'std::__ios_failure'\n  what():  basic_ios::clear: iostream error` | 134 (SIGABRT, core) |
| `-outfmt "7 std"`, `"6 std"`, `"0 std"` | byte-identical to `7`, `6`, `0` | none | 0 |
| `-outfmt abc` / `6x` / `6.0` / `''` | empty | `BLAST query/options error: '<x>' is not a valid output format\nPlease refer to the BLAST+ user manual.\n` | 1 |
| `-outfmt 22` / `-1` | empty | `Error: Formatting choice is out of range` | 255 |
| `-outfmt 17` / `19` / `21` | empty | `BLAST query/options error: SAM format is only applicable to blastn` / `AIRR rearrangement format is only applicable to igblastn` / `FASTA output format is only applicable to magicblast` (+ Please refer line) | 1 |
| `-outfmt ' 7'`, `'7 '`, `'+6'`, `'06'` | same as `7`, `7`, `6`, `6` | none | 0 |
| `-threshold 11.5` / `12.3456789` / `1e2` | epilog `Neighboring words threshold: 11.5` / `12.3457` / `100` | none | 0 |
| `-threshold 0` | empty | `BLAST query/options error: Non-zero threshold required` | 1 |
| `-window_size 0` | epilog has no `Window for multiple hits` line | none | 0 |
| `-num_threads 2` / `1000` with -subject | normal | `'num_threads' is currently ignored when 'subject' is specified.` / plus `Number of threads was reduced to 32 to match the number of available CPUs` (first) | 0 |
| `-query nonexist.fa` / `-subject nonexist.fa` / `-out /nonexistent/d/x` | empty | `Command line argument error: Argument "query" (resp. "subject", "out"). File is not accessible:  `<name>'` | 1 |
| no `-subject` | empty | `BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified` + Please refer line | 1 |
| empty `-subject` file | empty | `BLAST engine error: Empty CBlastQueryVector` | 3 |
| query with binary bytes | empty | `BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` | 1 |
| `-out -` | report on stdout | none | 0 |
| title `abcdefghijklmnopqrstuvwxyza` (27 chars) vs `...zab` (28) in an invalid batch | warning label `Query_1 abcdefghijklmnopqrstuvwxyza:` vs `Query_1 abcdefghijklmnopq.. :` | | 0 |
| 300 subjects with hits (`s300.fa`, `q_small.fa`), outfmt 0 default / `-max_target_seqs 280` | 300 description rows but alignments for 250 subjects / 280 and 280 | none (280 case: Examining warning not printed) | 0 |
| same, outfmt 6 default / 280 | 300 / 280 subjects | none | 0 |
| `-db_gencode 2` | search unchanged; outfmt 6 identity/mismatch and outfmt 0 Sbjct letters use code 2 | none | 0 |
| `-query_gencode 2` | search and display change | none | 0 |

## 3. LOSAT-base differences found while comparing (all TBLASTX)
- outfmt 0 and 7 are rejected by the clap value parser (`blastinput/value_parsers.rs:670`, default `-outfmt` is `0`, so a plain run fails with exit 2).
- No diagnostics channel: no `Query is Empty!`, no `Examining 5 or more...`, no invalid-query warnings; errors are plain `Error: ...` exit 1 (`main.rs:120`), not NCBI's CATCH_ALL texts/exit codes.
- All queries are searched in one call. With at least one valid query, invalid queries are skipped silently and the outfmt 6 bytes equal NCBI's (`mixed.fa`, `mixed_short.fa`, `mix_empty.fa`: identical). With no valid query in the whole input LOSAT-base fails (`Error: no valid query contexts`, rc 1; `two_allN.fa`, `allN.fa`, `short2.fa`) where NCBI exits 0 with the warnings and the not-searched reports; an all-N batch followed by valid queries (`b_seq2.fa`) is likewise not reproduced (NCBI warns for the batch that has no valid context).
- `-max_target_seqs` is parsed but never used by the TBLASTX engine or writer (`algorithm/tblastx/args.rs:57`): with two subjects with hits and `-max_target_seqs 1`, NCBI outputs 37 rows (best subject only) and LOSAT-base 95 rows (`m1_6.out` vs `m1_6.lout`).
- Query order/batches change NCBI's hits (see section 4); LOSAT-base always behaves like a single batch.
- A query with an empty or blank defline (`>` / `>   `) is `Query_1` in NCBI's tabular rows; LOSAT-base prints `unknown` (`emptytitle6.out` vs `.lout`).
- `-out -` creates a file named `-` (NCBI: stdout); `-out FILE` is not created for an empty query (NCBI creates/truncates it before reading anything); an empty `-subject` file is accepted (NCBI: rc 3).
- Approved differences (AGENTS.md): clap errors/`-help` (exit 2), thread warnings with `-subject`, non-zero exit on an outfmt 6/7 write failure, `-db_gencode` applied to the search.

## 4. Batches are not only about warnings: they change hits
NCBI default (10002-nt batches) vs one batch (`BATCH_SIZE=1000000`) for `b_seq1_valid.fa` (a 10003-nt query then a 3000-nt query) against `s2.fa`, outfmt 6: the rows of the 10003-nt query differ (e.g. its `swin40k` rows: `0.23`/`0.23` in the default run, `0.30`/`0.34` plus extra HSPs such as `7253 7212 2995 2954 0.34 24.0` in the one-batch run). With the queries in the opposite order they share a batch and NCBI equals LOSAT-base. LOSAT-base (which searches all queries together) also gives different rows for the long query when another query is present (`bs1v.lout` vs `bigonly.lout`). A tiny 120-nt companion query already changes the long query's rows in a shared batch; a 6000-nt companion did not (`big_tiny.fa`, `big_mid.fa`), so it is an engine detail of the concatenated query block (not query splitting: `split_query_cxx.cpp:60-61` disables splitting for ungapped searches). Queries below 10002 nt are batch-independent (`two.fa` with `BATCH_SIZE=1` is byte-identical). Cause not located (engine range). The port must therefore search batch by batch with NCBI's batch boundaries, not only print per batch.

## 5. What the port must do (with NCBI file:line)
1. Batching (`blast_input_aux.cpp:130-134`, `blast_input.cpp:138-171`, `tblastx_args.cpp:128-133`): tblastx batches hold 10002 nucleotides; the query that reaches/crosses the limit stays in the batch; no `CBatchSizeMixer`. Search and report batch by batch (hits depend on the batch, section 4). Numbering of local ids `Query_<n>` continues across batches.
2. Not-searched batch (`local_blast.cpp:177-225`, `blast_stat.c:2778-2823`, `prelim_stage.cpp:321-327`): a batch where no query has a context with Karlin parameters is not searched; every query gets the warning (label `Query_<n> <title>`, >35 bytes cut to 25 + `.. `; text ends with a space) and a report with `No hits found` and the `-1.00` blocks + `Gapped` heading; outfmt 7 prints no `# N hits found` line for such queries; outfmt 6 prints nothing. An invalid query inside a searched batch is silent, has `# 0 hits found` and a footer without Karlin blocks (`blast_format.cpp:451-478`).
3. Write order and flushes (`blast_format.cpp:395, 456, 2249`, `tabular.cpp:160-163`, `blast_app_util.hpp:252-255`): flush the prolog's first two newlines before reading the first batch; flush stdout before every stderr warning; failing write: outfmt 0 -> `BLAST failed to write output`, exit 6, stopping at the prolog flush (nothing else printed, not even query warnings); outfmt 6/7 -> abort-like non-zero exit (approved).
4. Option/formatting semantics (`blast_args.cpp:2794-2991`): `-outfmt` token parse and error texts/exit codes; `Examining 5 or more matches is recommended` once at option extraction when -max_target_seqs < 5; `-max_target_seqs` omitted (outfmt 0: 500 descriptions, 250 alignments, hit list 500) vs given (N, N, N); outfmt 6/7 prune to the hit-list size by subjects (`blast_format.cpp:813`).
5. Empty query, subject-before-query order, `-out` creation order and the CATCH_ALL mapping (`tblastx_app.cpp:132-135`, `blast_app_util.hpp:167-267`): `Query is Empty!` exit 0; `Command line argument error:` / `BLAST query/options error:` + `Please refer...` / `BLAST engine error:` / `Error:` with exit 1/1/3/255.
6. outfmt 7 (`blast_format.cpp:794-809`, `tabular.cpp:1266-1325`): per query `# TBLASTX 2.17.0+`, `# Query: <defline>`, `# Database: User specified sequence set (Input: <-subject text>)`, `# Fields: ...` only when hits exist, `# N hits found` only for searched queries, final `# BLAST processed N queries` (N = every query reaching PrintOneResultSet); outfmt 6 has neither header nor epilog.
7. outfmt 0 layout (`blast_format.cpp:386-441, 1491-1589, 2249-2288`): prolog/epilog with the nucleotide totals (`3,000 total letters`, `N sequences`), the `N` column header (sum statistics), `s_SetFlags` with `eTranslateNucToNucAlignment`, character middle line, protein alignment type, master/slave genetic codes, footer of the first valid context, epilog `Matrix`, `Neighboring words threshold` (C++ `%g`-like, 6 significant digits) and `Window for multiple hits` (omitted for 0); no `Gap Penalties` line.

## 6. Open questions / UNSURE
- `CBlastFormat::PrintOneResultSet` error branch (`blast_format.cpp:1446-1449`): no error-severity query message found for tblastx with in-scope options.
- `CIOException` other than eFlush (`blast_app_util.hpp:242-247`) leaves the status 0 silently; no trigger found.
- Origin of the 920-byte flush before the query footer (section 1).
- Cause of the batch dependence of hits (section 4).
- `PRE_FETCH_SEQS_LIMIT` warning `Error pre-fetching sequence data` (`blast_app_util.cpp:778-785`) not reachable without the environment variable.
- CBatchSizeMixer/`BATCH_SIZE` are out of scope for TBLASTX; `BATCH_SIZE=abc` makes NCBI die with an uncaught `CStringException`.
