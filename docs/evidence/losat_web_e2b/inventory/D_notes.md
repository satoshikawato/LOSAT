# S08 inventory, range D: tblastx `-outfmt 7` (and the shared parts of `-outfmt 6`)

Scope: `tabular.cpp` (`CBlastTabularInfo`), the tabular parts of `blast_format.cpp`, the driver in `tblastx_app.cpp`, and the pieces they call that decide the bytes (`align_format_util.cpp`, `blast_seqalign.cpp`, `local_blast.cpp`, `blast_stat.c`, `blast_setup_cxx.cpp`, `blast_results.cpp`). The table is `D.tsv` (60 rows). All NCBI line numbers are CR-stripped and refer to the pinned commit 598d8ae6.

LOSAT references are to the committed tree HEAD 4fc67f9ab (`LOSAT/src`). The working tree had uncommitted edits to `common.rs` and the TBLASTX engine files (`algorithm/tblastx/blast_engine/{mod,run_impl}.rs`, `blast_gapalign.rs`, `chaining.rs`, `sum_stats_linking/linking.rs`, ...) while this inventory was written (a port in progress); the line numbers I give for those files are HEAD's. Every other LOSAT file I cite is unchanged from HEAD. I did not build or run LOSAT (no binary exists and building is out of scope), so every "LOSAT does X" below is from reading the code; the TSV says so where it matters.

## 1. What NCBI prints for `tblastx -query Q -subject S -outfmt 7`

Call path: `CTblastxApp::Run` (`tblastx_app.cpp`) -> `PrintProlog` (nothing for outfmt 6/7, `blast_format.cpp:361`) -> per query batch `GetNextSeqBatch` (batch size 10002 nt, `blast_input_aux.cpp:130-134`; a batch ends with the query that makes the running length reach 10002, `blast_input.cpp:138-170`) -> `CLocalBlast::Run` -> `PrintOneResultSet` once per query -> `x_PrintTabularReport` -> `PrintEpilog` after the last batch. For `-subject` the adapter is in dbscan mode (`blast_app_util.cpp:205-210`, `BL2SEQ_LEGACY` unset), so there is one result per query (not per query/subject pair) and the Database line is printed, not a Subject line.

Block of one query (exact bytes, `cat -A`, from `tblastx -query q1.fa -subject s1.fa -outfmt 7`):

```
# TBLASTX 2.17.0+$
# Query: q1 first query window$
# Database: User specified sequence set (Input: s1.fa)$
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score$
# 176 hits found$
q1^Is1^I100.000^I344^I0^I0^I4412^I3381^I921^I1952^I0.0^I776$
...
```

and after the last query block: `# BLAST processed 3 queries$` (outfmt 7 only; outfmt 6 has no footer).

| Case | Bytes after the Database line | Source |
|---|---|---|
| searched, N >= 1 HSPs | `# Fields: ...` line, `# N hits found`, then N rows | `tabular.cpp:1277-1283` |
| searched, 0 HSPs (also an invalid query that shares its batch with a valid one) | `# 0 hits found` only (no Fields line) | `tabular.cpp:1279-1282`, `blast_seqalign.cpp:1559-1562` |
| unsearched batch (no valid context in the whole batch) | nothing: no Fields, no `# 0 hits found` | `if (align_set)` false, `local_blast.cpp:177-224` |

- `N` counts HSPs: the ungapped translated Seq-align of a subject is split into one Seq-align per HSP (`PrepareBlastUngappedSeqalign`, `showalign.cpp:3163-3215`, called at `blast_format.cpp:763-764`). Checked: the full LC738874 vs LC738875 run gives `# 4319 hits found` and 4319 rows.
- `# Query:` is `TruncateSpaces(idstring + " " + title)` (`align_format_util.cpp:732-734`). Without `-parse_deflines` the Bioseq id is local and its string is empty (`align_format_util.cpp:635`), so the line is the FASTA title: the whole defline (id token included), leading/trailing whitespace removed, internal whitespace kept, no HTML decoding, no `CDeflineGenerator` cleaning, cut at the first control character (tab). An empty defline gives `# Query: ` (trailing space kept) and the row id `Query_1`. The row id (first column) is the first token of the title (`s_ReplaceLocalId`, `tabular.cpp:474-508`).
- Rows are byte-identical to `-outfmt 6` (`cmp` on the full pair, on qq/ss, and with `-max_target_seqs 1`). Subject order and HSP order: hit-list order, then `Blast_HSPListSortByEvalue` per subject (`blast_seqalign.cpp:1571-1577`); LOSAT's `write_output_ncbi_order_evalue_hsp_order_to_writer` already does this (certified outfmt 6).
- Bit score below 10 keeps one leading space (`%4.1lf`): the last field is ` 9.3` (oracle: `-evalue 1000 -threshold 8`).

### Unsearched queries and batches (TBLASTX specifics)

- A batch is "unsearched" when `BlastScoreBlkCheck` fails, i.e. no query of the batch has a context with a computable ungapped Karlin block (`prelim_stage.cpp:321-327`, `blast_stat.c:2815-2823`).
- For translated programs the per-context warning is suppressed (`blast_stat.c:2787`), and only the batch-wide warning (`kBlastMessageNoContext`) is written, so:
  - an invalid query in a batch with a valid query: NO warning, block ends with `# 0 hits found` (the S07 BLASTN rule of one warning per invalid query does not apply to TBLASTX);
  - every query of an all-invalid batch: one stderr line each, block without any count line, exit status 0, outfmt 6 prints nothing at all for them.
- Invalid when alone (oracle, `-evalue 1e-100`): 300 N, `AC` (2 nt), poly-A 300, `ACGT` x 75, IUPAC-only `RYKMSWBDHVN` x 20. Valid: `ACG` (3 nt), 9 nt, 30 nt (searched, `# 0 hits found`). The non-obvious invalid ones are poly-A, the ACGT repeat and IUPAC-only.
- A record with no residues (`>x` only) in a batch with a valid query: warning `Query_2 x: Sequence contains no data `, block `# Query: x` + `# 0 hits found`. Alone in its batch it is fatal: stderr `BLAST engine error: Warning: Sequence contains no data `, exit 3, no output of that batch (earlier batches' blocks are already on stdout, no footer). With an all-N query in the same batch both are unsearched and the empty one's warning ends `...filtering options Sequence contains no data `.
- Warning text (one per invalid query of an unsearched batch): `Warning: [tblastx] Query_<n> <title>: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options ` + `\n`. `<n>` is the query's position in the whole file; the label is cut to 25 bytes + `.. ` when longer than 35 bytes (`blast_setup_cxx.cpp:534-543`, e.g. `Query_1 qLongName this is.. : `).

### Order of writes and flushes (`2>&1`)

`~CBlastTabularInfo` flushes stdout (`tabular.cpp:160-163`) at the end of each query's block, and `cerr` is tied to `cout`. A query's warning (`PrintOneResultSet`, `blast_format.cpp:1450-1452`) is therefore written after everything before it and before that query's header. Oracle `big_then_N.fa` (LC738874 full genome, then a 300-N query): lines 1-2449 are query 1's block, line 2450 is the warning for query 2, then query 2's header and `# BLAST processed 2 queries`. With two all-N queries in one later batch (`mb4.fa`) the warnings alternate with the headers: header of q3 ... `Warning Query_4`, header Query 4, `Warning Query_5`, header Query 5, footer. Option warnings at setup (`Examining 5 or more matches is recommended` for `-max_target_seqs` below 5) come before any stdout.

### Exit status and special inputs

| Input | stdout | stderr | rc |
|---|---|---|---|
| empty or whitespace-only query file | empty (no footer) | `Warning: [tblastx] Query is Empty!` | 0 |
| `-query -` with empty stdin (out of scope, noted) | `# BLAST processed 0 queries` | none | 0 |
| all-N / 2-nt query alone | header-only block + footer | warning above | 0 |
| `>x` only, alone | empty | `BLAST engine error: Warning: Sequence contains no data ` | 3 |
| write to a full device (`> /dev/full`, `-out /dev/full`), outfmt 6 or 7 | empty | `terminate called after throwing an instance of 'std::__ios_failure'` / `  what():  basic_ios::clear: iostream error` | 134 |

## 2. LOSAT: what exists and which writer TBLASTX can use

- TBLASTX accepts only `-outfmt 6` (`blastinput/value_parsers.rs:670`, `algorithm/tblastx/args.rs:90`). `write_tblastx_outputs` (`run_impl.rs:687`) always calls `write_output_ncbi_order_evalue_hsp_order_to_writer` (`common.rs:705`), which sorts (query, subject, HSP) in NCBI order and writes outfmt 6 rows with `write_hit_fields` (`report/outfmt6.rs:520`). That sorting and the row writer are reusable and certified.
- The generic outfmt 7 code is not usable for TBLASTX, and no program uses it any more:
  - `write_outfmt7_header` (`report/outfmt6.rs:659`): `# Fields: qaccver, saccver, pident, ...` (NCBI: `query acc.ver, subject acc.ver, % identity, ...`), program/version from the `ReportContext` (`BLAST 0.1.0` by default), Query/Database only when set, `# N hits found` always (so `# 0 hits found` for an unsearched query: the S07 audit item N1 is still true of this writer).
  - `write_outfmt7` (`:705`): one header for the whole hit list; `write_outfmt7_grouped` (`:779`): one block per query that has hits, titled with the query id (first token), no block for queries without hits, no footer.
  - `common.rs` `TabularWithComments` branches (`:630-633`, `:802-808`, `:997-1001`, `:1022-1031`): collect all queries' sorted hits and call `write_outfmt7` once; for no hits only a header.
- BLASTN's outfmt 7 (it is in `algorithm/blastn/hsp.rs`, not in `run.rs`): `write_output_blastn_hitlists_to_writer` (`:1143`) iterates all queries, calls the warnings object before each query, writes the header with `write_blastn_outfmt7_header` (`:179-198`), and the footer (`:1286-1288`) when `epilog` is true. The header takes `num_hits: Option<usize>`; `None` (unsearched query) prints only the three lines, so BLASTN no longer has item N1: it prints `# 0 hits found` only for searched queries. `unit test` at `hsp.rs:1595-1650`. It is private and hard-codes `# BLASTN 2.17.0+`.
- TBLASTN (`algorithm/tblastn/stage_e_report.rs:136-266`) and BLASTX (`algorithm/blastx/report.rs:1123-1330`, `:1629-1665`) have the same structure with their own literals; both use `search_skipped` for the missing count line, and BLASTX also has the `Sequence contains no data` warning (`:1762`), `TabularFlushError` (abort on write failure) and a nucleotide FASTA reader (`algorithm/blastx/input.rs:107`) that returns `title` and `internal_id`.
- Reusable helpers: `blastinput/query_batch.rs:19` (`query_batches(&lengths, 10_002)`), `report/query_warnings.rs` (`QueryWarnings::before_query` flushes the report, then writes the warning; `invalid_query_warning`; `few_matches_warning`).

Answer to "which writer": TBLASTX should get its own per-query writer shaped like `write_output_blastn_hitlists_to_writer` / TBLASTN's `write_tabular`, built on a generalisation of `write_blastn_outfmt7_header` (program label + version as parameters, `pub(crate)`; BLASTN passes its constants, so BLASTN bytes are unchanged), reusing the NCBI-ordered (query, subject, HSP) grouping of `common.rs` for the rows. Do not call `write_outfmt7*` or the `TabularWithComments` branches.

## 3. What the port must do (each item with the NCBI line it rests on)

1. Dispatch on the requested format in `write_tblastx_outputs` and let `tblastx_outfmt` accept `7` (still no custom fields): `blast_format.cpp:1454-1459`. The rows stay `write_hit_fields` rows (identical to outfmt 6, `tabular.cpp:1100-1108`).
2. Write one block per query for ALL queries of the file, in input order, also for queries without hits and unsearched ones (outfmt 6 prints nothing for them): `tblastx_app.cpp:203-205`, `blast_format.cpp:1430`.
3. Block lines, exactly: `# TBLASTX 2.17.0+`, `# Query: <title>`, `# Database: User specified sequence set (Input: <-subject text as typed>)`, then (a) unsearched: nothing; (b) 0 HSPs: `# 0 hits found`; (c) N > 0: `# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score` and `# N hits found`, N = HSPs printed (`tabular.cpp:1277-1283, 1111-1236, 1287-1320`; `blast_format.cpp:795-808`; `local_blast.cpp:177-224`).
4. Footer `# BLAST processed N queries` for outfmt 7 only, N = all queries of all batches, not written when a later batch fails and not written after "Query is Empty!" (`blast_format.cpp:2233-2238`, `tabular.cpp:1322-1325`, `tblastx_app.cpp:132-135, 209`).
5. Replace `anyhow::ensure!(any ctx.is_valid)` (`run_impl.rs:1266-1269`) by per-batch validity: batches from `query_batches(lengths, 10_002)`; a batch is unsearched when no query has a valid context; unsearched queries get the header-only block, the warning, and exit 0 (`blast_stat.c:2815-2823`, `prelim_stage.cpp:321-327`, `local_blast.cpp:177-224`). Validity must follow NCBI for poly-A, ACGT repeats, IUPAC-only and < 3 nt (see the UNSURE row).
6. Warnings only for an unsearched batch (every query of it) and for records without residues; none for an invalid query in a mixed batch; written through `QueryWarnings` so stdout is flushed first (`blast_stat.c:2787, 2819`, `blast_setup_cxx.cpp:534-543, 608-640`, `blast_results.cpp:276-293`, `tabular.cpp:160-163`).
7. `# Query:` text = reader title (`align_format_util.cpp:732-739`). TBLASTX reads queries with `bio` (`run_impl.rs:558`), which has no title: use `blastx::input::read_fasta(path, false, ...)` (`title`) or BLASTN's `fasta_defline` together with its control-character/leading-space rejection (`check_deflines`). Otherwise tab-separated deflines, empty deflines and `>` followed by spaces differ (rows too: `Query_1` vs `unknown`).
8. Apply `-max_target_seqs` in the TBLASTX engine (currently read by nobody): `# N hits found` and the rows count only the best N subjects (`blast_format.cpp:93, 813`, engine hit list), plus the `Examining 5 or more matches is recommended` warning for N < 5.
9. `Query is Empty!` for blank query files (`tblastx_app.cpp:132-135`, `blast_app_util.cpp:846-875`).
10. Write failure: abort with the `std::__ios_failure` text, exit 134 (as BLASTX does): `blast_format.cpp:119`, `tabular.cpp:160-163`.
11. Keep the `FormatProbe` hooks around each row (`probe.begin/end`) as BLASTN/TBLASTN do.

## 4. Oracle runs (all under `scratch_D/`, NCBI 2.17.0+ at `/home/kawato/micromamba/bin/tblastx`)

Inputs were cut from `LC738874.fasta` / `LC738875.fasta` (`mk.py`): `q1.fa` (4500 nt of LC738874 around 236000, hits against `s1.fa` = 4500 nt of LC738875), `q2.fa`/`s2.fa` (3000 nt), `qq.fa` = q1,q2,q3(1500 nt, few random hits), `ss.fa` = s1,s2, `qN.fa` (300 N), `qR.fa` (random 300 nt, no hits at `-evalue 1e-100`), plus the mixed files named in the TSV notes. Runner `run.sh`; results are the `*.out`, `*.err`, `*.both`, `*.rc` files.

| Run | Observation |
|---|---|
| full pair, `-outfmt 7` and `-outfmt 6` | `# 4319 hits found`, 4319 rows, rows `cmp`-identical, `# BLAST processed 1 queries`, rc 0, ~3.8 s |
| `qq.fa` vs `ss.fa` | 3 blocks (153, 22, 7 hits), q1 rows: 152 for s1 then 1 for s2; footer 3 |
| same, `-evalue 1e-100` | q1 block 42 hits; q2 and q3 `# 0 hits found` (no Fields line) |
| same, `-max_target_seqs 1` | 152 / 18 / 4 hits (only s1 kept), stderr `Warning: [tblastx] Examining 5 or more matches is recommended`, rows identical to outfmt 6 |
| `qN.fa` alone / `q2nt` / polyA / ACGT x 75 / IUPAC-only | header-only block + footer, warning on stderr, rc 0; outfmt 6 prints nothing |
| `q1N2`, `qN1N2`, `q1N`, `q1_2nt`, `q2nt_q1`, `q1_q2nt` | invalid query between/after valid ones: `# 0 hits found`, no warning, footer counts all |
| `a3`, `a9`, `a30` (3, 9, 30 nt) | valid: `# 0 hits found` |
| `mb1..mb5` (batches at 10002 nt) | all-N alone in its batch -> header-only; with a valid query in its batch -> `# 0 hits found`; two all-N queries in one later batch -> two warnings, two header-only blocks |
| `big_then_N.fa` with `2>&1` | warning between the blocks (section 1) |
| `e1.fa` / `e2.fa` / `donly.fa` / `emN.fa` / `Nem.fa` | record without sequence: see section 1 |
| `empty.fa`, `nl.fa`, `ws.fa` | `Query is Empty!`, no stdout, rc 0 |
| `printf '' | tblastx -query -` | `# BLAST processed 0 queries` |
| `/dev/full` (stdout and `-out`) | rc 134, terminate message |
| `-out o.txt` | same bytes in the file, stdout and stderr empty, rc 0 |
| `-num_threads 3` | stderr `Warning: [tblastx] 'num_threads' is currently ignored when 'subject' is specified.`, output unchanged |
| `-subject ./ss.fa` | `Input: ./ss.fa` |
| 19 deflines (`defl/d*.fa`) | `# Query:` and first column per section 1 |

## 5. Open questions / UNSURE

- UNSURE: whether LOSAT's per-context `is_valid` (`algorithm/tblastx/lookup/backbone.rs:384`, `compute_karlin_params_ungapped`) is false for exactly poly-A, the ACGT repeat, IUPAC-only and < 3 nt and true for 3, 9 and 30 nt, as NCBI is. Check with the six one-file oracle inputs of section 1 once the writer exists.
- UNSURE: tie behaviour when `-max_target_seqs` cuts between subjects with equal best e-value (the oracle only covers a clear winner). It belongs to the engine range.
- UNSURE: the NCBI exit status/text when the empty-sequence record is the first record of a later batch together with other invalid queries (only the single-record and the paired cases were run).
- The ABI v1/v2 `run_local` takes `bio` records (`run_impl.rs:657`), so the `# Query:` title must come from the record's id and description or the caller must pass titles; this is a design choice for the port.
- Several output formats in one run (`ReportOutputs.formats`): TBLASTN's comment says the warnings are written once; the same should hold for TBLASTX.
