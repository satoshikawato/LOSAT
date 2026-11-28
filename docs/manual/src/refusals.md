# What LOSAT refuses

LOSAT stops with an error for each option and value below. It never runs something else in its place. Use NCBI BLAST+ for these cases, except where an item says that NCBI BLAST+ fails there too.

There are five reasons for a refusal:

- NCBI BLAST+ 2.17.0 crashes, aborts, never finishes, or reads memory outside its buffers, so no output exists to match.
- NCBI's result depends on memory contents, an uninitialised value or the machine.
- The option needs NCBI's servers.
- The option is not implemented yet.
- The option is not part of what v0.2.0 provides.

A refusal comes before the first byte of output unless an item says otherwise. With `-outfmt` 0 to 4, every program writes the start of the report (the program and version line, the references and the database lines) before it reads the first query batch (with `-in_pssm`, after it reads the PSSM), as NCBI BLAST+ does, so a refusal found in a query batch, in the search or in the report comes after it. In several items LOSAT first writes what NCBI writes before it fails (apart from the version line, which names LOSAT), and then stops with the error.

## What you see

LOSAT names the option and the program. There are three forms.

An option refused by name, after NCBI's argument reader has accepted the command line (so NCBI's own argument errors come first, as in NCBI; exit status 1):

```text
Error: the NCBI BLAST+ option -remote is not supported by LOSAT's BLASTP
```

A value or setting refused after the options are read (exit status 1):

```text
Error: -use_sw_tback (the Smith-Waterman traceback) is not supported by LOSAT's PSIBLAST
```

An output format the program does not write (exit status 1):

```text
Error: output format 15 is not supported by LOSAT's PSIBLAST yet
```

For example, NCBI BLAST+ runs `psiblast` with the Smith-Waterman traceback, and LOSAT stops with the second message above:

<!-- example: refusals-1 -->
```bash
losat psiblast -query query.faa -subject proteins.faa -use_sw_tback
```

## How this page is organised

There is one section for each program, and one for the five search programs together: `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`. Common refusals appear directly in the tables below; rare and pathological inputs (where NCBI BLAST+ itself crashes, loops or reads memory it should not) are collected separately in [Rare and pathological inputs](refusals-edge-cases.md).

### At a glance: common CLI refusals

| Option or feature | Programs affected | Status in v0.2.0 |
| :--- | :--- | :--- |
| `-remote`, `-entrez_query` | All search programs | Out of scope (needs NCBI servers) |
| `-import_search_strategy` | `psiblast`, `rpsblast`, `rpstblastn`; strategy files in binary or XML form | Not yet supported (text strategy files run in the five search programs) |
| `-use_sw_tback` | `psiblast`, `rpstblastn` (with CBS), `blastp`, `blastx` and `tblastn` with `-ungapped` (at the first query batch that is searched), `tblastn -in_pssm` (when a match is found) | Refused |
| `-ungapped` | `rpstblastn` (without CBS), `tblastn -in_pssm` | Refused (`rpstblastn` with CBS gives NCBI error) |
| `-use_index true` | `blastn` | Index searches run, with indexes of only some volumes, two-hit (`-window_size` above 0) and `-db_hard_mask` too; refused for damaged or mismatched index volumes |
| `-phi_pattern` with more than one round, non-listed `-outfmt` | `psiblast` | Not yet supported (listed: `0` to `12`, `18` and `20`, with a FASTA query, `-in_pssm` or `-in_msa`; with `-phi_pattern`: `0`, `6`, `7`, `10` and `20`) |
| `-mk_counts -infmt seqids` | `windowmasker` | Out of scope (needs NCBI servers) |
| `deltablast` | Standalone executable | Not yet provided |

## All five search programs

These refusals apply to `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`.

| Option or value | What happens |
| :--- | :--- |
| `-import_search_strategy` with a strategy file in binary or XML form | Refused: not yet supported. Strategy files in text form, the form that `-export_search_strategy` writes, work. |
| `-import_search_strategy` with queries given as a list of locations, a GI list, a PHI-BLAST pattern, or a PSI-BLAST, RPS-BLAST or PSSM search | Refused: not yet supported. |
| `-import_search_strategy` with an option that the program does not have or refuses on its command line | Refused (for example composition-based statistics in `blastn`). The strategy's options that the program's command line runs run as there, the database mask, `-subject_besthit`, `-filtering_db` and the MegaBLAST index included. |
| `-import_search_strategy` with a MegaBLAST database index (`-use_index true`, in the strategy or on the command line) outside a `blastn` search of a database | Refused: not yet supported. |
| `-import_search_strategy` with the strategy of another program | Refused, except `blastn` and `blastx` running a `tblastx` strategy and `tblastx` running a `blastn` or `blastx` one, which run as NCBI runs them. Between the other programs NCBI stops with an error or reads the query as the wrong molecule. |
| `-import_search_strategy` with a negative `WindowSize` (a hand-edited strategy: `-window_size` takes 0 or more) | Runs as NCBI runs it, in the five programs: the one-hit search over NCBI's diagonal table, which passes from one subject to the next, so LOSAT writes NCBI's `-num_threads 1` result at every thread count, and the refusals of `-window_size 0` apply to it too. |
| `-import_search_strategy` with a value outside the range of the command line's option (a hand-edited strategy) | Runs as NCBI runs it, or stops with NCBI's error (for example `HSP Filtering options invalid` for a `BestHitOverhang` or `BestHitScoreEdge` outside 0 to 0.5, or with `Culling`). After an import the command line's rules do not apply: a strategy that sets a matrix or a word size without the threshold, the window or the gap costs keeps the task's values, as in NCBI. Refused: a genetic code that has no table (NCBI crashes), composition-based statistics outside 0 to 3 (NCBI crashes on 4; -1 turns some steps of NCBI's search on and others off), `UnifiedP` above 0 in `tblastn`, or without `CompositionBasedStats`, a `WordSize` of 1 in `blastp`, `blastx` and `tblastn` (NCBI builds a lookup table of one letter), a pair of gapped extension and traceback algorithms other than dynamic programming and megablast's greedy pair (NCBI crashes on some pairs and runs its jumper or Smith-Waterman extension on others), a `StrandOption` that contradicts `-strand`, and the values that the command line refuses too (an infinite `EvalueThreshold` in `blastp`, `blastx` and `tblastn`, a NaN one in `blastx`, a `MatchReward` below 1 with a penalty other than 0). |
| `-import_search_strategy` in `psiblast`, `rpsblast` and `rpstblastn` | Refused. |
| `-export_search_strategy -` (standard output) | Refused where NCBI writes its request into the report: at NCBI's first request (after its prolog and the first batch of queries is read, before the search). The report until there is NCBI's (its prolog, written before the refusal), and a run that NCBI ends before that point (an empty query, a database error, an error of an option or of `-outfmt`, a first query record it cannot read, a `-query_loc` past the end, `BATCH_SIZE=abc`, a PSSM run with `-subject`, an `-html` report whose `./.ncbirc` is malformed, a standard output that fails at the prolog, such as `/dev/full`) ends as in NCBI, with its error and exit status. Any file or device path works, including `/dev/stdout`, the same file as `-out`, symlinks, and failing devices such as `/dev/full` or a file at its size limit (matching NCBI's diagnostics, error archives, and exit status). NCBI's warnings from before its formatter come before the refusal. |
| `-export_search_strategy` with `-dbsize` and a list that keeps no sequence of the database | Refused. |
| `-export_search_strategy` in the library API and pylosat | Refused. The command line and LOSAT Web write the file. |
| `-export_search_strategy` in LOSAT Web with the standard output (`/dev/stdout`, `/dev/fd/1`), the standard error stream (`/dev/stderr`, `/dev/fd/2`), a device or a file of the system (`/dev/full`, `/proc/self/status`), the file of `-query`, `-subject`, `-in_pssm`, `-in_msa`, `-out_pssm`, `-out_ascii_pssm`, `-logfile` or `-window_masker_db`, a file that the search reads or checks in the run's folder, also where it is missing (a file of the database of `-db` or of BLASTN's `-filtering_db`, the list of `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`, `-seqidlist`, `-negative_seqidlist`, `-taxidlist` or `-negative_taxidlist`), or any other file of the run's files (such as `.ncbirc` or the file of `-conffile`) | Refused where NCBI writes its first request, as `-export_search_strategy -` is. The command line writes the requests into that stream, device or file. LOSAT Web gives the file of any other name among the run's files, under the command line's name of the file (`s.asn` for `./s.asn`), with the command line's bytes. A name in a folder that the run's files do not have (`out/s.asn` without a folder `out`) or that names a folder of them is NCBI's error, as on the command line run in that folder. |
| `-parse_deflines` with two subjects of one id | Refused. NCBI cannot report it correctly. |
| `-parse_deflines` with two queries of one id in one query batch, or in any batch with `-subject` or `-outfmt 5` | Refused, except in two cases without `-export_search_strategy`: (1) where every searched query of a batch fails query setup (`CSeqVector_CI: cannot locate segment` because its searched interval passes the earlier query's end, or `Sequence contains no data` for an empty interval, alongside any translated query below one codon), `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `rpsblast` and `rpstblastn` end with NCBI's error and output in every output format (`psiblast` stays refused because NCBI BLAST+ crashes there), and so does a `blastn` batch with DUST (the default) or WindowMasker filtering in which any searched query passes its earlier query's end (NCBI's masking stops the whole batch there); and (2) where NCBI stops first because its BLAST XML, XML2 or JSON report (`-outfmt 5`, `13`, `14`, `15`, `16`) does not know the `-matrix` name (such as `blosum62` in lower case), then a repeated query in any batch that NCBI reads up to that error ends with NCBI's error and output, which show nothing of the repeated query, including when a non-split batch's queries have no valid context (for example only `X` or `N` letters) and the subjects pass setup, and for a repeated query that `-query_loc` skips. It stays refused with `-export_search_strategy`; where a batch's subjects fail setup and NCBI's result would depend on the earlier query's letters (or the batch is split); where a batch mixes queries that pass query setup with a repeat past its earlier query's end (NCBI warns and searches the others) and is not such a `blastn` batch with DUST or WindowMasker filtering; and where a query's ids match those of an earlier query only in part, or those of two. `psiblast` refuses XML, XML2 and JSON output formats on its own (see psiblast below). |
| `-parse_deflines` with a query and a subject of one id that are not the same record | Refused. |
| `-db` with a sequence whose id NCBI cannot report correctly | Refused. The error names the sequence. |
| `-db` with a reported sequence whose stored title holds a byte outside `' '` to `'~'` | Refused in these cases: a format other than `-outfmt` 8, 9, 11 and 12, or 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647`; a search restricted by `-seqidlist`, `-negative_seqidlist`, `-gilist`, `-negative_gilist`, `-taxids`, `-negative_taxids`, `-taxidlist` or `-negative_taxidlist`, or a database whose alias file gives a list to the sequence's volume or whose volume has an OID mask type; a query batch whose headers also give NCBI's warning of a UTF8String. NCBI writes an error message for each such byte each time it reads the title, and these cases read it in other places. Otherwise LOSAT writes NCBI's messages and output, also for a database with parsed ids and for a query batch that ends at a header NCBI cannot read, where the messages of the titles read before that header come before NCBI's error in any format. |
| `-db` with a reported sequence with parsed ids whose stored Seq-id holds a byte outside `' '` to `'~'` | Refused. NCBI shows it as `#` and writes an error message each time it reads the definition line. |
| An alias file with `MEMB_BIT` above 32 on a sequence with membership bits | Refused. See Alias files below. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI` or `SEQ_ID_AVOID_GI` with a database that has parsed ids | Refused. |
| `SERIAL_VERIFY_DATA_WRITE` set to `ALWAYS` or `DEFVALUE_ALWAYS` with `-outfmt 5` | Refused, before any output. NCBI stops with an error in its XML report or adds an empty element to it. `psiblast`, `rpsblast` and `rpstblastn` refuse it too; the other output formats are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with a `-db` whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the search, after the same output. NCBI reports the member as missing and goes on. `psiblast`, `rpsblast`, `rpstblastn`, `blastdbcmd`, `blastdbcheck`, `makeblastdb -input_type blastdb`, `dustmasker`, `segmasker`, `windowmasker` and `makembindex` refuse it too. Other runs with these values are answered as NCBI does. |
| `-xmlhelp`, `-version`, `-version-full`, `-version-full-xml`, `-version-full-json` in the library API and pylosat | Refused. They work on the command line and in LOSAT Web, which prints the command line's text. |
| `-dryrun` | Reads and checks the arguments, then stops. It refuses what a normal run refuses. |
| `-logfile` with `tblastx`, `rpsblast` or `rpstblastn` in the WASI command-line build | Refused. That build cannot send its standard error to a file. |
| `-logfile` naming an input file read by the run (such as `-query`, `-subject`, `.ncbirc` / `ncbi.ini`, `<program>.ini`, `-conffile`, `-in`, `-ustat`, `-in_pssm`, `-in_msa`, `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`, `-seqidlist`, `-negative_seqidlist`, `-taxidlist`, `-negative_taxidlist`, `-window_masker_db`, `-window_masker_taxid` counts, `-import_search_strategy`, `-mask_data`, `-taxid_map`, `-entry_batch`, `-ipg_batch`, `-ids`, `-exclude_ids`, `-i`, `rpsbproc` data files, or a BLAST database alias, alias list, GI-mask, MegaBLAST index, or component file of `-db` / `-filtering_db`; including relative paths, symlinks, hard links, and missing inputs that `-logfile` would create) | Refused (exit status 1) in all programs, before opening or creating the file, leaving existing inputs unchanged. NCBI BLAST+ 2.17.0 appends diagnostics into the input file, which can corrupt it or enter an endless read loop. |
| `-logfile` naming another existing file of the run's directory (not read by the run) or a missing file the run checks for reading (database lookup candidates of `-db`, BLASTN's `-filtering_db` or `windowmasker -infmt blastdb` — alias `.nal`/`.pal`, index `.nin`/`.pin`, group alias `index.alx`, alias `DBLIST` targets and local-copy/search-path candidates — `taxdb.bti` or paired `taxdb.btd` when `-outfmt` consults the taxonomy database, `taxonomy4blast.sqlite3` when `-taxids`/`-taxidlist` expands against a file on the search path, entries read from `windowmasker -fa_list true`, or `-archive` of `blast_formatter`), or a `windowmasker` or `blast_formatter` `-out` naming an existing file of the run's directory or a missing input the command line reads after opening `-out` (`-in` of `-mk_counts`, including `-fa_list true` entries and `-infmt blastdb` lookup files, and `-ustat` of masking) in LOSAT Web | Refused (exit status 1). On the command line those files are appended to, created and read back, or overwritten; LOSAT Web reads the run's files and writes none of them. |
| `-remote`, `-entrez_query` | Refused. They need NCBI's servers. See the end of this page. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#all-five-search-programs).

## blastn

| Option or value | What happens |
| :--- | :--- |
| `-reward 0` with a penalty other than 0, and a negative reward (a search strategy's `MatchReward`; `-reward` takes 0 or more) | Refused (exit status 1): no query has Karlin-Altschul values (NCBI writes a warning for each query and no hits). The message is shown below the table. `-reward 0 -penalty 0` (the scores of `-task rmblastn`) runs. |
| `-reward 0 -penalty 0` with `-task megablast` and its greedy extension (without `-no_greedy`; also with `-ungapped` when the gap costs are 0 and 0) | Refused. NCBI divides by zero (exit status 136) or does not finish. |
| `-use_index true` with a new-style index volume of 36 to 39 bytes or one without read access, an empty first index volume, or an `-index_name` volume range whose loop NCBI never ends (a first volume of 0 or more with a negative last one, such as `oidx,,0,-1`, or a last one of -1) | Refused. NCBI reads such a header undefined, crashes on the empty volume, or never finishes the loop. |
| `-use_index true` with an index volume whose parts do not fit together (a list, offset or local id outside its sections, an Nmer length of 0, a stride of 32,768 or more), an index not made from the database (a sequence with more chunks than the index has), a seed of the index that starts at or past the last letter of its sequence, or an old-style index (`-index_name`) of fewer sequences than the database when the search runs more than one query batch or query chunk | Refused. NCBI reads past such an index's sections or past the sequence, crashes, or gives a result that depends on the order of its threads. With the old-style index NCBI stops with an error in its second search of the index (a later query batch or query chunk); LOSAT writes the same output of the batches before and refuses there. Other searches with an index give NCBI's output, two-hit searches (`-window_size` above 0) included. When NCBI finds no index, cannot read the header of the first index volume, or the word size is below the index's minimum, LOSAT writes NCBI's error lines and warning and searches without the index, as NCBI does. With an index of only some volumes of `-db`, LOSAT searches every query batch with the index (see [Differences from NCBI BLAST+](differences.md#blastn--use_index-true-with-an-index-of-only-some-volumes)); with `-db_hard_mask`, it reads the index as if the masked letters matched nothing (see [Differences from NCBI BLAST+](differences.md#blastn--use_index-true-with--db_hard_mask)). |
| `-xdrop_gap` or `-xdrop_gap_final` so large that NCBI's 32-bit arithmetic overflows | Refused. See [Rare and pathological inputs](refusals-edge-cases.md#blastn) for the limits. |
| A very large `-xdrop_ungap` (such as `1e10`) with `-subject` that leads to an alignment outside its query context whose query coordinates leave the query on the left | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). Where they pass the query's end, or when searching a `-db` database, LOSAT matches NCBI BLAST+ in every output format. |
| `-db_hard_mask` on a database whose sequences are masked to their end, with `-task blastn-short` | Refused. NCBI can crash. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order | Refused. |
| `-db` with a nucleotide index whose ambiguity offset passes the next sequence offset (`amb > next`), at `-num_threads` above 1 when traceback reads that sequence | Refused. NCBI `blastn` aborts (`SIGABRT`) with multiple threads. At `-num_threads 1` LOSAT reports NCBI's database error (`File error: could not get ambiguity data.`, exit status 2). |
| `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused. |
| `-window_masker_taxid` with a counts file whose unit size is 17 to 251 | Refused. NCBI's result is undefined. |
| Alignments merged across a subject or query chunk overlap whose density-scaled score leaves 32-bit integers (`Int4`, e.g. a very long high-scoring match over ~1.07 Gbp with reward 2) | Refused when the search reaches it, naming the sequence by OID (exit status 1). NCBI's float-to-int conversion is undefined there. |
| `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` (see [Environment](environment.md)) that split a query batch leaving a piece without a query (NCBI BLAST+ stops in its query splitter before masking or searching any piece), or a negative chunk size above a negative overlap whose pieces leave gaps between them (NCBI BLAST+ searches on and loses the hits there) | Refused after the output of the query batches before (and for the first batch, the prologs of `-outfmt` 0, 15, 16 and 20): for an empty piece before chunk masking, database setup errors and empty-database completion; for gaps where the batch's first piece with searchable query letters would read subjects in hand, after that piece's setup (with no subject to search, or no such piece, the batch is reported as NCBI BLAST+ reports it). Where every piece holds a query and NCBI BLAST+ would split a piece again, `blastn` searches each piece once (see [Differences from NCBI BLAST+](differences.md#query-chunks-that-ncbi-would-split-again)). |

The message for `-reward 0` with the default penalty of `-task megablast` is:

```text
Error: a reward of 0 with a penalty of -2 (NCBI BLAST+'s 16-bit values), for which no query has Karlin-Altschul values, is not supported by LOSAT's BLASTN; reward 0 needs penalty 0 (matrix-only scoring, as -task rmblastn)
```

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastn).

## blastp

| Option or value | What happens |
| :--- | :--- |
| `-use_sw_tback` with `-ungapped` (which only `-comp_based_stats 0` allows) | Refused at the first query batch that is searched, as in `blastx`. NCBI crashes there, whether the batch has hits or not. |
| `-ungapped` with a `-matrix` name that begins with `IDENTITY` but is not `IDENTITY` | Refused. NCBI aborts. |
| `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` (see [Environment](environment.md)) that split a query batch where NCBI BLAST+ stops (a piece it would split again, or one left without a query), or a negative chunk size above a negative overlap whose pieces leave gaps between them (NCBI BLAST+ searches on and loses the hits there) | Refused after the output of the query batches before: at the stop, before the database's errors and `-searchsp`'s average subject length error that NCBI BLAST+ meets later, and for the gaps where the batch's first piece with searchable query letters would read the subjects. With no subject to search, or no such piece, the batch is reported as NCBI BLAST+ reports it. `psiblast` and `tblastn` do the same. |
| `-ungapped` with a `-matrix` name that NCBI reads from a matrix file | Refused. NCBI aborts or fails. See [Environment](environment.md). |
| `-ungapped` with a protein subject of more than 5,000,000 residues, from `-subject` or from a database (without `-db_hard_mask`) | Refused. NCBI searches it in overlapping parts, which LOSAT does not do yet outside `-db_hard_mask`. With `-db_hard_mask`, LOSAT splits and merges overlapping parts. |
| `-comp_based_stats` with unified P-values (`1u`, `2u`, `3u`), when a query with a residue that can be scored is followed in the same batch of queries by one with none (all `X`, empty, or masked whole by `-seg yes`) and `-searchsp` is not given or is 0 | Refused. NCBI reports every e-value of the first of the two as 0. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7; 5 is the default of `-task blastp-fast`) | Refused. NCBI crashes. |
| `-threshold` with the compressed lookup table, when the neighboring words of a query batch fill NCBI's 1024 overflow banks | Refused. NCBI corrupts its heap. LOSAT writes the results of the query batches before and then stops with the error. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |
| An infinite `-evalue` (`+inf`, `1e999`) | Refused. |
| An `-evalue` of 1.7976931348623157e308 (`DBL_MAX`) with composition-based statistics (`-comp_based_stats` 1 to 3) | Refused. NCBI's `blastp` crashes on some inputs with it. With `-comp_based_stats 0` or `-matrix IDENTITY` it runs. |
| `-db_hard_mask` on a database with a part between overlapping masks that is shorter than the word, where NCBI's gapped search of the seed it saves there never ends | Refused when the search reaches it (exit status 1). Other such parts, ranges longer than 5,000,000 letters, and overlapping parts from out-of-order masks (except on the final chunk where NCBI reads the merge overlap past its list of ranges) are searched and merged across their overlap as NCBI searches them. |
| `-db_hard_mask` on a database with masks out of order whose final part starts before the end of an earlier part from another range (or where an earlier range's alignment from a short part meets the final part across that overlap) | Refused (before the search, or when the search reaches a short part's overlap; exit status 1). NCBI reads the merge overlap past its list of ranges there, so its result depends on prior sequences and threads. Merges within a final range longer than 5,000,000 letters are searched. |
| `-db_hard_mask` with `-ungapped -window_size 0`, when two one-hit alignments merged across a chunk overlap have query lengths summing to 0 | Refused when the search reaches it (exit status 1). NCBI divides by zero when computing their merged score. |
| `-db_hard_mask` on a database with a mask that runs to the end of its sequence after the part from its start, with masks stored after that mask | Refused. NCBI does not apply the later masks when it searches that part, but applies them to the letters of its alignments. |
| `-db_hard_mask` on a database with a mask that ends before it starts, or that ends past its sequence (`makeblastdb` never writes either) | Refused. |
| `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastp).

## blastx

| Option or value | What happens |
| :--- | :--- |
| `-query_gencode 32` | NCBI's argument error (exit status 1). Code 32 is not one of NCBI's codes. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7; 5 is the default of `-task blastx-fast`), or a name that begins with `IDENTITY` | Refused. NCBI aborts. |
| An `-evalue` that is not finite (`+inf`, `1e999`, `-nan`) | Refused. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. LOSAT writes the results of the query batches before and then stops with the error (exit status 1). |
| `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` (see [Environment](environment.md)) that split a query batch where NCBI BLAST+ stops (a piece it would split again, or one left without a query: every negative chunk size above a negative overlap leaves one), or an overlap that is not a multiple of 3 (NCBI BLAST+ searches pieces out of the query's reading frames) | Refused after the output of the query batches before, the stops before the database's errors that NCBI BLAST+ meets later. The overlap that is not a multiple of 3 is refused only where the batch has subjects to search whose sequences can be read; otherwise the batch is reported, or fails, as in NCBI BLAST+. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| `-query_loc` or `-subject_loc` with `-outfmt` 5, 8, 9 and 11 to 16 | Refused. Formats 0, 6, 7, 10 and 20 work. |
| A protein subject of more than 5,000,000 residues, from `-subject` or from a database (without `-db_hard_mask`) | Refused. NCBI searches it in overlapping parts, which LOSAT does not do yet outside `-db_hard_mask`. With `-db_hard_mask`, LOSAT splits and merges overlapping parts. |
| `-db_hard_mask` on a database with a part shorter than the word whose seed NCBI's gapped search never ends on, with masks out of order whose final part overlaps an earlier range across NCBI's out-of-bounds merge overlap read, with `-ungapped -window_size 0` alignments whose merged query lengths sum to 0, with a mask that runs to the end with masks stored after it, or with a mask that ends before it starts or past its sequence | Refused, as in `blastp`. Other databases whose masks overlap, are stored out of order, or leave more than 5,000,000 letters between masks are searched as NCBI searches them; one refused sequence refuses the whole database. |
| `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. |
| `-use_sw_tback` with `-ungapped` (which only `-comp_based_stats 0` allows) | Refused at the first query batch that is searched. NCBI crashes there, whether the batch has hits or not. |
| `-use_sw_tback` on a query batch that NCBI never finishes or crashes on | Refused, as in `blastp` and `tblastn`. See [Rare and pathological inputs](refusals-edge-cases.md#blastx). |

The error names the sequence and exits with status 1: `Error: subject sequence N has M residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX` for the N-th record of `-subject` (counted from 1), and the same sentence starting with `database sequence N` for a database sequence (N is its position in the database, counted from 0).

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastx).

## tblastn

| Option or value | What happens |
| :--- | :--- |
| `-use_sw_tback` with `-ungapped` (which only `-comp_based_stats 0` allows) | Refused at the first query batch that is searched, as in `blastx`. NCBI crashes there, whether the batch has hits or not. |
| `-ungapped` with a `-matrix` name that begins with `IDENTITY` but is not `IDENTITY`, or with a matrix file that NCBI BLAST+ finds (see [Environment](environment.md)) | Refused. NCBI aborts. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7) | Refused, as in `blastp`. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. LOSAT writes the results of the query batches before and then stops with the error. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| An infinite `-evalue` (`+inf`, `1e999`) | Refused. |
| An `-evalue` of 1.7976931348623157e308 (`DBL_MAX`) with composition-based statistics (`-comp_based_stats` 1 to 3) | Refused. NCBI's `tblastn` crashes on some inputs with it. With `-comp_based_stats 0`, `-matrix IDENTITY` or `-ungapped` it runs. |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback, or with `-ungapped` the ungapped re-evaluation | Refused. NCBI crashes. |
| `-use_sw_tback` on a query and subject that NCBI never finishes, or on a query batch that NCBI crashes on | Refused, as in `blastp`. See [Rare and pathological inputs](refusals-edge-cases.md#tblastn). |
| `-db_hard_mask` on a database whose masks are stored out of the order of their start (masks that overlap or nest in that order are searched), and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused (exit status 1). A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well, and so are the few sequences with overlapping masks listed in [Rare and pathological inputs](refusals-edge-cases.md#tblastn). |
| `-db_soft_mask` on a database sequence of more than about 15,000,000 bases whose masks are stored out of order, when a gapped search finds an alignment that starts in the letters before one of the 5,000,000-letter parts that NCBI BLAST+ searches a reading frame in (mostly with `-window_size 0`), and that alignment is longer than 11 letters or its middle lies before the part | Refused, naming the database sequence and the reading frame (exit status 1): `Error: TBLASTN: database sequence N, frame F: an initial alignment of the chunk from residue R starts before the chunk (...)`. NCBI BLAST+ crashes on these. Shorter alignments whose middle lies in the part are searched, with NCBI's rows. With `-ungapped` the search runs. |
| `-db` with a nucleotide index whose ambiguity offset is out of order (`amb <= start`, where NCBI reads the byte before the sequence, or `amb > next`, where NCBI reads into the next sequence's bytes) | Refused: not yet supported (`blastn`, `blastdbcmd`, and `blastdbcheck` handle in-mapping out-of-order ambiguity offsets). |
| `-ungapped -window_size 0` on a subject searched in overlapping parts (longer than 5,000,000 residues in a frame, or `-db_hard_mask`), when two one-hit alignments merged across a chunk overlap have query lengths summing to 0 or an out-of-range scaled score | Refused when the search reaches it, naming the sequence by OID and reading frame (exit status 1). NCBI divides by zero or converts an out-of-range float to a 32-bit integer there. |
| `-in_pssm` | See the next section. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastn).

## tblastx

| Option or value | What happens |
| :--- | :--- |
| `-query_gencode 32`, `-db_gencode 32` | NCBI's argument error (exit status 1). Code 32 is not one of NCBI's codes. |
| An `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62) when an extension passes into the next frame of the query | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the ungapped re-evaluation | Refused. NCBI crashes. |
| `-db_hard_mask` on a database whose masks are stored out of the order of their start (masks that overlap or nest in that order are searched), and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused (exit status 1). A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well, and so are the few sequences with overlapping masks listed in [Rare and pathological inputs](refusals-edge-cases.md#tblastx). |
| `-db` with a nucleotide index whose ambiguity offset is out of order (`amb <= start`, where NCBI reads the byte before the sequence, or `amb > next`, where NCBI reads into the next sequence's bytes) | Refused: not yet supported (`blastn`, `blastdbcmd`, and `blastdbcheck` handle in-mapping out-of-order ambiguity offsets). |
| `-outfmt 1` to `4` with `-parse_deflines` | See "All five search programs" above. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastx).

## Alias files and database arithmetic

These refusals apply to every search program that reads databases.

| Option or value | What happens |
| :--- | :--- |
| An alias file with `MEMB_BIT` above 32 when reading a definition line that carries membership bits | Refused where the line is read. NCBI shifts a 32-bit integer by `MEMB_BIT - 1`, which is undefined behaviour in C++. `MEMB_BIT` 1 to 32 (0 disables filtering), `GILIST`, `TILIST`, `SEQIDLIST` (version 4 and 5), `TAXIDLIST` (version 5; ignored on version 4 as in NCBI) and `MASKLIST` are read (see [Databases](databases.md)). Lines without membership bits are hidden deterministically, and commands that do not read filtered headers (`blastdbcmd -info`, `blastdbcheck`) run as in NCBI. |
| A `SEQIDLIST` file with a line that NCBI warns about, with an id that NCBI's accession lookup cannot simplify, or whose comparison with a search list posts a Seq-id warning | Refused. |
| A version 5 alias `SEQIDLIST` combined with `-seqidlist` or `-negative_seqidlist`, when the alias `SEQIDLIST` or the search's `-seqidlist` holds two differently spelled Seq-ids that NCBI treats as equivalent (such as `NP_000103` and `NP_000103.1`, or `LOC7` and `loc7`), or two unparseable Seq-ids that raise different exceptions, or an unparseable Seq-id alongside a Seq-id that posts a warning | Refused. NCBI's comparator for sorting these lists is not a strict weak ordering and parses ids during the sort, so its result or which error it reports depends on its C++ library's sort steps. Without a search seqid list, with a single unparseable exception text, or with equivalent pairs only inside `-negative_seqidlist`, the search runs as in NCBI. |
| An alias file whose `OID_MASK_TYPE` hides predicted RefSeq lines, over a sequence with an id that NCBI warns about when it reads the id's accession | Refused where the sequence's lines are read. NCBI writes the warning again at each read. |
| An `OIDLIST` file that is empty or shorter than 5 bytes (no byte after its 4-byte header) | Refused. NCBI crashes on an empty file, and its error for the other names a line of its source. A file whose bits NCBI would read past the memory page after its end is refused too. |
| A GI, TI or PIG index of a database (`.pni`/`.pnd`, `.pti`/`.ptd`, `.ppi`/`.ppd` and the `.n` files) whose data file is empty, or ends more than a memory page before the last term its index counts, where `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`, an alias `GILIST` or `TILIST`, or a `blastdbcmd` lookup opens it | Refused. NCBI reads the first and last terms of the data file when it opens the index, whatever the file's length, and there it reads past its mapping of the file. A data file of another length that ends inside that page is read as NCBI reads it: the index is not used, a GI or IPG list ends in NCBI's batch-mode error and a `TILIST` finds no sequence. |
| A volume's `.nog` (the GI of each sequence) that is empty, or that ends more than a memory page before the entry of a sequence, where the GI-indexed masks of that sequence are read (`blastdbcmd -outfmt %m` and `-mask_sequence_with`, `-db_soft_mask` and `-db_hard_mask`, the copy of `makeblastdb -input_type blastdb`, `makembindex -db_mask`) | Refused. NCBI crashes on an empty file, and past the page it reads whatever memory follows. A `.nog` that ends inside that page is read as NCBI reads it: the bytes past its end are zeros, and a sequence past its count takes the GI of its definition lines. |
| An alias file whose `FIRST_OID`-`LAST_OID` range lies outside the OIDs that the ranges of the alias files below it keep (it ends before them, or starts past the bytes NCBI keeps for them) | Refused where the database's OID list is first used. NCBI clears bits outside its array of those OIDs there: before the array it crashes, and past it it writes over its own memory, which crashes it unless the bytes happen to fall in unused memory. `blastdbcmd -info` does not use the list and works. |
| Databases of a `-db` list whose total lengths sum past 64 bits (`-outfmt` 0 to 5 and 13 to 16) | Refused. NCBI's arithmetic overflows. |
| Databases of a `-db` list whose numbers of sequences sum past 32 bits (`-outfmt` 5) | Refused. |
| `rpsblast` and `rpstblastn`: a reported hit of a profile that an alias file's `FIRST_OID`, `LAST_OID` and `OIDLIST` leave out | Refused. NCBI's search ignores these keys, and its report does not find the profile: it writes `Unknown` or stops with a null pointer exception. A profile whose id an included copy of its record holds is written as that copy, as NCBI finds it. |
| `rpsblast` and `rpstblastn` over an RPS database of several volumes: a reported hit of a profile whose definition line an alias file's `OID_MASK_TYPE`, id list (`GILIST`, `TILIST`, `SEQIDLIST` or `TAXIDLIST`) or `MEMB_BIT` hides | Refused. NCBI searches each volume without the alias files, and its report does not find the profile: it writes `Unknown` or stops with a null pointer exception. With one volume NCBI drops such a profile's hits, as LOSAT does. |

`TITLE`, `DBLIST`, `NSEQ`, `LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`, `FIRST_OID`, `LAST_OID` and `OIDLIST` work as in NCBI BLAST+, and `VDBLIST` and `MAX_SEQ_LENGTH` are ignored, as NCBI's reader ignores them. In the profile databases of `rpsblast` and `rpstblastn` they act as in NCBI's search of such a database: one volume is searched with all its keys, and each of several volumes with its own length, count and shortest length under the database's `LENGTH` or `STATS_TOTLEN` and `NSEQ` or `STATS_NSEQ` where these are above 0.

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#alias-files-and-database-arithmetic).

## psiblast

| Option or value | What happens |
| :--- | :--- |
| `-phi_pattern` with `-num_iterations` other than 1 | Refused: not yet supported. With a matrix and gap costs that have no pattern statistics, NCBI searches nothing and LOSAT gives its answer of one round. |
| `-phi_pattern` with `-outfmt` 1 to 5, 11, 13 to 16 or 18 | Refused: not yet supported. `-outfmt` 8, 9 and 12 get NCBI's own error (exit status 3) where an alignment is found, and no output otherwise. |
| `-phi_pattern` with `-html` | Refused: not yet supported. |
| `-phi_pattern` with `-sum_stats true` | Refused: not yet supported. |
| `-phi_pattern` with a subject longer than 5,000,000 residues | Refused: not yet supported. NCBI searches such a sequence in overlapping parts. |
| `-phi_pattern`: the order of these refusals | NCBI's own errors come first: a pattern that NCBI rejects, a query without the pattern, a last query record without letters, and with `-outfmt 11` NCBI's error archive. With a matrix and gap costs that have no pattern statistics, LOSAT gives NCBI's answer (its message and no hits) and refuses only `-outfmt` 5 and 11 and `-html`. |
| `phi_pattern` in the library and pylosat | Refused: not yet supported. Use `losat psiblast -phi_pattern`. |
| A `-phi_pattern` pattern with an unterminated `(`, `{` or `x(` | Refused. NCBI reads past the end of the pattern and crashes. |
| A `-phi_pattern` pattern with `x(n,m)` that does not give two numbers (`x(a,b)`, `x(3,)`) | Refused. NCBI reads values that were never set, and its result changes from run to run. |
| A `-phi_pattern` pattern with a count outside 32-bit integers | Refused. NCBI's count overflows. |
| A `-phi_pattern` pattern with positions past NCBI's 330 (`A(400)`, `M-x(330)`) | Refused. NCBI writes past its pattern array before it reports the pattern as too long. |
| A `-phi_pattern` pattern with a byte that is not ASCII outside `[]` or inside `{}` | Refused. NCBI reads past its residue table. |
| A `-phi_pattern` pattern that needs more than two machine words, with two variable wildcards in a row or more than 100 words | Refused. NCBI shifts by a negative count or writes past its arrays. |
| A `-phi_pattern` search where a sequence has more than 10,000 hits of a pattern of more than one machine word | Refused. NCBI writes past its hit array and stops (`stack smashing detected`). |
| A `-phi_pattern` alignment with more than 32 hits of one word of a pattern that needs more than two machine words | Refused. NCBI writes past its array of 64 entries for the hits of one word. |
| A `-phi_pattern` alignment with more than 20,000 word positions of a pattern that needs more than two machine words | Refused. NCBI writes past its array of 20,000 entries for the word positions. |
| `-in_msa` with a PHYLIP header of 0 sequences, or an empty command (a lone `;`) in a NEXUS file | Refused. NCBI crashes (a division by zero, a null pointer). |
| `-in_msa` with a FASTA alignment whose rows have different lengths, where a row that holds the query is shorter than the first row, or the first row has a letter past the end of a shorter row | Refused. NCBI swaps or reads the shorter row past its end, and its result depends on the memory allocator. Make the rows the same length. |
| `-in_msa` with a `secondary-accession` range that runs backwards (`AB000003-AB000001`) or has a number of more than 18 digits | Refused. NCBI loops without end on a backwards range and reads the numbers into a 32-bit integer that overflows; LOSAT reads 64 bits, so a number past 2147483647 is written with its digits, where NCBI's build wraps it. A range of any length is answered. |
| `-in_msa` with a master sequence from which no ungapped Karlin-Altschul parameters can be computed (X, O, `*` or lower-case x alone, or X with gaps: no residue that the parameters can be computed from), or with `-matrix IDENTITY` | Refused. NCBI crashes. |
| `-in_msa` with a master sequence that is the start of a longer sequence above it in the alignment (`-msa_master_idx` names it) | Refused. NCBI reads and writes past its alignment buffers there, and its result depends on the memory allocator. Put the master sequence above the longer one in the file. |
| `-use_sw_tback` | Refused. |
| `-outfmt` 13 to 16 | Refused: not yet supported, with a FASTA query, `-in_pssm` or `-in_msa`. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop that overflows NCBI's kernels | Refused, as in `blastp`. This includes later rounds, and `-in_pssm` with the PSSM's Lambda. |
| An infinite `-evalue`, or an `-evalue` of 1.7976931348623157e308 with composition-based statistics | Refused, as in `blastp`. |
| A `-subject` set whose sequences all have no letters | Refused. NCBI's `psiblast` crashes; `blastp` runs it. |
| `-outfmt` 1 to 4 with an `-in_pssm` PSSM whose query has the molecule `other`, at a round whose merge reads past a row's end, crashes, reads freed memory or wraps a segment length in NCBI | Refused after the report before the round's query-anchored section, which is NCBI's. NCBI reads such a query as a nucleotide sequence when it merges the alignments, three positions to a letter on three rows: where it reads past the last piece of a row while it orders the pieces, it stops with a C++ exception (`CObject is corrupted`, exit status 255); where a query letter against a gap of the subject falls on a row that no aligned piece opened, it crashes; where filling a row's unaligned stretch replaces a piece that it still lists, it reads memory it has freed; where a piece's walk finds the next piece before its own position, it computes a length below zero and then runs forever or goes on with lengths near 2^32. Rounds whose merge returns (an empty line for the section) or stops with NCBI's own merge exception (`CAlnException`, exit status 255) are answered as NCBI does. |
| `-outfmt 11` with a database and an `-in_pssm` file that holds frequency ratios only, without `params.rpsdbparams`, in the ASN.1 text archive (no `ARCHIVE_FORMAT`), `ARCHIVE_FORMAT=xml` or `ARCHIVE_FORMAT=bin` (NCBI stops with a C++ exception while it writes the request's PSSM, as it does for `-out_pssm`) | Refused. `-subject`, `-outfmt 6` and `-export_search_strategy` with such a file are answered as NCBI does. |
| `-outfmt 11` with an `-in_pssm` file whose author names (a Pubdesc's authors, which NCBI reads again only when it writes them in another format) are in another format than the archive (`ARCHIVE_FORMAT=xml`, `ARCHIVE_FORMAT=bin` or no variable for ASN.1 text) and hold a number too large for its type | Refused: not yet supported (NCBI stops with an exception after part of the archive). Every other PSSM file is written in all three archives, with NCBI's message for a byte outside the visible ASCII range in those names, and so are `-in_msa` and a FASTA query. A list element stated as `NULL` (`05 00` in ASN.1 binary) stops each archive at the element, as NCBI's writer stops, and the error archive follows. |
| `-out_pssm` or `-export_search_strategy` with an `-in_pssm` file in ASN.1 binary or XML whose author names (a Pubdesc's authors, which NCBI reads again when it writes them in ASN.1 text) hold a number too large for its type | Refused: not yet supported (NCBI writes part of the file and stops with an exception). With `-outfmt` 0 to 4 the refusal comes after the start of the report, which NCBI writes before that file. A list element stated as `NULL` is written as NCBI writes it: the file up to that element, then NCBI's exception (exit status 255). |
| `SERIAL_VERIFY_DATA_WRITE` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-in_pssm` file that holds a list element stated as `NULL`, where the run writes the PSSM (`-export_search_strategy` with `-db`, `-out_pssm`, `-outfmt 11` with `-db`) | Refused, before any output. NCBI leaves the element out there, where it otherwise stops with an error, as LOSAT does. `tblastn -in_pssm` refuses it too, for `-export_search_strategy` with `-db`. Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-in_pssm` file that lacks a member NCBI's reader requires | Refused, before any output, where NCBI's error would stop the run. NCBI reports the member as missing and goes on. `tblastn -in_pssm` refuses it too. Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_GET` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-in_pssm` file whose query has a feature table column without the `field-name` that NCBI reads | Refused where NCBI's error would stop the run (after the prolog of `-outfmt 0`). NCBI takes an empty name and goes on. `tblastn -in_pssm` refuses it too, before any output. Other runs with these values are answered as NCBI does. |
| `-out_pssm` or `-out_ascii_pssm` with a file that cannot be written | Refused. |
| `-out_pssm` or `-out_ascii_pssm` into the standard output, the standard error or the `-out` file | Refused. Where NCBI stops first with its own error of `-outfmt` (13 or 14 without `-out`, for example), you get that error. |
| `-out_pssm` into a named pipe, a block device or a symbolic link to a missing file | Refused. |
| The library calls `run` or `run_collect` with formats that take two searches (two hit list sizes under `-num_descriptions` or `-num_alignments`, or `psiblast`'s outfmt 5 or 11 with other formats) and PSSM files | Refused. One call per format writes the command line's output. |
| `-db` with a reported sequence whose stored title holds a byte outside `' '` to `'~'` | Refused in every format, including `-outfmt` 8, 9 and 11 (`psiblast`'s PSSM construction and unpruned hit lists read definition lines outside the reported Seq-aligns), except for a query batch that ends at a header NCBI cannot read (which writes NCBI's messages for the titles read before it, then NCBI's error, in any format). |

`-matrix IDENTITY` gives NCBI's own error.

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#psiblast).

## tblastn with `-in_pssm`

`tblastn -in_pssm` (PSI-TBLASTN) refuses the `tblastn` entries above and the PSSMs and PSSM queries that `psiblast -in_pssm` refuses. It refuses these too, where NCBI BLAST+ 2.17.0 crashes or aborts:

| Option or value | What happens |
| :--- | :--- |
| `-ungapped` | Refused. |
| `-use_sw_tback` | Refused when the search finds an alignment, after the report's header. Depending on the input, NCBI fails its traceback (`BLAST engine error: Traceback failed`), throws an exception, crashes (with `-subject` and with `-db`), never finishes, or writes rows that change from one run to the next. A search that finds nothing runs, as in NCBI. |
| `-outfmt 11` | Refused. NCBI searches, then stops with a C++ exception whose report (in its error archive too) changes from run to run. |
| `-mt_mode 1` with more than one thread and `-db` | Refused. |
| A PSSM query that `-seg` masks completely | Refused. LOSAT writes the report's header and then stops with the error. |
| `-db_gencode 32` | Refused. It is not one of NCBI's codes; LOSAT's `tblastn` takes it for its own translations. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastn-with--in_pssm).

## rpsblast

| Option or value | What happens |
| :--- | :--- |
| `-outfmt` 17, 19 and 21, and 13 and 14 without `-out` | NCBI's own errors. |
| `-entrez_query`, `-remote`, `-import_search_strategy` | Refused. |
| `-sum_stats true` | Refused. NCBI's output then depends on memory contents, or it crashes. |
| `-subject` | NCBI's argument error, as in NCBI (usage and `Unknown argument: "subject"`, exit status 1). It is not an `rpsblast` option. |
| `-outfmt 1` to `5` or `13` to `16`, `-outfmt 18` when the profile has a taxid, or `-outfmt 20` with a field that shows more of the profile than its id (`sallseqid`, `sallacc`, `stitle`, `salltitles`, `staxids`, `slen`, `sseq`, `btop`, `nident`, `pident`, `mismatch`, `positive`, `ppos`, as in the default fields, or a taxonomy field when the profile has a taxid), with a profile whose id is `lcl|Query_N`, the generated id of a query of its batch, in the hits | Refused. NCBI takes the query for the entry (the query's title in `Hit_def`, its letters in the identities, its taxonomy in the organism report). `-outfmt 8`, `9`, `11` and `12`, and `-outfmt 18` and `20` otherwise (such as `20 qseqid sseqid sacc evalue`), write NCBI's bytes (outfmt 11 apart from its `ident`: `LOSAT 0.2.0` where NCBI writes `2.17.0+`). |
| A query above 5,000,000 letters | Refused. |
| A database that lists one volume twice (`-db "x x"`) | Refused. The database layer refuses it for every program; NCBI's `rpsblast` searches the volume once. |
| A profile database whose volumes differ in matrix or gap costs | Refused. |
| An empty or inconsistent `.rps` or `.loo` file, or an inconsistent `.aux` file | Refused. NCBI reads past its buffers or crashes. |
| A hit on a profile whose consensus letters in `.psq` are more or fewer than its score columns in `.rps`, where NCBI's traceback then reads past the profile's rows, composition tables, copy of the consensus or score matrix, or starts an alignment past the consensus | Refused, for the query batch of that hit. NCBI crashes or reads memory it does not own. Other hits on such a profile are searched as NCBI searches them. `makeprofiledb` writes such a profile from a scored profile whose sequence is nucleotide-coded or of another length. |
| `rpsblast` in the library API and the Web engine | Refused. |
| The profile-database cases named in Alias files above | Refused. |
| An infinite `-evalue`, or an `-evalue` of 1.7976931348623157e308 with composition-based statistics (the default of `rpsblast`) | Refused, as in `blastp`. |
| `-db` with a reported profile whose stored title holds a byte outside `' '` to `'~'` | Refused in every format, including `-outfmt` 8, 9, 11 and 12 (NCBI's `CLocalRPSBlast` opens a fresh `CSeqDB` per query batch and per volume and merges volume alignments before `PruneSeqalign`), except for a query batch that ends at a header NCBI cannot read (which writes NCBI's messages for the titles read before it, then NCBI's error, in any format). |

NCBI's errors for a missing `.aux`, `.loo` or `.rps` file, a wrong magic number or an empty `.aux` come as in NCBI (exit status 3).

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#rpsblast).

## rpstblastn

| Option or value | What happens |
| :--- | :--- |
| `-outfmt` 17, 19 and 21, and 13 and 14 without `-out` | NCBI's own errors. |
| `-ungapped` without composition-based statistics | Refused. NCBI's `rpstblastn` crashes after the search. With composition-based statistics it gets NCBI's own error. |
| `-sum_stats true` | Refused. NCBI's output then depends on memory contents, or it crashes. |
| `-entrez_query`, `-remote`, `-import_search_strategy` | Refused. |
| `-subject` | NCBI's argument error, as in NCBI (usage and `Unknown argument: "subject"`, exit status 1). It is not an `rpstblastn` option. |
| `-outfmt 1` to `5` or `13` to `16`, `-outfmt 18` when the profile has a taxid, or `-outfmt 20` with a field of the alignment (such as `qstart`, `length` or the default fields) or a field that shows more of the profile than its id (`sallseqid`, `sallacc`, `stitle`, `salltitles`, `staxids`, `slen`, or a taxonomy field when the profile has a taxid), with a profile whose id is `lcl|Query_N`, the generated id of a query of its batch, in the hits | Refused. NCBI stops with an exception (after half a document in BLAST XML, after the header line in `-outfmt 20`) or takes the query for the entry. `-outfmt 8`, `9`, `11` and `12`, and `-outfmt 18` and `20` otherwise (such as `20 qseqid sseqid sacc evalue`), write NCBI's bytes (outfmt 11 apart from its `ident`: `LOSAT 0.2.0` where NCBI writes `2.17.0+`). |
| `-query_gencode 32` | NCBI's argument error (exit status 1). |
| `rpstblastn` in the library API and the Web engine | Refused. |
| `-use_sw_tback` with a profile whose Smith-Waterman redo NCBI repeats without end | Refused, for a window of a query frame: the e-value of a score of 0 passes `-evalue`, the window holds two or more of the profile's preliminary HSPs, and the profile's first position scores the window's first letter above 0. With `-strand minus` NCBI's redo reads a search space of 0, so every score passes `-evalue`. The refusal comes with the query batch that holds the profile, after the reports of the batches before it (in `-outfmt 0` after the header that NCBI writes before its first batch). |
| `-use_sw_tback` with a profile whose Smith-Waterman alignment starts with a gap in the query's row | Refused, at the same place. No search has been found that gives one, so there is no NCBI output to follow. An alignment that starts with a gap in the profile's row is written as NCBI writes it: the warning that the sequence data could not be pre-fetched, and `-outfmt 0` to `5` and `13` to `16` stop where NCBI's do. |
| The database refusals of `rpsblast` above, with its alias-file cases | Refused. |
| `-comp_based_stats 0` with a hit on a profile whose consensus letters in `.psq` are fewer than its score columns, where an alignment starts at or past those letters | Refused, for the query batch of that hit. NCBI prints alignments there whose positions, scores and identities disagree, or crashes when such an alignment then meets a letter of 28 or more. |
| `-db` with a reported profile whose stored title holds a byte outside `' '` to `'~'` | Refused in every format, including `-outfmt` 8, 9, 11 and 12 (NCBI's `CLocalRPSBlast` opens a fresh `CSeqDB` per query batch and per volume and merges volume alignments before `PruneSeqalign`), except for a query batch that ends at a header NCBI cannot read (which writes NCBI's messages for the titles read before it, then NCBI's error, in any format). |

NCBI's errors for a missing `.aux`, `.loo` or `.rps` file and for an empty `.aux` come as in NCBI (exit status 3). The refusals of options and databases come where `rpsblast`'s come (before the first byte of output, or with `-outfmt` 0 to 4 after the start of the report for one found in the search), also for queries without a valid frame, which NCBI does not search.

## rpsbproc

`rpsbproc` gives the domain hits, sites and motifs of protein and nucleotide queries.

| Input | What happens |
| :--- | :--- |
| A motif match whose result depends on bytes that NCBI reads beyond the residues of a site | Refused. |
| A site or motif at a position outside the query's residues (other than their end) | Refused. NCBI reads memory outside the string. |
| A `[datafiles]` entry of a registry file (`rpsbproc.ini` or the file of `-conffile`) | Refused. |
| A `.ncbirc` with a syntax error | Refused. NCBI reports it with the name of a source file of its own build. |
| The environment variables `OBJECTS_SEQ_DESCR_ALLOW_EMPTY`, and `SERIAL_VERIFY_DATA_READ` or `SERIAL_VERIFY_DATA_GET` set to turn the checks off | Refused. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI` or `SEQ_ID_AVOID_GI` with a query that has a gi | Refused. |
| An archive with parts LOSAT does not read | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#rpsbproc).

## makeblastdb

| Option or value | What happens |
| :--- | :--- |
| `-gi_mask` with two mask files whose GI-mask files are the same files (one `-gi_mask_name` twice, or two names of one path), when both would write one data file at the same time | Refused at the sequence where the second would start writing that file. NCBI's files then depend on how its C++ runtime buffers two streams on one file. The log head and the removal of an old database of the `-out` name come first; the build then ends as after an error: its volume files are removed, and the GI-mask files written so far and the alias file stay. Where one of the two writes nothing, or their data files are written one after the other, LOSAT writes NCBI's files. |
| `-input_type asn1_txt`, `asn1_bin` with titles outside the explicit, source-taxname, and default protein paths | Refused (`MolInfo` tech prefix/suffix, `Unverified`/`Unreviewed` title prefix, feature-derived protein names, `BioSource` modifiers, empty-title RefSeq `NC`/`NM`/`NR`, `WGS`/`TSA`/`TLS`, `PDB`, patent, segmented/map, and `AutoDef` title generation). Sequences with explicit titles, `ASN1_BlastDefLine` headers, `CFastaReader` user objects, source taxname titles, and default protein titles work. |
| `-input_type blastdb` with an empty `-in` and no `-title` | Refused. NCBI reads the title of a database that is not there and crashes. |
| `[BLAST] LONG_SEQID`, or a `[BLAST] MAX_SEQID_LENGTH` that `makeblastdb` would use with `-parse_seqids` | Refused. They change which ids are parsed. |
| `DIAG_*` variables and other `NCBI_CONFIG_*` entries outside the accepted list | Refused. |
| A non-empty `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` with `-input_type fasta` | Refused. It changes which titles are stored for FASTA input. With `-input_type asn1_txt` and `asn1_bin`, it is read as in NCBI BLAST+. |
| A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more: 4,294,967,293 to 4,294,967,295 bases of `A`, `C`, `G` and `T` in a row, or 4,294,967,295 bases in a row with other letters among them (the runs of `N` that NCBI makes of `>?` gap lines count) (`-dbtype nucl`) | Refused. NCBI crashes: it counts the bytes of the piece in 32 bits and gets 0. Other records of 4,294,967,296 letters or more are handled as NCBI handles them, with their length modulo 2^32. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with a `-mask_data` file, an `-input_type asn1_txt` / `asn1_bin` input, or an `-input_type blastdb` database that lacks a member NCBI's reader requires | Refused where NCBI's error would stop the build, after the same output; the build's files are removed. NCBI reports the member as missing and goes on (with a database, it then stops with an error writing the database). Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_GET` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-input_type asn1_txt` / `asn1_bin` input whose `Seq-inst.length` (or `Seq-literal.length`) is unset | Refused where NCBI's error would stop the build, after the same output; the build's files are removed. NCBI takes 0 and goes on. Other runs with these values are answered as NCBI does. |
| An index file (`.pin`, `.nin`) of a finished build that a failed write (a full disk) left with no bytes, or with a header field that reaches past the end of the file's last 4 KiB page | Refused. NCBI maps an empty file as a null pointer and crashes. What NCBI reads past its mapping depends on its memory. A file cut short inside that page is read as NCBI reads it, with zeros past its end. |

Example:

```text
losat makeblastdb: -gi_mask: two masking algorithms would write one GI-mask data file at the same time (the files of "masks"); NCBI makeblastdb's bytes then depend on how the C++ runtime buffers two streams on one file, and LOSAT does not reproduce that (a -gi_mask_name that names other files for one of them avoids it)
```

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#makeblastdb).

## blastdbcmd

| Option or value | What happens |
| :--- | :--- |
| A `-db` list in which an alias and another name give the same volume | Refused. |
| FASTA output of a sequence with a tax id, or `-outfmt` with `%S`, `%L`, `%B`, `%K`, `%C` or `%N`, or `-tax_info` with a name field, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. `-outfmt "%a %T"`, other formats that look nothing up, and `-info` work. |
| `-range ""`, `-taxids abc`, or a `-taxids` or `-taxidlist` id that is empty, blank, too large or not a number | Refused. NCBI crashes or prints a path. |
| `-mask_sequence_with` below -1 | Refused. NCBI crashes. |
| `-ipg` or `-ipg_batch` on a database without PIG files | Refused. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with a `-db` whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the run, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastdbcmd).

## blastdbcheck

`blastdbcheck` runs every option of its `-help`, including `-cdd_delta` (`-xmlhelp` on the command line only).

| Option or input | What happens |
| :--- | :--- |
| `-stride 0` or a negative `-stride` | Refused after the first OID passes. NCBI tests the first OID again and again and never finishes. |
| `-random` or the default sample of 200 OIDs on a database that includes none of its OIDs (an alias file's OID range or list leaves none) | Refused. NCBI draws OIDs without end. |
| An empty (0-byte) sequence file (`.psq`, `.nsq`) or CDD profile file (`.loo`, `.rps`, `.freq`, `.wcounts`, `.obsr`), a cut database or profile file whose offsets pass its last 4 KiB memory page, a nucleotide index whose ambiguity offset causes NCBI to read before the start of the sequence file mapping (such as offset 0), a cut string ISAM file (`.psd`, `.nsd`), or a file found by `-dir` that is not a database index | Refused. NCBI crashes or reads outside its memory mapping. Cuts that stay within the last 4 KiB page read zero-filled bytes as NCBI does, and out-of-order ambiguity offsets that stay within the mapping match NCBI's behavior. |
| `-must_have_taxids` with a sequence whose first definition line has no Seq-id (only a damaged database has one) | Refused. NCBI reads a null reference. |
| An error report at `-verbosity 3` or 4 whose text LOSAT does not know | Refused. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with a `-db` or `-dir` volume whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would fail the check of that OID, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

## dustmasker

`dustmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A title whose text encoding NCBI cannot tell | Refused. NCBI stops with a string error. |
| A database definition line with bytes NCBI reads as `#` | Refused. |
| A database sequence with a tax id, in any `-outfmt`, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. |
| A failed write of a `seqloc_*` or `maskinfo_*` output | Refused. NCBI's serial stream reports it; the other formats go on silently, as NCBI's. |
| A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more: 4,294,967,293 to 4,294,967,295 bases of `A`, `C`, `G` and `T` in a row, or 4,294,967,295 bases in a row with other letters among them (the runs of `N` that NCBI makes of `>?` gap lines count) | Refused. NCBI crashes: it counts the bytes of the piece in 32 bits and gets 0. Other records of 4,294,967,296 letters or more are read as NCBI reads them, with their length modulo 2^32. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-infmt blastdb` input whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the run, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#dustmasker).

## windowmasker

`windowmasker` runs NCBI's three modes, making counts files (`-mk_counts`, with `-checkdup`), converting them (`-convert`) and masking with them (`-ustat`), and every option of its `-help` (`-xmlhelp` on the command line only). Options of another mode are refused by NCBI's own argument check, as in NCBI's; `-meta` is accepted with `-ustat` and does nothing, as in NCBI's.

| Option or input | What happens |
| :--- | :--- |
| `-mk_counts` of an input without a unit of unambiguous letters (all `N`; or, with `-unit`, no sequence, as when `-ids` keeps none) | Refused. NCBI writes the percentages of the file as its C runtime prints 0/0 (`-nan` on Linux). |
| `-mk_counts -infmt seqids` | Refused. NCBI fetches each sequence from NCBI's servers. |
| `-mk_counts` or `-convert` of units of 1 to 3 letters to `oascii` or `obinary` with an `-smem` that is a multiple of 65536 | Refused. NCBI keeps the low 16 bits of `-smem` (0 megabytes); its search for a hash table is then undefined. With units of 4 letters or more NCBI's error is written as NCBI writes it. |
| A counts file whose `##pct:` line holds a number that C++ libraries read differently (`inf`, `nan`, a hexadecimal number such as `0x10`, a number with letters after it such as `5x`, or a number beyond the range of a double), with a `-t_*_pct` option or `-convert` | Refused. Other lines are read as NCBI reads them: a value that is missing, or a word that does not start with a number, is 0. |
| An `ascii` counts file with units of 17 letters or more | Refused. NCBI's result is undefined. |
| A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more: 4,294,967,293 to 4,294,967,295 bases of `A`, `C`, `G` and `T` in a row, or 4,294,967,295 bases in a row with other letters among them (the runs of `N` that NCBI makes of `>?` gap lines count) | Refused. NCBI crashes: it counts the bytes of the piece in 32 bits and gets 0. Other records of 4,294,967,296 letters or more are read as NCBI reads them, with their length modulo 2^32. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-infmt blastdb` input whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the run, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#windowmasker).

## segmasker

`segmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| Everything that `dustmasker` refuses: titles, database sequences and definition lines, the serial writer, and a failed write of a `seqloc_*` or `maskinfo_*` output | Refused, as for `dustmasker`. A FASTA record of 4,294,967,295 letters ends with NCBI's error. |
| A protein database whose residue bytes are above 127 | Refused. `makeblastdb` writes 0 to 27; NCBI's SEG reads larger bytes as negative characters. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-infmt blastdb` input whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the run, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#segmasker).

## convert2blastmask

`convert2blastmask` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A `Blast-db-mask-info` that LOSAT's serial writer cannot write as NCBI does | Refused. |
| A failed write of a `maskinfo_*` output | Refused. NCBI's serial stream reports it. `interval` goes on silently, as NCBI's, as for `dustmasker` and `segmasker`. |
| A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more: 4,294,967,293 to 4,294,967,295 bases of `A`, `C`, `G` and `T` in a row, or 4,294,967,295 bases in a row with other letters among them (the runs of `N` that NCBI makes of `>?` gap lines count) | Refused. NCBI crashes: it counts the bytes of the piece in 32 bits and gets 0. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#convert2blastmask).

## blast_formatter

`blast_formatter` runs every option of its `-help` (`-xmlhelp` on the command line only). It writes `-outfmt` `6` to `12` and `20` from archives in ASN.1 text, ASN.1 binary or XML, the tabular formats (`6`, `7`, `10` and `20`, with any field list) for archives of `blastn`, `blastp`, `psiblast`, `deltablast`, PHI-BLAST of a protein query (`psiblast -phi_pattern`) and `rpsblast`.

| Input | What happens |
| :--- | :--- |
| `-outfmt` `0`, `1`, `2`, `3`, `4`, `5`, `13`, `14`, `15`, `16`, `17` and `18` | Refused, with LOSAT's own error message. |
| The tabular output formats (`6`, `7`, `10` and `20`) of an archive of a translated search (`blastx`, `tblastn`, `tblastn -in_pssm`, `tblastx`, `rpstblastn`), of queries given as locations, of `-subject` sequences with `BL2SEQ_LEGACY` set in `-outfmt 7`, whose Seq-aligns name a sequence that NCBI BLAST+'s scope may find otherwise (a subject that the database holds under a header with a UTF8String), or whose Seq-aligns or Bioseqs take a form that BLAST+ does not write and whose result NCBI BLAST+ does not decide (each named in the message, such as a subject Bioseq with a BioSource while taxonomy fields are asked for) | Refused, with the case named. |
| A `-subject` sequence (`BL2SEQ_LEGACY` unset) in a form whose letters LOSAT does not read as NCBI BLAST+ does when it sets up the subject sequences: packed letters (`ncbi2na`, `ncbi4na`) that end inside their last byte before the length the sequence declares (NCBI BLAST+ reads past the data), segments of other sequences, segments that pass 4,294,967,294 letters | Refused, with the case named, in every output format. Its other errors (letters short of the declared length, a length of 4294967295, a residue outside NCBIstdaa's 28 codes, a Seq-data coding that does not convert) come as in NCBI (NCBI's engine error), and a sequence of 0 letters gets NCBI's warning. NCBI's own programs do not write such archives. |
| A database named in the archive that cannot be opened here as NCBI opens it (with the archive's taxid list) | Refused. NCBI BLAST+ then asks NCBI's servers for its sequences. |
| An archive on which NCBI BLAST+ has no defined result: a query or subject sequence without an id, a request without a query or subject, an element stated as `NULL` in a list, values that make NCBI's code abort, and in `-outfmt 11` a `FilterString` number that NCBI's `sscanf` does not read (such as `abc` in `D abc 64 1`) when no later option of the archive sets that value again | Refused. NCBI reads memory the archive does not decide, follows a null pointer or aborts; for the `FilterString` number, it writes the value of a variable it never set. NCBI's own programs do not write such archives. |
| An archive in a form that LOSAT does not read | Refused, with the form named. |
| A stored message of severity Fatal in the archive | Refused. NCBI BLAST+ posts it with a stack trace of its own binary and stops with exit status 255. |
| `PRE_FETCH_SEQS_LIMIT=2147483647` with `-outfmt` `8`, `9` or `12` | Refused. NCBI BLAST+ then loads every subject sequence ahead, from NCBI's servers where the database lacks one. |
| An archive of the mapper | Refused. NCBI's programs write none. |
| A value in an archive that LOSAT does not yet write in `-outfmt` `11` (request or results) or `12` (results): a NaN or an infinite real number in `12`, a UTF8String, a NULL value or an OCTET STRING | Refused, with the value and the format named. NCBI BLAST+ stops at a NaN or an infinity in `-outfmt` `12` with an error of its JSON stream. |
| `-outfmt 11` of an archive whose queries are a list of sequence locations (`seq-loc-list`) | Refused. NCBI BLAST+ fetches the sequences through its data loaders. |
| `-rid` | Refused. It names the results of a search that NCBI's servers ran, as `-remote` does for the search programs. |

An empty `-rid`, the archive's other errors (an archive of PHI-BLAST is NCBI's engine error in `-outfmt` `8`, `9` and `12`), the errors of `-outfmt` and a failed write of the output come as in NCBI. Without `BL2SEQ_LEGACY`, NCBI BLAST+ reads every `-subject` sequence when it sets up the formatter, in every output format, and LOSAT checks them there too; with it, neither does, and otherwise `BL2SEQ_LEGACY` changes nothing in `-outfmt` `6`, `8`, `9`, `10`, `11`, `12` and `20`, in NCBI as in LOSAT.

## makembindex

`makembindex` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| `-legacy false` with a `-stride` above `-ws_hint` − `-nmer` + 1 (`-stride 20` with the defaults) | Refused when the first sequence of an index volume is read. NCBI's `makembindex` never finishes such a volume. |
| `-legacy false` with a `-stride` of 32768 or more | Refused when the first sequence of an index volume is read. NCBI's index has no defined result. |
| The database sequences, definition lines and damaged taxonomy records that `dustmasker` refuses | Refused, as for `dustmasker`. |
| A database sequence without letters (only a damaged database has one), or, with `-old_style_index true`, one whose every definition line an alias file's `OID_MASK_TYPE` hides | Refused. NCBI's `makembindex` stops with a null-pointer report. |
| A nucleotide FASTA record that NCBI packs into a piece of 4,294,967,293 bases or more: 4,294,967,293 to 4,294,967,295 bases of `A`, `C`, `G` and `T` in a row, or 4,294,967,295 bases in a row with other letters among them (the runs of `N` that NCBI makes of `>?` gap lines count) | Refused. NCBI crashes: it counts the bytes of the piece in 32 bits and gets 0. |
| A FASTA record of 4,294,914,101 letters or more, counting the letters of its `>?` gap lines and modulo 2^32 (`-old_style_index true`) | Refused. NCBI's last chunk of 5,000,000 letters of it ends past 2^32; the end wraps at 32 bits and NCBI reads past the end of the sequence (it crashes). Other records of 4,294,967,296 letters or more are read as NCBI reads them, with their length modulo 2^32. |
| An index volume of more than 4 GiB of packed sequences, which needs a `-volsize` of about 18,000 or more and more than 16 Gbp of sequences in one index volume | Refused when the sequence that passes 4 GiB is read. NCBI writes an index whose sizes and positions wrap at 32 bits. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with an `-iformat blastdb` input whose definition lines lack a member that NCBI's reader requires (only a damaged database has such lines) | Refused where NCBI's error would stop the run, after the same output. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |

## makeprofiledb

`makeprofiledb` runs every option of its `-help` for `-dbtype rps` and `-dbtype delta` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` with a profile, or binary definition lines (a `user` object `ASN1_BlastDefLine`), that lacks a member NCBI's reader requires | Refused where NCBI's error would stop the build. NCBI reports the member as missing and goes on. Other runs with these values are answered as NCBI does. |
| `-dbtype cobalt` | Refused. |
| With `-dbtype delta`: a profile stored by row (`byRow TRUE`), a profile with `numRows < 28` and non-empty `weightedResFreqsPerPos`, an empty `numIndeptObsr` list, a `weightedResFreqsPerPos` list whose length is not a multiple of 28 when `len / 28 == numIndeptObsr.len()` and every full column of 28 sums within `0.0001` of `1.0`, a `weightedResFreqsPerPos` or `numIndeptObsr` value whose scaled conversion to `Uint4` is out of range or `NaN`, or a missing or unwritable temporary directory (`TMPDIR`) | Refused. NCBI's `makeprofiledb` copies into an empty `std::list::begin()`, dereferences or increments past `end()`, performs undefined floating-point-to-unsigned conversions, or stops with a `mkstemp()` error containing its process ID. |
| A query with `iupacna` data holding a byte `>= 0x80` in a pair position (`0..(len & !1)`, before any odd overhang byte) | Refused. NCBI's `makeprofiledb` indexes its `iupacna`-to-`ncbi4na` table with a signed `char` and reads before the table. Bytes `< 0x80`, the final odd overhang byte, and `mol na` (which stops with NCBI's molecule error first) match NCBI. |
| A query in a nucleotide coding with fewer bytes than its length needs: `ncbi2na` of fewer than a byte per 4 bases, `ncbi4na` of fewer than a byte per 2 bases, `iupacna` four letters or more short of an even length (three of an odd one) | Refused. NCBI's `makeprofiledb` reads past the query's data (its database then depends on memory outside the data, or it crashes). A query of enough bytes is built as NCBI builds it, and with frequency ratios only stops with NCBI's conversion error. |
| In a profile with frequency ratios only, a query whose `seq-data` is a `gap` | Refused. NCBI's `makeprofiledb` reads past its table of sequence codings. A query in any other coding gets NCBI's result: its conversion error for a nucleotide coding (`ncbi8na` too), the codes of `ncbi8aa`, NCBI's unknown exception for `ncbipna` and `ncbipaa`. |
| `-max_smp_vol 0` | Refused. NCBI's `makeprofiledb` crashes (division by zero). |
| `-in` naming a directory | Refused. NCBI's `makeprofiledb` never ends. |
| A profile stored by row (`byRow TRUE`) with more than 28 positions | Refused. NCBI's `makeprofiledb` reads its frequency ratios from outside its data. |
| A frequency ratio whose value times 10^9 does not fit a 64-bit integer (such as `1e300`), or a list of frequency ratios shorter than the profile for a profile with frequency ratios only | Refused. NCBI's conversion has no defined result, or it reads past the list. |
| A word score threshold (times the scale factor) that is not a number, a `-scale` beyond 32-bit integers | Refused. |
| A sequence with only a GI id between sequences with other ids, in a version 5 database with `-index true` | Refused with `makeblastdb`'s message. NCBI's `makeprofiledb` aborts; `-blastdb_version 4` builds it. |
| In a profile with frequency ratios only, a query residue that is not one of the 28 protein codes: in `ncbieaa` a character other than a capital letter, `*` or `-` (such as a lower-case letter), in `iupacaa` a character other than a letter, in `ncbistdaa` a code of 28 or more | Refused. NCBI's `makeprofiledb` writes outside its data, and then crashes or gives a result that depends on the residue. A query whose length does not match, or a matrix that NCBI does not know, gets NCBI's error first. |
| A query without Seq-ids and without descriptors (a version 4 database with `-index true` stops with NCBI's message) | Refused. NCBI's `makeprofiledb` aborts (version 5), or with `-index false` writes the `gnl\|BL_ORD_ID` id past the end of an empty list. |
| `-index false` with a query whose first binary definition line has no Seq-id | Refused, as above. |
| Binary definition lines (a `user` object `ASN1_BlastDefLine`) with a `Seq-id` element stated as `NULL` in a version 5 database with no prior sequence written yet and no active taxid map (including `-taxid` or an empty `-taxid_map`) | Refused as a writer-destructor abort after posting `(802.5)`. NCBI's version 5 taxonomy writer aborts on close. With `-taxid_map` (non-empty) or an `org` descriptor, `-blastdb_version 4`, `-index false` when the `NULL` is the first id of the first line, or version 5 after a prior sequence, LOSAT matches NCBI. Damaged binary definition lines and replaced visible characters are read and reported as in NCBI. |
| `:OtherType` pointers in ASN.1 text or `-binary` profiles and binary definition lines | Refused. NCBI's error depends on which ASN.1 classes the process has registered so far. |
| A profile stored by row without scores whose query length is 4,294,967,295 letters (`numColumns -1`) | Refused. NCBI's `makeprofiledb` writes before its matrix. |
| An error after which the database's writer fails again as it closes (a missing input file after the first with an empty `-out` in a version 5 database; a query without definition lines in a version 5 database) | Refused. NCBI's `makeprofiledb` aborts. |

## blastdb_aliastool

`blastdb_aliastool` runs every option of its `-help` (`-xmlhelp` on the command line only), with no program-specific refusals beyond the alias-file refusals above.

## cleanup-blastdb-volumes.py

`cleanup-blastdb-volumes.py` runs every option of NCBI BLAST+ 2.17.0's `cleanup-blastdb-volumes.py` (`-db`, `-dbtype`, `-dry-run`, `-verbose`, `-version`, `-h`/`--help`), matching NCBI's script running under Python 3.13. It is a port of NCBI's Python utility rather than an NCBI C++ Toolkit application, so it does not take C++ Toolkit flags (`-logfile`, `-conffile`, `-dryrun`, `-xmlhelp`, `-version-full`).

For the inputs below, NCBI BLAST+ 2.17.0's Python script exits with a Python traceback and status `1`; LOSAT exits with the same status `1` and writes a concise `Error:` message instead of a traceback:

| Input | What happens |
| :--- | :--- |
| An alias file (`.pal` or `.nal`) that cannot be opened or read, or contains invalid UTF-8 bytes | `Error: cannot read <path>: ...` or `Error: <path>: invalid UTF-8: ...` (exit status `1`), written after any buffered standard output from earlier chunks. `DBLIST` lines in preceding 8 KiB chunks are processed first, as Python's text reader does; if flushing buffered standard output fails and the pending bytes fit in stdout's block size (`st_blksize`), LOSAT also writes Python's `Exception ignored on flushing sys.stdout:` message and exits with status `120`. |
| A `.ncbirc` or `ncbi.ini` file with invalid UTF-8 or a `ConfigParser` syntax, missing-section-header, duplicate-section, duplicate-option, or interpolation error | `Error: <path>: ...` (exit status `1`), where NCBI raises a Python `UnicodeDecodeError` or `configparser.Error` traceback. |
| A matched volume file that cannot be unlinked (for example, a directory named `db.02.pxx` or a read-only directory) | `Error: cannot remove <path>: ...` (exit status `1`), written after the standard output already buffered (and, on a closed pipe where the pending bytes fit in stdout's block size, followed by Python's `Exception ignored on flushing sys.stdout:` message and exit status `120`), where NCBI raises a Python `OSError` traceback. |
| A removed or inaccessible working directory when the database is not found directly or through `BLASTDB` | `Error: cannot determine current working directory: ...` (exit status `1`), where NCBI fails in `os.getcwd()`. |
| A standard output write failure once output reaches the 8 KiB block buffer (or at a newline on a terminal), such as a closed pipe or full disk after 8 KiB of output | `Error: failed to write output: ...` (exit status `1`), stopping after the same file removal as Python's `BrokenPipeError`/`OSError` traceback. Smaller outputs fail on exit flush with Python's `Exception ignored on flushing sys.stdout:` message and exit status `120`. |

## Programs LOSAT does not provide

`blast_formatter_vdb`, `blast_vdb_cmd`, `blastn_vdb`, `tblastn_vdb`, `legacy_blast.pl`, `update_blastdb.pl` and `get_species_taxids.sh`.

## On the list for a later release

- The program `deltablast`.
- `makeprofiledb -dbtype cobalt`.
- `psiblast -phi_pattern` with more than one round, with `-html`, or with output formats other than `0`, `6`, `7`, `10` and `20`, and in
  the library and pylosat.
- The `psiblast` output formats `13` to `16`.
- `blast_formatter`: the output formats other than `6` to `12` and `20`, the tabular formats of the translated programs' archives, and the values of an archive that it does not yet write (in the request or results of `11`, the results of `12`).
- `-import_search_strategy` in `psiblast`, `rpsblast` and `rpstblastn`, and strategy files in binary or XML form.

## Out of scope: they need NCBI's servers

- `-remote` and `-entrez_query`: all programs.
- `-import_search_strategy`: `psiblast`, `rpsblast`, and `rpstblastn`.

## See also

- [Differences from NCBI BLAST+](differences.md)
- [Output formats](output-formats.md)
- [BLAST databases](databases.md)
