# What LOSAT refuses

LOSAT stops with an error for each option and value below. It never runs something else in its place. Use NCBI BLAST+ for these cases, except where an item says that NCBI BLAST+ fails there too.

There are five reasons for a refusal:

- NCBI BLAST+ 2.17.0 crashes, aborts, never finishes, or reads memory outside its buffers, so no output exists to match.
- NCBI's result depends on memory contents, an uninitialised value or the machine.
- The option needs NCBI's servers or its GI and IPG databases.
- The option is not implemented yet.
- The option is not part of what v0.2.0 provides.

A refusal comes before the first byte of output unless an item says otherwise. In several items LOSAT first writes what NCBI writes before it fails (apart from the version line, which names LOSAT), and then stops with the error.

## What you see

LOSAT names the option and the program. There are three forms.

An option refused by name, found while the command line is read (exit status 2):

```text
error: the NCBI BLAST+ option -html is not supported by LOSAT's BLASTP
```

A value or setting refused after the options are read (exit status 1):

```text
Error: -use_sw_tback (the Smith-Waterman traceback) is not supported by LOSAT's PSIBLAST
```

An output format the program does not write (exit status 1):

```text
Error: output format 5 is not supported by LOSAT's PSIBLAST
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
| `-gilist`, `-negative_gilist` | Five search programs, `psiblast` | Out of scope (GI lists) |
| `-ipglist`, `-negative_ipglist` | `blastp`, `blastx`, `psiblast` | Out of scope (IPG lists) |
| `-import_search_strategy` | All search programs | Out of scope (needs NCBI servers) |
| `-html` | All search programs | Not yet supported |
| `-use_sw_tback` | `blastx`, `psiblast`, `rpstblastn` (with CBS), `blastp` (with `-comp_based_stats 0` / `IDENTITY` / `-ungapped`), `tblastn -in_pssm` (when a match is found) | Refused |
| `-ungapped` | `tblastn`, `rpstblastn` (without CBS) | Refused (`rpstblastn` with CBS gives NCBI error) |
| `-use_index true` | `blastn` | Index searches run; refused for partial-volume indexes, two-hit (`-window_size` above 0), `-db_hard_mask`, or damaged/mismatched index volumes |
| `-phi_pattern`, non-listed `-outfmt` | `psiblast` | Not yet supported (`0, 6, 7, 10, 20`; `-in_pssm` / `-in_msa`: `0, 6, 7`) |
| `-mk_counts -infmt seqids` | `windowmasker` | Out of scope (needs NCBI servers) |
| `deltablast`, `makeprofiledb`, `blast_formatter` | Standalone executables | Not yet provided |

## All five search programs

These refusals apply to `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`.

| Option or value | What happens |
| :--- | :--- |
| `-html` | Refused: not yet supported. |
| `-import_search_strategy` | Refused. `psiblast`, `rpsblast` and `rpstblastn` refuse it too. |
| `-export_search_strategy -` (standard output) | Refused where NCBI writes its request into the report: LOSAT holds the report's bytes from NCBI's formatter until NCBI's first request (after its prolog and the first batch of queries is read, before the search), refuses there and writes none of them. A file name works. A run that NCBI ends before that point (an empty query, a database error, an error of an option or of `-outfmt`, a first query record it cannot read, a `-query_loc` past the end, `BATCH_SIZE=abc`, a PSSM run with `-subject`) ends as in NCBI: the held bytes (its prolog) are written, then its error and exit status. NCBI's warnings from before its formatter come before the refusal. |
| `-export_search_strategy` into the file of `-out`, or into standard output | Refused where NCBI writes its request into the report or over it, at the same point as the row above, with no byte written to the stream. A run that NCBI ends before its first request ends as in NCBI, including the held prolog. |
| `-export_search_strategy` with a file that cannot be written to the end | Refused. NCBI writes a message and exits with status 0. |
| `-export_search_strategy` with a file name that is not UTF-8 | Refused. |
| `-export_search_strategy` with `-dbsize` and a list that keeps no sequence of the database | Refused. |
| `-export_search_strategy` in the library API and pylosat | Refused. The command line writes the file. |
| `-parse_deflines` with two subjects of one id | Refused. NCBI cannot report it correctly. |
| `-parse_deflines` with two queries of one id in one query batch | Refused. |
| `-parse_deflines` with a query and a subject of one id that are not the same record | Refused. |
| `-db` with a sequence whose id NCBI cannot report correctly | Refused. The error names the sequence. |
| An alias file with `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `MEMB_BIT` or `MASKLIST` | Refused. See Alias files below. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI` or `SEQ_ID_AVOID_GI` with a database that has parsed ids | Refused. |
| `-xmlhelp`, `-version`, `-version-full`, `-version-full-xml`, `-version-full-json` in the library API and pylosat | Refused. They work on the command line and in LOSAT Web, which prints the command line's text. |
| `-dryrun` | Reads and checks the arguments, then stops. It refuses what a normal run refuses. |
| `-conffile` with a name that is not UTF-8 | Refused. |
| `-logfile` with `tblastx`, `rpsblast` or `rpstblastn` in the WASI command-line build | Refused. That build cannot send its standard error to a file. |
| `-remote`, `-entrez_query`, `-gilist`, `-negative_gilist` | Refused. They need NCBI's servers or GI databases. See the end of this page. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#all-five-search-programs).

## blastn

| Option or value | What happens |
| :--- | :--- |
| `-reward 0` with a penalty other than 0 | Refused (exit status 1). The message is shown below the table. `-reward 0 -penalty 0` (the scores of `-task rmblastn`) runs. |
| `-reward 0 -penalty 0` with `-task megablast` and its greedy extension (without `-no_greedy`; also with `-ungapped` when the gap costs are 0 and 0) | Refused. NCBI divides by zero (exit status 136) or does not finish. |
| `-gapopen` or `-gapextend` above 32767 with `-task megablast` | Refused. |
| `-filtering_db` in LOSAT Web | Refused. LOSAT Web does not read BLAST databases yet. |
| `-use_index true` with a MegaBLAST index (made by `makembindex`) that NCBI searches with and `-db_hard_mask`; or with an index of only some volumes of `-db`, in a later query batch or query chunk where the scan of a volume without an index saves a second initial hit into the list that NCBI frees (NCBI crashes there or reads freed memory) | Refused, after the output of the query batches before. Other searches with an index give NCBI's output, two-hit searches (`-window_size` above 0) and indexes of only some volumes included (there a later query batch finds no seeds in the last index volume searched before each volume without an index that the search reaches, as in NCBI's `-num_threads 1`, at every thread count). When NCBI finds no index, cannot read the header of the first index volume, or the word size is below the index's minimum, LOSAT writes NCBI's error lines and warning and searches without the index, as NCBI does. |
| `-use_index true` with a new-style index volume of 36 to 39 bytes or one without read access, an empty first index volume, or an `-index_name` volume range whose loop NCBI never ends (a first volume of 0 or more with a negative last one, such as `oidx,,0,-1`, or a last one of -1) | Refused. NCBI reads such a header undefined, crashes on the empty volume, or never finishes the loop. |
| `-use_index true` with an index volume whose parts do not fit together (a list, offset or local id outside its sections, an Nmer length of 0, a stride of 32,768 or more), an index not made from the database (a sequence with more chunks than the index has), a seed of the index that starts at or past the last letter of its sequence, or an old-style index (`-index_name`) of fewer sequences than the database when the search runs more than one query batch or query chunk | Refused. NCBI reads past such an index's sections or past the sequence, crashes, or gives a result that depends on the order of its threads. With the old-style index NCBI stops with an error in its second search of the index (a later query batch or query chunk); LOSAT writes the same output of the batches before and refuses there. |
| `-xdrop_gap` or `-xdrop_gap_final` so large that NCBI's 32-bit arithmetic overflows | Refused. See [Rare and pathological inputs](refusals-edge-cases.md#blastn) for the limits. |
| A very large `-xdrop_ungap` (such as `1e10`) that leads to an alignment outside its query context | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |
| `-db_hard_mask` on a database whose sequences are masked to their end, with `-task blastn-short` | Refused. NCBI can crash. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order | Refused. |
| `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused. |
| `-window_masker_taxid` with a counts file whose unit size is 17 to 251 | Refused. NCBI's result is undefined. |

The message for `-reward 0` with the default penalty of `-task megablast` is:

```text
Error: a reward of 0 with a penalty of -2 (NCBI BLAST+'s 16-bit values), for which no query has Karlin-Altschul values, is not supported by LOSAT's BLASTN; reward 0 needs penalty 0 (matrix-only scoring, as -task rmblastn)
```

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastn).

## blastp

| Option or value | What happens |
| :--- | :--- |
| `-use_sw_tback` without composition-based statistics (`-comp_based_stats 0`, and `-matrix IDENTITY`, which NCBI searches without them) | Refused. With `-comp_based_stats` 1 to 3 it works. |
| `-ungapped` with `-use_sw_tback` | Refused. NCBI crashes. |
| `-ungapped` with a `-matrix` name that begins with `IDENTITY` but is not `IDENTITY` | Refused. NCBI aborts. |
| `-ungapped` with a `-matrix` name that NCBI reads from a matrix file | Refused. NCBI aborts or fails. See [Environment](environment.md). |
| `-ungapped` with a protein subject of more than 5,000,000 residues, from `-subject` or from a database | Refused. NCBI searches it in overlapping parts, which LOSAT does not do yet. |
| `-comp_based_stats` with unified P-values (`1u`, `2u`, `3u`), when a query with a residue that can be scored is followed in the same batch of queries by one with none (all `X`, empty, or masked whole by `-seg yes`) and `-searchsp` is not given or is 0 | Refused. NCBI reports every e-value of the first of the two as 0. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7; 5 is the default of `-task blastp-fast`) | Refused. NCBI crashes. |
| `-threshold` with the compressed lookup table, when the neighboring words fill NCBI's 1024 overflow banks | Refused. NCBI corrupts its heap. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |
| `-db_hard_mask` on a database with more than 5,000,000 letters between two masks | Refused. NCBI searches them in overlapping parts, which LOSAT does not do yet. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastp).

## blastx

| Option or value | What happens |
| :--- | :--- |
| `-query_gencode 32` | Refused (exit status 2). Code 32 is not one of NCBI's codes. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7; 5 is the default of `-task blastx-fast`), or a name that begins with `IDENTITY` | Refused. NCBI aborts. |
| An `-evalue` that is not finite (`+inf`, `1e999`, `-nan`) | Refused. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| `-query_loc` or `-subject_loc` with `-outfmt` 5, 8, 9 and 11 to 16 | Refused. Formats 0, 6, 7, 10 and 20 work. |
| A protein subject of more than 5,000,000 residues, from `-subject` or from a database | Refused. NCBI searches it in overlapping parts, which LOSAT does not do yet. |
| `-db_hard_mask` on a database with more than 5,000,000 letters between two masks, and the mask checks of `blastn` | Refused, as in `blastp`. |
| `-use_sw_tback` with `-ungapped` (which only `-comp_based_stats 0` allows) | Refused at the first query batch that is searched. NCBI crashes there, whether the batch has hits or not. |
| `-use_sw_tback` on a query batch that NCBI never finishes or crashes on | Refused, as in `blastp` and `tblastn`. See [Rare and pathological inputs](refusals-edge-cases.md#blastx). |

The error names the sequence and exits with status 1: `Error: subject sequence N has M residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX` for the N-th record of `-subject` (counted from 1), and the same sentence starting with `database sequence N` for a database sequence (N is its position in the database, counted from 0).

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastx).

## tblastn

| Option or value | What happens |
| :--- | :--- |
| `-ungapped` | Refused. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7) | Refused, as in `blastp`. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| An `-evalue` of 1.7976931348623157e308 (`DBL_MAX`) or more | Refused. NCBI's `tblastn` crashes on some inputs with it. |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback | Refused. NCBI crashes. |
| `-use_sw_tback` on a query and subject that NCBI never finishes, or on a query batch that NCBI crashes on | Refused, as in `blastp`. See [Rare and pathological inputs](refusals-edge-cases.md#tblastn). |
| `-db_hard_mask` on a database with a sequence whose masks leave a reading frame an empty part at its start, followed by a part that NCBI searches as the whole frame (a mask within the sequence's first two bases, or within the one or two bases after its last whole codon, for example) | Refused. NCBI crashes. See [Rare and pathological inputs](refusals-edge-cases.md). |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |
| `-in_pssm` | See the next section. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastn).

## tblastx

| Option or value | What happens |
| :--- | :--- |
| `-query_gencode 32`, `-db_gencode 32` | Refused (exit status 2). Code 32 is not one of NCBI's codes. |
| An `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62) when an extension passes into the next frame of the query | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the ungapped re-evaluation | Refused. NCBI crashes. |
| `-db_hard_mask` on a database with a sequence whose masks leave a reading frame an empty part at its start, followed by a part that NCBI searches as the whole frame (a mask within the sequence's first two bases, or within the one or two bases after its last whole codon, for example) | Refused. NCBI crashes. See [Rare and pathological inputs](refusals-edge-cases.md). |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |
| `-outfmt 1` to `4` with `-parse_deflines` | See "All five search programs" above. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastx).

## Alias files and database arithmetic

These refusals apply to every search program that reads databases.

| Option or value | What happens |
| :--- | :--- |
| Alias-file keys `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `MEMB_BIT` and `MASKLIST` | Refused. They select sequences by a list, a mask or a membership bit. |
| An alias file whose `OID_MASK_TYPE` hides predicted RefSeq lines, over a sequence with an id that NCBI warns about when it reads the id's accession | Refused where the sequence's lines are read. NCBI writes the warning again at each read. |
| An `OIDLIST` file that is empty or shorter than 5 bytes (no byte after its 4-byte header) | Refused. NCBI crashes on an empty file, and its error for the other names a line of its source. A file whose bits NCBI would read past the memory page after its end is refused too. |
| An alias file whose `FIRST_OID`-`LAST_OID` range lies outside the OIDs that the ranges of the alias files below it keep (it ends before them, or starts past the bytes NCBI keeps for them) | Refused where the database's OID list is first used. NCBI clears bits outside its array of those OIDs there: before the array it crashes, and past it it writes over its own memory, which crashes it unless the bytes happen to fall in unused memory. `blastdbcmd -info` does not use the list and works. |
| Databases of a `-db` list whose total lengths sum past 64 bits (`-outfmt` 0 to 5 and 13 to 16) | Refused. NCBI's arithmetic overflows. |
| Databases of a `-db` list whose numbers of sequences sum past 32 bits (`-outfmt` 5) | Refused. |
| `rpsblast` and `rpstblastn`: a reported hit of a profile that an alias file's `FIRST_OID`, `LAST_OID` and `OIDLIST` leave out | Refused. NCBI's search ignores these keys, and its report does not find the profile: it writes `Unknown` or stops with a null pointer exception. A profile whose id an included copy of its record holds is written as that copy, as NCBI finds it. |
| `rpsblast` and `rpstblastn` over an RPS database of several volumes: a reported hit of a profile whose definition line an alias file's `OID_MASK_TYPE` hides | Refused. NCBI searches each volume without the alias files, and its report does not find the profile: it writes `Unknown` or stops with a null pointer exception. With one volume NCBI drops such a profile's hits, as LOSAT does. |

`TITLE`, `DBLIST`, `NSEQ`, `LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`, `FIRST_OID`, `LAST_OID` and `OIDLIST` work as in NCBI BLAST+, and `VDBLIST` and `MAX_SEQ_LENGTH` are ignored, as NCBI's reader ignores them. In the profile databases of `rpsblast` and `rpstblastn` they act as in NCBI's search of such a database: one volume is searched with all its keys, and each of several volumes with its own length, count and shortest length under the database's `LENGTH` or `STATS_TOTLEN` and `NSEQ` or `STATS_NSEQ` where these are above 0.

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#alias-files-and-database-arithmetic).

## psiblast

| Option or value | What happens |
| :--- | :--- |
| `-phi_pattern` | Refused: not yet supported. |
| `-in_msa` with a PHYLIP header of 0 sequences, or an empty command (a lone `;`) in a NEXUS file | Refused. NCBI crashes (a division by zero, a null pointer). |
| `-in_msa` with a FASTA alignment whose rows have different lengths, where a row that holds the query is shorter than the first row, or the first row has a letter past the end of a shorter row | Refused. NCBI swaps or reads the shorter row past its end, and its result depends on the memory allocator. Make the rows the same length. |
| `-in_msa` with a `secondary-accession` range that runs backwards (`AB000003-AB000001`) or has a number of more than 18 digits | Refused. NCBI loops without end on a backwards range and reads the numbers into a 32-bit integer that overflows; LOSAT reads 64 bits, so a number past 2147483647 is written with its digits, where NCBI's build wraps it. A range of any length is answered. |
| `-in_msa` with a master sequence from which no ungapped Karlin-Altschul parameters can be computed (X, O, `*` or lower-case x alone, or X with gaps: no residue that the parameters can be computed from), or with `-matrix IDENTITY` | Refused. NCBI crashes. |
| `-in_msa` with a master sequence that is the start of a longer sequence above it in the alignment (`-msa_master_idx` names it) | Refused. NCBI reads and writes past its alignment buffers there, and its result depends on the memory allocator. Put the master sequence above the longer one in the file. |
| `-html` | Refused: not yet supported. |
| `-use_sw_tback` | Refused. |
| `-outfmt` other than 0, 6, 7, 10 and 20 | Refused: not yet supported. With `-in_pssm` or `-in_msa`, only 0, 6 and 7 work. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop that overflows NCBI's kernels | Refused, as in `blastp`. This includes later rounds, and `-in_pssm` with the PSSM's Lambda. |
| A `-subject` set whose sequences all have no letters | Refused. NCBI's `psiblast` crashes; `blastp` runs it. |
| `-in_pssm` or `-in_msa` with a matrix or gap cost other than BLOSUM62 11/1 | Refused. |
| `-in_pssm` with an XML PSSM | Refused. |
| `-out_pssm` or `-out_ascii_pssm` with a file that cannot be written | Refused. |
| `-out_pssm` or `-out_ascii_pssm` into the standard output, the standard error or the `-out` file | Refused. Where NCBI stops first with its own error of `-outfmt` (13 or 14 without `-out`, for example), you get that error. |
| `-out_pssm` into a named pipe, a block device or a symbolic link to a missing file | Refused. |
| The library calls `run` or `run_collect` with formats that take two searches, `-num_descriptions` or `-num_alignments`, and PSSM files | Refused. One call per format writes the command line's output. |

`-matrix IDENTITY` gives NCBI's own error.

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#psiblast).

## tblastn with `-in_pssm`

`tblastn -in_pssm` (PSI-TBLASTN) refuses the `tblastn` entries above and the PSSMs and PSSM queries that `psiblast -in_pssm` refuses. It refuses these too, where NCBI BLAST+ 2.17.0 crashes or aborts:

| Option or value | What happens |
| :--- | :--- |
| `-ungapped` | Refused. |
| `-use_sw_tback` | Refused when the search finds an alignment, after the report's header. Depending on the input, NCBI fails its traceback (`BLAST engine error: Traceback failed`), throws an exception, crashes (with `-subject` and with `-db`), never finishes, or writes rows that change from one run to the next. A search that finds nothing runs, as in NCBI. |
| `-outfmt 11` | Refused. NCBI searches with it. |
| `-mt_mode 1` with more than one thread and `-db` | Refused. |
| A PSSM query that `-seg` masks completely | Refused. |
| `-outfmt` other than 0 to 4, 6, 7 and 10 | Refused. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors. |
| `-db_gencode 32` | Refused. It is not one of NCBI's codes; LOSAT's `tblastn` takes it for its own translations. |
| `-in_pssm` in LOSAT Web | Refused. The library API and pylosat take it as an input. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#tblastn-with--in_pssm).

## rpsblast

| Option or value | What happens |
| :--- | :--- |
| `-outfmt` 17, 19 and 21, and 13 and 14 without `-out` | NCBI's own errors. |
| `-html`, `-entrez_query`, `-remote`, `-import_search_strategy` | Refused. |
| `-sum_stats true` | Refused. NCBI's output then depends on memory contents, or it crashes. |
| `-subject` | Refused by the parser (exit status 2; NCBI prints its usage). It is not an `rpsblast` option. |
| `-outfmt 1` to `5` or `13` to `16`, `-outfmt 18` when the profile has a taxid, or `-outfmt 20` with a field that shows more of the profile than its id (`sallseqid`, `sallacc`, `stitle`, `salltitles`, `staxids`, `slen`, `sseq`, `btop`, `nident`, `pident`, `mismatch`, `positive`, `ppos`, as in the default fields, or a taxonomy field when the profile has a taxid), with a profile whose id is `lcl|Query_N`, the generated id of a query of its batch, in the hits | Refused. NCBI takes the query for the entry (the query's title in `Hit_def`, its letters in the identities, its taxonomy in the organism report). `-outfmt 8`, `9`, `11` and `12`, and `-outfmt 18` and `20` otherwise (such as `20 qseqid sseqid sacc evalue`), write NCBI's bytes (outfmt 11 apart from its `ident`: `LOSAT 0.2.0` where NCBI writes `2.17.0+`). |
| A query above 5,000,000 letters | Refused. |
| A database that lists one volume twice (`-db "x x"`) | Refused. The database layer refuses it for every program; NCBI's `rpsblast` searches the volume once. |
| A profile database whose volumes differ in matrix or gap costs | Refused. |
| An empty or inconsistent `.rps` or `.loo` file, or an inconsistent `.aux` file | Refused. NCBI reads past its buffers or crashes. |
| `rpsblast` in the library API and the Web engine | Refused. |
| The profile-database cases named in Alias files above | Refused. |

NCBI's errors for a missing `.aux`, `.loo` or `.rps` file, a wrong magic number or an empty `.aux` come as in NCBI (exit status 3).

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#rpsblast).

## rpstblastn

| Option or value | What happens |
| :--- | :--- |
| `-outfmt` 17, 19 and 21, and 13 and 14 without `-out` | NCBI's own errors. |
| `-use_sw_tback` with composition-based statistics | Refused. |
| `-ungapped` without composition-based statistics | Refused. NCBI's `rpstblastn` crashes after the search. With composition-based statistics it gets NCBI's own error. |
| `-sum_stats true` | Refused. NCBI's output then depends on memory contents, or it crashes. |
| `-html`, `-entrez_query`, `-remote`, `-import_search_strategy` | Refused. |
| `-subject` | Refused by the parser (exit status 2; NCBI prints its usage). It is not an `rpstblastn` option. |
| `-outfmt 1` to `5` or `13` to `16`, `-outfmt 18` when the profile has a taxid, or `-outfmt 20` with a field of the alignment (such as `qstart`, `length` or the default fields) or a field that shows more of the profile than its id (`sallseqid`, `sallacc`, `stitle`, `salltitles`, `staxids`, `slen`, or a taxonomy field when the profile has a taxid), with a profile whose id is `lcl|Query_N`, the generated id of a query of its batch, in the hits | Refused. NCBI stops with an exception (after half a document in BLAST XML, after the header line in `-outfmt 20`) or takes the query for the entry. `-outfmt 8`, `9`, `11` and `12`, and `-outfmt 18` and `20` otherwise (such as `20 qseqid sseqid sacc evalue`), write NCBI's bytes (outfmt 11 apart from its `ident`: `LOSAT 0.2.0` where NCBI writes `2.17.0+`). |
| `-query_gencode 32` | Refused (exit status 2). |
| `rpstblastn` in the library API and the Web engine | Refused. |
| The database refusals of `rpsblast` above, with its alias-file cases | Refused. |

NCBI's errors for a missing `.aux`, `.loo` or `.rps` file and for an empty `.aux` come as in NCBI (exit status 3). The refusals of options and databases come before the first byte of output, as `rpsblast`'s do, also for queries without a valid frame, which NCBI does not search.

## rpsbproc

`rpsbproc` gives the domain hits, sites and motifs of protein and nucleotide queries.

| Input | What happens |
| :--- | :--- |
| A motif match whose result depends on bytes that NCBI reads beyond the residues of a site | Refused. |
| A site or motif at a position outside the query's residues (other than their end) | Refused. NCBI reads memory outside the string. |
| A file name that is not UTF-8 | Refused. |
| A `[datafiles]` entry of a registry file (`rpsbproc.ini` or the file of `-conffile`) | Refused. |
| A `.ncbirc` with a syntax error | Refused. NCBI reports it with the name of a source file of its own build. |
| The environment variables `OBJECTS_SEQ_DESCR_ALLOW_EMPTY`, and `SERIAL_VERIFY_DATA_READ` or `SERIAL_VERIFY_DATA_GET` set to turn the checks off | Refused. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI` or `SEQ_ID_AVOID_GI` with a query that has a gi | Refused. |
| An archive with parts LOSAT does not read | Refused. See [Rare and pathological inputs](refusals-edge-cases.md). |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#rpsbproc).

## makeblastdb

| Option or value | What happens |
| :--- | :--- |
| `-gi_mask`, `-gi_mask_name` | Refused. The error names the option. |
| `-input_type asn1_txt`, `asn1_bin` | Refused. `fasta` and `blastdb` work. |
| `-input_type blastdb` with an empty `-in` and no `-title` | Refused. NCBI reads the title of a database that is not there and crashes. |
| An option value that is not UTF-8 (a `-title`, `-out` or `-metadata_output_prefix` of other bytes) | Refused. NCBI keeps the bytes. |
| A FASTA record of 4,294,967,296 letters or more, counting the letters of its gaps, whose sequence NCBI builds from all those letters | Refused. NCBI counts the length of a record in 32 bits, which wrap there. |
| `[BLAST] LONG_SEQID`, or a `[BLAST] MAX_SEQID_LENGTH` that `makeblastdb` would use with `-parse_seqids` | Refused. They change which ids are parsed. |
| `DIAG_*` variables and other `NCBI_CONFIG_*` entries outside the accepted list | Refused. |
| A non-empty `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` | Refused. It changes which titles are stored. |
| An index file (`.pin`, `.nin`) of a finished build that a failed write (a full disk) left with no bytes, or with a header field that reaches past the end of the file's last 4 KiB page | Refused. NCBI maps an empty file as a null pointer and crashes. What NCBI reads past its mapping depends on its memory. A file cut short inside that page is read as NCBI reads it, with zeros past its end. |

Example:

```text
error: the NCBI BLAST+ option -gi_mask is not supported by LOSAT's makeblastdb (it writes GI-indexed mask files, which LOSAT does not write yet; -mask_data, -mask_id and -mask_desc work)
```

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#makeblastdb).

## blastdbcmd

| Option or value | What happens |
| :--- | :--- |
| A `-db` list in which an alias and another name give the same volume | Refused. |
| FASTA output of a sequence with a tax id, or `-outfmt` with `%S`, `%L`, `%B`, `%K`, `%C` or `%N`, or `-tax_info` with a name field, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. `-outfmt "%a %T"`, other formats that look nothing up, and `-info` work. |
| A failed write of the `-metadata` output | Refused. NCBI's serial stream reports it. |
| `-range ""`, `-taxids abc`, or a `-taxids` or `-taxidlist` id that is empty, blank, too large or not a number | Refused. NCBI crashes or prints a path. |
| `-mask_sequence_with` below -1 | Refused. NCBI crashes. |
| `-ipg` or `-ipg_batch` on a database without PIG files | Refused. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#blastdbcmd).

## dustmasker

`dustmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A title whose text encoding NCBI cannot tell | Refused. NCBI stops with a string error. |
| A database definition line with bytes NCBI reads as `#` | Refused. |
| A database sequence with a tax id, in any `-outfmt`, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. |
| A failed write of a `seqloc_*` or `maskinfo_*` output | Refused. NCBI's serial stream reports it; the other formats go on silently, as NCBI's. |
| A FASTA record of 4,294,967,296 letters or more, counting the letters of its `>?` gap lines | Refused. NCBI counts the length of a record in 32 bits, which wrap there. |

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
| A FASTA record or database sequence of 4,294,967,296 letters or more | Refused. NCBI counts positions in 32 bits. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#windowmasker).

## segmasker

`segmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| Everything that `dustmasker` refuses: titles, database sequences and definition lines, the serial writer, a failed write of a `seqloc_*` or `maskinfo_*` output, and a FASTA record of 4,294,967,296 letters or more | Refused, as for `dustmasker`. |
| A protein database whose residue bytes are above 127 | Refused. `makeblastdb` writes 0 to 27; NCBI's SEG reads larger bytes as negative characters. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#segmasker).

## convert2blastmask

`convert2blastmask` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A command-line word that is not UTF-8 | Refused. |
| A `Blast-db-mask-info` that LOSAT's serial writer cannot write as NCBI does | Refused. |
| A failed write of a `maskinfo_*` output | Refused. NCBI's serial stream reports it. `interval` goes on silently, as NCBI's, as for `dustmasker` and `segmasker`. |

For rare and pathological inputs (where NCBI BLAST+ crashes, loops, overflows, or reads uninitialised memory), see [Rare and pathological inputs](refusals-edge-cases.md#convert2blastmask).

## makembindex

`makembindex` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| `-legacy false` with a `-stride` above `-ws_hint` − `-nmer` + 1 (`-stride 20` with the defaults) | Refused when the first sequence of an index volume is read. NCBI's `makembindex` never finishes such a volume. |
| `-legacy false` with a `-stride` of 32768 or more | Refused when the first sequence of an index volume is read. NCBI's index has no defined result. |
| A command-line word that is not UTF-8 | Refused. |
| The database sequences, definition lines and damaged taxonomy records that `dustmasker` refuses | Refused, as for `dustmasker`. |
| A database sequence without letters (only a damaged database has one), or, with `-old_style_index true`, one whose every definition line an alias file's `OID_MASK_TYPE` hides | Refused. NCBI's `makembindex` stops with a null-pointer report. |
| A FASTA record of 4,294,967,296 letters or more (`-old_style_index true`) | Refused. |
| An index volume of more than 4 GiB of packed sequences, which needs a `-volsize` of about 18,000 or more and more than 16 Gbp of sequences in one index volume | Refused when the sequence that passes 4 GiB is read. NCBI writes an index whose sizes and positions wrap at 32 bits. |

## Programs LOSAT does not provide

`blastdb_aliastool`, `blast_formatter_vdb`, `blast_vdb_cmd`, `blastn_vdb`, `tblastn_vdb`, `blastdbcheck`, `legacy_blast.pl`, `update_blastdb.pl`, `get_species_taxids.sh` and `cleanup-blastdb-volumes.py`.

## On the list for a later release

- The programs `deltablast`, `makeprofiledb` and `blast_formatter`.
- `psiblast -phi_pattern`.
- The `psiblast` output formats other than `0`, `6`, `7`, `10` and `20` (with `-in_pssm` or `-in_msa`: `0`, `6`, `7`).
- `-html`.

## Out of scope: they need NCBI's servers or its GI and IPG databases

- `-remote` and `-entrez_query`: all programs.
- `-gilist` and `-negative_gilist`: all five search programs and `psiblast`.
- `-ipglist` and `-negative_ipglist`: `blastp`, `blastx`, `psiblast`.
- `-import_search_strategy`: all programs that take it.

## See also

- [Differences from NCBI BLAST+](differences.md)
- [Output formats](output-formats.md)
- [BLAST databases](databases.md)
