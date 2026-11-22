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

There is one section for each program, and one for the five search programs together: `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`. Each section opens with a table of the common refusals in plain words. The rare inputs follow under "Rare inputs"; open them if your data is unusual. They are inputs on which NCBI BLAST+ itself crashes, loops or reads memory it should not.

## All five search programs

These refusals apply to `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx`.

| Option or value | What happens |
| :--- | :--- |
| `-html` | Refused: not yet supported. |
| `-best_hit_overhang`, `-best_hit_score_edge` with `blastx`, `tblastn` or `tblastx` | Refused: not yet supported. `blastn`, `blastp` and `psiblast` run them. |
| `-import_search_strategy` | Refused. `psiblast`, `rpsblast` and `rpstblastn` refuse it too. |
| `-export_search_strategy -` (standard output) | Refused where NCBI writes its request into the report: at NCBI's first request (after its prolog and the first batch of queries is read, before the search). The report until there is NCBI's (its prolog, written before the refusal), and a run that NCBI ends before that point (an empty query, a database error, an error of an option or of `-outfmt`, a first query record it cannot read, a `-query_loc` past the end, `BATCH_SIZE=abc`, a PSSM run with `-subject`, a standard output that fails at the prolog, such as `/dev/full`) ends as in NCBI, with its error and exit status. A file name works. NCBI's warnings from before its formatter come before the refusal. |
| `-export_search_strategy` into the file of `-out`, or into standard output | Refused where NCBI writes its request into the report or over it, at the same point as the row above, after NCBI's prolog. A run that NCBI ends before its first request ends as in NCBI (`-out /dev/full -export_search_strategy /dev/full`: `BLAST failed to write output`, exit status 6). |
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

<details>
<summary>Rare inputs</summary>

Output formats:

- With `-outfmt 1`, `2`, `3` or `4` and `-db`, LOSAT refuses a database sequence whose id is the id NCBI gives a query of the same batch. That id is `lcl|Query_1` for the first query, in any letter case. This applies when the queries are read without `-parse_deflines`.
- This covers the id in the alignment, and the id of the sequence's line in the description table (`lcl|Query_1` of `gi|777|lcl|Query_1`, from `makeblastdb -parse_seqids`).
- NCBI then shows the query in its place and can abort. LOSAT writes the earlier reports and stops.
- With `tblastn -in_pssm`, the id of the query is the local id of the PSSM query (`lcl|myprot`).
- The tabular fields `sacc`, `saccver` and `staxid` are refused for a subject whose title is in a text encoding that NCBI cannot tell. The default fields of formats 6, 7 and 10 contain `saccver`. NCBI makes the ids from that title.
- Such a subject is a `-subject` record whose first title word starts with `Subject_`, or whose id read with `-parse_deflines` is `lcl|Subject_N`. It is also a database sequence without parsed ids, or with an `lcl|Subject_N` id.
- For a database sequence, the title can also come from the taxonomy database: an empty title made from the organism name, or a protein title ending in `[ ]`.
- NCBI stops there with a null pointer exception (exit status 255), or prints an uninitialised number for `staxid`. LOSAT writes what NCBI writes before the batch's report and then stops. Other field lists and formats give NCBI's output (apart from the version line, which names LOSAT).
- `-outfmt 18` with a taxonomy database that names a taxid of a shown group with the text `<@bl_taxid@>` is refused. NCBI writes an uninitialised number. LOSAT stops before the first byte.
- `-outfmt 18` with a `./.ncbirc` that has a syntax error is refused. NCBI's report stops with an exception. LOSAT stops before the first byte.
- `-outfmt 18` with `[NCBI] .Inherits` in the `.ncbirc` is refused. NCBI reads the files it names.
- `-outfmt 18` with a query batch that NCBI does not search is refused. No query of it has a valid context ("Could not calculate ungapped Karlin-Altschul parameters", such as one query of only `N` or `X`).
- There NCBI stops with a null pointer exception (exit status 255). LOSAT writes what NCBI writes before it, then stops with the error.
- A query without a valid context in a batch that NCBI does search gets the report of a query without hits, as in NCBI.
- `-sorthits` is refused when the description table lists every subject and only the last subject has a total score above 99999 bits. NCBI pads that row with spaces without end and never finishes the report. LOSAT writes what NCBI writes before the table, then stops.

Two queries or two subjects with one id (`-parse_deflines`):

- Two queries with one id in one query batch are refused. With `-subject` or `-outfmt 5` this applies to any two batches, because NCBI keeps every batch's queries there. LOSAT first writes the results of the query batches before the second query's batch, as NCBI writes them when it stops at a batch, then stops with the error. An error that NCBI gives for an earlier batch comes first.
- With `-outfmt 11` LOSAT's output differs from NCBI's stop: NCBI ends the BLAST archive with an error record that holds its message (`request-id "Error"`). LOSAT's archive has the earlier batches' records without that record, because the refusal is LOSAT's message, not NCBI's.
- The query batches are NCBI's: like NCBI, LOSAT counts the second query into its batch by the first one's letters, so the second query's batch holds the queries of NCBI's batch, and an error of NCBI's reader comes in the batch where NCBI's reader meets it.
- In the second query's batch, NCBI's errors that do not depend on that query's letters come instead of the refusal, as NCBI gives them. Examples are subjects too short for `-searchsp` ("The average subject length is too short") and a database whose sequence file cannot be read ("Database memory map file error"). This holds where a query of that batch whose id is new in it has letters that NCBI can search. Otherwise the refusal comes.
- The refusal comes also where the second query's ids meet two earlier queries or only some of an earlier query's ids, and in `blastn` where the second query is longer than the first.
- With `-export_search_strategy`, LOSAT writes no search strategy for the refused batch. With one query batch, the file stays empty.
- A query and a subject with one id are refused unless they are one record: the same definition line and the same letters.
- `blastx` and `tblastn` refuse even one record that is both query and subject.
- With `-outfmt 0`, `blastn`, `blastp` and `tblastx` refuse such a record when it has masks.
- With `-outfmt 1` to `4`, `tblastx` writes NCBI's report (apart from the version line, which names LOSAT), whose sections show the record's alignments with itself on the query's rows.
- Where NCBI's merge of such a section throws, `tblastx` writes NCBI's report (apart from the version line, which names LOSAT) up to that query's description table and then NCBI's error, exit status 255. This happens for a repeat whose copies `-query_loc` and `-subject_loc` keep apart, a palindrome or a tandem repeat (`Start iterator already exists for seq2` or `Infinite loop detected`).
- A database sequence with the id of a query of the batch that reports it is refused. With `-outfmt 5` this applies to any batch.
- The same holds with `-outfmt 0` to `4` for `blastp`, `psiblast` and `rpsblast`, and with `-outfmt 1` to `4` for `blastn`, `blastx`, `tblastn` and `tblastx`. The id is then the id of the sequence's line in the description table (`lcl|abc` of `gi|779|lcl|abc`), where NCBI shows the query's title.
- `psiblast` searches one query per batch, so against a database it refuses only a database sequence with the id of the query whose rounds report it.
- With `-in_pssm` (`psiblast` and `tblastn`), a subject or database sequence whose id is the local id of the PSSM query (`lcl|Query_1`) is refused, even the PSSM query's own record.
- A PSSM query whose ids are not one local id is shown with its ids even without `-parse_deflines`. A database sequence with one of them (`sp|P69905.2|` against a PSSM of `swissprot P69905.2`) is refused with or without the option, because NCBI shows the PSSM query for it.
- A `-subject` record with one of those ids is refused under `-parse_deflines`, because NCBI's `tblastn` fails there.

Databases with parsed ids (the error names the sequence). No comparison with NCBI on its swissprot, pdbaa, 18S_fungal_sequences and 16S_ribosomal_RNA was refused. The full list is in [the decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md#databases-with-parsed-seq-ids-version-13-line-dbid).

- An id that NCBI's report code does not find, or finds on another sequence first, is refused in outfmt 0, 5, 6, 7, 10, 13 to 17 and 20. Examples are `gnl|PIG|` and `gim|` ids; in version 4 `bbs|` and `nat|` ids; `gb|X|` after `gb|X.1|`; in version 5 `lcl|300` after `bbs|300`.
- A stored `gnl|BL_ORD_ID|n` is refused when n is not the sequence's position in its volume.
- An id that another sequence also holds in any letter case is refused in the same formats (`lcl|555` and `gb|X.1|` with `lcl|555` in its second definition line). NCBI stops, writes `Unknown`, or writes the other sequence's names.
- The same record in two databases is not refused. With different tax ids it is refused only in outfmt 13 to 16 and in the tabular tax id and name fields, where NCBI writes the first copy's tax ids for both. It is also refused where a title or a taxonomy message depends on them.
- An id that NCBI treats as the id of a query read with `-parse_deflines` in the same batch, though it is spelled differently, is refused. Examples are `sp|P00001.1|NAME` against a query `sp||NAME`, and `gb|X|` against a query `gb|X.1|`.
- This applies in every format when the id is the id in the alignment. When it is the id of the sequence's line in the description table, it applies in outfmt 0 to 4 for `blastp`, `psiblast` and `rpsblast`, and in outfmt 1 to 4 for the other programs.
- NCBI sometimes loads database sequences one at a time while it writes the report: a batch with more alignments than its preload limit, or a limit of 0. There the same refusal applies to the id of any other sequence of the database.
- In version 5 this means any id with a name, or of `pir` or `prf`. It applies in outfmt 0, 5, 6, 7, 10, 13 to 17 and 20. NCBI writes the query's or the other sequence's names.
- Consecutive hits whose ids NCBI takes for one sequence, though they are not copies of one record, are refused (`pdb|1ABC|` and `pdb|1ABC|AA`, `gb|X|` and `gb|X.1|`). This applies in outfmt 0, 5 and 13 to 16 and to the tabular fields `qcovs` and `qcovus`.
- In every format of `psiblast`, the same hits are refused when they reach its `-max_target_seqs` cut.
- Copies of one record (the same definition lines but for their tax ids) that follow each other are shown in one hit, as in NCBI. The PSSM of `psiblast` puts such hits in one row, as in NCBI.
- `psiblast` compares sequences by NCBI's id handles in its inclusion sets, in a PSSM made by its rounds or by `-save_pssm_after_last_round` with `-out_pssm` or `-out_ascii_pssm`, and in the two-part table of outfmt 0.
- There it refuses two sequences whose ids the handles can take for one sequence, such as ids that differ only in letter case, or one accession under `gb`, `emb` and `dbj`. It refuses a `gim|` id in such a PSSM. It refuses two copies of one record with another hit between them, where NCBI crashes.
- Ids and titles that NCBI cannot write are refused, such as an outfmt 0 description id of 65 bytes or more.
- Titles that NCBI makes from taxonomy data are written as NCBI writes them: an empty title from the organism's name, and the end of a protein title ending in `]` rewritten from the first tax id (`alpha, partial [Bacillus sp.]` becomes `alpha [Bacillus sp.]`).
- A version 5 list of databases (`-db "a b"`), some with and some without the accession table of `makeblastdb -parse_seqids`, is refused together with a damaged `taxdb.bti` or `taxdb.btd` that writes text. NCBI's loads of other sequences then write such text or fail on their ids, in all five programs.
- With `-outfmt 5` it is also refused when a tax id is looked up after a batch whose XML report NCBI cuts at a query title or id outside US-ASCII (`blastp`, `blastn`, `tblastn`, `tblastx`).
- NCBI's stop here is `BLAST Database error: DB contains no accession info.` (exit status 2). Whether NCBI writes the texts before the error depends on the order of its loads. LOSAT refuses before the first byte.

Strategy files:

- `-export_search_strategy` is written as NCBI writes it for the five programs, `psiblast` (with a FASTA query or `-in_pssm`), `rpsblast`, `rpstblastn` and `tblastn -in_pssm`, apart from its `ident` (`LOSAT 0.2.0` where NCBI writes `2.17.0+`). For `tblastn -in_pssm` NCBI's request is a PSI-BLAST request that holds the PSSM as read, with `-db`. With `-subject` NCBI gives its error after the prolog.
- `rpsblast` and `rpstblastn` write the file of NCBI's `-num_threads 1` at every `-num_threads` and `-mt_mode`. With `-mt_mode 1` and more than one thread, NCBI's threads write their requests into the file at the same time, and for a database of several volumes they leave out the search space of the batch before.

</details>

## blastn

| Option or value | What happens |
| :--- | :--- |
| `-reward 0` with a penalty other than 0 | Refused (exit status 1). The message is shown below the table. `-reward 0 -penalty 0` (the scores of `-task rmblastn`) runs. |
| `-reward 0 -penalty 0` with `-task megablast` and its greedy extension (without `-no_greedy`; also with `-ungapped` when the gap costs are 0 and 0) | Refused. NCBI divides by zero (exit status 136) or does not finish. |
| `-gapopen` or `-gapextend` above 32767 with `-task megablast` | Refused. |
| `-filtering_db` in LOSAT Web | Refused. LOSAT Web does not read BLAST databases yet. |
| `-use_index true` with a megablast search that NCBI runs with a MegaBLAST index (made by `makembindex`): a new-style index of a volume of `-db`, an old-style index of `-index_name` or of the `-db` name, at a `-word_size` the index supports | Refused. When NCBI finds no index, cannot read the header of the first index volume, or the word size is below the index's minimum, LOSAT writes NCBI's error lines and warning and searches without the index, as NCBI does. |
| `-use_index true` with a new-style index volume of 36 to 39 bytes or one without read access, an empty first index volume, or an `-index_name` volume range whose loop NCBI never ends (a first volume of 0 or more with a negative last one, such as `oidx,,0,-1`, or a last one of -1) | Refused. NCBI reads such a header undefined, crashes on the empty volume, or never finishes the loop. |
| `-xdrop_gap` or `-xdrop_gap_final` so large that NCBI's 32-bit arithmetic overflows | Refused. See Rare inputs for the limits. |
| A very large `-xdrop_ungap` (such as `1e10`) that leads to an alignment outside its query context | Refused. See Rare inputs. |
| `-db_hard_mask` on a database whose sequences are masked to their end, with `-task blastn-short` | Refused. NCBI can crash. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order | Refused. |
| `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused. |
| `-window_masker_taxid` with a counts file whose unit size is 17 to 251 | Refused. NCBI's result is undefined. |

The message for `-reward 0` with the default penalty of `-task megablast` is:

```text
Error: a reward of 0 with a penalty of -2 (NCBI BLAST+'s 16-bit values), for which no query has Karlin-Altschul values, is not supported by LOSAT's BLASTN; reward 0 needs penalty 0 (matrix-only scoring, as -task rmblastn)
```

<details>
<summary>Rare inputs</summary>

- `-outfmt 1` to `4` with `-parse_deflines` and a query record that is also a subject (one record) is refused in one case, shared with `blastp`. The first of the longest aligned segments of the query's alignments (in their order) is then in an ungapped HSP of the record with itself that is not its first.
- NCBI leaves the record's subject rows out of its merged alignment, writes past its vectors and aborts (`double free or corruption`). LOSAT stops before the first byte of the report. For `blastn` it writes the prolog first, as NCBI does.
- In a run of many queries that LOSAT searches in blocks, it stops before the first byte of the block that holds the query, after the reports of the blocks before.
- A final gapped X-drop of dynamic programming (`-task blastn`, `blastn-short`, `dc-megablast`, `-no_greedy`) is refused when its raw score divided by the gap extension reaches 1,073,741,721. NCBI's arithmetic overflows there.
- A gapped X-drop of greedy extension (`-task megablast`) is refused where NCBI's greedy alignment reads or writes outside its arrays. NCBI then corrupts its heap, reads memory it did not write, crashes, or stops with `BLAST engine error: Unknown error code -1`.
- With the default scores (reward 1, penalty -2), `-xdrop_gap` or `-xdrop_gap_final` is refused from about -7.4 to -1.98e9 bits and from about 1.98e9 to 3.97e9 bits. It is also refused from about -1.98e9 to -3.97e9 bits, unless the other X-drop is as large in NCBI's doubled arithmetic.
- The limits move with the scores. With an even reward, a value past 32-bit integers is refused too.
- An alignment that NCBI keeps outside its query context and that the report reaches is refused. It comes from an ungapped extension that runs past an end of its context, such as with `-xdrop_ungap 1e10`, which NCBI turns into INT_MIN.
- The same refusal applies with `-ungapped` when the re-evaluation keeps part of an alignment outside the context, and with `-task megablast` when a greedy extension from the middle of such an extension stays outside it.
- No e-value, percent identity, query coverage, `-max_hsps`, `-max_target_seqs`, `-culling_limit` or `-subject_besthit` limit drops such an alignment first. LOSAT drops it where NCBI does.
- Where its query coordinates leave the query on the left (a start below 1), NCBI's report stops with an exception (exit status 255). Where they pass the query's end, NCBI's reports show them. LOSAT reproduces neither yet and stops before the rows of the query batch.
- Where the part outside the context is below the cutoff, or the greedy extension comes back into the context, LOSAT reports the alignment as NCBI does.
- `-db_hard_mask` is refused on a database that has a searched sequence of at least the word's letters whose first mask with letters runs to its end, when the word has fewer than 11 letters (`-task blastn-short`). That is a sequence that is all masked, or has a masked tail such as a poly-A tail.
- Further conditions for that refusal: at most 5,000,000 letters before the mask, and no other mask before it except masks of no letters at the start. NCBI then scans the masked letters past its first piece and can crash.
- With a word of 11 letters or more, such a sequence is searched before its mask; see [Differences](differences.md).
- `-db_hard_mask` is also refused for a database whose masks are stored out of order or overlap.
- It is refused too when the last mask is shorter than the 1 to 3 letters before its end that the last piece starts with, and that piece and an earlier piece that ends after its start are both long enough to hold a word.
- NCBI then merges their alignments by a value it reads past its list of pieces. This does not apply where the stretch after the mask is longer than 5,000,000 letters.
- LOSAT refuses before it writes anything, and the error names the first such sequence by its position in the database.
- `-db_soft_mask` or `-db_hard_mask` on a volume without the algorithm (a database whose volumes give the algorithm differently), or on a volume whose masks have only their big-endian file (`.nab` without `.nac`), is refused. NCBI reads masks of another algorithm or reversed numbers.
- `-window_masker_taxid` with an optimized ascii counts file is refused when the third line, which holds the hash parameters, leaves a value unread that NCBI then checks uninitialised.

</details>

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
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused. See Rare inputs. |
| `-db_hard_mask` on a database with more than 5,000,000 letters between two masks | Refused. NCBI searches them in overlapping parts, which LOSAT does not do yet. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |

<details>
<summary>Rare inputs</summary>

- `-outfmt 1` to `4` with `-parse_deflines` and a query record that is also a subject is refused in the case described under `blastn`.
- A volume without the algorithm can come from an alias or a `-db` list whose databases name the algorithm with different texts. NCBI numbers the algorithms over all volumes and can crash with composition-based statistics there.
- `-matrix` with gap costs: every matrix with any gap costs of NCBI's table for the matrix works. `IDENTITY` is refused with the compressed lookup table.
- A gapped X-drop whose raw score makes the sizes of NCBI's gapped alignment kernels leave 32-bit integers is refused. With BLOSUM62 11/1 this is from about 4.14e8 to 8.27e8 bits.
- This applies to `-xdrop_gap` at every `-comp_based_stats`, and to `-xdrop_gap_final` with `-comp_based_stats 0`. The redo of modes 1 to 3 runs with its own X-drop.
- A larger value leaves 32-bit integers itself and becomes NCBI's INT_MIN, which the kernels raise to the gap costs. Those searches run as NCBI's.
- A query batch whose search fails after a first batch that NCBI does not search (its queries give no Karlin-Altschul parameters) is refused in some cases. `-searchsp` or `OLD_FSC` with subjects of fewer letters than sequences, a database without sequences, and a database whose sequence file cannot be read, are examples.
- For `-outfmt` 0 to 4, 6, 7 and 10, `blastp` writes NCBI's report of that first batch (apart from the version line, which names LOSAT) and then the error. With the other formats it refuses before the first byte.
- With `-use_sw_tback`, a query and subject whose Smith-Waterman alignment NCBI never finishes is refused. Only a large `-evalue` or a small search space (`-searchsp`, or a small query and database) can reach it.
- With `-use_sw_tback`, a query batch on which NCBI crashes is refused: one whose last Smith-Waterman alignment fails after the batch's first query has kept one. Where the first query has kept none, the batch ends with NCBI's `BLAST engine error: Traceback failed`, and with LOSAT's `LOSAT engine error: Traceback failed`. A later match of a database sequence that `-negative_seqidlist` keeps without any of its definition lines (NCBI skips its alignment) clears that failure, as in NCBI, unless by then every query of the batch already holds `-max_target_seqs` alignments with e-values of 0.002 or less and the match's preliminary e-value is above 0.01.
- For both, LOSAT writes the results of the query batches before and then stops with the error.
- With `-ungapped -window_size 0` and a low `-threshold`, `-outfmt 1` to `4` is refused for a one-hit ungapped alignment of length 0 among the alignments shown, where NCBI's merge of the query-anchored alignments frees memory that it still reads. NCBI crashed in every such run we made but one; LOSAT writes that one. Every other format writes NCBI's report of such an alignment (apart from the version line, which names LOSAT).
- With `-ungapped` and a database, the 5,000,000-residue limit is checked against the longest sequence of the database, also when a sequence list or an alias file leaves that sequence out of the search.

</details>

## blastx

| Option or value | What happens |
| :--- | :--- |
| `-use_sw_tback` | Refused. |
| `-query_gencode 32` | Refused (exit status 2). Code 32 is not one of NCBI's codes. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7; 5 is the default of `-task blastx-fast`), or a name that begins with `IDENTITY` | Refused. NCBI aborts. |
| An `-evalue` that is not finite (`+inf`, `1e999`, `-nan`) | Refused. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| `-query_loc` or `-subject_loc` with `-outfmt` 5, 8, 9 and 11 to 16 | Refused. Formats 0, 6, 7, 10 and 20 work. |
| A protein subject of more than 5,000,000 residues, from `-subject` or from a database | Refused. NCBI searches it in overlapping parts, which LOSAT does not do yet. |
| `-db_hard_mask` on a database with more than 5,000,000 letters between two masks, and the mask checks of `blastn` | Refused, as in `blastp`. |

The error names the sequence and exits with status 1: `Error: subject sequence N has M residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX` for the N-th record of `-subject` (counted from 1), and the same sentence starting with `database sequence N` for a database sequence (N is its position in the database, counted from 0).

<details>
<summary>Rare inputs</summary>

- With `-ungapped -window_size 0`, a low `-threshold` and `-subject` (also with `-subject_loc`), LOSAT refuses in every format but `-outfmt 18` a one-hit ungapped alignment of length 0. It sits at the first letter of the searched part of a subject that has alignments of letters among the hits that the reports show.
- NCBI may stop earlier at another alignment of length 0. LOSAT refuses the batch before its first byte.
- NCBI's Seq-loc mapper cuts such an alignment into two alignments without scores. NCBI then writes rows of the query's id with an E-value of 0.0 and a bit score of -1.0, stops with an exception, or runs out of memory, depending on the frame, the fields and the memory limit.
- The same refusal applies to a subject whose every such alignment has length 0, one of them at the first or second letter of the searched part. The mapper, finding no sequence type, then makes alignments without scores of query nucleotides against the subject's letters.
- NCBI then writes more rows, throws exceptions or runs out of memory when there are two or more. The tabular reports of a subject's only one are written as NCBI writes them.
- Elsewhere LOSAT reports alignments of length 0 as NCBI does.
- With `-ungapped -window_size 0` and a low `-threshold`, `-outfmt 1` to `4` is refused for a one-hit ungapped alignment of length 0 among the alignments shown, where NCBI's merge of the query-anchored alignments frees memory that it still reads. NCBI crashed in every such run we made.
- One case is an alignment whose query row starts at the first nucleotide of a codon of an alignment in the query's reading frame +1. The other is an alignment followed on its query row by unaligned letters, whose segment the merge frees unless it is the last segment it made.
- LOSAT stops before the first byte of the report, unless the report of a query before it ends with NCBI's exception.

</details>

## tblastn

| Option or value | What happens |
| :--- | :--- |
| `-ungapped` | Refused. |
| `-use_sw_tback` | Refused. |
| `-matrix IDENTITY` with the compressed lookup table (word size 5 to 7) | Refused, as in `blastp`. |
| `-threshold` with the compressed lookup table | Refused, as in `blastp`. |
| A gapped X-drop so large that NCBI's alignment kernels overflow | Refused, as in `blastp`. |
| An `-evalue` of 1.7976931348623157e308 (`DBL_MAX`) or more | Refused. NCBI's `tblastn` crashes on some inputs with it. |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback | Refused. NCBI crashes. |
| `-db_hard_mask` on a database with a sequence whose masks leave a reading frame an empty part at its start, followed by a part that NCBI searches as the whole frame (a mask within the sequence's first two bases, or within the one or two bases after its last whole codon, for example) | Refused. NCBI crashes. See Rare inputs. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |
| `-in_pssm` | See the next section. |

<details>
<summary>Rare inputs</summary>

- A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback makes NCBI read the translation of the subject's frame, which it never makes for a subject this short. NCBI then crashes.
- In practice this happens with `-window_size 0`: no frame of such a subject, three letters at most, holds two words that do not overlap.
- With `-window_size 0` (the one-hit word finder) and a low `-threshold`, LOSAT refuses the cases below.
- Two gapped starts of one-hit alignments of negative length that lie at the end of their 5,000,000-letter part of a frame, or past it. NCBI never finishes on them.
- With `-in_pssm`, a one-hit alignment whose window of 11 letters passes the PSSM's last position. NCBI reads past its matrix.
- With `-in_pssm` and a PSSM whose `-` column scores above NCBI's minimum, a one-hit alignment whose gapped start lies past the end of NCBI's translation of the subject (after its last frame). NCBI's gapped alignment reads the byte after its end, which it never sets, or past it.
- A one-hit alignment whose gapped start lies after the end of its frame's translation and that reaches NCBI's traceback. NCBI reads a byte its translation does not set there.
- The report of a one-hit alignment whose traceback starts on the end marker. NCBI places its column past the subject's end, and LOSAT does not reproduce that result yet. Such alignments that NCBI drops before its reports give NCBI's output (apart from the version line, which names LOSAT).
- `-db_hard_mask` is refused on a database that has a searched sequence (of at least the word size times 3 plus 2 bases) whose masks leave a reading frame an empty part at its start, followed by a part that NCBI searches as the whole frame. NCBI searches each of the six reading frames in parts: it maps each stretch between the sequence's masks (or its ends) to the frame by dividing its positions by 3, counting the minus frames back from the end of frame -1, and searches each mapped stretch as a part, one that maps to no letter included. It searches the last part that starts at the frame's first letter as the whole frame when that part has at most 5,000,000 letters and no part follows it or the next part starts exactly at the frame's end. After an empty part at the frame's start, as after a mask within the sequence's first two bases or within the one or two bases after its last whole codon (for example [1200, 1202) of 1,202 bases), NCBI then frees memory that it did not allocate and crashes in every run. The message names the sequence, the mask and the frame.
- Where that part, from the frame's first letter, ends before the frame's last letter instead (masks at or near either end of a sequence), LOSAT searches the part alone; see [Differences](differences.md).

</details>

## tblastx

| Option or value | What happens |
| :--- | :--- |
| `-query_gencode 32`, `-db_gencode 32` | Refused (exit status 2). Code 32 is not one of NCBI's codes. |
| An `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62) when an extension passes into the next frame of the query | Refused. See Rare inputs. |
| A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the ungapped re-evaluation | Refused. NCBI crashes. |
| `-db_hard_mask` on a database with a sequence whose masks leave a reading frame an empty part at its start, followed by a part that NCBI searches as the whole frame (a mask within the sequence's first two bases, or within the one or two bases after its last whole codon, for example) | Refused. NCBI crashes. See Rare inputs. |
| `-db_hard_mask` on a database whose masks overlap or are stored out of order, and `-db_soft_mask` or `-db_hard_mask` on a volume with masks but without the algorithm, or with only the big-endian mask file | Refused, as in `blastn`. A mask that ends past its sequence (`makeblastdb` never writes one) is refused as well. |
| `-outfmt 1` to `4` with `-parse_deflines` | See "All five search programs" above. |

<details>
<summary>Rare inputs</summary>

- With `-window_size 0` and a low `-threshold`, a one-hit ungapped alignment of length 0 or less is refused when it meets another in the overlap of two 5,000,000-letter parts of a subject frame and their lengths add up to 0. NCBI merges them by dividing by that sum, and its result is undefined.
- An alignment on a subject (or its `-subject_loc` range) of fewer than 12 bases that reaches the ungapped re-evaluation makes NCBI read the translation of the subject's frame and crash.
- `tblastx` first writes in full the reports of the query batches before. NCBI writes them before it crashes, except what its output buffer still holds (with `-outfmt 0`, the end of the last one).
- An ungapped extension that passes the end of a query frame into the next frame of the query batch and scores higher there is refused, with an `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62; a raw X-drop above 32,768).
- In the cases checked, NCBI keeps such an alignment in the frame where it started, with query coordinates outside that frame. It writes coordinates past the query or leaves the alignment out of every report.
- LOSAT writes the results of the query batches before and then stops with the error.
- An extension that passes into the next frame without scoring higher there gives NCBI's output (apart from the version line, which names LOSAT).
- `-db_hard_mask` is refused for the sequences described under `tblastn`. A part that ends before its frame's end is searched alone, as in `tblastn`; see [Differences](differences.md).

</details>

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

<details>
<summary>Rare inputs</summary>

- `blastn` and `tblastx` refuse an alias file's `NSEQ 0`, and its `LENGTH` or `STATS_TOTLEN` without `NSEQ` or `STATS_NSEQ` over volumes without sequences. NCBI crashes on them where a query is searched.
- `blastp`, `blastx`, `tblastn` and `psiblast` refuse the same cases with `-searchsp` or `OLD_FSC`. `blastx` also refuses them with `-ungapped`.
- A run of such a database whose queries have no letters to search is reported as in NCBI by every program.
- `rpsblast` and `rpstblastn`: an alias file's `LENGTH` and `NSEQ` over one volume also set the length that NCBI's search covers, their sum, which can be shorter or longer than the profiles. Refused: a sum above 2,147,483,647 (NCBI frees the profiles twice and aborts), or one so long that NCBI's diagonal table passes 2^30 entries (its loop never ends).
- A sum shorter than the profiles is refused where two gapped alignments of a query (of a frame for `rpstblastn`) start past it plus one. NCBI's interval tree then splits without end. With one such alignment the search runs as in NCBI.
- A sum longer than the profiles is refused where an ungapped extension runs past the row after the last profile. This needs an X-drop above the profiles' delimiter scores and a sum of the profiles plus 2 or more (NCBI reads past its score matrix). At plus 1 the search runs as in NCBI where that row lies within the `.rps` file's last memory page, and is refused where it lies past it (NCBI reads past the file's mapping).
- `rpsblast` and `rpstblastn` with `-searchsp` or `OLD_FSC` refuse an `NSEQ 0` over one volume where a query is searched (NCBI crashes). A run of queries without a letter to search is reported as in NCBI.

</details>

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
| `-in_pssm` with a binary or XML PSSM | Refused. |
| `-out_pssm` or `-out_ascii_pssm` with a file that cannot be written | Refused. |
| `-out_pssm` or `-out_ascii_pssm` into the standard output, the standard error or the `-out` file | Refused. Where NCBI stops first with its own error of `-outfmt` (13 or 14 without `-out`, for example), you get that error. |
| `-out_pssm` into a named pipe, a block device or a symbolic link to a missing file | Refused. |
| The library calls `run` or `run_collect` with formats that take two searches, `-num_descriptions` or `-num_alignments`, and PSSM files | Refused. One call per format writes the command line's output. |

`-matrix IDENTITY` gives NCBI's own error.

<details>
<summary>Rare inputs</summary>

A PSSM given with `-in_pssm` is refused when its input or content is one of these. In each case NCBI crashes, does not finish, reads out of its buffers, goes on with a wrong value, or depends on the machine.

- A PSSM of frequency ratios only, written by row. NCBI reads the scores it computes for it in the wrong order.
- A PSSM whose position-based scaling gives score frequencies but no Karlin-Altschul parameters. NCBI goes on with a Lambda of -1.
- A PSSM with a negative gapped Karlin-Altschul value (`lambda`, `kappa` or `h`). NCBI goes on with a block it holds invalid. A negative ungapped value (`lambdaUngapped`, `kappaUngapped` or `hUngapped`) gets NCBI's error of a query without Karlin-Altschul parameters, as in NCBI.
- A PSSM of frequency ratios only from which NCBI's PSSM engine makes no scaled scores. NCBI does not finish or stops with an error.
- A query without a Seq-id, with a Geninfo import id (`gim|`, which NCBI cannot resolve), or with a database ordinal (`gnl|BL_ORD_ID|n`).
- A query with a local Seq-id among other Seq-ids (NCBI's reports then show it as no FASTA query is shown), with an empty local id, or with a Seq-id that has neither an accession nor a name.
- For the last kind, NCBI's run depends on the database and the report: a database error with some databases, `Failed to resolve SeqId` in the pairwise report, `Unknown` in tabular reports.
- A query with more than one title, or with white space at the ends of its title.
- A score outside -32,768 to 32,767.
- A PSSM whose columns' largest scores add up to 2^30 or more. Every PSSM of fewer than 19,800 columns passes.
- A query whose letters give no Karlin-Altschul parameters. NCBI crashes.
- A PSSM without a query, and query Seq-data other than `iupacaa`, `ncbieaa` and `ncbistdaa`.
- A query whose annotations (`annot`) give a case of NCBI's object manager whose output names a memory address or depends on the machine. These cases are an empty location (`Empty region in 0x...`), a NULL-pointer report, a Dense-seg that NCBI divides by zero, malformed BitMagic data of a Seq-table, reads past a Packed-seg's vectors, and integer overflows. BitMagic data is decoded as NCBI decodes it; refused are the data that NCBI reads past (such as empty or cut data), decodes with undefined behaviour (such as an XOR block or a write past BitMagic's blocks) or from memory it did not write, and a GAP block that is not one where NCBI's scalar, SSE2 and SSE4.2 builds decode it or answer a row of it differently.
- Also two Seq-ids of one location (a feature's, a Seq-table's or a Seq-table row's) with different messages that NCBI orders by address; a Seq-table of more than 10,000,000 rows; a Std-seg of more than 1,000,000 rows; and XML and JSON reports of an untitled query whose annotations hold features.
- An error about the query's annotations, or a Seq-id twice in the entry. Two Seq-ids count as one when NCBI's object manager takes them as one: equal ones, and ones that differ only where NCBI does not look, such as the case of letters, GenBank, EMBL or DDBJ for an accession alone, or PDB's `chain` and `chain-id`.
- A query Seq-entry set whose annotations NCBI's object manager takes, and which NCBI then stops at with its exception (`CSeq_entry::GetSeq()`).
- These last two are not refused where NCBI stops before the search and writes its prolog and error: `-export_search_strategy` with `-subject`, a PSSM that NCBI's checks reject, or a PSSM of frequency ratios only whose scores NCBI's PSSM engine cannot make.
- With `-subject`, a Seq-id twice in the query's entry where NCBI's exception may write a Seq-id as a subject record spelled it. That is a GenBank, EMBL or DDBJ accession alone, a name, a patent number or country with a letter, or more than 64 letters in a local or general id or in an accession.
- With `-out_pssm`, a PSSM of frequency ratios only without `rpsdbparams`. NCBI cannot write its matrix name.
- With `-out_ascii_pssm`, a PSSM whose information content comes without the other columns of NCBI's long ASCII form.
- With `-out_ascii_pssm`, a PSSM query with a residue that its coding does not have. NCBI prints its letter from memory it never set.
- A real number in a PSSM whose exponent overflows NCBI's 32-bit `int` (`1.55e-2147483647`) is read with the exact exponent, the smallest normal double; see [Differences](differences.md). The PSSM's ASN.1 text is read as NCBI reads it.

</details>

## tblastn with `-in_pssm`

`tblastn -in_pssm` (PSI-TBLASTN) refuses the `tblastn` entries above and the PSSMs and PSSM queries that `psiblast -in_pssm` refuses. It refuses these too, where NCBI BLAST+ 2.17.0 crashes or aborts:

| Option or value | What happens |
| :--- | :--- |
| `-ungapped` | Refused. |
| `-use_sw_tback` | Refused. |
| `-outfmt 11` | Refused. NCBI searches with it. |
| `-mt_mode 1` with more than one thread and `-db` | Refused. |
| A PSSM query that `-seg` masks completely | Refused. |
| `-outfmt` other than 0 to 4, 6, 7 and 10 | Refused. `-outfmt` 17, 19 and 21, and 13 and 14 without `-out`, get NCBI's own errors. |
| `-db_gencode 32` | Refused. It is not one of NCBI's codes; LOSAT's `tblastn` takes it for its own translations. |
| `-in_pssm` in LOSAT Web | Refused. The library API and pylosat take it as an input. |

<details>
<summary>Rare inputs</summary>

- Where NCBI stops before the search, as for `psiblast -in_pssm`, the output formats (every format `tblastn` writes, `-outfmt 11` included) and the query's titles are not refused. NCBI's prolog and error are written.
- NCBI's object manager can end the run as it adds the query's Seq-entry to its scope. LOSAT does not refuse what NCBI's `tblastn` makes after that point either: its formatter, its requests and its search.
- Examples of those later limits: an `-xdrop_gap` or `-xdrop_gap_final` that NCBI's kernels cannot size; `-dbsize` with `-outfmt 11` and a list that keeps no sequence; a database named `n/a` or `SRA` in a Seq-annot format; database lengths whose sums the report cannot write.
- NCBI's messages and error are written there, with `-outfmt 11` in its error archive.

</details>

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

<details>
<summary>Rare inputs</summary>

- An X-drop whose raw score, or its product with the database's scale, leaves 32-bit integers is refused. NCBI's conversion is undefined.
- A gapped X-drop whose raw score over the scaled gap extension reaches 1,073,741,721 is refused (about 4.1e8 bits, at scales 1 and 2 only). The sizes of NCBI's alignment kernels then overflow.
- A `.freq` file whose rows give a profile no Karlin-Altschul parameters is refused. NCBI goes on with a Lambda of -1.
- A profile database whose profile can score above 2,147,483,647 (the sum of its rows' best scores) is refused. NCBI's 32-bit scores overflow.
- A profile whose Smith-Waterman redo (`-use_sw_tback`) NCBI repeats without end is refused. This needs three conditions. The e-value of a score of 0 (the query's search space times K) passes `-evalue`. The profile has two or more preliminary HSPs. Its first position scores the query's first letter above 0. NCBI stops only when an allocation fails.
- A profile database whose matrix and gap costs are not BLOSUM62 11/1, or whose scale is not a whole number, is refused when a scaled gap cost, X-drop or score cutoff then leaves 32-bit integers (NCBI's products overflow).
- It is also refused when NCBI's Lambda changes from one e-value to the next at that scale (such as 531, where NCBI divides Lambda by the scale and multiplies it back).

</details>

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
| An archive with parts LOSAT does not read | Refused. See Rare inputs. |

<details>
<summary>Rare inputs</summary>

- A motif whose check makes NCBI hold 1 GiB or more at once is refused. This is its minimal map's string with the old buffer while it grows, or its match records with their copies, in allocations below 2^47 bytes.
- Whether NCBI's allocations succeed depends on the host's memory. A repeat that does not advance, such as `<*`, makes NCBI's records grow until the memory runs out. These refusals name the site and quote its motif.
- A motif match that may read beyond NCBI's buffer of the site's residues is refused.
- In an archive, a VisibleString with a byte outside the visible ASCII range is refused. NCBI replaces it by `#` and posts an error. The DefLine, a UTF8String, takes any byte.
- In an archive, a member named by its tag (`[0]`) is refused.
- In an archive, Seq-align segments other than `Dense-seg`, `Std-seg` and `Packed-seg` (`dendiag`, `disc`, `spliced`, `sparse`) are refused.
- A `Std-seg` alignment whose query id NCBI would take from a first location of the kind `packed-int`, `packed-pnt`, `mix`, `equiv`, `bond` or `feat` is refused. LOSAT does not look up those ids.
- Residues in NCBI8aa and NCBIstdaa codes from 28 are refused.
- A query without a title but with `annot` is refused.
- A mask location other than a `packed-int` of one Seq-id is refused.
- Data lines on which NCBI's loader reads outside the line are refused: `cdtrack.txt` lines that end without a root or after a field it skips.
- For nucleotide queries, a first `Std-seg` with fewer than 2 ids is refused.
- A `Std-seg` without locations, or whose only location is not `null`, is refused. One whose only location is `null` is skipped, as in NCBI.
- A minus-strand alignment with a subject location other than `int` or `pnt` is refused. NCBI reads past its segment lengths.
- A query whose hits carry sites or motifs and whose bases, or the IUPACna letters of their reverse complement, are fewer than 3 is refused. NCBI's `TranslateAll` erases from an empty string and crashes.
- A site at a residue outside its frame's translation, other than its end, is refused.
- In the XML input, these are refused (none of them occurs in the XML that NCBI's `rpsblast` and `rpstblastn` write): XML attributes, namespace prefixes, and control characters in character data. NCBI skips attributes, strips prefixes and replaces control characters with a warning.
- An `Hsp_hseq` shorter than its `Hsp_qseq` by more than one letter is refused. NCBI reads past it.
- An alignment that writes a residue other than NUL at or past the end of NCBI's string of its query's frame is refused.
- A query whose strings NCBI holds as 1 GiB or more at once is refused. Whether NCBI's allocation succeeds depends on the host's memory. A nucleotide query of 0 or 1 base gives frames of over a billion residues.
- An `Hsp_qseq` of 2^31 letters or more is refused. NCBI's check of the query's type reads past its end or loops over without end.

</details>

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

<details>
<summary>Rare inputs</summary>

- Some records of that size are handled as NCBI handles them, not refused: the error for a sequence longer than 2,147,483,647 bases, and the record skipped when the length modulo 2^32 is 0, when it has a gi above 4294967295 or when it starts with a gap of unknown length.
- A protein whose `>?` lines sum to a multiple of 2^32 gets NCBI's error.
- With `-input_type blastdb`, these input sequences are refused: one without a definition line, one whose definition lines have bytes that NCBI reads as `#`, a definition line with leaf tax ids and no id, and a first definition line whose `other-info` (where NCBI reads the PIG) is an empty list. NCBI reads past them or writes errors for each read.
- With `-input_type blastdb`, an input database with one masking algorithm is refused when one of its volumes does not list the algorithm and two of its sequences with masks are copied with no sequence with masks of another volume between them (in the input's order, any two; with `-randomize`, in the shuffled order, which depends on the second of the clock). For the second one NCBI reuses the algorithm id it kept from an earlier lookup and reads the masks of that id.
- With `-input_type blastdb`, an input database with masks is refused when the first column of one of its volumes is not the mask column (for example, a volume with another column and no masks). NCBI reads the first column of every volume as its masks. The error names that column's file in the volume.
- With `-input_type blastdb`, an input database that LOSAT's reader refuses (see [Alias files and database arithmetic](#alias-files-and-database-arithmetic) and [blastdbcmd](#blastdbcmd)) is refused.

</details>

## blastdbcmd

| Option or value | What happens |
| :--- | :--- |
| `-outfmt` with `%d` or `%b` (the ASN.1 texts) | Refused. |
| A `-db` list in which an alias and another name give the same volume | Refused. |
| FASTA output of a sequence with a tax id, or `-outfmt` with `%S`, `%L`, `%B`, `%K`, `%C` or `%N`, or `-tax_info` with a name field, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. `-outfmt "%a %T"`, other formats that look nothing up, and `-info` work. |
| A failed write of the `-metadata` output | Refused. NCBI's serial stream reports it. |
| `-range ""`, `-taxids abc`, or a `-taxids` or `-taxidlist` id that is empty, blank, too large or not a number | Refused. NCBI crashes or prints a path. |
| `-mask_sequence_with` below -1 | Refused. NCBI crashes. |
| `-ipg` or `-ipg_batch` on a database without PIG files | Refused. |

<details>
<summary>Rare inputs</summary>

- A name written from a record of a damaged `taxdb.btd` that reaches past NCBI's mapping of the file is refused: the super kingdom of `%K` and of `-tax_info`, and any name of a record whose third tab is past the mapping.
- NCBI reads memory past its mapping there. What it writes, or whether it crashes, depends on its memory. The other names of such a record are written as NCBI writes them, and so are NCBI's messages for the other damaged taxonomy databases.
- A taxonomy lookup that reads a record of an empty `taxdb.btd` is refused, because NCBI maps the file as a null pointer and crashes.
- After a failed write of the retrieved sequences the retrieval goes on, as NCBI's does. A closed pipe ends it at once, with exit status 0 where NCBI is ended by the signal.
- A database date that NCBI does not read for `-metadata`, an empty index file of `-list` with `%t`, and damaged column files are refused. NCBI crashes or prints a path.
- A database whose volume of an alias has masks but not the algorithm is refused.
- `%m` digits that overflow NCBI's `int` are refused.
- LOSAT stops with an error where NCBI crashes or writes uninitialised values. Details are in [the decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

</details>

## dustmasker

`dustmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A title whose text encoding NCBI cannot tell | Refused. NCBI stops with a string error. |
| A database definition line with bytes NCBI reads as `#` | Refused. |
| A database sequence with a tax id, in any `-outfmt`, when `taxdb.btd` is empty | Refused. NCBI maps the file as a null pointer and crashes. |
| A failed write of a `seqloc_*` or `maskinfo_*` output | Refused. NCBI's serial stream reports it; the other formats go on silently, as NCBI's. |
| A FASTA record of 4,294,967,296 letters or more, counting the letters of its `>?` gap lines | Refused. NCBI counts the length of a record in 32 bits, which wrap there. |

<details>
<summary>Rare inputs</summary>

- A title on which NCBI's definition line generator reads past its end is refused.
- A title made from the scientific name of a tax id whose record in a damaged `taxdb.btd` has its third tab past NCBI's mapping of the file is refused.
- A sequence whose title is stored, or whose name lies inside the mapping, is written as NCBI writes it. NCBI's messages for the other damaged taxonomy databases, such as a wrong magic number, are written as NCBI writes them, where it writes them.
- A `Seq-loc` or `Blast-db-mask-info` that LOSAT's serial writer cannot write as NCBI does is refused.

</details>

## windowmasker

`windowmasker` masks with a counts file (`-ustat`) and runs every option of its `-help` for that (`-xmlhelp` on the command line only). The options that only make or convert counts files (`-checkdup`, `-fa_list`, `-mem`, `-unit`, `-genome_size`, `-sformat`, `-smem`) are refused with `-ustat` by NCBI's own argument check, as in NCBI's; `-meta` is accepted and does nothing, as in NCBI's.

| Option or input | What happens |
| :--- | :--- |
| `-mk_counts`, `-convert` | Refused. LOSAT does not make or convert counts files yet. |
| An `-ids` or `-exclude_ids` file matched as text (without `-text_match F`) with an empty id (a line that starts with a space or a tab, `>`, `\|` or `>\|`) | Refused. NCBI writes outside its list of ids and crashes. |
| A counts file whose `##pct:` line holds a word that NCBI reads as a number with the C++ runtime, such as `inf` or `0x10`, with a `-t_*_pct` option | Refused. Plain decimal numbers are read as NCBI reads them. |
| An `ascii` counts file with units of 17 letters or more | Refused. NCBI's result is undefined. |
| A FASTA record or database sequence of 4,294,967,296 letters or more | Refused. NCBI counts positions in 32 bits. |

<details>
<summary>Rare inputs</summary>

- An optimized ascii counts file is refused when the third line, which holds the hash parameters, leaves a value unread that NCBI then checks uninitialised.
- An `obinary` counts file with units of 1 or 2 letters is refused: NCBI reads past an empty array.
- The database and serial-writer refusals of `dustmasker` apply to `windowmasker` too.

</details>

## segmasker

`segmasker` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| Everything that `dustmasker` refuses: titles, database sequences and definition lines, the serial writer, a failed write of a `seqloc_*` or `maskinfo_*` output, and a FASTA record of 4,294,967,296 letters or more | Refused, as for `dustmasker`. |
| A protein database whose residue bytes are above 127 | Refused. `makeblastdb` writes 0 to 27; NCBI's SEG reads larger bytes as negative characters. |

<details>
<summary>Rare inputs</summary>

- Some database errors are reported without the place in NCBI's source that NCBI's report names. These are a missing database or volume, an alias file without `DBLIST`, an empty database name, and an alias value that NCBI cannot convert. LOSAT reports them as NCBI does, apart from that place.

</details>

## convert2blastmask

`convert2blastmask` runs every option of its `-help` (`-xmlhelp` on the command line only).

| Input | What happens |
| :--- | :--- |
| A command-line word that is not UTF-8 | Refused. |
| A `Blast-db-mask-info` that LOSAT's serial writer cannot write as NCBI does | Refused. |
| A failed write of a `maskinfo_*` output | Refused. NCBI's serial stream reports it. `interval` goes on silently, as NCBI's, as for `dustmasker` and `segmasker`. |

<details>
<summary>Rare inputs</summary>

- A byte outside the visible ASCII range in `-masking_options`, or a DEL in a `-parse_seqids` id, is written as NCBI writes it. This is `#` with NCBI's error message for each such byte in `maskinfo_asn1_text` and `maskinfo_asn1_bin`, and for each control character other than TAB, LF and CR in `maskinfo_xml`. It is not a refusal.

</details>

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
- `windowmasker -mk_counts` and `-convert`, which make and convert counts files.
- `psiblast -phi_pattern`.
- The `psiblast` output formats other than `0`, `6`, `7`, `10` and `20` (with `-in_pssm`: `0`, `6`, `7`).
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
