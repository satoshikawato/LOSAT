# Rare and pathological inputs

This companion page collects rare and pathological input refusals separated from [What LOSAT refuses](refusals.md). In almost all of these cases, NCBI BLAST+ 2.17.0 itself crashes, overflows 32-bit integers, enters an endless loop, or reads uninitialised memory outside its buffers.

## All five search programs

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

Database header records (the definition lines of a sequence):

- A header record that NCBI cannot read (damaged, cut short, or with a tag, length or tax id that NCBI's reader does not take) is not refused. Each program writes NCBI's error where NCBI does, after the reports of the query batches before it. A record cut inside its first tag or length ends the run there with exit status 0, as in NCBI. A record that no report reads fails nothing.
- A UTF8String where a VisibleString is expected is read as NCBI reads it, with NCBI's warning `(803.10)` once per run, where NCBI first reads such a record (none in `segmasker` and `dustmasker`). A tag in the long form is read with no message.

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
- Where NCBI's lookup of a reported sequence's id first reads the header record of another sequence that has a UTF8String (a database list built with `-parse_seqids`), the warning `(803.10)` comes where NCBI loads that id: in the batch's prefetch, among the taxonomy database's texts of its other loads, or inside the report where `PRE_FETCH_SEQS_LIMIT` turns the prefetch off.
- Refused there: a prefetch in which NCBI orders two such loads whose warnings name different members, or such a load and one whose taxonomy database writes a text, by the memory addresses of its id records (two ids with one packed number, such as `gb|X00005.1|` and `emb|Y00005.1|`, or two ids that it does not pack), and a prefetch of more than 200 ids, whose order LOSAT does not know. NCBI writes such loads in another order from one run to the next.
- Titles that NCBI makes from taxonomy data are written as NCBI writes them: an empty title from the organism's name, and the end of a protein title ending in `]` rewritten from the first tax id (`alpha, partial [Bacillus sp.]` becomes `alpha [Bacillus sp.]`).
- A version 5 list of databases (`-db "a b"`), some with and some without the accession table of `makeblastdb -parse_seqids`, is refused together with a damaged `taxdb.bti` or `taxdb.btd` that writes text. NCBI's loads of other sequences then write such text or fail on their ids, in all five programs.
- With `-outfmt 5` it is also refused when a tax id is looked up after a batch whose XML report NCBI cuts at a query title or id outside US-ASCII (`blastp`, `blastn`, `tblastn`, `tblastx`).
- NCBI's stop here is `BLAST Database error: DB contains no accession info.` (exit status 2). Whether NCBI writes the texts before the error depends on the order of its loads. LOSAT refuses before the first byte.

Strategy files:

- `-export_search_strategy` is written as NCBI writes it for the five programs, `psiblast` (with a FASTA query or `-in_pssm`), `rpsblast`, `rpstblastn` and `tblastn -in_pssm`, apart from its `ident` (`LOSAT 0.2.0` where NCBI writes `2.17.0+`). For `tblastn -in_pssm` NCBI's request is a PSI-BLAST request that holds the PSSM as read, with `-db`. With `-subject` NCBI gives its error after the prolog.
- `rpsblast` and `rpstblastn` write the file of NCBI's `-num_threads 1` at every `-num_threads` and `-mt_mode`. With `-mt_mode 1` and more than one thread, NCBI's threads write their requests into the file at the same time, and for a database of several volumes they leave out the search space of the batch before.

## blastn

- `-outfmt 1` to `4` with `-parse_deflines` and a query record that is also a subject (one record) is refused in one case, shared with `blastp`. The first of the longest aligned segments of the query's alignments (in their order) is then in an ungapped HSP of the record with itself that is not its first.
- NCBI leaves the record's subject rows out of its merged alignment, writes past its vectors and aborts (`double free or corruption`). LOSAT stops before the first byte of the report. For `blastn` it writes the prolog first, as NCBI does.
- In a run of many queries that LOSAT searches in blocks, it stops before the first byte of the block that holds the query, after the reports of the blocks before.
- A final gapped X-drop of dynamic programming (`-task blastn`, `blastn-short`, `dc-megablast`, `-no_greedy`) is refused when its raw score divided by the gap extension reaches 1,073,741,721. NCBI's arithmetic overflows there.
- A gapped X-drop of greedy extension (`-task megablast`) is refused where NCBI's greedy alignment reads or writes outside its arrays. NCBI then corrupts its heap, reads memory it did not write, crashes, or stops with `BLAST engine error: Unknown error code -1`.
- With the default scores (reward 1, penalty -2), `-xdrop_gap` or `-xdrop_gap_final` is refused from about -7.4 to -1.98e9 bits and from about 1.98e9 to 3.97e9 bits. It is also refused from about -1.98e9 to -3.97e9 bits, unless the other X-drop is as large in NCBI's doubled arithmetic.
- The limits move with the scores. With an even reward, a value past 32-bit integers is refused too.
- An alignment that NCBI keeps outside its query context can reach the report. It comes from an ungapped extension that runs past an end of its context, such as with `-xdrop_ungap 1e10`, which NCBI turns into INT_MIN.
- The same happens with `-ungapped` when the re-evaluation keeps part of an alignment outside the context, and with `-task megablast` when a greedy extension from the middle of such an extension stays outside it.
- No e-value, percent identity, query coverage, `-max_hsps`, `-max_target_seqs`, `-culling_limit` or `-subject_besthit` limit drops such an alignment first. LOSAT drops it where NCBI does.
- Where its query coordinates pass the query's end (a plus-strand alignment right of its context, a minus-strand one left of it), LOSAT shows the alignment as NCBI does in every report. Its query letters outside the context are those of the query's record (with `-query_loc`, also outside the interval), and past the record's end the gap character.
- With `-query_loc`, an ungapped alignment that the interval's start moves into the query's record is shown too, with the record's letters before the interval.
- Where its query coordinates leave the query's record on the left (a start below 1: a start left of the context on the plus strand, an end past it on the minus strand; with `-query_loc` and a gapped search, a start before the interval's first letter), the alignment is refused.
- With `-subject`, NCBI then stops with an exception in the search of the query batch (exit status 255). LOSAT refuses there, after what NCBI writes before: the prologs and the reports of the query batches before.
- From a database, NCBI's reports read such an alignment each their own way (some stop with an exception, others leave it out). LOSAT does not reproduce that yet and stops before the rows of the query batch.
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

## blastp

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
- With `-use_sw_tback` and `-ungapped`, NCBI's traceback still sets up the Smith-Waterman alignment, reads gapped statistics that an ungapped search never computes, and crashes at every query batch it searches, with hits or without. LOSAT refuses at the first query batch that it searches, as `blastx` does. Batches before it whose queries give no Karlin-Altschul parameters are reported first, and a run without a searched batch ends as NCBI's does.
- With `-ungapped -window_size 0` and a low `-threshold`, `-outfmt 1` to `4` is refused for a one-hit ungapped alignment of length 0 among the alignments shown, where NCBI's merge of the query-anchored alignments frees memory that it still reads. NCBI crashed in every such run we made but one; LOSAT writes that one. Every other format writes NCBI's report of such an alignment (apart from the version line, which names LOSAT).
- With `-ungapped` and a database, the 5,000,000-residue limit is checked against the longest sequence of the database, also when a sequence list or an alias file leaves that sequence out of the search.
- `-db_hard_mask` searches a database whose masks overlap, nest or are stored out of order (`segmasker` writes such masks) as NCBI searches it, with the rows NCBI gives. NCBI cuts each sequence into parts at its masks in the order stored. Where masks overlap, a part ends before it starts, and NCBI's table of diagonals moves back by it. A part shorter than the word gets one word from NCBI's two-hit scan, read from the stored letters at its start and past it, and LOSAT gives the rows NCBI makes from that word. Where that word runs past the end of the sequence into the next sequence's letters, it finds nothing: no query word can hit a word that holds the end of a sequence, since the FASTA reader drops the gap letter `-` from a query (with NCBI's warning).
- A sequence with overlapping masks whose word hits, after those masks, reach entries of NCBI's table of diagonals that belong to the sequences searched before it, or leave entries that the next sequence reaches, is searched on its own: LOSAT moves the table past those entries, so the sequence's rows do not depend on the other sequences' entries or on `-num_threads`. Like NCBI's 1-thread rows, they can depend on the lengths of the sequences searched before it: where NCBI's 1-thread search places the sequence so that its word hits reach the part of the table that no sequence has written yet, LOSAT starts it on such a table at that place. NCBI keeps one table per thread. With one thread its rows for such a sequence depend on the neighbouring sequences, and with several threads they can change from run to run (in a search of 200 copies of one pair of sequences with 4 threads, three runs found the second sequence's alignment in 5, 6 and 3 copies, where one thread finds it in none). See [Differences](differences.md#-db_hard_mask-with-overlapping-masks-that-reach-a-neighbouring-sequence).
- With NCBI's defaults for `segmasker`, A. thaliana TAIR10 has 118 such proteins among the 35,642 that `segmasker` masks at the default window, and 578 with `-window_size 0`; the first 60,000 sequences of UniRef50 have 242 among their 59,484 masked ones, and 1,438 with `-window_size 0`. `-db_hard_mask` searches both databases.
- The cases below are refused. One refused sequence refuses the whole database.
- A part between overlapping masks shorter than the word whose seed NCBI's gapped search never ends on: the part's length, negative after overlapping masks, sets the bounds of the tree that NCBI keeps the seed's alignments in, and the tree splits its nodes without end while its memory grows. Also refused, until LOSAT does that merge: such a part whose alignments meet another part's in NCBI's merge across their overlap. LOSAT stops with exit status 1 when the search reaches the sequence, and the error names it by its position (OID) in the database.
- Masks stored out of order whose part after the later-stored mask starts before the end of a part before the earlier one, both parts able to hold a word: NCBI joins the two parts' alignments by a merge across that overlap, which LOSAT does not do yet. Masks out of order whose parts do not meet so are searched.
- A mask that runs to the end after the part from the start, with masks stored after it: NCBI searches that part past its end without the later masks and applies them only to the letters of its alignments. Without later masks such a sequence is searched over that part (see [Differences](differences.md)).
- A mask that ends before it starts, or ends past its sequence.
- The refusals other than those found during the search come before the search, so before the first byte of output, and name the first such sequence by its position in the database.
- `tblastn` and `tblastx` search masks that overlap or nest in the order of their start as NCBI BLAST+ does, and refuse the cases in their sections.

## blastx

- With `-ungapped -window_size 0`, a low `-threshold` and `-subject` (also with `-subject_loc`), LOSAT refuses in every format but `-outfmt 18` a one-hit ungapped alignment of length 0. It sits at the first letter of the searched part of a subject that has alignments of letters among the hits that the reports show.
- NCBI may stop earlier at another alignment of length 0. LOSAT refuses the batch before its first byte.
- NCBI's Seq-loc mapper cuts such an alignment into two alignments without scores. NCBI then writes rows of the query's id with an E-value of 0.0 and a bit score of -1.0, stops with an exception, or runs out of memory, depending on the frame, the fields and the memory limit.
- The same refusal applies to a subject whose every such alignment has length 0, one of them at the first or second letter of the searched part. The mapper, finding no sequence type, then makes alignments without scores of query nucleotides against the subject's letters.
- NCBI then writes more rows, throws exceptions or runs out of memory when there are two or more. The tabular reports of a subject's only one are written as NCBI writes them.
- Elsewhere LOSAT reports alignments of length 0 as NCBI does.
- With `-ungapped -window_size 0` and a low `-threshold`, `-outfmt 1` to `4` is refused for a one-hit ungapped alignment of length 0 among the alignments shown, where NCBI's merge of the query-anchored alignments frees memory that it still reads. NCBI crashed in every such run we made.
- One case is an alignment whose query row starts at the first nucleotide of a codon of an alignment in the query's reading frame +1. The other is an alignment followed on its query row by unaligned letters, whose segment the merge frees unless it is the last segment it made.
- LOSAT stops before the first byte of the report, unless the report of a query before it ends with NCBI's exception.
- With `-use_sw_tback`, a query and subject whose Smith-Waterman alignment NCBI never finishes is refused. Only a large `-evalue` or a small search space (`-searchsp`, `-dbsize`, or a small query and subject) can reach it.
- With `-use_sw_tback`, a query batch on which NCBI crashes is refused: one whose last Smith-Waterman alignment fails after the batch's first query has kept one. Where the first query has kept none, the batch ends with NCBI's `BLAST engine error: Traceback failed`, and with LOSAT's `LOSAT engine error: Traceback failed`.
- For these two, LOSAT writes the results of the query batches before and then stops with the error.

## tblastn

- A subject (or `-subject_loc` range) of fewer than 12 bases whose alignment reaches the traceback (with `-ungapped`, the ungapped re-evaluation) makes NCBI read the translation of the subject's frame, which it never makes for a subject this short. NCBI then crashes. With `-ungapped`, a one-hit alignment of length 0 or less (`-window_size 0`) is deleted there without that read, so a subject whose only alignments are such ones is searched as NCBI searches it.
- In practice this happens with `-window_size 0`: no frame of such a subject, three letters at most, holds two words that do not overlap.
- LOSAT writes the results of the query batches before and then stops with the error.
- With `-window_size 0` (the one-hit word finder) and a low `-threshold`, LOSAT refuses the cases below.
- Two gapped starts of one-hit alignments of negative length that lie at the end of their 5,000,000-letter part of a frame, or past it. NCBI never finishes on them.
- With `-in_pssm`, a one-hit alignment whose window of 11 letters passes the PSSM's last position. NCBI reads past its matrix.
- With `-in_pssm` and a PSSM whose `-` column scores above NCBI's minimum, a one-hit alignment whose gapped start lies past the end of NCBI's translation of the subject (after its last frame). NCBI's gapped alignment reads the byte after its end, which it never sets, or past it.
- A one-hit alignment whose gapped start lies after the end of its frame's translation and that reaches NCBI's traceback. NCBI reads a byte its translation does not set there.
- The report of a one-hit alignment whose traceback starts on the end marker. NCBI places its column past the subject's end, and LOSAT does not reproduce that result yet. Such alignments that NCBI drops before its reports give NCBI's output (apart from the version line, which names LOSAT).
- For each of these `-window_size 0` refusals but the last, LOSAT writes the results of the query batches before (with `-in_pssm`, one batch: the report's header) and then stops with the error.
- With `-db_hard_mask`, where a sequence's masks leave a reading frame an empty part at its start (a mask within the sequence's first two bases, or within the one or two bases after its last whole codon, for example [1200, 1202) of 1,202 bases, or an empty mask at its start), NCBI searches the next part as the whole frame, then frees memory that it did not allocate and crashes in every run. LOSAT searches the frame as outside the mask; see [Differences](differences.md).
- Where a frame's first part, from the frame's first letter, ends at a mask with no other part of the frame after it (masks at or near either end of a sequence), LOSAT searches the part alone; see [Differences](differences.md).
- With `-db_hard_mask` on a database whose masks overlap or nest in the order of their start, `tblastn` searches the database as NCBI BLAST+ does, searching any sequence whose hits would reach another sequence's diagonal entries on its own as `blastp` and `blastx` do (see [Differences](differences.md#-db_hard_mask-with-overlapping-masks-that-reach-a-neighbouring-sequence)), with the same rows at every `-num_threads`, and refuses the cases below (one refused sequence refuses the whole database, exit status 1):
  - Masks stored out of the order of their start, or a mask that ends before it starts or past its sequence (refused before the search).
  - A sequence with overlapping masks one of whose reading frames NCBI searches whole past the later masks (a part from the frame's start ending at a mask with no later part after it, or an empty first part before overlapping masks; refused before the search).
  - A part between overlapping masks shorter than the word whose seed NCBI's gapped search never ends on (its negative length sets the bounds of NCBI's interval tree, which splits nodes without end while memory grows), or whose alignments meet another part's in NCBI's merge across an overlap read past the list of ranges (refused during the search, naming the sequence by OID). Also refused before the search for a query whose lookup holds residue 0 (no FASTA query does): a short part whose one word reads past a reading frame's end.
- With `-use_sw_tback`, a query and subject whose Smith-Waterman alignment NCBI never finishes is refused. Only a large `-evalue` or a small search space (`-searchsp`, `-dbsize`, or a small query and subject) can reach it.
- With `-use_sw_tback`, a query batch on which NCBI crashes is refused: one whose last Smith-Waterman alignment fails after the batch's first query has kept one. Where the first query has kept none, the batch ends with NCBI's `BLAST engine error: Traceback failed`, and with LOSAT's `LOSAT engine error: Traceback failed`.
- With `-use_sw_tback` and `-ungapped`, NCBI's traceback still sets up the Smith-Waterman alignment, reads gapped statistics that an ungapped search never computes, and crashes at every query batch it searches, with hits or without. LOSAT refuses at the first query batch that it searches, as `blastx` does.
- For these `tblastn` cases, LOSAT writes the results of the query batches before and then stops with the error.

## tblastx

- With `-window_size 0` and a low `-threshold`, a one-hit ungapped alignment of length 0 or less is refused when it meets another in the overlap of two 5,000,000-letter parts of a subject frame and their lengths add up to 0. NCBI merges them by dividing by that sum, and its result is undefined.
- An alignment on a subject (or its `-subject_loc` range) of fewer than 12 bases that reaches the ungapped re-evaluation makes NCBI read the translation of the subject's frame and crash. A one-hit alignment of length 0 or less (`-window_size 0`) is deleted there without that read.
- `tblastx` first writes in full the reports of the query batches before. NCBI writes them before it crashes, except what its output buffer still holds (with `-outfmt 0`, the end of the last one).
- An ungapped extension that passes the end of a query frame into the next frame of the query batch and scores higher there is refused, with an `-xdrop_ungap` of about 15,000 bits or more (BLOSUM62; a raw X-drop above 32,768).
- In the cases checked, NCBI keeps such an alignment in the frame where it started, with query coordinates outside that frame. It writes coordinates past the query or leaves the alignment out of every report.
- LOSAT writes the results of the query batches before and then stops with the error.
- An extension that passes into the next frame without scoring higher there gives NCBI's output (apart from the version line, which names LOSAT).
- `-db_hard_mask` on the sequences described under `tblastn` is searched as in `tblastn`: a frame after an empty part at its start whole (its masked bases stay masked), and a part that ends before its frame's end alone; see [Differences](differences.md).
- With `-db_hard_mask` on a database whose masks overlap or nest in the order of their start, `tblastx` searches the database as NCBI BLAST+ does, searching any sequence whose hits would reach another sequence's diagonal entries on its own as `blastp` and `blastx` do (see [Differences](differences.md#-db_hard_mask-with-overlapping-masks-that-reach-a-neighbouring-sequence)), with the same rows at every `-num_threads`, and refuses the same shapes as `tblastn` before the search (masks stored out of start order or ending before they start or past the sequence, and a reading frame searched whole past later overlapping masks) as well as, during the search, a short part whose alignments meet another part's in NCBI's merge across an overlap read past its list of ranges. (`tblastx` is ungapped, so NCBI's endless gapped interval-tree split does not occur.)

## Alias files and database arithmetic

- `blastn` and `tblastx` refuse an alias file's `NSEQ 0`, and its `LENGTH` or `STATS_TOTLEN` without `NSEQ` or `STATS_NSEQ` over volumes without sequences. NCBI crashes on them where a query is searched.
- `blastp`, `blastx`, `tblastn` and `psiblast` refuse the same cases with `-searchsp` or `OLD_FSC`. `blastx` also refuses them with `-ungapped`.
- A run of such a database whose queries have no letters to search is reported as in NCBI by every program.
- `rpsblast` and `rpstblastn`: an alias file's `LENGTH` and `NSEQ` over one volume also set the length that NCBI's search covers, their sum, which can be shorter or longer than the profiles. Refused: a sum above 2,147,483,647 (NCBI frees the profiles twice and aborts), or one so long that NCBI's diagonal table passes 2^30 entries (its loop never ends).
- A sum shorter than the profiles is refused where two gapped alignments of a query (of a frame for `rpstblastn`) start past it plus one. NCBI's interval tree then splits without end. With one such alignment the search runs as in NCBI.
- A sum longer than the profiles is refused where an ungapped extension runs past the row after the last profile. This needs an X-drop above the profiles' delimiter scores and a sum of the profiles plus 2 or more (NCBI reads past its score matrix). At plus 1 the search runs as in NCBI where that row lies within the `.rps` file's last memory page, and is refused where it lies past it (NCBI reads past the file's mapping).
- `rpsblast` and `rpstblastn` with `-searchsp` or `OLD_FSC` refuse an `NSEQ 0` over one volume where a query is searched (NCBI crashes). A run of queries without a letter to search is reported as in NCBI.

## psiblast

A PSSM given with `-in_pssm` is refused when its input or content is one of these. In each case NCBI crashes, does not finish, reads out of its buffers, goes on with a wrong value, or depends on the machine.

- A PSSM of frequency ratios only, written by row. NCBI reads the scores it computes for it in the wrong order.
- A PSSM whose position-based scaling gives score frequencies but no Karlin-Altschul parameters. NCBI goes on with a Lambda of -1.
- A PSSM with a negative ungapped Karlin-Altschul value (`lambdaUngapped`, `kappaUngapped` or `hUngapped`) or an ungapped value that is not a number gets NCBI's error of a query without Karlin-Altschul parameters, as in NCBI. Gapped values that are negative, infinite or not a number, and `lambdaUngapped` of `PLUS-INFINITY`, in ASN.1 text or binary, are searched as in NCBI: the hits, E-values (`0.0`), bit scores (`inf`, `-1.0`) and report footers are NCBI's. A gapped `lambda` below 0 (a finite one), or one so small that NCBI's search for the score of the hit-saving cutoff passes 32-bit integers, gives a report without hits, where NCBI never ends or writes the same report (see [Differences from NCBI BLAST+](differences.md#-in_pssm-with-a-very-small-or-negative-gapped-lambda)).
- A PSSM of frequency ratios only from which NCBI's PSSM engine makes no scaled scores. NCBI does not finish or stops with an error.
- A query without a Seq-id, with a Geninfo import id (`gim|`, which NCBI cannot resolve), or with a database ordinal (`gnl|BL_ORD_ID|n`).
- A query with a local Seq-id among other Seq-ids (NCBI's reports then show it as no FASTA query is shown), with an empty local id, or with a Seq-id that has neither an accession nor a name.
- For the last kind, NCBI's run depends on the database and the report: a database error with some databases, `Failed to resolve SeqId` in the pairwise report, `Unknown` in tabular reports.
- A query with more than one title, or with white space at the ends of its title.
- A score outside -32,768 to 32,767.
- A PSSM whose columns' largest scores add up to 2^30 or more. Every PSSM of fewer than 19,800 columns passes.
- A query whose letters give no Karlin-Altschul parameters. NCBI crashes.
- A PSSM query whose Seq-inst `length` is not its number of letters, except a `length` that is one less than the letters: that query is searched as its first `length` letters, as in NCBI. With other lengths (two or more fewer, more, 0 to 10, 255, 2^32 - 1) NCBI writes outside its query buffer or crashes. A PSSM of frequency ratios only whose `length` is one less than its letters or more, below 1 GiB, gets NCBI's PSSM engine error (`Number of columns returned by IPssmInputFreqRatiosFreqRatios does not match query length`, exit status 3); with two or more letters past its `length` NCBI writes outside its query buffer, and from 1 GiB on whether NCBI allocates the buffer depends on the host's memory: both are refused.
- Query Seq-data other than `iupacaa`, `ncbieaa` and `ncbistdaa`. (A PSSM without a query, or with `query NULL`, ends with NCBI's exception of its empty Seq-entry, as in NCBI.)
- Two of NCBI's pointer forms, in the cases where NCBI's result is not decided by the PSSM alone. One is a value stated as `:Type` of another type than its place declares, where what NCBI does depends on the types its process has met before: in the author names that NCBI skips, and where reading the value would post a message, meet another such value or stop with an error that NCBI does not report as a PSSM of an unrecognized format (a number that is not one). Elsewhere such a value (`:Object-id` where `params` belongs) ends with NCBI's `Unrecognized format for PSSM`, as in NCBI. The other is an element stated as `NULL` (`05 00` in ASN.1 binary) of a list that NCBI's object manager dereferences when the query joins its scope, which ends with a report of NCBI's memory addresses or a crash: a Seq-id of the query, a member of a Bioseq-set or a Seq-id of a member, an annotation, a feature, an alignment or a graph, an annotation's id or description, the parts of a location, a feature's ids or cross-references, a Seq-table column and a list inside an alignment. An element stated as `NULL` in any other list is read as NCBI reads it: in a descriptor of the query it ends with NCBI's engine error, the lists that NCBI drops or does not read (such as a variation's population data or a feature's database cross-references) are searched, and `-out_pssm`, `-export_search_strategy` and the archive of `-outfmt 11` stop at the element as NCBI's writer does. `NULL` for an optional member and `:Type` of the declared type are read as NCBI reads them, and `@n` ends with NCBI's error, in ASN.1 text and binary.
- A query without letters, or whose number of letters is not the `length` of its Seq-inst, other than the lengths above. For these NCBI crashes or writes past its buffer (seen with 0, 1, 2, 10, 28, 31, 60, 255 and 4,294,967,295 for 30 letters of a PSSM with scores).
- A PSSM in ASN.1 binary form that ends inside a string or a byte block (OCTET STRING or BIT STRING) whose length asks for 1 GiB or more and less than 2^47 bytes. Whether NCBI's allocation succeeds depends on the host's memory; at 2^47 bytes or more, NCBI's `BLAST ran out of memory` (exit status 4) is written, as in NCBI.
- A query whose annotations (`annot`) give a case of NCBI's object manager whose output names a memory address or depends on the machine. These cases are an empty location (`Empty region in 0x...`), a NULL-pointer report, a Dense-seg that NCBI divides by zero, malformed BitMagic data of a Seq-table, reads past a Packed-seg's vectors, and integer overflows. BitMagic data is decoded as NCBI decodes it; refused are the data that NCBI reads past (such as empty or cut data), decodes with undefined behaviour (such as an XOR block or a write past BitMagic's blocks) or from memory it did not write, and a GAP block that is not one where NCBI's scalar, SSE2 and SSE4.2 builds decode it or answer a row of it differently. In each refusal of the query's entry (these cases, the next two and the elements stated as `NULL` above), psiblast first writes what NCBI writes before its object manager stops: the start of the report of `-outfmt 0` and the header line of `-outfmt 20` (tblastn writes nothing before it).
- Also two Seq-ids of one location (a feature's, a Seq-table's or a Seq-table row's) with different messages that NCBI orders by address; a Seq-table of more than 10,000,000 rows; a Std-seg of more than 1,000,000 rows; and XML and JSON reports of an untitled query whose annotations hold features.
- An error about the query's annotations, or a Seq-id twice in the entry. Two Seq-ids count as one when NCBI's object manager takes them as one: equal ones, and ones that differ only where NCBI does not look, such as the case of letters, GenBank, EMBL or DDBJ for an accession alone, or PDB's `chain` and `chain-id`.
- A query Seq-entry set whose annotations NCBI's object manager takes, and which NCBI then stops at with its exception (`CSeq_entry::GetSeq()`).
- These last two are not refused where NCBI stops before the search and writes its prolog and error: `-export_search_strategy` with `-subject`, a PSSM that NCBI's checks reject, or a PSSM of frequency ratios only whose scores NCBI's PSSM engine cannot make.
- With `-subject`, a Seq-id twice in the query's entry where NCBI's exception may write a Seq-id as a subject record spelled it. That is a GenBank, EMBL or DDBJ accession alone, a name, a patent number or country with a letter, or more than 64 letters in a local or general id or in an accession.
- With `-out_pssm`, and with `-outfmt 11` and a database (for each `ARCHIVE_FORMAT` that writes: no variable, `xml` or `bin`), a PSSM of frequency ratios only without `rpsdbparams`. NCBI cannot write its matrix name and stops with a C++ exception while it writes the file or the archive's request. With `-subject`, with `-outfmt 6` and with `-export_search_strategy` NCBI writes no such object, and LOSAT answers as NCBI does.
- With `-out_ascii_pssm`, a PSSM whose information content comes without the other columns of NCBI's long ASCII form.
- With `-out_ascii_pssm`, a PSSM query with a residue that its coding does not have. NCBI prints its letter from memory it never set.
- With `-out_pssm`, a PSSM in ASN.1 binary or XML form whose author names (`Auth-list` `names`, which NCBI reads again only when it writes them in ASN.1 text) give an error when they are read again, such as a number too large for its type: NCBI writes part of the file and stops with an exception. A byte outside the visible ASCII range there is not one of them: LOSAT gives NCBI's message for it and writes `#`, as NCBI does. Nor is the warning of a UTF8String tag: NCBI gives it once, when it first reads the file, and writes the names.
- A real number in a PSSM whose exponent overflows NCBI's 32-bit `int` (`1.55e-2147483647`) is read with the exact exponent, the smallest normal double; see [Differences](differences.md). The PSSM's ASN.1 text, ASN.1 binary and XML are read as NCBI reads them.

## tblastn with `-in_pssm`

- Where NCBI stops before the search, as for `psiblast -in_pssm`, the output formats (every format `tblastn` writes, `-outfmt 11` included) and the query's titles are not refused. NCBI's prolog and error are written.
- NCBI's object manager can end the run as it adds the query's Seq-entry to its scope. LOSAT does not refuse what NCBI's `tblastn` makes after that point either: its formatter, its requests and its search.
- Examples of those later limits: an `-xdrop_gap` or `-xdrop_gap_final` that NCBI's kernels cannot size; `-dbsize` with `-outfmt 11` and a list that keeps no sequence; a database named `n/a` or `SRA` in a Seq-annot format; database lengths whose sums the report cannot write.
- NCBI's messages and error are written there, with `-outfmt 11` in its error archive.

## rpsblast

- An X-drop whose raw score, or its product with the database's scale, leaves 32-bit integers is refused. NCBI's conversion is undefined.
- A gapped X-drop whose raw score over the scaled gap extension reaches 1,073,741,721 is refused (about 4.1e8 bits, at scales 1 and 2 only). The sizes of NCBI's alignment kernels then overflow.
- A `.freq` file whose rows give a profile no Karlin-Altschul parameters is refused. NCBI goes on with a Lambda of -1.
- A profile database whose profile can score above 2,147,483,647 (the sum of its rows' best scores) is refused. NCBI's 32-bit scores overflow.
- A profile whose Smith-Waterman redo (`-use_sw_tback`) NCBI repeats without end is refused. This needs three conditions. The e-value of a score of 0 (the query's search space times K) passes `-evalue`. The profile has two or more preliminary HSPs. Its first position scores the query's first letter above 0. NCBI stops only when an allocation fails.
- A profile database whose matrix and gap costs are not BLOSUM62 11/1, or whose scale is not a whole number, is refused when a scaled gap cost, X-drop or score cutoff then leaves 32-bit integers (NCBI's products overflow).
- It is also refused when NCBI's Lambda changes from one e-value to the next at that scale (such as 531, where NCBI divides Lambda by the scale and multiplies it back).

## rpsbproc

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

## makeblastdb

- Some records of that size are handled as NCBI handles them, not refused: the error for a sequence longer than 2,147,483,647 bases, and the record skipped when the length modulo 2^32 is 0, when it has a gi above 4294967295 or when it starts with a gap of unknown length.
- A protein whose `>?` lines sum to a multiple of 2^32 gets NCBI's error.
- With `-input_type blastdb`, these input sequences are refused: one without a definition line, one whose definition lines have bytes that NCBI reads as `#`, a definition line with leaf tax ids and no id, and a first definition line whose `other-info` (where NCBI reads the PIG) is an empty list. NCBI reads past them or writes errors for each read.
- With `-input_type blastdb`, an input database with one masking algorithm is refused when one of its volumes does not list the algorithm and two of its sequences with masks are copied with no sequence with masks of another volume between them (in the input's order, any two; with `-randomize`, in the shuffled order, which depends on the second of the clock). For the second one NCBI reuses the algorithm id it kept from an earlier lookup and reads the masks of that id.
- With `-input_type blastdb`, an input database with masks is refused when the first column of one of its volumes is not the mask column (for example, a volume with another column and no masks). NCBI reads the first column of every volume as its masks. The error names that column's file in the volume.
- With `-input_type blastdb`, an input database that LOSAT's reader refuses (see [Alias files and database arithmetic](#alias-files-and-database-arithmetic) and [blastdbcmd](#blastdbcmd)) is refused.

## blastdbcmd

- A name written from a record of a damaged `taxdb.btd` that reaches past NCBI's mapping of the file is refused: the super kingdom of `%K` and of `-tax_info`, and any name of a record whose third tab is past the mapping.
- NCBI reads memory past its mapping there. What it writes, or whether it crashes, depends on its memory. The other names of such a record are written as NCBI writes them, and so are NCBI's messages for the other damaged taxonomy databases.
- A taxonomy lookup that reads a record of an empty `taxdb.btd` is refused, because NCBI maps the file as a null pointer and crashes.
- After a failed write of the retrieved sequences the retrieval goes on, as NCBI's does. A closed pipe ends it at once, with exit status 0 where NCBI is ended by the signal.
- A database date that NCBI does not read for `-metadata`, an empty index file of `-list` with `%t`, and damaged column files are refused. NCBI crashes or prints a path.
- A database whose volume of an alias has masks but not the algorithm is refused.
- `%m` digits that overflow NCBI's `int` are refused.
- LOSAT stops with an error where NCBI crashes or writes uninitialised values. Details are in [the decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

## dustmasker

- A title on which NCBI's definition line generator reads past its end is refused.
- A title made from the scientific name of a tax id whose record in a damaged `taxdb.btd` has its third tab past NCBI's mapping of the file is refused.
- A sequence whose title is stored, or whose name lies inside the mapping, is written as NCBI writes it. NCBI's messages for the other damaged taxonomy databases, such as a wrong magic number, are written as NCBI writes them, where it writes them.
- A `Seq-loc` or `Blast-db-mask-info` that LOSAT's serial writer cannot write as NCBI does is refused.

## windowmasker

- An optimized ascii counts file is refused when the third line, which holds the hash parameters, leaves a value unread that NCBI then checks uninitialised.
- An `obinary` counts file with units of 1 or 2 letters is refused: NCBI reads past an empty array.
- `-mk_counts` is refused when its counters cannot be allocated on the computer (a large `-mem` with long units), and an `oascii` or `obinary` output when its cache bit array cannot: NCBI's result then depends on the computer's memory.
- `-convert` of an `ascii` counts file with units of 17 letters or more is refused. NCBI converts it without an error, but its reverse complement of each unit shifts the 32-bit unit past its width, which is undefined.
- `-convert` of an `ascii` counts file whose unit size line is 0 (or a multiple of 256, which NCBI reads as 0) to `oascii` or `obinary` is refused. NCBI's search for the size of its hash table never ends. The `ascii` and `binary` outputs of such a file are written as NCBI writes them.
- A `binary` or `obinary` counts file whose `##pct:` line names 4294967295 or more as its largest count is refused. NCBI then reads past its table of percentages.
- The database and serial-writer refusals of `dustmasker` apply to `windowmasker` too.

## segmasker

- Some database errors are reported without the place in NCBI's source that NCBI's report names. These are a missing database or volume, an alias file without `DBLIST`, an empty database name, and an alias value that NCBI cannot convert. LOSAT reports them as NCBI does, apart from that place.

## convert2blastmask

- A byte outside the visible ASCII range in `-masking_options`, or a DEL in a `-parse_seqids` id, is written as NCBI writes it. This is `#` with NCBI's error message for each such byte in `maskinfo_asn1_text` and `maskinfo_asn1_bin`, and for each control character other than TAB, LF and CR in `maskinfo_xml`. It is not a refusal.
