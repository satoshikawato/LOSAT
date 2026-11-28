# Product Decision: database search (`-db`)

- Decision ID: `PD-LOSAT-DATABASE-SEARCH`
- Version: 1.3
- Date: 2026-11-02 (version 1.0: 2026-10-12; 1.1: 2026-10-13; 1.2: 2026-11-01, several databases
  in one `-db`, P4, line DI-8f, decided under the maintainer's standing delegation; its PR asks
  the maintainer to accept it; 1.3: databases with parsed Seq-ids, line DBID, on the
  maintainer's request of 2026-10-21, with the refusals under the same delegation, which the
  maintainer accepted on 2026-11-02; the same day line DBID N2a refused a stored
  `gnl|BL_ORD_ID|n` of another number, which version 1.3 listed as a known difference, and
  line DBID N2b left every kind of Seq-align id to the lookup, which ported refusals that
  version 1.3 listed as wider than NCBI's failure; on 2026-11-03 line DBID N2c refused the weak
  matches of NCBI's scope, one of which version 1.3 listed as a known difference, and line DBID
  N1b compared the same record with other tax ids as the reports show it, which version 1.3 listed
  as a refusal wider than NCBI's failure, and line DBID N3 showed copies of one record that follow
  each other in one hit, as NCBI shows them, which version 1.3 refused; on 2026-11-04 line DBID
  N1b L1 refused such copies for a title made from the tax ids only where a report shows that
  title and an empty title only where the taxonomy database names its tax id, which version 1.3
  listed as refusals wider than NCBI's failure, and refused the tabular id that NCBI takes from
  a nucleotide title made from the tax ids, and line DBID H3 refused the outfmt 0 of an empty
  title whose tax id the taxonomy database names in a database without parsed Seq-ids, which
  LOSAT wrote without the organism, and line DBID P1 wrote psiblast's messages of the taxonomy
  database's prefetch, which version 1.3 listed as a known difference, and line DBID P2 wrote
  these messages in BLASTN, TBLASTN, BLASTX and TBLASTX where `PRE_FETCH_SEQS_LIMIT` turns the
  prefetch off, after each query's warnings as NCBI's formatters fetch, which line DBID P1 listed
  as a known difference, and line DBID P3 wrote them in outfmt 8, 9 and 12 with
  `PRE_FETCH_SEQS_LIMIT=2147483647`, at which NCBI prefetches for these formats too, which version
  1.3 listed as a known difference, and refused the subjects with parsed Seq-ids that this
  prefetch cannot read as NCBI does; on 2026-11-05 line DBID L1 looked a subject up for a
  tabular report whose only field of the alignment block is `nident` only where one of its shown
  HSPs has identities, as NCBI's formatter reads it with the prefetch off, which line DBID P2
  listed as a known difference; on 2026-11-07 line DBID pfbatch refused what the prefetch reads
  for a report that reads nothing of the subject only in a batch whose Seq-aligns do not
  outnumber a positive `PRE_FETCH_SEQS_LIMIT`, in every program, which line DBID tabfields listed
  as a known difference; on 2026-11-08 line DBID pfwarn wrote NCBI's `Error pre-fetching sequence
  data` where a batch's prefetch fails at a sequence that no report reads, which lines DBID P3 and
  pfbatch refused, and took the strand of BLASTN's ungapped Seq-align of a sequence from its first
  HSP, which line DBID pfbatch listed as a known difference, and line DBID weakbatch refused a
  weak match with another sequence of the database only in a batch whose prefetch
  `PRE_FETCH_SEQS_LIMIT` turns off, which line DBID N2c listed as a known difference, and in psiblast
  checked the sequences that each round's PSSM reads, and refused following hits that its PSSM
  joins, where LOSAT wrote other later rounds than NCBI, and line GZ read a `-in` file named
  `.gz` as makeblastdb's decompressor reads it, where a gzip trailer whose CRC-32 or length is
  wrong, or that is cut short, stopped the build, and line DBID pssmrows made psiblast's PSSM
  count its rows by the Seq-align ids' handles and join following hits whose ids `Match`, as NCBI
  does, which lifted psiblast's refusal of the hits that its PSSM joins and showed copies of one
  record as NCBI's PSSM does, and refused the ids that NCBI's handles can take for one and copies
  with another hit between them; the same day `losat makeblastdb` wrote makeblastdb's exception
  reports (`Error: NCBI C++ Exception:`) for the errors of its general handler and read
  `-max_file_sz` where makeblastdb reads it, after the head of the log, which version 1.0 listed as
  known differences; on 2026-11-11 line SW (SW8 follow-up 2) searched a volume with column files,
  the masks of `makeblastdb -mask_data` or a user column, damaged ones included, without reading
  them, as NCBI's searches do without a database mask named by `-db_soft_mask`, `-db_hard_mask` or
  an imported search strategy, all three of which LOSAT still refuses; the reader of stage DB-1
  refused such a volume; the same day line DBID pssmread made psiblast's PSSM read each alignment's
  letters by its Seq-align id as NCBI's scope gets them, another sequence's or none with NCBI's
  warning, which line DBID pssmrows refused; on 2026-11-14 line DBID ordname wrote the tabular
  subject id of a sequence without parsed Seq-ids whose title gives no word (a title of spaces, a
  nucleotide's empty title) as NCBI does, the general id `gnl|BL_ORD_ID|n` and `BL_ORD_ID:n` with
  its OID, which version 1.0 refused, and refused it where the taxonomy database names the tax id
  of a nucleotide's empty title, from which NCBI makes the id, as line DBID H3 refused outfmt 0,
  only for a table that prints a subject id, which version 1.0 refused for every table; the same
  day line RG item T5 made an empty title of a sequence whose tax id the taxonomy database names
  as NCBI makes it from the source of its Bioseq (`unnamed protein product [Homo sapiens]`,
  `Homo sapiens`, the first letter of a nucleotide's name upper-cased, a TPA or MAG prefix dropped
  in the heading) in outfmt 0, 18 and the tabular subject ids (`Homo` for a nucleotide's), in
  databases with and without parsed Seq-ids, which lines DBID H3, N1b L1 and ordname refused
  (also of a RefSeq chromosome, mRNA or ncRNA sequence, whose own NCBI functions make nothing of
  a database Bioseq), and kept refused one whose tax id's lookup reads past the end of
  `taxdb.btd` (NCBI reads past its mapping of the file; whether it crashes depends on its memory
  map); the same
  day line RG item T3a refused in BLASTP's outfmt 1 to 4 a reported sequence whose Seq-align id is
  the id `lcl|Query_N` of a query of its batch read without `-parse_deflines`, and also its
  display id in outfmt 1 to 4 (outfmt 0 stays line RG item X7's), and in outfmt 0 to 4 of
  BLASTP, psiblast and rpsblast a sequence of the description table whose display id is the id
  of a query read with `-parse_deflines`, where that table looks the title up again by the
  display id, which LOSAT wrote as the sequence's own; on 2026-11-15 line RG item T3b refused the
  same in BLASTN's and TBLASTN's outfmt 1 to 4, which it opened, and line RG item T3c in BLASTX's
  and TBLASTX's outfmt 1 to 4, which it opened; on 2026-11-14 line RG item T5c
  wrote the end of a protein title ending in `]` that NCBI rewrites from the first tax id
  (`x_AdjustProteinTitleSuffix`, with or without a taxonomy database) in the outfmt 0 heading,
  outfmt 18 and the tabular subject ids, which version 1.3 refused in databases with parsed
  Seq-ids and LOSAT wrote unchanged in the tabular ids of a database without them, and where the
  rewritten title has no first word (`MAG: (plastid) [ ]` cut to `MAG:` from the tax id) the
  general id of line DBID ordname, in every program; on 2026-11-15 line RG item T5b refused a
  lookup of the taxonomy database that reads past NCBI's mapping of `taxdb.btd` only where NCBI
  faults for certain (a record of an empty `taxdb.btd`, which NCBI maps as a null pointer) and
  where a report writes a name whose bytes come from past the mapping, and wrote the others as NCBI
  writes them where it completes, which lines DBID and RG item T5 refused at every lookup
  (approved exception 14 of PD-LOSAT-NCBI-DEFECTS, which the Owner decided on 2026-11-15); on
  2026-11-16 line RG item T5d wrote NCBI's exception where NCBI's read of a header record fails,
  which LOSAT refused, and on 2026-11-17 its fix 1 read the header records with a UTF8String or a
  tag in the long form as NCBI reads them, which it had refused, its fix 2 wrote NCBI's warning of
  such a record also where NCBI's fetches or its version check of an accession read it, and its
  fix 3 where NCBI's scope reads it for another sequence that holds a reported sequence's id; on
  2026-11-17 line DBID item v4late looked a Seq-align id up in a version 4 search's ISAM files in
  the volumes' order up to the first one that holds it, never a later one, as NCBI's lookup
  returns there, where LOSAT stopped with NCBI's `Error: Could not open input file` at a later
  volume's ISAM file that does not open, and line DBID item pfoff wrote outfmt 13 to 16 where a
  batch's prefetch is off and a query's report reaches a `Hit` whose Seq-align id's lookup throws
  (`DB contains no accession info.`, a version 5 database without the accession table in the
  search's list) as NCBI writes them, with the taxonomy database's texts at the loads that NCBI
  makes (none of the throwing `Hit` or of the `Hit`s after it in that report; a later `Hit` or a
  later batch's prefetch makes them), in every program, which version 1.3 refused; on 2026-11-18
  line DBID item pfoff-3 wrote outfmt 0 to 4 there as NCBI writes them: the description line of
  the Seq-align id and `Unknown`, then the stop at the first alignment of outfmt 0 that loads such
  a sequence under a Seq-align id that is not local (`GetDbType` does not look a local id up) or at
  any alignment of outfmt 1 to 4 that loads such a sequence, or outfmt 0's skipped line for each of
  its HSPs; the description lines load the sequences that they show, which in psiblast's rounds
  after the first are those of the two-part table of the first `-num_descriptions` HSPs (at least the first; psiblast ignores `-sorthits` and `-sorthsps`), and a report that stops the run ends it there, with no later batch's prefetch; in a call of LOSAT's library the refusal of outfmt 0 to 4 with a format of another writer is the call's, also where it searches once for each hit list size; on 2026-11-23 line DBID item pfoff-4 wrote outfmt 5 there as NCBI writes it: each `Hit` of such a sequence is the Seq-align id's FASTA form, `Unknown` and length 0, its HSPs the scores, the coordinates 0, the Seq-align's identities and no letters, and the report goes on (on 2026-11-24, at its merge of line PS item PS15b-1, which brought psiblast's outfmt 5, 8, 9 and 11, also psiblast's outfmt 5, each round's `Iteration`, and the refusal of a sequence that NCBI gives for such a Seq-align id without the lookup also at psiblast's PSSM read); on 2026-11-18 line NP item M4a-1c read the mask files with a string byte outside 0x20-0x7e and a
  `UTF8String` tag for a `VisibleString`, with NCBI's posts, which version 1.3 refused, and wrote
  makeblastdb's report for a damaged mask file, where LOSAT wrote its own text, keeping the
  refusals of a nesting deeper than 400 levels and of a string length past the end of the file
  that the host's memory decides, and its review the same day read NCBI's pointer syntax in mask
  files (`NULL` for an optional member, `@n` with NCBI's report, `:Type value` of the declared
  type), refusing a `:Type` of another type and a list element stated as `NULL` that makeblastdb
  dereferences; on 2026-11-20 line LIB read the local time of the creation date with the crate
  jiff (the maintainer's 「TZ=A」 of that day), which gives glibc's time in zone files with leap
  seconds, where version 1.3 listed a time 27 seconds later, and gives UTC for the forms of `TZ`
  that glibc reads outside POSIX and jiff's time for rules whose change falls in another year in
  UTC (tables of them; daylight time without a rule followed glibc's `posixrules`), which wait for
  the maintainer's decision on its pull request; on 2026-11-23 line NP item M (b) read the
  Seq-locs of a mask file nested to any depth, which LOSAT refused past 400 levels, as makeblastdb
  reads them with a stack large enough (PD-LOSAT-NCBI-DEFECTS approved exception 28, the Owner's
  「Q-MDEPTH=B」 of 2026-11-24); on 2026-11-24 line RG item D1c-1c-b wrote the
  sequences with a definition line without a Seq-id as NCBI writes them, which version 1.3 listed
  as not compared: those whose lines all lack one get no Seq-align, and beside lines with Seq-ids
  each format writes its `catch`, refused only at `sallacc` (NCBI's null pointer), in the organism
  report (outfmt 18, which waits for the maintainer's answer to question Q-D1c1cb-18) and in
  `rpsblast` and `rpstblastn`, not compared, and made LOSAT's refusal of a tabular row where NCBI
  stops or prints an uninitialised `staxid` (these, and the `sacc`, `saccver` and `staxid` rows of
  the tables) come after the reports of the query batches before in every program, as NCBI writes
  them (BLASTP, psiblast, TBLASTN and TBLASTX wrote none of them before, BLASTN none where it writes
  its reports after the loop), and after the failing batch's messages that NCBI writes to stderr
  before that row (its reader messages, posts, and for each of its queries up to the row's its own
  warnings, its formatter's taxdb warning and its subject reads up to the row; TBLASTN, TBLASTX and
  BLASTN dropped its reader messages, every program the rest); on 2026-11-27 line RG (series D) wrote
  such sequences of a profile database in `rpsblast` and `rpstblastn` as NCBI writes them, which
  item D1c-1c-b had left refused and not compared (refused at `sallacc` and in the organism report
  as in the other programs); the rest of this history
  under the maintainer's standing delegation (Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07))
- Status: Accepted by the maintainer on 2026-10-12, after the probe and design of the
  database search campaign (decisions P1 to P4 below, each as recommended). The decision
  that version 1.0 left open, the thread count of the parity target, was taken by the
  maintainer on 2026-10-13 (P1.1).
- Implementation status: stage DB-1 added the database reader (`LOSAT/src/seqdb`). Stage
  DB-2 added the subcommand `losat makeblastdb` (`LOSAT/src/writedb`, section "What `losat makeblastdb`
  of stage DB-2 writes and rejects"). Stage DB-3 adds `-db` and `-dbsize` to BLASTP (section
  "What BLASTP `-db` of stage DB-3 does and rejects"). Stage DB-4 adds `-db` and `-dbsize` to
  BLASTN (section "What BLASTN `-db` of stage DB-4 does and rejects"). Stage DB-5 adds `-db`
  and `-dbsize` to TBLASTN first (section "What TBLASTN `-db` of stage DB-5 does and rejects"),
  then to BLASTX (section "What BLASTX `-db` of stage DB-5 does and rejects") and to TBLASTX
  (section "What TBLASTX `-db` of stage DB-5 does and rejects"); stage DB-5 is complete (TBLASTN,
  BLASTX and TBLASTX), and all five programs search a database. Line DBID added databases with
  parsed Seq-ids (section "Databases with parsed Seq-ids"). The stage sections below describe
  what each stage added; later lines accept `-searchsp` (line DI-7) and the id and taxonomy
  filters (line DI-8), which those sections still list as rejected. The README lists what is
  refused today.

## Scope

How LOSAT searches a sequence database instead of a FASTA subject set: what the search
must reproduce, what a database is, how it is named on the command line, and what the
first release covers. The `-query` / `-subject` path is not changed by this decision.

## Background

Measured on 2026-10-12 with NCBI BLAST+ 2.17.0 (Swiss-Prot, 575,748 proteins; RefSeq 16S,
28,036 sequences):

- `-subject` already computes E-values from the size of the whole subject set, as NCBI
  2.17.0 does for `-subject`. LOSAT `-subject` gave the bytes of NCBI `-db` on 20 of 21
  probe cells; the 21st was a thread-dependent defect of BLASTP that is fixed.
- NCBI `-db` and NCBI `-subject` of the same sequences still differ in defined ways: the
  shortest subject length behind the hit-saving cutoff of blastp, blastx and tblastn (a
  database always answers 10, so `-db` keeps more weak hits), the bases a database stores
  in place of ambiguous nucleotides, and the `Database:` lines of the reports.
- Reading and encoding a large FASTA subject set costs a constant 0.3 to 0.7 s per run and
  2.7 times NCBI's memory on Swiss-Prot. A database read in blocks removes both.

## Decisions

### P1. Parity target

`losat <program> -db NAME` must give the output of NCBI BLAST+ 2.17.0 `-db NAME` on the same
database files, byte for byte, for outfmt 0, 6 and 7 (version 1.3: every format the program
writes, with and without parsed Seq-ids). `-db` and `-subject` of one FASTA
therefore differ exactly as NCBI's do. This replaces the earlier wording "E-values are
computed with the size of the whole database", which `-subject` already meets.

### P1.1. Thread count of the parity target (2026-10-13)

NCBI's own `-db` output can depend on `-num_threads`. With composition-based statistics,
blastp ends the redo of a query early by the hit lists of its query batch, and with more
than one thread NCBI keeps one set of hit lists per thread. Measured: a random
10,500-residue query followed by PUR4_HUMAN against Swiss-Prot gives 42 different output
rows at 8 threads and at 1 thread. P1 names "NCBI `-db` on the same database" without a
thread count.

Decision: LOSAT's `-db` reproduces NCBI `-db` run with `-num_threads 1`, and LOSAT gives
the same output at every thread count.

- NCBI's threaded redo (one set of hit lists per thread) is not reproduced.
- Comparison runs of NCBI `-db` use `-num_threads 1` and state the thread count.
- "LOSAT's output at 1, 2, 4 and 8 threads is the same" is a gate of every stage that adds
  `-db` to a program.

### P2. Database format

A LOSAT database is an NCBI BLAST database: the volume layout that makeblastdb writes.
LOSAT reads format versions 4 and 5 directly and writes both with `losat makeblastdb`
(version 5 by default, as makeblastdb; `-blastdb_version 4` for the older layout). There is no LOSAT-only format.

Reasons: the comparison oracle searches the very same files, so ordinal ids, the stored
2-bit bases and the title are shared by construction; existing databases are usable as they
are; the writer can be checked byte for byte against `makeblastdb -blastdb_version 4` and
`-blastdb_version 5`.

### P3. Command-line surface

- `-db`, `-dbsize` and `-searchsp`, under NCBI's names and with NCBI's meaning.
- A subcommand `losat makeblastdb` with makeblastdb's option names (`-in`, `-dbtype`, `-out`,
  `-title`, `-max_file_sz`).

### P4. Scope of the first release

- Databases without parsed Seq-ids: makeblastdb without `-parse_seqids`, and `losat makeblastdb`.
- One database per search.
- No id or taxonomy filters.

Databases distributed by NCBI (they store parsed Seq-ids), id and taxonomy filters, alias
files that filter or replace the statistics, and several databases in one `-db` are a
second step. Until then LOSAT rejects them explicitly; it never searches such a database
with part of its data ignored.

Version 1.2 (line DI-8f): several databases in one `-db` (`-db "a b"`, `-db '"a b" c'`) are
searched as one database, as NCBI's top alias node does: the names are resolved in order
(the first that is not found is NCBI's error), their volumes are opened in NCBI's order (by
base name, then path; independent of the order of the names), and a volume that two names
give is opened and searched once. The reports show one record per name, each the database
of that name opened alone (NCBI's `GetBlastDbInfo`): the outfmt 0 header joins the titles
with `; ` and sums the counts, the epilog has one block per name, a name given twice shows
twice. A version 5 list that mixes databases with and without the accession table of
`makeblastdb -parse_seqids` is searched as NCBI searches it, and a batch whose prefetch looks a
reported parsed Seq-id up in a file without the table stops as NCBI's does, with `BLAST Database
error: DB contains no accession info.` (exit 2) after the prologs and the reports of the batches
before (line DBID item mixdb2), and with the texts of the taxonomy database of the loads that NCBI's
prefetch makes before the stopping lookup, in the order of the Seq-id handles of the batch's Seq-align
ids (line DBID item loadord); refused is a run whose reports read such a Seq-id with the prefetch
off, a stopping batch where a load whose order to the stop NCBI takes from the addresses of
its Seq-id handles writes such a text, or where a lookup that throws on its id may come before the
stop, one where a load that NCBI makes only where it skips an earlier stopping lookup writes such a
text, one whose stopping ids NCBI may all never look up (another sequence holds an id of its
handle that the batch may load first or the run loaded before, in a batch, query or psiblast
round before, or in a call of the library's batches before; or a query that NCBI's scope holds at
the prefetch has it: the batch's records as read, those that `-query_loc` drops included, by
their `lcl|Query_N` or the ids of `-parse_deflines`, or the query of `-in_pssm`, and with outfmt 5
those of the batches before; of several stopping ids, the first that NCBI surely looks up ends
the loads), and one where a load whose order to the stop LOSAT does not know would make a lookup
of a tax id that a BLAST XML report before the batch left (outfmt 5: NCBI's XML stream cuts the
report before its hits at a query title or id outside US-ASCII, and it makes none of its loads,
line RG item K4h; the texts of such lookups come at the first load after the cut report that
makes the same lookup: a load of the reports or prefetches of the batches between the cut report
and the stopping batch, else a load of the stopping batch's prefetch before the stopping lookup,
or not at all, in every search program since review 5 of line DBID item loadord; the cut report
loads nothing, so the sequences that it shows count as loaded for no stop, since review 6)
(`PD-LOSAT-NCBI-DEFECTS`, rule 3; before
line DBID item mixdb, every such
list where NCBI starts the search; before item mixdb2, every run whose reports or prefetch
looked such a Seq-id up). `losat blastdbcmd` keeps one
database per `-db`. Since line RG (series D, D1a, 2026-11-12) a volume reached through an
alias file and again is read once and counted for each node that names it, as NCBI's walks
count it.

Line RG, D1a (2026-11-12): the alias keys `NSEQ`, `LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`,
`MIN_SEQ_LENGTH`, `FIRST_OID` and `LAST_OID` are read as NCBI's `CSeqDB` reads them: the counts
are NCBI's walks of the alias tree (a node with the key gives its value and is not descended
into), the statistics read `STATS_*` before `GetNumSeqs` and `GetTotalLength`, the reports
print the walks, and the OID ranges select the OIDs of the volume set (nested ranges ANDed,
two paths to one volume ORed, a taxid or seqid list ANDed with them). A value that NCBI's
`NStr` functions do not convert is NCBI's `Error: <what()>`, exit status 255, in NCBI's order of
the walks. `VDBLIST` is ignored (NCBI's `CSeqDB` does not read it). Since item D1a-2
(2026-11-13): a `FIRST_OID` or `LAST_OID` of 2^31 or more is NCBI's negative `int`; the OID
range of the only alias file of a `-db` is NCBI's iteration range, also with another name of the
same `-db` (a taxid or seqid list's search reads every listed OID, as NCBI's list forgets the
range); a `LENGTH` or `STATS_TOTLEN` over volumes without sequences is searched as NCBI searches
it (no hit, the empty-database warning where no count is given); `blastdbcmd -ipg` and
`-ipg_batch` write a PIG whose OID a range leaves out, as NCBI's formatter reads by OID; and a
database's errors come in NCBI's order (a list, the alias files' values, the first header, the
accession tables). Since item D1a-3 (2026-11-14): an `NSEQ` negative as NCBI's `int` or of
2^31 or more, and a `LENGTH` of 2^63 or more (a negative `Int8`), go through the statistics and
the reports as NCBI's `int` and `Int8` (`SeqDb::statistics`, `stats::setup`): the effective
lengths with the negative count, NCBI's test `(Uint4)(Int4)(total / count) == 0` of "The average
subject length is too short" (it fails for a count of 2^31 or a length of 2^32 times the count),
`Unknown error code 1` for BLASTX and TBLASTN without the Gumbel block and a negative average
(NCBI's hit-saving parameters), the Gumbel block's negative length (negative e-values, printed
0), no sort of the results at the end of the traceback where `GetTotalLength` is not positive (the
order of NCBI's one-thread traceback: by OID, then the hit list's heap), and the reports'
signed counts with NCBI's commas, which lose digits from 4,294,967,296,000,000,000. Refused,
where NCBI's arithmetic overflows: the total lengths of a `-db` list's databases summed past
`Int8` (outfmt 0 to 5 and 13 to 16) and their numbers of sequences summed past `int` (outfmt 5).
A summed length of exactly 2^63 (NCBI's Int8 -2^63), which NCBI's serial writer prints as `-0` in
BLAST XML and XML2 (outfmt 5, 14, 16), the strategy and the archive, is written as
`-9223372036854775808` (PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits,
the Owner's answer Q1 = A of
2026-11-17 to D1a-3's question, DECISIONS OQ-1117 of line RG; refused before). The requests of `-outfmt 11` and
`-export_search_strategy` compute their `EffectiveSearchSpace` from the negative length as NCBI does
(the query length of context 0; probe run2 of D1a-3). The preliminary stage of every search reads
NCBI's option spaces of the statistics before their fallback to OID 0
(`SplitQuery_SetEffectiveSearchSpace`, prelim_stage.cpp:219-220; `stats::setup::
preliminary_option_space`): for a negative length each context's query length, for a length of 0
none (with `-dbsize` the spaces of `-dbsize` and the walk's count), while the traceback and the
reports read the setup's. Where `GetTotalLength` is 0 (an alias file's `LENGTH 0`, also with
`STATS_*`, or a list whose count visits no sequence) NCBI's preliminary search computes the
statistics of each subject it searches (`BLAST_OneSubjectUpdateParameters`, blast_engine.c:1434-1443;
TBLASTN and TBLASTX skip the shorter subjects first), and so does LOSAT in the six programs: the
effective lengths, the hit-saving cutoffs without the composition-based stretch of the e-value, the
word cutoffs and the link cutoffs of the subject's length; the e-values of the traceback and the
reports' counts and spaces are the setup's (Spouge's database factor 1), and the order of the
results is the one-thread traceback's. The overflowing sums of BLAST XML, XML2 and JSON are refused
where NCBI formats them (after the search: its own errors of the search come first). Line RG item
D1a-4 part A (2026-11-17): RPSBLAST and RPSTBLASTN read these keys as NCBI's RPS search does
(`CLocalRPSBlast`, rpsblast_local.cpp:331-400; one owner, `rps_search::volumes::RpsUnits`): a
database of one volume is one search unit, the `-db` name with all its alias keys; one of several
volumes searches each volume bare (its own length and count, no `STATS_*` key, the shortest length
10, no OID range) under the database-wide `DbLength` (`STATS_TOTLEN`, else `GetTotalLength`) and
`DbSeqNum` (`STATS_NSEQ`, else `GetNumSeqs`) that `x_AdjustDbSize` sets without `-searchsp` and
`-dbsize`, each read only where above 0, while the Gumbel block reads `DbLength` where it is not 0.
The footers print the first unit's search space, the search strategy's requests after the first hold
that `DbLength` (0 and negatives too) and the last unit's search space, and the traceback sorts the
results before its filters whatever the total length and after them (`-max_hsps`, `-qcov_hsp_perc`)
only where the unit's `GetTotalLength` is above 0 (blast_traceback.c:1306, 1768-1771; part B, review A
MF2). The preliminary stage's hit cutoff without the Gumbel block reads the space that
`SplitQuery_SetEffectiveSearchSpace` sets (prelim_stage.cpp:219-221: the unit's total before the
fallback to OID 0, its count without it), the traceback the setup's: they differ over one volume
whose total is negative, and with `-dbsize` over one volume whose total is 0 or negative
(SplitQuery's count is the unit's, the setup's 1; part B, review A MF1, review B N2). The setup of
every unit fails as NCBI's where its average subject length is below one letter without the Gumbel
block (an alias file's `NSEQ` above its `LENGTH` over one volume: `BLAST engine error: The average
subject length is too short`, exit status 3, both programs; fix round 1, review A MF5). The search
ignores `FIRST_OID` and `LAST_OID` (NCBI's lookup table holds every profile of a volume); a hit of a
profile that they leave out among the profiles the reports show (those of the first
`-max_target_seqs` subjects of a query's hit list, following profiles of one Seq-id counting as one
subject, as NCBI's formatter counts them) is written as the included copy of its record that holds
its Seq-id where there is one (NCBI's lookup by Seq-id finds it; fix round 2), and is refused where
NCBI's report does not find it (`Unknown`, or a null pointer exception: the D1a-4 entry of
PD-LOSAT-NCBI-DEFECTS)
or finds another record, RPSBLAST before any output and RPSTBLASTN after the reports of
the query batches before; one past that cut of several volumes' merged lists is NCBI's run, its
HSPs counted by outfmt 7 (fix round 1, review A MF3). A unit's concatenated length `Lc` (`LENGTH`
plus `NSEQ` over one volume, NCBI's `Int4` of the `Int8` sum; a bare volume's rows) sizes the
diagonal table, bounds the ungapped extensions to the right and ends the gapped stage's interval
tree, below the profiles, above them or negative, as in NCBI (blast_engine.c:1129-1149; part B); at
the rows plus 1 the extensions read the row after the profiles as NCBI does (the zeros of the
`.rps` file's last memory page: DECISIONS D1a4-9). Refused (the D1a-4 entry of PD-LOSAT-NCBI-DEFECTS,
rule 3): a
sum `LENGTH` plus `NSEQ` above 2,147,483,647 (NCBI frees the profiles twice and aborts), an `Lc` whose
diagonal table passes 2^30 entries (`Lc` plus the window above 2^30: NCBI's loop never ends), both
where NCBI's engine meets them (the unit's search of the first query batch NCBI searches: an empty
query file or a first batch without letters ends as NCBI's; fix round 1, review B MF1), a query
batch where two gapped alignments of one tree start past `Lc + 1` (NCBI's interval tree splits
without end) and one where `Lc` is above the rows and an ungapped extension reads the row after the
profiles past the `.rps` file's last page or goes on past that row (NCBI reads past its score
matrix), RPSBLAST before any output and RPSTBLASTN after the reports of the query batches before, a
setup count of 0 without the Gumbel block (`NSEQ 0` over one volume with `-searchsp` or `OLD_FSC`
where a query is searched: SIGFPE; a run without a valid query is reported as NCBI reports it, fix
round 3, review C MF1; a bare volume counts its profiles). With `-export_search_strategy` the requests after the
first (a run of two query batches or more) hold `DbLength` = D, also -2^63 (`LENGTH` or
`STATS_TOTLEN` 2^63 over several volumes), which NCBI writes as `-0` and LOSAT as
`-9223372036854775808` (approved exception 17, Int8 -2^63 written as digits; with one query batch
NCBI's only request holds no
`DbLength`: fix round 1, review A MF4; refused until fix round 2). Refused for
good: a statistics count of 0 (an alias file's `NSEQ 0`, or its `LENGTH` or `STATS_TOTLEN` over
volumes without sequences without `NSEQ` or `STATS_NSEQ`) in a search without the Gumbel block
(BLASTN, TBLASTX, and BLASTP, BLASTX, TBLASTN, PSI-BLAST with `-searchsp` or `OLD_FSC`, BLASTX
also with `-ungapped`), where NCBI divides by zero (SIGFPE) as a query is searched. A run of
queries without a valid one, which NCBI reports without dividing, is reported as NCBI's by
BLASTX and, since fix round 3 of D1a-4, by BLASTP and
PSI-BLAST (rows `db.ak.bp.nseq0.oldfsc.allx.7`, `db.ak.bp.nseq0.ssp.allx.7` and
`db.ak.ps.nseq0.oldfsc.allx.7`) and, since D1a-5 of line RG, by BLASTN, TBLASTN and TBLASTX (rows
`db.ak.bn.nseq0.alln.*`, `db.ak.tbx.nseq0.alln.*`, `db.ak.tbn.nseq0.oldfsc.allx.*` and
`db.ak.tbn.nseq0.ssp.allx.7`; they refused it at the database's open, probe e8 of fix round 3 of
D1a-4). Every program makes the refusal with one test (`common::database::check_setup_count`) in
the setup of the first query batch (BLASTN: or query chunk) with a valid context, where NCBI
divides, before it reads a sequence (also over a volume without its `.nsq`): BLASTN, TBLASTX and
TBLASTN after the reports of the batches before (NCBI's output before its crash in the next batch
in the tabular formats: 150 and 152 bytes of outfmt 7 with `BATCH_SIZE=1`, two queries without
letters to search then a valid one; the crash cuts NCBI's buffered output of the other formats
short), BLASTP and PSI-BLAST before any output. With the Gumbel block, a TBLASTN database whose `GetMinLength` is
below 3 (an alias file's `MIN_SEQ_LENGTH` 1 or 2) fails as NCBI's does: `BLAST engine error: The
average subject length is too short`, exit status 3.

Version 1.3 (line DBID): databases with parsed Seq-ids are searched, NCBI's distributed
databases included (section "Databases with parsed Seq-ids").

## Not part of this document

- LOSAT Web taking a database as engine input (accepted on the same day) is written into
  [`PD-LOSAT-WEB-APP-BOUNDARY`](PD-LOSAT-WEB-APP-BOUNDARY.md) with the stage that
  implements it (DB-7).
- The README statement on large database searches is decided when the benchmark stage
  (DB-6) has numbers.
- The BLASTP seed scan on large databases is a performance package of its own; its order
  relative to the other work is the maintainer's.

## What the reader of stage DB-1 rejects

These follow from P4 and from the files the first stages do not read; they are
implementation scope, listed here so that a rejection is not mistaken for a defect. Each
is an explicit error of `LOSAT/src/seqdb`:

| Rejected | Why |
| --- | --- |
| A header with a Seq-id other than `gnl|BL_ORD_ID|n` of the sequence's number n in its volume | parsed Seq-ids (P4); searched since version 1.3, section "Databases with parsed Seq-ids" |
| Alias keys `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `MEMB_BIT`, `MASKLIST` (all read) | filters (P4); `GILIST`, `TILIST` and a version 4 database's `SEQIDLIST` read since line RG D1c-1 (2026-11-22; the listed OIDs searched and the definition lines that no list names hidden, as NCBI's; a TI list finds no OID, as NCBI BLAST+ 2.17.0's; NCBI's exception where a `ti` general id's tag is no number; a `SEQIDLIST` token NCBI warns about, an id NCBI's simplification throws on or changes a second time over two volumes refused), a version 4 database's `TAXIDLIST` read nowhere as NCBI's, a version 5 database's `SEQIDLIST` and `TAXIDLIST` read through LMDB since line RG D1c-2 (2026-11-25; Pass 2 of `CSeqDBOIDList::x_Setup`, combining with Pass 1 filters and user search lists, rescuing definition lines by alias `TAXIDLIST` before negative user taxids, refusing a `SEQIDLIST` comparison that posts a Seq-id warning, and refusing an alias `SEQIDLIST` or search `-seqidlist` with two differently spelled Seq-ids that `s_CompareSeqId` treats as equivalent when combined with a search seqid list, `PD-LOSAT-NCBI-DEFECTS` rule 3); `MEMB_BIT` read since line RG D1d (2026-11-26; parsed as NCBI's unsigned integer between `SEQIDLIST` and `TAXIDLIST`, applied per volume with NCBI's `MEMB_BIT error: conflicting bit found.` on conflicting bits, and filtering each definition line's first membership bitmask before id lists, taxid lists and `OID_MASK_TYPE`, with a sequence left with none searched and not reported as NCBI's and an RPS profile of a database of several volumes refused where reported; values past 32 refused when testing a line with memberships because C++ shifts past a 32-bit integer); `FIRST_OID`, `LAST_OID` read since line RG D1a, `OID_MASK_TYPE` since D1b-2 (2026-11-21; its bit 0 hides the definition lines with a predicted RefSeq id, a sequence left with none is searched and not reported, as NCBI's (an RPS profile of a database of several volumes refused where reported, NCBI's report not finding it); a sequence id that NCBI's accession guide warns about under it refused, NCBI repeating the warning at each read), `OIDLIST` since D1b-1 (2026-11-21; an empty file or one of its header alone refused, LOSAT's words where NCBI crashes or names a line of its source), `VDBLIST` ignored as NCBI's reader ignores it; `MASKLIST` names the GI-indexed masks that `makeblastdb -gi_mask` writes (line NP item M4c-1, 2026-11-20), which every reader reads as NCBI's does since item M4c-2: the masks of a database named alone whose own alias file has it (`blastdbcmd`, the `-db_soft_mask` / `-db_hard_mask` of blastn, blastp, blastx, tblastn and tblastx, `makeblastdb -input_type blastdb`, `makembindex -db_mask` and `-show_filters`, whose new-style index opens each volume alone and so finds no `MASKLIST` in a database of several volumes, `unrecognized filter algorithm id` as NCBI's; `LOSAT/src/seqdb/gi_mask.rs`), passed over elsewhere as NCBI's reader passes it over; a `MASKLIST` name that is not UTF-8 read by its bytes as NCBI's (line RG item K10a, 2026-11-24, which withdrew the refusal of review N8); makeblastdb's own reads of an alias file (the removal of the old database, the metadata after the build) pass it as NCBI's do |
| Alias keys `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`; `NSEQ` or `LENGTH` whose value is not the number of sequences or the total length of the alias file's databases | read since line RG D1a (2026-11-12) as NCBI's walks, negative and wrapped counts since D1a-3 (2026-11-14); the values listed above as refused |
| A name that lists several databases | one database per search (P4); searched since version 1.2 |
| A volume listed more than once, a group alias file (`index.alx`) in a directory, two volumes of one name in different directories | NCBI derives counts or the volume order from them in ways the reader does not reproduce; since version 1.2 two volumes of one name in different directories are read (NCBI orders them by base name, then path); since line RG D1a a volume listed more than once is read once and counted per alias node, as NCBI does |
| Column files of a volume (mask data of `makeblastdb -mask_data`) | not read; searched since 2026-11-11 (line SW, SW8 follow-up 2): NCBI's searches open no column file without a database mask named by `-db_soft_mask`, `-db_hard_mask` or an imported search strategy (`-import_search_strategy`), so LOSAT searches the volume without reading them, whatever they hold; since 2026-11-14 BLASTN reads the mask column for `-db_soft_mask` and `-db_hard_mask` (line NP, M4a-2a and M4a-2b; LOSAT refuses the third) |

Two errors are NCBI's own and are reproduced: volumes of version 4 and of version 5 in one
database, and a version 5 volume whose LMDB file is missing.

Left unread without an error: the id and taxonomy files of a version 5 database, which serve
lookups by id or taxonomy only, the LMDB file's volume list, which NCBI compares with the
volumes, and (since 2026-11-11) the column files of a volume, which NCBI's searches read only for
a database mask named by `-db_soft_mask`, `-db_hard_mask` or an imported search strategy (since
2026-11-14 BLASTN reads the mask column for `-db_soft_mask` and `-db_hard_mask`; LOSAT refuses the third). The LMDB file is opened when the
database is opened. Since version 1.3 a search looks
the Seq-align id of a reported sequence with parsed Seq-ids up in the LMDB file or the ISAM
files, as NCBI's formatters do. Since stage DB-3 a `-db`
name is looked up along NCBI's search path (the current directory, then `BLASTDB`, then
`[BLAST] BLASTDB` of `.ncbirc`), with NCBI's texts when it is not found. Not ported, each an
error or a difference that a test records: group alias files (`index.alx`), several names
in one `-db` (searched since version 1.2), index files cut short, NCBI's threshold of open files and its warning at
1,024 links, escapes in registry values.

## What `losat makeblastdb` of stage DB-2 writes and rejects

`losat makeblastdb -in FASTA[ FASTA...] -dbtype prot|nucl -out NAME [-title T] [-max_file_sz S]
[-blastdb_version 4|5] [-logfile F] [-parse_seqids] [-hash_index] [-taxid N | -taxid_map F]
[-mask_data F[,F...] [-mask_id I[,I...] [-mask_desc D[,D...]]]]`
writes the files of `makeblastdb` (format version 5 unless `-blastdb_version 4` is given): sequence files
(`.psq`, `.nsq`) and header files (`.phr`, `.nhr`) byte for byte, the index file (`.pin`,
`.nin`) but for the creation date, for several volumes the alias file (`.pal`, `.nal`)
but for its date line, and the JSON metadata file (`.pjs`, `.njs`) but for the date
it records. Version 5 adds, once for the whole database, the two LMDB files (`.pdb` with
the volume tables and `.ptf` with the tax ids; `.ndb`, `.ntf`) and the two tax-id lookup
files (`.pot`, `.pto`; `.not`, `.nto`), page for page. With `-parse_seqids` the header files hold
the Seq-ids of the definition lines as makeblastdb parses them (its accession guide included),
duplicate ids stop the build with makeblastdb's message, and the id indexes are written: of
version 4 the string and GI indexes (`.psi`/`.psd`, `.pni`/`.pnd`), of both versions the OID-to-GI
file (`.pog`), of version 5 the accession keys in the LMDB file (`acc2oid` of `.pdb`) and the
OID-to-ids file (`.pos`). `-hash_index` adds the sequence hash index (`.phi`/`.phd`, `.pog`),
`-taxid` and `-taxid_map` (which needs `-parse_seqids`) the tax ids of the definition lines and of
the version 5 tax-id files. As in makeblastdb, `-taxid_map -` reads the map from the standard
input (before the FASTA input, so with `-in -` the map takes all of it), a map file whose reading
fails (a directory) keeps the lines read before the failure without a message, and a plain FASTA
id with a control byte or a byte from 0x80 up is not a valid local id: the definition line keeps
the id up to that byte. `-mask_data` (with `-mask_id` and `-mask_desc`) reads the masks of
dustmasker, segmasker, windowmasker and convert2blastmask, in ASN.1 text or binary, and writes the
mask column of each volume (`.?aa`, `.?ab`, `.?ac`, version 4 and 5 alike) as makeblastdb does;
a mask reaches a sequence by the same Seq-id handle as in NCBI (so `-parse_seqids` is needed for
masks that name accessions, and the log, the "Mask file:" lines, and the messages of the masks
are makeblastdb's). Several input files can be named in
`-in` (as makeblastdb splits the value), and `-logfile` redirects the log and the messages
to a file (appended); as in makeblastdb, `-logfile -` keeps them on the standard streams,
and `-logfile /dev/null` and `-logfile ""` discard the messages (the empty name then ends
the build with exit status 1, silently). The FASTA file is read as makeblastdb reads it (its reader settings
differ from those of `-subject`). The log on the standard output and the messages on the
standard error stream are makeblastdb's, in makeblastdb's order, and so is the exit status
of a build and of makeblastdb's own errors (0; 1 for an input or option error; 2 when the
finished database does not resolve; 4 when a sequence does not fit in memory; 255
otherwise). The differences are listed below.

Platforms. The expected files of the tests and gates are makeblastdb 2.17.0's on Linux
x86-64. Where NCBI's source asks the platform, LOSAT asks it too and does what that source
does there (`LOSAT/src/writedb/platform.rs`): the page size of the LMDB files
(`GET_PAGESIZE`: 4 KiB on Linux x86-64 and Windows x64, 16 KiB on macOS on Apple silicon,
the kernel's on Linux arm64, at most 32 KiB), their default map sizes (Windows: 500000 bytes
for both files, elsewhere 700000000 and 300000000), the length of an LMDB file on Windows
(LMDB sets it to the map size), and `sizeof(long)` in the header size of an index file (8 on
Linux and macOS, 4 on Windows), which moves the padding of the date and the cut of a volume
under `-max_file_sz`. The Web app's builds (wasm32) write what makeblastdb writes on Linux
x86-64. The output on Windows and macOS has not been compared with a makeblastdb run on those
platforms; the tests there check LOSAT's platform choice and, on the Linux x86-64 settings,
the same bytes as on Linux. Errors that makeblastdb reports through its general handler
(`Error: NCBI C++ Exception:` and one line per exception, `T0 "<file>", line <n>: Error:
(<class>::<code>) <function>() - <message>`) name the source files in the build tree of the
official Linux x64 release of makeblastdb 2.17.0, on every platform, as LOSAT's other NCBI
exception reports do; a makeblastdb built elsewhere (bioconda, macOS, Windows) names its own tree.

| Not written or rejected | Why |
| --- | --- |
| `-gi_mask` with two algorithms whose GI-mask files are the same files (one `-gi_mask_name` twice, or two names of one path) | written since line NP item M4c-1 (2026-11-21) as makeblastdb 2.17.0 writes them, but for the alias file of a single-volume build, which names the volume by its file name where NCBI's names the database itself (`PD-LOSAT-NCBI-DEFECTS` approved exception 25, 1.81), and where both algorithms would write one data file at the same time: an explicit rejection of a deterministic result (`PD-LOSAT-NCBI-DEFECTS` 1.81, rule 2), whose bytes depend on how the C++ runtime buffers two streams on one file; files written one after the other (an algorithm that reaches no GI writes none) are NCBI's |
| `[BLAST] MAX_SEQID_LENGTH` (registry or `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`) with `-parse_seqids`, when makeblastdb uses the value (the first layer that has the entry gives a number `NStr::StringToUInt` converts, not 0; an empty, unconvertible or 0 value is ignored by makeblastdb and accepted by LOSAT); `[BLAST] LONG_SEQID` | they change which FASTA ids makeblastdb parses; LOSAT keeps NCBI's defaults and refuses to run |
| `-input_type` other than `fasta` and `blastdb` | FASTA and existing databases (`blastdb`, line RG item D3c) only (a FASTA input plain, or named `.gz`, `.bz2`, `.zst`: inflated while it is read, in any case of the suffix, as makeblastdb decides it; `-in -` reads the standard input as it comes). Every gzip member, bzip2 stream and zstd frame of a concatenated file is read, where makeblastdb reads every gzip member but only the first bzip2 stream or zstd frame and drops the rest silently (`PD-LOSAT-NCBI-DEFECTS`, the makeblastdb exception). A `.gz` file is read as makeblastdb's decompressor (`CZipDecompressor`) reads it (`GzipMembers` in `LOSAT/src/blastinput/compressed.rs`, which inflates with zlib-rs, a Rust port of zlib's inflate, where makeblastdb calls the host's zlib: the official Linux build links `libz.so.1`, and a host whose `libz.so.1` is zlib-ng's zlib-compatible build could in principle differ on damaged data): neither a member's trailer (CRC-32 and length) nor its header CRC is checked, a trailer cut short reads as the whole member, a member's header must end within the next 4 KiB, bytes after the last member are read only where makeblastdb reads them (fewer than 4 KiB that inflate to no byte, such as up to four zero bytes or a line feed), and a file that does not start with a gzip header is read as a zlib stream, up to that stream's end. A file that is not what its name says, or is cut short (for `.gz`: where makeblastdb's decompressor fails too), stops the build with `losat makeblastdb: -in: 'x.gz' is named .gz but is not gzip data (...)` (exit status 1, the files removed as for another failed input); makeblastdb warns about the decompressor's error, says the file does not appear to be FASTA, and ends with `BLAST failed to write output` (exit status 6) |
| NCBI settings that change makeblastdb's messages or stored titles (`DIAG_*`, `NCBI_CONFIG_*`, registry entries outside the accepted list, a non-empty `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS`) | as for the search programs: LOSAT refuses to run instead of giving a different result |
| Replacing a database of the same name whose alias file the reader of stage DB-1 rejects | the files to remove are found by resolving the alias file; LOSAT stops before it removes anything |

### Where `losat makeblastdb` knowingly differs from makeblastdb

Removal of an existing database and the paths of the output:

- An output path with a space or a double quote in its absolute form (for example a build
  inside a directory `my data`). makeblastdb reads the path as a list of names when it looks
  for the volumes of an existing database: it removes only that database's LMDB files
  (`.pdb`, `.pos`, `.pot`, `.ptf`, `.pto`, and logs "Deleted existing ..."), keeps its volume
  files, writes the new files over or beside them and ends with exit status 2 when it opens
  the result. LOSAT
  takes the path as one name, removes the whole database (a directory at one of its names with
  its contents), replaces it and ends with 0; a build
  that fails after the removal (a `-max_file_sz` that is no size) leaves none of the old
  database, where makeblastdb leaves its volume files. (`-out` itself may not contain a space,
  in both programs.)
- An existing `sub/db.pal` while the current directory holds volumes of the same names
  (`./db.00.*`). NCBI looks for a volume in the current directory first, so makeblastdb
  removes the files of the current directory. LOSAT removes the volumes beside the alias
  file (a directory at one of those names with its contents). Other behaviour of makeblastdb is kept: the existing database is removed before the
  input is read, so a failing build also removes it; a directory at one of the names removed
  is removed with its contents, read-only parts included, and a symbolic link is removed, not
  followed (`CDirEntry::Remove`); a volume file's final name is cleared the same way before
  the volume is renamed to it (a directory `NAME.pin` that is no database is removed then).
  Where the removal fails (the directory holds a directory that cannot be read), the volume
  file keeps its numbered name and the build goes on, as in makeblastdb, which does not look at
  the result of the rename: it ends where the database is read (`No alias or index file found
  ...` or `Could not find volume or alias file ...`, exit status 2) when the file is the
  volume's `.pin`, `.nin`, and ends normally when it is another file (`.phr`, `.psq`, an
  index file).
- An output file that cannot be created (a directory of its name that the removal of the old
  database could not delete) or written (a full device): LOSAT writes it as makeblastdb does,
  through a model of the `ofstream` of makeblastdb 2.17.0 (the official Linux x64 release,
  which links libstdc++ 13.2 statically: a put area of 8191 bytes, writes that do not fit go out
  at once, a failed write leaves the stream failed and every later output a no-op; the C++
  libraries of makeblastdb on macOS and Windows are not compared), and does not test it, as
  makeblastdb does not. The build goes on;
  a volume file under its numbered name is renamed whatever is there (a directory, a link as
  the link); the new database is opened as makeblastdb opens it, so a missing index or alias
  file ends with `No alias or index file found ...` or `Could not find volume or alias file
  ...` (exit status 2) and any other missing file is left out of the metadata file. A device
  that fills during the build leaves regular files cut short: the totals of the metadata file
  are those of the volumes the new database resolves to (a cut alias file names fewer than the
  writer made), which LOSAT takes from its writer (below): they are the counts the header of
  each index file holds (compared with makeblastdb on files cut by `RLIMIT_FSIZE` at 100 bytes
  to 100 KB: equal). The metadata file's serializer reports what makeblastdb's reports, at the same byte: `Error:
  (802.5) [makeblastdb] CObjectOStream: error at byte N: stack is empty: write fault` (exit
  status 0) when a write of its 4096-byte buffer fails, and, when the stream's flush fails,
  that error with `cannot close output stream`, the report of the failed close, and `BLAST
  failed to write output: COStreamBuffer::Flush: failed` (exit status 6). The errors of
  `close(2)` are not seen (makeblastdb does not test the close either). The data files of LMDB
  (`.pdb`, `.ptf`) are not `ofstream`s: a write of one that fails (a device that fills) ends
  LOSAT with `Error: <path>: <system error>` (exit status 255), where makeblastdb prints `BLAST
  Database creation error: mdb_txn_commit: Input/output error` (or terminates with that text)
  and then aborts or crashes (an NCBI crash, not reproduced). A duplicate id that the close finds in
  the last sequence removes the files of the build and then runs the writer's destructors in
  makeblastdb's order (the `.pot` and `.pto` files and the commit of `.ptf`, the `.pos` file and
  the commit into the removed `.pdb`, the close of the last volume, whose ISAM and OID to GI
  files are made under its numbered name). Where one of these destructors throws, makeblastdb
  aborts and the steps after it do not run; LOSAT ends with the build's error there (the lock
  files: below).
- A volume that an existing alias file names and that is found only through `BLASTDB` or
  `.ncbirc` is not removed (`losat makeblastdb` does not use the search path of `-db`).
- Directories that makeblastdb creates get mode 0755; LOSAT does the same on Unix only.

Version 5 and the options added to the writer:

- The `dbname` of the metadata file is the file name of the output path, which LOSAT takes
  as one name (makeblastdb splits a path with a space and ends with exit status 2 before it
  writes the file). A file `NAME.pjs` left by a replaced database of several volumes is
  listed and counted, as in makeblastdb.
- The totals of the metadata file come from the volumes the writer made that the new database
  resolves to; the volumes are not opened afterwards (opening them would hold the index of
  every sequence in memory). A volume the resolution finds that the writer did not make (one of
  the same name in the current directory, which NCBI's reader finds first) adds nothing.
- `-logfile` is read as makeblastdb's option, not wherever it stands in the command line
  (makeblastdb's corelib also finds it as the value of another option).
- LMDB's lock files (`NAME.pdb-lock`, `NAME.ptf-lock`) are made as makeblastdb makes them:
  LMDB creates each when it opens its environment, before the data file, with 8192 bytes of a
  lock region that no transaction has used, and the writer removes it when it is destroyed,
  also when the build ends with an error that makeblastdb reports. Where makeblastdb aborts
  instead (the crash inputs listed below, and a write that fails while its LMDB writer is
  destroyed, for example on a full disk) it leaves both lock files, with the `.pot` and `.pto`
  files; LOSAT ends those builds with an error and removes the lock files. A lock file stays
  where makeblastdb leaves it: when LMDB has made the lock file and cannot finish opening the
  environment (a directory of the data file's name that the removal of the old database could
  not delete: `BLAST Database error: LMDB runtime error: mdb_env_open: Is a directory`, exit
  status 2; or a data file that cannot be written, `No space left on device`), the lock file of
  that environment is left, and with it the `NAME.pdb` of the other environment when it is the
  `NAME.ptf` that cannot be made. A lock file that LMDB cannot create stops LOSAT where it
  stops makeblastdb: an output name of 247 bytes or more, whose `NAME.pdb-lock` passes the
  file-name limit of 255 bytes, ends with "BLAST Database error: LMDB runtime error:
  mdb_env_open: File name too long" (exit status 2) when the second sequence is added, or that
  creation error (exit status 1) for one sequence, with no file written, in both programs; a
  directory of a lock file's name that the removal of the old database could not remove (it
  holds a directory that cannot be read) ends the same way with "Is a directory". The bytes of
  the lock file are those of the official Linux x64 release, which keeps robust
  process-shared POSIX mutexes in it; the lock file of makeblastdb on macOS (named
  semaphores) and on Windows (named mutexes) holds other data, which LOSAT does not
  reproduce, and neither platform has been compared with a makeblastdb run: LMDB's
  `CreateFileW` gives a Windows error code, where LOSAT reports the text of the system's error
  number, and LOSAT's paths are not bound by `MAX_PATH`, which binds LMDB's `CreateFileW`.

Errors that makeblastdb does not detect, or detects elsewhere:

- On Windows, an output directory that exists but is not writable: makeblastdb stops before
  its log (it reads the directory's access rights). LOSAT does not read them: a version 4
  build writes its files to nowhere, as any file that cannot be created, and ends where the
  database is opened (`No alias or index file found ...`, exit status 2); a version 5 build
  stops where LMDB cannot create its lock file (`LMDB runtime error: mdb_env_open: ...`, as
  above). From the code; not run on Windows.
- Usage errors (unknown option, missing `-dbtype`, `-taxid -1`, `-taxid_map` without
  `-parse_seqids` or with `-taxid`, `-mask_id` without `-mask_data`, `-mask_desc` without
  `-mask_id`) end with exit status 1 and NCBI's usage and error lines, as makeblastdb does
  (since line RG item K4i-2b, 2026-11-23; before, LOSAT's option parser's text and exit status 2).
- With more volumes than the process may hold open files, makeblastdb writes the database
  and ends with exit status 2 ("Too many open files"). LOSAT resolves the finished
  database without opening its volumes and ends with 0.
- An allocation failure outside the reading of the input ends LOSAT abnormally; makeblastdb
  prints "BLAST ran out of memory" (exit status 4), as LOSAT does when a length written in
  the input does not fit.

Inputs on which makeblastdb 2.17.0 crashes, refused with an error (exit status 1):

- Version 5 with `-parse_seqids`: a sequence whose ids are all GIs between two sequences with
  other ids (makeblastdb aborts when it writes the `.pos` file); LOSAT names the sequence.
  GI-only sequences before the first sequence with other ids or after the last one are written
  as makeblastdb writes them (its `.pos` file then lists the ids from OID 0 on, a deterministic
  NCBI result).
- Version 5 with `-parse_seqids`: a duplicate id inside the first sequence. makeblastdb prints
  the duplicate message and then aborts in its tax-id writer; LOSAT prints the same message.
- Version 5 with `-parse_seqids`: an accession key longer than 511 bytes (a general id with a
  long database name), which LMDB refuses (`MDB_BAD_VALSIZE`); makeblastdb aborts after its
  `.pos` file.
- Version 4 with `-parse_seqids`: a string-index key longer than 256 bytes (makeblastdb stops
  with "buffer overflow detected").
- A `-taxid_map` line of white space only (makeblastdb crashes).
- Version 5 with `-parse_seqids`: `LMDB_SPLIT_CHUNK_SIZE` or `LMDB_MIN_SPLIT_SIZE` that is not
  a number of the range of `unsigned int` (makeblastdb aborts in its LMDB writer, after the
  volumes are written, and leaves `NAME.pdb`, `NAME.pos`, `NAME.phr`, `NAME.pin`, `NAME.pog`,
  `NAME.psq`, `NAME.ptf` and the two lock files, without `NAME.pot`, `NAME.pto` and
  `NAME.pjs`); LOSAT stops with
  `BLAST Database creation error: Cannot convert string '<value>' to ...` (exit status 1),
  removes the volumes and leaves `NAME.pos` and `NAME.ptf`.
- `-in ""` (no name): makeblastdb reads the first name of an empty list after the log line
  "Maximum file size" and `-taxid_map`, and crashes (a segmentation fault). LOSAT refuses it
  there, so that an existing database of the output name has been removed first, as
  makeblastdb removes it, and the errors of `-max_file_sz` and `-taxid_map` come first.

Reproduced although it looks wrong (deterministic NCBI behaviour): a file that is already at the
numbered name of a file the writer never created (`NAME.00.pnd`, `NAME.00.ppi` of a protein
database with `-parse_seqids`, an ISAM file without entries) is given the final name as well
(`NAME.pnd`, `NAME.ppi`), and is counted among the files of the metadata file; a file of the
names makeblastdb does not make for the build (`NAME.00.pog` without `-parse_seqids`, `NAME.00.psd`
of a version 5 database, trace index files) is left alone. With `-hash_index`, no
`-parse_seqids` and several volumes, the first sequence of a later volume is indexed under the
ordinal id it would have had in the earlier volume, so makeblastdb and LOSAT both stop with
"Duplicate seq_ids are found: GNL|BL_ORD_ID|n" when that volume holds more than n sequences.
`-taxid_map` lines end at a line feed on every platform (makeblastdb on macOS also ends them at a
carriage return).

Mask data (`-mask_data`): every Seq-loc alternative of a mask is read (`null`, `empty`, `whole`,
`int`, `packed-int`, `pnt`, `packed-pnt`, `mix`, `equiv`, `bond`; a `feat` or a Seq-loc of several
Seq-ids ends the run with NCBI's "Cannot get masked range Seq-id"), several Seq-locs of one
sequence are combined as makeblastdb combines them, empty and whole ranges included, and a `giim`
Seq-id finds the masks that makeblastdb gives it (which follows the order of the mask files). What
makeblastdb reads and LOSAT refuses (`losat makeblastdb: ...`, exit status 1, without
makeblastdb's closing "No sequences matched any of the masks provided."):

- Ranges whose order NCBI takes from memory addresses: ranges of two different Seq-id handles
  with the same packed number (`gi|5` and `gb|AY000005.1|`; a local id and a `general` string
  tag; two `giim` handles of one number) written on one strand of one sequence and algorithm,
  when the two orders write different ranges. That order shows only in the files of a build that
  succeeds, so the error comes at the end of such a build and names the Seq-ids and the sequence;
  an error of makeblastdb that comes later in the run (a mask of the same sequence past its end, a
  later sequence or input file that fails) is reported as makeblastdb reports it. The files are
  left as makeblastdb leaves them after an error of its build at the last sequence (none in
  version 4; `NAME.nos`, `NAME.not`, `NAME.ntf`, `NAME.nto` in version 5, `.p..` for proteins).
- A string whose length passes the end of the file and is 1 GiB to 2^47 bytes (`mask data that
  LOSAT does not read: ...`), before the build: makeblastdb reserves the length before it reads
  (`CIStreamBuffer::GetChars`, `util/strbuffer.cpp:496-523`), which succeeds or throws
  `std::bad_alloc` by the host's memory (the bound of line PS's DECISIONS 40, which `losat
  rpsbproc` uses).
- A value stated as `:Type` (NCBI's pointer syntax) of a type other than the one its place
  declares (`mask data that LOSAT does not read: ...`), before the build: makeblastdb reads it,
  or stops with "class not found", by the types its process has met before
  (`CClassTypeInfoBase::GetClassInfoByName`, `serial/classinfob.cpp:244-252`).
- An element of `masks`, `mix`, `equiv` or `packed-int` stated as `NULL` (`05 00` in binary)
  where makeblastdb dereferences it (`mask data that LOSAT does not read: ...`), before the
  build: in the `GetId` of a mask, after the masks before it and after any error of NCBI's reader
  in the rest of its object. NCBI's report of the null pointer prints the run's memory addresses
  (rule 3 of PD-LOSAT-NCBI-DEFECTS). A `NULL` that `GetId` does not reach (after two Seq-ids or a
  `feat`) gives NCBI's "Cannot get masked range Seq-id".

A string with a byte outside 0x20-0x7e is read as makeblastdb reads it: the byte is stored as `#`
and NCBI's error is written (`Error: (803.7) [makeblastdb] <the object's path>`, `Bad char [0x..]
in string at line N` or `at byte N`, and the string as NCBI's reader holds it;
`ReplaceVisibleChar`, `serial/objistr.cpp:1855-1897`); a `UTF8String` tag where a
`VisibleString` belongs is read with NCBI's warning (803.10), once per run. NCBI's pointer
syntax is read as makeblastdb reads it where NCBI reads an empty `CRef`: `NULL` leaves an
optional member (`fuzz-from`, `fuzz-to`, `fuzz`, `Seq-bond.b`, `PDB-seq-id.rel`) unset, `@n`
ends the run with NCBI's exception report (`invalid object index: NO_COLLECT defined`), and
`:Type value` of the declared type is that value. A mask file that
NCBI's reader rejects as damaged ends the run as makeblastdb ends it, with its exception report
(`Error: NCBI C++ Exception:` and the trace that names NCBI's build tree, exit status 255), and for
a string longer than the file: "BLAST ran out of memory" (exit status 4) from 2^47 bytes (no host
can allocate them) to 2^63 - 1 bytes, `Error: basic_string::_M_create` (exit status 255) past the
largest `std::string`, and NCBI's end-of-file report below 1 GiB.

Undefined behaviour of NCBI that is not reproduced:

- A volume file past 4 GiB is an error (makeblastdb's 32-bit offsets wrap).
- A version 5 database of more than 2^31 sequences is an error at the sequence that would get
  OID 2^31 (makeblastdb's `int` OID overflows). The number of sequences in the metadata file
  is an `Int8`, the type of the field, summed without NCBI's `int`.
- Either of these two errors, raised while the close writes the last sequence, ends the build
  as a duplicate id found there does (above): the files of the build are removed, then the
  writer's destructors write the `.pos`, `.pot` and `.pto` files and a `.ptf` of 16 KiB, and
  the last volume's ISAM and OID to GI files under its numbered name. A database whose first
  sequence is the one stopped leaves only the `.ptf` of 8 KiB that LMDB made in version 5 (the
  tax id writer holds nothing: "No tax info for any oid" ends the destructors) and nothing in
  version 4.
- A map size that the platform's `size_t` cannot hold (a 32-bit build) is refused with LMDB's
  error for a map it cannot make (`Cannot allocate memory`), not cut to 32 bits.
- A sequence longer than 2,147,483,647 residues is refused with NCBI's message, also where
  makeblastdb's 32-bit length wraps and it stores a short sequence. A nucleotide sequence of
  exactly 2,147,483,647 bases with ambiguity is written in full (makeblastdb's size
  computation overflows and it writes a shorter file).
- The test whether an input looks like FASTA does not read past the bytes it sampled.

The creation date is the local time, found without the C library. On Linux, macOS and WASI the
zone is the one that glibc finds for NCBI's Linux build (observed with glibc 2.39): with `TZ`
unset, the zone file `/etc/localtime`; otherwise `TZ` less one leading `:` (UTC where nothing is
left) names the zone file itself when it starts with `/`, else the file of that name under
`TZDIR` (`/usr/share/zoneinfo` when `TZDIR` is unset or empty), and where that is no zone file it
is a POSIX rule. Under WASI the environment and the files are those that the host gives the
module. The crate jiff reads the zone file and the rule (line LIB, 2026-11-20, on the
maintainer's 「TZ=A」 of that day); LOSAT reads the leap-second records of a zone file itself, so
that tzdata's `right/` zones give glibc's time, second 60 included. On Windows the zone comes
from `TZ` and the system zone as Microsoft's C runtime (UCRT) reads them for NCBI's Windows
build. A `TZ` offset beyond 596,523 hours overflows UCRT's 32-bit arithmetic (undefined in C++);
LOSAT takes the mathematical value. When no zone can be read, the date is UTC; so it is on
Windows where UCRT's `localtime` fails (a `TZ` offset that puts the local time before 1969-12-31
12:00 or after 3001-01-19 21:59:59, where NCBI's `CTime` throws).

Where glibc reads a `TZ` that POSIX does not write, or a zone file that jiff rejects, the date of
NCBI's Linux build is glibc's reading and LOSAT's is UTC (8,588 rules and every zone of tzdata
2026e in four builds compared with glibc, 1970 to 2100; the forms below are in
LOSAT/tests/unit/local_time_glibc_vectors.tsv, which CI compares with LOSAT and with glibc). The
local time at 2026-01-01 00:00 UTC and at 2026-07-25 17:20 UTC, as glibc gives it and as LOSAT
does:

| `TZ` | glibc | LOSAT |
| --- | --- | --- |
| `XXX+25` (hours past 24) | UTC-24:00, UTC-24:00 | UTC |
| `XXX-168` | UTC+24:00, UTC+24:00 | UTC |
| `XXX99999` | UTC-24:00, UTC-24:00 | UTC |
| `XXX5:70` (minutes past 59) | UTC-05:59, UTC-05:59 | UTC |
| `EST5EDT,M13.2.0,M11.1.0` (a month past 12) | UTC-04:00, UTC-05:00 | UTC |
| `EST5EDT,M3.6.0,M11.1.0` (a week past 5) | UTC-04:00, UTC-04:00 | UTC |
| `EST5EDT,M3.2.7,M11.1.0` (a weekday past 6) | UTC-04:00, UTC-04:00 | UTC |
| `EST5EDT,J0,J300` (a Julian day 0) | UTC-04:00, UTC-05:00 | UTC |
| `EST5EDT,366,300` (a day past 365) | UTC-04:00, UTC-04:00 | UTC |
| `EST5EDT,M3.2.0/168,M11.1.0` (a rule time of 168 hours) | UTC-05:00, UTC-04:00 | UTC |
| `AAA--8`, `JST--9` (two signs) | UTC+24:00, UTC+24:00 | UTC |
| `JST+-9` | UTC-24:00, UTC-24:00 | UTC |
| `EST5EDT,M3.2.0` (one date of the rule) | UTC-05:00, UTC-04:00 | UTC |
| `EST5EDT,` (an empty rule) | UTC-05:00, UTC-04:00 | UTC |
| `EST5EDT,M3,M11.1.0`, `EST5EDT,M3.2,M11.1.0`, `EST5EDT,,M11.1.0` (a date cut short) | UTC-04:00, UTC-04:00 | UTC |
| `EST5EDT,m3.2.0,m11.1.0` (a lowercase `m`) | UTC-04:00, UTC-04:00 | UTC |
| `EST5EDT4`, `XXX5YYY`, `XXX5YYY4` (daylight time without a rule: glibc takes the rules of the zone file `posixrules` under `TZDIR`, New York's on Debian and Ubuntu, or without it the United States' since 2007) | UTC-05:00, UTC-04:00 | UTC |
| `XXX-9YYY`, `JST-9junk` | UTC+09:00, UTC+10:00 | UTC |
| `XXX1:2` (minutes of one digit) | UTC-01:02, UTC-01:02 | UTC |
| `XXX-5:3` | UTC+05:03, UTC+05:03 | UTC |
| `XXX5YYY4:3,M3.2.0,M11.1.0` | UTC-05:00, UTC-04:03 | UTC |
| `JST-9 `, `JST-9:` (text after the rule) | UTC+09:00, UTC+09:00 | UTC |
| `JST-9,` | UTC+09:00, UTC+00:00 | UTC |
| `EST5EDT,M3.2.0,M11.1.0 `, `EST5EDT,M3.2.0,M11.1.0,M12.1.0` | UTC-05:00, UTC-04:00 | UTC |
| `America/Ojinaga` of a zone file built with `zic -b slim` (its footer disagrees with its last transition; the distributions' fat file is read) | UTC-06:00, UTC-05:00 | UTC |

A rule whose change falls in another year in UTC than in local time (a valid POSIX rule) gives
LOSAT another local time than glibc's for up to a day at the turn of the year: jiff keeps the
change in its year in UTC, and glibc puts it where the rule says (231 of 5,000 such rules
generated differ; no zone of tzdata 2026e does). The local time at 2025-12-31 23:59:59 UTC:

| `TZ` | glibc | LOSAT |
| --- | --- | --- |
| `AAA-9BBB,0,365`, `AAA-9BBB,0/0,J365/25`, `AAA-9BBB,J60,365` | UTC+09:00 | UTC+10:00 |
| `AAA1BBB,J60,365/-3` | UTC-01:00 | UTC+00:00 |
| `XXX-24YYY,0/0,J365/25` | UTC+24:00 | UTC+25:00 |

A zone file that breaks RFC 8536 (transitions or leap seconds out of order, a leap-second step
other than one second, a footer that disagrees with the only local time type of a file without
transitions) can give another time than glibc's, or UTC; zic writes no such file.

A database of 1,000 volumes or more is not in FASTA order when it is read, with makeblastdb
and with `losat makeblastdb` alike: NCBI opens volumes in the sorted order of their names, and
the names makeblastdb gives then (`db.0010`, `db.00100`, `db.00101`, ...) do not sort by
number. `losat makeblastdb` writes the same names and LOSAT's reader opens them in NCBI's order.

## What BLASTP `-db` of stage DB-3 does and rejects

`losat blastp -db NAME` searches a protein database of the reader above and gives the output
of NCBI BLAST+ 2.17.0 `blastp -db NAME -num_threads 1` for outfmt 0, 6 and 7, the same bytes
at every `-num_threads` (P1, P1.1). The `-subject` path is unchanged; both paths run the same
search code with a different subject source.

- Statistics: the number of sequences and the total length of the database, and 10 as the
  shortest subject length behind the hit-saving cutoff (a database always answers 10;
  `-subject` answers the shortest subject of the set, at least 10). This is why `-db` and
  `-subject` of the same sequences can differ near the cutoff, in NCBI and in LOSAT alike.
- `-dbsize N` (NCBI's `Int8`) for `-db` and `-subject`, with NCBI's meaning: any value but 0
  replaces the database length in the Gumbel parameters, a value above 0 also in the
  effective lengths.
- The database is read in blocks of about 1 MiB. With one thread the blocks are searched in
  order; with more, each worker starts a block from the state the one-thread order has
  there, and the results are merged in database order. During the search the memory taken
  by the database is its index (12 bytes per sequence) and the blocks in work; afterwards
  only the sequences of the subjects in the hit lists and the titles of the reported ones
  are read.
- Residues: a stored `O` is read as `X` for the search and as `O` for the composition-based
  redo (NCBI compares the stored bytes there); `U`, `B`, `Z`, `J`, `*` as stored.
- Reports: the `Database:` lines of outfmt 0 (title, counts, `Posted date`), `>` and the
  stored title, subject ids made from the title as NCBI makes them (the general id
  `gnl|BL_ORD_ID|n`, `BL_ORD_ID:n` for the accessions, with n the sequence's OID, only for a
  title that gives no word), `# Database:` of outfmt 7 with the `-db` value as given.
- Errors: NCBI's texts and exit statuses for a database that is not found or cannot be
  read, at the point where NCBI reports them (a header file is first read when the reports
  of a query batch are made). A header record that NCBI's object stream fails on (a damaged
  record, or one that another writer encoded so: a taxid of another type or outside `int`, a
  tag or a length NCBI's reader does not take, a record cut short) gives NCBI's
  `CSerialException`, made where NCBI throws it (position, message, code, the frames of NCBI's
  object stack and the backlog `error while reading <member>` of each member with a set flag),
  as each program reports it: the search programs `Error: <what()>` with exit status 255, after
  the outfmt 0 prolog and the reports of the batches before the batch whose reports read the
  record (a record cut inside its first tag or length, the opening bytes of its set of
  definition lines, throws NCBI's `CEofException`, which `CATCH_ALL` takes as a `CIOException`
  and does not post: the run ends there with exit status 0, fix 1 of item T5d; a record cut
  after them, before its first definition line, gives the exception of NCBI's `HandleEOF` with
  exit status 255 as any other cut record, fix 2); blastdbcmd `Error: [blastdbcmd] <GetMsg()>` with exit status 1 after the records
  before; segmasker `<what()>` (exit 1); dustmasker NCBI's report of each exception and its
  `Application's execution failed` line (exit 3). A record that no report reads fails nothing
  (NCBI reads no header when it opens a database). Line RG item T5d. Two encodings that
  makeblastdb does not write are read as NCBI reads them (fix 1 of item T5d): a UTF8String where
  a VisibleString is expected, as the VisibleString (`SERIAL_READ_ANY_VISIBLESTRING_TAG`), with
  NCBI's warning `Warning: (803.10) [<program>] CObjectIStreamAsnBinary: UTF8String data for
  VisibleString member <member>, ASN.1 specification may need an update` once per run
  (`ERR_POST_X_ONCE`) where NCBI first reads such a record, whichever read it is: in the search
  programs where a batch's search or Seq-aligns read it (after the batch's reader messages,
  before its prefetch and its reports; a message of the BLAST archive). Under a negative seqid
  list with ids NCBI's traceback and composition-based redo read the header of each subject they
  fetch (the redo only those its early termination test leaves), also of one that no report shows
  (fix 2 of item T5d); a seqid list's Seq-aligns read every final hit, also one that gets no
  Seq-align; the version check of an accession reads its candidates. NCBI's scope reads the
  header of the sequence that the lookup of a reported sequence's Seq-align id gives first, where
  that is another one (the data loader's `x_GetOid`, fix 3 of item T5d): in every format that
  reads the subject's Bioseq and where the batch's prefetch loads it (also a table of fields that
  read nothing of it, and outfmt 8 with `PRE_FETCH_SEQS_LIMIT=2147483647`), not in the archive nor
  the Seq-annots without the prefetch; another holder after the subject is not read. A sequence
  that differs from the subject only by such a string is the subject to every report. The
  warning comes at the load of that id (fix 4 of item T5d): in the batch's prefetch after the
  taxonomy database's texts of the loads before it in NCBI's order of the Seq-id handles (by packed
  number), before a failed lookup's `Error pre-fetching sequence data` and a stopping prefetch's
  `BLAST Database error`; where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off, inside the
  report at its first read of the subject (outfmt 0's description line, a tabular row, a BLAST XML,
  XML2 or JSON `Hit`; psiblast's PSSM after the round's report where the report reads nothing of
  it), and in the prefetch of a batch that the limit leaves on with a format that reads nothing of
  the subject. Refused there: two such loads whose warnings name different members (two of one
  member write the same bytes in either order, fix 5), or that load and one that asks the taxonomy
  database where the batch writes a text of it, whose order NCBI takes from the addresses of its
  Seq-id infos (one packed number, or neither packed: NCBI 2.17.0 wrote the warning before the texts of
  opening in one run and after them in the next, blastx probe `gt.bad_magic`, fix 4) or that LOSAT
  does not know (a batch of more than 200 Seq-align ids, as line DBID loadord's stops). A
  prefetch whose loads end at a lookup that throws on its id (`Error pre-fetching sequence data`)
  before that load does not make it, so it writes no warning (fix 5 of item T5d): the warning comes
  at the run's next load of the id (a later batch's prefetch) or not at all; a report that reads the
  subjects reads the one whose lookup throws too, which LOSAT refuses. In blastdbcmd before the
  record's text, and before `Entry not found` where the version check of `-entry` or
  `-entry_batch` read it; segmasker and dustmasker write none. A tag in the long form is read
  with no message (`PeekLongTag`), or fails as NCBI's read does.
- `-db` with `-subject` or `-subject_loc` is a usage error: NCBI's usage and error lines, exit
  status 1, as in NCBI (since line RG item K4i-2b, 2026-11-23; before, LOSAT's parser text and exit
  status 2). LOSAT Web rejects `-db` until stage DB-7.

Rejected explicitly (NCBI runs these, or crashes, or LOSAT cannot reproduce what it does):

| Rejected | Why |
| --- | --- |
| `-searchsp` | NCBI then uses a different statistics path (`OLD_FSC`, no Gumbel block); a stage of its own |
| A database sequence with a gap residue (`-`, stored by makeblastdb), a byte above the protein alphabet, or a sequence not between two NUL bytes | NCBI's `-db` never aligns a stored gap and its formatter counts it as a gap; not reproduced |
| A database sequence without parsed Seq-ids with more or fewer than one definition line | makeblastdb without `-parse_seqids` writes one; a sequence with parsed Seq-ids may have several (version 1.3) |
| A reported subject whose title NCBI's cleaner reads past (titles like `, ,`: NCBI crashes) | tabular formats and outfmt 0 |
| A shown subject with an empty title, or a protein title ending in `]` (line RG item T5c), and a tax id (makeblastdb `-taxid`) whose scientific name NCBI reads past its mapping of `taxdb.btd`: a record of an empty `taxdb.btd`, or a record that reaches past the mapping's last 4096-byte page without its third tab inside it (line RG items T5 and T5b; outfmt 0, 18 and the tabular subject ids, all six programs) | An empty `taxdb.btd` NCBI maps as a null pointer and crashes on; past the mapping of another it reads memory whose bytes, or whether it crashes, depend on its memory map (PD-LOSAT-NCBI-DEFECTS: what the output writes from past the mapping is refused, rule 3; the rest is written as NCBI writes it where it completes, approved exception 14). A title made from a name the taxonomy database gives is written as NCBI makes it (line RG item T5; line DBID H3 refused it in outfmt 0, line DBID ordname in a nucleotide's tabular subject id) |
| Any other lookup in the taxonomy database that reads a record of an empty `taxdb.btd` (the prefetch's lookups of a reported subject's tax ids, psiblast's PSSM, the name fields), and a name a report writes whose bytes NCBI reads past its mapping of a damaged `taxdb.btd`: a super kingdom (`sskingdom`, `sskingdoms`) of a record past the mapping, any name of a record without its third tab inside the mapping (`sscinames` and the other name fields, the `sciname` of outfmt 13 to 16, outfmt 18) (line RG item T5b) | As the row above. A lookup whose names no report writes (the prefetch's, psiblast's PSSM) and the names inside the mapping are written as NCBI writes them where it completes |
| A volume of the database that holds a tax id whose record begins past the end of a damaged `taxdb.btd` (NCBI writes `Offset error at end of taxdb file.` at each lookup of it that asks the taxonomy database) and another positive tax id of its slot of NCBI's 256-slot taxonomy cache of a volume (`CSeqDBVol::m_TaxCache`, the slot `taxid & 255`), in any definition line of the volume (all six programs): a query batch's prefetch where the order of its loads decides how often NCBI writes the offset error, and the run's other fetches. NCBI's prefetch loads a batch's sequences in the order of their Seq-id handles (`CSeq_id_Handle::operator<`): the ids that its mapper packs (GIs, accessions and the others of the prefetch's stop in PD-LOSAT-NCBI-DEFECTS) by their packed value less 1 as an unsigned number (`gnl|BL_ORD_ID|0` after every positive one), then the others, and the ids of one packed value, by the addresses of their handles' records (`CSeq_id_Info`); a batch of more than 200 ids it loads in pieces of 100. LOSAT follows that order (line NP item O4b2-3) and writes NCBI's texts where the loads of one key give the same count in any order (rows `db.nz.bp.c150.pf.*`, `db.nz.bp.c150.ung.lim3.*` and `db.tm.bp.cachegi.*` of BLASTP (gi ids: six offset errors in the order of their gi numbers, one in that of the OIDs), `db.nz.cq.*.c150.b1.x2.*` of BLASTN, TBLASTN and TBLASTX, `db.nz.cq.x.bn.rep.b1.6`, `db.nz.cq.x.bx.reppo.b1.qx.0`, `db.tm.psi.c150.qbig.6`). It takes all of a volume's loads as one key where the batch has more than 200 ids or a Seq-align id of two sequences, and refuses: a batch whose loads of one key give another count in another order (two loads of that tax id and one of the other among ids that its mapper does not pack: `db.tm.bp.cache.6.reject`, `db.tm.bx.cache.6.reject`); a later load whose texts depend on what such a batch left in the slot; a load of a sequence by a Seq-align id that loaded it before (NCBI's data loader names a loaded sequence by its ordinal id with the Seq-id handle that asked for it), after NCBI resets its scope (`ResetScopeHistory`, between query batches in formats other than BLAST XML outfmt 5, and between queries in psiblast) on a database of more than 10 sequences (or without a pinning prefetch on a database of more than 20 sequences or with `OBJMGR_*`), or by another Seq-align id that the Bioseq of a blob still live in NCBI's scope holds (outside a report read whose exact Seq-id handle is held by a kept live-scope Bioseq, or the prefetch that made it), where NCBI may keep it or load it again (where making the load and not making it give other texts); a load by a Seq-align id whose Seq-id handle LOSAT does not know (a text id with neither accession nor name, for which NCBI makes a new handle at each read, or a giim id); a batch with a Seq-align id whose lookup fails and a load of the slot whose place to the failed lookup LOSAT does not know (a Seq-id handle of the failed one's packed value, or neither packed; a Seq-align id of two sequences; a batch of more than 200 ids or with a handle that LOSAT does not know), which NCBI may make before the failed lookup or not at all, where making it and not making it give other texts; a load of the slot before the stop of a prefetch whose place among the loads LOSAT does not know, and one that a stopping id whose lookup NCBI may not make is given instead; a report's read of the slot on a database of more than 20 sequences, or where an environment variable of NCBI's object manager is set (`OBJMGR_*`, whose sizes decide what it keeps; such a variable also makes every read again a load that NCBI may make); a BLAST XML report that NCBI's XML stream cuts before a read of the slot, leaving it out; and a secondary XML2 or JSON report after a non-XML2 first format in a multi-format library call on a colliding volume. Since line NP item O4b2-5 LOSAT also follows the slot through the reports' reads in a batch whose prefetch is off (`PRE_FETCH_SEQS_LIMIT`; rows `db.nz.bp.c150.lim0.*`, `db.nz.bp.c150.lim1.nident.0`, `db.tm.bp.cachegi.lim0.*`, `db.nz.*.c150.off0.*`, `db.nz.bx.c150.lim1.nident.6` and `*.lim1.*.0` of BLASTN, TBLASTN, BLASTX and TBLASTX, and PSIBLAST's PSSM reads of a later query, `db.id.prd3.tc.lastoff.pf0.q2` and `db.id.prd4.m2.lastoff.pf0.q2`) and through the loads that surely come before a prefetch's stop, in NCBI's order of a prefetch, and takes a read again of a sequence of a database of at most 10 sequences by the Seq-align id that read it before as no load (`db.nz.bp.c150.pf.b1.0`, `db.nz.bp.c150.lim0.b1.0`). Since line NP item O4b2-6 the XML2 and JSON reports follow the slot by the same rules, over the `Hit`s that their writers make (a copy of one record that follows another reads nothing; where the prefetch is off, a report reads up to its `Hit` whose lookup throws), with the texts where NCBI writes them: a batch's prefetch where the batch begins, after the posts before its first query's report (the texts of opening first), a read where the prefetch is off at the read's `Hit`, before the `Hit`'s own lookups (rows `db.nz.bp.c150.pf.{15,16,b1.16}`, `db.nz.bp.c150.lim0.{15,b1.16}`, `db.nz.cq.x.bp.c150.off0.{15,16}`, `db.nz.cq.x.bn.rep.off0.16`, `db.nz.cq.x.bn.repo.{b1.15,qd.16}`, `db.nz.cq.x.tn.rep.off0.15`, `db.nz.cq.x.b{p,x}.reppo.b1.qx.16`, `db.nz.cq.x.tx.rep.off0.b1.16`). Since line NP item O4b2-7 a load of a sequence by another Seq-align id than its loads before, which no Bioseq of those loads holds, is a load of its own that follows the slot, and in one prefetch always (it asks its scope only for ids resolved before it; rows `db.tm.hd.*`: the second volume's Seq-align id that NCBI's data loader gives as the first volume's sequence looks that sequence's slot up again after its own load, also where that id is in the same definition line, `db.tm.hd.g.*`), and the loads of one key that NCBI may make or not are counted both ways, refused only where the two give other texts. Since line NP item O4b2-8 LOSAT tracks NCBI's per-batch scope lifetime (`CBlastFormat::ResetScopeHistory` at the end of each query batch in formats other than BLAST XML outfmt 5, and between queries in psiblast while keeping scope across rounds of a query): blobs whose scope was reset no longer block another Seq-align id's fresh load (`db.tm.hd.g2.bp.qe.q2b1.6`), outfmt 5 keeps its scope across batches so succeeding prefetches stay pinned (`SetLinkUsedTSE`) regardless of database size and non-prefetched reads stay live up to 20 sequences (`db.tm.bp.cachegi.q2b1.{5,lim0.5}`, `db.tm.bp.gi21.q2b1.5`, `db.tm.psi.cachegi.{it2.6,q2.5}`, `db.tm.{bx.cachegi,bn.cachengi,tn.cachengi,tx.cachengi}.q2b1.5`), report reads matching an exact handle of a kept live-scope Bioseq load nothing, and `prefetch_stop` sure loads on colliding databases (`db.tm.bp.stopgi.{6,0,16,q2b1.6}`) plus the 11-20, >20, and BLAST XML stream-cut cases (`db.tm.bp.cachegi12a.q2b1.lim0.6`, `db.tm.bp.gi21.*.reject`, `db.nz.b{p,x}.c150.lim0.xmlcut.5.reject`) are covered. Since line NP item O4b2-9 XML2 and JSON (outfmt 13 to 16) share the unified batch-prefetch path with the other formats, placing opening and volume-lookup texts at batch start for non-colliding and slot-silent damaged-taxdb runs (`db.o4b2_9.{bp,bn}.r1.*`), emitting `Error pre-fetching sequence data` warnings on failing lookups (`db.o4b2_9.bp.r2.*`), interleaving UTF8String (803.10) warnings with taxonomy opening and colliding slot texts in NCBI's Seq-id handle order while refusing packed-key address-order ties (`db.o4b2_9.bp.u8.*`), and following colliding slot state across batches before a stopping batch (`db.o4b2_9.bp.{prestop1,prestop2,stopgi.q2b1}.{15,16}`), while truncating one-by-one reads at a throwing `Hit` and collapsing adjacent joined copies (`db.o4b2_9.bp.h2.lim0.{15,16}`). Since line NP item O4b2-10 a batch's prefetch with a Seq-align id whose lookup throws on the id follows NCBI's order of its Seq-id handles in all six programs: the loads before the first such lookup write their texts and take their steps through the slot, the loads after it are not made (a report that reads such a sequence loads it then), and a load whose place to it LOSAT does not know is refused only where it writes a text, where it reads a header with a UTF8String and no load that surely comes before it read one (NCBI writes that warning once in a run), or where it takes the slot of a colliding volume and making it and not making it give other texts (rows `db.o4b2_10.*` on `ids/ntagx_t5`, `ids/ntagc_t5` (a colliding volume), `ids/ptagx_t5`, `ids/ptagc_t5`, `ids/ntag_t5` and `ids/pshare_t5`, `db.id.pfw.bn.tag.off.6` and `db.id.pfw.bp.tag.bad.6`; `db.o4b2_10.bp.pshare.q2.mts2.6qe.b1.reject`: `lcl|555` and the text tag, both unpacked, where NCBI wrote 0, 1, 1 and 0 offset errors in four heap layouts). Since its second round the loads that such a prefetch made go to NCBI's data source's cache of 10, not to the scope: its map releases them the largest key first, so the loads of the 10 smallest keys stay there and a report's read of another loads it again and looks its slot up again (rows `db.o4b2_10.bp.pslot{11,12,10,c,r,m}.*`); since its third round every such load counts among those 10, whatever its taxids, where the scope holds no read of the run since its reset (where it may, as in BLAST XML, which keeps its scope, a load that the scope resolves releases nothing, and LOSAT refuses where that decides the texts: `db.o4b2_10.bp.pslotscope.b1.5.reject`; `db.o4b2_10.bp.pslotmix.*`, `db.o4b2_10.bp.pslotmix12.*`, and `db.o4b2_10.bp.pslotscope.b1.0` and `.b1.6sci`, whose scope a report reset), and a load of the failing lookup's key that NCBI may make or not, of a sequence that NCBI holds nowhere before, is loaded again at its read either way after 10 such loads (`db.o4b2_10.bp.pslottie.*`, refused before); a load of the slot that NCBI may make before the failure or not at all is followed both ways (`db.o4b2_10.bp.ptie.*`, `db.o4b2_10.bp.ptie150.*`, and `db.o4b2_9.bp.r2.mix.off.16`, a refusal before); the 200 ids of NCBI's map count every Seq-align id of the batch, also one that a query resolves; BLASTN's Seq-aligns on the minus strand read no length where a lookup of the batch fails (`db.o4b2_10.bn.rc.*`, refused before); and rows on two volumes (`db.o4b2_10.bp.pv{ab,cd}.*`) and of psiblast's later rounds (`db.o4b2_10.psi.it{2,3}.*`) | Where the addresses decide, NCBI's count changes with its memory layout (a relative or absolute `BLASTDB`, the length of an unused environment variable: 1, 1, 2 and 1 offset errors in four layouts of `db.tm.bp.cache.6.reject`). In a batch of more than 200 ids NCBI's count changes from run to run too (a volume of 573 proteins with gi ids whose taxids alternate, 200 queries: 563 and 566 offset errors from one command line in two runs). NCBI's object manager keeps a sequence that it loaded while a handle holds it, and a batch's prefetch links what it loads to a Bioseq it adds to the scope, which holds them to the scope's reset after the batch; released, a sequence waits in its scope's queue of 10, then in its data source's cache of 10, and is dropped after that, so its next read loads it again and looks the slot up again (`OBJMGR_SCOPE_AUTORELEASE_SIZE` and `OBJMGR_BLOB_CACHE`, 10 each by default), and its data loader names a loaded sequence by its ordinal id with the Seq-id handle that asked for it: a database of at most 10 sequences never loses one, the reads of one of at most 20 lose none, and above that a read that LOSAT does not model (a report's second read of a subject) can write the offset error again. The XML2 and JSON lookups come in an order that NCBI fixes (the same bytes in the four layouts of the probe of line NP item O4b2-3, and of the 12 rows of line NP item O4b2-6); where a batch's prefetch is ordered by addresses, they are refused as the other formats are (`db.nz.cq.x.bp.c150.15.reject`) |
| A reported title with a raw byte outside ' '..'~' in the header file, in the cases narrowed by `PD-LOSAT-NCBI-DEFECTS` 1.120 and 1.127 (`psiblast`, `rpsblast` and `rpstblastn` except a query batch that ends at a header NCBI cannot read; `-outfmt` 0 to 5, 13 to 16 and 18; with `-outfmt` 6, 7, 10, 17 and 20, or 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647`, a subject whose blob NCBI may still hold or not from an earlier batch or whose header cache slot such a subject takes (not determined, rule 3), and, determined but not modelled yet (line RG-KB item 9b-3), a prefetch that may be off, parsed Seq-ids, `OBJMGR_*` settings, a run whose prefetch's loads write texts of the taxonomy database; a search restricted by a seqid, GI or taxid list, an alias file with a list or an OID mask, or a volume with an OID mask type; a query batch whose headers also give NCBI's warning of a UTF8String; line RG-KB items 9(b) and 9b-2) | NCBI prints `#` and an error message (`Error: (803.7)`) each time it reads the definition line; in the other runs of the five search programs with `-outfmt` 8, 9, 11 and 12, 6, 7, 10, 17 and 20 (the messages of the Seq-aligns, of the prefetch's loads and of the table rows that read the definition line) (and in any format of all eight search programs for a query batch that ends at a header NCBI cannot read), LOSAT writes NCBI's messages and output |
| `-dbsize` whose effective search space overflows NCBI's `Int8` | NCBI's arithmetic overflows |
| A database of more letters than the target's `usize` (32-bit targets) | |
| A database that cannot be searched, when the first query batch holds invalid queries only | NCBI would report that batch first; not reproduced |
| A registry value `[BLAST] BLASTDB` that LOSAT does not read as NCBI does | |

Known differences that remain: the reader's (section above), and a version 5 database
copied under another name, which NCBI refuses (`Input db vol does not match lmdb vol`) and
LOSAT searches.

## What BLASTN `-db` of stage DB-4 does and rejects

`losat blastn -db NAME` (megablast, blastn, dc-megablast) searches a nucleotide database of the
reader above and gives the output of NCBI BLAST+ 2.17.0 `blastn -db NAME -num_threads 1` for
outfmt 0, 6 and 7, the same bytes at every `-num_threads` (P1, P1.1). The `-subject` path is
unchanged; both paths run the same per-subject search, the same traceback and the same report
writers with a different subject source. The database opener, the reader's error texts and the
title rules are those of BLASTP's `-db`.

- Bases: the search reads the stored 2-bit bases, and the alignments read the letters rebuilt
  from the 2-bit bases and the stored ambiguity codes, as NCBI's search and traceback do.
  makeblastdb (and `losat makeblastdb`) draw the 2-bit base of each ambiguous letter from `CRandom`
  seeded with the sequence length, as `-subject` draws it at search time (writedb_convert.cpp,
  blast_objmgr_tools.cpp), so `-db` and `-subject` of the same FASTA search the same bases, in
  NCBI and in LOSAT alike, except where the FASTA has hyphens (makeblastdb stores them as `N`,
  `-subject` drops them), with `-subject_loc` (the bases are drawn over the range), or for a
  database written otherwise.
- Statistics: the number of sequences and the total length of the database, kept in 64 bits as
  NCBI keeps them (the 2^31 - 1 limit of `-subject`, NCBI's `int` there, does not apply). BLASTN
  has no shortest-subject cutoff (BLASTP's 10 does not apply). `-dbsize N` (NCBI's `Int8`), for
  `-db` and `-subject`: a value above 0 replaces the database length in the length adjustment
  and the effective search space; the counts, the reports and the query batching keep the real
  numbers.
- The database is read in blocks of about 1 MiB of 2-bit data. With more threads the blocks are
  searched in parallel and their results merged in database order. After the preliminary stage
  only the sequences of the hit lists are read again, one at a time, for the alignments; the
  definition lines are read for the final hit lists of each query batch, where NCBI reads them.
  A sequence longer than 5,000,000 bases is searched in NCBI's overlapping chunks, as with
  `-subject`. The alignments decode the whole sequence in hand (NCBI decodes the ranges it
  needs and refetches the whole sequence when an alignment leaves them, with the same result),
  and outfmt 0 holds the letters of the reported sequences until the reports are written. On
  2026-11-14 line NP M4a-2a searched BLASTN with `-db_soft_mask` (section "BLASTN database masks")
  and read `-db_hard_mask` up to the search, which M4a-2b searches the same day.
- Reports: the `Database:` lines of outfmt 0 (title, counts, `Posted date`), `>` and the stored
  title, subject ids made from the title as for BLASTP, `# Database:` of outfmt 7 with the `-db`
  value as given. A title of punctuation that NCBI's cleaner reads past (`, ,`; NCBI BLAST+
  2.17.0 crashes on it in every format) is written in outfmt 0 as BLASTN writes it for
  `-subject` (PD-LOSAT-NCBI-DEFECTS) and rejected in the tabular formats, as for BLASTP.
- Errors: NCBI's texts and exit statuses for a database that is not found or cannot be read,
  at the point where NCBI reports them, as for BLASTP. A database without sequences fails in the
  first searched query batch, as in NCBI.
- `-db` with `-subject` or `-subject_loc` is a usage error, as for BLASTP. LOSAT Web rejects
  `-db` until stage DB-7.

Rejected explicitly, besides what BLASTP `-db` and the reader reject (the same rules apply):

| Rejected | Why |
| --- | --- |
| `-searchsp`, `-db_soft_mask`, `-db_hard_mask` | a stage of its own; masks are not read (P4). Since 2026-11-14 `-db_soft_mask` and `-db_hard_mask` are searched (line NP, M4a-2a and M4a-2b; section "BLASTN database masks") |
| In the tabular formats that print a subject id (the default fields, `sseqid`, `sacc`, `saccver`, `sallseqid`, `sallacc`), a reported subject with an empty title and a tax id that the taxonomy database names (line DBID ordname) | NCBI makes the title from the tax id, as in outfmt 0 (line DBID H3, the row of BLASTP), and prints its first word, the organism's (`Homo` for 9606), in the place of the id; the other fields show the stored title (`N/A`) and the taxonomy as LOSAT does. A title of spaces, and an empty title without such a tax id, give no word: NCBI prints the general id, `gnl|BL_ORD_ID|n` for `sseqid` and `BL_ORD_ID:n` for `sacc` and `saccver` with n the sequence's OID, as LOSAT does (for a protein an empty title gives `unnamed`) |
| `-dbsize` whose effective search space overflows NCBI's `Int8` | NCBI's arithmetic overflows |
| A database of more letters than the target's `usize` (32-bit targets) | |

## What TBLASTN `-db` of stage DB-5 does and rejects

`losat tblastn -db NAME` searches a nucleotide database of the reader above and gives the output
of NCBI BLAST+ 2.17.0 `tblastn -db NAME -num_threads 1` for outfmt 0, 6 and 7, the same bytes at
every `-num_threads` (P1, P1.1). The `-subject` path is unchanged; both paths run the same
search, the same traceback and composition-based redo, and the same report writers. The database
opener, the reader's error texts and the title rules are those of BLASTP's and BLASTN's `-db`.

- Bases: the preliminary stage translates the stored 2-bit bases, and the alignments, the
  composition-based redo and the reports read the letters with the stored ambiguity codes, as
  NCBI's search and traceback do. makeblastdb (and `losat makeblastdb`) draw the 2-bit base of each
  ambiguous letter from `CRandom` seeded with the sequence length, as `-subject` draws it at
  search time, so the bases of `-db` and `-subject` differ only where the FASTA has hyphens
  (stored as `N`, dropped by `-subject`), with `-subject_loc`, or for a database written
  otherwise.
- Statistics: the number of sequences and the total length of the database, in 64 bits as NCBI
  keeps them. The shortest subject length behind the hit-saving cutoff is the database's 10
  (3 residues), not the shortest sequence, so `-db` keeps the weak hits NCBI `-db` keeps.
  `-dbsize N` (NCBI's `Int8`), for `-db` and `-subject`: a value above 0 replaces the database
  length in the length adjustment and the effective search space, and a value other than 0 in
  the database length of the E-value statistics (NCBI's Gumbel block); the counts and the
  reports keep the real numbers. `-db_gencode` translates the database as it translates
  `-subject`.
- Memory: the database's sequences are read once, in blocks of about 1 MiB, and kept packed
  (two bits per base and the runs of ambiguity codes, about the size of the `.nsq` file); NCBI
  maps the same file and its search touches all of it. The search, its parallel split and its
  translations then run as for `-subject` with packed subjects. A database larger than memory
  is not searched in blocks yet.
- Reports: the definition lines are read for the final hit lists of each query batch, where
  NCBI reads them; the `Database:` lines of outfmt 0, `>` and the stored title, subject ids made
  from the title as for BLASTN, `# Database:` of outfmt 7 with the `-db` value as given.
- Errors: NCBI's texts and exit statuses for a database that is not found or cannot be read, at
  the point where NCBI reports them. A database without sequences or a sequence file that cannot
  be read fails the first query batch with a valid context (NCBI does not search a batch without
  one), after the reports of the batches before; a run without such a batch succeeds. Damaged
  files are LOSAT's refusal, with no report: ambiguity data that the reader cannot read (NCBI
  reads them only for the sequences it aligns), a sequence file shorter than its index (NCBI
  reads past its end), and a definition line that cannot be decoded (NCBI throws after the
  reports of the batches before).
- `-db` with `-subject` or `-subject_loc` is a usage error, as for BLASTP. LOSAT Web rejects
  `-db` until stage DB-7.

Rejected explicitly, besides what BLASTN `-db` and the reader reject (the same rules apply):
`-searchsp`, `-db_soft_mask`, `-db_hard_mask`, and the options TBLASTN rejects with `-subject`.

## BLASTN database masks (line NP, M4a-2a and M4a-2b, 2026-11-14; M4a-2c, 2026-11-15)

`losat blastn -db NAME -db_soft_mask ID` searches with the masks of the filtering algorithm `ID`
that `makeblastdb -mask_data` stored in the database, as NCBI BLAST+ 2.17.0 does, byte for byte:

- The option is an id when it is a 32-bit integer in NCBI's form (`+11` and `011` are 11) and a
  name otherwise (looked up as stored, letter case included: the description of `-mask_id`'s name).
  An id above 0 that the database does not have is NCBI's options error when the database is
  opened (exit 1); after the test for an empty query, a name that the database does not have and
  the id -1 give NCBI's warning and a search without masks, and any other id must be known to the
  database of each name of the `-db` list, else NCBI's database error with that database's
  algorithms (exit 2).
- Every database sequence is searched as a masked subject: words inside its masks start no hit;
  the extensions and the alignments read the real letters. The ranges between the masks are
  taken in the order the masks are stored, so unsorted or overlapping masks give NCBI's results.
- Outfmt 0 shows `Masked using:` after the database title (the first name's algorithm) and in
  each database block of the end (each name's), and the database letters inside a mask that meets
  an HSP in lower case. The BLAST archive and the search strategy write
  `DbFilteringAlgorithmKey` for a name or `DbFilteringAlgorithmId`, then `SubjectMaskingType`.
  The other formats change only through the results.

`losat blastn -db NAME -db_hard_mask ID` (M4a-2b) reads the option the same way and searches
with the masked letters cut out, as NCBI BLAST+ 2.17.0 does, byte for byte:

- The preliminary search takes each database sequence in chunks, one per stretch between its
  masks (in the order the masks are stored; a stretch longer than 5,000,000 letters in NCBI's
  overlapping chunks), each starting up to three masked letters early (at a multiple of four),
  so no word inside a mask starts a hit and no preliminary extension crosses a mask.
- The traceback (and `-ungapped`'s re-evaluation) reads the masked letters as `B`, so an
  alignment may still cross a mask at the cost of its mismatches.
- The identities that the reports count from the rows they show (outfmt 0, 5, 13 to 16, tabular
  rows with `qseq`, `sseq`, `positive`, `ppos` or `btop`) count the stored letters; the
  Seq-align's (outfmt 8, 9, 11, 12, the standard tabular columns, SAM's `PI`, `-sorthits` and
  "Max ident") count the masked letters as mismatches: 594 and 570 of 600 for an HSP across 24
  masked letters. Outfmt 0 shows the masks in lower case as with `-db_soft_mask`.
- Before anything is written, the masks of the searched sequences (after `-seqidlist`,
  `-taxids`, alias ranges and the like) and their lengths are read once, each sequence is cut
  into the chunks the search cuts (one chunker for both), and the search is refused at the first
  sequence whose masks NCBI's chunks do not handle; the error names it by its OID and volume:
  masks stored out of order or overlapping (hand-made mask data), whose result NCBI ties to the
  slices of its partial fetch; at a scanned word (the word size, the template's length for
  dc-megablast) below 11 letters (blastn-short's 7), a chunk from the sequence's start that a mask running to
  its end follows (the first mask with letters, with at most 5,000,000 letters and no other mask
  before it than masks of no letters at 0: all masked, or one masked tail such as a poly-A tail),
  which NCBI scans over the whole sequence past the chunk, in a sequence of at least the
  word's letters (NCBI scans every BLASTN subject, however short): rows that differ from those
  of the input with the mask one letter shorter, or a crash. At a scanned word
  of 11 or more (megablast, blastn, dc-megablast) such a sequence is searched over that chunk
  alone, the range before the mask, with the masked letters read as `B` in the traceback (the
  approved NCBI defect PD-LOSAT-NCBI-DEFECTS exception 13, the Owner's decision of 2026-11-15,
  line NP M4a-2c): NCBI's own result where NCBI ends normally and its scan past the chunk adds
  no alignment; where that scan adds one (a word across the mask's start, in an alignment with
  too few letters before the mask for a word), leaves diagonal state that drops seeds of the
  sequence searched next in the same thread (which sequence depends on NCBI's threads), or ends in a
  crash, the result of NCBI's normal chunks, as for the input with the mask one letter shorter; the
  sequences after it are searched as without it. After a
  mask of no letters elsewhere, or after more than 5,000,000 letters, the chunk from 0 is one of
  NCBI's normal chunks and the sequence is searched; and a last mask
  after which the final chunk starts before the end of an earlier chunk of another stretch (the
  mask is shorter than the 1 to 3 letters before its end that the final chunk starts with),
  where both chunks hold a word, so that the overlap that NCBI reads past its list of ranges
  decides which of their alignments merge across the mask. The chunk of a stretch of 0 to 3
  letters between two masks holds no word and so no alignment to merge, and NCBI does not merge
  a final chunk without one: such sequences are searched (review 2 of M4a-2b, 2026-11-15).
- That overlap, read one range past the list after a sequence's final chunk (another sequence's
  mask end, or memory never set), is `dbseq_chunk_overlap` (100) in LOSAT, as the search without
  masks already takes it after every sequence over 5,000,000 letters (NCBI reads
  `(&full_range)[1]` there). LOSAT refuses it only where it can decide a merge across a mask
  (above), a shape that only hand-made mask data give; where the last stretch between masks is
  longer than 5,000,000 letters, its final chunk merges with the full chunk before it in the
  same stretch, as in a sequence without masks, and is searched.
- These refusals come before the first line of output, where every limit of LOSAT is checked,
  so they apply also where NCBI would stop earlier at an error of the query batch or find no
  word to start a hit (a `-query_loc` outside the query, a query of `N` only).

Refused: a search that reads the masks of a volume that has masks but not the algorithm
(NCBI reads the masks of another volume's algorithm through a cache, or an unset value; M/DECISIONS
9 of line NP) or of a volume whose mask column has its big-endian file alone (NCBI reads its
numbers reversed). Such a volume is searched when no sequence that the search reads from it
(after `-seqidlist`, `-taxids` and the like) has mask data, as NCBI then reads none of it.
LOSAT Web refuses both options until it searches databases.

## What BLASTX `-db` of stage DB-5 does and rejects

`losat blastx -db NAME` searches a protein database of the reader above and gives the output of
NCBI BLAST+ 2.17.0 `blastx -db NAME -num_threads 1` for outfmt 0, 6 and 7, the same bytes at
every `-num_threads` (P1, P1.1). The `-subject` path is unchanged; both paths run the same
search, the same composition-based redo and traceback, and the same report writers. The database
opener, the reader's error texts and the title rules are those of BLASTP's `-db`.

- Residues: the search, the composition-based redo, the traceback and the reports read the stored
  bytes as NCBI's `-db` does. A stored `O` scores as `X` (NCBI's matrix) and is shown as `O`;
  `U`, `B`, `Z`, `J`, `*` as stored.
- Statistics: the number of sequences and the total length of the database. The shortest subject
  length behind the hit-saving cutoff is the database's 10 (BLASTX translates the query, so it is
  not divided), not the shortest sequence, so `-db` keeps the weak hits NCBI `-db` keeps.
  `-dbsize N` (NCBI's `Int8`), for `-db` and `-subject`: a value above 0 replaces the database
  length in the length adjustment and the effective search space, and a value other than 0 in
  the database length of the E-value statistics (NCBI's Gumbel block); the link cutoffs of an
  ungapped search, the counts and the reports keep the real numbers.
- Memory: the database's sequences are read once, in blocks of about 1 MiB, and kept as stored
  (about the size of the `.psq` file, and 8 bytes per sequence); NCBI maps the same file and its
  search touches all of it. A database larger than memory is not searched in blocks yet.
- Reports: the definition lines of the sequences a query batch reports are read before the
  batch's reports, where NCBI reads them; the `Database:` lines of outfmt 0, `>` and the stored
  title, the description lines made from the stored title, subject ids made from the title as
  for BLASTP, `stitle` as stored, `# Database:` of outfmt 7 with the `-db` value as given.
- Errors: NCBI's texts and exit statuses for a database that is not found or cannot be read, at
  the point where NCBI reports them. A database without sequences or a sequence file that cannot
  be read fails the first query batch with a valid context (NCBI does not search a batch without
  one), after the reports of the batches before; a run without such a batch succeeds. A header
  file that cannot be read stops the run at the batch that reports from it, after the reports of
  the batches before.
- `-db` with `-subject` is a usage error, as for BLASTP. NCBI's application settings (its
  environment and registry files) are read for `-db` as for BLASTP's `-db`. LOSAT Web rejects
  `-db` until stage DB-7.

Rejected explicitly, besides what BLASTP `-db` and the reader reject (the same rules apply to the
database's sequences and titles): a database sequence longer than 5,000,000 residues, as with
`-subject` a record (NCBI searches it in overlapping parts; refused at the first query batch with a
valid context, as the reader's errors, in one form, `database sequence N has M residues; ...` or
`subject sequence N has M residues; ...`, exit status 1),
`-searchsp`, `-db_soft_mask`, `-db_hard_mask`, the id and taxonomy lists, and the options BLASTX
rejects with `-subject`.

## What TBLASTX `-db` of stage DB-5 does and rejects

`losat tblastx -db NAME` searches a nucleotide database of the reader above and gives the output
of NCBI BLAST+ 2.17.0 `tblastx -db NAME -num_threads 1` for outfmt 0, 6 and 7, the same bytes at
every `-num_threads` (P1, P1.1). The `-subject` path is unchanged; both paths run the same
search, the same re-evaluation and linking, and the same report writers. The database opener,
the reader's error texts and the title rules are those of TBLASTN's `-db`.

- Bases: the seed and ungapped stages translate the stored 2-bit bases, and the ungapped
  re-evaluation, the identities and the reports read the letters with the stored ambiguity codes,
  as NCBI's search does. makeblastdb (and `losat makeblastdb`) draw the 2-bit base of each ambiguous
  letter from `CRandom` seeded with the sequence length, as `-subject` draws it at search time
  (writedb_convert.cpp, blast_objmgr_tools.cpp), so the bases of `-db` and `-subject` differ only
  where the FASTA has hyphens (stored as `N`, dropped by `-subject`), with `-subject_loc`, or for
  a database written otherwise; a frozen row searches such a database.
- Statistics: the number of sequences and the total length of the database, in 64 bits as NCBI
  keeps them. TBLASTX has no shortest-length cutoff (an ungapped search: the hit-saving length is
  the average length), so `-db` and `-subject` of sequences without ambiguous letters give the
  same hits. `-dbsize N` (NCBI's `Int8`), for `-db` and `-subject`: a value above 0 replaces the
  database length in the length adjustment and the effective search space; the average length,
  the linking cutoffs, the counts and the reports keep the real numbers. `-db_gencode` translates
  the database as it translates `-subject`.
- Memory: the database's sequences are read once, in blocks of about 1 MiB, and kept packed
  (about the size of the `.nsq` file); each thread makes the letters and the frames of the
  sequence it searches, as for `-subject`, and the reports read the letters of each reported
  HSP's codons from the packed sequences. A long sequence with ambiguity codes is searched with
  both its frame sets (about 6 bytes a base, above NCBI's peak; a later change). A database
  larger than memory is not searched in blocks yet.
- Reports: the definition lines are read for the reported sequences of each query batch, where
  NCBI reads them; the `Database:` lines of outfmt 0, `>` and the stored title, subject ids made
  from the title as for TBLASTN, `# Database:` of outfmt 7 with the `-db` value as given.
- Errors: as TBLASTN's `-db` (a database without sequences or a sequence file that cannot be read
  fails the first query batch with a valid context, after the reports of the batches before;
  damaged files are LOSAT's refusal, with no report).
- `-db` with `-subject` or `-subject_loc` is a usage error, as for BLASTP. LOSAT Web rejects
  `-db` until stage DB-7.

Rejected explicitly, besides what TBLASTN `-db` and the reader reject (the same rules apply):
`-searchsp`, `-db_soft_mask`, `-db_hard_mask`, and the options TBLASTX rejects with `-subject`.

## Databases with parsed Seq-ids (version 1.3, line DBID)

`blastp`, `blastn` (megablast, blastn, dc-megablast), `tblastn`, `blastx`, `tblastx` and `psiblast`
search a database of `makeblastdb -parse_seqids` with `-db`, of format version 4 or 5, NCBI's
distributed databases included. A sequence may have several definition lines, GIs, tax ids and
memberships. The output is that of NCBI BLAST+ 2.17.0 run with `-num_threads 1` in every format the
program writes, the same bytes at every `-num_threads` (P1, P1.1), apart from what this section lists.

- Ids: `LOSAT/src/seqid` holds NCBI's Seq-id types and its choice of id and label for each format,
  and `LOSAT/src/report/db_subject.rs` applies them to a reported sequence. A database without
  parsed Seq-ids keeps the ordinal-id path and its reports.
- Reading: a sequence's definition lines are decoded when it is reported. BLASTN and TBLASTX, which
  write their reports batch by batch, first read every header of a database of up to 2^20 sequences
  to learn whether a refusal below can occur. The Seq-align id of a reported sequence is looked up
  as NCBI's formatters look it up (`CSeqDB::SeqidToOid`: the LMDB keys of version 5, the ISAM files
  of version 4, under the search's id or taxonomy list; an ISAM with a file that is not a regular
  file, such as a directory, is none, as NCBI's `CSeqDBIsam::IndexExists` says, so the lookup
  finds nothing in that volume, since line DBID item loadord; the ISAMs of the volumes in their
  order up to the first one that holds the id, never a later one, since line DBID item v4late),
  whatever its kind: version 5 finds `bbs|`,
  `nat|` and `lcl|<integer>` by their keys, and `gnl|BL_ORD_ID|n` of the sequence's own number in
  its volume is found also beside another id. NCBI's scope first looks among the sequences of the
  batch it has already read, without case, so another reported sequence that holds the id can be
  shown in its place (refused below).
- Following hits whose Seq-align ids `Match` (`CSeq_id::Match`): outfmt 0, 5 and 13 to 16 show the
  alignments of copies of one record (the same definition lines but for their tax ids) that follow
  each other in one hit, with the description line, heading and names of the first and its HSPs
  then the others' (`Hsp_num` going on); `-num_descriptions`, `-num_alignments` and `psiblast`'s cut
  to `-max_target_seqs` count them as one subject. The other formats keep a row per copy (the
  tabular `qcovs` and `qcovus` count the group, which for copies is each copy's own coverage). Other
  sequences whose Seq-align ids `Match` are refused below.
- Approved difference: a `gb`, `emb` or `dbj` id prints the stored type
  ([`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md), version 1.10).

Rejected explicitly, with the sequence named. The last column is the rule of
[`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md) that applies: 3 where NCBI fails, 2 where it gives
a result that LOSAT does not reproduce.

| Rejected, for a reported sequence | Formats | NCBI BLAST+ 2.17.0 | Rule |
| --- | --- | --- | --- |
| A Seq-align id that the lookup does not find: `gnl|PIG|` (makeblastdb fills the PIG files from no FASTA id), in version 4 `bbs|` and `nat|` (makeblastdb writes no ISAM key for them), `gnl|BL_ORD_ID|n` where n plus the first ordinal id of each volume is no sequence (makeblastdb stores such an id only when the FASTA names it), any other id of a version 4 database without ISAM files; and a local id that alone decides the definition line (`lcl|777` before `gi|777`), which the lookup finds but the formatters do not; a `giim` id (`gim|123`), which the lookup finds but NCBI's Seq-id handles never take for one with an equal id, so the scope never finds its Bioseq | 0, 5, 13 to 17; 6, 7, 10, 20 with a field that reads the sequence (below) | The default tabular fields and any list with `sacc`, `saccver` or a field that reads the letters stop (exit 255). Other lists exit 0: `sseqid` writes `Unknown`, `sallseqid` the stored id. Outfmt 0 and 5 write `Unknown`, 13 to 16 and SAM the id | 2, 3 |
| A Seq-align id that the lookup finds first on another sequence, for example `gb|X|` after `gb|X.1|`, `emb|X.1|` after `gb|X.1|`, `gb|X|` after `emb|X.1|`, `lcl|X` after `gb|X.1|`, in version 5 `lcl|300` after `bbs|300` and `bbs|301` after `lcl|301` (one LMDB key), `gnl|BL_ORD_ID|n` where n plus the first ordinal id of a volume is another sequence. The lookup does not refuse the same record in two volumes or in two databases of an alias (with other tax ids: the row below) | 0, 5, 13 to 17; 6, 7, 10, 20 with a field that reads the sequence (below) | Writes the other sequence's names, in every program but BLASTN also its identities, where one of its ids equals the Seq-align id (`gb|X.1|` for `gb|X|`, a `gnl|BL_ORD_ID|n` that is its own number); otherwise as in the row above | 2, 3 |
| A Seq-align id that another sequence also holds, in any letter case, where the lookup finds the sequence itself (or the same record) first: `lcl|555` and `gb|X.1|` with `lcl|555` in a second definition line, `lcl|abc` and `lcl|ABC`, `bbs|123` in two databases of an alias, `pir||PIRNAME` and `pir||pirname`, `prf||PRFNAME` and `prf||prfname` or the names alone `gb||locusy` and `gb||LOCUSY` in two databases (version 5 keys these in their own letter case; LOSAT looks up the other keys of their letters and keeps the sequences that hold an id of the same handle: not `gb||locusy` and `emb||LOCUSY`, nor `sp||NAMEY` and `tr||namey`). Not the same record in two volumes or databases (with other tax ids: the row below) | 0, 5, 13 to 17; 6, 7, 10, 20 with a field that reads the sequence (below) | Where both are reported (also for different queries), its scope, which reads them in the order of its prefetch and compares ids without case, can write the other sequence's names, in every program but BLASTN also its identities; for `pdb|1ABC|A` and `pdb|1abc|A` in two databases it stops (exit 2) | 2, 3 |
| A Seq-align id for which LOSAT cannot list the other sequences that hold it, where NCBI's lookup does not read what fails: in version 4 an ISAM file of a volume after the first one that holds the id in the search's list that does not open (a regular file of mode 000; NCBI's lookup returns at the first volume that holds the id and never opens the later ones, since line DBID item v4late: LOSAT stopped there with NCBI's `Error: Could not open input file` before), in version 5 a page of the other keys of its letters (the row above) | 0, 5, 13 to 17; 6, 7, 10, 20 with a field that reads the sequence (below); where only the prefetch reads it, the formats whose taxonomy database's texts differ between the sequence that the lookup gives and none (another holder's Bioseq that its scope read before stands for it without a load) | Goes on (exit 0) with the names and letters that its lookup or its scope gives, as in the two rows above, or stops where the sequence the lookup gives is another one whose definition lines `GetBioseq` rejects for the id (`oid headers do not contain target gi/seq_id`); an id that only such a volume can hold, or any id where it comes first in the order of the volumes' names, opens it and stops (exit 2, `BLAST Database error: Error: Could not open input file`), which LOSAT writes where it knows what the batch's prefetch writes before the stop | 2 |
| The same record in another volume or database of an alias, or of `-db "a b"`, with other tax ids, which the lookup gives first or which also holds the Seq-align id. In the formats that show a title made from the tax ids (outfmt 0 for an empty title whose tax id the taxonomy database names or a protein title ending in `]` where a definition line has no title; the tabular formats where a printed id is the first word of the title that replaces an `lcl|Subject_N` id, for a nucleotide or a protein title ending in `]`), and in every format where the taxonomy database's texts depend on the copy: texts of opening it where only one copy has a tax id, the offset error of a tax id that a copy holds (also for the same tax ids in another volume) | 13 to 16 and the tabular `staxid`, `staxids` and name fields; for such a title the formats named; for the texts 0, 5, 6, 7, 10, 13 to 17, 20, and 8, 9, 12 with `PRE_FETCH_SEQS_LIMIT=2147483647` | Shows for both the copy that the lookup gives first: its tax ids and their names, the title made from them, and the taxonomy database's texts of the lookups that the volume of that copy makes (the offset error of `taxdb.btd` once for that volume); the other formats write the same bytes | 2 |
| A Seq-align id whose lookup fails (an LMDB file without the accession table where the batch's prefetch is off, or in a batch whose prefetch loads other sequences that throw on their ids, or, at a lookup that throws on its id, other sequences whose places to it in NCBI's order of the batch's Seq-id handles LOSAT does not know and whose loads write texts of the taxonomy database, above; `gnl|BL_ORD_ID|0`, which makeblastdb stores with a text tag) | 17; 0 to 7, 10, 13 to 16, 18 and 20 where the run writes another format too (formats that NCBI's runs of each stop at different places after the sequence, in a call of LOSAT's library: outfmt 0, outfmt 1 to 4, the formats that end the run there, outfmt 13 and 14 and the field lists below that stop, each at their rows, outfmt 15 and 16, which post the error on standard error, and the other formats, outfmt 5 and those that read nothing of the sequence among them, since reviews 1 and 2 of line DBID item pfoff-4, M1 and M3; also where the call makes one search for each hit list size of `-num_descriptions` and `-num_alignments`: the refusal is the call's, over all of its formats; NCBI writes each format in a run of its own, and those runs stop at different places with different exit statuses and texts, where the call has one result); 6, 7, 10 and 20 with the prefetch off at the first row of a query's report with `sacc` or `saccver` (the default fields) or `staxid`, and for the text tag at its row with `sacc` or `saccver` (since line DBID item pfoff-7, every other format of the text tag is written as NCBI writes it, with the prefetch on or off; the prefetch's failure there makes the batch's reports read one by one, so the refusals below of the prefetch off apply to it, also a weak match of another sequence that the reports read in such a batch); every format with the prefetch off where the run loaded before the sequence, in NCBI's order (a report of a query, batch, psiblast round or library call before, or a row, `Hit` or description line before it in its report; for psiblast's PSSM, which reads each sequence that it includes by its Seq-align id after the round's reports, in the order of the hit list, also with outfmt 8, 9 and 11, which read nothing of the sequences, those reports and the reads before it; of equal scores NCBI shows the later sequence of the database first, and it orders a `-db` list by name) and while NCBI keeps it, another sequence of the database whose own lookup does not fail and that NCBI gives for the Seq-align id (listed in the LMDB files that have the accession table): for the whole run one with an id of the same type, accession and version, in any letter case, in any of its definition lines, until the end of the query batch (of the run in outfmt 5) also one of another type of the accession and version or a weak match in the definition line that its load made its Bioseq of (the line of its own Seq-align id: its first line with a gi, else its first line; since review 3 of line DBID item pfoff-4, M1) (the versions of an accession without one, an id of the same name; in a version 5 database LOSAT cannot list the weak matches of a name or of a `pir` or `prf` id and takes any sequence); after a library call before this one, of which LOSAT has no record, any such sequence, each call a batch of its own; with the words of a sequence that also holds the Seq-align id, or of holders that LOSAT cannot list (since reviews 1 and 2 of line DBID item pfoff-4, H1 and H2, H1 and M1; a report after one whose first row ends NCBI's run at the null pointer of `sacc` or `saccver` is not reached, and that refusal's words stay; this also refuses runs whose output with the holder equals the one without it, M2: 534 runs of review 2's probe that LOSAT wrote as NCBI does before) | Where the run loaded before another sequence that NCBI gives for the Seq-align id (`gi|5001|ref|NP_100001.1|` in the first query's report, for `ref|NP_100001.1|` in the second query's), NCBI gives that sequence for it without the lookup, and every format writes that sequence's record and goes on (exit 0): its data loader keeps the ids of every definition line of every sequence that it loaded for the whole run and gives one for an id of the same type, accession and version in any letter case; its scope keeps the sequences that it loaded, each with the ids of the definition line that its load targets (the first that holds the id it loads by: for a hit its Seq-align id, so `gi|7051|ref|XP_900001.1|` then `gi|7052|ref|NP_100001.1|` is held by its first line only), until the end of the query batch (of the run in outfmt 5, `ResetScopeHistory`) and also gives one for another type of the accession and version (`emb|NP_100001.1|` for `ref|NP_100001.1|`) or a weak match (`ref|NP_100001.1|` for `ref|NP_100001|`, `gb|AAA09999.1|LOCUSN` for `gb||LOCUSN`); a holder that loads only after the throw, or one that the scope alone gives and an earlier batch loaded (outside outfmt 5), leaves it thrown, as below. Otherwise stops (exit 2); for the text tag (a lookup that throws on the id itself: a string tag of `PIG` or `BL_ORD_ID`, a `ti` tag that is not a number) its prefetch writes "Error pre-fetching sequence data" and the batch's formats load its sequences one by one, as with the prefetch off below, where a format that stops ends with `Error: ` and the exception (`CObject_id::GetId(): Invalid choice selection`, exit 255) in place of `BLAST Database error` (written by LOSAT since line DBID item pfoff-7); with the prefetch off writes `Unknown` for it, then outfmt 17 stops (exit 2) and 5 goes on (exit 0). Outfmt 5 with the prefetch off, which LOSAT writes as NCBI does since line DBID item pfoff-4 in BLASTP, psiblast (each round's `Iteration`), BLASTN (every task), BLASTX, TBLASTN and TBLASTX: each `Hit` looks its subject up by its first Seq-align's id, and where that throws (the `catch` of `s_SeqAlignToXMLHit`) its `Hit_id` is the Seq-align id's FASTA form (`ref|NP_100001.1|`, a local id's too), `Hit_def` and `Hit_accession` `Unknown` and `Hit_len` 0, and each of its HSPs (the `catch` of `s_SeqAlignSetToXMLHsps`) writes `Hsp_num`, the bit score, the score and the E-value, `Hsp_query-from`, `Hsp_query-to`, `Hsp_hit-from` and `Hsp_hit-to` 0, `Hsp_identity` the Seq-align's identities (`GetAlnScores`: -1 where it has none) and empty `Hsp_qseq` and `Hsp_hseq`, nothing else; nothing of it is loaded or written on standard error, and the report goes on with the next `Hit` (exit 0). Outfmt 0 to 4 with the prefetch off, which LOSAT writes as NCBI does since line DBID item pfoff-3 in BLASTP, psiblast (outfmt 0, each round's report), BLASTN (every task), BLASTX, TBLASTN and TBLASTX: the description line of such a sequence is its Seq-align id and `Unknown` (`NP_100001.1 Unknown`; with `-show_gis` the id's FASTA form and no space, `ref|NP_100001.1|Unknown`) with its scores, and loads nothing (the `catch` of `x_GetDeflineInfo`); outfmt 0 stops with `BLAST Database error` (exit 2) after the query's description lines where the first alignment that it shows (after `-sorthits` and `-num_alignments`; psiblast ignores `-sorthits` and `-sorthsps`, as NCBI's psiblast does, which never sets the hits' sort) is of such a sequence whose Seq-align id is not local (`GetDbType`, which has no `catch`), and otherwise writes `Sequence with id <accession without version, or the whole local id> no longer exists in database...alignment skipped` for each HSP of such a sequence and goes on (exit 0); outfmt 1 to 4 stop there where any alignment that they show is of such a sequence, and else write as before; with `-html`, outfmt 0 stops already at the description line of the first such sequence (exit 2), so a holder after it is not reached (review 3 of line DBID item pfoff-4, M2); `-num_alignments 0` shows no alignment and stops nowhere. The description lines load the sequences that they show in their order (in psiblast's rounds after the first, the two-part table of the first `-num_descriptions` HSPs, not hits, and at least the first, which NCBI adds before it counts, so `-num_descriptions 0` shows its sequence too; the part of the previous round's set first: `x_SplitSeqAlign`), then the alignments those that they show up to the stop; after the stop nothing loads, so the taxonomy database's texts are those of these loads. A report that stops the run (outfmt 0 to 4 here, 13 and 14, and the field lists below that stop) ends NCBI's application there: no later query's report runs and no later batch is prefetched, so a later batch's prefetch writes no text and does not stop (BLASTP, psiblast and TBLASTN, which make every batch's prefetch before the reports, since review 3 of line DBID item pfoff-3). Outfmt 13 to 16 with the prefetch off, which LOSAT writes as NCBI does since line DBID item pfoff, in BLASTP, BLASTN, BLASTX, TBLASTN and TBLASTX: the query's report throws at that `Hit` (`DB contains no accession info.`) and loads nothing of it or of the `Hit`s after it, so the taxonomy database's texts are those of the loads before it, and a text that NCBI writes once (opening the taxonomy database, a volume's first lookup of a tax id) comes at the next load that NCBI makes, a later query's `Hit` or a later batch's prefetch (and its stop); 15 and 16 post `Error: [<program>] ` with the exception, leave the query's report out and go on (exit 0), 13 and 14 leave the query's file empty and stop with `BLAST Database error` (exit 2). Outfmt 6, 7, 10, 20 (10 after the line of the field names) and 18 with the prefetch off, which LOSAT writes as NCBI does since line DBID item pfoff-2 in the same programs (in psiblast the tables; its outfmt 18 is not written): a field list with a field of the subject block (an id, title or taxonomy field, `sall*`, `sstrand`; the default fields) writes the row that the block's `catch` leaves (`slen`, `length`, the counts and coordinates 0, `pident` 0.000, `ppos` 0.00, `frames` 0/0, `sstrand` N/A, the letters empty; the scores and `qcovs` as usual; `qlen` only with a field of the query block; the ids, titles and taxonomy of the query's last row that loaded its subject, else `Unknown`, `sgi` 0 and `N/A`; `sall*` that row's lists and the Seq-align id of each row since) and goes on; a field list that reads the sequence only after that block (`slen`; the letters and the counts that read them: `positive`, `ppos`, and outside BLASTN `mismatch`, `pident` and `nident`; in BLASTX, TBLASTN and TBLASTX `sframe` and `frames` with a field of the alignment block: a coordinate, `length`, `gaps`, `gapopen`, `qframe`) writes the rows before and stops with `BLAST Database error` (exit 2), but `nident` without another such field reads the sequence only at a row with identities, so a row of none prints 0 and goes on; a row of a throwing load loads nothing, so the taxonomy database writes no text for it and the names and tax ids are looked up only for the subjects that a report loads (a damaged taxonomy database's text of opening comes at the first row that loads its subject, and `taxdb.btd` is not read for a throwing subject's tax ids). The organism report (18) leaves out each HSP whose load throws (the `catch` of `x_InitTaxInfoMap`) and goes on (exit 0). At the first row of a query's report NCBI stops at `sacc` and `saccver` (null pointer, exit 255) and prints an uninitialised `staxid` | 3 |
| Following hits of a query whose Seq-align ids `Match` but that are not copies of one record: `pdb|1ABC|` and `pdb|1ABC|AA`, `gb|X|` and `gb|X.1|`, one accession and version under `gb`, `emb` or `dbj` | 0, 5, 13 to 16; tabular `qcovs` and `qcovus` | Shows them as one hit; writes the coverage of both on the first row and `N/A` on the other | 2 |
| In `psiblast`, two sequences whose Seq-align ids differ but that NCBI's Seq-id handles can take for one, where NCBI compares them: one with an alignment below `-inclusion_ethresh` and one of the previous round's inclusion set (convergence); two such, where the round makes a PSSM (`-save_pssm_after_last_round` with `-out_pssm` or `-out_ascii_pssm`, a round that makes the next round's PSSM); in outfmt 0, a reported sequence and one of the previous round's set (LOSAT checks every reported sequence, a superset of the `-num_descriptions` first that NCBI's two-part table looks up). Such ids: ids that differ only in letter case (`gb|ABC1.1|` and `gb|abc1.1|`), one accession under `gb`, `emb` and `dbj` (`gb|ABC1.1|` and `emb|ABC1.1|`, one handle where the accession packs), and ids that LOSAT does not tell apart as the handles do (an integer tag and the same digits as a string, a patent's number and application number, PDB ids that differ only in their release or the case of the molecule's first letter). Also a `giim` id (`gim|123`) with an alignment below `-inclusion_ethresh` where the round makes a PSSM | every format | Takes them for one sequence in its inclusion sets (convergence, the two-part description table) and counts them as one row of the PSSM's multiple alignment (`CSeq_id_Handle`); never takes a `giim` id for one handle with an equal id (its tree's search for an equal id is inverted), so it counts each of its alignments as a row and cannot read its letters | 2 |
| In `psiblast` where the run makes a PSSM, two copies of one record with another hit between them in the hit list before the cut, both with an alignment below `-inclusion_ethresh` | every format | Counts one row for the copies (`x_CountAndSelectQualifyingAlignments`) but moves to a new row at each (`x_ExtractAlignmentData`) and writes past the last row: it crashes (SIGSEGV, exit 139) | 3 |
| In `psiblast`, such hits that its cut of a longer list to `-max_target_seqs` counts (copies of one record count as one) | every format | `PruneSeqalign` counts them as one subject and keeps one more | 2 |
| In `psiblast`, a sequence of its longer list beyond `-max_target_seqs`, which no report shows: one whose Seq-align id's lookup throws on the id (`gnl|BL_ORD_ID|0` stored with a text tag) where another sequence of the longer list, which the failing prefetch may have loaded, holds an id that NCBI's scope takes for it, or where a shown sequence that the formatter then reads one by one can match an id weakly (otherwise written as NCBI writes it since line DBID item pfoff-7); one for which the prefetch can read another sequence whose definition lines do not hold the id (`pdb|1abc|A` for `pdb|1ABC|A` in two databases of an alias; line DBID P3), or, where the taxonomy database's texts can differ between what NCBI's prefetch can read for it: another sequence that also holds the id (as the blobs read before in the run decide), or the same record in another volume (an id found nowhere is not read, as by NCBI; another sequence that the lookup finds first with an id that compares equal is read, its tax ids looked up, as by NCBI, also for a subject whose definition line the formatters do not find) | 0, 6, 7, 10, 20 | For the failed lookup its prefetch writes "Error pre-fetching sequence data" and leaves the shown sequences it has not loaded yet to the formatter, which reads them one by one (its scope can match an id of a sequence read before weakly, as with `PRE_FETCH_SEQS_LIMIT`), and its PSSM drops that sequence's letters ("Failed to retrieve sequence"); otherwise it stops (exit 2, `oid headers do not contain target gi/seq_id`) or looks up the other sequence's tax ids | 2, 3 |
| Where NCBI's prefetch reads a reported sequence by its Seq-align id and no report reads it (the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`, at which the prefetch reads every reported sequence also for them; a tabular field list that reads nothing of the sequence, below; psiblast's sequences beyond `-max_target_seqs`), in a batch whose prefetch is on (a positive `PRE_FETCH_SEQS_LIMIT` turns off the prefetch of a batch whose Seq-aligns outnumber it; each batch counts its own): another sequence that the prefetch can read for it whose definition lines do not hold the id (`pdb|1ABC|A` and `pdb|1abc|A` in two databases of an alias); where the taxonomy database's texts can differ between what the prefetch can read, as in psiblast's row above; a Seq-align id whose lookup fails on a file of the database (an LMDB page that does not read: NCBI's prefetch rethrows the `CSeqDBException`); and a Seq-align id whose lookup throws on the id in a batch whose other reported sequences' loads make the taxonomy database write a text (its texts of opening, the offset error of a tax id) where LOSAT does not know whether NCBI's prefetch makes them before that lookup (a Seq-id handle of its packed value, or neither packed; a Seq-align id of two sequences; a batch of more than 200 ids or with a handle that LOSAT does not know; since line NP item O4b2-10, the loads of the handles before it are written and those after it are not made). Without a taxonomy database on the search path too, but for the texts | 8, 9, 12; 6, 7, 10, 20 without a field that reads the sequence | Stops (exit 2, `oid headers do not contain target gi/seq_id`, or `BLAST Database error`); looks up the other sequence's tax ids; the prefetch loads the batch's sequences in the order of its Seq-id handles (the ids its mapper packs, such as GIs and accessions, by their packed value, then the others by the address of each id's record, which LOSAT does not model) and ends at the failed lookup with "Error pre-fetching sequence data" | 2, 3 |
| A Seq-align id or display id that is the id of a query read with `-parse_deflines` in the query batch that reports it (with outfmt 5, in any batch: NCBI's formatter empties its scope after each batch but there); in BLASTP's, BLASTN's, BLASTX's, TBLASTN's and TBLASTX's outfmt 1 to 4 also the id `lcl|Query_N` of a query read without it (any letter case). The display id is the id of the sequence's line in the description table (`GetDisplayIds`, `lcl|abc` of `gi|779|lcl|abc`), by which that line looks the title up again (`showdefline.cpp:498`) | the Seq-align id: the formats in which NCBI's scope looks the sequence up by it, a report that reads its Bioseq or the batch's prefetch (`ParsedShown::resolves_subject`): every format but 8, 9, 11 and 12, which read only the Seq-align (they too with `PRE_FETCH_SEQS_LIMIT=2147483647`), and with `PRE_FETCH_SEQS_LIMIT=0` but a table whose fields read nothing of the sequence; psiblast every format where the call makes a PSSM from the round's Seq-aligns (another round follows or the last PSSM is saved, `psiblast_app.cpp:296-310`), whose PSSM reads the aligned letters from that scope (`psi_pssm_input.cpp:409-445`; narrowed from every format on 2026-11-24, line RG-KB item 12 review 2, N1, where item 12's key of a pdb id by its handle had extended it to such formats; to the calls that make a PSSM on 2026-11-25, review 3 of item 12, F1); the display id: 0 to 4 in BLASTP, psiblast, rpsblast, BLASTN, BLASTX, TBLASTN and TBLASTX, where the table lists the sequence (one of the first `-num_descriptions` subjects; in a psiblast round with the previous round's set, a subject of the first `-num_descriptions` HSPs, at least one); for `lcl|Query_N` both in 1 to 4 only (outfmt 0 is line RG item X7's). The Seq-align id is not refused where NCBI's scope finds one query of the batch for the sequence by it (by its own keys, `QueryScopeKeys::own`, not a weak match) and every report of the run shows NCBI's bytes from that query's Bioseq (2026-11-25, line RG-KB item 12, decision KB-12-1 and review 3, F1 and F2, `report::query_bioseq`, one owner of the run's query Bioseqs that the tables and the organism report read): BLASTP's (psiblast's where the call makes no PSSM and searches no `-in_pssm`), BLASTN's, BLASTX's, TBLASTN's and TBLASTX's tables whose fields do not read the found Bioseq's letters or hit BLASTX's `CreateTranslatedDensegFromNADenseg` exception (not `sseq`, `btop`, `positive` or `ppos`; outside BLASTN, whose formatter fetches no sequence for them, not `mismatch`, `pident`, `nident` or the default fields; in BLASTN not `qseq` combined with `mismatch`, `pident` or `nident`; in TBLASTN not `qseq`; in TBLASTX not `qseq` when an HSP reaches past the found query's length; in BLASTX not any field of `opens_counts`, i.e., `qstart`, `qend`, `sstart`, `send`, `length`, `gapopen`, `gaps`, `qframe`, `sstrand`, `qseq`, `sseq`, `pident`, `nident`, `mismatch`, `positive`, `ppos` or `btop`; NCBI throws in most of them), so also a table that reads nothing of the sequence where only the batch's prefetch looks it up (which finds the query and loads nothing of the sequence), the organism report and outfmt 8, 9, 11 and 12, without `-query_loc` and without a taxonomy database whose loads write a text or throw; the tables then show the query's ids, `N/A` for `stitle` and `salltitles`, the query record's length for `slen` and no taxonomy (`sframe`/`frames` in TBLASTX and TBLASTN from the found query's length, TBLASTN `sstrand` `N/A`), and the organism report the query's line without a taxid (`tabular.cpp:836-866`) | Shows the query's title and letters for it; for the display id the query's title in that line | 2 |
| A Seq-align id or display id (row above) by which NCBI's scope finds a query read with `-parse_deflines` in the query batch that reports it (with outfmt 5, in any batch) as a weak match (`x_GetSeqMatch`): `sp||NAME` for `sp|P00001.1|NAME` and the other way, `gb|X.1|` for `gb|X|`, `pir||ABCD` for `pir|ABCD|` | the Seq-align id: the formats of the row above; the display id: 0 to 4 in BLASTP, psiblast, rpsblast, BLASTN, BLASTX, TBLASTN and TBLASTX, where the table lists the sequence (as in the row above) | Shows the query's title and letters for it; for the display id the query's title in that line (not `gb|X|` for `gb|X.1|`, nor `pir|ABCD|` for `pir||ABCD`) | 2 |
| In a query batch whose prefetch `PRE_FETCH_SEQS_LIMIT` turns off (0 or below, or a value below the batch's number of Seq-aligns): a Seq-align id that an id of another sequence matches weakly, as in the row above, in any letter case and across the alternatives (`tpg|X.1|NAME` and `prf||NAME`). In version 4 where the lookup of the Seq-align id or of its name alone finds such a sequence; in version 5 an accession without a version where the lookup finds a version of it, and any Seq-align id with a name or of `pir` or `prf` (version 5 keys a name only in its own letter case and only without an accession, and a `pir` or `prf` id by its FASTA text) | 0, 5, 6, 7, 10, 13 to 17, 20 | It reads the batch's reported sequences one by one (the prefetch of a batch loads each by its own id), and where it has read that other sequence first, it writes that sequence's names (in BLASTP also its identities; psiblast's PSSM also reads its letters) | 2 |
| A stored Seq-id with a byte outside ' ' to '~', or a stored title with such a byte in the remaining refused cases of `PD-LOSAT-NCBI-DEFECTS` 1.120 and 1.127 (makeblastdb stores `#` for such a byte; line RG-KB items 9(b) and 9b-2) | stored Seq-id: every format; stored title (except a query batch that ends at a header NCBI cannot read, which writes NCBI's `Error: (803.7)` messages then NCBI's error in any format of all eight search programs): every format for `psiblast`, `rpsblast` and `rpstblastn`, `-outfmt` 0 to 5, 13 to 16 and 18, with 6, 7, 10, 17 and 20 (or 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647`) a subject whose blob NCBI may still hold or not or whose header cache slot such a subject takes (rule 3), and, not modelled yet (line RG-KB item 9b-3), a prefetch that may be off, parsed Seq-ids, `OBJMGR_*` settings, prefetch loads that write texts of the taxonomy database, searches with a seqid, GI or taxid list or an alias/volume list or OID mask, or a query batch with a UTF8String header warning | Prints `#` and an error message (`Error: (803.7)`) each time it reads the line, exit 0 | 2 |
| A definition line without a Seq-id (makeblastdb writes none) that no list of the search erases, in a sequence with a Seq-id in another line (the paragraph below the table): the `sallacc` of each row of the sequence | 6, 7, 10, 20 with `sallacc` | Stops at the row's `sallacc` after its earlier fields (null pointer, exit 255) | 3 |
| The same line in the organism report | 18 | Stops at that line for each HSP of the sequence and shows the lines before it once per HSP (Owner question Q-D1c1cb-18 of line RG item D1c-1c-b) | 2 |
| An outfmt 0 description id of 65 bytes or more | 0 | Pads the line without end, from 65 bytes with a title of 3 bytes or more, from 67 with any title | 3 |
| A protein title of punctuation that `x_CleanAndCompress` reads past | 0; tabular formats for an `lcl|Subject_N` id | Crashes | 3 |
| A title in no encoding that NCBI recognises (or the name in the taxonomy database from which NCBI makes an empty title, or which `x_AdjustProteinTitleSuffix` writes at the end of a protein title ending in ` [ ]`), for an `lcl|Subject_N` id: the rows that `SetSubjectId` cannot make (the subject Bioseq's ids hold that id), and the first row of a query's report that either `SetSubjectId` or `x_SetSubjectIds` cannot make (`sallseqid`, `sallgi` or `sallacc` of a sequence whose definition line with that id follows one with another id): `x_SetTaxInfo` comes after them in the same `try` and never runs (`tabular.cpp:839-854`) | 6, 7, 10, 20 with `sacc` or `saccver` (the default fields of 6, 7, 10); `staxid` in such a first row | Stops at the row's `sacc` or `saccver` (null pointer, exit 255); prints an uninitialised `staxid` | 3 |

A definition line without a Seq-id (a hand-edited database: makeblastdb writes none; under any list,
also one left without ids or taxids, NCBI's header read erases such a line and LOSAT does the same:
rows `db.il.*.ln.*`) is NCBI's in BLASTP, psiblast, BLASTN, BLASTX, TBLASTN and TBLASTX, since line
RG item D1c-1c-b, and in `rpsblast` and `rpstblastn` (makeprofiledb writes none either; rows
`db.il.rps.ln.*`, `db.il.rtbn.ln.*`) since 2026-11-27. A sequence whose lines all lack one gets no
Seq-align: its hits leave the report after the cut of the hit list, which they keep counting (rows
`db.il.*.ln.q2*`, `db.il.*.lnn.*`, `db.il.rps.ln.c.*`; a profile of an RPS database of several volumes
leaves in its own volume's search, before the volumes' lists are merged, and takes no place there:
`db.il.rps.ln.3vc.*`). Beside
lines with Seq-ids, each format's walk of the lines takes the best id of such a line, null, and its
`catch` decides (rows `db.il.*.x1.*`, `db.il.bp.xw.*`, `db.il.*.ln.qa.*`, `db.il.*.pn.qa.*`): the
outfmt 0 to 4 description line is the Seq-align id and `Unknown` (its GI with `-show_gis`; with
`-html` an id link without a score link) where `GetDisplayIds` labels such a line before the line of
the display id or the walk after the title meets it before the text passes 4096 bytes (a text of
4096 bytes walks on, one of 4097 stops), and the heading of each HSP is the skipped id; outfmt 5's
`Hit` is the Seq-align id's FASTA form, `Unknown`, `Unknown` and length 0 where a titled line without
one (an empty title counts) follows the first titled line, its HSPs as always; XML2 and JSON keep
the `HitDescr` of the lines before the first such line, with their taxonomy, and then the Seq-align
id alone, with the sequence's length; the tabular formats print `Unknown` in `sallseqid`, 0 in
`sallgi` and `N/A` for an untitled line in `salltitles`, and are refused only at `sallacc` (the
table above), after the reports of the query batches before, which NCBI writes before it stops
there (rows `db.il.*.qab1.rej.6`, `db.il.ps.ln.qa*.rej.6`, `db.il.bp.x1.q01b1.rej.6`; so does every
refused tabular row where NCBI stops or prints an uninitialised `staxid`, in every program, BLASTN's
reports written after the loop too: rows `db.il.bn.pn.qwb1.{dbs,chk,chkl}.rej.6`), and after what NCBI
writes to stderr of the failing batch before that row: the batch's reader messages, its search
strategy posts, its UTF8String warning, its prefetch texts, and for each of its queries up to the
row's the query's own warnings, the taxdb warning of its formatter (a taxonomy name field without
`taxdb.bti`) and, with the prefetch off, its subject reads up to the row (rows `db.il.*.qwb1.*.6`,
`db.il.ps.ln.qwq.*.6`, `db.il.*.qo*.wrn.6`, `db.il.*.tx.rej.6`, `db.il.*.qwa.pf0.rej.6`; the rows
marked `ncbi-stderr` compare that whole stderr with NCBI's).

A sequence without parsed Seq-ids or with an `lcl|Subject_N` id whose title, or the name in the
taxonomy database from which NCBI makes an empty title or the end of a protein title ending in ` [ ]`
(line RG item T5c, rows `db.t5.sfu.*`), is in no encoding that NCBI recognises is NCBI's in every format: the
tabular rows of NCBI's catch (`Unknown`, the Seq-align id in `sallseqid`, `gnl|BL_ORD_ID|n` or
`Subject_N`, after the ids of the definition lines before the first that holds `lcl|Subject_N`; the
titles and taxonomy of the row before) and refused only where NCBI stops or prints an uninitialised
number (the row above; [`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md) version 1.37, line RG
item T7).

A tabular field list reads the sequence that NCBI's formatters find by the Seq-align id only with
the default fields, a subject id (`sseqid`, `sgi`, `sacc`, `saccver`, `sallseqid`, `sallgi`,
`sallacc`), a title, tax id or name, `sstrand`, `slen`, or a field of the alignment's letters
(`qseq`, `sseq`, `positive`, `ppos`, `btop`, and in every program but BLASTN `pident`, `nident` and
`mismatch`; in BLASTX, TBLASTN and TBLASTX `sframe` or `frames` with a field of the alignment). Another
list (`6 qseqid evalue`, the coordinates, `length`, `gaps`, `gapopen`, `qcovhsp`, `qcovs`) writes
the same bytes whatever the lookup finds, so the rows above that name the tabular formats do not
refuse it, nor do the ids of a stored title that it does not print; the prefetch of a batch still
reads such a sequence (the row of the prefetch), unless the batch's Seq-aligns outnumber a positive
`PRE_FETCH_SEQS_LIMIT`. Where that prefetch fails, at a Seq-align id whose lookup fails or, in BLASTN
and TBLASTX, at an alignment on the minus strand of the query and the plus strand of a sequence it
gets no Bioseq of (an id found nowhere; the sequence that the lookup gives first, the sequence itself
or another one, has no id that compares equal, or the definition line that NCBI's reader keeps for
the id does not hold it, as where a local id alone decides that line; BLASTN's minus-strand HSPs, in
an ungapped search the sequence's first HSP; in TBLASTX the first HSP of the sequence's alignment),
LOSAT writes NCBI's `Warning: [<program>] Error pre-fetching sequence data` once for the batch, after
the taxonomy database's texts of its loads (at a failed lookup, of the loads that NCBI's order of the
batch's Seq-id handles puts before it) and before the batch's query warnings, and the rows;
where the data loader reads no sequence for the id and another sequence also holds it, whose
Bioseq the scope can hold from an earlier read (NCBI warns in some runs and not in others), it
refuses the run. Its `qcovs` and `qcovus` count following hits whose
Seq-align ids `Match` as one (the row of such hits). psiblast with a PSSM file (`-out_pssm`,
`-out_ascii_pssm`), and each of its rounds that makes the PSSM of the next round
(`-num_iterations` other than 1), reads the letters of each aligned sequence by its Seq-align id for
the PSSM, also for the sequences of its longer list beyond `-max_target_seqs`: those with an
alignment below `-inclusion_ethresh`, the only ones the PSSM reads (in a round whose prefetch
`PRE_FETCH_SEQS_LIMIT` turns off, one by one). Each alignment reads its range of the sequence that
NCBI's scope gets by that id: the sequence itself, or the other sequence that the lookup gives first
where one of its ids compares equal to the Seq-align id (one accession under `gb`, `emb`, `dbj`, `tpg`
or `tpe`, another version) and the definition line that NCBI's reader keeps for that id holds it;
where it gets none (a `ref` or local id of the same text, an id found nowhere, a lookup that throws on
the id, a giim id, a sequence, itself or another, whose definition line that NCBI's reader keeps for
the id does not hold it, as `lcl|777` before `gi|777`) or the range ends beyond that sequence's
letters, LOSAT writes NCBI's `Warning: [psiblast] Failed to retrieve sequence` with the Seq-align id
for the alignment, after the round's report (with one round, after the query's report and before
the next query's warnings, the last query's before the epilog), and the alignment's row of the
multiple alignment becomes the query's letters, which the purge removes. In a batch whose prefetch is
off, a read that loads a sequence no fetch of the run has loaded before (also a sequence without
parsed Seq-ids, which the PSSM reads by `gnl|BL_ORD_ID|n`) looks each taxid of its definition lines up
in its volume's 256-slot taxonomy cache, which keeps a taxid for the run's later queries and rounds and
asks the taxonomy database at the volume's first lookup of a taxid and where the taxid's slot
(`taxid & 255`) holds another one. LOSAT writes the taxonomy database's texts of those lookups (a wrong
magic number or count at the run's first lookup, the offset error of the last entry at each lookup of
its taxid that asks) where NCBI writes them: after the query's report, with the PSSM's warnings in the
order of its reads (twice for one sequence whose definition lines hold that taxid before and after
another taxid of its slot, again for a later sequence after a load of such a taxid), also where a PSSM
error stops the run after that query (where the query's report read the sequence first, or an earlier
lookup of the run wrote the texts of opening, there is nothing to write). NCBI's scope does not keep a
query's sequences for its later rounds and queries: a read of a sequence loaded before loads it again
or not (after 40 other included subjects it does, after none it does not), which writes nothing where
the slot holds that taxid either way. Refused: an id that
another sequence also holds besides the one the lookup gives (NCBI's scope can read that one where it
has read it before, and with the prefetch on as its memory layout falls), a lookup that fails on a
file, the same record whose read writes taxonomy texts, in a batch whose prefetch is off an id that
another sequence's id matches weakly, a read of a sequence loaded before whose load again would write
the offset error again (another taxid of its slot was looked up since), and a load after such a read
that may have put another taxid of the slot there (whether NCBI writes the text there follows what its
scope keeps).
Its multiple alignment has a row for each handle of these Seq-align ids, and a Seq-align whose id
`Match`es the previous Seq-align's id goes into the previous row, as NCBI's `CPsiBlastInputData`
makes them: copies of one record that follow each other share a row, and so do two other sequences
whose ids `Match` (`gb|X|` and `gb|X.1|`, two handles, which leave the last row empty), each with the
letters read by its own id.

A nucleotide title of punctuation is not refused: it stops at the end of the string, in outfmt 0 and
in the tabular first word of an `lcl|Subject_N` id (approved exception 2 of
[`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md)).

The environment variables `SEQ_ID_PREFER_ACCESSION_OVER_GI` and `SEQ_ID_AVOID_GI`, which change NCBI's
choice of ids, are refused when the database has parsed Seq-ids (any value). `[BLAST] LONG_SEQID` is
refused in every search.

Known differences that remain:

- Refusals wider than NCBI's failure, which rule 4 of `PD-LOSAT-NCBI-DEFECTS` would port:
  - A sole `gnl|BL_ORD_ID|n` of another number whose n names a sequence with the same definition
    lines and letters, which NCBI shows as the subject itself: that sequence is read without parsed
    Seq-ids, so the lookup does not take it for the same record.
  - A Seq-align id that another sequence also holds is refused also where NCBI writes the subject's
    own report: where that sequence is not reported, or its scope reads the subject first; in
    outfmt 8, 9 and 12 with `PRE_FETCH_SEQS_LIMIT=2147483647` (and in psiblast's list beyond
    `-max_target_seqs`) also where the prefetch reads the subject itself, when that sequence's
    definition lines do not hold the id.
  - In a batch whose prefetch `PRE_FETCH_SEQS_LIMIT` turns off, a Seq-align id that an id of another
    sequence matches weakly is refused also where NCBI writes the subject's own report: where that
    sequence is not reported or not read first, and where that id is in another definition line
    than the one NCBI reads that sequence by. In version 5 a Seq-align id with a name, or of `pir` or `prf`, is
    refused whatever the database holds.
  - The weak matches ignore the release of an id (`sp` and `tr` ids carry `reviewed` and
    `unreviewed`; NCBI matches a name only within one release), in both refusals of a weak match.
  - A weak match with a query read with `-parse_deflines` is refused also in outfmt 8, 9, 11 and 12,
    which write the Seq-align id alone (as the exact match).
  - A tabular field list that reads the sequence without the fields on which NCBI stops
    (`6 sseqid`, `6 sallseqid`, `6 stitle`), where NCBI exits 0: on the rows of a sequence that it
    does not find it writes `Unknown` (`sseqid`), the query's earlier rows' ids and titles
    (`sallseqid`, `stitle`, which it does not reset), and the alignment's fields as reset.
  - A title ending in `]` without `taxdb.bti` on the search path, where NCBI keeps the title.
  - The same record with other tax ids is refused in outfmt 13 to 16 and the tabular tax id and name
    fields also where NCBI writes the subject's own tax ids (the copy it shows is the subject, or is
    not reported), in the formats that show a title made from the tax ids also where that title
    does not change (the copy it shows is the subject; a protein title ending in `]` without
    `taxdb.bti` on the search path), and in every format where the taxonomy database's texts depend
    on the copy, also where the bytes do not change (the copy with a tax id fetched first, which
    makes the same first lookup).
- A damaged LMDB file can stop NCBI and LOSAT with other exit statuses and texts.

## Evidence

Design and measurements of 2026-10-12 are kept outside the repository by the campaign
(`dbsearch-20261012`: design, probe tables, decision record). The reader's tests are in
`LOSAT/src/seqdb/tests.rs` with fixtures written by NCBI BLAST+ 2.17.0 under
`LOSAT/tests/seqdb`. The writer's tests are in `LOSAT/src/writedb` (`tests.rs`,
`fasta_source/tests.rs`) with runs of makeblastdb 2.17.0 recorded in the generated tables
beside them. The tests of `-db` (BLASTP, BLASTN, TBLASTN, BLASTX and TBLASTX) are frozen NCBI outputs replayed by
`LOSAT/tests/db_search_fixtures.py` (rows in `LOSAT/tests/db_search_manifest/`, inputs in
`LOSAT/tests/db_search`, outputs in `LOSAT/tests/fixtures/db_search`; run by
`ci_fast_regressions.py` at 1 and 4 threads), the search-path cases recorded from NCBI in
`LOSAT/tests/seqdb/expected/search_path`, the unit tests of
`LOSAT/src/algorithm/blastp/blast_engine/database.rs` and
`LOSAT/src/algorithm/blastn/blast_engine/database.rs` and
`LOSAT/src/algorithm/blastx/database.rs`, `LOSAT/tests/blastn_db_search.rs` (BLASTN `-db` at 1,
2 and 4 threads) and `LOSAT/tests/blastx_db_search.rs` (BLASTX `-db` at 1, 2 and 4 threads, the
shortest-length cutoff, the first batch with a valid context); for TBLASTN the unit tests of
`LOSAT/src/sequence/packed_subject.rs` (`from_stored` against the reader's letters), `LOSAT/tests/tblastn_db_search.rs` (TBLASTN
`-db` at 1, 2 and 4 threads, the shortest length of a database, the sequence file read by the
first batch with a valid context) and the TBLASTN rows of the same manifest; for TBLASTX
`LOSAT/tests/tblastx_db_search.rs` (TBLASTX `-db` at 1, 2 and 4 threads and against `-subject`,
`-dbsize` on both paths, the sequence file read by the first batch with a valid context) and the
TBLASTX rows of the same manifest. For parsed Seq-ids: the `db.id.*` rows of the same manifest (formats 0, 5 to 9,
11 to 17, and most refusals above; the failed lookup is a unit test of `db_subject.rs`), with their
databases in `LOSAT/tests/db_search/ids` and `LOSAT/tests/seqid`; the ID4 sweep (BLASTP and BLASTX
on swissprot and pdbaa, BLASTN on 18S_fungal_sequences and 16S_ribosomal_RNA, TBLASTN and TBLASTX
on 18S_fungal_sequences, at 1 and 4 threads: 256 rows equal to NCBI's, the 2 SAM rows without their
`@PG` line, last on 2026-11-01, comment on PR #683); outfmt 10 and 20 on the fixtures by a probe
of the line (`blastdb-ids-20261021/R/id9`: 48 cases, also at 4 threads); outfmt 10 and 20 on the same
pairs of programs and databases, and PSIBLAST on swissprot and pdbaa (outfmt 0, 6, 7, 10, 20), at 1
and 4 threads (39 cases equal to NCBI's, 2026-11-02). A scan of those four databases found no two
sequences with one id and no `pdb` id without a chain.
