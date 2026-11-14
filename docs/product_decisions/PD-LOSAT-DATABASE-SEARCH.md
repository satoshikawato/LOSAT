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
  only for a table that prints a subject id, which version 1.0 refused for every table)
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
`makeblastdb -parse_seqids` is searched as NCBI searches it; refused is a run whose reports or
prefetch look a parsed Seq-id up in a file without the table, where NCBI stops with `DB contains
no accession info.` (`PD-LOSAT-NCBI-DEFECTS`, rule 3; before line DBID item mixdb, every such
list where NCBI starts the search). `losat blastdbcmd` keeps one
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
accession tables). Refused until item D1a-3: `NSEQ` negative as NCBI's `int`, `LENGTH` 0 or at
least 2^63 while sequences are included (NCBI computes the statistics of each subject, or reads
a negative `Int8`), `LENGTH` of 10^18 or more (NCBI's commas lose digits); until item D1a-4:
OID ranges, `STATS_NSEQ`, `STATS_TOTLEN` and `MIN_SEQ_LENGTH` in RPSBLAST databases. Refused for
good: a statistics count of 0 (an alias file's `NSEQ 0`, or its `LENGTH` or `STATS_TOTLEN` over
volumes without sequences without `NSEQ` or `STATS_NSEQ`) in a search without the Gumbel block
(BLASTN, TBLASTX, and BLASTP, BLASTX, TBLASTN, PSI-BLAST with `-searchsp` or `OLD_FSC`, BLASTX
also with `-ungapped`), where NCBI divides by zero (SIGFPE). With the Gumbel block, a TBLASTN database whose `GetMinLength` is
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
| Alias keys `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `OIDLIST`, `OID_MASK_TYPE`, `MEMB_BIT`, `MASKLIST` | filters (P4); `FIRST_OID`, `LAST_OID` read since line RG D1a, `VDBLIST` ignored as NCBI's reader ignores it |
| Alias keys `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`; `NSEQ` or `LENGTH` whose value is not the number of sequences or the total length of the alias file's databases | read since line RG D1a (2026-11-12) as NCBI's walks; the values listed above as refused until D1a-3 and D1a-4 |
| A name that lists several databases | one database per search (P4); searched since version 1.2 |
| A volume listed more than once, a group alias file (`index.alx`) in a directory, two volumes of one name in different directories | NCBI derives counts or the volume order from them in ways the reader does not reproduce; since version 1.2 two volumes of one name in different directories are read (NCBI orders them by base name, then path); since line RG D1a a volume listed more than once is read once and counted per alias node, as NCBI does |
| Column files of a volume (mask data of `makeblastdb -mask_data`) | not read; searched since 2026-11-11 (line SW, SW8 follow-up 2): NCBI's searches open no column file without a database mask named by `-db_soft_mask`, `-db_hard_mask` or an imported search strategy (`-import_search_strategy`), so LOSAT searches the volume without reading them, whatever they hold; since 2026-11-14 BLASTN reads the mask column for `-db_soft_mask` (line NP, M4a-2a; LOSAT refuses the other two) |

Two errors are NCBI's own and are reproduced: volumes of version 4 and of version 5 in one
database, and a version 5 volume whose LMDB file is missing.

Left unread without an error: the id and taxonomy files of a version 5 database, which serve
lookups by id or taxonomy only, the LMDB file's volume list, which NCBI compares with the
volumes, and (since 2026-11-11) the column files of a volume, which NCBI's searches read only for
a database mask named by `-db_soft_mask`, `-db_hard_mask` or an imported search strategy (since
2026-11-14 BLASTN reads the mask column for `-db_soft_mask`; LOSAT refuses the other two). The LMDB file is opened when the
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
| `-randomize`, `-gi_mask`, `-gi_mask_name`, `-oid_masks`, `-metadata_output_prefix` | not implemented; each is an explicit error that names the option (`-gi_mask` and `-gi_mask_name` write GI-indexed mask files, a later item of line NP) |
| `[BLAST] MAX_SEQID_LENGTH` (registry or `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`) with `-parse_seqids`, when makeblastdb uses the value (the first layer that has the entry gives a number `NStr::StringToUInt` converts, not 0; an empty, unconvertible or 0 value is ignored by makeblastdb and accepted by LOSAT); `[BLAST] LONG_SEQID` | they change which FASTA ids makeblastdb parses; LOSAT keeps NCBI's defaults and refuses to run |
| `-input_type` other than `fasta` | FASTA only (plain, or named `.gz`, `.bz2`, `.zst`: inflated while it is read, in any case of the suffix, as makeblastdb decides it; `-in -` reads the standard input as it comes). Every gzip member, bzip2 stream and zstd frame of a concatenated file is read, where makeblastdb reads every gzip member but only the first bzip2 stream or zstd frame and drops the rest silently (`PD-LOSAT-NCBI-DEFECTS`, the makeblastdb exception). A `.gz` file is read as makeblastdb's decompressor (`CZipDecompressor`) reads it (`GzipMembers` in `LOSAT/src/blastinput/compressed.rs`, which inflates with zlib-rs, a Rust port of zlib's inflate, where makeblastdb calls the host's zlib: the official Linux build links `libz.so.1`, and a host whose `libz.so.1` is zlib-ng's zlib-compatible build could in principle differ on damaged data): neither a member's trailer (CRC-32 and length) nor its header CRC is checked, a trailer cut short reads as the whole member, a member's header must end within the next 4 KiB, bytes after the last member are read only where makeblastdb reads them (fewer than 4 KiB that inflate to no byte, such as up to four zero bytes or a line feed), and a file that does not start with a gzip header is read as a zlib stream, up to that stream's end. A file that is not what its name says, or is cut short (for `.gz`: where makeblastdb's decompressor fails too), stops the build with `losat makeblastdb: -in: 'x.gz' is named .gz but is not gzip data (...)` (exit status 1, the files removed as for another failed input); makeblastdb warns about the decompressor's error, says the file does not appear to be FASTA, and ends with `BLAST failed to write output` (exit status 6) |
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
  `-mask_id`) end with exit status 2 and the text of LOSAT's option
  parser, as in every LOSAT subcommand; makeblastdb uses 1.
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

- A string with a byte outside 0x20-0x7e, which NCBI repairs and logs, and a `UTF8String` tag
  where a `VisibleString` belongs (`mask data that LOSAT does not read: ...`): LOSAT does not
  write NCBI's diagnostics for them.
- Ranges whose order NCBI takes from memory addresses: ranges of two different Seq-id handles
  with the same packed number (`gi|5` and `gb|AY000005.1|`; a local id and a `general` string
  tag; two `giim` handles of one number) written on one strand of one sequence and algorithm,
  when the two orders write different ranges. That order shows only in the files of a build that
  succeeds, so the error comes at the end of such a build and names the Seq-ids and the sequence;
  an error of makeblastdb that comes later in the run (a mask of the same sequence past its end, a
  later sequence or input file that fails) is reported as makeblastdb reports it. The files are
  left as makeblastdb leaves them after an error of its build at the last sequence (none in
  version 4; `NAME.nos`, `NAME.not`, `NAME.ntf`, `NAME.nto` in version 5, `.p..` for proteins),
  where the refusal of a mask string comes before the build and leaves none.

A mask file that NCBI's reader rejects as damaged is no refusal: both programs stop, LOSAT with
exit status 255 and `Error: mask data in '<file>': <LOSAT's text>` where NCBI prints an
exception report that names its build tree.

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

The creation date is the local time, found without the C library: the zone comes from `TZ`,
`TZDIR` and `/etc/localtime` as glibc reads them. On Windows and under WASI, and when no
zone can be read, the date is UTC. Zone files with leap-second records give a time 27
seconds later than the C library's, and a `TZ` value with a daylight-saving name but no
rule follows the United States rules unless `posixrules` says otherwise.

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
  of a query batch are made).
- `-db` with `-subject` or `-subject_loc` is a usage error (LOSAT's parser text and exit
  status 2; NCBI prints its usage and exits 1). LOSAT Web rejects `-db` until stage DB-7.

Rejected explicitly (NCBI runs these, or crashes, or LOSAT cannot reproduce what it does):

| Rejected | Why |
| --- | --- |
| `-searchsp` | NCBI then uses a different statistics path (`OLD_FSC`, no Gumbel block); a stage of its own |
| A database sequence with a gap residue (`-`, stored by makeblastdb), a byte above the protein alphabet, or a sequence not between two NUL bytes | NCBI's `-db` never aligns a stored gap and its formatter counts it as a gap; not reproduced |
| A database sequence without parsed Seq-ids with more or fewer than one definition line | makeblastdb without `-parse_seqids` writes one; a sequence with parsed Seq-ids may have several (version 1.3) |
| A reported subject whose title NCBI's cleaner reads past (titles like `, ,`: NCBI crashes) | tabular formats and outfmt 0 |
| In outfmt 0, a shown subject with an empty title and a tax id (makeblastdb `-taxid`) that the taxonomy database names (line DBID H3; all six programs) | NCBI makes the title from the tax id and the taxonomy database, as for a sequence with parsed Seq-ids: `unnamed protein product [Homo sapiens]`, for a nucleotide `Homo sapiens`. For a tax id that it does not name (not in it, no `taxdb.btd`, a database that does not open) NCBI writes the stored title, as LOSAT does; the other formats show the stored title or its first word, but a nucleotide's tabular subject id is the first word of the title made from the tax id (BLASTN's row below) |
| A reported title with a raw byte outside ' '..'~' in the header file | NCBI prints `#` and a warning a number of times that depends on the format and the field |
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
  and read `-db_hard_mask` up to the search, which it refuses.
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
| `-searchsp`, `-db_soft_mask`, `-db_hard_mask` | a stage of its own; masks are not read (P4). Since 2026-11-14 `-db_soft_mask` is searched and `-db_hard_mask` is refused where it would mask the search (line NP, M4a-2a; section "BLASTN database masks") |
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

## BLASTN database masks (line NP, M4a-2a, 2026-11-14)

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

Refused: `-db_hard_mask` where it would mask the search (its errors and its warning are NCBI's;
line NP, M4a-2b), a search that reads the masks of a volume that has masks but not the algorithm
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
database's sequences and titles): a database sequence longer than 5,000,000 residues (NCBI searches it
in overlapping parts; refused at the first query batch with a valid context, as the reader's errors),
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
  of version 4, under the search's id or taxonomy list), whatever its kind: version 5 finds `bbs|`,
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
| The same record in another volume or database of an alias, or of `-db "a b"`, with other tax ids, which the lookup gives first or which also holds the Seq-align id. In the formats that show a title made from the tax ids (outfmt 0 for an empty title whose tax id the taxonomy database names or a protein title ending in `]` where a definition line has no title; the tabular formats where a printed id is the first word of the title that replaces an `lcl|Subject_N` id, for a nucleotide or a protein title ending in `]`), and in every format where the taxonomy database's texts depend on the copy: texts of opening it where only one copy has a tax id, the offset error of a tax id that a copy holds (also for the same tax ids in another volume) | 13 to 16 and the tabular `staxid`, `staxids` and name fields; for such a title the formats named; for the texts 0, 5, 6, 7, 10, 13 to 17, 20, and 8, 9, 12 with `PRE_FETCH_SEQS_LIMIT=2147483647` | Shows for both the copy that the lookup gives first: its tax ids and their names, the title made from them, and the taxonomy database's texts of the lookups that the volume of that copy makes (the offset error of `taxdb.btd` once for that volume); the other formats write the same bytes | 2 |
| A Seq-align id whose lookup fails (an LMDB file without the accession table; `gnl|BL_ORD_ID|0`, which makeblastdb stores with a text tag) | 0, 5, 6, 7, 10, 13 to 17, 20 | Stops (exit 2; for the text tag "Error pre-fetching sequence data", exit 255) | 3 |
| Following hits of a query whose Seq-align ids `Match` but that are not copies of one record: `pdb|1ABC|` and `pdb|1ABC|AA`, `gb|X|` and `gb|X.1|`, one accession and version under `gb`, `emb` or `dbj` | 0, 5, 13 to 16; tabular `qcovs` and `qcovus` | Shows them as one hit; writes the coverage of both on the first row and `N/A` on the other | 2 |
| In `psiblast`, two sequences whose Seq-align ids differ but that NCBI's Seq-id handles can take for one, where NCBI compares them: one with an alignment below `-inclusion_ethresh` and one of the previous round's inclusion set (convergence); two such, where the round makes a PSSM (`-save_pssm_after_last_round` with `-out_pssm` or `-out_ascii_pssm`, a round that makes the next round's PSSM); in outfmt 0, a reported sequence and one of the previous round's set (LOSAT checks every reported sequence, a superset of the `-num_descriptions` first that NCBI's two-part table looks up). Such ids: ids that differ only in letter case (`gb|ABC1.1|` and `gb|abc1.1|`), one accession under `gb`, `emb` and `dbj` (`gb|ABC1.1|` and `emb|ABC1.1|`, one handle where the accession packs), and ids that LOSAT does not tell apart as the handles do (an integer tag and the same digits as a string, a patent's number and application number, PDB ids that differ only in their release or the case of the molecule's first letter). Also a `giim` id (`gim|123`) with an alignment below `-inclusion_ethresh` where the round makes a PSSM | every format | Takes them for one sequence in its inclusion sets (convergence, the two-part description table) and counts them as one row of the PSSM's multiple alignment (`CSeq_id_Handle`); never takes a `giim` id for one handle with an equal id (its tree's search for an equal id is inverted), so it counts each of its alignments as a row and cannot read its letters | 2 |
| In `psiblast` where the run makes a PSSM, two copies of one record with another hit between them in the hit list before the cut, both with an alignment below `-inclusion_ethresh` | every format | Counts one row for the copies (`x_CountAndSelectQualifyingAlignments`) but moves to a new row at each (`x_ExtractAlignmentData`) and writes past the last row: it crashes (SIGSEGV, exit 139) | 3 |
| In `psiblast`, such hits that its cut of a longer list to `-max_target_seqs` counts (copies of one record count as one) | every format | `PruneSeqalign` counts them as one subject and keeps one more | 2 |
| In `psiblast`, a sequence of its longer list beyond `-max_target_seqs`, which no report shows: one whose Seq-align id's lookup throws on the id (`gnl|BL_ORD_ID|0` stored with a text tag) where the report reads the sequences it shows or the run writes a PSSM; one for which the prefetch can read another sequence whose definition lines do not hold the id (`pdb|1abc|A` for `pdb|1ABC|A` in two databases of an alias; line DBID P3), or, where the taxonomy database's texts can differ between what NCBI's prefetch can read for it: another sequence that also holds the id (as the blobs read before in the run decide), or the same record in another volume (an id found nowhere is not read, as by NCBI; another sequence that the lookup finds first with an id that compares equal is read, its tax ids looked up, as by NCBI, also for a subject whose definition line the formatters do not find) | 0, 6, 7, 10, 20 | For the failed lookup its prefetch writes "Error pre-fetching sequence data" and leaves the shown sequences it has not loaded yet to the formatter, which reads them one by one (its scope can match an id of a sequence read before weakly, as with `PRE_FETCH_SEQS_LIMIT`), and its PSSM drops that sequence's letters ("Failed to retrieve sequence"); otherwise it stops (exit 2, `oid headers do not contain target gi/seq_id`) or looks up the other sequence's tax ids | 2, 3 |
| Where NCBI's prefetch reads a reported sequence by its Seq-align id and no report reads it (the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`, at which the prefetch reads every reported sequence also for them; a tabular field list that reads nothing of the sequence, below; psiblast's sequences beyond `-max_target_seqs`), in a batch whose prefetch is on (a positive `PRE_FETCH_SEQS_LIMIT` turns off the prefetch of a batch whose Seq-aligns outnumber it; each batch counts its own): another sequence that the prefetch can read for it whose definition lines do not hold the id (`pdb|1ABC|A` and `pdb|1abc|A` in two databases of an alias); where the taxonomy database's texts can differ between what the prefetch can read, as in psiblast's row above; a Seq-align id whose lookup fails on a file of the database (an LMDB page that does not read: NCBI's prefetch rethrows the `CSeqDBException`); and a Seq-align id whose lookup throws on the id in a batch whose other reported sequences' loads make the taxonomy database write a text (its texts of opening, the offset error of a tax id). Without a taxonomy database on the search path too, but for the texts | 8, 9, 12; 6, 7, 10, 20 without a field that reads the sequence | Stops (exit 2, `oid headers do not contain target gi/seq_id`, or `BLAST Database error`); looks up the other sequence's tax ids; the prefetch loads the batch's sequences in the order of its Seq-id handles (the ids its mapper packs, such as GIs and accessions, by their packed value, then the others by the address of each id's record, which LOSAT does not model) and ends at the failed lookup with "Error pre-fetching sequence data" | 2, 3 |
| A Seq-align id that is the id of a query read with `-parse_deflines` in the query batch that reports it (with outfmt 5, in any batch: NCBI's formatter empties its scope after each batch but there) | every format | Shows the query's title and letters for it | 2 |
| A Seq-align id by which NCBI's scope finds a query read with `-parse_deflines` in the query batch that reports it (with outfmt 5, in any batch) as a weak match (`x_GetSeqMatch`): `sp||NAME` for `sp|P00001.1|NAME` and the other way, `gb|X.1|` for `gb|X|`, `pir||ABCD` for `pir|ABCD|` | every format | Shows the query's title and letters for it (not `gb|X|` for `gb|X.1|`, nor `pir|ABCD|` for `pir||ABCD`) | 2 |
| In a query batch whose prefetch `PRE_FETCH_SEQS_LIMIT` turns off (0 or below, or a value below the batch's number of Seq-aligns): a Seq-align id that an id of another sequence matches weakly, as in the row above, in any letter case and across the alternatives (`tpg|X.1|NAME` and `prf||NAME`). In version 4 where the lookup of the Seq-align id or of its name alone finds such a sequence; in version 5 an accession without a version where the lookup finds a version of it, and any Seq-align id with a name or of `pir` or `prf` (version 5 keys a name only in its own letter case and only without an accession, and a `pir` or `prf` id by its FASTA text) | 0, 5, 6, 7, 10, 13 to 17, 20 | It reads the batch's reported sequences one by one (the prefetch of a batch loads each by its own id), and where it has read that other sequence first, it writes that sequence's names (in BLASTP also its identities; psiblast's PSSM also reads its letters) | 2 |
| A stored title or Seq-id with a byte outside ' ' to '~' (makeblastdb stores `#` for such a byte) | every format | Prints `#` and an error message each time it reads the line, exit 0 | 2 |
| A definition line without a Seq-id (makeblastdb writes none) | every format | Not compared | |
| An outfmt 0 description id of 65 bytes or more | 0 | Pads the line without end, from 65 bytes with a title of 3 bytes or more, from 67 with any title | 3 |
| A protein title of punctuation that `x_CleanAndCompress` reads past | 0; tabular formats for an `lcl|Subject_N` id | Crashes | 3 |
| A protein with a tax id and a title ending in `]`, where a definition line has no title (0) or the id is `lcl|Subject_N` (tabular formats) | 0, tabular formats | With `taxdb.bti` on the search path `x_AdjustProteinTitleSuffix` rewrites the title from the tax ids; without it NCBI keeps the title | 2 |
| An empty title of a sequence with a tax id that the taxonomy database names (not a `pdb` id with a chain or a patent) | 0; tabular formats that print a nucleotide `lcl|Subject_N` id | Makes the title from the tax ids and the taxonomy database (the tabular id is its first word, the organism's; for a protein `unnamed`, as LOSAT writes it); for a tax id the database does not name, the title of a sequence without one | 2 |
| A title in no encoding that NCBI recognises, for an `lcl|Subject_N` id | tabular formats | Not compared | |

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
the taxonomy database's texts of its loads and before the batch's query warnings, and the rows;
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
