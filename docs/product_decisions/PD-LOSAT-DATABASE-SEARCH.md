# Product Decision: database search (`-db`)

- Decision ID: `PD-LOSAT-DATABASE-SEARCH`
- Version: 1.1
- Date: 2026-10-13 (version 1.0: 2026-10-12)
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
  BLASTX and TBLASTX), and all five programs search a database. `-searchsp` stays an explicit
  rejection until a stage of its own.

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
database files, byte for byte, for outfmt 0, 6 and 7. `-db` and `-subject` of one FASTA
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
| A header with a Seq-id other than `gnl|BL_ORD_ID|n` | parsed Seq-ids (P4); checked on the first sequence of every volume when the database is opened, and at every header read after that |
| Alias keys `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `OIDLIST`, `OID_MASK_TYPE`, `FIRST_OID`, `LAST_OID`, `MEMB_BIT`, `MASKLIST`, `VDBLIST` | filters (P4) |
| Alias keys `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`; `NSEQ` or `LENGTH` whose value is not the number of sequences or the total length of the alias file's databases | they replace the values behind the statistics (second step). `NSEQ` and `LENGTH` equal to those sums, as `blastdb_aliastool -dblist` writes them, change nothing in NCBI and are read (since stage DB-3) |
| A name that lists several databases | one database per search (P4) |
| A volume listed more than once, a group alias file (`index.alx`) in a directory, two volumes of one name in different directories | NCBI derives counts or the volume order from them in ways the reader does not reproduce |
| Column files of a volume (mask data of `makeblastdb -mask_data`) | not read |

Two errors are NCBI's own and are reproduced: volumes of version 4 and of version 5 in one
database, and a version 5 volume whose LMDB file is missing.

Left unread without an error: the contents of the LMDB file and the id and taxonomy files of
a version 5 database, which serve lookups by id or taxonomy only (NCBI also compares the
LMDB file's volume list with the volumes; the reader does not). Since stage DB-3 a `-db`
name is looked up along NCBI's search path (the current directory, then `BLASTDB`, then
`[BLAST] BLASTDB` of `.ncbirc`), with NCBI's texts when it is not found. Not ported, each an
error or a difference that a test records: group alias files (`index.alx`), several names
in one `-db`, index files cut short, NCBI's threshold of open files and its warning at
1,024 links, escapes in registry values.

## What `losat makeblastdb` of stage DB-2 writes and rejects

`losat makeblastdb -in FASTA[ FASTA...] -dbtype prot|nucl -out NAME [-title T] [-max_file_sz S]
[-blastdb_version 4|5] [-logfile F] [-parse_seqids] [-hash_index] [-taxid N | -taxid_map F]`
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
the version 5 tax-id files. Several input files can be named in
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
the same bytes as on Linux.

| Not written or rejected | Why |
| --- | --- |
| `-randomize`, `-mask_data`, `-mask_id`, `-mask_desc`, `-gi_mask`, `-gi_mask_name`, `-oid_masks`, `-metadata_output_prefix` | not implemented; each is an explicit error that names the option |
| `[BLAST] MAX_SEQID_LENGTH` (registry or `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`) with `-parse_seqids`; `[BLAST] LONG_SEQID` | they change which FASTA ids makeblastdb parses; LOSAT keeps NCBI's defaults and refuses to run |
| `-input_type` other than `fasta`, compressed input (`.gz`, `.bz2`, `.zst`) | uncompressed FASTA only; `-in -` reads the standard input |
| NCBI settings that change makeblastdb's messages or stored titles (`DIAG_*`, `NCBI_CONFIG_*`, registry entries outside the accepted list, a non-empty `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS`) | as for the search programs: LOSAT refuses to run instead of giving a different result |
| Replacing a database of the same name whose alias file the reader of stage DB-1 rejects | the files to remove are found by resolving the alias file; LOSAT stops before it removes anything |

### Where `losat makeblastdb` knowingly differs from makeblastdb

Removal of an existing database and the paths of the output:

- An output path with a space or a double quote in its absolute form (for example a build
  inside a directory `my data`). makeblastdb reads the path as a list of names: it removes
  nothing of an existing database, writes the new files beside the old ones and ends with
  exit status 2. LOSAT takes the path as one name, replaces the database and ends with 0.
  (`-out` itself may not contain a space, in both programs.)
- An existing `sub/db.pal` while the current directory holds volumes of the same names
  (`./db.00.*`). NCBI looks for a volume in the current directory first, so makeblastdb
  removes the files of the current directory. LOSAT removes the volumes beside the alias
  file. Other behaviour of makeblastdb is kept: the existing database is removed before the
  input is read, so a failing build also removes it.
- A volume that an existing alias file names and that is found only through `BLASTDB` or
  `.ncbirc` is not removed (`losat makeblastdb` does not use the search path of `-db`).
- Directories that makeblastdb creates get mode 0755; LOSAT does the same on Unix only.

Version 5 and the options added to the writer:

- The `dbname` of the metadata file is the file name of the output path, which LOSAT takes
  as one name (makeblastdb splits a path with a space and ends with exit status 2 before it
  writes the file). A file `NAME.pjs` left by a replaced database of several volumes is
  listed and counted, as in makeblastdb.
- The totals of the metadata file come from the volumes the writer made; the volumes are
  not opened afterwards (opening them would hold the index of every sequence in memory).
- `-logfile` is read as makeblastdb's option, not wherever it stands in the command line
  (makeblastdb's corelib also finds it as the value of another option).
- Of several `-in` files, one that exists and cannot be opened is an error (exit status 255)
  that removes what the earlier files added; makeblastdb skips it without a word. `-in ""`
  (no name) is refused; makeblastdb reads an empty list.
- LMDB's lock files (`NAME.pdb-lock`, `NAME.ptf-lock`) are not created (they change no byte
  of a database; one left at the path is removed, makeblastdb's final state). So an output
  name of 247 to 251 bytes, whose `NAME.pdb-lock` passes the file-name limit of 255 bytes,
  stops makeblastdb ("BLAST Database error: LMDB runtime error: mdb_env_open: File name too
  long", exit status 2, no volume written) while LOSAT builds the database (exit status 0;
  every file is valid and its own name fits).

Errors that makeblastdb does not detect, or detects elsewhere:

- A file that cannot be created, written or renamed, and an input file that cannot be
  opened: makeblastdb does not test its streams and ends with exit status 0 and a damaged
  or missing database. LOSAT removes the files of the build and ends with
  `Error: <path>: <system error>` and exit status 255.
- An output directory that exists but is not writable: makeblastdb stops before its log;
  LOSAT stops when it creates the first file, after the first lines of the log.
- `-max_file_sz` that is no size or is 4 GiB or more: LOSAT checks it before it removes an
  existing database or writes the log; makeblastdb after. The exit status is the same.
- Usage errors (unknown option, missing `-dbtype`, `-taxid -1`, `-taxid_map` without
  `-parse_seqids` or with `-taxid`) end with exit status 2 and the text of LOSAT's option
  parser, as in every LOSAT subcommand; makeblastdb uses 1.
- Errors of makeblastdb's general handler carry the source path of its binary
  (`Error: NCBI C++ Exception: ...`); LOSAT prints `Error:` and the sentence.
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

Reproduced although it looks wrong (deterministic NCBI behaviour): with `-hash_index`, no
`-parse_seqids` and several volumes, the first sequence of a later volume is indexed under the
ordinal id it would have had in the earlier volume, so makeblastdb and LOSAT both stop with
"Duplicate seq_ids are found: GNL|BL_ORD_ID|n" when that volume holds more than n sequences.
`-taxid_map` lines end at a line feed on every platform (makeblastdb on macOS also ends them at a
carriage return).

Undefined behaviour of NCBI that is not reproduced:

- A volume file past 4 GiB is an error (makeblastdb's 32-bit offsets wrap).
- A version 5 database of more than 2^31 sequences is an error at the sequence that would get
  OID 2^31 (makeblastdb's `int` OID overflows). The number of sequences in the metadata file
  is an `Int8`, the type of the field, summed without NCBI's `int`.
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
  stored title, subject ids made from the title as NCBI makes them (`gnl|BL_ORD_ID|n` is
  never shown), `# Database:` of outfmt 7 with the `-db` value as given.
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
| A database sequence with more or fewer than one definition line | makeblastdb without `-parse_seqids` writes one; several need the parsed-Seq-id rules (P4) |
| A reported subject whose title gives no subject id (a title of spaces: NCBI shows `gnl|BL_ORD_ID|n`), or whose title NCBI's cleaner reads past (titles like `, ,`: NCBI crashes) | tabular formats, and outfmt 0 for the second |
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
  and outfmt 0 holds the letters of the reported sequences until the reports are written.
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
| `-searchsp`, `-db_soft_mask`, `-db_hard_mask` | a stage of its own; masks are not read (P4) |
| A reported subject whose title gives no subject id: a title of spaces, or an empty title (NCBI writes `gnl|BL_ORD_ID|n` and `BL_ORD_ID:n`; for a protein an empty title gives `unnamed`, which BLASTP writes) | as BLASTP: ids of the general form are not written |
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

## Evidence

Design and measurements of 2026-10-12 are kept outside the repository by the campaign
(`dbsearch-20261012`: design, probe tables, decision record). The reader's tests are in
`LOSAT/src/seqdb/tests.rs` with fixtures written by NCBI BLAST+ 2.17.0 under
`LOSAT/tests/seqdb`. The writer's tests are in `LOSAT/src/writedb` (`tests.rs`,
`fasta_source/tests.rs`) with runs of makeblastdb 2.17.0 recorded in the generated tables
beside them. The tests of `-db` (BLASTP, BLASTN, TBLASTN, BLASTX and TBLASTX) are frozen NCBI outputs replayed by
`LOSAT/tests/db_search_fixtures.py` (rows in `LOSAT/tests/db_search_manifest.tsv`, inputs in
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
TBLASTX rows of the same manifest.
