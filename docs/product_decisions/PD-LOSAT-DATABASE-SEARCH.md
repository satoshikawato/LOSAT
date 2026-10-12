# Product Decision: database search (`-db`)

- Decision ID: `PD-LOSAT-DATABASE-SEARCH`
- Version: 1.0
- Date: 2026-10-12
- Status: Accepted by the maintainer on 2026-10-12, after the probe and design of the
  database search campaign (decisions P1 to P4 below, each as recommended). One decision
  is open (section "Open decision"); the maintainer takes it before stage DB-3.
- Implementation status: stage DB-1 added the database reader (`LOSAT/src/seqdb`); no
  search program calls it yet. Stage DB-2 adds the subcommand `LOSAT makedb`
  (`LOSAT/src/writedb`, section "What `LOSAT makedb` of stage DB-2 writes and rejects").
  Until a stage adds them, `-db`, `-dbsize` and `-searchsp` stay explicit rejections.

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

`LOSAT <program> -db NAME` must give the output of NCBI BLAST+ 2.17.0 `-db NAME` on the same
database files, byte for byte, for outfmt 0, 6 and 7. `-db` and `-subject` of one FASTA
therefore differ exactly as NCBI's do. This replaces the earlier wording "E-values are
computed with the size of the whole database", which `-subject` already meets.

### P2. Database format

A LOSAT database is an NCBI BLAST database: the volume layout that makeblastdb writes.
LOSAT reads format versions 4 and 5 directly and writes the version 4 layout with
`LOSAT makedb`. There is no LOSAT-only format.

Reasons: the comparison oracle searches the very same files, so ordinal ids, the stored
2-bit bases and the title are shared by construction; existing databases are usable as they
are; the writer can be checked byte for byte against `makeblastdb -blastdb_version 4`.

### P3. Command-line surface

- `-db`, `-dbsize` and `-searchsp`, under NCBI's names and with NCBI's meaning.
- A subcommand `LOSAT makedb` with makeblastdb's option names (`-in`, `-dbtype`, `-out`,
  `-title`, `-max_file_sz`).

### P4. Scope of the first release

- Databases without parsed Seq-ids: makeblastdb without `-parse_seqids`, and `LOSAT makedb`.
- One database per search.
- No id or taxonomy filters.

Databases distributed by NCBI (they store parsed Seq-ids), id and taxonomy filters, alias
files that filter or replace the statistics, and several databases in one `-db` are a
second step. Until then LOSAT rejects them explicitly; it never searches such a database
with part of its data ignored.

## Open decision (the maintainer's, before stage DB-3)

NCBI's own `-db` output can depend on `-num_threads`. With composition-based statistics,
blastp ends the redo of a query early by the hit lists of its query batch, and with more
than one thread NCBI keeps one set of hit lists per thread. Measured: a random
10,500-residue query followed by PUR4_HUMAN against Swiss-Prot gives 42 different output
rows at 8 threads and at 1 thread. P1 names "NCBI `-db` on the same database" without a
thread count.

To decide: which NCBI thread count LOSAT's `-db` reproduces. Recommendation of the
campaign (not a decision): NCBI with one thread, at every LOSAT thread count, as `-subject`
does today. Comparisons with NCBI state the thread count in the meantime.

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
| A header with a Seq-id other than `gnl|BL_ORD_ID|n` | parsed Seq-ids (P4); checked on the first sequence of every volume when the database is opened |
| Alias keys `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`, `OIDLIST`, `OID_MASK_TYPE`, `FIRST_OID`, `LAST_OID`, `MEMB_BIT`, `MASKLIST`, `VDBLIST` | filters (P4) |
| Alias keys `NSEQ`, `LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH` | they replace the values behind the statistics (second step) |
| A name that lists several databases | one database per search (P4) |
| A volume listed more than once, a group alias file (`index.alx`) in a directory, two volumes of one name in different directories | NCBI derives counts or the volume order from them in ways the reader does not reproduce |
| Column files of a volume (mask data of `makeblastdb -mask_data`) | not read |

Two errors are NCBI's own and are reproduced: volumes of version 4 and of version 5 in one
database, and a version 5 volume whose LMDB file is missing.

Left unread without an error: the contents of the LMDB file and the id and taxonomy files of
a version 5 database, which serve lookups by id or taxonomy only (NCBI also compares the
LMDB file's volume list with the volumes; the reader does not). The reader looks for a
database at the given path; NCBI's search path (`BLASTDB`, `.ncbirc`) belongs to the stage
that adds `-db`.

## What `LOSAT makedb` of stage DB-2 writes and rejects

`LOSAT makedb -in FASTA -dbtype prot|nucl -out NAME [-title T] [-max_file_sz S]` writes the
files of `makeblastdb -blastdb_version 4` without `-parse_seqids`: sequence files
(`.psq`, `.nsq`) and header files (`.phr`, `.nhr`) byte for byte, the index file (`.pin`,
`.nin`) but for the creation date, and for several volumes the alias file (`.pal`, `.nal`)
but for its date line. The FASTA file is read as makeblastdb reads it (its reader settings
differ from those of `-subject`), the log on the standard output, the messages on the
standard error stream and the exit status are makeblastdb's, and an existing database of the
same name is replaced as makeblastdb replaces it.

| Not written or rejected | Why |
| --- | --- |
| The JSON metadata file (`.pjs`, `.njs`) | NCBI's search programs do not read it; their output is the same without it |
| `-parse_seqids`, `-hash_index`, `-randomize`, `-mask_data`, `-mask_id`, `-mask_desc`, `-gi_mask`, `-gi_mask_name`, `-taxid`, `-taxid_map`, `-oid_masks`, `-metadata_output_prefix`, `-logfile`, `-version` | parsed Seq-ids, masks and taxonomy are a second step (P4), the others are not implemented; each is an explicit error that names the option |
| `-blastdb_version 5` | the LMDB file is not written; version 4 is the default of `LOSAT makedb` (makeblastdb's default is 5) |
| `-input_type` other than `fasta`, several input files, compressed input (`.gz`, `.bz2`, `.zst`) | one uncompressed FASTA input; `-in -` reads the standard input |
| Replacing a database of the same name whose alias file the reader of stage DB-1 rejects | the files to remove are found by resolving the alias file as NCBI does; LOSAT stops before it removes anything |

Where NCBI's behaviour is undefined, LOSAT does not reproduce it: a volume file past 4 GiB
(makeblastdb's 32-bit offsets wrap) is an error, and the test whether an input looks like
FASTA does not read past the bytes it sampled.

A database of 1,000 volumes or more is not in FASTA order when it is read, with makeblastdb
and with `LOSAT makedb` alike: NCBI opens volumes in the sorted order of their names, and
the names makeblastdb gives then (`db.0010`, `db.00100`, `db.00101`, ...) do not sort by
number. `LOSAT makedb` writes the same names and LOSAT's reader opens them in NCBI's order.

## Evidence

Design and measurements of 2026-10-12 are kept outside the repository by the campaign
(`dbsearch-20261012`: design, probe tables, decision record). The reader's tests are in
`LOSAT/src/seqdb/tests.rs` with fixtures written by NCBI BLAST+ 2.17.0 under
`LOSAT/tests/seqdb`. The writer's tests are in `LOSAT/src/writedb` (`tests.rs`,
`fasta_source/tests.rs`) with runs of makeblastdb 2.17.0 recorded in the generated tables
beside them.
