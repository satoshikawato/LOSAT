# Product Decision: Definition line of a record in the BLAST archive (outfmt 11)

- Decision ID: `PD-LOSAT-ARCHIVE-DEFLINE`
- Version: 1.0
- Date: 2026-10-25
- Status: Recorded for the v0.2.0 release. Owner-delegated (standing instruction 2026-09-29,
  restated 2026-10-07): the difference was found by the release inventory and was written
  only in code comments; the Owner confirms it when the pull request is made ready.

## Scope

`-outfmt 11` (the BLAST archive, ASN.1) of `blastn`, `blastp` and `tblastn`, the three
programs whose `run_local` entry takes parsed records and says so in a code comment.
`blastx` writes the same fallback (see below) but has no such entry; `tblastx` does not write
outfmt 11 (README.md:64). This record covers one thing: the definition line that the archive
stores for a query or subject record. Everything else in
the archive is under the root [`AGENTS.md`](../../AGENTS.md) bit-perfect rule.

## NCBI behaviour

`CFastaReader` stores the raw definition line of each record in the Bioseq as a `CFastaReader`
user object, `DefLine`, through `NStr::PrintableString`
(`objtools/readers/fasta.cpp:597-601`). The archive therefore holds the line as read, with its
extra spaces and control bytes.

## LOSAT behaviour

The reader keeps the line as read when it is configured to (`ReaderConfig::with_defline(true)`,
`LOSAT/src/blastinput/fasta_reader/mod.rs:94`; `FastaRecord::archive_defline`, `:259`), and
the archive writes it (`LOSAT/src/report/archive/bioseq.rs:3-11,166`). Every command-line
run of `-outfmt 11` reads its files that way (`LOSAT/src/algorithm/blastn/blast_engine/run.rs:6967,7136`,
`LOSAT/src/algorithm/blastp/blast_engine.rs:5493,5636`,
`LOSAT/src/algorithm/tblastn/args.rs:902,1224`; BLASTX: `LOSAT/src/algorithm/blastx/input.rs:697`), so the command line equals NCBI, which the
frozen archive rows check (`LOSAT/tests/outfmt0_manifest.tsv` rows `sa.arch.bn.*` and
`sa.arch.tbn.*`; `LOSAT/tests/db_search_manifest.tsv` rows `db.sa.arch.bn.*`,
`db.sa.arch.tbn.*`).

A record that has no kept definition line gets `>` followed by its title
(`report::archive::defline_of`, `LOSAT/src/report/archive/mod.rs:64-74`; BLASTX has its own
copy, `LOSAT/src/algorithm/blastx/archive.rs:219-234`). The title is the line after `>` up to
the first control byte, trimmed (`bioseq.rs:26`), so the archive differs from NCBI's for a
definition line with extra spaces or control bytes. The code comments say so at
`LOSAT/src/algorithm/blastn/blast_engine/run.rs:7439-7441` and
`LOSAT/src/algorithm/blastp/blast_engine.rs:5765-5766`; `LOSAT/src/algorithm/tblastn/args.rs:1390-1392`
states the fallback without the comparison.

### Exact condition

Both must hold:

1. `-outfmt 11`, and
2. the records reach the search without a kept definition line. The command line never does
   this. The entries that take already parsed records can: the library functions
   `run_local_blastn`, `run_local_blastp` and `run_local_tblastn`
   (`LOSAT/src/api/local_blast.rs:40,52,65`), which LOSAT Web's adapter uses, and records made
   by `FastaRecord::new` (`fasta_reader/mod.rs:182`). A caller that reads its records with
   `with_defline(true)` gets NCBI's line (`LOSAT/tests/run_local_blastn.rs:123`). ABI v1
   reads only formats 0, 6 and 7 (`LOSAT/src/web_api/v1_blastn.rs:84-92`), so it never writes
   an archive.

The code comments name only the stored definition line as different, and this record
claims no other difference: for a definition line without extra spaces or control bytes,
`>` and the title are the line as read.

## Release-notes sentence

In `-outfmt 11`, a record that the caller handed to the library without a kept definition line
is written with `>` and its title, which differs from NCBI BLAST+ 2.17.0 for a definition line
with extra spaces or control bytes; the command line keeps the line as read and matches NCBI.
