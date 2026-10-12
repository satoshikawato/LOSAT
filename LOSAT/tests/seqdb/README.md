# Fixtures of the database reader (`LOSAT/src/seqdb`)

Small NCBI BLAST databases for the unit tests in `LOSAT/src/seqdb/tests.rs`.

- `prot.faa`, `nucl.fna`: the FASTA files the databases were built from.
- `db/`: databases written by makeblastdb 2.17.0 (versions 4 and 5, single and in several
  volumes, with and without `-parse_seqids`), and alias files written by hand or by
  blastdb_aliastool. Only the files LOSAT reads are kept (`.pin .psq .phr`, `.nin .nsq .nhr`,
  `.pal`, `.nal`).
- `expected/NAME.info.tsv`, `expected/NAME.entries.tsv`: what NCBI's blastdbcmd reads back
  from each database (`-info`; `-entry all -outfmt "%o %l %i %s %t"`).
- `expected/date_pairs.tsv`: the date blastdbcmd reports for a database of two volumes with
  the given date strings (first, second, reported or `ERROR`).

`make_fixtures.py` rebuilds everything. NCBI BLAST+ is the comparison oracle only; the tests
read the committed files and run no NCBI program. A rebuild changes the dates stored in the
index files and in `expected/`.
