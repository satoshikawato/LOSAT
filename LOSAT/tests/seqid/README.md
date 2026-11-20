# Fixtures of the Seq-id module (`LOSAT/src/seqid`)

Small NCBI BLAST databases with parsed Seq-ids and what NCBI BLAST+ 2.17.0 prints from them, for
the unit tests in `LOSAT/src/seqid/tests.rs`. `make_fixtures.py` (its docstring has the usage)
wrote every file outside `fasta/`; NCBI is the comparison oracle only and the tests never run it.

- `fasta/`: the inputs: `prot_ids.fa`, `nucl_ids.fa` and their taxid maps (every Seq-id form that
  makeblastdb writes, several definition lines joined by Ctrl-A, GIs, `lcl|Subject_` ids),
  `prot_extra.fa` (a GI in the second definition line only, 12 titles longer than 60 bytes, an
  empty title), `variants.fa` (one record per further form; the gibbsq, `gnl|PIG|` and `nat|`
  records hold a sequence that no query hits, because NCBI's formatters cannot fetch their
  Bioseqs: `Unknown` titles, and blastp's tabular output stops with an exception), `qprot.fa`
  and `qnucl.fa` (the queries of the reports).
- `db/`: the databases (`makeblastdb -parse_seqids -blastdb_version 4`; the files LOSAT reads).
- `expected/NAME.deflines.tsv`: `blastdbcmd -entry all`, one line per definition line (ordinal id,
  FASTA form and bare id of the shown id, first GI, taxid, sum of the memberships, title).
- `expected/NAME.asn.tsv`: `blastdbcmd -outfmt %d` (the Blast-def-line-set in ASN.1 text), one
  line per stored Seq-id (ordinal id, definition line, position, the Seq-id).
- `expected/18S_sample.ber`, `expected/18S_sample.deflines.tsv`: 26 header records of NCBI's
  distributed database `18S_fungal_sequences` (version 5; memberships, GIs above 2^31, up to 23
  definition lines) and blastdbcmd's lines for them. The test of the whole database is ignored
  unless `LOSAT_SEQID_18S_DB` and `LOSAT_SEQID_18S_DEFLINES` name the database and its dump.
- `expected/CASE.fmt{6,5,8,12,0}.tsv`: the subject id and title fields of blastp and blastn
  reports (`-max_hsps 1`) on `prot_ids`, `nucl_ids`, `prot_extra` and `variants` (tabular
  `sseqid ... staxids`; XML `Hit_id`, `Hit_def`, `Hit_accession`; the Seq-align subject id in
  ASN.1 text and in JSON; outfmt 0 description lines and definition-line blocks, with and
  without `-show_gis`).

Fields use backslash escapes (`\\`, `\t`, `\n`).
