- `losat blastdbcmd`: `-entry`, `-entry_batch`, `-get_dups` and `-target_only` look ids up in a version 5
  database built with `-parse_seqids` (NCBI's distributed databases such as swissprot and pdbaa) as NCBI BLAST+
  2.17.0 blastdbcmd does, byte for byte: accessions, names and local ids in the LMDB table `acc2oid` (all the
  sequences of a key in LMDB's order, several volumes, alias files of several LMDB files, partial aliases), GIs
  in the numeric ISAM, with NCBI's `Entry not found`, `Skipped` and `Entry or entries not found` lines. An accession
  on more than about 500 sequences, on which NCBI stops with `Accessions to Oids lookup error`, is looked up like
  any other (`PD-LOSAT-NCBI-DEFECTS` exception 6). The lookup in a version 4 database built with `-parse_seqids`
  stays refused.
- Every program that opens a version 5 database now reads its LMDB file there, as NCBI does: a file that is not an
  LMDB file, an empty one, or one that lists other volumes gives NCBI's `BLAST Database error` and exit status 2, and
  a truncated one (on which NCBI crashes) is refused.
