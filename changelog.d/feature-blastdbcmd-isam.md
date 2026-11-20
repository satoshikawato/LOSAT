- `losat blastdbcmd`: `-entry`, `-entry_batch`, `-get_dups` and `-target_only` look ids up in a version 4 database
  built with `-parse_seqids` as NCBI BLAST+ 2.17.0 blastdbcmd does, byte for byte: `SeqDB_SimplifyAccession`, the
  string ISAM of each volume (`gb|` forms, the id itself, a version stripped and then checked against the sequence's
  ids, the FASTA form of the id), GIs and PIGs in their numeric ISAMs, trace ids through the string ISAM, all the
  sequences of a key in the order of the volumes and of the ISAM, with NCBI's `Entry not found`, `Skipped` and
  exception lines (`ID overflows range of specified type.`, the `CObject_id` and `Int8` conversion errors). A string
  ISAM whose header does not fit its files, on which NCBI crashes, is refused.
