- `blastn`, `blastp` and `blastx -db` with `-num_threads` above 1: a volume whose sequence file is empty or missing, after the
  point where NCBI BLAST+'s scan stops (a sequence of negative length in a damaged index, or a sequence at offset 0 of an
  empty `.nsq`), no longer makes LOSAT refuse or fail; the search completes as NCBI BLAST+'s does.
- `blastp` and `blastx -db` with a list (such as `-gilist` or `-taxids`) and `-num_threads` above 1: a later volume whose
  `.psq` is missing, which NCBI BLAST+'s scan reaches after the last sequence the list includes, now ends the search with
  NCBI BLAST+'s error (`Database memory map file error`, exit status 3) where LOSAT completed; where the list includes no
  sequence, that error replaces the warning of an empty database, as in NCBI BLAST+.
- An empty (0-byte) sequence file: NCBI BLAST+ reads nothing from it for a protein whose length in the index is 0 or less
  or for a nucleotide that starts at offset 0, and LOSAT no longer refuses such sequences: `blastp` and `blastx` complete
  where the scan stops at such a protein or the list leaves out the others, `blastn -num_threads 1` skips such a
  nucleotide, and `blastdbcmd -entry` writes nothing for such a sequence (`-outfmt %s`: NCBI BLAST+'s error `could not get
  sequence or range.`), as NCBI BLAST+ does. `-outfmt %l` of such a nucleotide stays refused (NCBI BLAST+ crashes).
- A database with a list or an OID range whose nucleotide index gives a sequence a negative approximate length, or a
  sequence of 2^31 bases or more: its total, longest and shortest length (`blastdbcmd -info` and the search's statistics)
  now follow NCBI BLAST+'s 32-bit arithmetic.
- `blastn -filtering_db`: the search of the repeats runs as NCBI BLAST+'s does, without threads, at any `-num_threads`. A
  sequence of negative length in a damaged database of repeats no longer stops it at `-num_threads` above 1, so the repeats
  after that sequence mask the queries as in NCBI BLAST+; a database error in it (a sequence whose ambiguity data cannot be
  read) now ends the query batch with NCBI BLAST+'s error and exit status 3, where LOSAT refused or failed.
