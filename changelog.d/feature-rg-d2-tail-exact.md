- `blastn`, `blastp` and `blastx -db` with `-num_threads` above 1: a volume whose sequence file is empty or missing, after the
  point where NCBI BLAST+'s scan stops (a sequence of negative length in a damaged index, or a first sequence at offset 0 of
  an empty `.nsq`), no longer makes LOSAT refuse or fail; the search completes as NCBI BLAST+'s does.
- A database with a list or an OID range whose damaged nucleotide index gives a sequence a negative approximate length: its
  total and shortest length now follow NCBI BLAST+'s arithmetic.
- `blastn -filtering_db`: the search of the repeats runs as NCBI BLAST+'s does, without threads, at any `-num_threads`. A
  sequence of negative length in a damaged database of repeats no longer stops it at `-num_threads` above 1, so the repeats
  after that sequence mask the queries as in NCBI BLAST+; a database error in it (a sequence whose ambiguity data cannot be
  read) now ends the query batch with NCBI BLAST+'s error and exit status 3, where LOSAT refused or failed.
