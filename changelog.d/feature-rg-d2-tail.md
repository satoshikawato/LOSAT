- `blastn -db` with `-num_threads` above 1, on an alias file whose `FIRST_OID`-`LAST_OID` range ends inside a later volume of
  a damaged database (an ambiguity offset past the sequence file): LOSAT now refuses exactly where NCBI BLAST+ crashes reading
  past the range's end. How far NCBI BLAST+ reads there depends on the largest database file it has read so far, which LOSAT
  now follows through earlier query batches, `-db_soft_mask` and `-db_hard_mask`, `-filtering_db`, an alias file's id list
  and an alias file's `LENGTH 0`.
- `blastn`, `blastp` and `blastx -db` with `-num_threads` above 1: a volume whose sequence file is empty or missing, after the
  point where NCBI BLAST+'s scan stops (a sequence of negative length in a damaged index, or a first sequence at offset 0 of
  an empty `.nsq`), no longer makes LOSAT refuse or fail; the search completes as NCBI BLAST+'s does.
- A database with a list or an OID range whose damaged nucleotide index gives a sequence a negative approximate length: its
  total and shortest length now follow NCBI BLAST+'s arithmetic.
