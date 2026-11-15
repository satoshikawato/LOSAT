- `-seqidlist` and `-negative_seqidlist` with `gnl|BL_ORD_ID|n` on a database list or a multi-volume database whose accession
  table holds such ids (made with `makeblastdb -input_type blastdb` from a database without `-parse_seqids`): NCBI BLAST+
  finds the sequence by the number stored in the table, then tests its definition line with the first OID of its volume
  added, so a sequence after the first volume keeps no line, and NCBI leaves its hits out of the report without a word.
  LOSAT wrote those hits (or, with a negative list, the hits NCBI leaves out). It now tests the line as NCBI does and
  refuses such a search, naming the sequence, as it does for other seqid lists that keep a sequence without a line; the
  sequences of the first volume are searched as before.
