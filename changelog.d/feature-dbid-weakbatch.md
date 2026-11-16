- `-db` with parsed Seq-ids (all six programs): with `PRE_FETCH_SEQS_LIMIT`, a sequence whose
  Seq-align id another sequence's id can match weakly (`gb|X|` and `gb|X.1|`, `sp||NAME` and
  `sp|P00001.1|NAME`) is now refused only where a query batch that reports it has its prefetch
  off (0, or a value below the batch's number of alignments): NCBI BLAST+ 2.17.0 reads such a
  batch's sequences one by one and can show the other sequence for it, while a batch whose
  prefetch is on loads each sequence by its own id. A positive value at or above a batch's
  alignments, such as `PRE_FETCH_SEQS_LIMIT=1000`, now writes those sequences as NCBI does
  instead of refusing the run. psiblast's PSSM follows the batch of its round.
- psiblast `-db` with parsed Seq-ids and `-num_iterations` other than 1: each round that makes the
  next round's PSSM now checks the sequences beyond `-max_target_seqs` that the PSSM reads by their
  Seq-align ids (those with an alignment below `-inclusion_ethresh`), as `-out_pssm` did, which now
  also skips the others. Before, a sequence there that NCBI BLAST+ 2.17.0 reads
  as another one (an accession without a version that the lookup finds as a versioned sequence
  first, or a weak match with the prefetch off) gave later rounds that differ from NCBI's.
