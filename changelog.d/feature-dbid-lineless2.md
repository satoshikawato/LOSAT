- `-negative_seqidlist` that keeps a sequence while removing all its definition lines (every line of the sequence has
  only a gi or a ti, or a `gnl|BL_ORD_ID|n` of a later volume): every search program now reports such a run as NCBI
  BLAST+ does. The sequence is searched and counts in the statistics, but NCBI excludes it where the search reads it
  again for the alignments (the traceback, the composition-based statistics, the ungapped re-evaluation of a
  nucleotide sequence), so it takes no place among the `-max_target_seqs` hits: the next sequences take them
  (blastx `-ungapped`, which reads no protein sequence again, drops it after the cut, as for `-seqidlist`). LOSAT
  refused such a run where that sequence was to be reported.
