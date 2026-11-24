- Database searches: a sequence whose definition lines hold no Seq-id at all (a hand-edited database;
  makeblastdb writes none) leaves the reports as in NCBI BLAST+ 2.17.0, also without a sequence list: it is
  searched, counts in the statistics and takes its place before the `-max_target_seqs` cut, and its hits are
  then dropped instead of the run being refused.
