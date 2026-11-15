- `-seqidlist` that keeps a sequence while matching none of its definition lines (a pdb molecule named without its
  chain, a local id in other letters, a Swiss-Prot name on a version 4 database, a `gnl|BL_ORD_ID|n` of a later volume
  or a later database of a `-db` list): every search program now reports such a run as NCBI BLAST+ does. The sequence
  is searched and counts in the statistics, and its hits take their places in the hit list before the cut to
  `-max_target_seqs` (not in psiblast, whose alignments are made of the longer preliminary list and cut after), but
  NCBI makes no alignment of them, so they leave every report, the prefetch, the taxonomy database's lookups and
  psiblast's PSSM without a word. LOSAT refused such a list before the
  search, whether or not that sequence would be hit. A `-negative_seqidlist` that leaves a sequence without a line is
  still refused where that sequence is to be reported.
