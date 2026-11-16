- `psiblast -db` searches databases built with `makeblastdb -parse_seqids` (NCBI's distributed protein
  databases such as swissprot and pdbaa), with NCBI's subject ids and titles in outfmt 0, 6, 7, 10 and the
  tabular fields, and blastp's refusals of the subjects NCBI does not write correctly. psiblast also refuses
  a cut to `-max_target_seqs` where two following hits have Seq-align ids that NCBI counts as one subject
  (NCBI then keeps one more). The refusals of psiblast's database name PSIBLAST, not BLASTP.
