- Match NCBI BLAST+ 32-bit signed protein sequence length handling across `seqdb`, `blastdbcmd`,
  `segmasker`, `blastp`, `blastx`, `psiblast`, and `makeblastdb -input_type blastdb` for protein
  databases containing sequences of `2^31` or more residues or `0` residues, while keeping
  nucleotide lengths widened.
