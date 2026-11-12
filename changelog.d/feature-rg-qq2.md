- psiblast runs every matrix of NCBI's tables (BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90,
  PAM30, PAM70, PAM250) with each of its table's gap costs in all rounds of a FASTA query, as NCBI
  BLAST+ 2.17.0: the PSSM and the later rounds' search, statistics and positives take the matrix of
  the search (line RG, Q2). A PSSM query (`-in_pssm`) keeps BLOSUM62 11/1.
