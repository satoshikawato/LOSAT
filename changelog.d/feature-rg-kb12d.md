- Pairwise output (outfmt 0, also with `-html`) of blastn, blastp, blastx, tblastn and tblastx:
  where the query's id and the subject's id name the same sequence (for example a query
  `gb|KCW91010|` and a database sequence `gb|KCW91010.1|`, with `-parse_deflines`), a query mask
  that lies outside the query's aligned range but inside the subject's now shows in lower case at
  the columns of the subject's range, and a subject mask likewise at the columns of the query's
  range, as NCBI BLAST+ 2.17.0 shows them. Such runs used to print these letters in upper case.
