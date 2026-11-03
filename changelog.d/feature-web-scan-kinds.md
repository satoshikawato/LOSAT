- LOSAT Web: the input's record table is read with the program's own FASTA reader (NCBI's, with the nucleotide or
  the protein flags), so it counts what the search reads: a protein-like record in a BLASTN query no longer stops the
  run, a sequence without a defline is one record without a title, and a query of white space only runs as the
  CLI's does ("Query is Empty!"). A gap line (`>?`) is refused when the input is read.
