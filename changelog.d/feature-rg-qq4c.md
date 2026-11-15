- psiblast runs `-num_iterations` other than 1 on a query that NCBI BLAST+ 2.17.0 splits into chunks
  (19,800 letters and more by default) as NCBI does: every round's preliminary search runs the chunks
  with the matrix, its traceback the whole query with the PSSM (line RG, Q4c). It was refused.
- psiblast and tblastn `-in_pssm` search a PSSM of 19,800 letters and more whole, as NCBI BLAST+ 2.17.0
  does (NCBI never splits a PSSM query). It was refused.
