- psiblast runs `-query_loc` with `-num_iterations` other than 1 as NCBI BLAST+ 2.17.0 does: every
  round searches the range, the PSSM of a later round is made of the whole record, and the search
  reads its rows from the range's first letter on (the range's first letter with the record's first
  row), so the later rounds' alignments, identities and PSSM files are NCBI's (line RG, Q4a). It
  was refused.
