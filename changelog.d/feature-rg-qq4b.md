- psiblast runs `-seg` (other than `no`) and `-soft_masking true` with `-num_iterations` other than 1
  and with `-in_pssm` as NCBI BLAST+ 2.17.0 does, with one query for all the rounds: round 1 masks
  the record's letters (SEG and, with `-lcase_masking`, the lower case), each later round runs SEG
  again on the letters that the rounds before masked (X, none with `-soft_masking true`) without the
  lower-case mask, and every round's PSSM lookup table leaves out the words of its masked letters;
  the reports show the record's letters and count the identities on them (line RG, Q4b). They were
  refused.
