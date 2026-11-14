- psiblast runs `-seg` (other than `no`) and `-soft_masking true` with `-num_iterations` other than 1
  and with `-in_pssm` as NCBI BLAST+ 2.17.0 does: every round's PSSM search leaves the masked letters'
  words out of its lookup table (also with `-lcase_masking`), and a hard mask's X letters reach the
  search's Karlin-Altschul values (`kbp_std`: the gap trigger, the ASCII PSSM's `Standard Ungapped`
  line) while the PSSM is made of the query's letters unmasked, so the later rounds' alignments and
  PSSM files are NCBI's (line RG, Q4b). They were refused.
