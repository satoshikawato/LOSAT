- BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX write outfmt 8 and 9 (Seq-annot, ASN.1 text and
  binary) with a `-db` title or a sequence id holding a byte outside the visible ASCII range (such
  as `café` or a `-parse_deflines` id with DEL) as NCBI BLAST+ does: `#` for each such byte and
  NCBI's error message for each byte of each query's document, in the same place among the output
  when standard output and standard error go to one file; LOSAT refused these before.
- psiblast `-out_pssm` (with `-save_pssm_after_last_round`, `-save_each_pssm` and rounds) writes
  a query title with such a byte as NCBI does (`#` in the PSSM file and NCBI's error message at
  each write of the file); LOSAT refused it before.
