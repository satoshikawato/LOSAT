- Engine: BLASTN, TBLASTN and TBLASTX delete HSPs by percent identity and minimal hit length
  (NCBI's `Blast_HSPTest`) with the one function that BLASTX uses, instead of three copies of
  their own. No output changes.
