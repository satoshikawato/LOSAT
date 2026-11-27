- A search of a database with tax ids and an empty `taxdb.btd` (where NCBI BLAST+ crashes) is refused at the same place
  in every program: after the results of the query batches before, as NCBI's own error of that batch is placed.
  `blastp`, `psiblast`, `tblastn` and `tblastx` refused it before writing anything, and `blastn` and `blastx` lost the
  results before it in `-outfmt` 5, 15 and 16. The words of the refusal and its exit status (1) are unchanged.
