- A search of a database with tax ids and an empty `taxdb.btd` (where NCBI BLAST+ crashes) is refused at the same place
  in every program: after the results of the query batches before, as NCBI's own error of that batch is placed.
  `blastp`, `tblastn`, `tblastx`, `rpsblast` and `psiblast` with one round stopped before writing those results,
  `blastn`, `blastx` and `rpstblastn` lost them in `-outfmt` 5, 15 and 16, `blastp`, `blastn`, `tblastn`, `tblastx`,
  `rpsblast` and `psiblast` lost the field line of `-outfmt` 20 when the first query batch failed, and every program but
  `psiblast` lost the files of the query batches before with `-outfmt` 13 and 14 and `-out`; where the first query batch
  fails, `blastp`, `blastn`, `blastx`, `tblastn`, `tblastx`, `rpsblast` and `rpstblastn` now write the opening lines of
  `-outfmt` 15 and 16 before the error, the lines that NCBI BLAST+ writes before its first batch (its crash there loses
  them; these programs wrote nothing). The exit status (1) is unchanged.
