- A search of a database with tax ids and an empty `taxdb.btd` (where NCBI BLAST+ crashes) is refused at the same place
  in every program: after the results of the query batches before, as NCBI's own error of that batch is placed.
  `blastp`, `tblastn`, `rpsblast` and `psiblast` with one round refused it before writing anything, `tblastx` wrote only
  the heading of `-outfmt` 0 to 4, `blastn`, `blastx` and `rpstblastn` lost the results before it in `-outfmt` 5, 15
  and 16, `blastn` lost the field line of `-outfmt` 20 when the first query batch failed, and every program but
  `psiblast` lost the files of the query batches before with `-outfmt` 13 and 14 and `-out`. The exit status (1) is
  unchanged.
