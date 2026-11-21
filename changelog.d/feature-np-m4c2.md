- `losat blastdbcmd` and `losat blastn` read the GI-indexed masks of a database built with `makeblastdb -gi_mask`
  (the alias file's `MASKLIST`), as NCBI BLAST+ 2.17.0 reads them: `-info`, `-outfmt %m`, `-mask_sequence_with` by id
  or name, `-db_soft_mask` and `-db_hard_mask`, with NCBI's output, messages and exit status. Every search program
  searches such a database (until now LOSAT refused its alias file).
  `blastn` refuses `-db_soft_mask` and `-db_hard_mask` on such a database with `-num_threads` above 1, where NCBI's
  threads stop or crash in some runs. A protein in FASTA with such masks is masked over its whole length, as for masks
  in columns.
