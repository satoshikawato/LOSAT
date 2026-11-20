- `losat blastdbcmd` and `losat blastn` read the GI-indexed masks of a database built with `makeblastdb -gi_mask`
  (the alias file's `MASKLIST`), as NCBI BLAST+ 2.17.0 reads them: `-info`, `-outfmt %m`, `-mask_sequence_with` by id
  or name, `-db_soft_mask` and `-db_hard_mask`, with NCBI's output, messages and exit status. Every search program
  searches such a database (until now LOSAT refused its alias file).
