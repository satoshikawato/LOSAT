- `losat blastdbcmd`, `losat blastn`, `losat blastp`, `losat blastx`, `losat tblastn` and `losat tblastx` read the
  GI-indexed masks of a database built with `makeblastdb -gi_mask` (the alias file's `MASKLIST`), as NCBI BLAST+
  2.17.0 reads them: `-info`, `-outfmt %m`, `-mask_sequence_with` by id or name, `-db_soft_mask` and `-db_hard_mask`,
  with NCBI's output, messages and exit status. Every search program searches such a database (until now LOSAT refused
  its alias file). With `-db_soft_mask` or `-db_hard_mask` on such a database, the search programs write NCBI's
  `-num_threads 1` output at every thread count, where NCBI's own search programs stop or crash in some runs on more
  than one thread. A protein in FASTA with such masks is masked over its whole length, as for masks in columns.
- `losat makeblastdb -input_type blastdb` with such a database as its one input copies its GI-indexed masks into the
  copy's mask column, with NCBI's algorithm ids, as NCBI's `makeblastdb` does.
- `losat makembindex -db_mask` and `-show_filters` read such masks too, by id or by the names of `-gi_mask_name`, as
  NCBI's makembindex does (on a database of several volumes a new-style index stops with NCBI's
  `unrecognized filter algorithm id`; `-old_style_index true -iformat blastdb` reads them).
