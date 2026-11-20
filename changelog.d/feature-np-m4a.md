- `losat makeblastdb` reads `-mask_data`, `-mask_id` and `-mask_desc`: the masks of dustmasker, segmasker, windowmasker and
  convert2blastmask (ASN.1 text or binary) become the mask column of the database, in format version 5 and 4, with or
  without `-parse_seqids`, so that `blastdbcmd -mask_sequence_with` and `%m` of NCBI BLAST+ 2.17.0 and of `losat blastdbcmd`
  read them. Its log, messages and exit statuses are makeblastdb's, and the library's `run_makeblastdb` takes the same options.
  A mask reaches the sequences that NCBI's Seq-id handle matching gives it to, every Seq-loc alternative is read, and
  several Seq-locs of one sequence are combined as makeblastdb combines them (empty and whole ranges included). Only
  ranges whose order NCBI takes from memory addresses (two Seq-ids of one sequence with the same packed value writing
  different ranges on one strand) are refused, by name, at the end of a build that would otherwise succeed (makeblastdb's
  later errors come first).
  `-gi_mask` and `-gi_mask_name` are still refused.
