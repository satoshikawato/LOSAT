- `losat makeblastdb` reads `-mask_data`, `-mask_id` and `-mask_desc`: the masks of dustmasker, segmasker, windowmasker and
  convert2blastmask (ASN.1 text or binary) become the mask column of the database, in format version 5 and 4, with or
  without `-parse_seqids`, so that `blastdbcmd -mask_sequence_with` and `%m` of NCBI BLAST+ 2.17.0 and of `losat blastdbcmd`
  read them. Its log, messages and exit statuses are makeblastdb's, and the library's `run_makeblastdb` takes the same options.
  A mask reaches the sequences that NCBI's Seq-id handle matching gives it to; the two cases whose answer depends on
  memory addresses in NCBI (a `giim` id, two ids of one sequence with the same packed value) are refused by name.
  `-gi_mask` and `-gi_mask_name` are still refused.
