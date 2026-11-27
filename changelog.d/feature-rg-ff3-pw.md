- `losat blast_formatter`: writes the pairwise report (`-outfmt 0`, its default, so
  `blast_formatter -archive results.asn` works without `-outfmt`), also with `-html`, from BLAST
  archives of `blastn` (all tasks), `blastp` (all tasks, gapped and ungapped), `psiblast` (all
  rounds, `-in_msa`, `-in_pssm`), `deltablast`, PHI-BLAST of a protein query and `rpsblast`, for
  `-db` and `-subject` searches and batch archives, with `-num_descriptions`, `-num_alignments`,
  `-max_target_seqs`, `-line_length`, `-sorthits`, `-sorthsps`, `-show_gis`, `-parse_deflines` and
  `-out`, matching NCBI BLAST+ 2.17.0 bytes, messages and exit statuses (apart from the LOSAT
  version line and error prefixes). Archives of translated searches and `-subject` archives with
  `BL2SEQ_LEGACY` set remain refused by name in this format.
