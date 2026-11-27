- `losat blast_formatter`: writes the pairwise report (`-outfmt 0`, its default, so
  `blast_formatter -archive results.asn` works without `-outfmt`) and the query-anchored reports
  (`-outfmt` 1 to 4), also with `-html`, from BLAST
  archives of `blastn` (all tasks), `blastp` (all tasks, gapped and ungapped), `psiblast` (all
  rounds, `-in_msa`, `-in_pssm`), `deltablast`, PHI-BLAST of a protein query and `rpsblast`, for
  `-db` and `-subject` searches and batch archives, with `-num_descriptions`, `-num_alignments`,
  `-max_target_seqs`, `-line_length`, `-sorthits`, `-sorthsps`, `-show_gis`, `-parse_deflines` and
  `-out`, matching NCBI BLAST+ 2.17.0 bytes, messages and exit statuses (apart from the LOSAT
  version line and error prefixes). Archives of translated searches and `-subject` archives with
  `BL2SEQ_LEGACY` set remain refused by name in these formats.
- LOSAT Web: the **Format Archive** tab offers `-outfmt` 0 (now the selected default, as in
  `blast_formatter`) to 4 beside 8, 9 and 12.
