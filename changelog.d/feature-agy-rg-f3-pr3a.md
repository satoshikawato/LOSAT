- `losat blast_formatter`: writes tabular reports (`-outfmt 6`, `7`, `10` and `20`, with default
  or custom field lists and delimiters) from ASN.1 text BLAST archives of `blastn` (all tasks),
  `blastp` (gapped and ungapped), `psiblast` (all iterations, `-in_msa`, `-in_pssm`), `deltablast`,
  PHI-BLAST of a protein query (`psiblast -phi_pattern`) and `rpsblast`, for both `-db` and
  `-subject` searches, batch archives, and `-max_target_seqs`, `-parse_deflines` and `-out`,
  matching NCBI BLAST+ 2.17.0 bytes, messages and exit statuses (apart from the LOSAT version line
  and error prefixes). Archives of translated searches (`blastx`, `tblastn`, `tblastx`,
  `rpstblastn`), queries stored as `seq-loc-list`, and `-subject` with `BL2SEQ_LEGACY` in
  `-outfmt 7` remain refused by name.
