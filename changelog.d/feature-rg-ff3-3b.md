- `losat blast_formatter`: the tabular reports (`-outfmt 6`, `7`, `10` and `20`) now also cover
  archives of the translated searches (`blastx`, `tblastn`, `tblastx`, `rpstblastn`, and
  PSI-TBLASTN as NCBI BLAST+'s `blast_formatter -rid` saves it), gapped and ungapped, with their
  frames, translated letters and genetic codes; queries that a `-parse_deflines` search keeps as
  locations (read from the archive's database); and `-subject` archives with `BL2SEQ_LEGACY` in
  `-outfmt 7` (the `# Subject:` lines). Bytes, messages and exit statuses match NCBI BLAST+ 2.17.0
  (apart from the LOSAT version line and error prefixes). Queries given as locations that the
  database does not hold, which NCBI BLAST+ asks its servers for, stay refused by name.
