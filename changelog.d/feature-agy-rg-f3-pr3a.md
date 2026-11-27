- `losat blast_formatter`: writes tabular reports (`-outfmt 6`, `7`, `10` and `20`, with default
  or custom field lists and delimiters) from BLAST archives (ASN.1 text, ASN.1 binary or XML) of
  `blastn` (all tasks),
  `blastp` (gapped and ungapped), `psiblast` (all iterations, `-in_msa`, `-in_pssm`), `deltablast`,
  PHI-BLAST of a protein query (`psiblast -phi_pattern`) and `rpsblast`, for both `-db` and
  `-subject` searches, batch archives, and `-max_target_seqs`, `-parse_deflines` and `-out`,
  matching NCBI BLAST+ 2.17.0 bytes, messages and exit statuses (apart from the LOSAT version line
  and error prefixes). Archives of translated searches (`blastx`, `tblastn`, `tblastx`,
  `rpstblastn`), queries stored as `seq-loc-list`, and `-subject` with `BL2SEQ_LEGACY` in
  `-outfmt 7` remain refused by name.
- `losat blast_formatter`: in every output format, an archive of `-subject` sequences (without
  `BL2SEQ_LEGACY`) whose subject sequence does not read (more letters declared than held, a protein
  residue outside NCBIstdaa's 28 codes, a sequence coding that does not convert), or reads but
  declares 4294967295 letters, ends with NCBI BLAST+'s engine error before any output, as NCBI's
  sets up the subject sequences; a subject sequence of 0 letters gets NCBI's warning and the report
  goes on.
- `losat blast_formatter`: in the tabular formats, a query or subject sequence that holds fewer
  letters than it declares (or whose segments end before its length) ends the run with NCBI
  BLAST+'s error at the row that reads past its letters, with the position NCBI's read names.
- `losat blast_formatter`: in the tabular formats, `-subject` sequences with `BL2SEQ_LEGACY` and a
  subject (or query) that declares 4294967295 letters give NCBI BLAST+'s rows: the length as NCBI
  reads it at each row (the declared one, or that of the letters once read), and where NCBI's
  prefetch of the subjects meets such a subject, its rows with only the ids and the scores.
