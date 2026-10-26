- `blastp -db` and `blastx -db` search databases built with `makeblastdb -parse_seqids` (versions 4 and 5,
  several definition lines, GIs, taxids; NCBI's distributed protein databases) and write NCBI's subject ids
  and titles in every format (0, 5, 6, 7, 8, 9, 11, 12, 13 to 16). LOSAT refuses the subjects whose
  Seq-align ids NCBI's formatters do not find again, matching Seq-align ids that follow each other in a hit
  list (NCBI joins them in outfmt 0, 5 and 13 to 16), names from the taxonomy database, and the variables
  `SEQ_ID_PREFER_ACCESSION_OVER_GI` and `SEQ_ID_AVOID_GI` with such a database.
