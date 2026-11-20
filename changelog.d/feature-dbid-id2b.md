- `blastn -db` (megablast, blastn, dc-megablast), `tblastn -db` and `tblastx -db` search databases built
  with `makeblastdb -parse_seqids` (versions 4 and 5, several definition lines, GIs, taxids) and write NCBI's
  subject ids and titles in every format (0, 5, 6, 7, 8, 9, 11, 12, 13 to 16, and SAM for blastn), with the
  same refusals as `blastp -db` and `blastx -db`. psiblast still refuses such a database.
