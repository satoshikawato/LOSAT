- Database search: a database sequence whose id is the id NCBI BLAST+ gives a query of the same
  batch (`lcl|Query_1` for the first query read without `-parse_deflines`, in any letter case,
  or the id of a `-in_pssm` query) is now refused in every format where NCBI shows that query in
  its place: `-outfmt 0` to `5` and `13` to `16`, the tables whose fields show what NCBI finds by
  the id (the default fields outside `blastn`), the organism report of a sequence with a taxid
  that the taxonomy database names, and `psiblast` rounds that make the next PSSM, in all eight
  search programs. LOSAT wrote the sequence's own output there, which differed from NCBI's
  (before: only `-outfmt 1` to `4` of five programs and some `rpsblast` and `rpstblastn`
  formats were refused). `-outfmt 8`, `9`, `11`, `12` and `17`, and tables that read nothing of
  the sequence, still give NCBI's output.
