- XML2 and JSON (`-outfmt 13` to `16`) on a `-db` list of version 5 databases made with and without `makeblastdb
  -parse_seqids`: where an earlier query batch opened the taxonomy database (a damaged `taxdb`, whose message NCBI
  writes with that batch) and a later batch stops at a parsed Seq-id that a database of the list cannot look up,
  every search program but psiblast now writes the earlier batches' documents and NCBI's message in LOSAT's words (`LOSAT database error: DB
  contains no accession info.`) (exit status 2) as NCBI BLAST+ does; LOSAT refused such a run. A stopping batch whose
  other sequences make the taxonomy database write a message after the batches before stays refused, as in the
  other formats.
