- `-db` on a database without parsed Seq-ids, in every search program: the tabular formats (outfmt 6, 7, 10 and
  20) of a sequence whose title gives no word for the subject id, a title of spaces or a nucleotide's empty title,
  now print NCBI BLAST+ 2.17.0's general id, `gnl|BL_ORD_ID|n` for `sseqid` and `sallseqid` and `BL_ORD_ID:n`
  for `sacc`, `saccver` and `sallacc`, with the sequence's OID n (also under `-taxids` or an alias file's
  `FIRST_OID`), where LOSAT refused the run; the library's reported sequences carry `BL_ORD_ID:n` as their id. A
  nucleotide's empty title whose tax id the taxonomy database names prints the organism's first word (`Homo`)
  there, as NCBI makes that title from the tax id (line RG, T5). A table of other fields (`6 qseqid stitle staxid evalue`) is
  written as NCBI writes it, where LOSAT refused it.
