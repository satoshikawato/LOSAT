- `-db` searches of `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` with `-outfmt` 8, 9, 11 and 12
  now run on a database whose stored titles hold bytes outside the printable ASCII range, also with parsed
  ids: the output and NCBI BLAST+'s error message for each such byte are NCBI's, where LOSAT refused them
  before. In any format of all eight search programs (including `psiblast`, `rpsblast` and `rpstblastn`),
  a query batch that ends at a header NCBI cannot read now gives NCBI's messages for the titles read
  before it, then NCBI's error. Otherwise `psiblast`, `rpsblast`, `rpstblastn`, the other formats, searches
  with a seqid, GI or taxid list, alias files with a list or an OID mask, stored Seq-ids with such a byte,
  and a query batch that also gives NCBI's warning of a UTF8String header are still refused.
