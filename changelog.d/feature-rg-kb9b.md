- `-db` searches with `-outfmt` 8, 9, 11 and 12 now run on a database whose stored titles hold bytes
  outside the printable ASCII range: the output and NCBI BLAST+'s error message for each such byte
  are NCBI's, where LOSAT refused them before. The other formats still refuse such a title.
